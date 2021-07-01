TITLE Activity-dependent weight adjuster using Van Rossum et al. (2000)

COMMENT

Homeostatic plasticity (activity-dependent weight scaling) according to 
equation (3) in article Van Rossum et al. (2000) - Stable hebbian Learning 
from Spike Timing-Dependent Plasticity.

Usage
-----

- The default tau and beta are very slow as they approximate the biological 
  rates (Van Rossum et al. 2000). Lower this if you want to use the mechanism
  for calibrating synaptic weights rather than study a biologically
  plausible scenario.


Example
-------

Python example:

>> hpwa = h.VanRossumHPWA(soma(0.5))
>> hpwa.scaling = 0
>> hpwa.sensing = 0
>> hpwa.activitytau = 1e3 # sense activity over window of last ~ 1 second
>> hpwa.activitybeta = 1e-1
>> hpwa.set_target_rate(10.0)
>> 
>> h.setpointer(soma(0.5)._ref_v, 'vm', hpwa)
>> for syn in synapse_list:
>>     h.setpointer(syn._ref_gmax_AMPA, 'temp_wref', hpwa)
>>     hpwa.add_wref(1)
>>     h.setpointer(syn._ref_gmax_NMDA, 'temp_wref', hpwa)
>>     hpwa.add_wref(1)
>> 
>> spike_rec = h.NetCon(soma(0.5)._ref_v, hpwa, sec=soma)
>> spike_rec.threshold = 0.0
>> spike_rec.delay = 1
>> spike_rec.weight[0] = 1
>> 
>> 
>> cvode = h.CVode()
>> def enable_sensing():
>>     hpwa.sensing = 1
>>     print("Start sensing at time {}".format(h.t))
>> def enable_scaling():
>>     hpwa.scaling = 1
>>     print("Start scaling at time {}".format(h.t))
>> def add_events():
>>     cvode.event(10.0, enable_sensing)
>>     cvode.event(500.0, enable_scaling)
>> fih = h.FInitializeHandler(add_events)

ENDCOMMENT


VERBATIM
// Definitions for synaptic scaling procs
// void raise_activity_sensor(double time);
// void decay_activity_sensor(double time);
// void update_scale_factor(double time);
// double get_avg_activity();


// Linked list node for storing refs to observed hoc variables
typedef struct node {
    double id;          // user-specified identifier
    double* hoc_ref;    // pointer to hoc variable (weight)
    double initial_val; // initial value of hoc variable
    struct node* next;  // next node in linked list
} WeightRef;

#define WREFS_INH (*((WeightRef**)&(_p_inh_weight_refs)))
#define WREFS_EXC (*((WeightRef**)&(_p_exc_weight_refs)))

static WeightRef* cwref; // temporary pointer

ENDVERBATIM


NEURON {
    POINT_PROCESS VanRossumHPWA2

    POINTER temp_wref       : temporary variable for passing weight reference
    POINTER exc_weight_refs
    POINTER inh_weight_refs
    POINTER vm

    : Parameters
    RANGE scaling, sensing
    RANGE goal_activity     : set using set_target_rate(hz)
    RANGE scaleinhib        : Set to TRUE (1) for I-cell scaling in addition to 
                            : E-cell scaling. Default is off (0).
    RANGE activitytau       : Activity time constant (ms^-1)  
    RANGE activitybeta      : Scaling strength constant (s^-1 Hz^-1)
    RANGE activitygamma     : Scaling update constant (s^-2 Hz^-1)

    RANGE vm_block_th       : lower threshold for depolarization block
    RANGE vm_spike_th       : upper threshold for depolarization block
    RANGE detect_block      : whether mechanism should try to detect depolarization block
    RANGE vm_buffer         : buffer for moving average Vm values
    RANGE vblock_activity   : (kHz) surrogate activity value when in depolarization block

    : Recordable variables
    RANGE scalefactor       : scale factor for initial values of weights
    RANGE activity          : activity value in kHz (due to ms timing)
    RANGE spkcnt
}


DEFINE VBUFSIZE 50 : max number of samples availabe to the ring buffer


PARAMETER {
    : Local storage
    temp_wref = 0        : temporary variable for passing weight reference
    exc_weight_refs = 0  : refs to excitatory weights
    inh_weight_refs = 0  : refs to inhibitory weights

    vm
    vm_block_th = -35   : threshold for depolarization block
    vm_spike_th = -10   : upper threshold for depolarization block
    detect_block = 1
    vblock_activity = 300e-3 : 300 Hz to kHz

    scaling = 0         : enable weight scaling
    sensing = 1         : enable monitoring of activity
    timeout = 5         : (ms) self-update scheduled after this waiting period
    scaleinhib = 1      : Whether or not we should scale I cells as well as E cells
                        : 0 should be default but because we add weight refs explicitly,
                        : you expect them to be scaled.
    
    goal_activity = 10      : target firing rate
    activitytau = 100.0e3     : Activity sensor time constant (ms) (100s in van Rossum et al., 2000)
    activityoneovertau      : Store 1 / tau for faster calculation
    activitybeta = 4.0e-8   : Scaling strength constant (s^-1 Hz^-1) (4e-5 in van Rossum et al., 2000)
    activitygamma = 1.0e-10 : Scaling update constant (s^-2 Hz^-1) (1e-7 in van Rossum et al., 2000)
}

ASSIGNED {
    : State variables
    activity    : activity sensor value ('a' in Van Rossum et al. (2000))
    activity_integral_err : intergral of cell's activity divergence from target
    scalefactor : default scaling factor for this cell's excitatory synapses
    spkcnt      : spike count of post-synaptic cell
    lastupdate  : time of last activity sensor decay / spike udpate
    max_err     : max error value
    max_scale   : max scaling factor
    : Depolarization block detection
    vm_buffer[VBUFSIZE]
    ibuf        : index for ring buffer
    block_detected  : whether depolarization block is dectected
}


INITIAL {
    activity = 0    : Sensor for this cell's recent activity (default 0kHz i.e. cycles per ms)
    spkcnt = 0      : Total spike count
    max_err = 0     : Max error value
    max_scale = 4   : Max scaling factor
    lastupdate = 0  : Time of last activity sensor decay / spike update
    scalefactor = 1.0 : Default scaling factor for this cell's AMPA synapses
    : goal_activity = -1 : Cell's target activity (kHz i.e. cycles per ms)
    activity_integral_err = 0.0 : Integral of cell's activity divergence from target activity
    activityoneovertau = 1.0 / activitytau

    FROM i=0 TO VBUFSIZE-1 {
        vm_buffer[i] = -68
    }
    ibuf = 0
    block_detected = 0

    : schedule an update in case no spikes received
    net_send(timeout, 1)
}

: Constructor, called only once
CONSTRUCTOR {
VERBATIM {

    //_p_exc_weight_refs = (void*)ecalloc(1, sizeof(WeightRef));
    WREFS_EXC = emalloc(sizeof(WeightRef));
    WREFS_EXC->id = 0.0;
    WREFS_EXC->hoc_ref = NULL;
    WREFS_EXC->next = NULL;

    //_p_inh_weight_refs = (void*)ecalloc(1, sizeof(WeightRef));
    WREFS_INH = emalloc(sizeof(WeightRef));
    WREFS_INH->id = 0.0;
    WREFS_INH->hoc_ref = NULL;
    WREFS_INH->next = NULL;
}
ENDVERBATIM
}


DESTRUCTOR {
VERBATIM {
    // Free linked lists
    WeightRef* current = WREFS_EXC;
    WeightRef* next_node;
    while (current != NULL) {
        next_node = current->next;
        free(current);
        current = next_node;
    }

    current = WREFS_INH;
    while (current != NULL) {
        next_node = current->next;
        free(current);
        current = next_node;
    }

}
ENDVERBATIM
}


: Receive spikes from post-synaptic cell or self-scheduled updates.
: 
: @param    w : double
:           NetCon weight, w > 0 indicates spikes from post-synaptic cell.
:
: @param    flag : int
:           Implicit parameter: flag supplied to `net_send(delay, flag)`
:           or `NetCon.event(t_deliver, flag)`.
:               0 = spike event
:               1 = scheduled update
:               2 = turn on weight scaling
:               3 = turn off weight scaling
NET_RECEIVE(w) { LOCAL spike_received, scheduled_update, vm_avg, jj
    INITIAL { w=w }

    spike_received = (flag==0 && w==1)
    scheduled_update = (flag==1)

    : schedule update in case cell doesn't fire
    if (flag==1) {
      net_send(timeout, 1)
    }

    : WARNING - testing weight causes segmentation fault: need to test with ifarg()
:VERBATIM
:    if (ifarg(1)) {
:        fprintf(stderr, "event with weight arg");
:        if (_args[0]==2.0) {
:            scaling = 1;
:            fprintf(stderr, "turned on scaling");
:        }
:        if (_args[0]==3.0) {
:            scaling = 0;
:        }
:    }
:ENDVERBATIM
:    if (flag==2 || w==2) {
:        scaling = 1
:    }
:    if (flag==3 || w==3) {
:        scaling = 0
:    }
 
    : only update on spike or timeout
    if (spike_received || (scheduled_update && (t-lastupdate >= timeout))) {

        : very rudimentary depolarization block detector
        : NOTE: ideally this should be changed to a low-pass filter of Vm
        :       Current approach counts on the fact that events are received in AHP
        :       of spike or during self-scheduled update.
        block_detected = 0
        if (detect_block) {
            : Add vm to moving average buffer if we are not in spike
            if (vm < vm_spike_th) {
                vm_buffer[ibuf] = vm
                ibuf = ibuf + 1
                if (ibuf >= VBUFSIZE) {
                    ibuf = 0
                }
            }
            : Calculate moving average and decide whether we are in depolarization block
            vm_avg = 0
            FROM i=0 TO VBUFSIZE-1 {
              vm_avg = vm_avg + vm_buffer[i] / VBUFSIZE
            }
            if ((vm_avg >= vm_block_th) && (vm_avg < vm_spike_th)) {
                block_detected = 1
                VERBATIM
                fprintf(stderr, "block detect at t = %f ms", t);
                ENDVERBATIM
            }
        }
        
        if (sensing) {
            decay_activity_sensor(t) : Allow activity sensor to decay on every update
            if (spike_received) {
                raise_activity_sensor(t) : Update activity sensor
                spkcnt = spkcnt + 1
            }
            if (block_detected) { : Also update sensor in depolarization block
                vblock_period = 1/vblock_activity
                jj = 0
                while (((t - lastupdate) - (jj * vblock_period)) > 0) {
                    : simulate spiking at rate vblock_activity
                    raise_activity_sensor(t)
                    jj = jj + 1
                }
            }
            lastupdate = t : Store time of last update
        }

        
        
        if (scaling) {
            : Only update if cell is not inhib OR we are scaling all I+E cells
            update_scale_factor(t) : Run synaptic scaling procedure to find scalefactor
            
            if (goal_activity < 0) {
                : If scaling has just been turned on, set goal activity to historical average firing rate This is only meaningful if sensor has had a chance to measure correct activity over a relatively long period of time, so don't call setscaling(1) until at least ~800s.
                
                : goal_activity = get_avg_activity()
                goal_activity = activity : Take current activity sensor value
                : max_err = goal_activity * 0.5; : Error value saturates at +- 50% of goal activity rate
            }
            : Multiply all excitatory wref by scalefactor
            : Scalefactor remains constant when goal_activity reached.
            : This means we have to multiply initial weight with scale factor
            : or we get exponential increase in weight.
 VERBATIM
            WeightRef* current = WREFS_EXC;
            while (current != NULL) {
                if (current->hoc_ref != NULL) {
                    // *(current->hoc_ref) *= scalefactor;
                    *(current->hoc_ref) = current->initial_val * scalefactor;
                }
                current = current->next;
            }
            current = WREFS_INH;
            while (current != NULL) {
                if (current->hoc_ref != NULL) {
                    // *(current->hoc_ref) *= 1 / scalefactor;
                    *(current->hoc_ref) = current->initial_val / scalefactor;
                }
                current = current->next;
            }
 ENDVERBATIM
        }
    }
}


: Set the target spike rate
PROCEDURE set_target_rate(rate) {
    goal_activity = rate * 1e-3 : convert to kHz
}

: Get average spike rate of the cell since t = 0.
FUNCTION get_avg_activity () {
  get_avg_activity = spkcnt / t
}

: Update the cell's activity sensor value, assuming this function has been 
: called at the same time as a spike at time t.
:
: @param  time : double
:         time of current spike in ms
PROCEDURE raise_activity_sensor(time) {
  : max activity of 1 = 1 kHz, rate of increase saturates for higher activity
  activity = activity + (-activity + 1.0) / activitytau
}

: Decay the cell's activity sensor value according to the time since last decay update.
: In algorithm described in van Rossum et al. (2000), this is called every discrete timestep t
: But this procedure is only called on NET_RECEIVE events, so we need to decay
: taking into account the time since the last decay operation.
PROCEDURE decay_activity_sensor(time) {
    activity = activity * exp(-activityoneovertau * (time - lastupdate))
}


: Implements weight scaling according to van Rossum et al. (2000), eq. (3).
: Calculates dw/dt according to PI controller and target activity.
PROCEDURE update_scale_factor(time) {
  LOCAL err, timecorrection
  : Get difference between goal and current activity
  err = goal_activity - activity
  if (block_detected) {
    err = goal_activity - vblock_activity
  }

  : Set scalefactor
  scalefactor = scalefactor + (activitybeta * scalefactor * err + activitygamma * scalefactor * activity_integral_err)

  : Bound scalefactor to max_scale to prevent Inf values
  if (scalefactor > max_scale) {
    scalefactor = max_scale
  }

  : Calculate integral error term between sensor and target activity for next time (t')
  timecorrection = time - lastupdate
  : e.g. If last update was 1ms ago, then the time correction = 1
  : If last update was 0.1ms ago correction = 0.1, so the accumulated error will be much smaller
  : If it's been a long time since the last update, the error will be correspondingly much larger

  activity_integral_err = activity_integral_err + (err * timecorrection)
}


FUNCTION add_wref(excitatory) {
VERBATIM
    // Look for end of linked list and append controlled variable
    WeightRef* current;
    if (_lexcitatory) {
        current = WREFS_EXC;
    } else {
        current = WREFS_INH;
    }
    double i = 0;
    while (current->next != NULL) {
        current = current->next;
        i++;
    }

    current->next = emalloc(sizeof(WeightRef));
    current->next->id = i;
    current->next->hoc_ref = _p_temp_wref;
    current->next->initial_val = (double)*_p_temp_wref;
    current->next->next = NULL;

    _ladd_wref = i; // return identifier so wref can be removed
ENDVERBATIM
}
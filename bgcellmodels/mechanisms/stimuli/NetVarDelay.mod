COMMENT
NetVarDelay is a NetStim that passes on a spike event after a variable delay.
The delays for subsequent spikes are read from a Vector.

One use case is to sample the ISI of a neuron uniformly by feeding back its 
own spikes via a NetVarDelay. This can be used to construct a Phase Response
Curve for the neuron.

Python Example:
--------------

from neuron import h

cell = MyCell()

delays = h.Vector(range(10))
net_delay = h.NetVarDelay()
net_delay.set_delays(delays)

nc_in = h.NetCon(net_delay, cell)
nc_in.delay = 0.05

nc_out = h.NetCon(cell, net_delay, sec=sec)
nc_out.threshold = 0.0
nc_out.delay = 0.05
nc_out.weight = 1.0


cvode = h.CVode()
h.finitialize()
h.run(100)

ENDCOMMENT

NEURON {
    ARTIFICIAL_CELL NetVarDelay
    RANGE tstart
    POINTER ptr
}

PARAMETER {
    tstart = 0 (ms)
}

ASSIGNED {
    index
    delay (ms)
    ptr
}

INITIAL {
    index = 0
}

NET_RECEIVE (w) {

    : flag is zero for external events
    if (flag == 0 && t >= tstart) {
        element() : read next delay from vector and assign to delay
        if (index > 0) {
            net_event(t + delay)
        }
    }

}

DESTRUCTOR {
VERBATIM
    void* vv = (void*)(_p_ptr);  
        if (vv) {
        hoc_obj_unref(*vector_pobj(vv));
    }
ENDVERBATIM
}

: Read next element from vector and assign it to delay
PROCEDURE element() {
VERBATIM    
  { void* vv; int i, size; double* px;
    i = (int)index;
    if (i >= 0) {
        vv = (void*)(_p_ptr);
        if (vv) {
            size = vector_capacity(vv);
            px = vector_vec(vv);
            if (i < size) {
                delay = px[i];
                index += 1.;
            }else{
                index = -1.;
            }
        }else{
            index = -1.;
        }
    }
  }
ENDVERBATIM
}

: Set delays for incoming spikes
PROCEDURE set_delays() {
VERBATIM
    void** pv;
    void* ptmp = NULL;
    if (ifarg(1)) {
        ptmp = vector_arg(1);
        hoc_obj_ref(*vector_pobj(ptmp));
    }
    pv = (void**)(&_p_ptr);
    if (*pv) {
        hoc_obj_unref(*vector_pobj(*pv));
    }
    *pv = ptmp;
ENDVERBATIM
}

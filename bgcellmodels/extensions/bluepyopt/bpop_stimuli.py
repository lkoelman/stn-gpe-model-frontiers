"""
Extensions to bluepyopt.ephys.stimuli
"""

import bluepyopt.ephys as ephys

from bgcellmodels.mechanisms import stimuli # loads VecStim.mod

import logging
logger = logging.getLogger('bpop_ext')

class NrnSpaceClamp(ephys.stimuli.Stimulus):

    """Square pulse current clamp injection"""

    def __init__(self,
                 step_amplitudes=None,
                 step_durations=None,
                 total_duration=None,
                 location=None):
        """
        Constructor
        
        Args:
            step_amplitudes (float):    amplitude (nA)
            step_durations (float):     duration (ms)
            total_duration (float):     total duration of stimulus and its effects (ms)
            location (Location):        stimulus Location
        """

        super(NrnSpaceClamp, self).__init__()
        self.step_amplitudes = step_amplitudes
        self.step_durations = step_durations
        self.location = location
        self.total_duration = total_duration
        self.seclamp = None


    def instantiate(self, sim=None, icell=None):
        """Run stimulus"""

        icomp = self.location.instantiate(sim=sim, icell=icell)
        logger.debug(
            'Adding space clamp to {} with '
            'durations {}, and amplitudes {}'.format(
            str(self.location),
            self.step_durations,
            self.step_amplitudes))

        # Make SEClamp (NEURON space clamp)
        self.seclamp = sim.neuron.h.SEClamp(icomp.x, sec=icomp.sec)
        for i in range(3):
            setattr(self.seclamp, 'amp%d' % (i+1), self.step_amplitudes[i])
            setattr(self.seclamp, 'dur%d' % (i+1), self.step_durations[i])


    def destroy(self, sim=None):
        """Destroy stimulus"""

        self.seclamp = None


    def __str__(self):
        """String representation"""

        return "Square pulse amps {} durations {} totdur {} at {}".format(
            self.step_amplitudes,
            self.step_durations,
            self.total_duration,
            self.location)


class NrnVecStimStimulus(ephys.stimuli.Stimulus):
    """
    NetStimStimulus equivalent with spike times read from Vector.
    """

    def __init__(self,
                 locations=None,
                 total_duration=None,
                 times=None,
                 weight=1):
        """Constructor
        Args:
            location: synapse point process location to connect to
            times (list[float]) : stimulus times
        """

        super(NrnVecStimStimulus, self).__init__()
        if total_duration is None:
            raise ValueError(
                'NrnNetStimStimulus: Need to specify a total duration')
        else:
            self.total_duration = total_duration

        self.locations = locations
        self.times = times
        self.times_vec = None
        self.weight = weight
        self.connections = {}

    def instantiate(self, sim=None, icell=None):
        """Run stimulus"""

        self.times_vec = sim.neuron.h.Vector(self.times)

        for location in self.locations:
            self.connections[location.name] = []
            for synapse in location.instantiate(sim=sim, icell=icell):
                netstim = sim.neuron.h.VecStim()
                netstim.play(self.times_vec)
                netcon = sim.neuron.h.NetCon(netstim, synapse)
                netcon.weight[0] = self.weight

                self.connections[location.name].append((netcon, netstim))

    def destroy(self, sim=None):
        """Destroy stimulus"""

        self.connections = None
        self.times = None
        self.times_vec = None

    def __str__(self):
        """String representation"""

        return "VecStim at %s" % ','.join(
            location
            for location in self.locations) \
            if self.locations is not None else "Netstim"


class NetVarDelayStimulus(ephys.stimuli.Stimulus):
    """
    NetStimStimulus equivalent that delays incoming spike times by a variable
    delay, read from a Vector.

    @see    NetVarDelay.mod
    """

    def __init__(self,
                 delays=None,
                 start_time=None,
                 target_locations=None,
                 source_location=None,
                 source_threshold=None,
                 source_delay=0.1,
                 target_delay=0.1,
                 target_weight=1.0,
                 total_duration=None):
        """
        Constructor
        
        @param  location: 
                synapse point process location to connect to

        @param  source_location : NrnSectionCompLocation
                NEURON segment where voltage will be monitored for spikes

        @param  delays : (list[float])
                Delays for incoming spikes
        """

        super(NetVarDelayStimulus, self).__init__()
        if total_duration is None:
            raise ValueError(
                'NrnNetStimStimulus: Need to specify a total duration')
        else:
            self.total_duration = total_duration

        self.target_locations = target_locations
        self.locations = target_locations # compatibility with NetStimStimulus
        self.source_location = source_location
        self.source_threshold = source_threshold
        self.source_delay = source_delay
        self.target_delay = target_delay
        self.delays = delays
        self.start_time = start_time
        self.delays_vec = None
        self.target_weight = target_weight
        self.connections = {}

    def instantiate(self, sim=None, icell=None):
        """Run stimulus"""

        self.delays_vec = sim.neuron.h.Vector(self.delays)

        # Stimulator object
        self.stim = net_delay = sim.neuron.h.NetVarDelay()
        net_delay.tstart = self.start_time
        net_delay.set_delays(self.delays_vec)

        # Connection to source of events
        src_seg = self.source_location.instantiate(sim=sim, icell=icell)
        self.source_netcon = nc = sim.neuron.h.NetCon(
            src_seg._ref_v, net_delay, sec=src_seg.sec)
        nc.threshold    = self.source_threshold
        nc.delay        = self.source_delay
        nc.weight[0]    = 1.0 # not used

        # Connection to targets (event consumers, e.g. synapses)
        for location in self.target_locations:
            self.connections[location.name] = []
            for synapse in location.instantiate(sim=sim, icell=icell):
                netcon = sim.neuron.h.NetCon(self.stim, synapse)
                netcon.delay     = self.target_delay
                netcon.weight[0] = self.target_weight

                self.connections[location.name].append((netcon, self.stim))

    def destroy(self, sim=None):
        """Destroy stimulus"""

        self.connections = None
        self.delays = None
        self.delays_vec = None

    def __str__(self):
        """String representation"""

        return "NetVarDelayStimulus at %s" % ','.join(
            location for location in self.target_locations) \
            if self.target_locations is not None else "Netstim"
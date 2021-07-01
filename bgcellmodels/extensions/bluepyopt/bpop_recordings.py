"""
Extensions to bluepyopt.ephys.recordings
"""

import bluepyopt.ephys as ephys
from bluepyopt.ephys.responses import TimeVoltageResponse

import logging
logger = logging.getLogger('bpop_ext')


class NetStimRecording(ephys.recordings.Recording):

    """Response to stimulus"""

    def __init__(
            self,
            name=None,
            netstim=None):
        """
        Constructor
        Args:
            name (str): name of this object
            netstim (NrnNetStimStimulus) : NetStim to record
        """

        super(NetStimRecording, self).__init__(name=name)

        # assert len(netstim.locations) == 1

        self.netstim = netstim
        self.location = netstim.locations[0]
        self.tvector = None
        self.instantiated = False


    @property
    def response(self):
        """Return recording response"""

        if not self.instantiated:
            return None

        spike_times = self.tvector.to_python()

        return TimeVoltageResponse(self.name,
                                   spike_times,
                                   [1.0 for t in spike_times])

    def instantiate(self, sim=None, icell=None):
        """
        Instantiate recording

        @pre    NetStims are instantiated and thus have attribute 'connections'
        """

        self.tvector = sim.neuron.h.Vector()

        # Alternative
        # netstim = self.netstim.connections[self.location.name][1]
        # self.rec_netcon = sim.neuron.h.NetCon(netstim, None)
        # self.rec_netcon.record(self.tvector)

        # Re-use existing netcon for recording (ensures delay is same)
        netcon = self.netstim.connections[self.location.name][0][0]
        netcon.record(self.tvector)

        self.instantiated = True

    def destroy(self, sim=None):
        """Destroy recording"""
        self.tvector = None
        self.instantiated = False

    def __str__(self):
        """String representation"""

        return '%s: recording %s' % (self.name, self.netstim)
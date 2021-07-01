"""
Test and run Ephys protocols for Gunay (2008) GPe cell model and reduced version

@author Lucas Koelman

@date   12/09/2017
"""

# BluePyOpt modules
import bluepyopt.ephys as ephys

# Our custom BluePyOpt modules
from bgcellmodels.extensions.bluepyopt import (
    bpop_cellmodels, feature_factory, protocol_analysis as anls_proto)

from bgcellmodels.extensions.bluepyopt.bpop_protocol_ext import (
    BpopProtocolWrapper, PROTOCOL_WRAPPERS)

from bgcellmodels.models.GPe.Gunay2008 import mechanisms, gunay_model
from bgcellmodels.models.GPe.Gunay2008.optimization import (
    gpe_protocols_prc as prc_protos) # load to register

from bgcellmodels.common.stimprotocols import StimProtocol


################################################################################
# OPTIMIZATION EXPERIMENTS
################################################################################


def test_protocol(stim_proto, model, export_locals=True):
    """
    Test stimulation protocol applied to full cell model.

    @param  stim_proto : StimProtocol
            Enum instance identifying stimulation protocol

    @param  model_type : StnModel or str
            Type of STN model to use ('full', 'reduced', StnModel.Gillies2005)

    USAGE
    -----

    >>> test_protocol(StimProtocol.MIN_SYN_BURST, 'full')
    
    """
    

    # instantiate protocol
    proto = BpopProtocolWrapper.make(stim_proto)

    # Get protocol mechanisms that need to be isntantiated
    proto_mechs, proto_params = proto.get_mechs_params()


    if isinstance(model, str) and model.endswith('.pkl'):
        cellmodel = bpop_cellmodels.PickledCellModel(
                        name='GpeFromPickle',
                        cell_pkl_file=model,
                        mechs=proto_mechs,
                        params=proto_params)

        nrnsim = ephys.simulators.NrnSimulator(dt=0.025, cvode_active=False)
    
    elif model == 'full':
        model, nrnsim = gunay_model.create_cell(
            model=gunay_model.MODEL_GUNAY2008_AXONLESS)
    
    else:
        model, nrnsim = gunay_model.create_cell(model=model)

    # Apply protocol and simulate
    responses = proto.ephys_protocol.run(
                    cell_model      = cellmodel, 
                    param_values    = {},
                    sim             = nrnsim,
                    isolate         = False) # allows us to query cell with h.allsec()

    protos_responses = {proto.ephys_protocol.name: responses}

    # Get features for each stimulation protocol
    # - features are defined based on the response to each stimulation protocol
    #   that we want to capture
    protos_features = feature_factory.make_opt_features([proto])

    # Calculate target values from full model responses
    feature_factory.calc_feature_targets(protos_features, protos_responses,
                                         raise_check=False)

    # Plot protocol responses
    # anls_proto.plot_proto_responses(protos_responses)

    ## Plot protocol responses
    anls_proto.plot_proto_responses(protos_responses)

    if export_locals:
        globals().update(locals())


if __name__ == '__main__':

    pkl_file = ('/home/luye/cloudstore_m/simdata/GunayEdgerton2008_reduced'
            '/gpe-cell_gunay2008_reduce-BushSejnowski_dL-1.pkl')

    test_protocol(StimProtocol.PRC_SYN_EXC_DIST, pkl_file)
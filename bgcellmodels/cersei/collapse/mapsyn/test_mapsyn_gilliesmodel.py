"""
Test synapse mapping procedures on Gillies STN Model.
"""
# print date and time of script execution
import os
gillies_model_dir = '/home/luye/workspace/bgcellmodels/GilliesWillshaw'
os.chdir(gillies_model_dir)

# Python standard library
import pprint
pp = pprint.PrettyPrinter(indent=2)
from bgcellmodels.common import logutils

# BluePyOpt
import bluepyopt.ephys as ephys

# Custom BluePyOpt modules
from cersei_cellmodel import StnCellReduced
from optimize.bpop_protocols_stn import BpopProtocolWrapper

# Physiology parameters
from evalmodel.proto_common import StimProtocol
SP = StimProtocol

# Adjust verbosity of loggers
logutils.setLogLevel('quiet', ['marasco', 'folding', 'redops', 
                               'bluepyopt.ephys.parameters', 
                               'bluepyopt.ephys.recordings'])


def test_mapsyn_excitatory_inhibitory(export_locals=False):
    # ## Create Protocols

    # Protocols to use for optimisation
    opt_proto = SP.SYN_BACKGROUND_LOW
    proto_kwargs = { # SETPARAM: extra keyword arguments for validation protocol
        'impl_proto': opt_proto,
        'base_seed': 8,
        'num_syn_gpe': 12,
    }

    stimprotos_wrappers = {
        SP.SYN_BACKGROUND_LOW: BpopProtocolWrapper.make(opt_proto, **proto_kwargs)
    }

    proto_wrappers = stimprotos_wrappers.values()
    opt_stim_protocols = stimprotos_wrappers.keys()
    ephys_protos = [p.ephys_protocol for p in proto_wrappers]

    # Collect al frozen mechanisms and parameters required for protocols to work
    proto_mechs, proto_params = BpopProtocolWrapper.all_mechs_params(proto_wrappers)


    # ## Run Reduced Model

    # Create reduced model and get parameters
    red_model = StnCellReduced(
                    reduction_method='BushSejnowski',
                    name='StnFolded',
                    mechs=proto_mechs,
                    params=proto_params)

    nrnsim = ephys.simulators.NrnSimulator(dt=0.025, cvode_active=False)

    # Simulate protocols
    red_responses = {}
    for e_proto in ephys_protos:
        
        # Make sure recording functions are executes
        e_proto.record_contained_traces = True
        
        # NOTE: isolate=False only if model not previously build
        red_responses[e_proto.name] = e_proto.run(
                                            cell_model      = red_model, 
                                            param_values    = {},
                                            sim             = nrnsim,
                                            isolate         = False)

    if export_locals:
        globals().update(locals())


if __name__ == '__main__':
    test_mapsyn_excitatory_inhibitory(export_locals=True)
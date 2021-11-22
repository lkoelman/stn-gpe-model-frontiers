"""
Analysis of STN model optimization using BluePyOpt.

@author Lucas Koelman

@date   6/10/2017

@see    based on scripts:
            https://github.com/BlueBrain/BluePyOpt/blob/master/examples/l5pc/l5pc_analysis.py
            https://github.com/BlueBrain/BluePyOpt/blob/master/examples/l5pc/opt_l5pc.py
"""

import pickle

# Scipy
from matplotlib import pyplot as plt

import bluepyopt.ephys as ephys


def plot_responses(responses, show=True, plot_kws=None):
    """
    Plot responses of a single protocol

    @param  responses : dict[str, responses.TimeVoltageResponse]
            Responses returned by run_protocol(...)
    """
        
    fig, axes = plt.subplots(len(responses))
    try:
        iter(axes)
    except TypeError:
        axes = [axes]

    if plot_kws is None:
        plot_kws = {}

    for index, (resp_name, response) in enumerate(sorted(responses.items())):

        if resp_name.endswith('_times'):
            axes[index].plot(response['time'], response['voltage'],
                             marker='|', linestyle='', snap=True)
        else:
            axes[index].plot(response['time'], response['voltage'],
                             label=resp_name, **plot_kws)
            axes[index].set_title(resp_name)
    
    fig.tight_layout()
    if show:
        plt.show(block=False)
    
    return fig, axes

def plot_proto_responses(proto_responses, show=True):
    """
    Plot responses of multiple protocols

    @param  proto_responses : dict[str: protocol_name, dict: responses]
    """
    figs = []
    for proto_name, responses in proto_responses.items():
        fig, axes = plot_responses(responses, show=False)
        plt.suptitle("Responses of {}".format(proto_name))
        figs.append(fig)

    # n_total_resp = sum((len(responses) for responses in proto_responses.values()))    
    
    # fig, axes = plt.subplots(n_total_resp)
    # try:
    #   iter(axes)
    # except TypeError:
    #   axes = [axes]
    
    # index = 0
    # for proto_name, responses in proto_responses.items():
    #   for resp_name, response in sorted(responses.items()):
    #       axes[index].plot(response['time'], response['voltage'], label=resp_name)
    #       axes[index].set_title(resp_name)
    #       index += 1
        
    #   fig.tight_layout()
    
    if show:
      plt.show(block=False)

    return figs


def save_proto_responses(responses, filepath):
    """
    Save protocol responses to file.
    """

    # Save to file
    with open(filepath, 'w') as recfile:
        pickle.dump(responses, recfile)

    print("Saved responses to file {}".format(filepath))


def load_proto_responses(filepath):
    """
    Load protocol responses from pickle file.

    @return     dictionary {protocol: responses}
    """
    with open(filepath, 'r') as recfile:
        responses = pickle.load(recfile)
        return responses


def run_proto_responses(cell_model, ephys_protocols, param_values=None, isolate=True):
    """
    Run protocols using given cell model and return responses,
    indexed by protocol.name.

    @return     dict<str, TimeVoltageResponse> {proto_name: proto_response}
    """
    nrnsim = ephys.simulators.NrnSimulator(dt=0.025, cvode_active=False)

    if param_values is None:
        param_values = {}

    # Run each protocol and get its responses
    all_responses = {}
    for e_proto in ephys_protocols:

        response = e_proto.run(
                        cell_model      = cell_model, 
                        param_values    = param_values,
                        sim             = nrnsim,
                        isolate         = isolate)

        all_responses[e_proto.name] = response

    return all_responses
# -*- coding: utf-8 -*-
"""
Basal Ganglia network model consisting of morphologically detailed
cell models for the major cell types.

Parameterized model construction based on configuration file / dictionary.

@author     Lucas Koelman

@date       31/05/2018

@see        PyNN manual for building networks:
                http://neuralensemble.org/docs/PyNN/building_networks.html
            PyNN examples of networks:
                https://github.com/NeuralEnsemble/PyNN/tree/master/examples

USAGE
-----

To run using MPI, use the `mpirun` or `mpiexec` command like:

`mpirun -n 6 python model_parameterized.py --scale 0.5 --dur 3000 --dnorm --no-lfp --seed 888 --transient-period 0.0 --write-interval 3000 -id CALDNORM --outdir ~/storage/BBB_LuNetStnGpe --config configs/StnGpe_template_syn-V8.json 2>&1 | tee CALDNORM.log`

To run from an IPython shell, use the %run magic function like:

`%run model_parameterized.py --scale 0.5 --dd --dur 100 --seed 888 --transient-period 0.0 --write-interval 1000 --no-gui -id test1 --outdir ~/storage --config configs/syn-V18_template.json`


NOTES
-----

- It may look like some imports are not used but they may be called dynamically
  using eval() based on the config file.

"""
from __future__ import print_function
import time
import os
from datetime import datetime

import numpy as np

# MPI support
from mpi4py import MPI
comm = MPI.COMM_WORLD
mpi_size = comm.Get_size() # number of processes
mpi_rank = comm.Get_rank() # rank of current process
WITH_MPI = mpi_size > 1

# PyNN library
import pyNN.neuron as sim
from pyNN import space

from pyNN.random import RandomDistribution
from pyNN.utility import init_logging # connection_plot is bugged
import neo.io

# Custom PyNN extensions
from bgcellmodels.extensions.pynn.synapses import (
    GluSynapse, GabaSynapse, GabaSynTmHill, GABAAsynTM,
    GabaSynTm2, NativeMultiSynapse)
from bgcellmodels.extensions.pynn.utility import connection_plot
from bgcellmodels.extensions.pynn.populations import Population
import bgcellmodels.extensions.pynn.spiketrains as spikegen

# Monkey-patching of pyNN.neuron.Population class
# from bgcellmodels.extensions.pynn.recording import TraceSpecRecorder
# sim.Population._recorder_class = TraceSpecRecorder

# Custom NEURON mechanisms
from bgcellmodels.mechanisms import synapses, noise # loads MOD files

# Custom cell models
import bgcellmodels.models.STN.GilliesWillshaw.gillies_pynn_model as gillies
import bgcellmodels.models.GPe.Gunay2008.gunay_pynn_model as gunay

import bgcellmodels.cellpopdata.connectivity as connectivity # for use in config files
ConnectivityPattern = connectivity.ConnectivityPattern
make_connection_list = connectivity.make_connection_list
make_divergent_pattern = connectivity.make_divergent_pattern

# Our physiological parameters
# from bgcellmodels.cellpopdata.physiotypes import Populations as PopID
#from bgcellmodels.cellpopdata.physiotypes import ParameterSource as ParamSrc
# from bgcellmodels.cellpopdata.cellpopdata import CellConnector

from bgcellmodels.common.configutil import eval_params
from bgcellmodels.common.stdutil import getdictvals
from bgcellmodels.common import logutils, fileutils

# Debug messages
logutils.setLogLevel('quiet', [
    'Neo',
    'bpop_ext',
    'bluepyopt.ephys.parameters',
    'bluepyopt.ephys.mechanisms',
    'bluepyopt.ephys.morphologies'])


def nprint(*args, **kwargs):
    """
    Print only on host with rank 0.
    """
    if mpi_rank == 0:
        print(*args, **kwargs)


def make_stn_lateral_connlist(pop_size, num_adjacent, fraction, rng):
    """
    Make connection list for STN lateral connectivity.

    @param  pop_size : int
            Number of cells in STN population.

    @param  fraction : float
            Fraction of STN neurons that project laterally to neighbors.

    @param  num_adjacent : int
            Number of neighbors on each side to project to.

    @param  rng : numpy.Random
            Random number generator.

    @return conn_list : list(list[int, int])
            Connection list with [source, target] pairs.
    """
    if fraction == 0 or num_adjacent == 0:
        return []
    source_ids = rng.choice(range(pop_size), int(fraction*pop_size), replace=False)
    targets_relative = range(1, num_adjacent+1) + range(-1, -num_adjacent-1, -1)
    return make_divergent_pattern(source_ids, targets_relative, pop_size)



def write_population_data(pop, output, suffix, gather=True, clear=True):
    """
    Write recorded data for Population to file.

    @param  output : str
            Output path including asterisk as placeholder: "/path/to/*.ext"

    @note   gathers data from MPI nodes so should be executed on all ranks.
    """
    if output is None:
        return
    outdir, extension = output.split('*')
    # Get Neo IO writer for file format associated with extension
    if extension.endswith('h5'):
        IOClass = neo.io.NixIO
    elif extension.endswith('mat'):
        IOClass = neo.io.NeoMatlabIO
    elif extension.endswith('npz'):
        IOClass = neo.io.PyNNNumpyIO
    else:
        IOClass = str # let PyNN guess from extension
    outfile =  "{dir}{label}{suffix}{ext}".format(dir=outdir,
                    label=pop.label, suffix=suffix, ext=extension)
    io = IOClass(outfile)
    pop.write_data(io, variables='all', gather=gather, clear=clear,
                       annotations={'script_name': __file__})


def run_simple_net(
        pop_scale       = 1.0,
        sim_dur         = 500.0,
        export_locals   = True,
        with_gui        = True,
        output          = None,
        report_progress = None,
        config          = None,
        seed            = None,
        calculate_lfp   = None,
        dopamine_depleted = None,
        transient_period = None,
        max_write_interval = None):
    """
    Run a simple network consisting of an STN and GPe cell population
    that are reciprocally connected.

    @param      output : str (optional)
                File path to save recordings at in following format:
                '~/storage/*.mat'


    @param      config : dict
                Dictionary with one entry per population label and one
                key 'simulation' for simulation parameters.
    """

    ############################################################################
    # SIMULATOR SETUP
    ############################################################################

    sim.setup(timestep=0.025, min_delay=0.1, max_delay=10.0, use_cvode=False)
    if mpi_rank == 0:
        init_logging(logfile=None, debug=True)


    print("""\nRunning net on MPI rank {} with following settings:
    - sim_dur = {}
    - output = {}""".format(mpi_rank, sim_dur, output))

    print("\nThis is node {} ({} of {})\n".format(
          sim.rank(), sim.rank() + 1, sim.num_processes()))

    h = sim.h
    sim.state.duration = sim_dur # not used by PyNN, only by our custom funcs
    sim.state.rec_dt = 0.05
    sim.state.mcellran4_rng_indices = {} # Keep track of MCellRan4 indices for independent random streams.
    finit_handlers = []

    # Make one random generator that is shared and should yield same results
    # for each MPI rank, and one with unique results.
    # - The shared (parallel-safe) RNGs should be used in functions that are
    #   executed on all ranks, like instantiating Population and Projection
    #   objects.
    # - The default RNG for Connectors is NumpyRNG(seed=151985012)
    if seed is None:
        seed = config['simulation']['shared_rng_seed']
    sim.state.shared_rng_seed = shared_seed = seed # original: 151985012
    sim.state.rank_rng_seed = rank_seed = sim.state.native_rng_baseseed + sim.state.mpi_rank
    # RNGs that can be passed to PyNN objects like Connector subclasses
    # Store them on simulator.state so we can access from other custom classes
    sim.state.shared_rng = shared_rng_pynn = sim.NumpyRNG(seed=shared_seed)
    sim.state.rank_rng = rank_rng_pynn = sim.NumpyRNG(seed=rank_seed)
    # Raw Numpy RNGs (numpy.random.RandomState) to be used in our own functions
    shared_rng = shared_rng_pynn.rng
    rank_rng = rank_rng_pynn.rng

    # Global physiological conditions
    DD = dopamine_depleted
    if DD is None:
        DD = dopamine_depleted = config['simulation'].get('DD', None)
    if DD is None:
        raise ValueError("Dopamine depleted condition not specified "
                         "in config file nor as simulation argument.")
    nprint("Dopamine state is " + "DEPLETED" if DD else "NORMAL")

    ############################################################################
    # LOCAL FUNCTIONS
    ############################################################################

    params_global_context = globals()
    params_local_context = locals() # capture function arguments

    def get_pop_parameters(pop, *param_names):
        """
        Get population parameters from config and evaluate them.
        """
        config_locals = config[pop].get('local_context', {})
        param_specs = getdictvals(config[pop], *param_names, as_dict=True)
        pvals = eval_params(param_specs, params_global_context,
                            [params_local_context, config_locals])
        return getdictvals(pvals, *param_names)

    def get_param_group(pop, group_name=None, mapping=None):
        """
        Get a group of parameters for a population as dictionary.
        """
        config_locals = config[pop].get('local_context', {})
        if group_name is None:
            param_specs = config[pop]
        else:
            param_specs = config[pop][group_name]
        if mapping is not None:
            param_specs = {v: param_specs[k] for k,v in mapping.items()}
        return eval_params(param_specs, params_global_context,
                           [params_local_context, config_locals])

    def get_cell_parameters(pop):
        """
        Get PyNN cell parameters as dictionary of numerical values.
        """
        config_locals = config[pop].get('local_context', {})
        param_specs = config[pop].get('PyNN_cell_parameters', {})
        return eval_params(param_specs, params_global_context,
                           [params_local_context, config_locals])

    def synapse_from_config(pre, post):
        """
        Make Synapse object from config dict
        """
        config_locals = config[post].get('local_context', {})
        syn_type, syn_params = getdictvals(config[post][pre]['synapse'],
                                           'name', 'parameters')
        syn_class = synapse_types[syn_type]
        syn_pvals = eval_params(syn_params, params_global_context,
                                [params_local_context, config_locals])
        num_contacts = config[post][pre].get('num_contacts', 1)
        syntype_obj = syn_class(**syn_pvals)
        syntype_obj.num_contacts = num_contacts
        return syntype_obj

    def connector_from_config(pre, post, rng=None):
        """
        Make Connector object from config dict
        """
        config_locals = config[post].get('local_context', {})
        con_type, con_params = getdictvals(config[post][pre]['connector'],
                                           'name', 'parameters')
        connector_class = getattr(sim, con_type)
        con_pvals = eval_params(con_params, params_global_context,
                               [params_local_context, config_locals])
        connector = connector_class(**con_pvals)
        if rng is not None:
            connector.rng = rng
        return connector


    # LFP calculation: command line args get priority over config file
    if calculate_lfp is None:
        calculate_lfp, = get_pop_parameters('STN', 'calculate_lfp')

    # Set NEURON integrator/solver options
    if calculate_lfp:
        sim.state.cvode.use_fast_imem(True)
    sim.state.cvode.cache_efficient(True) # necessary for lfp, also 33% reduction in simulation time

    ############################################################################
    # POPULATIONS
    ############################################################################
    # Define each cell population with its cell type, number of cells
    # NOTE:
    # - to query cell model attributes, use population[i]._cell
    print("rank {}: starting phase POPULATIONS.".format(mpi_rank))

    config_pop_labels = [k for k in config.keys() if not k in ('simulation',)]

    #===========================================================================
    # CTX POPULATION

    # CTX spike sources
    ctx_pop_size, = get_pop_parameters('CTX', 'base_population_size')
    ctx_burst_params = get_param_group('CTX', 'spiking_pattern')
    spikegen_name = ctx_burst_params.pop('algorithm')
    spikegen_func = getattr(spikegen, spikegen_name)

    ctx_spike_generator = spikegen_func(duration=sim_dur,
                                        rng=rank_rng,
                                        **ctx_burst_params)

    pop_ctx = Population(
        int(ctx_pop_size * pop_scale),
        cellclass=sim.SpikeSourceArray(spike_times=ctx_spike_generator),
        label='CTX')

    #===========================================================================
    # STR.MSN POPULATION

    # STR.MSN spike sources
    msn_pop_size, = get_pop_parameters(
        'STR.MSN', 'base_population_size')
    # msn_pop_size, synchronous, bursting_fraction, = get_pop_parameters(
    #     'STR.MSN', 'base_population_size', 'synchronous', 'bursting_fraction')
    # T_burst, dur_burst, f_intra, f_inter, f_background = get_pop_parameters(
    #     'STR.MSN', 'T_burst', 'dur_burst', 'f_intra', 'f_inter', 'f_background')
    # burst_intervals, randomize_bursting = get_pop_parameters(
    #     'STR.MSN', 'burst_intervals', 'randomize_bursting_cells')

    # if burst_intervals is not None:
    #     msn_spike_generator = spikegen.bursty_spiketrains_during(
    #                             burst_intervals, bursting_fraction,
    #                             T_burst, dur_burst, f_intra, f_inter, f_background,
    #                             sim_dur, randomize_bursting, rank_rng)
    # else:
    #     msn_spike_generator = spikegen.make_bursty_spike_generator(
    #                             bursting_fraction=bursting_fraction,
    #                             synchronous=synchronous, rng=rank_rng,
    #                             T_burst=T_burst, dur_burst=dur_burst,
    #                             f_intra=f_intra, f_inter=f_inter,
    #                             f_background=f_background, duration=sim_dur)

    msn_burst_params = get_param_group('STR.MSN', 'spiking_pattern')
    spikegen_name = msn_burst_params.pop('algorithm')
    spikegen_func = getattr(spikegen, spikegen_name)
    msn_spike_generator = spikegen_func(duration=sim_dur,
                                        rng=rank_rng,
                                        **msn_burst_params)

    pop_msn = Population(
        int(msn_pop_size * pop_scale),
        cellclass=sim.SpikeSourceArray(spike_times=msn_spike_generator),
        label='STR.MSN')

    #===========================================================================
    # STN POPULATION
    stn_ncell_base, stn_dx, = get_pop_parameters('STN',
                                'base_population_size', 'grid_dx')
    stn_grid = space.Line(x0=0.0, dx=stn_dx, y=0.0, z=0.0)
    stn_ncell_biophys = int(stn_ncell_base * pop_scale)

    # FIXME: set electrode coordinates
    stn_cell_params = get_cell_parameters('STN')
    stn_type = gillies.StnCellType(
                        calculate_lfp=calculate_lfp,
                        **stn_cell_params)

    vinit = stn_type.default_initial_values['v']
    initial_values={
        'v': RandomDistribution('uniform', (vinit-5, vinit+5), rng=shared_rng_pynn)
    }

    pop_stn = Population(stn_ncell_biophys,
                         cellclass=stn_type,
                         label='STN',
                         structure=stn_grid,
                         initial_values=initial_values)

    #---------------------------------------------------------------------------
    # STN Surrogate spike sources

    frac_surrogate, surr_rate = get_pop_parameters('STN',
        'surrogate_fraction', 'surrogate_rate')

    ncell_surrogate = int(stn_ncell_biophys * frac_surrogate)
    if ncell_surrogate > 0:
        pop_stn_surrogate = Population(ncell_surrogate,
                                       sim.SpikeSourcePoisson(rate=surr_rate),
                                       label='STN.surrogate')
    else:
        pop_stn_surrogate = None

    #---------------------------------------------------------------------------
    # STN Assembly (Biophys + Surrogate)

    if pop_stn_surrogate is None:
        asm_stn = sim.Assembly(pop_stn, pop_stn_surrogate,
                               label='STN.all')
    else:
        asm_stn = sim.Assembly(pop_stn, label='STN.all')
    stn_pop_size = asm_stn.size

    #===========================================================================
    # GPE POPULATION (prototypic)

    # Get common parameters for GPE cells
    gpe_dx, gpe_ncell_base, frac_proto, frac_arky = get_pop_parameters(
        'GPE.all', 'grid_dx', 'base_population_size',
        'prototypic_fraction', 'arkypallidal_fraction')

    gpe_common_params = get_cell_parameters('GPE.all')

    gpe_grid = space.Line(x0=0.0, dx=gpe_dx,
                          y=1e6, z=0.0)

    proto_type = gunay.GpeProtoCellType(**gpe_common_params)

    vinit = proto_type.default_initial_values['v']
    initial_values={
        'v': RandomDistribution('uniform', (vinit-5, vinit+5), rng=shared_rng_pynn)
    }

    ncell_proto = int(gpe_ncell_base * pop_scale * frac_proto)
    pop_gpe_proto = Population(ncell_proto,
                               cellclass=proto_type,
                               label='GPE.proto',
                               structure=gpe_grid,
                               initial_values=initial_values)


    #---------------------------------------------------------------------------
    # GPE Surrogate spike sources

    frac_surrogate, surr_rate = get_pop_parameters('GPE.all',
        'surrogate_fraction', 'surrogate_rate')

    ncell_surrogate = int(gpe_ncell_base * pop_scale * frac_surrogate)
    if ncell_surrogate > 0:
        pop_gpe_surrogate = Population(ncell_surrogate,
                                       sim.SpikeSourcePoisson(rate=surr_rate),
                                       label='GPE.surrogate')
    else:
        pop_gpe_surrogate = None

    #---------------------------------------------------------------------------
    # GPE Assembly (Proto + Arky + Surrogate)

    if pop_gpe_surrogate is None:
        asm_gpe = sim.Assembly(pop_gpe_proto, pop_gpe_surrogate,
                               label='GPE.all')
    else:
        asm_gpe = sim.Assembly(pop_gpe_proto, label='GPE.all')
    gpe_pop_size = asm_gpe.size


    #===========================================================================

    # Make all Population and Projection objects accessible by label
    all_pops = {pop.label : pop for pop in Population.all_populations}
    all_asm = {asm.label: asm for asm in (asm_gpe,)}
    all_proj = {pop.label : {} for pop in Population.all_populations}
    all_proj[asm_gpe.label] = {} # add Assembly projections manually

    # Make distinction between 'real' and surrogate subpopulations
    # (note: NativeCellType is common base class for all NEURON cells)
    biophysical_pops = [pop for pop in Population.all_populations if isinstance(
                        pop.celltype, sim.cells.NativeCellType)]
    artificial_pops = [pop for pop in Population.all_populations if not isinstance(
                        pop.celltype, sim.cells.NativeCellType)]

    # Update local context for eval() statements
    params_local_context.update(locals())


    ############################################################################
    # CONNECTIONS
    ############################################################################

    # Allowed synapse types (for creation from config file)
    synapse_types = {
        "GluSynapse": GluSynapse,
        "GABAAsynTM": GABAAsynTM,
        "GabaSynTm2": GabaSynTm2,
        "GabaSynTmHill" : GabaSynTmHill, # Desthexhe-like signaling pathway
        "NativeMultiSynapse" : NativeMultiSynapse,
    }

    # Make all Projections directly from (pre, post) pairs in config
    for post_label, pop_config in config.items():
        # get post-synaptic Population
        if post_label in all_pops.keys():
            post_pop = all_pops[post_label]
        elif post_label in all_asm.keys():
            post_pop = all_asm[post_label]
        else:
            continue
        print("rank {}: starting phase {} AFFERENTS.".format(mpi_rank, post_label))

        for pre_label in pop_config.keys():
            # get pre-synaptic Population
            if pre_label in all_pops.keys():
                pre_pop = all_pops[pre_label]
            elif pre_label in all_asm.keys():
                pre_pop = all_asm[pre_label]
            else:
                continue
            proj_config = pop_config[pre_label]

            # make PyNN Projection
            all_proj[pre_label][post_label] = sim.Projection(pre_pop, post_pop,
                connector=connector_from_config(pre_label, post_label, rng=shared_rng_pynn),
                synapse_type=synapse_from_config(pre_label, post_label),
                receptor_type=proj_config['receptor_type'])

    #---------------------------------------------------------------------------
    # Post-constructional modifications

    # Reduce dendritic branching and number of GLU synapses in DD
    num_prune = config['STN'].get('prune_dendritic_GLUR', 0)
    if DD and num_prune > 0:
        # PD: dendritic AMPA & NMDA-NR2B/D afferents pruned
        num_disabled = np.zeros(pop_stn.size)
        for conn in all_proj['CTX']['STN'].connections:
            if num_disabled[conn.postsynaptic_index] < num_prune:
                conn.GLUsyn_gmax_AMPA = 0.0
                conn.GLUsyn_gmax_NMDA = 0.0
                num_disabled[conn.postsynaptic_index] += 1

    # Disable somatic/proximal fast NMDA subunits
    if config['STN'].get('disable_somatic_NR2A', False):
        # NOTE: config uses a separate NMDAsyn point process for somatic NMDAR
        all_proj['CTX']['STN'].set(NMDAsynTM_gmax_NMDA=0.0)

    # Only allow GABA-B currents on reported fraction of cells
    num_without_GABAB = config['STN'].get('num_cell_without_GABAB', 0)
    if num_without_GABAB > 0:
        # Pick subset of cells with GABA-B disabled
        pop_sample = pop_stn.sample(num_without_GABAB, rng=shared_rng_pynn)
        stn_ids = pop_sample.all_cells  # global ids
        for pre in 'GPE.all', 'GPE.proto', 'GPE.surrogate':
            if pre in all_proj and 'STN' in all_proj[pre]:
                for conn in all_proj[pre]['STN'].connections:
                    if conn.postsynaptic_cell in stn_ids:
                        conn.gmax_GABAB = 0.0
                        print('Disabled GABAB on STN cell with id {}'.format(conn.postsynaptic_cell))
        # TODO: for cells without GABAB, create new Projection with only GABA-A synapses
        #       - either from surrogate only or whole population (choose)

    #---------------------------------------------------------------------------
    # Sanity check: make sure all populations and projections are instantiated

    undefined_pops = [cpop for cpop in config_pop_labels if (
                        cpop not in all_pops and cpop not in all_asm)]
    undefined_proj = [(pre, post) for (post, pre) in config.items() if (
                        (pre in config_pop_labels and post in config_pop_labels)
                        and (pre not in all_proj or post not in all_proj[pre]))]

    err_msg = ''
    if len(undefined_pops) > 0:
        err_msg += ("\nFollowing populations in config file were not "
                    "instantiated in simulator: {}".format(undefined_pops))

    if len(undefined_proj) > 0:
        err_msg += ("\nFollowing projections in config file were not "
                    "instantiated in simulator: {}".format(undefined_proj))

    if err_msg:
        raise Exception(err_msg)

    ############################################################################
    # RECORDING
    ############################################################################
    print("rank {}: starting phase RECORDING.".format(mpi_rank))

    # Default traces
    traces_biophys = {
        'Vm':       {'sec':'soma[0]', 'loc':0.5, 'var':'v'},
        # Can use SynMech[a:b:c], key must be formattable with index.
        # 'gAMPA{:d}': {'syn':'GLUsyn[0]', 'var':'g_AMPA'},
        # 'gNMDA{:d}': {'syn':'GLUsyn[::2]', 'var':'g_NMDA'},
        # 'gGABAA{:d}': {'syn':'GABAsyn[1]', 'var':'g_GABAA'},
        # 'gGABAB{:d}': {'syn':'GABAsyn[1]', 'var':'g_GABAB'},
    }

    for pop in biophysical_pops:
        pop.record(traces_biophys.items(), sampling_interval=.05)

    for pop in list(all_pops.values()):
        pop.record(['spikes'], sampling_interval=.05)

    if calculate_lfp:
        pop_stn.record(['lfp'], sampling_interval=.05)

    # Traces defined in config file
    for pop_label, pop_config in config.items():
        if 'traces' not in pop_config:
            continue
        if pop_label in all_pops:
            target_pop = all_pops[pop_label]
        elif pop_label in all_asm:
            target_pop = all_asm[pop_label]
        else:
            raise ValueError("Unknown population to record from: {}".format(pop_label))

        # Translate trace group specifier to Population.record() call
        for trace_group in pop_config['traces']:
            pop_sample = trace_group['cells']
            if isinstance(pop_sample, int):
                target_cells = target_pop.sample(pop_sample, rng=shared_rng_pynn)
            elif isinstance(pop_sample, (str, unicode)):
                slice_args = [int(i) if i!='' else None for i in pop_sample.split(':')]
                target_cells = target_pop[slice(*slice_args)]
            elif isinstance(pop_sample, list):
                target_cells = target_pop[pop_sample]
            else:
                raise ValueError("Cannot interpret cell indices '{}'".format(pop_sample))
            target_cells.record(trace_group['specs'].items(),
                                sampling_interval=trace_group['sampling_period'])


    ############################################################################
    # INITIALIZE & SIMULATE
    ############################################################################
    print("rank {}: starting phase SIMULATE.".format(mpi_rank))

    # Set physiological conditions
    h.celsius = 36.0
    h.set_aCSF(4) # Hoc function defined in Gillies code

    # Simulation statistics
    num_segments = sum((sec.nseg for sec in h.allsec()))
    num_cell = sum((1 for sec in h.allsec()))
    each_num_segments = comm.gather(num_segments, root=0)
    if mpi_rank == 0:
        # only rank 0 receives broadcast result
        total_num_segments = sum(each_num_segments)
        print("Entire network consists of {} segments (compartments)".format(
              total_num_segments))

    print("MPI rank {} will simulate {} segments ({} sections) for {} ms.".format(
            mpi_rank, num_segments, num_cell, sim_dur))

    tstart = time.time()
    outdir, filespec = os.path.split(output)
    progress_file = os.path.join(outdir, '{}_sim_progress.log'.format(
        datetime.fromtimestamp(tstart).strftime('%Y.%m.%d-%H.%M.%S')))
    report_interval = 50.0 # (ms) in simulator time
    last_report_time = tstart

    # Times for writing out data to file
    if transient_period is None:
        transient_period = 0.0 # (ms)
    steady_period = sim_dur - transient_period
    if max_write_interval is None:
        max_write_interval = 10e3 # (ms)
    homogenize_intervals = False
    if homogenize_intervals:
        write_interval = steady_period / (steady_period // max_write_interval + 1)
    else:
        write_interval = max_write_interval
    if transient_period == 0:
        first_write_time = write_interval
    else:
        first_write_time = transient_period
    write_times = list(np.arange(first_write_time, sim_dur, write_interval)) + [sim_dur]
    last_write_time = 0.0

    # SIMULATE
    while sim.state.t < sim_dur:
        sim.run(report_interval)

        # Report simulation progress
        if mpi_rank == 0:
            tnow = time.time()
            t_elapsed = tnow - tstart
            t_stepdur = tnow - last_report_time
            last_report_time = tnow
            # ! No newlines in progress report - passed to shell
            progress = ("Simulation time is {} of {} ms. "
                        "CPU time elapsed is {} s, last step took {} s".format(
                        sim.state.t, sim_dur, t_elapsed, t_stepdur))
            print(progress)

            if report_progress:
                stamp = datetime.fromtimestamp(tnow).strftime('%Y-%m-%d@%H:%M:%S')
                os.system("echo [{}]: {} >> {}".format(stamp, progress, progress_file))

        # Write recorded data
        if len(write_times) > 0 and abs(sim.state.t - write_times[0]) <= 5.0:
            suffix = "_{:.0f}ms-{:.0f}ms".format(last_write_time, sim.state.t)
            for pop in list(all_pops.values()):
                write_population_data(pop, output, suffix, gather=True, clear=True)
            write_times.pop(0)
            last_write_time = sim.state.t

    # Report simulation statistics
    tstop = time.time()
    cputime = tstop - tstart
    each_num_segments = comm.gather(num_segments, root=0)
    if mpi_rank == 0:
        # only rank 0 receives broadcast result
        total_num_segments = sum(each_num_segments)
        print("Simulated {} segments for {} ms in {} ms CPU time".format(
                total_num_segments, sim.state.tstop, cputime))


    ############################################################################
    # WRITE PARAMETERS
    ############################################################################
    print("rank {}: starting phase INTEGRITY CHECK.".format(mpi_rank))

    # NOTE: - any call to Population.get() Projection.get() does a ParallelContext.gather()
    #       - cannot perform any gather() operations before initializing MPI transfer
    #       - must do gather() operations on all nodes
    saved_params = {'dopamine_depleted': DD}

    # Save cell information
    for pop in list(all_pops.values()) + list(all_asm.values()):
        saved_params.setdefault(pop.label, {})['gids'] = pop.all_cells.astype(int)

    # Save connection information
    for pre_pop, post_pops in all_proj.items():
        saved_params.setdefault(pre_pop, {})
        for post_pop, proj in post_pops.items():

            # Plot connectivity matrix ('O' is connection, ' ' is no connection)
            utf_matrix, float_matrix = connection_plot(proj)
            nprint("{}->{} connectivity matrix (dim[0,1] = [src,target]: \n".format(
                proj.pre.label, proj.post.label) + utf_matrix)

            # This does an mpi gather() on all the parameters
            conn_params = ["delay", "weight"]
            gsyn_params = ['gmax_AMPA', 'gmax_NMDA', 'gmax_GABAA', 'gmax_GABAB']
            conn_params.extend(
                [p for p in gsyn_params if p in proj.synapse_type.default_parameters])
            pre_post_params = np.array(proj.get(conn_params, format="list",
                                       gather='all', multiple_synapses='sum'))

            # Sanity check: minimum and maximum delays and weights
            mind = min(pre_post_params[:,2])
            maxd = max(pre_post_params[:,2])
            minw = min(pre_post_params[:,3])
            maxw = max(pre_post_params[:,3])
            nprint("Error check for projection {pre}->{post}:\n"
                  "    - delay  [min, max] = [{mind}, {maxd}]\n"
                  "    - weight [min, max] = [{minw}, {maxw}]\n".format(
                    pre=pre_pop, post=post_pop, mind=mind, maxd=maxd,
                    minw=minw, maxw=maxw))

            # Make (gid, gid) connectivity pairs
            pop_idx_pairs = [tuple(pair) for pair in pre_post_params[:, 0:2].astype(int)]
            cell_gid_pairs = [(int(proj.pre[a]), int(proj.post[b])) for a, b in pop_idx_pairs]

            # Append to saved dictionary
            proj_params = saved_params[pre_pop].setdefault(post_pop, {})
            proj_params['conn_matrix'] = float_matrix
            proj_params['conpair_pop_indices'] = pop_idx_pairs
            proj_params['conpair_gids'] = cell_gid_pairs
            proj_params['conpair_pvals'] = pre_post_params
            proj_params['conpair_pnames'] = conn_params


    # Write model parameters
    print("rank {}: starting phase WRITE PARAMETERS.".format(mpi_rank))

    if mpi_rank==0 and output is not None:
        outdir, extension = output.split('*')

        # Save projection parameters
        import pickle
        extension = extension[:-4] + '.pkl'
        params_outfile = "{dir}pop-parameters{ext}".format(dir=outdir, ext=extension)
        with open(params_outfile, 'wb') as fout:
            pickle.dump(saved_params, fout)

    ############################################################################
    # PLOT DATA
    ############################################################################

    if mpi_rank==0 and with_gui:
        # Only plot on one process, and if GUI available
        import analysis
        pop_neo_data = {
            pop.label: pop.get_data().segments[0] for pop in list(all_pops.values())
        }
        analysis.plot_population_signals(pop_neo_data)

    if export_locals:
        globals().update(locals())

    print("rank {}: SIMULATION FINISHED.".format(mpi_rank))


if __name__ == '__main__':
    # Parse arguments passed to `python model.py [args]`
    import argparse

    parser = argparse.ArgumentParser(description='Run basal ganglia network simulation')

    parser.add_argument('-d', '--dur', nargs='?', type=float, default=500.0,
                        dest='sim_dur', help='Simulation duration')

    parser.add_argument('--scale', nargs='?', type=float, default=1.0,
                        dest='pop_scale', help='Scale for population sizes')

    parser.add_argument('--seed', nargs='?', type=int, default=None,
                        dest='seed', help='Seed for random number generator')

    parser.add_argument('-wi', '--writeinterval', nargs='?', type=float, default=None,
                        dest='max_write_interval',
                        help='Interval between successive write out of recording data')

    parser.add_argument('-tp', '--transientperiod', nargs='?', type=float, default=None,
                        dest='transient_period',
                        help=('Duration of transient period at start of simulation. '
                              'First data write-out is after transient period'))

    parser.add_argument('--lfp',
                        dest='calculate_lfp', action='store_true',
                        help='Calculate Local Field Potential.')
    parser.add_argument('--nolfp',
                        dest='calculate_lfp', action='store_false',
                        help='Calculate Local Field Potential.')
    parser.set_defaults(calculate_lfp=None)

    parser.add_argument('--dd',
                        dest='dopamine_depleted', action='store_true',
                        help='Set dopamine depleted condition.')
    parser.add_argument('--dnorm',
                        dest='dopamine_depleted', action='store_false',
                        help='Set dopamine normal condition.')
    parser.set_defaults(dopamine_depleted=None)

    parser.add_argument('-g', '--gui',
                        dest='with_gui',
                        action='store_true',
                        help='Enable graphical output')
    parser.add_argument('-ng', '--nogui',
                        dest='with_gui',
                        action='store_false',
                        help='Enable graphical output')
    parser.set_defaults(with_gui=False)

    parser.add_argument('-o', '--outdir', nargs='?', type=str,
                        default='~/storage/',
                        dest='output',
                        help='Output destination in format \'/outdir/*.ext\''
                             ' or /path/to/outdir/ with trailing slash')

    parser.add_argument('-p', '--progress',
                        dest='report_progress', action='store_true',
                        help='Report progress periodically to progress file')
    parser.set_defaults(report_progress=False)

    parser.add_argument('-c', '--config', nargs=1, type=str,
                        metavar='/path/to/config.json',
                        dest='config_file',
                        help='Configuration file in JSON format')

    parser.add_argument('-id', '--identifier', nargs=1, type=str,
                        metavar='<job identifer>',
                        dest='job_id',
                        help='Job identifier to tag the simulation')

    args = parser.parse_args() # Namespace object
    parsed_dict = vars(args) # Namespace to dict

    # Parse config JSON file to dict
    config_file = os.path.expanduser(parsed_dict.pop('config_file')[0])
    config_name, ext = os.path.splitext(os.path.basename(config_file))
    sim_config = fileutils.parse_json_file(config_file, nonstrict=True)
    parsed_dict['config'] = sim_config

    # Post process output specifier
    out_basedir = parsed_dict['output']
    if out_basedir is None or out_basedir == '': # shell can pass empty string
        out_basedir = '~/storage'
    job_id = parsed_dict.pop('job_id')[0]
    time_now = time.time()
    timestamp = datetime.fromtimestamp(time_now).strftime('%Y.%m.%d_%H.%M.%S')

    # Default output directory
    # NOTE: don't use timestamp -> mpi ranks will make different filenames
    out_subdir = 'LuNetStnGpe_{stamp}_job-{job_id}_{config_name}'.format(
                                            stamp=timestamp,
                                            job_id=job_id,
                                            config_name=config_name)

    # File names for data files
    # Default output format is hdf5 / NIX io
    filespec = '*_{stamp}_scale-{scale}_dur-{dur}_job-{job_id}.mat'.format(
                                            scale=parsed_dict['pop_scale'],
                                            dur=parsed_dict['sim_dur'],
                                            stamp=timestamp,
                                            job_id=job_id)

    # Make output directory if non-existing, but only on one host
    out_basedir = os.path.expanduser(out_basedir)
    if not os.path.isdir(out_basedir) and mpi_rank == 0:
        os.mkdir(out_basedir)

    # Don't make directory with variable timestamp -> mpi ranks will make different
    out_fulldir = os.path.join(out_basedir, out_subdir)
    if not os.path.isdir(out_fulldir) and mpi_rank == 0:
        os.mkdir(out_fulldir)
    parsed_dict['output'] = os.path.join(out_fulldir, filespec)

    # Copy config file to output directory
    if mpi_rank == 0:
        import shutil
        shutil.copy2(config_file, os.path.join(out_fulldir, 'sim_config.json'))

    # Run the simulation
    run_simple_net(**parsed_dict)

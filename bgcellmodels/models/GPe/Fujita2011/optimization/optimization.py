"""
Optimization of Gunay et al. (2008) GPe neuron model.


@author     Lucas Koelman
@date       24/06/2019
"""

# Python standard library
import os.path, re
import pickle
from collections import OrderedDict

# Third party libraries
from neuron import h
import numpy as np
import matplotlib.pyplot as plt

# Optimization libraries
# import pyelectro, neurotune
# from pyelectro import analysis
from neurotune import optimizers, evaluators

# Our custom libraries
from bgcellmodels.common import analysis as recording
from bgcellmodels.models.GPe.Fujita2011 import fujita_pynn_model as fujita
from bgcellmodels.models.axon.foust2011 import AxonFoust2011


class GpeController(object):
    """ 
    The Controller class maps candidate parameter sets to raw data
    by setting up the model for the candidate parameter set and 
    simulating it

    The only interface it must adhere to is providing a run()
    method with the following signature:

    run(candidates: list(list(float)), parameters:list(string)) -> (t_trace, v_trace)
    """

    def __init__(
            self,
            plot_traces     = False,
            with_axon       = True,
            axon_filepath   = None,
            axon_identifier = None,
            protocol_func   = None,
        ):
        """
        Make new controller.
        """
        self.plot_traces = plot_traces

        if axon_filepath is None:
            axon_filepath = (
                '/home/luye/workspace/bgcellmodels/bgcellmodels/models/network/'
                'LuNetDBS/configs/axons/axon_coordinates_full.pkl')
        if axon_identifier is None:
            axon_identifier = 'axon.GPe-STN.nurbs.5'

        # Load axon trajectory
        if with_axon:
            with open(axon_filepath, 'rb') as axon_file:
                all_axons = pickle.load(axon_file)
            axon_coords = np.array(all_axons[axon_identifier]) * 1e-3
        else:
            axon_coords = []

        # Cell parameters
        cell_params = dict(fujita.FujitaGpePrototypic.default_parameters)
        cell_params['transform'] = np.eye(4)
        cell_params['streamline_coordinates_mm'] = axon_coords
        cell_params['axon_class'] = AxonFoust2011
        cell_params['with_extracellular'] = False
        cell_params['owning_gid'] = 1

        # instantiate cell
        self.model = fujita.FujitaGpePrototypic.model(**cell_params)
        self.soma = list(self.model.icell.somatic)[0]
        self.ais = list(self.model.icell.axonal)[0]

        # Modify cell
        self.soma.diam = 25.0

        # Save default parameters
        self.set_base_parameters()
        self.set_recordings()

        # Set up protocol once (cell model not rebuilt)
        getattr(self, protocol_func)()


    def set_base_parameters(self):
        """
        Set base parametes so we can reset them.

        @post   self.base_parameters is dict[str, float] containing all parameter
                names and their base values
        """
        self.base_parameters = {}
        self.target_sections = {'soma': self.soma}

        for mech, varnames in self.model._mechs_params_dict.items():
            for name in varnames:
                rangevar = name + '_' + mech
                for secname, sec in self.target_sections.items():
                    if h.ismembrane(mech, sec=sec):
                        opt_param_name = secname + '_' + rangevar
                        self.base_parameters[opt_param_name] = getattr(
                            sec(0.5), rangevar)


    def set_recordings(self):
        """
        Set up all recordings.
        """
        self.rec_targets = {
            'soma': self.soma,
            'ais': self.ais,
        }

        self.trace_specs = OrderedDict([
            ('V_soma', {'var':'v', 'sec':'soma', 'loc': 0.5}),
            ('V_ais',  {'var':'v', 'sec':'ais', 'loc': 0.5}),
            ('t_global', {'var':'t'})
        ])

        # Record
        rec_dt = 0.05
        vec_dict, markers = recording.recordTraces(
                                self.rec_targets, self.trace_specs, rec_dt)

        self.rec_data = vec_dict
        self.rec_makers = markers


    def run(self, candidates, parameter_names):
        """
        Generate traces for each candidate by simulating it.

        @param  candidates: list[list[float]]
                List of candidates. Each candidate is a list of parameter values.

        @param  parameters: list[str]
                Names of parameters making up a candidate.

        @return candidate_results: list[tuple[np.array, np.array]]
                List of sample times and signal values, e.g. (t, v) for each
                candidate
        """
        traces = []
        for candidate_vals in candidates:
            param_names_vals = dict(zip(parameter_names, candidate_vals))
           
            t, v = self.run_individual(param_names_vals)
            traces.append([t,v])

        return traces


    def run_individual(self, candidate):
        """
        Generate simulation results for a candidate parameter set.

        @param  candidate : dict[str, float]
                Candidate as mapping of paramer names to values.
        """
        # Adjust cell model using candidate parametres
        for param_name, param_value in candidate.items():
            if param_name in self.base_parameters.keys():
                self.apply_parameter(param_name, param_value)
            else:
                raise Exception("Unrecognized parameter '{}' with value <{}>".format(
                                param_name, param_value))

        h.run()

        v_rec = np.array(self.rec_data['V_soma'])
        t_rec = np.array(self.rec_data['t_global'])

        if self.plot_traces:
            plt.plot(t_rec, v_rec)

        return t_rec, v_rec


    def apply_parameter(self, param_name, param_value):
        """
        Apply parameter specified as '<secname>_<param_name>''
        """
        matches = re.search(r'^([a-zA-Z0-9]+)_(.*)', param_name)
        if matches is None:
            raise Exception("Unrecognized parameter '{}' with value <{}>".format(
                                param_name, param_value))
        
        sec_name, sec_param_name = matches.groups() # we want error if len(groups) != 2
        sec = self.target_sections[sec_name]
        
        if sec_param_name in ('Ra', 'L'):
            setattr(sec, sec_param_name, param_value)
        else:
            for seg in sec:
                setattr(seg, sec_param_name, param_value)


    def save_latest_trace(self, trace_name, file_path):
        """
        Save trace to file in PyElectro-compatible CSV format.
        """
        V_soma = np.array(self.rec_data[trace_name], ndmin=2)
        T_soma = np.array(self.rec_data['t_global'], ndmin=2)
        TV_soma = np.concatenate((T_soma, V_soma), axis=0) * 1e-3 # pyelectro expects SI units: seconds, Volts
        np.savetxt(file_path, TV_soma.T, delimiter=',', fmt=['%.3E', '%.7E'])
        print("Wrote trace to " + os.path.abspath(file_path))


    def setup_protocol_SPONT(self):
        """
        Edgerton 2010, Fig. 2

        Spontaneous firing for Arkypallidal cells in shown in:

            Abdi, Mallet et al (2015), Fig. 7 : f = 3 Hz
            Bogacz, Moraud, et al (2016), Fig. 3 : f = 2 Hz
        """
        h.dt = 0.025
        h.tstop = 1000.0
        
        h.celsius = 35.0
        h.v_init = -68.0
        
        h.init()
        # h.run()
        # nrnsim.run(h.tstop, h.dt)
        
        self.protocol_vars = {}


    def setup_protocol_POSPULSE(self):
        """
        Stimulation with 100 pA

        See article Gunay (2008), Fig. 1 and Fig. 2.
        """
        # Amplitude adjustment: soma surface was changed by factor 1 / 13.4^2 == pi*1^2 / pi*13.4^2
        # However: remaining compartments were changed by smaller factor, so this is not good adjustment
        surf_factor = 0.05 # TODO: adjust surface factor for final diameter
        stim = h.IClamp(self.soma(0.5))
        stim.delay = 1000
        stim.dur = 1000
        stim.amp = 0.1 * surf_factor # 100 pA = 0.1 nA

        h.dt = 0.025
        h.tstop = 3000.0
        
        h.celsius = 35.0
        h.v_init = -68.0
        
        h.init()
        
        self.protocol_vars = {'electrodes': [stim]}


    def setup_protocol_NEGPULSE(self):
        """
        Stimulation with -100 pA

        See article Gunay (2008), Fig. 1 and Fig. 2.
        """
        # Amplitude adjustment: soma surface was changed by factor 1 / 13.4^2 == pi*1^2 / pi*13.4^2
        # However: remaining compartments were changed by smaller factor, so this is not good adjustment
        surf_factor = 0.05
        stim = h.IClamp(self.soma(0.5))
        stim.delay = 1000
        stim.dur = 1000
        stim.amp = -0.1 * surf_factor # 100 pA = 0.1 nA
        
        h.dt = 0.025
        h.tstop = 3000.0
        
        h.celsius = 35.0
        h.v_init = -68.0
        
        h.init()
        
        self.protocol_vars = {'electrodes': [stim]}


# Global optimization variables
protocol_intervals = {
    'setup_protocol_SPONT': (0.0, 1000.0),
    'setup_protocol_POSPULSE': (0.0, 3000.0), # pulse 1000-2000 ms
    'setup_protocol_NEGPULSE': (0.0, 3000.0), # pulse 1000-2000 ms
}

# TODO: re-generate trace files for Fujita model
protocol_trace_files = {
    'setup_protocol_SPONT': '/home/luye/workspace/bgcellmodels/bgcellmodels/models/GPe/Gunay2008/cellvalidation/v_soma_SPONT.csv',
    'setup_protocol_POSPULSE': '/home/luye/workspace/bgcellmodels/bgcellmodels/models/GPe/Gunay2008/cellvalidation/v_soma_POSPULSE.csv',
    'setup_protocol_NEGPULSE': '/home/luye/workspace/bgcellmodels/bgcellmodels/models/GPe/Gunay2008/cellvalidation/v_soma_NEGPULSE.csv',
}

opt_param_defaults = OrderedDict([
    # TODO: plot m_tau and m_inf variables to see which ones to optimize
    ('soma_gmax_NaF',    0.050000),
    ('soma_gmax_NaP',    0.000750),
    ('soma_gmax_Kv2',    0.000100),
    ('soma_gmax_Kv3',    0.010000),
    ('soma_gmax_Kv4f',   0.002000),
    ('soma_gmax_Kv4s',   0.001000),
    ('soma_gmax_KCNQ',   0.000150),
    ('soma_gmax_SK',     0.000400),
    ('soma_gmax_CaH',    0.000300),
    ('soma_gmax_HCN',    0.000100),
])


def optimize(export_locals=False):
    """
    Main method for the optimization routine

    HOWTO adjust for optimization:
        - Simulate protocol with full model and save trace with same time step
        - Adjust protocol for reduced model (stim timing + adjust stim current to polarize to same level)
        - Adjust optimization parameters (see ADJUSTME tags):
            - adjust parameters to optimize ('chromosome')
            - adjust seed candidates
            - adjust targets and their weights based on protocol
    """
    
    protocol_name = 'setup_protocol_SPONT'
    protocol_interval = protocol_intervals[protocol_name]
    protocol_trace_file = protocol_trace_files[protocol_name]


    # Parameters that constitute a candidate ('DNA') and their bounds
    gfac = 10.0
    

    opt_param_names = opt_param_defaults.keys()

    opt_param_bounds = OrderedDict([
        (param_name, (val/gfac, val*gfac))
            for param_name, val in opt_param_defaults.items()
    ])


    # Seeds (initial candidates), in same order as parameter names
    # SETPARAM: Add candidates in final population of past optimization run
    opt_param_seeds = [
        list(opt_param_defaults.values()),
    ]

    # Parameters for calculation of voltage trace metrics
    trace_analysis_params = {
        'peak_delta':       1e-4, # the value by which a peak or trough has to exceed its neighbours to be considered outside of the noise
        'baseline':         0.0, # !!! units of volt (V): voltage at which AP width is measured
        'dvdt_threshold':   0, # used in PPTD method described by Van Geit 2007
    }

    # Weight for components of cost function / error measure (i.e. targets)
    # NOTE: see metrics in http://pyelectro.readthedocs.io/en/latest/pyelectro.html
    # NOTE: default targets are all keys of in IClampAnalysis.analysis_results dict
    error_weights = {
        # Spike timing/frequency related
        'first_spike_time': 1.0,        # time of first AP
        'max_peak_no': 2.0,             # number of AP peaks
        'min_peak_no': 1.0,             # number of AP throughs
        'spike_frequency_adaptation': 1.0,  # slope of exp fit to initial & final frequency
        'trough_phase_adaptation': 1.0, # slope of exp fit to phase of first and last through
        'mean_spike_frequency': 2.0,    # mean AP frequency
        'interspike_time_covar': 1.0,   # coefficient of variation of ISIs
        'average_maximum': 1.0,         # average value of AP peaks
        'average_minimum': 1.0,         # average value of AP throughs
        # Spike shape related
        'spike_broadening': 1.0,        # ratio first AP width to avg of following AP widths
        'spike_width_adaptation': 1.0,  # slope of exp fit to first & last AP width
        'peak_decay_exponent': 1.0,     # Decay of AP peak height
        'trough_decay_exponent': 1.0,   # Decay of AP through height
        'pptd_error': 2.0               # Phase-plane trajectory density (see Van Geit (2008))
    }

    # Create controller to run simulations
    controller = GpeController(protocol_func=protocol_name)

    # Make evaluator to map candidate parameter sets to fitness values
    evaluator = evaluators.IClampEvaluator(
                        controller = controller,
                        analysis_start_time = protocol_interval[0],
                        analysis_end_time = protocol_interval[1],
                        target_data_path = protocol_trace_file,
                        parameters = opt_param_names,
                        analysis_var = trace_analysis_params,
                        weights = error_weights,
                        targets = None, # if not None: provide self-computed metrics
                        automatic = True # if automatic: metrics computed from target trace
                    )

    # make an optimizer
    min_constraints = [opt_param_bounds[par][0] for par in opt_param_names]
    max_constraints = [opt_param_bounds[par][1] for par in opt_param_names]

    # The optimizer creates an inspyred.ec.EvolutionaryComputation() algorithm
    # and calls its evolve() method (see docs at http://pythonhosted.org/inspyred/reference.html#inspyred.ec.EvolutionaryComputation)
    my_optimizer = optimizers.CustomOptimizerA(
                        max_constraints,
                        min_constraints,
                        evaluator,
                        population_size = 15, # initial number of individuals/candidates
                        max_evaluations = 400,
                        num_selected = 3, # how many individuals should become parents
                        num_offspring = 6, # total number of offspring (default = pop size)
                        num_elites = 1,
                        seeds = opt_param_seeds #  iterable collection of candidate solutions to include in the initial population
                    )

    #run the optimizer
    my_optimizer.optimize(summary_dir='.')

    # Make all variables accessible in interpreter
    if export_locals:
        globals().update(locals())


def plot_candidate():
    controller = GpeController(
        protocol_func='setup_protocol_SPONT',
        plot_traces=True)

    # dna = [1.088976095416787, 0.0102, 0.009997311369070006, 0.07472870427860101, 0.13623172224458147, 0.21694343549992162, 0.019683786548542716, 0.215855161863178, 0.00040706088365998237, 0.009355204839748313, 0.000301543803557004, 15.836359620045412, 0.4, 4.0885947477682105, 3.007863143713057, 5.964701100254942, 5.971654624191065, 0.002689233205754577]
    dna = [1.3671088282410035, 0.0035913145065284277, 2e-06, 2e-05, 4e-05, 0.11699895070177069, 0.00301167684885779, 0.19884371844046614, 0.000882161848342077, 0.002312673450659118, 0.00845191706794941, 23.82284173591584, 0.09818379230204187]
    dna = opt_param_defaults.values()
    
    candidate_descr = dict(zip(opt_param_defaults.keys()[:-5], dna))

    controller.run_individual(candidate_descr)

    plt.show(block=False)


if __name__ == '__main__':
    # optimize(export_locals=True)
    plot_candidate()

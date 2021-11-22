
# coding: utf-8

# # Optimize parameters of reduced morphology model 
# 
# Reduction method is Marasco method, 7 folding passes.

# In[1]:

# Enable connecting with ipyton console --existing
# %connect_info

# print code version (hash of checked out version)
get_ipython().system(u'git log -1 --format="%H"')

# print date and time of script execution
import datetime
print("\nNotebook executed at at {} in following directory:".format(datetime.datetime.now()))
get_ipython().magic(u'pwd')


# Import our optimization modules
from optimize.bpop_cellmodels import StnReducedModel
from optimize.bpop_optimize_stn import *

import pickle
import pprint
pp = pprint.PrettyPrinter(indent=2)


# ## Define Cell Model & Parameters

# In[3]:

# Choose model we want to optimize
fold_method, num_fold_passes = 'marasco', 7
red_model = StnReducedModel(
                name        = 'StnFolded',
                fold_method = fold_method,
                num_passes  = num_fold_passes)

nrnsim = ephys.simulators.NrnSimulator(dt=0.025, cvode_active=False)

# Make parameters
gleak_orig = 7.84112e-05
gleak_fit = 12.43169e-5 # fit to match Zin_DC (see praxis_passive.py)
dend_gl_param = ephys.parameters.NrnSectionParameter(
                    name		= 'gleak_dend_param',
                    param_name	= gleak_name,
                    locations	= [StnParameters.dendritic_region],
                    bounds		= [gleak_orig*1e-1, gleak_orig*1e1],
                    value		= gleak_fit, # SETPARAM: use fitted gl value
                    frozen		= True)

cm_orig = 1.0
cm_fit = cm_orig * (gleak_fit / gleak_orig) # preserve membrane time constant
dend_cm_param = ephys.parameters.NrnSectionParameter(
                    name		= 'cm_dend_param',
                    param_name	= 'cm',
                    locations	= [StnParameters.dendritic_region],
                    bounds		= [cm_orig*1e-1, cm_orig*1e1],
                    value		= cm_fit, # SETPARAM: use fitted cm value
                    frozen		= True)

# FROZEN PARAMETERS are passive parameters fit previously in passive model
frozen_params = [dend_gl_param, dend_cm_param] # SETPARAM: frozen params from previous optimisations

# FREE PARAMETERS are active conductances with large impact on response
free_params = StnParameters.dend_active_params # SETPARAM: parameters that are optimised (must be not frozen)


# ## Protocols for optimisation

# In[4]:

stn_model_type = StnModel.Gillies_FoldMarasco # SETPARAM: model type to optimise

# Protocols to use for optimisation
opt_stim_protocols = [CLAMP_REBOUND, MIN_SYN_BURST]

# Make all protocol data
red_protos = {stim_proto: BpopProtocolWrapper.make(stim_proto, stn_model_type) 
                for stim_proto in opt_stim_protocols}

# Collect al frozen mechanisms and parameters required for protocols to work
proto_mechs, proto_params = BpopProtocolWrapper.all_mechs_params(red_protos.values())

# Distinguish between sets of parameters (used, frozen, free/optimised)
frozen_params += proto_params
used_params = frozen_params + free_params
for param in frozen_params: assert param.frozen
for param in free_params: assert (not param.frozen)

# Assign parameters to reduced model
red_model.set_mechs(proto_mechs)
red_model.set_params(used_params)


# ## Features for optimisation

# In[5]:

# Get protocol responses for full model
if PROTO_RESPONSES_FILE is not None:
    full_responses = load_proto_responses(PROTO_RESPONSES_FILE)
else:
    full_protos = [BpopProtocolWrapper.make(stim_proto, stn_model_type) for stim_proto in opt_stim_protocols]
    full_mechs, full_params = BpopProtocolWrapper.all_mechs_params(full_protos)
    full_model = StnFullModel(
                    name		= 'StnGillies',
                    mechs		= full_mechs,
                    params		= full_params)
    full_responses = run_proto_responses(full_model, full_protos)

# Make EFEL feature objects
stimprotos_feats = StnFeatures.make_opt_features(red_protos.values())

# Calculate target values from full model responses
StnFeatures.calc_feature_targets(stimprotos_feats, full_responses)


# ## Objectives, Fitness Calculator, Cell Evaluator

# In[6]:

# Collect characteristic features for all protocols used in evaluation
all_opt_features, all_opt_weights = StnFeatures.all_features_weights(stimprotos_feats.values())

# # Make final objective function based on selected set of features
# total_objective = ephys.objectives.WeightedSumObjective(
#                                         name	= 'optimise_all',
#                                         features= all_opt_features,
#                                         weights	= all_opt_weights)
# all_objectives = [total_objective]

# ALTERNATIVE: set weights using 'exp_std'
all_objectives = [ephys.objectives.SingletonObjective(f.name, f) for f in all_opt_features]

# Calculator maps model responses to scores
fitcalc = ephys.objectivescalculators.ObjectivesCalculator(all_objectives)

# Make evaluator to evaluate model using objective calculator
opt_ephys_protos = {k.name: v.ephys_protocol for k,v in red_protos.items()}
opt_params_names = [param.name for param in free_params]

cell_evaluator = ephys.evaluators.CellEvaluator(
                    cell_model			= red_model,
                    param_names			= opt_params_names, # fitted parameters
                    fitness_protocols	= opt_ephys_protos,
                    fitness_calculator	= fitcalc,
                    sim					= nrnsim,
                    isolate_protocols	= True) # SETPARAM: enable multiprocessing


# ## Evaluate initial model - adjust objective scales

# In[7]:

used_param_names = [p.name for p in used_params]
default_params = {k:v for k,v in StnParameters.default_params.items() if k in used_param_names}

# Evaluate initial model
scores = cell_evaluator.evaluate_with_dicts(default_params)
pp.pprint(scores)

# NOTE: efeature objects are not copied, just references, so can change these
for efeat, weight in zip(all_opt_features, all_opt_weights):
    logger.debug('Scaling EFeature {} : exp_std / weight = {} / {}'.format(efeat.name, scores[efeat.name], weight))
    efeat.exp_std = scores[efeat.name] / weight
    
# Verify: all scores schould be 1.0 * weight
# init_scores = cell_evaluator.evaluate_with_dicts(default_params)
# pp.pprint(init_scores)


################################################################################
# Optimisation
################################################################################

for deap_seed in range(3, 20):

    # Optimisation parameters
    num_generations = 10 # Maximum number of generations in optimization
    parallel = True

    # Checkpoints: for each generation save [population, generation, parents, halloffame, history, logbook, rndstate]
    import uuid
    uuid_head = str(uuid.uuid1())[0:8]
    checkpoints_file = '/home/luye/cloudstore_m/simdata/marasco_folding/opt_checkpoints_' + uuid_head + '.pkl'

    # Make optimisation using the model evaluator
    optimisation = bpop.optimisations.DEAPOptimisation(
                        evaluator		= cell_evaluator,
                        offspring_size	= 15,
                        map_function = get_map_function(parallel),
                        seed = deap_seed)


    # In[ ]:

    # Run optimisation
    final_pop, hall_of_fame, logs, hist = optimisation.run(
                                            max_ngen = num_generations,
                                            cp_filename = checkpoints_file)


    # ## Save optimisation data

    # Save dict {StimProtol: {feat_name : {exp_mean/std : value } } }
    proto_feat_info = {}
    for stim_proto, feat_dict in stimprotos_feats.items():
        proto_feat_info[stim_proto.name] = {
            feat_name: {'weight': feat_data[1], 'exp_mean': feat_data[0].exp_mean, 'exp_std': feat_data[0].exp_std} 
                for feat_name, feat_data in feat_dict.items()
        }

    # Module info
    import sys
    code_version_info = {}
    modulenames = set(sys.modules) & set(('bluepyopt', 'efel', 'elephant')) # & set(globals())
    for module_name in modulenames:
        code_version_info[module_name] = getattr(sys.modules[module_name], '__version__', 'unknown')
    head_SHA1 = get_ipython().magic(u'sx git log -1 --format="%H"')
    code_version_info['bgcellmodels'] = head_SHA1


    info_opt = {
        'opt_param_names': cell_evaluator.param_names, # in same order as individual params
        'free_params': {p.name: p.value for p in free_params},
        'frozen_params': {p.name: p.value for p in frozen_params},
        'deap_seed': deap_seed,
        'protos_feats': proto_feat_info,
        'code_version_info': code_version_info,
    }

    # pp.pprint(info_opt)
    with open(checkpoints_file, 'ab') as f: # appends to file stream, load using second 'pickle.load()'
        pickle.dump(info_opt, f)


    # ## Analyze optimization results

    # Import our analysis modules
    from optimize.bpop_analysis_stn import plot_log

    # Plot evolution of fitness values
    plot_log(logs)


    # Print fittest individual (best parameter set)
    best_params = optimisation.evaluator.param_dict(hall_of_fame[0])
    print(pp.pprint(best_params))

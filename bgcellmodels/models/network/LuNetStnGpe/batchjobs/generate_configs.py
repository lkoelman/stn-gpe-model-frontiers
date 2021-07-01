"""
Generate simulation configuration files from a template.

USAGE
-----

The script does not work with input arguments. Just modify the parameters
directly in the script, by looking for lines marked with 'SETPARAM'.

>>> python generate_configs.py

To generate a 2D sweep of two independent variables:
- generate configs for sweep of one variable
- use these as inputs for the next sweep: paste in 'template_paths'
"""
import json, os.path
from bgcellmodels.common import fileutils
import numpy as np

# SETPARAM: template file and output directory
template_paths = """
/home/luye/workspace/bgcellmodels/bgcellmodels/models/network/LuNetStnGpe/configs/syn-V18_severe-gpe-to-stn.json
""".strip().split()

for template_path in template_paths:

    template_dir, template_name = os.path.split(template_path)
    outdir = "../configs/sweeps_f-burst-input" # SETPARAM: output dir
    config = fileutils.parse_json_file(template_path, nonstrict=True, ordered=True)

    # SETPARAM: substitutions
    factors = f_beta = np.arange(3.0, 63.0, 3.0)

    # Partial derivatives: correct changes in average population spike rate
    dfavg_dfbeta = 0.1 * 0.050 * (4.0/0.050 - 10.0)
    dfavg_from_dfbeta = 0.35 * (f_beta - 20.0) # d(f_spk)/d(f_beta) * d(f_beta) : change in average rate because of reduction in beta, increase for f_beta > 20
    dfavg_dfbg = 1 - (0.1 * 0.050 * f_beta) # d(f_spk)/d(f_spk_bg) : 
    dfbg = dfavg_from_dfbeta / dfavg_dfbg
    f_bg = 10.0 - dfbg

    f_bg_min = 2.0
    f_bg[f_bg < f_bg_min] = f_bg_min
    dfavg_dpcell = f_beta * 0.050 * (4.0/0.050 - f_bg_min) # fill in f_bg at this point
    dfavg_from_dfbeta_dfbg = (dfavg_from_dfbeta - dfavg_dfbg * (10.0 - f_bg_min))
    dpcell = dfavg_from_dfbeta_dfbg / dfavg_dpcell
    dpcell[dpcell < 0] = 0.0
    pcell = 0.1 - dpcell

    
    # gs_gabaa = config['STN']['GPE.all']['synapse']['parameters'][
    #                   'gmax_GABAA']['locals']['gmax_base']
    # gs_gabab = config['STN']['GPE.all']['synapse']['parameters'][
    #                   'gmax_GABAB']['locals']['gmax_base']
    # cs_ampa = config['STN']['CTX']['synapse']['parameters'][
    #                   'GLUsyn_gmax_AMPA']['locals']['gmax_base']
    # cs_nmda_dend = config['STN']['CTX']['synapse']['parameters'][
    #                   'GLUsyn_gmax_NMDA']['locals']['gmax_base']
    # cs_nmda_soma = config['STN']['CTX']['synapse']['parameters'][
    #                   'NMDAsynTM_gmax_NMDA']['locals']['gmax_base']
    # gg_gabaa = config['GPE.proto']['GPE.all']['synapse']['parameters'][
    #                   'gmax_GABAA']['locals']['gmax_base']
    # gg_gabab = config['GPE.proto']['GPE.all']['synapse']['parameters'][
    #                   'gmax_GABAB']['locals']['gmax_base']
    # sg_ampa  = config['GPE.proto']['STN']['synapse']['parameters'][
    #                   'gmax_AMPA']['locals']['gmax_base']
    # mg_gabaa = config['GPE.proto']['STR.MSN']['synapse']['parameters'][
    #                   'gmax_GABAA']['locals']['gmax_base']
    # ratio of scale factors for MSN-GPE and GPE-GPE to preserve EXC/INH ratio
    # mg_gg_scale_ratio = 0.0265/0.024
    # gg_factors = [1.0 - i*0.1 for i in range(1,10)]
    # mg_factors = [1.0 + i*0.1*mg_gg_scale_ratio for i in range(1,10)]

    # Replace all occurrences of format keywords
    substitutions = {
        # ('STN', 'GPE.all', 'synapse', 'parameters', 'gmax_GABAA', 'locals', 
        #     'gmax_base'): [f*gs_gabaa for f in factors],
        # ('STN', 'GPE.all', 'synapse', 'parameters', 'gmax_GABAB', 'locals', 
        #     'gmax_base'): [0.2*gs_gabab for f in factors],
        # ('STN', 'CTX', 'synapse', 'parameters', 'GLUsyn_gmax_AMPA', 'locals', 
        #     'gmax_base'): [f*cs_ampa for f in factors],
        # ('STN', 'CTX', 'synapse', 'parameters', 'GLUsyn_gmax_NMDA', 'locals', 
        #     'gmax_base'): [f*cs_nmda_dend for f in factors],
        # ('STN', 'CTX', 'synapse', 'parameters', 'NMDAsynTM_gmax_NMDA', 'locals', 
        #     'gmax_base'): [f*cs_nmda_soma for f in factors],
        # ('GPE.proto', 'GPE.all', 'synapse', 'parameters', 'gmax_GABAA', 'locals', 
        #     'gmax_base'): [1.333*gg_gabaa for f in factors],
        # ('GPE.proto', 'GPE.all', 'synapse', 'parameters', 'gmax_GABAB', 'locals', 
        #     'gmax_base'): [1.333*gg_gabab for f in factors],
        # ('GPE.proto', 'STN', 'synapse', 'parameters', 'gmax_AMPA', 'locals', 
        #     'gmax_base'): [f*sg_ampa for f in factors],
        # ('GPE.proto', 'STR.MSN', 'synapse', 'parameters', 'gmax_GABAA', 'locals', 
        #     'gmax_base'): [f*mg_gabaa for f in factors],
        ('CTX', 'spiking_pattern', 'T_burst'): [1e3/f for f in factors],
        ('CTX', 'spiking_pattern', 'f_background'): f_bg,
        ('CTX', 'spiking_pattern', 'bursting_fraction'): pcell,
    }
    suffix_format = 'severe_gpe2stn_f-burst-{:.1f}-Hz' # SETPARAM: format string for json filename
    suffix_substitutions = factors
    sweep_length = len(suffix_substitutions)

    for i in range(sweep_length):

        # Update each config parameter (nested dictionary) with sweep variable
        for nested_keys, sweep_vals in substitutions.items():
            assert len(sweep_vals) == sweep_length
            parent_dict = config
            for k in nested_keys[:-1]:
                if k not in parent_dict:
                    parent_dict[k] = {}
                parent_dict = parent_dict[k]
            sweep_param = nested_keys[-1]

            # SETPARAM: uncomment line if you want to check the modified entry exists
            if sweep_param not in parent_dict:
                raise ValueError("Key {} not present at nesting level {}".format(
                    sweep_param, nested_keys))
            # SETPARAM: substitution or multiplication of target value
            parent_dict[sweep_param] = sweep_vals[i]
            # parent_dict[sweep_param] *= sweep_vals[i]
            print("Updated key {} for sweep {}".format(sweep_param, i))

        # Write config after doing all substitutions for current sweep value
        # SETPARAM: config filename substitution
        outname = template_name.replace('severe-gpe-to-stn.json',
                        suffix_format.format(suffix_substitutions[i]) + '.json')
        outfile = os.path.join(outdir, outname)
        
        with open(outfile, 'w') as f:
            json.dump(config, f, indent=4)

        # Write unicode (encoding='utf-8')
        # import io
        # with io.open(filename, 'w', encoding=encoding) as f:
        #     f.write(json.dumps(morph_dicts, ensure_ascii=False))

        print("Wrote config file {}".format(outfile))

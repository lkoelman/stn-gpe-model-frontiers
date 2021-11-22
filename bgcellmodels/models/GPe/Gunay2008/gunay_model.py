"""
Globus pallidus multi-compartmental model by Gunay C, Edgerton JR, Jaeger D (2008).
Using channel mechanisms ported to NEURON by Kitano (2011)

@author     Lucas Koelman
@date       07-01-2018

@reference  Gunay C, Edgerton JR, Jaeger D (2008) Channel density distributions 
            explain spiking variability in the globus pallidus: a combined physiology 
            and computer simulation database approach. J Neurosci 28:7476-91

            https://senselab.med.yale.edu/modeldb/ShowModel.cshtml?model=114639


@reference  Fujita T, Fukai T, Kitano K (2012) Influences of membrane properties 
            on phase response curve and synchronization stability in a model 
            globus pallidus neuron. J Comput Neurosci 32(3):539-53
            
            https:#senselab.med.yale.edu/modeldb/ShowModel.cshtml?model=143100

@reference  Based on following example templates:
            https:#github.com/BlueBrain/BluePyOpt/blob/master/examples/l5pc/l5pc_model.py
"""

import os
from math import pi as PI

import neuron
import bluepyopt.ephys as ephys

from bgcellmodels.common import units, fileutils

units.set_units_module('pint')
h = neuron.h

# Load NEURON libraries, mechanisms
script_dir = os.path.dirname(__file__)
NRN_MECH_PATH = os.path.join(script_dir, 'mechanisms')
neuron.load_mechanisms(NRN_MECH_PATH)
h.load_file("stdlib.hoc")
h.load_file("stdrun.hoc")

from bgcellmodels.common.electrotonic import calc_min_nseg_hines, calc_lambda_AC
import bgcellmodels.extensions.bluepyopt.bpop_locations as ext_locations
from bgcellmodels.extensions.bluepyopt.bpop_parameters import NrnSegmentParameter

# Debug messages
from bgcellmodels.common import logutils
logger = logutils.getBasicLogger('gunay')
logutils.setLogLevel('quiet', ['gunay'])

# Channel mechanisms (key = suffix of mod mechanism) : max conductance parameters
gbar_dict = {
    # Nonspecific channels
    'HCN':      ['gmax'],
    'HCN2':     ['gmax'],
    'pas':      ['g'],
    # Na channels
    'NaF':      ['gmax'],
    'NaP':      ['gmax'],
    # K-channels
    'Kv2':      ['gmax'],
    'Kv3':      ['gmax'],
    'Kv4f':     ['gmax'],
    'Kv4s':     ['gmax'],
    'KCNQ':     ['gmax'],
    'SK':       ['gmax'],
    # Calcium channels / buffering
    'CaHVA':      ['gmax'],
}
gleak_name = 'g_pas'

# Mechanism parameters that are changed from default values in original model code
mechs_params_dict = {
    # Nonspecific channels
    'HCN':      ['gmax', 'e'], # HCN channel
    'HCN2':     ['gmax', 'e'], # HCN channel
    'pas':      ['g', 'e'],
    # Na channels
    'NaF':      ['gmax'],
    'NaP':      ['gmax'],
    # K-channels
    'Kv2':      ['gmax'],
    'Kv3':      ['gmax'],
    'Kv4f':     ['gmax'],
    'Kv4s':     ['gmax'],
    'KCNQ':     ['gmax'],
    'SK':       ['gmax'],
    # Calcium channels / buffering
    'Calcium':  [''],
    'CaHVA':    ['gmax', 'e'], # high-voltage-activated calcium channel
    # Built-in RANGE params
    '':         ['ena', 'ek'],
}

# All mechanism parameters that are not conductances
mechs_params_nogbar = dict(mechs_params_dict)
for mech, params in mechs_params_nogbar.items():
    for gbar_param in gbar_dict.get(mech, []):
        try:
            params.remove(gbar_param)
        except ValueError:
            pass


# List of mechanisms, max conductance params, active conductances
mechs_list = [k for k in mechs_params_dict.keys() if k!=''] # all mechanisms
gbar_list = [gname+'_'+mech for mech,chans in gbar_dict.items() for gname in chans]
active_gbar_names = [gname for gname in gbar_list if gname != gleak_name]


# Different models from Edgerton lab:
MODEL_GUNAY2008_AXONLESS = "GUNAY2008_AXONLESS"
MODEL_GUNAY2008_FULL = "GUNAY2008_FULL"
MODEL_GUNAY2008_REBOUNDBURSTING = "GUNAY2008_t1768"
MODEL_HENDRICKSON2011_AXONLESS = "HENDRICKSON2011_AXONLESS"
MODEL_HENDRICKSON2011_FULL = "HENDRICKSON2011_FULL"


def define_mechanisms(filename, exclude_mechs=None):
    """
    Create list of mechanism descriptions that link NEURON mechanisms to 
    specific regions in the cell, identified by named section lists.

    @param      filename : str
                Filename of json file containing MOD mechanisms for each
                region (section list)

    @param      exclude : list(str)
                List of mechanism names to exclude.
    
    @return     mechanisms: list(ephys.mechanisms.NrnModMechanism)
                List of NEURON mechanism descriptions as Ephys objects.
    """
    if exclude_mechs is None:
        exclude_mechs = []

    full_filename = os.path.join(script_dir, filename)
    mech_definitions = fileutils.parse_json_file(full_filename, nonstrict=True)

    mechanisms = []
    for seclist_name, mod_names in mech_definitions.items():
        
        seclist_loc = ephys.locations.NrnSeclistLocation(
                                        seclist_name,
                                        seclist_name=seclist_name)
        
        for channel in mod_names:
            if channel in exclude_mechs:
                continue
            mechanisms.append(ephys.mechanisms.NrnMODMechanism(
                            name='{}.{}'.format(channel, seclist_name),
                            mod_path=None,
                            suffix=channel,
                            locations=[seclist_loc],
                            preloaded=True))

    return mechanisms


def define_locations(locations_file):
    """
    Create named locations based on Ephys Location definitions in JSON file.

    @param      locations_file : str
                Relative path to location definitions JSON file.

    @return     locations : dict<str, ephys.Location>
    """

    fullfile = os.path.join(script_dir, locations_file)
    location_specs = fileutils.parse_json_file(fullfile, nonstrict=True)

    locations = {}
    for spec in location_specs:
        loc_name = spec['loc_name']
        locations[loc_name] = ext_locations.SomaDistanceDiamLocation(
            loc_name,
            seclist_name=spec['sectionlist'],
            distance_range=parse_values(
                spec, float, tuple, 'lower_distance', 'upper_distance'),
            diameter_range=parse_values(
                spec, float, tuple, 'lower_diameter', 'upper_diameter'))

    return locations


def parse_values(src_dict, target_type, container_type, *keys):
    return container_type(target_type(src_dict[k]) for k in keys)


def define_parameters(
        genesis_params_file,
        params_mapping_file,
        named_locations,
        exclude_mechs=None,
        genesis_params_overrides=None,
    ):
    """
    Create list of parameter descriptions that link (distributions of)
    mechanism parameters to specific regions in the cell, 
    identified by named section lists.

    Arguments
    ---------

    @param      named_locations: dict<str, ephys.Location>
                Predefined locations on the cell.

    @param      genesis_params_overrides : dict
                Dictionary containing GENESIS model parameter names and values.

    Returns
    -------
    
    @return     parameters: list(ephys.parameters.NrnParameter)
                List of NEURON parameter descriptions as Ephys objects.
    """
    if exclude_mechs is None:
        exclude_mechs = []

    fullfile = os.path.join(script_dir, genesis_params_file)
    genesis_params = fileutils.parse_json_file(fullfile, nonstrict=True)
    if genesis_params_overrides is not None:
        genesis_params.update(genesis_params_overrides)
    
    fullfile = os.path.join(script_dir, params_mapping_file)
    param_specs = fileutils.parse_json_file(fullfile, nonstrict=True)

    parameters = []

    # Create dummy section so we can query all units
    h('create dummy')
    dummysec = h.dummy
    for mech_name in mechs_list:
        dummysec.insert(mech_name)

    for param_spec in param_specs:

        # Check if parameter should be excluded
        if 'mech' in param_spec and param_spec['mech'] in exclude_mechs:
            logger.debug('Skipping parameter {} because its mechanism '
                        'is in excluded mechanisms list'.format(param_spec))
            continue

        # Get param name in NEURON
        if 'param_name' in param_spec:
            param_name = param_spec['param_name']
        elif 'mech' in param_spec and 'mech_param' in param_spec:
            param_name = '{}_{}'.format(param_spec['mech_param'], param_spec['mech'])
        else:
            raise ValueError(
                'Not enough information to resolve NEURON parameter name: {}'.format(
                    param_spec))
        
        # 'value' is an expression that maps from original GENESIS parameter name
        # to NEURON parameter value
        if 'value' in param_spec:
            frozen = False # TODO: any reason this should True/False?

            # Interpret spec: can be expression or value
            spec = param_spec['value']
            if isinstance(spec, (float, int)):
                value = spec
            elif isinstance(spec, str):
                value = eval(spec.format(**genesis_params))
            else:
                raise ValueError(
                    "Unexpected value {} for parameter '{}'".format(
                        spec, param_name))
            bounds = None

        elif 'bounds' in param_spec:
            frozen = False
            bounds = param_spec['bounds']
            value = None
        
        else:
            raise Exception(
                'Parameter config has to have bounds or value: {}'.format(
                param_spec))

        # Correct units
        if 'units' in param_spec:
            if value is not None:
                quantity = units.Quantity(value, param_spec['units'])
                converted_quantity = units.to_nrn_units(quantity, h, param_name)
                value = converted_quantity.magnitude
            if bounds is not None:
                quantity = units.Quantity(bounds, param_spec['units'])
                converted_quantity = units.to_nrn_units(quantity, h, param_name)
                bounds = converted_quantity.magnitude

        # Make Ephys description of parameter
        if param_spec['type'] == 'global':
            parameters.append(
                ephys.parameters.NrnGlobalParameter(
                    name=param_name,
                    param_name=param_name,
                    frozen=frozen,
                    bounds=bounds,
                    value=value))
        
        elif param_spec['type'] in ['section', 'segment', 'range']:
            
            # Spatial distribution of parameter
            if param_spec['dist_type'] == 'uniform':
                scaler = ephys.parameterscalers.NrnSegmentLinearScaler()
            
            elif param_spec['dist_type'] == 'exp':
                scaler = ephys.parameterscalers.NrnSegmentSomaDistanceScaler(
                    distribution=param_spec['dist'])
            
            if 'location' in param_spec:
                param_loc = named_locations[param_spec['location']]
            else:
                param_loc = ephys.locations.NrnSeclistLocation(
                                param_spec['sectionlist'],
                                seclist_name=param_spec['sectionlist'])

            name = '{}.{}'.format(param_name, param_loc.name)

            # Section parameter is uniform in a Section
            if param_spec['type'] == 'section':
                parameters.append(
                    ephys.parameters.NrnSectionParameter(
                        name=name,
                        param_name=param_name,
                        value_scaler=scaler,
                        value=value,
                        frozen=frozen,
                        bounds=bounds,
                        locations=[param_loc]))

            elif param_spec['type'] == 'segment':
                parameters.append(
                    NrnSegmentParameter(
                        name=name,
                        param_name=param_name,
                        value=value,
                        frozen=frozen,
                        locations=[param_loc]))
            
            # RANGE parameters vary over x-loc of Section
            elif param_spec['type'] == 'range':
                parameters.append(
                    ephys.parameters.NrnRangeParameter(
                        name=name,
                        param_name=param_name,
                        value_scaler=scaler,
                        value=value,
                        frozen=frozen,
                        bounds=bounds,
                        locations=[param_loc]))
        else:
            raise Exception(
                'Param config type has to be global, section, segment, or range: {}'.format(
                param_spec))

        p = parameters[-1]
        logger.debug("Created parameter description:\n" + "\n\t".join(
            ["{} : {}".format(k,getattr(p,k)) for k in ['name', 'value', 'bounds']]))

    # delete dummy section
    h.delete_section(sec=dummysec)
    del dummysec

    return parameters


def define_morphology(filename, replace_axon):
    """
    Define morphology (don't instantiate yet).

    The morphology determines the named SecArray and SectionLists 
    available as cell attributes.

    @note   The morphology is instantiated when cell.instantiate() is called.
    """

    return ephys.morphologies.NrnFileMorphology(
                os.path.join(script_dir, filename),
                do_replace_axon=False)


def fix_num_compartments(cell, f_lambda, sim, report=False):
    """
    Correct the discretization (minimum number of compartments per section)
    based on Hines' rule of thumb: L/lambda ~= 0.1
    """

    f_lambda = 100.0

    # Calculate total nseg
    all_sec = list(cell.icell.all)
    tot_nseg = sum((sec.nseg for sec in all_sec))
    if report:
        print("Number of segments before adjustment: nseg = {}".format(tot_nseg))

    # Adjust minimum number of segments
    for sec in cell.icell.all:
        # Check if we need to re-discretize
        lamb = calc_lambda_AC(f_lambda, sec.diam, sec.Ra, sec.cm)
        Lseg = sec.L/sec.nseg
        min_nseg = calc_min_nseg_hines(f_lambda, sec.L, sec.diam, sec.Ra, sec.cm, round_up=False)

        if min_nseg > sec.nseg:
            if report:
                print("Discretization too coarse:\n"
                      "\tL/lambda = {} -- nseg = {}\n"
                      "\t=> new nseg = {}".format(Lseg/lamb, sec.nseg, min_nseg))
            sec.nseg = min_nseg

    # Recalculate total nseg
    new_nseg = sum((sec.nseg for sec in all_sec))
    nseg_fraction = float(new_nseg) / tot_nseg
    if report:
        print("Number of segments after adjustment: nseg = {}".format(new_nseg))
        print("Change is {} percent".format(nseg_fraction))

    # If number of segments has changed we need to reset properties
    for param in cell.params.values():
        param.instantiate(sim=sim, icell=cell.icell)


def fix_comp_dimensions(cell):
    """
    Fix dimensions for each compartment after instantiating cell
    with specific quantities so that cable equation corresponds to
    the one solved by GENESIS.

    @see    Use of variable 'dia' in files names 'GPcomps.g'
            - Gunay (2008) code : https://senselab.med.yale.edu/modeldb/ShowModel.cshtml?model=114639&file=/GunayEdgertonJaeger2008/common/GPcomps.g#tabs-2
            - Hendrickson (2010) code : https://senselab.med.yale.edu/modeldb/ShowModel.cshtml?model=127728&file=%2farticleCode%2fcommonGPFull%2fGPcomps.g#tabs-2
    """
    # Everywhere L=1 and diam=1.
    # In soma they use sphere with diam=1, so we get L from following equality:
    # pi * diam^2 = pi*diam*L  <=>  L = diam
    for sec in cell.icell.all:
        sec.diam = 1.0
        sec.L = 1.0
        # h.define_shape() # not necessary


def define_cell(model, exclude_mechs=None):
    """
    Create GPe cell model

    @param  model : str
            Model identifier.
    """

    morphology = define_morphology(
                        'morphology/bg0121b_axonless_GENESIS_import.swc',
                        replace_axon=False)

    locations = define_locations(
                        'config/locations.json')

    mechanisms = define_mechanisms(
                        'config/mechanisms.json',
                        exclude_mechs=exclude_mechs)

    if model == MODEL_HENDRICKSON2011_AXONLESS:
        parameters = define_parameters(
                            'config/params_hendrickson2011_GENESIS.json',
                            'config/map_params_hendrickson2011.json',
                            locations,
                            exclude_mechs=exclude_mechs)
    elif model == MODEL_GUNAY2008_AXONLESS:
        parameters = define_parameters(
                            'config/params_gunay2008_GENESIS.json',
                            'config/map_params_gunay2008_v2.json',
                            locations,
                            exclude_mechs=exclude_mechs)
    elif model == MODEL_GUNAY2008_REBOUNDBURSTING:
        parameters = define_parameters(
                            'config/params_gunay2008_GENESIS_model-t1768.json',
                            'config/map_params_gunay2008_v2.json',
                            locations,
                            exclude_mechs=exclude_mechs)
    else:
        raise ValueError("Unknown model '{}'".format(model))

    cell = ephys.models.CellModel(
                        'GPe',
                        morph=morphology,
                        mechs=mechanisms,
                        params=parameters)

    return cell


def create_cell(model=MODEL_GUNAY2008_AXONLESS):
    """
    Instantiate GPe cell in NEURON simulator.
    """
    cell = define_cell(model)
    nrnsim = ephys.simulators.NrnSimulator(dt=0.025, cvode_active=False)

    cell.instantiate(sim=nrnsim)
    fix_comp_dimensions(cell)
    
    return cell, nrnsim


if __name__ == '__main__':
    # Make GPe cell
    ephys_cell, nrnsim = create_cell(model=MODEL_GUNAY2008_AXONLESS)
    icell = ephys_cell.icell

    # Save cell
    # from bgcellmodels.morphology import morph_io
    # import pickle
    # cell_data = morph_io.cell_to_dict(
    #                     section=icell.soma[0],
    #                     descr='Gunay (2008) GPe cell, axonless',
    #                     icell=icell)
    # pkl_path = 'gpe-cell_gunay2008-axonless_full.pkl'
    # with open(pkl_path, 'wb') as file:
    #         pickle.dump(cell_data, file)

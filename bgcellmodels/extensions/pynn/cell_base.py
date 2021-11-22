"""
Base classes for multi-compartmental PyNN cell models.

@author     Lucas Koelman

@date       14/02/2018

"""

# Standard library
import re
import math
import logging

# Third party libraries
import numpy as np
import quantities as pq
import bluepyopt.ephys as ephys
from neuron import h

# PyNN imports
import pyNN.neuron
from pyNN.parameters import ArrayParameter
from pyNN.neuron.cells import NativeCellType

# Our own imports
from bgcellmodels.common import nrnutil, configutil
from bgcellmodels.common import treeutils
from bgcellmodels.emfield import xtra_utils
from bgcellmodels.extensions.pynn import connection as ext_conn

# Global variables
logger = logging.getLogger('ext.pynn.cell_base')
ephys_nrn_sim = None

# Define additional units not recognized by quantities module
pq.UnitQuantity('nanovolt', pq.V * 1e-9, symbol='nV')


def ephys_sim_from_pynn():
    """
    Get Ephys NrnSimulator without changing parameters of wrapped
    NEURON simulator from those set by pyNN.
    """
    global ephys_nrn_sim

    if ephys_nrn_sim is None:
        cvode_active = pyNN.neuron.state.cvode.active()
        cvode_minstep = pyNN.neuron.state.cvode.minstep()
        ephys_nrn_sim = ephys.simulators.NrnSimulator(
                            dt=pyNN.neuron.state.dt,
                            cvode_active=cvode_active,
                            cvode_minstep=cvode_minstep)
    return ephys_nrn_sim


class UnitFetcherPlaceHolder(dict):
    """
    At the moment there isn't really a robust way to get the correct
    units for all variables.

    This can be set explicity for each possible variable that is recordable
    from the cell type, but with our custom TraceSpecRecorder the variable name
    can be anything.
    """
    def __getitem__(self, key):
        """
        Return units according to trace naming convention.

        @return     units : str
                    Unit string accepted by 'Quantities' python package
        """
        if key.lower().startswith('v'):
            return 'mV'
        elif key.lower().startswith('isyn'):
            return 'nA'
        elif key.lower().startswith('gsyn'):
            return 'uS' # can also be nS but this is for our GABAsyn/GLUsyn
        elif key.lower().startswith('i'):
            return 'mA/cm^2' # membrane currents of density mechanisms
        elif key.lower().startswith('g'):
            return 'S/cm^2' # membrane conductances of density mechanisms
        elif key.lower().startswith('lfp'):
            return 'nV' # LFP calculator mechanism 'LfpSumStep'
        elif key == 'spikes':
            return 'ms'
        else:
            return 'dimensionless'


################################################################################
# INTERFACE METHODS
################################################################################

# NOTE: can assign to class or put in mixin object representing interface


def irec_resolve_section(self, spec, multiple=False):
    """
    Resolve a section specification in the following format:

        "seclist:secname[index]"

    where seclist is the SectionList name on the NEURON cell object, secname
    is a string that mathces part of the section name, and index is an index
    in the resulting collection of matches.

    Provides compabitibility with TraceSpecRecorder.

    @return     nrn.Section
                The section intended by the given section specifier.


    @param      spec : str

                Section specifier in the format 'section_container[index]',
                where 'section_container' is the name of a secarray, Section,
                or SectionList that is a public attribute of the icell,
                and the index part is optional.


    @note       The default Section arrays for Ephys.CellModel are:
                'soma', 'dend', 'apic', 'axon', 'myelin'.

                The default SectionLists are:
                'somatic', 'basal', 'apical', 'axonal', 'myelinated', 'all'

                If 'section_container' matches any of these names, it will
                be treated as such. If not, it will be treated as an indexable
                attribute of the icell instance.
    """

    matches = re.search(
        r'^(?P<seclist>\w+)(:(?P<secname>\w+))?(?P<slice>\[(?P<index>.+)\])?', spec)
    sec_list = matches.group('seclist')
    sec_name = matches.group('secname')
    sec_slice = matches.group('slice')

    # No index: spec refers to single section
    if sec_slice is None:
        assert sec_name is None
        return getattr(self.icell, sec_list)

    seclist = list(getattr(self.icell, sec_list))
    if sec_name is not None:
        seclist = [sec for sec in seclist if sec_name in sec.name()]

    try:
        secs = configutil.index_with_str(seclist, sec_slice)
    except IndexError:
        raise IndexError("Invalid index expression {} for seclist {} in cell {}".format(
            spec, seclist, self))

    if not multiple and isinstance(secs, list):
        return secs[0]
    else:
        return secs


def init_extracellular_stim_rec(
        self, all_seclist, stim_seclist, rec_seclist, tracker_seg):
    """
    Set up extracellular stimulation and recording

    @param  seclist : neuron.SectionList

            Sections that will be the target/source for extracellular
            stimulation/recording IF they have the mechanism 'extracellular'
            inserted.
    """
    # if stim_seclist is rec_seclist:
    #     all_sections = stim_seclist
    # else:
    #     all_sections = nrnutil.join_seclists(stim_seclist, rec_seclist)
    all_sections = all_seclist

    # Insert mechanism that mediates between extracellular variables and
    # recording & stimulation routines.
    for sec in all_sections:
        if h.ismembrane('extracellular', sec=sec):
            sec.insert('xtra')

    # Store coordinates of each compartment's center in 'xtra'
    h.xtra_segment_coords_from3d(all_sections)
    h.xtra_setpointers(all_sections)

    # Set transfer impedances between electrode and each compartment's center
    if (self.transfer_impedance_matrix_um is None
            or len(self.transfer_impedance_matrix_um) == 0):
        # Set transfer impedance analytically
        x_elec, y_elec, z_elec = self.electrode_coordinates_um
        h.xtra_set_impedances_pointsource(all_sections,
            self.rho_extracellular_ohm_cm,
            x_elec, y_elec, z_elec)
    else:
        # Set using look-up table / nearest-neigbhors
        Z_coords = self.transfer_impedance_matrix_um[:, :3]
        Z_values = self.transfer_impedance_matrix_um[:, -1]
        xtra_utils.set_transfer_impedances_nearest(
            all_sections, Z_coords, Z_values, max_dist=5.0, warn_dist=0.1,
            min_electrode_dist=10.0, electrode_coords=self.electrode_coordinates_um,
            Z_intersect=1e12)

    # Disable stimulation or recording if desired
    if not (all_seclist is stim_seclist and all_seclist is rec_seclist):
        # Use python list instead of Hoc SectionList list
        # -> does not build Section stack, no stack overflow for large cells
        stim_pylist = list(stim_seclist)
        rec_pylist = list(rec_seclist)
        for sec in list(all_sections):
            if sec not in stim_pylist:
                for seg in sec: seg.scale_stim_xtra = 0.0
            if sec not in rec_pylist:
                for seg in sec: seg.scale_rec_xtra = 0.0


    # Set up LFP calculation
    # NOTE: Recorder class records lfp_tracker.summator._ref_summed
    # if logger.level <= logging.WARNING:
    #     h.XTRA_VERBOSITY = 1
    self.lfp_summator = h.xtra_sum(tracker_seg)
    self.lfp_tracker = h.ImembTracker(self.lfp_summator, rec_seclist, "xtra")


################################################################################
# BASE CLASSES
################################################################################

class MorphCellType(NativeCellType):
    """
    PyNN native cell type that has Ephys model as 'model' attribute.

    Attributes
    ----------

    @attr   units : UnitFetcherPlaceHolder
            Required by PyNN, celltype must have method `units(varname)` that
            returns the units of recorded variables


    @attr   receptor_types : list(str)
            Required by PyNN: receptor types accepted by Projection constructor.
            This attribute is created dynamically by combining
            celltype.model.region with the celltype.receptor_types declared
            in the subclass
    """

    # Population.find_units() queries this for units
    units = UnitFetcherPlaceHolder()

    # Parameter sets to be used by subclasses
    _axon_parameters = {
        'streamline_coordinates_mm': ArrayParameter([]),
        'termination_method': np.array('terminal_sequence'),
        'netcon_source_spec': np.array('main_branch:-1'),
        'collateral_branch_points_um': ArrayParameter([]),
        'collateral_target_points_um': ArrayParameter([]),
        'collateral_lvl_lengths_um': ArrayParameter([]),
        'collateral_lvl_num_branches': ArrayParameter([]),
    }

    _emf_parameters = {
        'with_extracellular_stim': False,
        'with_extracellular_rec': False,
        'seclists_with_dbs': np.array('all'), # separated by ';', e.g.:
        'seclists_with_lfp': np.array('all'), # somatic;basal;axonal
        'electrode_coordinates_um' : ArrayParameter([]),
        'rho_extracellular_ohm_cm' : 0.03,
        'transfer_impedance_matrix_um': ArrayParameter([]),

    }


    def __init__(self, **kwargs):
        """
        The instantated cell type is passed to Population.

        @param      extra_receptors : iterable(str)

                    Synaptic mechanism names of synapses that should be allowed
                    on this cell type.

        @post       The list of receptors in self.receptor_types will be
                    updated with each receptor in 'with_receptors' for every
                    cell location/region.
        """
        extra_receptors = kwargs.pop('extra_receptors', None)
        super(MorphCellType, self).__init__(**kwargs)

        # Combine receptors defined on the cell type with regions
        # defined on the model class
        celltype_receptors = type(self).receptor_types
        if extra_receptors is None:
            all_receptors = celltype_receptors
        else:
            all_receptors = celltype_receptors + list(extra_receptors)

        region_receptors = []
        for region in self.model.regions:
            for receptor in all_receptors:
                region_receptors.append(region + "." + receptor)

        self.receptor_types = region_receptors


    def get_schema(self):
        """
        Get mapping of parameter names to allowed parameter types.
        """
        # Avoids specifying default values for each scale parameter and
        # thereby calling the property setter for each of them
        schema = super(MorphCellType, self).get_schema()
        if hasattr(self.model, 'rangevar_names'):
            schema.update({
                varname + '_scale': float for varname in self.model.rangevar_names
            })
        return schema


class MorphModelBase(object):
    """
    Base functionality for our custom PyNN cell models.

    This class implements all methods required to work with our custom
    Connector and Recorder classes or marks them as abstract methods
    for implementation in the subclass.


    ATTRIBUTES
    ----------

    @attr   spike_threshold_source_sec : float
            Spike threshold (mV) in self.source_section

    @attr   synapse_spacing : float
            Synapse spacing (um) for placing synapses

    @attr   regions : list[str]
            Names of subcellular regions

    @attr   region_boundaries : dict[str, tuple[float, float]]
            Boundaries of subcellular regions measured in distance from soma.

    @attr   rangevar_names : list[str]
            List of NEURON RANGE variables that can be scaled

    @attr   rangevar_scaled_seclists : dict[str, list[str]]
            Mapping of RANGE variable name to a section list name
            that's available as a self.icell attribute. By default, scaling
            is applied to icell.all.


    DEVNOTES
    --------

    Called by ID._build_cell() defined in module pyNN/neuron/simulator.py

    - Population._create_cells()
        - Population.cell_type.model(**cell_parameters)
            - MorphModelBase.__init__()
                - instantiate()
    """

    # Combined with celltype.receptors in MorphCellType constructor
    # to make celltype.receptor_types in format 'region.receptor'
    regions = []
    multicompartment = True

    def __init__(self, *args, **kwargs):
        """
        As opposed to the original CellModel class,
        this class instantiates the cell in its __init__ method.

        ARGUMENTS
        ---------

        @param      **kwargs : dict(str, object)
                    Parameter name, value pairs

        CODE CONTRACT
        -------------

        @pre        self.parameter_names is enumerable containing parameter
                    names passed to constructor

        @post       self.icell contains the instantiated Hoc cell model

        """
        # Pass param names from cell_type.(default_parameters + extra_parameters)
        # self.parameter_names must be defined in the subclass body.
        self.parameter_names.append('owning_gid') # gid of owning ID object

        # Make parameters accessible as attributes
        for param_name in list(kwargs.keys()):

            # Check if we should handle the parameter
            if param_name not in self.parameter_names:
                raise ValueError('Unknown parameter "{}":'.format(
                    param_name), kwargs[param_name])

            elif param_name.endswith('_scale'):
                continue # handle after instantiation

            # Set the parameter
            param_value = kwargs.pop(param_name)
            if isinstance(param_value, ArrayParameter):
                param_value = param_value.value

            setattr(self, param_name, param_value)

        # Support for multi-compartment connections
        self.region_to_gid = {}
        self.region_to_source = {}
        self.region_to_section = {}

        # Random number generation for cell
        self.rng_seed = self.owning_gid
        self.rng_numpy = np.random.RandomState(self.rng_seed)

        # Instantiate cell and synapses in NEURON
        self.axon = None
        self.instantiate(sim=ephys_sim_from_pynn())
        self._init_synapses()

        # Handle post-instantiation parameters
        for param_name, param_value in kwargs.items():
            if isinstance(param_value, ArrayParameter):
                param_value = param_value.value
            setattr(self, param_name, param_value)

        # Attributes required by PyNN
        self.source_section = self.icell.soma[0]
        self.source = self.icell.soma[0](0.5)._ref_v
        self.rec = h.NetCon(self.source, None,
                            self.get_threshold(), 0.0, 0.0,
                            sec=self.source_section)

        # Update sources for connections
        self.region_to_gid['soma'] = self.owning_gid
        self.region_to_source['soma'] = self.source
        self.region_to_section['soma'] = self.source_section

        # See pyNN.neuron.recording.Recorder._record()
        self.spike_times = h.Vector(0)
        self.traces = {}
        self.recording_time = False


    def instantiate(self, sim=None):
        """
        Instantiate cell in simulator.

        Virtual method. Implement your own instantiate.

        @pre    self.template_name is a Hoc template that is loaded
                by Hoc interpreter

        @note   The call order is:

                - MyCellModel.__init__()
                    `- cell_base.MorphModelBase.__init__()
                     `- MyCellModel.instantiate()
        """
        raise NotImplementedError("Virtual method. Implement your own instantiate")


    def _post_build(self, population, pop_index):
        """
        Hook called after Population object is instantiated. Used to perform additional actions based on Population properties.

        @see        Population._create_cells() -> ID._build_cell()
        """
        self._init_memb_noise(population, pop_index)


    def memb_init(self):
        """
        Set initial values for all variables in this cell.

        @override   memb_init() required by PyNN interface for cell models.
        """
        for sec in self.icell.all:
            for seg in sec:
                seg.v = self.v_init # set using pop.init(v=v_init) or default_initial_values

        noise_init = getattr(self, 'noise_rng_init', None) # function added dynamically
        if noise_init is not None:
            noise_init()


    def get_threshold(self, region='soma'):
        """
        Get spike threshold for self.source variable (usually points to membrane
        potential). This threshold is used when creating NetCon connections.

        @override   get_threshold() required by pyNN interface for cell models.

        @return     threshold : float
        """
        return self.spike_threshold[region]

    # Properties defined for synapse position
    synapse_spacing = None          # (um)
    region_boundaries = {
        'proximal': (None, None),   # (um)
        'distal':   (None, None),   # (um)
    }


    def _init_synapses(self):
        """
        Initialize synapse map.

        By default it just creates an empty synapse map, but the subclass
        can override this to create synapses upon cell creation.

        @post       self._synapses is a nested dict with depth=2 and following
                    structure:
                    { region: {receptors: list({ 'synapse':syn, 'used':int }) } }
        """

        # Create data structure for synapses
        self._synapses = {region: {} for region in self.regions}

        # Indicate to Connector that we don't allow multiple NetCon per synapse
        self.allow_synapse_reuse = False

        # Sample each region uniformly and place synapses there
        synapse_spacing = self.synapse_spacing # [um]

        target_secs = list(self.icell.somatic) + list(self.icell.basal) + \
                      list(self.icell.apical)

        self._cached_region_segments = {}
        for region in self.regions:
            self._cached_region_segments[region] = []

        for sec in target_secs:
            nsyn = np.ceil(sec.L / synapse_spacing)
            for i in range(int(nsyn)):
                seg = sec((i+1.0)/nsyn)
                for region in self.regions:
                    if self.segment_in_region(seg, region):
                        self._cached_region_segments[region].append(seg)


    def segment_in_region(self, segment, region):
        """
        Test if segment is in subcellular region.
        """
        h.distance(0, 0.5, sec=self.icell.soma[0]) # reference for distance measurement
        segment_dist = h.distance(1, segment.x, sec=segment.sec)
        lower, upper = self.region_boundaries[region]
        return lower <= segment_dist <= upper


    # @abstractmethod
    def get_synapses(self, region, receptors, num_contacts, **kwargs):
        """
        Get synapse in subcellular region for given receptors.
        Called by Connector object to get synapse for new connection.

        Parameters
        ---------

        @param      region : str
                    Region descriptor.

        @param      receptors : enumerable(str)
                    Receptor descriptors.

        @param      mark_used : bool
                    Whether synapse is 'consumed' by caller and should be
                    marked as used.

        @param      **kwargs
                    (Unused) extra keyword arguments for use when overriding
                    this method in a subclass.

        Returns
        -------

        @return     synapse, num_used : tuple(nrn.POINT_PROCESS, int)

                    NEURON point process object that can serve as the
                    target of a NetCon object, and the amount of times
                    the synapse has been used before as the target
                    of a NetCon.
        """
        syns = self.make_synapses_cached_region(region, receptors,
                                                num_contacts, **kwargs)
        synmap_key = tuple(sorted(receptors))
        self._synapses[region].setdefault(synmap_key, []).extend(syns)
        return syns


    def make_synapses_cached_region(self, region, receptors, num_synapses, **kwargs):
        """
        Make a new synapse by sampling a cached region.

        @pre    model must have attribute '_cached_region_segments' containing a
                dict[str, list(nrn.Segment)] that maps the region name
                to eligible segments in that region.
        """
        rng = kwargs.get('rng', np.random)
        region_segs = self._cached_region_segments[region]
        seg_ids = rng.choice(len(region_segs),
                             num_synapses, replace=False)
        if 'mechanism' in kwargs:
            mech_name = kwargs['mechanism']
        elif all((rec in ('AMPA', 'NMDA') for rec in receptors)):
            mech_name = self.default_GLU_mechanism
        elif all((rec in ('GABAA', 'GABAB') for rec in receptors)):
            mech_name = self.default_GABA_mechanism

        return [getattr(h, mech_name)(region_segs[i]) for i in seg_ids]


    def make_new_synapse(self, receptors, segment, mechanism=None):
        """
        Make a new synapse that implements given receptors in the
        given segment.

        @see    get_synapse() for documentation
        """
        if mechanism is not None:
            syn = getattr(h, mechanism)(segment)
        elif all((rec in ('AMPA', 'NMDA') for rec in receptors)):
            syn = getattr(h, self.default_GLU_mechanism)(segment)
        elif all((rec in ('GABAA', 'GABAB') for rec in receptors)):
            syn = getattr(h, self.default_GABA_mechanism)(segment)
        else:
            raise ValueError("No synapse mechanism found that implements all "
                             "receptors {}".format(receptors))
        return syn


    # For method ideas look in git history of ephys_models.py
    # def make_synmap_tree():
    # def make_synapse_balanced(self, region, receptors, mechanism=None):
    # def get_existing_synapse(self, region, receptors, mark_used):


    def get_synapse_list(self, region, receptors):
        """
        Get synapses in region that implements all given receptors.

        @param      region : str
                    Region descriptor.

        @param      receptors : enumerable(str)
                    Receptor descriptors.

        @return     synapse_list : list(dict())
                    List of synapse containers
        """
        return next((syns for recs, syns in self._synapses[region].items() if all(
                        (ntr in recs for ntr in receptors))), [])


    def get_synapses_by_mechanism(self, mechanism):
        """
        Get synapses of given NEURON mechanism.

        @param      mechanism : str
                    NEURON mechanism name

        @return     synapse_list : list(synapse)
                    List of synapses in region
        """
        all_syns = sum((synlist for region in self._synapses.values() for recs, synlist in region.items()), [])
        return [syn for syn in all_syns if nrnutil.get_mod_name(syn) == mechanism]


    def resolve_synapses(self, spec):
        """
        Resolve string definition of a synapse or point process
        on this cell (used by Recorder).

        @param      spec : str

                    Synapse specifier in format "mechanism_name[slice]" where
                    'slice' as a slice expression like '::2' or integer.

        @return     list(nrn.POINT_PROCESS)
                    List of NEURON synapse objects.
        """
        matches = re.search(r'^(?P<mechname>\w+)', spec)
        mechname = matches.group('mechname')

        synlist = self.get_synapses_by_mechanism(mechname)
        if len(synlist) == 0:
            return []
        else:
            return configutil.index_with_str(synlist, spec)


    # Bind interface method
    resolve_section = irec_resolve_section


    def __getattr__(self, name):
        """
        Override __getattr__ to support dynamic properties that are not
        explicitly declared.

        The following dynamic attributes are supported:

            - scale factors for NEURON RAMGE variables

        @pre    Subclass must define attribute 'rangevar_names'

        @note   __getattr__ is called as last resort and the default behavior
                is to raise AttributeError
        """
        attr_name_matches = re.search(r'^(\w+)_scale$', name)
        handle_attribute = (attr_name_matches is not None)

        if handle_attribute and (any(
            [v.startswith(attr_name_matches.groups()[0]) for v in self.rangevar_names])):

            # search for NEURON RANGE variable to be scaled
            varname = attr_name_matches.groups()[0]
            scale_factor_name = '_' + varname + '_scale'

            # Only if it has been scaled before does the scale attribute exist
            if not hasattr(self, scale_factor_name):
                return 1.0 # not scaled
            else:
                # Call will bypass this method if attribute exists
                return getattr(self, scale_factor_name)
        else:
            # return super(MorphModelBase, self).__getattr__(name)
            # raise AttributeError
            return self.__getattribute__(name)


    def __setattr__(self, name, value):
        """
        Override __setattr__ to support dynamic properties.

        The following dynamic attributes are supported:

            - scale factors NEURON RANGE variables, by assigning
              '<rangevar>_scale' where <rangevar> is RANGE variable that is
              a member of self.rangevar_names

        @pre    Subclass must define attribute 'rangevar_names'
        """
        # Check if attribute is of form '<rangevar>_scale' where rangevar
        # is a member of self.rangevar_names
        matches = re.search(r'^(\w+)_scale$', name)
        if (matches is None) or (not any(
                [v.startswith(matches.groups()[0]) for v in self.rangevar_names])):
            return super(MorphModelBase, self).__setattr__(name, value)
        varname = matches.groups()[0]
        scale_factor_name = '_' + varname + '_scale'

        # Calculate actual scale factor based on previous scale
        if not hasattr(self, scale_factor_name):
            old_scale = 1.0
        else:
            old_scale = getattr(self, scale_factor_name)
        if value == old_scale:
            return # scale is already correct
        relative_scale = value / old_scale

        # Find NEURON mechanism name
        matches = re.search(r'.+_([a-zA-Z0-9]+)$', varname)
        if matches is None:
            raise ValueError('Could not extract NEURON mechanism name '
                             'from RANGE variable name {}'.format(varname))
        mechname = matches.groups()[0]

        # Scale neuron RANGE variable on sections in target region
        if hasattr(self, 'rangevar_scaled_seclists'):
            scaled_seclists = self.rangevar_scaled_seclists.get(varname, ['all'])
        else:
            scaled_seclists = ['all']

        for seclist_name in scaled_seclists:
            for sec in getattr(self.icell, seclist_name):
                if not h.ismembrane(mechname, sec=sec):
                    continue
                for seg in sec:
                    setattr(seg, varname, relative_scale * getattr(seg, varname))

        # Save scale factor as self._<varname>_scale
        setattr(self, scale_factor_name, value)


    def get_all_sections(self):
        """
        Get all neuron.Section objects that make up this cell.

        @return     neuron.SectionList containing all sections
        """
        return self.icell.all


    def _init_axon(self, axon_class, with_ais_compartment=True):
        """
        Create and append axon.

        @post   self.axon contains references to axonal sections.
        """

        # Build axon
        axonal_secs = list(self.icell.axonal)
        if len(axonal_secs) > 0:
            # Attach axon to axon stub/AIS if present
            axon_terminal_secs = treeutils.leaf_sections(axonal_secs[0], subtree=True)
            assert len(axon_terminal_secs) == 1
            axon_parent_sec = axon_terminal_secs[0]
        else:
            # Attach axon directly to soma
            axon_parent_sec = self.icell.soma[0]

        # Parameters for axon building
        axon_builder = axon_class(
            self.streamline_coordinates_mm,
            termination_method=self.termination_method,
            interp_method='arclength',
            parent_cell=self.icell,
            parent_sec=axon_parent_sec,
            connection_method='translate_axon_closest',
            tolerance_mm=1e-4,
            without_extracellular=not (self.with_extracellular_rec or self.with_extracellular_stim),
            rng=self.rng_numpy)

        if not with_ais_compartment:
            # If cell model already has AIS, remove it from axon model
            axon_builder.initial_comp_sequence.pop(0) # = ['aismyelin']

        # Axon collaterals are specified as Nx3 matrix of branch points
        self.with_collaterals = (hasattr(self, 'collateral_branch_points_um') and
            not all([n[0] == 0 for n in self.collateral_lvl_num_branches]) and
            not all(np.isnan(self.collateral_branch_points_um).flatten())) # NaN when None is specified
        if self.with_collaterals:
            assert self.collateral_branch_points_um.ndim == 2
            assert self.collateral_branch_points_um.shape[1] == 3
            for i, branch_point_um in enumerate(self.collateral_branch_points_um):
                target_point_um = self.collateral_target_points_um[i]
                level_lengths_um = self.collateral_lvl_lengths_um[i]
                levels_num_branches = self.collateral_lvl_num_branches[i]
                step_length_um = 10.0
                step_angles_deg = np.zeros_like(level_lengths_um) + 30.0
                axon_builder.add_collateral_definition(
                                branch_point_um, target_point_um, level_lengths_um,
                                step_length_um, levels_num_branches, step_angles_deg)
        self.axon = axon_builder.build_axon()

        # Change source for NetCons (see pyNN.neuron.simulator code)
        terminal_sec = axon_builder.get_terminal_section(self.netcon_source_spec)
        terminal_source = terminal_sec(0.5)._ref_v # source for connections
        self.source_section = terminal_sec
        self.source = terminal_source

        # Support for multicompartment connections
        region = 'axon_terminal'
        self.region_to_gid[region] = None # mark as unset
        self.region_to_source[region] = terminal_source
        self.region_to_section[region] = terminal_sec


    _init_extracellular_stim_rec = init_extracellular_stim_rec


    def _init_emfield(self):
        """
        Set up extracelullar stimulation and recording.

        @pre    mechanism 'extracullular' must be inserted and parameters set
                in all compartments that should contribute to the LFP and are
                targets for stimulation
        """
        # Check if we need extracellular layers for recording/stimulation
        if not (self.with_extracellular_rec or self.with_extracellular_stim):
            return

        # Get section lists with recorded / stimulated sections
        seclists_with_dbs = set(getattr(self, 'seclists_with_dbs', '').split(';'))
        seclists_with_dbs.difference_update({''})
        if self.with_extracellular_stim and len(seclists_with_dbs) == 0:
            seclists_with_dbs = {'all'}

        seclists_with_lfp = set(getattr(self, 'seclists_with_lfp', '').split(';'))
        seclists_with_lfp.difference_update({''})
        if self.with_extracellular_rec and len('seclists_with_lfp') == 0:
            seclists_with_lfp = {'all'}

        seclists_with_xtra = seclists_with_dbs.union(seclists_with_lfp)

        # Make Hoc.SectionList with stimulated / recorded sections
        if not self.with_extracellular_stim:
            stim_seclist = h.SectionList()
        elif 'all' in seclists_with_dbs:
            stim_seclist = self.icell.all
        else:
            stim_seclist = nrnutil.join_seclists(
                *(getattr(self.icell, sl) for sl in seclists_with_dbs))

        if not self.with_extracellular_rec:
            rec_seclist = h.SectionList()
        elif 'all' in seclists_with_lfp:
            rec_seclist = self.icell.all
        elif seclists_with_dbs == seclists_with_lfp:
            rec_seclist = stim_seclist
        else:
            rec_seclist = nrnutil.join_seclists(
                *(getattr(self.icell, sl) for sl in seclists_with_lfp))

        # Insert mechanism 'extracellular' providing extracellular layers
        extra_layer_indices = range(2)
        for seclist_name in seclists_with_xtra:
            seclist = getattr(self.icell, seclist_name)
            for sec in seclist:
                if h.ismembrane('extracellular', sec=sec):
                    continue # assume already configured

                sec.insert('extracellular')
                for i in extra_layer_indices:
                    sec.xraxial[i] = 1e9
                    sec.xg[i] = 1e9
                    sec.xc[i] = 0.0

        self._init_extracellular_stim_rec(self.icell.all, stim_seclist, rec_seclist,
                                          self.icell.soma[0](0.5))


    def _init_memb_noise(self, population, pop_index):
        # Insert membrane noise
        if self.membrane_noise_std > 0:
            # Configure RNG to generate independent stream of random numbers.
            num_picks = int(pyNN.neuron.state.duration / pyNN.neuron.state.dt)
            seed = (1e4 * population.pop_gid) + pop_index
            rng, init_rng = nrnutil.independent_random_stream(
                                    num_picks, pyNN.neuron.state.mcellran4_rng_indices,
                                    force_low_index=seed)
            rng.normal(0, 1)
            self.noise_rng = rng
            self.noise_rng_init = init_rng

            soma = self.icell.soma[0]
            self.noise_stim = stim = h.ingauss2(soma(0.5))
            std_scale =  1e-2 * sum((seg.area() for seg in soma)) # [mA/cm2] to [nA]
            stim.mean = 0.0
            stim.stdev = self.membrane_noise_std * std_scale
            stim.noiseFromRandom(rng)
        else:
            def fdummy():
                pass
            self.noise_rng = None
            self.noise_rng_init = fdummy

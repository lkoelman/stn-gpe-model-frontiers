"""
Reduction of Balbi et al. motoneuron model using CERSEI folding tools
"""

import re
import sys

pkg_root = ".." # root dir for our packages
sys.path.append(pkg_root)

import cersei.collapse.redutils as redutils
import cersei.collapse.cluster as cluster
import bgcellmodels.common.logutils as logutils
import bgcellmodels.common.treeutils as treeutils
from bgcellmodels.common.nrnutil import ExtSecRef, getsecref
from cersei.collapse.fold_reduction import ReductionMethod, FoldReduction

import balbi_model
from neuron import h

import logging
logger = logging.getLogger('balbi')
logutils.setLogLevel('verbose', ['balbi', 'marasco', 'folding'])

################################################################################
# Cell model-specific implementations of reduction functions
################################################################################

class BalbiFoldReduction(FoldReduction):
    """
    Model-specific functions for folding reduction of Balbi et al. (2015)
    motoneuron model.
    """

    def __init__(self, **kwargs):
        """
        Initialize new Balbi cell model reduction.

        @param  balbi_motocell_id   (int) morphology file identifier
        """
        self.balbi_motocell_id = kwargs.pop('balbi_motocell_id')
        named_secs = balbi_model.get_named_sec_lists()

        # Group subtrees by trunk
        somatic = list(named_secs['soma'])
        # Dendritic sections have other subtrees (trunks = first level branchpoints)
        dendritic = list(named_secs['dend'])
        # Axonic sections are first subtree (trunk = AH)
        axonic = sum((list(named_secs[name]) for 
                        name in ('AH', 'IS', 'node', 'MYSA', 'FLUT', 'STIN')), [])
        
        nonsomatic = axonic + dendritic
        kwargs['soma_secs'] = somatic
        kwargs['dend_secs'] = nonsomatic

        # Whole cell
        wholecell = treeutils.wholetree_secs(somatic[0]) # sum(named_secs.values(), [])
        kwargs['fold_root_secs'] = self._get_fold_roots(named_secs, wholecell)

        # Set all parameters
        kwargs['gleak_name'] = balbi_model.gleak_name
        kwargs['mechs_gbars_dict'] = balbi_model.gbar_dict
        kwargs['mechs_params_nogbar'] = balbi_model.mechs_params_nogbar
        super(BalbiFoldReduction, self).__init__(**kwargs)


    def _get_fold_roots(reduction, named_secs, all_secs):
        """
        Get folding roots.
        """
        # Get the folding branchpoints for the dendritic subtrees
        all_refs = [ExtSecRef(sec=sec) for sec in all_secs]
        root_sec = all_refs[0].root; h.pop_section() # true root, pushes CAS
        root_ref = getsecref(root_sec, all_refs)

        # Choose root sections for folding operations
        cluster.assign_topology_attrs(root_ref, all_refs)
        trunk_refs = redutils.find_roots_at_level(2, all_refs) # level 2 roots
        
        # At which nodes should tree be folded?
        fold_roots = [] # [named_secs['AH'][0]]
        fold_roots.extend([ref.sec for ref in trunk_refs])
        logger.debug("Folding roots:\n" + "\n".join((str(ref) for ref in fold_roots)))

        return fold_roots


    def assign_region_label(reduction, secref):
        """
        Assign region labels to sections.
        """
        arrayname = re.sub(r"\[\d+\]", "", secref.sec.name())
        if secref.is_original:
            if arrayname == 'soma':
                secref.region_label = 'somatic'
            elif arrayname == 'dend':
                secref.region_label = 'dendritic'
            elif arrayname in ['AH', 'IS', 'node', 'MYSA', 'FLUT', 'STIN']:
                secref.region_label = 'axonic'
            else:
                raise Exception("Unrecognized original section {}".format(secref.sec))
        
        elif hasattr(secref, 'merged_region_labels'):
            secref.region_label = '-'.join(sorted(secref.merged_region_labels))

        elif hasattr(secref, 'orig_props'):
            secref.region_label = '-'.join(sorted(secref.orig_props.merged_region_labels))


    def fix_section_properties(self, new_sec_refs):
        """
        Fix properties of newly created sections.

        @override   FoldReduction.fix_section_properties
        """
        # Set global mechanism parameters again
        h.fix_global_params() # defined in 3_ins_ch.hoc

        # Set region-specific properties
        for ref in new_sec_refs:
            if not 'axonic' in ref.region_label:
                ref.sec.insert('extracellular')
                nlayer = 2 # see https://neuron.yale.edu/neuron/static/new_doc/modelspec/programmatic/mechanisms/mech.html -> 'extracellular'
                for i in range(nlayer):
                    ref.sec.xraxial[i] = 1e9
                    ref.sec.xg[i] = 1e10
                    ref.sec.xc[i] = 0
            else:
                raise Exception('Axonic sections should not be collapsed. '
                                '(section {}'.format(ref.sec))


    @staticmethod
    def init_cell_steadystate():
        """
        Initialize cell for analyzing electrical properties.
        
        @note   Ideally, we should restore the following, like SaveState.restore():
                - all mechanism STATE variables
                - voltage for all segments (seg.v)
                - ion concentrations (nao,nai,ko,ki, ...)
                - reversal potentials (ena,ek, ...)
        """
        h.celsius = 37

        # NOTE: cannot use h.load_steadystate(). This uses SaveState.restore() which means 
        #       between a save and a restore, you cannot create or delete sections, 
        #       NetCon objects, or point processes, nor change the number of segments, 
        #       insert or delete mechanisms, or change the location of point processes.
        
        # Uniform ion concentrations, verified in all sections after h.load_steadystate()
        h.nai0_na_ion = 10.0
        h.nao0_na_ion = 140.0
        h.ki0_k_ion = 54.0
        h.ko0_k_ion = 2.5
        h.cai0_ca_ion = 5e-5
        h.cao0_ca_ion = 2.0

        h.init()


################################################################################
# Gillies Model Reduction Experiments
################################################################################

def make_fold_reduction():
    """
    Make FoldReduction object with Marasco method.
    """
    # Instantiate NEURON model
    BALBI_CELL_ID = 1 # morphology file to load
    balbi_model.make_cell_balbi(model_no=BALBI_CELL_ID)
    

    # Reduce model
    reduction = BalbiFoldReduction(method=ReductionMethod.BushSejnowski,
                                   balbi_motocell_id=BALBI_CELL_ID)

    # Extra Reduction parameters
    reduction.set_reduction_params({
            'Z_freq' :              25.,
            'Z_init_func' :         reduction.init_cell_steadystate,
            'Z_linearize_gating' :  False,
            'f_lambda':             100.0,
            'gbar_scaling' :        'area',
            'syn_map_method' :      'Ztransfer',
            'post_tweak_funcs' :    [],
        })

    return reduction


def reduce_model_folding(export_locals=True):
    """
    Reduce cell model using given folding algorithm
    
    @param  export_locals       if True, local variables will be exported to the global
                                namespace for easy inspection
    """
    from timeit import default_timer as timer

    # Make reduction object
    tstart = timer()
    reduction = make_fold_reduction()
    tstop = timer()
    logger.debug("Building model took {}".format(tstop-tstart))
    
    # Do reduction
    reduction.reduce_model(num_passes=1) # SETPARAM: numbe of folding passes
    tstart, tstop = tstop, timer()
    logger.debug("Reducing model took {}".format(tstop-tstart))

    if export_locals:
        globals().update(locals())

    return reduction._soma_refs, reduction._dend_refs


if __name__ == '__main__':
    reduce_model_folding()

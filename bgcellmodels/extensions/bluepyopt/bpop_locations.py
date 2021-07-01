"""
Extensions to Ephys.locations module

@author Lucas Koelman

@date   8/05/2018
"""

from neuron import h
from bgcellmodels.common import treeutils

from bluepyopt.ephys.locations import Location
from bluepyopt.ephys.serializer import DictMixin

import re
import logging
logger = logging.getLogger('bpop_ext')


class NrnNamedSecLocation(Location, DictMixin):
    """
    Location in a specific section, identified by name.
    """

    SERIALIZED_FIELDS = (
        'name',
        'comment',
        'seclist_name',
        'sec_index',
        'comp_x',
    )

    def __init__(
            self,
            name,
            sec_name=None,
            comp_x=None,
            comment=''):
        """
        Constructor
        
        Args:
            name (str): name of the object
            seclist_name (str): name of Neuron section list (ex: 'somatic')
            sec_index (int): index of the section in the section list
            comp_x (float): segx (0..1) of segment inside section
        """

        super(NrnNamedSecLocation, self).__init__(name, comment)
        self.sec_name = sec_name
        self.comp_x = comp_x

    def instantiate(self, sim=None, icell=None):  # pylint: disable=W0613
        """
        Find the instantiate compartment
        """
        iseclist = icell.all
        
        try:
            isection = next((sec for sec in iseclist if sec.name()==self.sec_name))
        except StopIteration:
            raise ValueError("Section with name {} not found on this cell".format(self.sec_name))
        icomp = isection(self.comp_x)
        return icomp

    def __str__(self):
        """String representation"""
        return '%s(%s)' % (self.sec_name, self.comp_x)


class NrnSeclistLocationExt(Location, DictMixin):
    """
    Section in a sectionlist
    """

    SERIALIZED_FIELDS = (
        'name', 
        'comment', 
        'seclist_name',
        'section_filter'
    )

    def __init__(
            self,
            name,
            seclist_name=None,
            comment='',
            secname_filter=None):
        """
        Constructor
        
        Args:
            name (str): name of the object
            seclist_name (str): name of NEURON section list (ex: 'somatic')
            secname_filter: regular expression to match section name

        NOTE: can't serialize functions, unless we change DixtMixon to make 
              use of external library for serialization, e.g. Dill or Pyro.
        """

        super(NrnSeclistLocationExt, self).__init__(name, comment)
        self.seclist_name = seclist_name
        
        if secname_filter is None:
            secname_filter = r'' # matches any string
        self.secname_filter = secname_filter

    def instantiate(self, sim=None, icell=None):  # pylint: disable=W0613
        """
        Find the instantiate compartment
        """

        isectionlist = getattr(icell, self.seclist_name)

        return (isection for isection in isectionlist if re.search(self.secname_filter, isection.name()))

    def __str__(self):
        """
        String representation
        """

        return '%s' % (self.seclist_name)


class SomaDistDiamLocSimple(Location, DictMixin):
    """
    All sections in a SectionList with diameter in a given range
    that lie in a given distance range from the soma compartment.
    """
    
    SERIALIZED_FIELDS = (
        'name',
        'comment',
        'seclist_name',
        'distance_range',
        'diameter_range')

    def __init__(
            self,
            name,
            seclist_name=None,
            distance_range=None,
            diameter_range=None,
            comment=''):
        """
        Constructor

        @param      name : str
                    Name of this object

        @param      seclist_name : str
                    Name of Neuron section list (ex: 'apical')

        @param      distance_range : tuple(float, float)

        @param      diameter_range : tuple(float, float)
        """

        super(SomaDistDiamLocSimple, self).__init__(
            name,
            comment=comment)

        self.distance_range = distance_range
        self.diameter_range = diameter_range
        self.seclist_name = seclist_name


    def instantiate(self, sim=None, icell=None):
        """
        Get all Sections in named SectionList within requested distance range
        and diameter range.
        """

        isegments = []
        target_seclist = list(getattr(icell, self.seclist_name))

        # Initialize distance measurement at soma
        somatic = list(icell.somatic)
        # trunk_sections = [ch for sec in somatic for ch in sec.children() if ch not in somatic]
        somaref = h.SectionRef(sec=somatic[0])
        soma = somaref.root; h.pop_section() # changes the cas
        assert soma in somatic, "Root Section of tree is not somatic"
        h.distance(0, 0.5, sec=soma)
        soma_offset = soma.L / 2.0

        # do tree traversal
        for seg in treeutils.ascend_segmentwise_dfs(soma):
            if seg.sec not in target_seclist:
                continue

            # first check distance constraint
            sec_dist = max(0.0, h.distance(seg.x, sec=seg.sec) - soma_offset)
            lower, upper = self.distance_range
            if not (lower < sec_dist <= upper):
                continue

            # then check diameter constraint
            lower, upper = self.diameter_range
            if not (lower < seg.diam <= upper):
                continue

            isegments.append(seg)
            print("Add {} to distance {} and diam {}".format(
                seg, self.distance_range, self.diameter_range))


        return (seg for seg in isegments)


    def __str__(self):
        """
        String representation of location.
        """
        return 'Distance range {}-{} micron from soma in {} with diameter {}-{} um'.format(
            self.distance_range[0], self.distance_range[1],
            self.seclist_name,
            self.diameter_range[0], self.diameter_range[1])


class SomaDistanceDiamLocation(Location, DictMixin):
    """
    All sections in a SectionList with diameter in a given range
    that lie in a given distance range from the soma compartment.
    """
    
    SERIALIZED_FIELDS = (
        'name',
        'comment',
        'seclist_name',
        'distance_range',
        'diameter_range')

    def __init__(
            self,
            name,
            seclist_name=None,
            distance_range=None,
            diameter_range=None,
            comment=''):
        """
        Constructor

        @param      name : str
                    Name of this object

        @param      seclist_name : str
                    Name of Neuron section list (ex: 'apical')

        @param      distance_range : tuple(float, float)

        @param      diameter_range : tuple(float, float)
        """

        super(SomaDistanceDiamLocation, self).__init__(
            name,
            comment=comment)

        self.distance_range = distance_range
        self.diameter_range = diameter_range
        self.seclist_name = seclist_name


    def instantiate(self, sim=None, icell=None):
        """
        Get all Sections in named SectionList within requested distance range
        and diameter range.

        Whether a section falls within the diameter range depends on the
        transition rules defined by J.Edgerton.

        ALGORITHM
        ---------

        - start at root
        - assign initial rating
        - ascend tree: go to next section
            - assign temporary rating based on section's diam
            - look ahead for 4 sections, along each possible path (depth-first)
                - if each path has <= 1 higher rating, keep temporary rating
                - if a path has > 1 higher rating, use parent rating
        """

        isegments = []
        target_seclist = list(getattr(icell, self.seclist_name))

        # Initialize distance measurement at soma
        somatic = list(icell.somatic)
        # trunk_sections = [ch for sec in somatic for ch in sec.children() if ch not in somatic]
        assert len(somatic) == 1, "more than one somatic section"
        soma = somatic[0]
        h.distance(0, 0.5, sec=soma)
        soma_offset = soma.L / 2.0

        # do tree traversal
        for root_sec in soma.children():
            parent_rating = 1 # initialize to above diameter range
            for seg in treeutils.ascend_segmentwise_dfs(root_sec):
                # if sec in trunk_sections or sec in somatic:
                #   parent_rating = self.diameter_rating(sec.diam)

                # check diameter constraint
                temp_rating = self.diameter_rating(seg.diam)
                if parent_rating < 0 or temp_rating < 0:
                    # Monotonically decreasing constraint
                    diam_rating = -1
                elif parent_rating == 0 and temp_rating >= 0:
                    # Monotonically decreasing constraint
                    diam_rating = 0
                elif parent_rating > 0 and temp_rating == 0:
                    # transition
                    if self.children_within_range(seg, 4, 1):
                        diam_rating = 0
                    else:
                        diam_rating = parent_rating
                else: # parent_rating > 0 and temp_rating > 0
                    diam_rating = parent_rating
                parent_rating = diam_rating

                if seg.sec not in target_seclist:
                    # print("Not in seclist: {}".format(sec.name()))
                    continue

                # check distance constraint
                sec_dist = max(0.0, h.distance(seg.x, sec=seg.sec) - soma_offset)
                lower, upper = self.distance_range
                if (lower < sec_dist <= upper) and (diam_rating == 0):
                    isegments.append(seg)


        return (sec for sec in isegments)


    def diameter_rating(self, diam):
        """
        Check whether given diameter is within range, above it, or below it.

        @return     relative_range : int
                    +1 if diameter is above range
                    0 if diamtere is within range
                    -1 if diameter is below range
        """
        lower, upper = self.diameter_range
        if diam <= lower:
            return -1
        elif diam <= upper:
            return 0
        else:
            return 1


    def children_within_range(self, seg, max_depth, max_greater, depth=0, num_greater=0):
        """
        Check that all paths along child sections up to max_depth have
        maximum <max_greater> diameters that are above the given diameter range.

        @param  sec : nrn.Section
                starting section

        @pram   max_depth : int
                Maximum depth of ascent / number of sections from starting section.
        """
        if depth > max_depth:
            return True
        if num_greater > max_greater:
            return False

        if depth != 0 and self.diameter_rating(seg.diam) > 0:
            num_greater += 1

        for child_seg in treeutils.next_segs(seg):
            if not self.children_within_range(child_seg, max_depth, max_greater, 
                                              depth+1, num_greater):
                return False
        return True


    def __str__(self):
        """
        String representation of location.
        """
        return 'Distance range {}-{} micron from soma in {} with diameter {}-{} um'.format(
            self.distance_range[0], self.distance_range[1],
            self.seclist_name,
            self.diameter_range[0], self.diameter_range[1])
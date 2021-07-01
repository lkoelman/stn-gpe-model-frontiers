"""
Classes describing the ontology of Basal Ganglia physiology and histology.

@author     Lucas Koelman
@date       12/07/2017
"""

from enum import unique
from bgcellmodels.common.stdutil import IntEnumDescriptor


@unique
class PhysioState(IntEnumDescriptor):
    """
    Physiological state of the cell
    """

    # bit flags for states
    PARKINSONIAN =      2**0
    DBS =               2**1
    AWAKE =             2**2
    SLEEP_SWS =         2**3

    # Combinations
    NORMAL = 0
    PARK_DBS = PARKINSONIAN | DBS
    # ALL = 2**32 - 1

    def __contains__(self, item):
        """
        Bitwise operation to detect if item is subset/member of this state

        NOTE: this only works when item is an Enum item, e.g. check using
              expression (item_A in item_B)
        """
        return  (self.value & item.value) == item.value

    def is_subset(self, int_item):
        """
        Check if this item is a subset of item described by bitflags
        in given integer 'int_item'

        NOTE: this can be used for checking bitflags in an Enum item
              againt bitflags in an aribrary integer, e.g.: item.is_subset(8-1)
        """
        return (self.value & int_item) == self.value


@unique
class Populations(IntEnumDescriptor):
    """
    Physiological state of the cell
    """
    STN = 0
    CTX = 1
    GPE = 2
    THA = 4
    PPN = 5
    STR = 6


@unique
class NTReceptors(IntEnumDescriptor):
    """
    NTReceptors used in synaptic connections
    """
    AMPA = 0
    NMDA = 1
    GABAA = 2
    GABAB = 3


@unique
class ParameterSource(IntEnumDescriptor):
    """
    A source (e.g. scientific article) for physiological parameters.
    """
    Default = 0
    CommonUse = 1       # Widely used or commonly accepted values
    Chu2015 = 2
    Baufreton2009 = 3
    Fan2012 = 4
    Atherton2013 = 5
    Kumaravelu2016 = 6
    Custom = 7          # Custom parameters, e.g. for testing
    Gradinaru2009 = 8
    Mallet2016 = 9 # Mallet (2016), Neuron 89
    Nambu2014 = 10
    Bergman2015RetiCh3 = 11
    Hendrickson2011 = 12 # Paper limitations of reduced compartmental models
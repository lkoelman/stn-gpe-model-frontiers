"""
Common functionality for folding algorithms.

@author     Lucas Koelman

@date       04-01-2018
"""

from enum import Enum, unique
from abc import ABCMeta, abstractmethod # abstract base class

@unique
class ReductionMethod(Enum):
    Rall = 0
    Stratford = 1           # Stratford, K., Mason, A., Larkman, A., Major, G., and Jack, J. J. B. (1989) - The modelling of pyramidal neurones in the visual cortex
    BushSejnowski = 2       # Bush, P. C. & Sejnowski, T. J. Reduced compartmental models of neocortical pyramidal cells. Journal of Neuroscience Methods 46, 159-166 (1993).
    Marasco = 3             # Marasco, A., Limongiello, A. & Migliore, M. Fast and accurate low-dimensional reduction of biophysically detailed neuron models. Scientific Reports 2, (2012).

    @classmethod
    def from_str(cls, descr):
        method_str = descr.lower()
        if method_str in ['rall']:
            return ReductionMethod.Rall
        elif method_str in ['stratford', 'stratfordmason', 'stratfordmasonlarkman']:
            return ReductionMethod.Stratford
        elif method_str in ['bush', 'bushsejnowski']:
            return ReductionMethod.BushSejnowski
        elif method_str in ['marasco', 'marascolimongiello', 'marascolimongiellomigliore']:
            return ReductionMethod.Marasco
        else:
            return ValueError('Unrecognized reduction method {}'.format(descr))


class FoldingAlgorithm(object):
    """
    Abstract base class for folding algorithms.
    """

    __metaclass__ = ABCMeta # Python2 way


    @abstractmethod
    def preprocess_reduction(self):
        """
        Preprocess cell for folding reduction.

        Calculates properties needed during reduction and saves them on reduction
        object and individual SectionRef instances.

        @param  reduction   FoldReduction object

        @effect             - assign identifiers to each sections
                            - compute electrotonic properties in each segment of cell
                            - determine interpolation path for channel distributions and save them
        """
        raise NotImplementedError(
                "Virtual method of abstract base class FoldingAlgorithm not implemented.")


    @abstractmethod
    def fold_one_pass(self, i_pass, Y_criterion):
        """
        Collapse branches at branch points identified by given criterion.
        """
        raise NotImplementedError(
                "Virtual method of abstract base class FoldingAlgorithm not implemented.")


    @abstractmethod
    def postprocess_reduction(self):
        """
        Postprocess cell after reduction procedure. Executed once.
        """
        raise NotImplementedError(
                "Virtual method of abstract base class FoldingAlgorithm not implemented.")
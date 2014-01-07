"""
GeneTree acessibility/creation functions for CodingBlockGraphs
"""

# Module metadata
__authors__ = "Ate van der Burgt"
__license__ = "MIT"

# graphAbgp Imports
import graphAbgp.conversion

# Abgp imports

# Python Imports

# Global Variable Import


class GeneTreeOfCodingBlockFunctions:
    """
    Subclass with GeneTreeGraph functions for CodingBlockGraphs
    
    @attention: not functional on its own, required to be inherited from
    """
    def genetree(self):
        """
        Get GeneTreeGraph of this CodingBlockGraph

        @attention: alias for self.get_genetree()

        @rtype:  GeneTreeGraph
        @return: GeneTreeGraph instance
        """
        if self._GENETREE:
            return self._GENETREE
        else:
            return self.set_genetree()

    # end of function genetree


    def get_genetree(self):
        """
        Get GeneTreeGraph of this CodingBlockGraph

        @attention: alias for self.genetree()

        @rtype:  GeneTreeGraph
        @return: GeneTreeGraph instance
        """
        if self._GENETREE:
            return self._GENETREE
        else:
            return self.set_genetree()

    # end of function get_genetree


    def set_genetree(self):
        """
        Convert CodingBlockGraph 2 GeneTree

        @rtype:  GeneTreeGraph
        @return: GeneTreeGraph instance
        """
        self._GENETREE = conversion.CodingBlockGraph2GeneTreeGraph(self)
        return self._GENETREE

    # end of function set_genetree

# end of class GeneTreeOfCodingBlockFunctions

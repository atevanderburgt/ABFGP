"""
GeneTree acessibility/creation functions for CodingBlockGraph like classes
"""

# graphAbgp imports
import conversion

# Python Imports

class GeneTreeOfCodingBlockFunctions:
    """
    """
    def genetree(self):
        """
        Get GeneTreeGraph of this CodingBlockGraph

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

        @attention: THIS FUNCTION MUST BE OVERWRITTEN IN NON-CodingBlockGraph classes!
        """
        self._GENETREE = conversion.CodingBlockGraph2GeneTreeGraph(self)
        return self._GENETREE

    # end of function set_genetree

# end of class GeneTreeOfCodingBlockFunctions

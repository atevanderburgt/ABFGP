"""
Accessing EMBOSS-TCODE data in Graphs joined in a class that is inherited in graph classes for Aligment Based Gene Predictions

In __init__() function of graph objects, some specific attributes must be set:
    self._tcode5pscore = {}
    self._tcode3pscore = {}
    self._TCODE_5P_WINDOWSIZE = tcode_5p_windowsize
    self._TCODE_3P_WINDOWSIZE = tcode_3p_windowsize

"""

# Module metadata
__authors__ = "Ate van der Burgt"
__license__ = "MIT"


# Python Imports

class GraphTcodeDataAccesFunctions:
    """
    """
    def average_5p_tcode_score(self):
        """ """
        if not self.get_nodes():
            return 0.0 # empty graph
        else:
            return sum(self._tcode5pscore.values()) / self.organism_set_size()

    # end of function average_5p_tcode_score


    def average_3p_tcode_score(self):
        """ """
        if not self.get_nodes():
            return 0.0 # empty graph
        else:
            return sum(self._tcode3pscore.values()) / self.organism_set_size()

    # end of function average_3p_tcode_score

# end of class GraphTcodeDataAccesFunctions

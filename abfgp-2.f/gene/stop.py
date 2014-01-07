"""
"""
# Module metadata
__authors__ = "Ate van der Burgt"
__license__ = "MIT"

# Imports
from gene_gff import BasicGFF
from gene_exceptions import InproperlyAppliedArgument

# Import Global variables


class StopCodon(BasicGFF):
    def __init__(self,pos,gff={}):
        """ """
        BasicGFF.__init__(self)
        self._gff.update(gff)
        self.pos        = pos
        self.start      = self.pos
        self.end        = self.start+3
        self.pssm_score = 1.0  # default, dummy value
        self.phase      = 0
    # end of function __init__

# end of class StopCodon


def _handle_stopcodon(stop):
    """ """
    if stop.__class__.__name__ in ['StopCodon']:
        stopcodon = stop
    elif type(stop) == type(int()):
        stopcodon = StopCodon(stop)
    else:
        raise InproperlyAppliedArgument, "``stop`` not an expected object or integer"
    # return stop object
    return stopcodon

# end of function _handle_stopcodon

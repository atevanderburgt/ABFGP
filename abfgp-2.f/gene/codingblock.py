"""
"""
# Module metadata
__authors__ = "Ate van der Burgt"
__license__ = "MIT"

# Imports
from gene_gff import BasicGFF

# Import Global variables


class CodingBlockStart(BasicGFF):
    def __init__(self,pos,gff={}):
        """ """
        BasicGFF.__init__(self)
        self._gff.update(gff)       
        self.pos        = pos
        self.start      = self.pos
        self.end        = self.start+2
        self.pssm_score = 1.0  # default, dummy value
        self.phase      = None # no phase for this feature!

    # end of function __init__

    def __str__(self):
        """ """
        return "<%s pos=%s>" % (self.__class__.__name__,self.pos)
    # end of function __str__

# end of class CodingBlockStart


class CodingBlockEnd(BasicGFF):
    def __init__(self,pos,gff={}):
        """ """
        BasicGFF.__init__(self)
        self._gff.update(gff)
        self.pos        = pos
        self.start      = self.pos-2
        self.end        = self.pos
        self.pssm_score = 1.0  # default, dummy value
        self.phase      = None # no phase for this feature!

    # end of function __init__

    def __str__(self):
        """ """
        return "<%s pos=%s>" % (self.__class__.__name__,self.pos)
    # end of function __str__

# end of class CodingBlockEnd

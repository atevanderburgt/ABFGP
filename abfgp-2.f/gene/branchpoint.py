"""
"""
# Module metadata
__authors__ = "Ate van der Burgt"
__license__ = "MIT"

# Python Imports
from re import finditer, compile, IGNORECASE
BP_PATTERN = compile("[CT]T[AG]A[CT]",IGNORECASE)
BP_PATTERN_SINGLE_MISMATCH = compile("([^CT]T[AG]A[CT]|[CT]T[^AG]A[CT]|[CT]T[AG]A[^CT])",IGNORECASE)


# Imports
from gene_gff import BasicGFF
from gene_exceptions import InproperlyAppliedArgument

# Import Global variables

class BranchPoint(BasicGFF):
    """ """
    def __init__(self,start,end,pattern="",gff={}):
        """ """
        BasicGFF.__init__(self)
        self._gff.update(gff)
        self.start      = start
        self.end        = end
        self.pos        = self.start
        self.pattern    = pattern

    # end of function __init__

    #def __str__(self):
    #    return BasicGFF.__str__(self).replace(">"," [%s]>" % self._gff['fscore'])
    ## end of function __str__

    def __str__(self):
        """ print one-liner with BranchPoint data """
        _pattern = ""
        _score   = ""
        if hasattr(self,"pattern") and self.pattern != None:
            _pattern = "%s" % self.pattern
            _score   = " "
        if self._gff['fscore'] not in ("",".",None):
            _score = " (%2.1f) " % self._gff['fscore']
        return "<%s %s-%s%s%s>" % (
            self.__class__.__name__,
            self.start, self.end,
            _score, _pattern
            )
    # end of function __str__

# end of class BranchPoint


def scan_regex_branchpoint(seq,allow_single_mismatch=False):
    """
    Find BranchPoint sequences by regular expression of the input sequence

    @type seq:  string
    @param seq: DNA sequence

    @rtype:  list
    @return: list with BranchPoint objects
    """
    base2score = { 'T':1.0, 'G':0.5, 'C':0.0, 'A': -2.0 }
    sites = []
    for m in finditer(BP_PATTERN,seq):
        # instantiate BranchPoint object
        bp = BranchPoint(m.start(),m.end(),pattern=m.group())
        bp._gff['fscore'] = 1
        sites.append(bp)
    if allow_single_mismatch:
        for m in finditer(BP_PATTERN_SINGLE_MISMATCH,seq):
            # instantiate BranchPoint object
            seq = list(m.group())
            if seq[0] not in ['C','T']:
                seq[0] = seq[0].lower()
            elif seq[2] not in ['A','G']:
                seq[2] = seq[2].lower()
            elif seq[4] not in ['C','T']:
                seq[4] = seq[4].lower()
            bp = BranchPoint(m.start(),m.end(),pattern="".join(seq))
            bp._gff['fscore'] = 0 
            sites.append(bp)
        sites.sort()

    # and return
    return sites

# end of function scan_regex_branchpoint


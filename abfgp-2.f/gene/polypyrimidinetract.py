"""
"""
# Module metadata
__authors__ = "Ate van der Burgt"
__license__ = "MIT"

# Python Imports
from re import finditer, compile, IGNORECASE
PPT_PATTERN = compile("[^AN]{5,}",IGNORECASE)

# Imports
from gene_gff import BasicGFF
from gene_exceptions import InproperlyAppliedArgument

# Import Global variables


class PolyPirimidineTract(BasicGFF):
    """ """
    def __init__(self,start,end,pattern="",gff={}):
        """ """
        BasicGFF.__init__(self)
        self._gff.update(gff)
        self.start      = start
        self.end        = end
        self.pos        = self.start
        self.pattern    = pattern
        self.length     = self.end - self.start
        if not self._gff.has_key('fscore'):
            self._gff['fscore'] = self.length
    # end of function __init__


    def __str__(self):
        """ print one-liner with PolyPirimidineTract data """
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


# end of class PolyPirimidineTract


def scan_regex_polypyrimidinetract(seq):
    """
    Find Polypyrimidine sequences by regular expression of the input sequence

    @type seq:  string
    @param seq: DNA sequence

    @rtype:  list
    @return: list with Polypyrimidine objects
    """
    base2score = { 'T':1.0, 'G':0.5, 'C':0.0, 'A': -2.0 }
    sites = []
    for m in finditer(PPT_PATTERN,seq):
        if m.group().upper().find("T") == -1:
            continue
        score = sum([ base2score[base]+1.0 for base in list(m.group().upper()) ])
        if score >= 8.0:
            # instantiate PolyPirimidineTract object
            ppt = PolyPirimidineTract(m.start(),m.end(),pattern=m.group(),gff={'fscore':score})
            sites.append(ppt)

    # and return
    return sites

# end of function scan_regex_polypyrimidinetract
 

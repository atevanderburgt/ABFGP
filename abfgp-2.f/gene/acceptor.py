"""
"""
# Module metadata
__authors__ = "Ate van der Burgt"
__license__ = "MIT"

# Imports
from splicesite import SpliceSiteBase, _projected_site_total_score

# Import Global variables
from settings.splicesites import *


class SpliceAcceptorAG(SpliceSiteBase):
    def __init__(self,start,phase=None,strand=None,gff={}):
        # initialize object from parental SliceSiteBase class
        SpliceSiteBase.__init__(self,start,phase=phase,strand=strand,pattern=None,
            pattern_offset=(0,0),pssm_score=None,gff=gff)
        self.end        = self.start+2
        self.pos        = self.end
        self.acceptor   = "AG"

    # end of function __init__

# end of class SpliceAcceptorAG


class SpliceAcceptor(SpliceSiteBase):
    def __init__(self,start,pattern,acceptor="AG",phase=None,strand=None,
        pattern_offset=IC_ACCEPTOR_PATTERN_OFFSET,pssm_score=None,gff={}):
        # initialize object from parental SliceSiteBase class
        SpliceSiteBase.__init__(self,start,phase=phase,strand=strand,gff=gff,
            pattern=pattern,pattern_offset=pattern_offset,pssm_score=pssm_score)
        self.end        = self.start + 2 + self._offset_5p + self._offset_3p
        self.pos        = self.end - self._offset_3p
        self.acceptor   = acceptor

    # end of function __init__

    def is_canonical(self):
        """ is splicesite canonical or not """
        if self.canonical_acceptor == self.acceptor:
            return True
        else:
            return False

    # end of function is_canonical

    def is_nagnag_acceptor(self):
        """ is this site a NAGNAG acceptor? """
        if self.pattern:
            if self.pattern[self._offset_5p-3:self._offset_5p-1].upper() == "AG":
                return True
            elif self.pattern[-self._offset_3p+1:-self._offset_3p+3].upper() == "AG":
                return True
            else:
                return False
        else:
            return None

    # end of function is_nagnag_acceptor 

    def togff(self,gff={}):
        """
        Overrides SpliceSiteBase.togff()
        Takes care for proper fstart and fstop setting
        """
        self._gff['fstart']  = self.pos-1
        self._gff['fstop']   = self.pos
        return SpliceSiteBase.togff(self,gff=gff)

    # end of function togff


# end of class SpliceAcceptor


class ProjectedSpliceAcceptor(SpliceAcceptor):
    def __init__(self,start,pattern,spliceacceptor,distance=None,
        entropy=None,pssm_score=None,cases=1):
        # initialize object from parental SpliceAcceptor class
        SpliceAcceptor.__init__(self,start,phase=spliceacceptor.phase,
            strand=spliceacceptor.strand, gff=spliceacceptor._gff,
            pattern=pattern,pattern_offset=(0,0) )
        self.end        = self.start + 2 + self._offset_5p + self._offset_3p
        self.pos        = self.start + self._offset_5p
        self.acceptor   = pattern.upper()
        self.distance   = distance
        self.cases      = cases
        # needed for compatibilty with SpliceSiteBase object
        self.pssm_score = pssm_score
        self.entropy    = entropy
        # special thingies for Projected flavour
        self.total_entropy      = entropy
        self.total_pssm_score   = pssm_score
        self.average_entropy    = entropy / float(cases)
        self.average_pssm_score = pssm_score / float(cases)

    # end of function __init__
    
    def is_canonical(self):
        """ is splicesite canonical or not """
        return False

    # end of function is_canonical

    def total_score(self):
        """ total score of this pojected splice site"""
        return _projected_site_total_score(self)

    # end of function total_score

    def __str__(self):
        """
        Overrides SpliceAcceptor.__str__()
        """
        retstr = SpliceAcceptor.__str__(self)
        retstr = retstr.replace(", score=", " [%sx], b.e.t.w. %1.2f, score=" % (self.cases,self.total_entropy) )
        return retstr

    # end of function __str__

    def togff(self,gff={}):
        """
        Overrides SpliceSiteBase.gff()
        Takes care for proper fstart and fstop setting
        """
        self._gff['column9data'] = {
                'Entropy'       : self.total_entropy,
                'Alignedsites'  : self.cases,
                'Distance'      : self.distance,
            }
        return SpliceAcceptor.togff(self,gff=gff)

    # end of function togff

# end of class ProjectedSpliceAcceptor


def find_AG_acceptor_sites(seq):
    """
    Find canonical AG acceptor sites on input sequence

    @type seq:  string
    @param seq: DNA sequence

    @rtype:  list
    @return: list with AG =end= coordinates (0-based Python string coords)
    """
    sites = []
    for item in finditer( "AG", seq.upper() ):
        sites.append( item.start()+2 )
    return sites

# end of function find_AG_acceptor_sites

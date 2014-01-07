"""
"""
# Module metadata
__authors__ = "Ate van der Burgt"
__license__ = "MIT"

# Imports
from splicesite import SpliceSiteBase, _projected_site_total_score

# Import Global variables
from settings.splicesites import *

class SpliceDonorGT(SpliceSiteBase):
    def __init__(self,start,phase=None,strand=None,gff={}):
        # initialize object from parental SliceSiteBase class
        SpliceSiteBase.__init__(self,start,phase=phase,strand=strand,pattern=None,
            pattern_offset=(0,0),pssm_score=None,gff=gff)
        self.end        = self.start+2
        self.pos        = self.start
        self.donor      = "GT"

    # end of function __init__

# end of class SpliceDonorGT


class SpliceDonor(SpliceSiteBase):
    def __init__(self,start,pattern,donor="GT",phase=None,strand=None,
        pattern_offset=IC_DONOR_PATTERN_OFFSET,pssm_score=None,gff={}):
        # initialize object from parental SliceSiteBase class
        SpliceSiteBase.__init__(self,start,phase=phase,strand=strand,gff=gff,
            pattern=pattern,pattern_offset=pattern_offset,pssm_score=pssm_score)
        self.end        = self.start + 2 + self._offset_5p + self._offset_3p
        self.pos        = self.start + self._offset_5p
        self.donor      = donor.upper()

    # end of function __init__

    def is_canonical(self):
        """ is splicesite canonical or not """
        if self.canonical_donor == self.donor:
            return True
        else:
            return False

    # end of function is_canonical

    def togff(self,gff={}):
        """
        Overrides SpliceSiteBase.togff()
        Takes care for proper fstart and fstop setting
        """
        self._gff['fstart']  = self.pos+1
        self._gff['fstop']   = self.pos+2
        if not self.is_canonical():
            if not self._gff.has_key('column9data'):
                self._gff['column9data'] = {}
            self._gff['column9data']['Non-canonical'] = self.donor
        # return basal pliceSiteBase.togff() function
        return SpliceSiteBase.togff(self,gff=gff)

    # end of function togff

# end of class SpliceDonor


class ProjectedSpliceDonor(SpliceDonor):
    def __init__(self,start,pattern,splicedonor,distance=None,
        entropy=None,pssm_score=None,cases=1):
        # initialize object from parental SpliceDonor class
        SpliceDonor.__init__(self,start,phase=splicedonor.phase,
            strand=splicedonor.strand, gff=splicedonor._gff,
            pattern=pattern,pattern_offset=(0,0) )
        self.end        = self.start + 2 + self._offset_5p + self._offset_3p
        self.pos        = self.start + self._offset_5p
        self.donor      = pattern.upper()
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
        Overrides SpliceDonor.__str__()
        """
        retstr = SpliceDonor.__str__(self)
        retstr = retstr.replace(", score=", " [%sx], b.e.t.w. %1.2f, score=" % (self.cases,self.total_entropy) )
        return retstr

    # end of function __str__

    def togff(self,gff={}):
        """
        Overrides SpliceSiteBase.togff()
        Takes care for proper fstart and fstop setting
        """
        self._gff['column9data'] = {
                'Entropy'       : self.total_entropy,
                'Alignedsites'  : self.cases,
                'Distance'      : self.distance,
            }
        return SpliceDonor.togff(self,gff=gff)

    # end of function togff

    
# end of class ProjectedSpliceDonor


def find_GT_donor_sites(seq):
    """
    Find canonical GT donor sites on input sequence

    @type seq:  string
    @param seq: DNA sequence

    @rtype:  list
    @return: list with GT =start= coordinates (0-based Python string coords)
    """
    sites = []
    for item in finditer( "GT", seq.upper().replace("U","T") ):
        sites.append( item.start() )
    sites.reverse()
    return sites

# end of function find_GT_donor_sites

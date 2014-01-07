"""
"""
# Module metadata
__authors__ = "Ate van der Burgt"
__license__ = "MIT"

# Imports
from gene_gff import BasicGFF
from validators import *
from gene_exceptions import (
    WrongPatternLengthApplied,
    InproperlyAppliedArgument,
    UnexpectedSpliceSitePhase,
    IncompatibleSpliceSitePhases,
    )
from pssm import parse_ic_file, Pssm, pssmscore


# Import Global variables
from settings.splicesites import (
    IC_DONOR_PATTERN_OFFSET,
    IC_DONOR_DATA_FILE,
    IC_DONOR_NCGC_DATA_FILE,
    IC_ACCEPTOR_PATTERN_OFFSET,
    IC_ACCEPTOR_DATA_FILE,
    )

# parse IC PSSM files of cannonical sites
IC_ACCEPTOR    = parse_ic_file(IC_ACCEPTOR_DATA_FILE)
IC_DONOR       = parse_ic_file(IC_DONOR_DATA_FILE)
IC_NC_GC_DONOR = parse_ic_file(IC_DONOR_NCGC_DATA_FILE)

class SpliceSiteBase(BasicGFF):
    """ """
    def __init__(self,start,phase=None,strand='+',pattern=None,
        pattern_offset=(0,0),pssm_score=None,gff={}):
        """
        Initialization function of Basal SpliceSite logic
        Recommended is to use only one of the inheriting classes

        @type  start: number
    	@param start: start coord of site (e.g GT) or pattern (e.g. tgtGTcgat)

        @type  phase: number
    	@param phase: [0,1,2] or None

        @type  phase: number
    	@param phase: [0,1,2] or None, default None

        @type  strand: string
    	@param strand: ['+','-'] or None, default '+'

        @type  pattern: string
    	@param pattern: splice site pattern sequence (e.g. tgtGTcgat)

        @type  pattern_offset: tuple
    	@param pattern_offset: tuple with 2 integers, specifying 5p and 3p
                               offset towards actual splice site.
                               e.g. (3,4) for tgtGTcgat

        @type  pssm_score: float
    	@param pssm_score: PSSM score of splice site pattern

        @type  gff: dictionary
    	@param gff: dictionary with gff data. Possible keys are:
                    'fref', 'fmethod', 'fstrand', 'fsource', 'gclass', 'gname'

        """
        BasicGFF.__init__(self)
        self.start      = start
        self.strand     = strand
        self.phase      = phase
        self.end        = None  # to be filled in inheriting classes
        self.pos        = None  # to be filled in inheriting classes
        self.canonical_donor = "GT"
        self.canonical_acceptor = "AG"
        # pattern, pattern_offset and pssm_score are only truely
        # functional in inheriting classes SpliceDonor and SpliceAcceptor
        self.pattern    = pattern
        self.pssm_score = pssm_score
        self._offset_5p = pattern_offset[0]
        self._offset_3p = pattern_offset[1]
        # and some GFF data thingies
        self._gff['fscore'] = self.pssm_score
        self._gff.update(gff)
        # and error-check the data
        self.check()

    # end of function __init__

    def check(self):
        """
        Basic error check of object's validity
        """
        IsProperIntronPhase(self.phase)
        if self.pattern:
            if len(self.pattern) != self._offset_5p + self._offset_3p + 2:
                raise WrongPatternLengthApplied
        else:
            if ( self._offset_5p, self._offset_3p ) != (0,0):
                raise WrongPatternLengthApplied

    # end of function check

    def togff(self,gff={}):
        """ Overwrites BasicGFF.togff() """
        if not gff.has_key('gname'): gff['gname'] = str(self.pos)
        return BasicGFF.togff(self,gff=gff)

    # end of function togff

    def __str__(self):
        """ print one-liner with SpliceSite data """
        checkA = self.pssm_score != None
        checkB = self.pattern != None
        _pattern = ""
        _score   = ""
        if checkA:
            _pattern = ", %s%s%s" % (
                self.pattern[0:self._offset_5p].lower(),
                self.pattern[self._offset_5p:-self._offset_3p].upper(),
                self.pattern[-self._offset_3p:].lower(),
                )       
        if checkB:
            _score = ", score=%2.3f" % self.pssm_score
        if self.pssm_score != None:
            return "<%s %s-%s (%s)%s%s>" % (
                self.__class__.__name__,
                self.start, self.end, self.phase, 
                _score, _pattern
                )
        else:
            return "<%s %s-%s (%s)>" % (
                self.__class__.__name__,
                self.start, self.end, self.phase, 
                )

    # end of function __str__


    def is_canonical(self):
        """ is splicesite canonical or not """
        return True

    # end of function is_canonical

# end of class SpliceSiteBase


def _projected_site_total_score(self):
    """
    Function in use as total_score() in 
    ProjectedSpliceDonor
    ProjectedSpliceAcceptor
    The total score of a projected site is a
    balanced average of PSSM, entropy, distance and
    number of cases this projection is supported by
    """
    case_correction = max([ 0.75*float(self.cases) , 1.0 ])
    pssm_part    = self.total_pssm_score / case_correction
    entropy_part = self.total_entropy / case_correction

    weight_pssm_part    = 1.0  # a single pssm score is in range -3.0..9.0
    weight_entropy_part = 1.0  # a single entropy score has range 0.0..1.0
                               # ATTENTION!! this lower range is corrected
                               # for in the ABFGP code at a later stage
                               # so, leave here as 1.0

    return ( ( pssm_part * weight_pssm_part ) +\
             ( entropy_part * weight_entropy_part ) ) /\
             float( abs(self.distance)/3 + 1 )

# end of function _projected_site_total_score


def _score_splice_site(seq,splicetype="donor"):
    """
    Return PSSM splice site score

    @type seq:  string
    @param seq: DNA sequence of EXACT length of the PSSM

    @type splicetype:   string
    @param splicetype:  'donor' or 'acceptor'

    @rtype:     float
    @return:    PSSM score for this splice site
    """
    if splicetype=='acceptor':
        return pssmscore(seq,IC_ACCEPTOR)
    else:
        if seq[IC_DONOR_PATTERN_OFFSET[0]:-IC_DONOR_PATTERN_OFFSET[-1]].upper() == "GC":
            return pssmscore(seq,IC_NC_GC_DONOR)
        else:
            return pssmscore(seq,IC_DONOR)

# end of _score_splice_site


class DonorScorer(Pssm):
    def __init__(self,ic=[],ignore_unambiguity=True,relativescore=True):
        """ """
        if not ic: ic=IC_DONOR
        Pssm.__init__(self,ic=ic,ignore_unambiguity=ignore_unambiguity,relativescore=relativescore)
    # end of function __init__
# end of class DonorScorer


class AcceptorScorer(Pssm):
    def __init__(self,ic=[],ignore_unambiguity=True,relativescore=True):
        """ """
        if not ic: ic=IC_ACCEPTOR
        Pssm.__init__(self,ic=ic,ignore_unambiguity=ignore_unambiguity,relativescore=relativescore)
    # end of function __init__
# end of class AcceptorScorer


def scan_pssm_donor(seq,**kwargs):
    """
    @attention: see scan_pssm_splice_site for **kwargs
    """
    return scan_pssm_splice_site(seq,splicetype="donor",**kwargs)

# end of function scan_pssm_donor


def scan_pssm_acceptor(seq,**kwargs):
    """
    @attention: see scan_pssm_splice_site for **kwargs
    """
    return scan_pssm_splice_site(seq,splicetype="acceptor",**kwargs)

# end of function scan_pssm_donor


def scan_pssm_splice_site(seq,splicetype="donor",
    override_pattern_offset=(),min_pssm_score=None,
    allow_non_canonical=False,non_canonical_min_pssm_score=0.0,
    ignore_unambiguity=False,relativescore=False,):
    """
    Find splice sites by a PSSM on input sequence

    @type  seq:  string
    @param seq: DNA sequence of EXACT length of the PSSM

    @type  splicetype:   string
    @param splicetype:  'donor' or 'acceptor'

    @type  min_pssm_score:   float
    @param min_pssm_score:

    @type  allow_non_canonical:  boolean
    @param allow_non_canonical: True of False

    @type  non_canonical_min_pssm_score:   float
    @param non_canonical_min_pssm_score:

    @type  override_pattern_offset:  tuple
    @param override_pattern_offset: tuple with 2 integers; use cautiously!!

    @rtype:  list
    @return: list with SpliceDonors or SpliceAcceptors
    """
    if splicetype == 'acceptor':
        PSSM_MATRIX     = IC_ACCEPTOR
        pattern_offset  = IC_ACCEPTOR_PATTERN_OFFSET 
        canonical       = "AG"
        # initialize Psmm (scoring) class
        canssPssm = Pssm(ic=IC_ACCEPTOR,ignore_unambiguity=ignore_unambiguity,relativescore=relativescore)
        # import output SpliceAcceptor object
        from acceptor import SpliceAcceptor
    elif splicetype == 'donor':
        PSSM_MATRIX     = IC_DONOR
        pattern_offset  = IC_DONOR_PATTERN_OFFSET 
        canonical       = "GT"
        # initialize Psmm (scoring) class
        canssPssm = Pssm(ic=IC_DONOR,ignore_unambiguity=ignore_unambiguity,relativescore=relativescore)
        # import output SpliceDonor object
        from donor import SpliceDonor
    else:
        message = "'splicetype' (%s) not in [donor,acceptor]" % splicetype
        raise InproperlyAppliedArgument, message

    if allow_non_canonical:
        # obtain PSSM_IC for non-canonical (GC) donors
        IC_NCGC_DONOR = parse_ic_file(IC_DONOR_NCGC_DATA_FILE)
        noncanonical = ["GC"]
        # initialize Psmm (scoring) class
        noncanssPssm = Pssm(ic=IC_NCGC_DONOR,ignore_unambiguity=ignore_unambiguity,relativescore=relativescore)


    # hmm... somebody knows what he or she is doing ;-)
    if override_pattern_offset:
        pattern_offset = override_pattern_offset

    pssmlength = len(PSSM_MATRIX)
    sites = []
    for offset in range(0, len(seq) - pssmlength + 1 ):
        # get sequence slice of pattern and actual splice site
        seqpart = seq[offset:offset+pssmlength].upper()
        splicesite = seqpart[pattern_offset[0]:-pattern_offset[1]]

        # continue if non-canonical sites if not requested for
        if not allow_non_canonical and splicesite != canonical:
            continue
        elif splicesite == canonical:
            # score this splicesite
            #score = _score_splice_site(seqpart,splicetype=splicetype)
            score = canssPssm.score(seqpart)
            # check if site must be stored
            if min_pssm_score or min_pssm_score == 0.0:
                if score < min_pssm_score:
                    continue
        elif splicesite != canonical and splicetype == 'donor' and splicesite in noncanonical:
            # score non-canonical donor site
            #score = pssmscore(seqpart,IC_NCGC_DONOR) 
            score = noncanssPssm.score(seqpart)
            # check if site must be stored
            if non_canonical_min_pssm_score or non_canonical_min_pssm_score == 0.0:
                if score < non_canonical_min_pssm_score:
                    continue
            ####print seqpart, score, offset
        else:
            continue

        if splicetype=='acceptor':
            a = SpliceAcceptor(offset,seqpart,acceptor=splicesite,pssm_score=score)
            sites.append(a)
        else:
            d = SpliceDonor(offset,seqpart,donor=splicesite,pssm_score=score)
            sites.append(d)

    # return sites for Donor
    if splicetype == 'donor':
        sites.reverse()

    # and return
    return sites

# end of function scan_pssm_splice_site


def scan_orf_for_pssm_splice_sites(orf,splicetype="donor",min_pssm_score=None,
    allow_non_canonical=False,non_canonical_min_pssm_score=0.0):
    """
    """
    if splicetype=='acceptor':
        pattern_offset  = IC_ACCEPTOR_PATTERN_OFFSET 
        offset_5p       = pattern_offset[0]+4 
        offset_3p       = 0
    else:
        pattern_offset  = IC_DONOR_PATTERN_OFFSET 
        offset_5p       = 0
        offset_3p       = pattern_offset[1]+4 

    # get (elongated) sequence of this orf
    seqslice = orf.nucleotidesequence(extra_start=offset_5p,extra_end=offset_3p)

    # scan for occurrences
    sites = scan_pssm_splice_site(seqslice,splicetype=splicetype,
        min_pssm_score=min_pssm_score,allow_non_canonical=allow_non_canonical,
        non_canonical_min_pssm_score=non_canonical_min_pssm_score)

    # correct site positions to absolute coords
    # and set correct phase of the splice site
    for site in sites:
        site.start = site.start + orf.startPY - offset_5p
        site.end   = site.end   + orf.startPY - offset_5p
        site.pos   = site.pos   + orf.startPY - offset_5p
        site.phase = (site.pos - orf.startPY) % 3

    # and return the splice sites on this orf
    return sites

# end of function scan_orf_for_pssm_splice_sites


def get_shared_nucleotides_at_splicesite(orfA,orfD,acceptor,donor):
    """
    """
    # do splice site phase compatibility check
    if acceptor.phase != donor.phase:
        raise IncompatibleSpliceSitePhases

    if donor.phase == 0:
        shared_nts = ""
    elif donor.phase == 1:
        shared_nts = "%s%s" % (
            orfD.inputgenomicsequence[donor.pos-1:donor.pos],
            orfA.inputgenomicsequence[acceptor.pos:acceptor.pos+2]
            )
    elif donor.phase == 2:
        shared_nts = "%s%s" % (
            orfD.inputgenomicsequence[donor.pos-2:donor.pos],
            orfA.inputgenomicsequence[acceptor.pos:acceptor.pos+1]
            )
    else:
        # what else !? Unexpected error!
        raise UnexpectedSpliceSitePhase

    # return shared_nucleotides
    return shared_nts

# end of function get_shared_nucleotides_at_splicesite

"""
"""
# Module metadata
__authors__ = "Ate van der Burgt"
__license__ = "MIT"

# Imports
from gene_gff import BasicGFF
from gene_exceptions import InproperlyAppliedArgument
from validators import IsProperPhase
from pssm import parse_ic_file, Pssm, pssmscore

# Python Imports
from re import finditer, compile
from copy import deepcopy

# Import Global variables
from settings.translationalstartsites import (
    TSS_MIN_PSSM_SCORE,
    TSS_ALLOW_NON_CANONICAL,
    TSS_NON_CANONICAL_MIN_PSSM_SCORE,
    IC_TSS_DATA_FILE,
    IC_TSS_PATTERN_OFFSET,
    )

# parse IC PSSM file of TSS
IC_TSS = parse_ic_file(IC_TSS_DATA_FILE)



class StartCodon(BasicGFF):
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

# end of class StartCodon


def _handle_startcodon(start):
    """ """
    start_site_classes = ['TranslationalStartSiteATG','TranslationalStartSite']
    if start.__class__.__name__ in start_site_classes:
        startcodon = start
    elif start.__class__.__name__ in ['StartCodon']:
        startcodon = start
    elif type(start) == type(int()):
        startcodon = StartCodon(start)
    else:
        message = "``start`` not an expected object or integer"
        raise InproperlyAppliedArgument, message
    # return start object
    return startcodon

# end of function _handle_startcodon


class TranslationalStartSiteBase(BasicGFF):
    """
    Class describing the position of a putative TranslationalStartSite (ATG)
    """
    def __init__(self,start,phase='.',strand='+',pattern=None,
        pattern_offset=(0,0),pssm_score=None,gff={}):
        """
        Initialization function of Basal TranslationalStartSiteBase logic
        Recommended is to use only one of the inheriting classes

        @type  start: number
    	@param start: start coord of site (ATG) or pattern (e.g. nncancATGgcnnn)

        @type  phase: number
    	@param phase: [0,1,2] or None; only relevant when sequence ISA orf

        @type  strand: string
    	@param strand: ['+','-'] or None, default '+'

        @type  pattern: string
    	@param pattern: splice site pattern sequence (e.g. nncancATGgcnnn)

        @type  pattern_offset: tuple
    	@param pattern_offset: tuple with 2 integers, specifying 5p and 3p
                               offset towards actual splice site.
                               e.g. (0,0) for the base class (no PSSM)
                               e.g. (6,5) for nncancATGgcnnn
                               e.g. (9,7) for another implementation

        @type  pssm_score: float
    	@param pssm_score: PSSM score of splice site pattern

        @type  gff: dictionary
    	@param gff: dictionary with gff data. Possible keys are:
                    'fref', 'fmethod', 'fstrand', 'fsource', 'gclass', 'gname'

        """
        # initialize basal BasicGff
        BasicGFF.__init__(self)
        self.start      = start
        self.strand     = strand
        self.phase      = phase
        self.end        = None  # to be filled in inheriting classes
        self.pos        = None  # to be filled in inheriting classes
        self.tss        = None  # to be filled in inheriting classes
        self.canonical_tss = "ATG"
        # pattern, pattern_offset and pssm_score are only truely
        # functional in inheriting class TranslationalStartSite,
        # not in TranslationalStartSiteATG
        self.pattern    = pattern
        self.pssm_score = pssm_score
        self._offset_5p = pattern_offset[0]
        self._offset_3p = pattern_offset[1]
        # update gff dict
        self._gff.update(gff)
        # and error-check the data
        IsProperPhase(self.phase)
        self.check()

    # end of function __init__

    def check(self):
        """
        Basic error check of object's validity
        """
        if self.pattern:
            if len(self.pattern) != self._offset_5p + self._offset_3p + 3:
                raise WrongPatternLengthApplied
        else:
            if ( self._offset_5p, self._offset_3p ) != (0,0):
                raise WrongPatternLengthApplied

    # end of function check

    def __str__(self):
        """ print one-liner with TranslationStartSite data """
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
            _score = ", score=%2.2f" % self.pssm_score
        return "<%s %s (%s)%s%s (%s-%s)>" % (
                self.__class__.__name__,
                self.pos, self.phase, 
                _score, _pattern, self.start, self.end
                )

    # end of function __str__


    def is_canonical(self):
        """ is TSS canonical or not """
        return True

    # end of function is_canonical

# end of class TranslationalStartSiteBase


class TranslationalStartSiteATG(TranslationalStartSiteBase):
    """ class describing a(ny) ATG as a (putative) TranslationalStartSite """
    def __init__(self,start,phase='.',strand='+',gff={}):
        # initialize object from parental TranslationalStartSiteBase class
        TranslationalStartSiteBase.__init__(self,start,phase=phase,strand=strand,pattern=None,
            pattern_offset=(0,0),pssm_score=None,gff=gff)
        self.end        = self.start+3
        self.pos        = self.start
        self.tss        = "ATG"

    # end of function __init__

# end of class TranslationalStartSiteATG


class TranslationalStartSite(TranslationalStartSiteBase):
    """ class describing a pssm around a(ny) ATG as a (putative) TranslationalStartSite """
    def __init__(self,start,pattern,tss="ATG",phase='.',strand='+',
        pattern_offset=IC_TSS_PATTERN_OFFSET,pssm_score=None,gff={}):
        # initialize object from parental TranslationalStartSiteBase class
        TranslationalStartSiteBase.__init__(self,start,phase=phase,strand=strand,gff=gff,
            pattern=pattern,pattern_offset=pattern_offset,pssm_score=pssm_score)
        # +2 and not +3 for length(ATG) because of GFF/Python off-by-one coordinate
        self.end        = self.start + 2 + self._offset_5p + self._offset_3p
        self.pos        = self.start + self._offset_5p
        self.tss        = tss.upper()

    # end of function __init__

    def is_canonical(self):
        """ is TranslationalStartSite canonical or not """
        if self.canonical_tss == self.tss:
            return True
        else:
            return False

    # end of function is_canonical

    def togff(self,gff={}):
        """
        Overrides TranslationalStartSiteBase.togff()
        Takes care for proper fstart and fstop setting
        """
        self._gff['fstart']  = self.pos+1
        self._gff['fstop']   = self.pos+3
        self._gff['gname']   = "TSS%s(%s)" % (self.pos+1,self.phase)
        return TranslationalStartSiteBase.togff(self,gff=gff)

    # end of function togff

# end of class TranslationalStartSite


def find_ATG_tss(seq):
    """
    Find canonical ATG TSS's on input sequence

    @type seq:  string
    @param seq: DNA sequence

    @rtype:  list
    @return: list with ATG =start= coordinates (0-based Python string coords)
    """
    sites = []
    for item in finditer( "ATG", seq.upper().replace("U","T") ):
        sites.append( item.start() )
    sites.reverse()
    return sites

# end of function find_ATG_tss


def score_tss(seq,ignore_unambiguity=False):
    """
    Return PSSM TSS score

    @type seq:  string
    @param seq: DNA sequence of EXACT length of the PSSM

    @rtype:     float
    @return:    PSSM score for this TSS (translational start site)
    """
    return pssmscore(seq,IC_TSS,ignore_unambiguity=ignore_unambiguity)

# end of score_tss


def scan_pssm_tss(seq,override_pattern_offset=(),min_pssm_score=None,
    allow_non_canonical=False,non_canonical_min_pssm_score=0.0,
    ignore_unambiguity=False,relativescore=False,):
    """
    Find TSS's by a PSSM on input sequence

    @type seq:  string
    @param seq: DNA sequence of EXACT length of the PSSM

    @type min_pssm_score:   float
    @param min_pssm_score:

    @type allow_non_canonical:  boolean
    @param allow_non_canonical: True of False

    @type non_canonical_min_pssm_score:   float
    @param non_canonical_min_pssm_score:

    @type override_pattern_offset:  tuple
    @param override_pattern_offset: tuple with 2 integers; use cautiously!!

    @rtype:  list
    @return: list with TranslationalStartSites
    """
    pattern_offset  = IC_TSS_PATTERN_OFFSET
    canonical       = "ATG"
    regex_pattern   = "ATG"
    # obtain Psmm (scoring) class
    tssPssm = Pssm(ic=IC_TSS,ignore_unambiguity=ignore_unambiguity,relativescore=relativescore)

    # hmm... somebody knows what he or she is doing ;-)
    if override_pattern_offset:
        pattern_offset = override_pattern_offset

    if allow_non_canonical:
        # http://en.wikipedia.org/wiki/Start_codon
        # Blattner, F. R.; Plunkett, G.; Bloch, C. A.; Perna, N. T.; Burland, V.; Riley, M.; Collado-vides, J.; Glasner, J. D. et al (1997). "The Complete Genome Sequence of Escherichia coli K-12". Science 277 (5331): 1453-62
        noncanonicalGTG = "GTG"
        noncanonicalTTG = "TTG"
        regex_pattern   = "[ATG]TG"
        tssPssmGTG = deepcopy(tssPssm)
        tssPssmGTG.ic[pattern_offset[0]]['A'],tssPssmGTG.ic[pattern_offset[0]]['G'] = tssPssmGTG.ic[pattern_offset[0]]['G'],tssPssmGTG.ic[pattern_offset[0]]['A'] 
        tssPssmTTG = deepcopy(tssPssm)
        tssPssmTTG.ic[pattern_offset[0]]['A'],tssPssmGTG.ic[pattern_offset[0]]['T'] = tssPssmTTG.ic[pattern_offset[0]]['T'],tssPssmTTG.ic[pattern_offset[0]]['A']  

    pssmlength = len(IC_TSS)
    sites = []
    
    for atg in finditer(compile(regex_pattern),seq.upper()[0:-pattern_offset[1]]):
        if atg.start() < pattern_offset[0]: continue
        seqpart = seq[atg.start()-pattern_offset[0]:atg.end()+pattern_offset[1]]
        tss = atg.group()
        offset = atg.start()-pattern_offset[0]

        # continue if non-canonical sites if not requested for
        if tss == canonical:
            # score this TSS
            score = tssPssm.score(seqpart)

            # continue if score is to low
            if min_pssm_score or min_pssm_score == 0.0:
                if score < min_pssm_score:
                    continue
        else:
            if atg.group() == noncanonicalGTG:
                score = tssPssmGTG.score(seqpart)
            elif atg.group() == noncanonicalTTG: 
                score = tssPssmTTG.score(seqpart)
            else:
                continue

            # continue if score is to low
            if non_canonical_min_pssm_score or non_canonical_min_pssm_score == 0.0:
                if score < non_canonical_min_pssm_score: 
                    continue


        # phase of this tss
        phase = (offset+pattern_offset[0]) % 3

        # instantiate tss object
        t = TranslationalStartSite(offset,seqpart,pssm_score=score,phase=phase)
        sites.append(t)

    # and return
    return sites

# end of function scan_pssm_tss

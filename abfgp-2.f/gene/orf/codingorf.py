"""
"""

# Module metadata
__authors__ = "Ate van der Burgt"
__license__ = "MIT"

# Python Imports
from re import finditer as refinditer

# Orf Imports
from basicorf import BasicOrf

# Orf Imports
from pssm import (
    scan_orf_for_pssm_splice_sites,
    scan_orf_for_pssm_tss,
    )


# Global variable Imports
from settings.genestructure import MIN_EXON_NT_LENGTH
from settings.splicesites import (
    MIN_DONOR_PSSM_SCORE,
    ALLOW_NON_CANONICAL_DONOR,
    NON_CANONICAL_MIN_DONOR_PSSM_SCORE,
    MIN_ACCEPTOR_PSSM_SCORE,
    ALLOW_NON_CANONICAL_ACCEPTOR,
    NON_CANONICAL_MIN_ACCEPTOR_PSSM_SCORE,
    )

DONOR_KWARGS = {
    'min_pssm_score':                MIN_DONOR_PSSM_SCORE,
    'allow_non_canonical':           ALLOW_NON_CANONICAL_DONOR,
    'non_canonical_min_pssm_score':  NON_CANONICAL_MIN_DONOR_PSSM_SCORE,
    }

ACCEPTOR_KWARGS = {
    'min_pssm_score':                MIN_ACCEPTOR_PSSM_SCORE,
    'allow_non_canonical':           ALLOW_NON_CANONICAL_ACCEPTOR,
    'non_canonical_min_pssm_score':  NON_CANONICAL_MIN_ACCEPTOR_PSSM_SCORE,
}

def _update_kwargs(kwargs,defaults):
    """
    Update **kwargs dict with defaults for non-existing keys
    """
    for k,v in defaults.iteritems():
        if not kwargs.has_key(k):
            kwargs[k] = v
# end of function_update_kwargs


class CodingOrf(BasicOrf):
    """
    """
    def __init__(self,*args,**kwargs):
        """
        """
        BasicOrf.__init__(self,*args,**kwargs)

        # data attributes for splice sites and TSS on this orf
        self._donor_sites    = []
        self._acceptor_sites = []
        self._has_donor_sites_predicted    = False
        self._has_acceptor_sites_predicted = False
        self._projected_donor_sites    = []
        self._projected_acceptor_sites = []
        self._tss_sites = []
        self._has_tss_predicted = False
        self._signalp_sites = []
        self._has_signalp_sites_predicted = False

        # attributes that hold the settings with which
        # pssm sites are searched for
        self.TSS_MIN_PSSM_SCORE      = None
        self.MIN_ACCEPTOR_PSSM_SCORE = None
        self.MIN_DONOR_PSSM_SCORE    = None

    # end of function __init__


    def scan_orf_for_pssm_tss(self,min_pssm_score=None,forced=False):
        """
        Scan Orf object for Translational Start Sites by a PSSM

        @type  min_pssm_score: float or None
    	@param min_pssm_score: minimal splice site pssm score

        @type  forced: Boolean
    	@param forced: if True, redo the scan even if _has_tss_predicted is True
        """
        if self._has_tss_predicted and not forced:
            # start sites are already predicted
            pass
        else:
            self._tss_sites = scan_orf_for_pssm_tss(
                self, min_pssm_score=min_pssm_score)
            # set _has_tss_predicted to True 
            self._has_tss_predicted = True
            self.TSS_MIN_PSSM_SCORE = min_pssm_score 

    # end of scan_orf_for_pssm_tss


    def scan_orf_for_pssm_donor_sites(self,forced=False,**kwargs):
	"""
	Scan Orf object for donor splice sites by a PSSM

        @type  forced: Boolean
    	@param forced: re-scan even if self._has_donor_sites_predicted = True

        @type  min_pssm_score: float or None
    	@param min_pssm_score: minimal splice site pssm score

        @type  allow_non_canonical: Boolean
    	@param allow_non_canonical: allow non-canonical (non GT) splice site

        @type  non_canonical_min_pssm_score: float or None
    	@param non_canonical_min_pssm_score: minimal non-canonical pssm score
	"""
	if self._has_donor_sites_predicted and not forced:
	    # splice sites are already predicted
	    pass
	else:
	    # set defaults for **kwargs when not applied
	    _update_kwargs(kwargs,DONOR_KWARGS)
	    # predict donor sites
	    self._donor_sites = scan_orf_for_pssm_splice_sites(
		    self,splicetype = 'donor',
		    min_pssm_score = kwargs['min_pssm_score'],
		    allow_non_canonical = kwargs['allow_non_canonical'],
		    non_canonical_min_pssm_score =
			kwargs['non_canonical_min_pssm_score']
		    )
	    # set _has_donor_sites_predicted to True 
	    self._has_donor_sites_predicted = True
	    self.MIN_DONOR_PSSM_SCORE = kwargs['min_pssm_score'] 

    # end of function scan_orf_for_pssm_donor_sites


    def scan_orf_for_pssm_acceptor_sites(self,forced=False,**kwargs):
	"""
	Scan Orf object for acceptor splice sites by a PSSM

        @type  forced: Boolean
    	@param forced: re-scan even if self._has_acceptor_sites_predicted = True

        @type  min_pssm_score: float or None
    	@param min_pssm_score: minimal splice site pssm score

        @type  allow_non_canonical: Boolean
    	@param allow_non_canonical: allow non-canonical (non AG) splice site

        @type  non_canonical_min_pssm_score: float or None
    	@param non_canonical_min_pssm_score: minimal non-canonical pssm score
	"""
	if self._has_acceptor_sites_predicted and not forced:
	    # splice sites are already predicted
	    pass
	else:
	    # set defaults for **kwargs when not applied
	    _update_kwargs(kwargs,ACCEPTOR_KWARGS)
	    # predict donor sites
	    self._acceptor_sites = scan_orf_for_pssm_splice_sites(
		    self,splicetype = 'acceptor',
		    min_pssm_score = kwargs['min_pssm_score'],
		    allow_non_canonical = kwargs['allow_non_canonical'],
		    non_canonical_min_pssm_score =
			kwargs['non_canonical_min_pssm_score']
		    )
            # set _has_acceptor_sites_predicted to True 
            self._has_acceptor_sites_predicted = True
            self.MIN_ACCEPTOR_PSSM_SCORE = kwargs['min_pssm_score']

    # end of function scan_orf_for_pssm_acceptor_sites


    def scan_orf_for_pssm_splice_sites(self,**kwargs):
        """
        Scan Orf object for donor & acceptor splice sites by a PSSM

	@attention: see scan_orf_for_pssm_donor_sites for **kwargs
	@attention: see scan_orf_for_pssm_acceptor_sites for **kwargs
        """
	self.scan_orf_for_pssm_donor_sites(**kwargs)
	self.scan_orf_for_pssm_acceptor_sites(**kwargs)

    # end of scan_orf_for_pssm_splice_sites


    def has_methionine(self):
        """
        Does Orf has a Methionine codon?

        @rtype:  Boolean
        @return: True or False
        """
        if self.protein_sequence.upper().find("M") == -1:
            return False
        else:
            return True

    # end of function has_methionine


    def has_start(self):
        """
        Does Orf has a Methionine codon?

	@attention: alias for has_methionine()

        @rtype:  Boolean
        @return: True or False
        """
	return self.has_methionine()

    # end of function has_start


    def has_translationalstartsite(self):
        """
        Does Orf has a TranslationalStartSite assigned?

        @rtype:  Boolean
        @return: True or False
        """
        if self._tss_sites: return True
        else:               return False

    # end of function has_translationalstartsite


    def has_tss_upstream_of(self,pos):
        """
        Does Orf has a TranslationalStartSite upstream of a certain position?

        @type  pos: positive integer
        @param pos: absolute AA coordinate

        @rtype:  Boolean
        @return: True or False
        """
        if not self._tss_sites: return False
        for site in self._tss_sites:
            if (site.pos / 3) <= pos:
                return True
        else:
            return False

    # end of function has_tss_upstream_of
 

    def has_start_upstream_of(self,pos):
        """
        Does Orf has a start (Methionine) codon upstream of a certain position?

        @type  pos: positive integer
    	@param pos: absolute AA coordinate

        @rtype:  Boolean
    	@return: True or False
        """
        offset = pos - self.protein_startPY
        if self.protein_sequence[0:offset].upper().find("M") == -1:
            return False
        else:
            return True

    # end of function has_start_upstream_of


    def potential_start_aa_positions(self,COORDINATES='ABSOLUTE'):
        """
        Return AA coordinates of start positions
        """
        positions = []
        for mm in refinditer("M",self.protein_sequence):
            positions.append( mm.start() + self.protein_startPY )
        return positions

    # end of function potential_start_aa_positions


    def potential_start_nt_positions(self,COORDINATES='ABSOLUTE'):
        """
        Return NT coordinates of start positions
        """
        return [ self.aapos2dnapos(pos) for pos in self.potential_start_aa_positions() ]

    # end of function potential_start_nt_positions


    def can_orf_encode_internalexon(self,min_exon_nt_length=MIN_EXON_NT_LENGTH,**kwargs):
	"""
	Can this OpenReadingFrame potentially encode an internal coding exon?

	@type  min_exon_nt_length: integer
	@param min_exon_nt_length: minimal exon length (nt)

	@rtype:  Boolean
	@return: True or False

	@attention: see _filter_exon_start_element for additional **kwargs
	@attention: see _filter_exon_end_element for additional **kwargs
	"""
	self.scan_orf_for_pssm_splice_sites(splicetype="acceptor")
	if not self._acceptor_sites:
	    return False
	self.scan_orf_for_pssm_splice_sites(splicetype="donor")
	if not self._donor_sites:
	    return False
	if self._donor_sites[0].pos - self._acceptor_sites[0].pos >= min_exon_nt_length:
	    acceptors = self._filter_exon_start_element(self._acceptor_sites,**kwargs)
	    if not acceptors:
		return False
	    donors = self._filter_exon_end_element(self._donor_sites,**kwargs)
	    if not donors:
		return False
	    if donors[0].pos - acceptors[0].pos >= min_exon_nt_length:
		return True
	    else:
		return False
	else:
	    return False

    # end of function can_orf_encode_internalexon

	
    def can_orf_encode_firstexon(self,min_exon_nt_length=MIN_EXON_NT_LENGTH,**kwargs):
	"""
	Can this OpenReadingFrame potentially encode a first coding exon?

	@type  min_exon_nt_length: integer
	@param min_exon_nt_length: minimal exon length (nt)

	@rtype:  Boolean
	@return: True or False

	@attention: see _filter_exon_start_element for additional **kwargs
	@attention: see _filter_exon_end_element for additional **kwargs
	"""
	self.scan_orf_for_pssm_tss()
	if not self._tss_sites:
	    return False
	self.scan_orf_for_pssm_splice_sites(splicetype="donor")
	if not self._donor_sites:
	    return False
	if self._donor_sites[0].pos - self._tss_sites[0].pos >= min_exon_nt_length:
	    tss = self._filter_exon_start_element(self._tss_sites,**kwargs)
	    if not tss:
		return False
	    donors = self._filter_exon_end_element(self._donor_sites,**kwargs)
	    if not donors:
		return False
	    if donors[0].pos - tss[0].pos >= min_exon_nt_length:
		return True
	    else:
		return False
	else:
	    return False

    # end of function can_orf_encode_firstexon


    def can_orf_encode_finalexon(self,min_exon_nt_length=MIN_EXON_NT_LENGTH,**kwargs):
	"""
	Can this OpenReadingFrame potentially encode a final coding exon?

	@type  min_exon_nt_length: integer
	@param min_exon_nt_length: minimal exon length (nt)

	@rtype:  Boolean
	@return: True or False

	@attention: see _filter_exon_start_element for additional **kwargs
	@attention: see _filter_exon_end_element for additional **kwargs
	"""
	self.scan_orf_for_pssm_splice_sites(splicetype="acceptor")
	if not self._acceptor_sites:
	    return False
	if self.endPY - self._acceptor_sites[0].pos >= min_exon_nt_length:
	    acceptors = self._filter_exon_start_element(self._acceptor_sites,**kwargs)
	    if not acceptors:
		return False
	    if self.endPY - acceptors[0].pos >= min_exon_nt_length:
		return True
	    else:
		return False
	else:
	    return False

    # end of function can_orf_encode_finalexon


    def can_orf_encode_singleexon(self,min_exon_nt_length=MIN_EXON_NT_LENGTH,**kwargs):
	"""
	Can this OpenReadingFrame potentially encode a single exon (gene)?

	@type  min_exon_nt_length: integer
	@param min_exon_nt_length: minimal exon length (nt)

	@rtype:  Boolean
	@return: True or False

	@attention: see _filter_exon_start_element for additional **kwargs
	@attention: see _filter_exon_end_element for additional **kwargs
	"""
	self.scan_orf_for_pssm_tss()
	if not self._tss_sites:
	    return False
	if self.endPY - self._tss_sites[0].pos >= min_exon_nt_length:
	    tss = self._filter_exon_start_element(self._tss_sites,**kwargs)
	    if not tss:
		return False
	    if self.endPY - tss[0].pos >= min_exon_nt_length:
		return True
	    else:
		return False
	else:
	    return False

    # end of function can_orf_encode_singleexon


    def _filter_exon_start_element(self,elementlist,
	exon_dna_start=None,exon_protein_start=None,**kwargs):
	"""
	"""
	if exon_dna_start:
	    offset = int(exon_dna_start)
	elif exon_protein_start:
	    offset = int(exon_protein_start)*3
	elif exon_dna_start and exon_protein_start:
	    offset = max( [ int(exon_protein_start)*3, int(exon_dna_start) ] )
	else:
	    return elementlist
	passedelements = []
	for elem in elementlist:
	    if elem.pos < offset:
		continue
	    else:
		passedelements.append(elem)
	# return list with passed elements (TSS,acceptors,etc)
	return passedelements
    
    # end of function _filter_exon_start_element


    def _filter_exon_end_element(self,elementlist,
	exon_dna_end=None,exon_protein_end=None,**kwargs):
	"""
	"""
	if exon_dna_end:
	    offset = int(exon_dna_end)
	elif exon_protein_end:
	    offset = int(exon_protein_end)*3
	elif exon_dna_end and exon_protein_end:
	    offset = min( [ int(exon_protein_end)*3, int(exon_dna_end) ] )
	else:
	    return elementlist
	passedelements = []
	for elem in elementlist:
	    if elem.pos > offset:
		continue
	    else:
		passedelements.append(elem)
	passedelements.reverse()
	# return list with passed elements (donors,etc)
	return passedelements
    
    # end of function _filter_exon_end_element

# end of class CodingOrf

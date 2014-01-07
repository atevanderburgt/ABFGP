"""


"""

# Module metadata
__authors__ = "Ate van der Burgt"
__license__ = "MIT"

# Imports
from lib_intron import *

# Global variables
try:
    # import global variables from settings file
    from settings.genestructure import *
except:
    # when using as stand-alone library (not in abfgp), specify these variables
    # length range of short leading exon and first intron
    SHORT_LEADINGEXON_MIN_NT_LENGTH                         = 10 
    SHORT_LEADINGEXON_MAX_NT_LENGTH                         = 300
    SHORT_LEADINGEXON_MAX_INTRON_NT_LENGTH                  = 300
    SHORT_LEADINGEXON_MIN_INTRON_NT_LENGTH                  = MIN_INTRON_NT_LENGTH
    # thresholds for TSS of short leading exon
    SHORT_LEADINGEXON_TSS_MIN_PSSM_SCORE                    = float(0)
    SHORT_LEADINGEXON_TSS_ALLOW_NON_CANONICAL               = False
    SHORT_LEADINGEXON_TSS_NON_CANONICAL_MIN_PSSM_SCORE      = float(0)
    # thresholds for splice site signals of short leading exon
    SHORT_LEADINGEXON_MIN_TOTAL_PSSM_SCORE                  = None
    SHORT_LEADINGEXON_MIN_DONOR_PSSM_SCORE                  = 0.0
    SHORT_LEADINGEXON_MIN_ACCEPTOR_PSSM_SCORE               = 0.0
    SHORT_LEADINGEXON_ALLOW_NON_CANONICAL_DONOR             = False     # not implemented yet
    SHORT_LEADINGEXON_ALLOW_NON_CANONICAL_ACCEPTOR          = False     # not implemented yet
    SHORT_LEADINGEXON_NON_CANONICAL_MIN_PSSM_SCORE          = None      # DEPRECATED
    SHORT_LEADINGEXON_NON_CANONICAL_MIN_DONOR_PSSM_SCORE    = float(0)  # not implemented yet
    SHORT_LEADINGEXON_NON_CANONICAL_MIN_ACCEPTOR_PSSM_SCORE = float(0)  # not implemented yet

# Logging function for debugging
def logf(*args):
    """ Logging function for debugging purposes. """
    x = ", ".join([str(item) for item in args])
    if x[0] == "#": print x
    else:           pass

# InproperlyAppliedArgument exception
class InproperlyAppliedArgument(Exception):
    """ """
    def __init__(self, message):
        self.message = message
    def __str__(self):
        return repr(self.message)



def find_leading_exon_on_orf(orfD,orfA,order_by='total_pssm',
    max_leadingexon_nt_length=SHORT_LEADINGEXON_MAX_NT_LENGTH,
    min_leadingexon_nt_length=SHORT_LEADINGEXON_MIN_NT_LENGTH,
    max_leadingexon_intron_nt_length=SHORT_LEADINGEXON_MAX_INTRON_NT_LENGTH,
    min_leadingexon_intron_nt_length=SHORT_LEADINGEXON_MIN_INTRON_NT_LENGTH,
    min_total_pssm_score=SHORT_LEADINGEXON_MIN_TOTAL_PSSM_SCORE,
    min_donor_pssm_score=SHORT_LEADINGEXON_MIN_DONOR_PSSM_SCORE,
    min_acceptor_pssm_score=SHORT_LEADINGEXON_MIN_ACCEPTOR_PSSM_SCORE,
    min_tss_pssm_score=SHORT_LEADINGEXON_TSS_MIN_PSSM_SCORE,
    allow_non_canonical_donor=SHORT_LEADINGEXON_ALLOW_NON_CANONICAL_DONOR,
    allow_non_canonical_acceptor=SHORT_LEADINGEXON_ALLOW_NON_CANONICAL_ACCEPTOR,
    non_canonical_min_donor_pssm_score=SHORT_LEADINGEXON_NON_CANONICAL_MIN_DONOR_PSSM_SCORE,
    non_canonical_min_acceptor_pssm_score=SHORT_LEADINGEXON_NON_CANONICAL_MIN_ACCEPTOR_PSSM_SCORE,
    subsequent_acceptor=None,
    subsequent_acceptor_pos=None,
    verbose=False):
    """

    @type  orfD: Orf object
	@param orfD: Orf object to scan for a leadingexon

    @type  orfA: Orf object
	@param orfA: Orf object from where to start the search

    @type  subsequent_acceptor: object
	@param subsequent_acceptor: SpliceAcceptorAG or SpliceAcceptor object

    @type  subsequent_acceptor_pos: number
	@param subsequent_acceptor_pos: (estimated) position of preceding acceptor in nt coordinates

    @type  max_leadingexon_nt_length: integer
	@param max_leadingexon_nt_length: positive integer, largest length of leadingexon in nt

    @type  min_leadingexon_nt_length: integer
	@param min_leadingexon_nt_length: positive integer, smallest length of leadingexon in nt

    @type  max_leadingexon_intron_nt_length: integer
	@param max_leadingexon_intron_nt_length: positive integer, largest length of bridging intron in nt

    @type  min_leadingexon_intron_nt_length: integer
	@param min_leadingexon_intron_nt_length: positive integer, smallest length of bridging intron in nt

    @type  min_total_pssm_score: float or None
	@param min_total_pssm_score: minimal sum of donor and acceptor pssm score

    @type  min_donor_pssm_score: float or None
	@param min_donor_pssm_score: minimal donor pssm score

    @type  min_acceptor_pssm_score: float or None
	@param min_acceptor_pssm_score: minimal acceptor pssm score

    @type  allow_non_canonical_donor: Boolean
    @param allow_non_canonical_donor: Boolean (default False)

    @type  allow_non_canonical_acceptor: Boolean
    @param allow_non_canonical_acceptor: Boolean (default False)

    @type  non_canonical_min_donor_pssm_score: float or None
        @param non_canonical_min_donor_pssm_score: minimal donor pssm score

    @type  non_canonical_min_acceptor_pssm_score: float or None
        @param non_canonical_min_acceptor_pssm_score: minimal acceptor pssm score

    @type  min_tss_pssm_score: float or None
	@param min_tss_pssm_score: minimal translational start site pssm score

    @type  order_by: TODO
	@param order_by: TODO

    @rtype:  list
	@return: list of elegiable short leading exons

    @attention: Global vars that have to be set upon usage:
        # DEPRECATED!! ELIGABLE_ALIGNED_SPLICE_SITES_AA_OFFSET
        # and all SHORT_LEADINGEXON variable named
        SHORT_LEADINGEXON_MIN_NT_LENGTH                         
        SHORT_LEADINGEXON_MAX_NT_LENGTH                         
        SHORT_LEADINGEXON_MAX_INTRON_NT_LENGTH                  
        SHORT_LEADINGEXON_MIN_INTRON_NT_LENGTH                  
        SHORT_LEADINGEXON_TSS_MIN_PSSM_SCORE                    
        SHORT_LEADINGEXON_TSS_ALLOW_NON_CANONICAL               
        SHORT_LEADINGEXON_TSS_NON_CANONICAL_MIN_PSSM_SCORE      
        SHORT_LEADINGEXON_MIN_TOTAL_PSSM_SCORE
        SHORT_LEADINGEXON_MIN_DONOR_PSSM_SCORE                  
        SHORT_LEADINGEXON_MIN_ACCEPTOR_PSSM_SCORE               
        SHORT_LEADINGEXON_ALLOW_NON_CANONICAL_DONOR             
        SHORT_LEADINGEXON_ALLOW_NON_CANONICAL_ACCEPTOR          
        SHORT_LEADINGEXON_NON_CANONICAL_MIN_PSSM_SCORE          
        SHORT_LEADINGEXON_NON_CANONICAL_MIN_DONOR_PSSM_SCORE    
        SHORT_LEADINGEXON_NON_CANONICAL_MIN_ACCEPTOR_PSSM_SCORE 

    """

    # check if start of orfD is before end of orfA; if so, return []
    if orfD.startPY > orfA.endPY:
        return []

    # do some input data processing on subsequent_acceptor
    if subsequent_acceptor == None:
        pass
    elif subsequent_acceptor.__class__.__name__ in ['SpliceAcceptorAG','SpliceAcceptor']:
        # set/override variable `subsequent_acceptor_pos` from splice site object
        subsequent_acceptor_pos   = subsequent_acceptor.pos
    else:
        message = "subsequent_acceptor is not a `SpliceAcceptorAG` or `SpliceAcceptor` object"
        raise InproperlyAppliedArgument, message

    # check if subsequent_acceptor_pos is applied/set/overriden
    if not subsequent_acceptor_pos:
        message = "subsequent_acceptor_pos must be specified by either subsequent_acceptor_pos or subsequent_acceptor SpliceAcceptor object"
        raise InproperlyAppliedArgument, message

    # check if subsequent_acceptor has indeed high enough pssm score
    if subsequent_acceptor and (min_acceptor_pssm_score or min_acceptor_pssm_score == 0.0) and\
    subsequent_acceptor.pssm_score < min_acceptor_pssm_score:
        # applied acceptor splicesite doesn't have a high enough pssm_score
        return []

    # check if this orf is correctly positiond towards subsequent_acceptor_pos
    if orfD.startPY > subsequent_acceptor_pos:
        if verbose: print "orfD.startPY > subsequent_acceptor_pos -> return []"
        return []

    # scan orfD for splice sites
    orfD.scan_orf_for_pssm_splice_sites(splicetype="donor",
                forced=True,
                min_pssm_score=min_donor_pssm_score,
                allow_non_canonical=allow_non_canonical_donor,
                non_canonical_min_pssm_score=non_canonical_min_donor_pssm_score)

    # check if there are donor sites found
    if not orfD._donor_sites:
        if verbose: print "orfD has no donor sites -> return []"
        return []

    # scan orfD for TSS pssm sites;
    # set forced=True to redo the search in case done with other settings
    orfD.scan_orf_for_pssm_tss(forced=True,
        min_pssm_score=min_tss_pssm_score)

    # check if there are TSS sites found
    if not orfD._tss_sites:
        if verbose: print "orfD has no TSS sites -> return []"
        return []

    # scan orfA for acceptor splice sites,
    # but only if no `subsequent_acceptor` applied
    eligable_acceptor_sites = []
    if not subsequent_acceptor:
        orfA.scan_orf_for_pssm_splice_sites(splicetype="acceptor",forced=True,
                    min_pssm_score=min_acceptor_pssm_score,
                    allow_non_canonical=SHORT_LEADINGEXON_ALLOW_NON_CANONICAL_ACCEPTOR,
                    non_canonical_min_pssm_score=SHORT_LEADINGEXON_NON_CANONICAL_MIN_ACCEPTOR_PSSM_SCORE)

        # check if there are acceptor sites found
        if not orfA._acceptor_sites:
            if verbose: print "orfA has no acceptor sites -> return []"
            return []

        # gather eligable acceptor sites only
        for acceptor in orfA._acceptor_sites:
            if acceptor.pos <= subsequent_acceptor_pos:
                eligable_acceptor_sites.append( acceptor )
    else:
        eligable_acceptor_sites = [ subsequent_acceptor ]

    if verbose: print "orfD:%s" % orfD.id, [ (d.pos,d.phase,"%1.1f" % d.pssm_score) for d in orfD._donor_sites ]
    if verbose: print "orfD:%s" % orfD.id, [ (t.pos,"%1.1f" % t.pssm_score) for t in orfD._tss_sites ]
    if verbose: print "orfA:%s" % orfA.id, [ (a.pos,a.phase,"%1.1f" % a.pssm_score) for a in orfA._acceptor_sites ]
    if verbose: print "orfA:%s" % orfA.id, subsequent_acceptor_pos, "ELEG", [ (a.pos,a.phase,"%1.1f" % a.pssm_score) for a in eligable_acceptor_sites ]
    if verbose: print "min-max: %s-%s" % (min_leadingexon_nt_length,max_leadingexon_nt_length), "D", min_donor_pssm_score,
    if verbose: print  "A", min_acceptor_pssm_score, "I", min_leadingexon_intron_nt_length, max_leadingexon_intron_nt_length

    # return list with short leading exons and bridging introns
    leadingexons = []

    for tss in orfD._tss_sites:
        for donor in orfD._donor_sites:
            if donor.pos < tss.pos:
                continue
            if (min_donor_pssm_score or min_donor_pssm_score == 0.0) and\
            donor.pssm_score < min_donor_pssm_score:
                # splicesite doesn't have a high enough pssm_score
                continue
            if (donor.pos - tss.pos) < min_leadingexon_nt_length:
                # short atg-exon to short
                continue
            if (donor.pos - tss.pos) > max_leadingexon_nt_length:
                # short atg-exon to long
                continue

            if verbose: print "okay:", tss.pos, donor.pos

            for acceptor in eligable_acceptor_sites:
                if acceptor.phase != donor.phase:
                    # incompatible acceptor and donor phases
                    continue
                if (min_acceptor_pssm_score or min_acceptor_pssm_score == 0.0) and\
                acceptor.pssm_score < min_acceptor_pssm_score:
                    # splicesite doesn't have a high enough pssm_score
                    continue
                if acceptor.pos - donor.pos > max_leadingexon_intron_nt_length:
                    # intron to long
                    continue
                if acceptor.pos - donor.pos < min_leadingexon_intron_nt_length:
                    # intron to short
                    continue
                if (min_total_pssm_score or min_total_pssm_score == 0.0) and\
                    acceptor.pssm_score + donor.pssm_score < min_total_pssm_score:
                    # splicesites together do not have a high enough pssm_score
                    continue

                # make LeadingExonOnOrf object
                exon       = FirstExonOnOrf(tss,donor,orfD)
                shared_nts = get_shared_nucleotides_at_splicesite(orfA,orfD,acceptor,donor)
                intron     = IntronConnectingOrfs(donor,acceptor,shared_nts,orfD,orfA)
                # and append exon and intron to leadingexons list
                leadingexons.append( (exon,intron) )

                # some printing messages
                protstart = tss.pos / 3
                protend   = donor.pos / 3
                seq = orfD.getaas(abs_pos_start=protstart,abs_pos_end=protend)
                logf( orfD.id, orfD, orfD.frame, len(orfD._donor_sites) )
                logf("%s-%s" % (exon.start, exon.end), "\t", "%s-%s" % (intron.start, intron.end), "\t",exon.length, intron.length, "\t", exon.donor.phase, intron.donor.pssm_score, intron.acceptor.pssm_score, seq )

    # return the found exons
    return leadingexons

# end of function find_leading_exon_on_orf

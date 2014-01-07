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
    # length range of short tailing exon and first intron
    SHORT_TAILINGEXON_MIN_NT_LENGTH                         = 10 
    SHORT_TAILINGEXON_MAX_NT_LENGTH                         = 300
    SHORT_TAILINGEXON_MAX_INTRON_NT_LENGTH                  = 300
    SHORT_TAILINGEXON_MIN_INTRON_NT_LENGTH                  = MIN_INTRON_NT_LENGTH
    # thresholds for splice site signals of short tailing exon
    SHORT_TAILINGEXON_MIN_TOTAL_PSSM_SCORE                  = None
    SHORT_TAILINGEXON_MIN_DONOR_PSSM_SCORE                  = 0.0
    SHORT_TAILINGEXON_MIN_ACCEPTOR_PSSM_SCORE               = 0.0
    SHORT_TAILINGEXON_ALLOW_NON_CANONICAL_DONOR             = False     # not implemented yet
    SHORT_TAILINGEXON_ALLOW_NON_CANONICAL_ACCEPTOR          = False     # not implemented yet
    SHORT_TAILINGEXON_NON_CANONICAL_MIN_PSSM_SCORE          = None      # DEPRECATED
    SHORT_TAILINGEXON_NON_CANONICAL_MIN_DONOR_PSSM_SCORE    = float(0)  # not implemented yet
    SHORT_TAILINGEXON_NON_CANONICAL_MIN_ACCEPTOR_PSSM_SCORE = float(0)  # not implemented yet

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



def find_tailing_exon_on_orf(orfD,orfA,order_by='total_pssm',
    max_tailingexon_nt_length=SHORT_TAILINGEXON_MAX_NT_LENGTH,
    min_tailingexon_nt_length=SHORT_TAILINGEXON_MIN_NT_LENGTH,
    max_tailingexon_intron_nt_length=SHORT_TAILINGEXON_MAX_INTRON_NT_LENGTH,
    min_tailingexon_intron_nt_length=SHORT_TAILINGEXON_MIN_INTRON_NT_LENGTH,
    min_total_pssm_score=SHORT_TAILINGEXON_MIN_TOTAL_PSSM_SCORE,
    min_donor_pssm_score=SHORT_TAILINGEXON_MIN_DONOR_PSSM_SCORE,
    min_acceptor_pssm_score=SHORT_TAILINGEXON_MIN_ACCEPTOR_PSSM_SCORE,
    current_donor=None,
    current_donor_pos=None,
    verbose=False):
    """

    @type  orfD: Orf object
	@param orfD: Orf object from where to start the search

    @type  orfA: Orf object
	@param orfA: Orf object to scan for a tailingexon

    @type  current_donor: object
	@param current_donor: SpliceDonorGT or SpliceDonor object

    @type  current_donor_pos: number
	@param current_donor_pos: (estimated) position of preceding donor in nt coordinates

    @type  max_tailingexon_nt_length: integer
	@param max_tailingexon_nt_length: positive integer, largest length of tailingexon in nt

    @type  min_tailingexon_nt_length: integer
	@param min_tailingexon_nt_length: positive integer, smallest length of tailingexon in nt

    @type  max_tailingexon_intron_nt_length: integer
	@param max_tailingexon_intron_nt_length: positive integer, largest length of bridging intron in nt

    @type  min_tailingexon_intron_nt_length: integer
	@param min_tailingexon_intron_nt_length: positive integer, smallest length of bridging intron in nt

    @type  min_total_pssm_score: float or None
	@param min_total_pssm_score: minimal sum of donor and acceptor pssm score

    @type  min_donor_pssm_score: float or None
	@param min_donor_pssm_score: minimal donor pssm score

    @type  min_acceptor_pssm_score: float or None
	@param min_acceptor_pssm_score: minimal acceptor pssm score

    @type  order_by: TODO
	@param order_by: TODO

    @rtype:  list
	@return: list of elegiable short tailing exons

    @attention: Global vars that have to be set upon usage:
        # DEPRECATED!! ELIGABLE_ALIGNED_SPLICE_SITES_AA_OFFSET
        # and all SHORT_TAILINGEXON variable named
        SHORT_TAILINGEXON_MIN_NT_LENGTH                         
        SHORT_TAILINGEXON_MAX_NT_LENGTH                         
        SHORT_TAILINGEXON_MAX_INTRON_NT_LENGTH                  
        SHORT_TAILINGEXON_MIN_INTRON_NT_LENGTH                  
        SHORT_TAILINGEXON_MIN_TOTAL_PSSM_SCORE
        SHORT_TAILINGEXON_MIN_DONOR_PSSM_SCORE                  
        SHORT_TAILINGEXON_MIN_ACCEPTOR_PSSM_SCORE               
        SHORT_TAILINGEXON_ALLOW_NON_CANONICAL_DONOR             
        SHORT_TAILINGEXON_ALLOW_NON_CANONICAL_ACCEPTOR          
        SHORT_TAILINGEXON_NON_CANONICAL_MIN_PSSM_SCORE          
        SHORT_TAILINGEXON_NON_CANONICAL_MIN_DONOR_PSSM_SCORE    
        SHORT_TAILINGEXON_NON_CANONICAL_MIN_ACCEPTOR_PSSM_SCORE 

    """

    # check if end of orfA is before start of orfD; if so, return []
    if orfA.endPY < orfD.startPY:
        return []

    # do some input data processing on current_donor
    if current_donor == None:
        pass
    elif current_donor.__class__.__name__ in ['SpliceDonorGT','SpliceDonor']:
        # set/override variable `current_donor_pos` from splice site object
        current_donor_pos   = current_donor.pos
    else:
        message = "current_donor is not a `SpliceDonorGT` or `SpliceDonor` object"
        raise InproperlyAppliedArgument, message

    # check if current_donor_pos is applied/set/overriden
    if not current_donor_pos:
        message = "current_donor_pos must be specified by either current_donor_pos or current_donor SpliceDonor object"
        raise InproperlyAppliedArgument, message

    # check if current_donor has indeed high enough pssm score
    if current_donor and (min_donor_pssm_score or min_donor_pssm_score == 0.0) and\
    current_donor.pssm_score < min_donor_pssm_score:
        # applied donor splicesite doesn't have a high enough pssm_score
        return []

    # check if this orf is correctly positiond towards subsequent_acceptor_pos
    if orfA.endPY < current_donor_pos:
        if verbose: print "orfA.endPY (%s) < current_donor_pos (%s) -> return []" % (orfA.endPY, current_donor_pos)
        return []

    # scan orfA for splice sites
    orfA.scan_orf_for_pssm_splice_sites(splicetype="acceptor",
                min_pssm_score=min_acceptor_pssm_score,
                allow_non_canonical=SHORT_TAILINGEXON_ALLOW_NON_CANONICAL_ACCEPTOR,
                non_canonical_min_pssm_score=SHORT_TAILINGEXON_NON_CANONICAL_MIN_ACCEPTOR_PSSM_SCORE)

    # check if there are acceptor sites found
    if not orfA._acceptor_sites:
        if verbose: print "orfA no acceptor sites (%s,%s,%s)" % (
                min_acceptor_pssm_score,
                SHORT_TAILINGEXON_ALLOW_NON_CANONICAL_ACCEPTOR,
                SHORT_TAILINGEXON_NON_CANONICAL_MIN_ACCEPTOR_PSSM_SCORE
                )
        return []


    # scan orfD for donor splice sites,
    # but only if no `current_donor` applied
    eligable_donor_sites = []
    if not current_donor:
        orfD.scan_orf_for_pssm_splice_sites(splicetype="donor",forced=True,
                    min_pssm_score=min_donor_pssm_score,
                    allow_non_canonical=SHORT_TAILINGEXON_ALLOW_NON_CANONICAL_DONOR,
                    non_canonical_min_pssm_score=SHORT_TAILINGEXON_NON_CANONICAL_MIN_DONOR_PSSM_SCORE)

        # check if there are acceptor sites found
        if not orfD._donor_sites:
            if verbose: print "orfD no donor sites (%s,%s,%s)" % (
                    min_donor_pssm_score,
                    SHORT_TAILINGEXON_ALLOW_NON_CANONICAL_DONOR,
                    SHORT_TAILINGEXON_NON_CANONICAL_MIN_DONOR_PSSM_SCORE
                    )
            return []

        # gather eligable donor sites only
        for donor in orfD._donor_sites:
            if donor.pos >= current_donor_pos:
                eligable_donor_sites.append( donor )
    else:
        eligable_donor_sites = [ current_donor ]

    if verbose: print "orfD:%s" % orfD.id, [ (d.pos,d.phase,"%1.1f" % d.pssm_score) for d in orfD._donor_sites ]
    if verbose: print "orfD:%s" % orfD.id, current_donor_pos, "ELEG", [ (d.pos,d.phase,"%1.1f" % d.pssm_score) for d in eligable_donor_sites ]
    if verbose: print "orfA:%s" % orfA.id, [ (a.pos,a.phase,"%1.1f" % a.pssm_score) for a in orfA._acceptor_sites ]
    if verbose: print "min-max: %s-%s" % (min_tailingexon_nt_length,max_tailingexon_nt_length), "D", min_donor_pssm_score,
    if verbose: print  "A", min_acceptor_pssm_score, "I", min_tailingexon_intron_nt_length, max_tailingexon_intron_nt_length



    # return list with short tailing exons and bridging introns
    tailingexons = []

    for donor in eligable_donor_sites:
        if (min_donor_pssm_score or min_donor_pssm_score == 0.0) and\
        donor.pssm_score < min_donor_pssm_score:
            # splicesite doesn't have a high enough pssm_score
            continue

        for acceptor in orfA._acceptor_sites:
            if acceptor.phase != donor.phase:
                # incompatible acceptor and donor phases
                continue
            if acceptor.pos < donor.pos + min_tailingexon_intron_nt_length:
                # intron to short
                continue
            if acceptor.pos > donor.pos + max_tailingexon_intron_nt_length:
                # intron to long
                continue
            if orfA.endPY - acceptor.pos < min_tailingexon_nt_length:
                # short exon to short
                continue
            if orfA.endPY - acceptor.pos > max_tailingexon_nt_length:
                # short exon to long
                continue
            if (min_acceptor_pssm_score or min_acceptor_pssm_score == 0.0) and\
            acceptor.pssm_score < min_acceptor_pssm_score:
                # splicesite doesn't have a high enough pssm_score
                continue
            if (min_total_pssm_score or min_total_pssm_score == 0.0) and\
                acceptor.pssm_score + donor.pssm_score < min_total_pssm_score:
                # splicesites together do not have a high enough pssm_score
                continue

            # make LeadingExonOnOrf object
            exon       = FinalExonOnOrf(acceptor,orfA.endPY,orfA)
            shared_nts = get_shared_nucleotides_at_splicesite(orfA,orfD,acceptor,donor)
            intron     = IntronConnectingOrfs(donor,acceptor,shared_nts,orfD,orfA)
            # and append exon and intron to tailingexons list
            tailingexons.append( (exon,intron) )

    # return the found exons
    return tailingexons

# end of function find_tailing_exon_on_orf

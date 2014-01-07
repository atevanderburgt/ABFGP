"""


"""

# Module metadata
__authors__ = "Ate van der Burgt"
__license__ = "MIT"

# Imports
from lib_intron import *
from lib_intron import _order_intron_list
from listofcodingblockgraphs import ListOfCodingBlockGraphs
from dna2prot import dna2proteinbyframe 
from lib_crossdatafunctions import create_pacbpcollectiongraph_from_crossdata
from lib_genestructure import make_consensus_genestructure_from_compatible_pacb_graphs
import pacb


# Global variables
try:
    # import global variables from settings file
    from settings.genestructure import *
except:
    # when using as stand-alone library (not in abfgp), specify these variables
    TINYEXON_MAX_NT_LENGTH                          = 36
    TINYEXON_MIN_NT_LENGTH                          = 1
    TINYEXON_MAX_INTRON_NT_LENGTH                   = 1000
    TINYEXON_MIN_INTRON_NT_LENGTH                   = 3
    TINYEXON_MIN_PSSM_SCORE                         = None      # DEPRECATED
    TINYEXON_MIN_DONOR_PSSM_SCORE                   = 0.0
    TINYEXON_MIN_ACCEPTOR_PSSM_SCORE                = 0.0
    TINYEXON_ALLOW_NON_CANONICAL_DONOR              = False     # not implemented yet 
    TINYEXON_ALLOW_NON_CANONICAL_ACCEPTOR           = False     # not implemented yet
    TINYEXON_NON_CANONICAL_MIN_PSSM_SCORE           = None      # DEPRECATED
    TINYEXON_NON_CANONICAL_MIN_DONOR_PSSM_SCORE     = float(0)  # not implemented yet 
    TINYEXON_NON_CANONICAL_MIN_ACCEPTOR_PSSM_SCORE  = float(0)  # not implemented yet 

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


def scan_orf_for_tiny_exon(orfX,order_by='total_pssm',
    max_tinyexon_nt_length=TINYEXON_MAX_NT_LENGTH,
    min_tinyexon_nt_length=TINYEXON_MIN_NT_LENGTH,
    min_donor_pssm_score=TINYEXON_MIN_DONOR_PSSM_SCORE,
    min_acceptor_pssm_score=TINYEXON_MIN_ACCEPTOR_PSSM_SCORE,
    min_total_pssm_score=TINYEXON_MIN_TOTAL_PSSM_SCORE,
    allow_non_canonical_donor=TINYEXON_ALLOW_NON_CANONICAL_DONOR,
    allow_non_canonical_acceptor=TINYEXON_ALLOW_NON_CANONICAL_ACCEPTOR,
    min_intron_nt_length=None,
    max_intron_nt_length=None,
    donor_phase=None,
    acceptor_phase=None,
    preceeding_donor_site=None,
    subsequent_acceptor_site=None,
    min_acceptor_pos=None,
    max_donor_pos=None):
    """
    Find tiny exons on an orf by length range

    @type  orfX: Orf object
	@param orfX: Orf object to scan for a tinyexon

    @type  max_tinyexon_nt_length: integer
	@param max_tinyexon_nt_length: positive integer, largest length of tinyexon in nt

    @type  min_tinyexon_nt_length: integer
	@param min_tinyexon_nt_length: positive integer, smallest length of tinyexon in nt

    @type  min_total_pssm_score: float or None
	@param min_total_pssm_score: minimal sum of donor - acceptor pssm score pair of tinyexon

    @type  min_donor_pssm_score: float or None
	@param min_donor_pssm_score: minimal donor pssm score of tinyexon

    @type  min_acceptor_pssm_score: float or None
	@param min_acceptor_pssm_score: minimal acceptor pssm score of tinyexon

    @type  max_donor_pos: integer or None
	@param max_donor_pos: maximal elegiable donor position

    @type  min_acceptor_pos: integer or None
	@param min_acceptor_pos: minimal elegiable acceptor position

    @type  order_by: TODO
	@param order_by: TODO
    """

    # scan for splice sites on this (tiny) orf
    orfX.scan_orf_for_pssm_splice_sites(
            splicetype="donor",
            min_pssm_score=min_donor_pssm_score,
            allow_non_canonical=TINYEXON_ALLOW_NON_CANONICAL_DONOR,
            non_canonical_min_pssm_score=TINYEXON_NON_CANONICAL_MIN_DONOR_PSSM_SCORE)
    orfX.scan_orf_for_pssm_splice_sites(
            splicetype="acceptor",
            min_pssm_score=min_acceptor_pssm_score,
            allow_non_canonical=TINYEXON_ALLOW_NON_CANONICAL_ACCEPTOR,
            non_canonical_min_pssm_score=TINYEXON_NON_CANONICAL_MIN_ACCEPTOR_PSSM_SCORE)

    # return list with exons
    exons = []

    # most quickest scan possible: are there donors & acceptors?
    if orfX._donor_sites == [] or orfX._acceptor_sites == []:
        # no exons possible because splice sites are missing
        return exons

    # make a list of compatible_acceptor_sites
    compatible_acceptor_sites = []
    for acceptor in orfX._acceptor_sites:
        if acceptor_phase in [0,1,2] and acceptor.phase != acceptor_phase:
            continue
        if acceptor.pssm_score < min_acceptor_pssm_score:
            continue
        if min_acceptor_pos and acceptor.pos < min_acceptor_pos:
            continue
        if preceeding_donor_site:
            if preceeding_donor_site.phase != acceptor.phase:
                continue
            if min_intron_nt_length and acceptor.pos - preceeding_donor_site.pos < min_intron_nt_length:
                continue
            if max_intron_nt_length and acceptor.pos - preceeding_donor_site.pos > max_intron_nt_length:
                continue

        # if we reach this point, compatible site!
        compatible_acceptor_sites.append( acceptor )

    # make a list of compatible_donor_sites
    compatible_donor_sites = []
    for donor in orfX._donor_sites:
        if donor_phase in [0,1,2] and donor.phase != donor_phase:
            continue
        if donor.pssm_score < min_donor_pssm_score:
            continue
        if max_donor_pos and donor.pos > max_donor_pos:
            continue
        if subsequent_acceptor_site:
            if subsequent_acceptor_site.phase != donor.phase:
                continue
            if min_intron_nt_length and subsequent_acceptor_site.pos - donor.pos < min_intron_nt_length:
                continue
            if max_intron_nt_length and subsequent_acceptor_site.pos - donor.pos > max_intron_nt_length:
                continue
        # if we reach this point, compatible site!
        compatible_donor_sites.append( donor )

    ###print "lib_tinyexon, comp d & a:", len(compatible_donor_sites), len(compatible_acceptor_sites), "orf:", orfX.id, min_donor_pssm_score, min_acceptor_pssm_score


    # and combine sites to exons!
    for acceptor in compatible_acceptor_sites:
        for donor in compatible_donor_sites:
            # length of exon
            exon_length = donor.pos - acceptor.pos
            # continue if exon to short
            if exon_length < min_tinyexon_nt_length: continue
            # continue if exon to long
            if exon_length > max_tinyexon_nt_length: continue

            # check sum of donor and acceptor pssm score
            if (min_total_pssm_score or min_total_pssm_score==0.0) and\
            donor.pssm_score + acceptor.pssm_score < min_total_pssm_score:
                continue

            # make a Exon object
            exon = ExonOnOrf(acceptor,donor,orfX)
            exons.append(exon)

    # return ordered exon list
    return _order_intron_list(exons,order_by=order_by)

# end of function scan_orf_for_tiny_exon



def find_tiny_exon_on_orf(orfX,order_by='total_pssm',
    max_tinyexon_nt_length=TINYEXON_MAX_NT_LENGTH,
    min_tinyexon_nt_length=TINYEXON_MIN_NT_LENGTH,
    max_tinyexon_intron_nt_length=TINYEXON_MAX_INTRON_NT_LENGTH,
    min_tinyexon_intron_nt_length=TINYEXON_MIN_INTRON_NT_LENGTH,
    min_donor_pssm_score=TINYEXON_MIN_DONOR_PSSM_SCORE,
    min_acceptor_pssm_score=TINYEXON_MIN_ACCEPTOR_PSSM_SCORE,
    min_total_pssm_score=TINYEXON_MIN_TOTAL_PSSM_SCORE,
    preceding_donor=None,
    subsequent_acceptor=None,
    preceding_donor_pos=None,
    subsequent_acceptor_pos=None):
    """
    Find a tiny exon on an orf by a leading donor and a trailing acceptor site.

    @type  orfX: Orf object
	@param orfX: Orf object to scan for a tinyexon

    @type  preceding_donor: object
	@param preceding_donor: SpliceDonorGT or SpliceDonor object

    @type  subsequent_acceptor: object
	@param subsequent_acceptor: SpliceAcceptorAG or SpliceAcceptor object

    @type  max_tinyexon_nt_length: integer
	@param max_tinyexon_nt_length: positive integer, largest length of tinyexon in nt

    @type  min_tinyexon_nt_length: integer
	@param min_tinyexon_nt_length: positive integer, smallest length of tinyexon in nt

    @type  max_tinyexon_intron_nt_length: integer
    @param max_tinyexon_intron_nt_length: positive integer, largest length of intron around tinyexon in nt

    @type  min_tinyexon_intron_nt_length: integer
    @param min_tinyexon_intron_nt_length: positive integer, smallest length of intron around tinyexon in nt

    @type  min_total_pssm_score: float or None
	@param min_total_pssm_score: minimal sum of donor - acceptor pssm score pair of tinyexon

    @type  min_donor_pssm_score: float or None
	@param min_donor_pssm_score: minimal donor pssm score of tinyexon

    @type  min_acceptor_pssm_score: float or None
	@param min_acceptor_pssm_score: minimal acceptor pssm score of tinyexon

    @type  order_by: TODO
	@param order_by: TODO

    @attention: Global vars that have to be set upon usage:
        MIN_DONOR_PSSM_SCORE
        MIN_ACCEPTOR_PSSM_SCORE
        # and all TINYEXON variable named
        TINYEXON_MAX_NT_LENGTH                          
        TINYEXON_MIN_NT_LENGTH                          
        TINYEXON_MAX_INTRON_NT_LENGTH                   
        TINYEXON_MIN_INTRON_NT_LENGTH                   
        TINYEXON_MIN_PSSM_SCORE                         
        TINYEXON_MIN_DONOR_PSSM_SCORE                   
        TINYEXON_MIN_ACCEPTOR_PSSM_SCORE                
        TINYEXON_ALLOW_NON_CANONICAL_DONOR              
        TINYEXON_ALLOW_NON_CANONICAL_ACCEPTOR           
        TINYEXON_NON_CANONICAL_MIN_PSSM_SCORE           
        TINYEXON_NON_CANONICAL_MIN_DONOR_PSSM_SCORE     
        TINYEXON_NON_CANONICAL_MIN_ACCEPTOR_PSSM_SCORE  

    """

    # scan for splice sites on this (tiny) orf
    orfX.scan_orf_for_pssm_splice_sites(
            splicetype="donor",
            min_pssm_score=TINYEXON_MIN_DONOR_PSSM_SCORE,
            allow_non_canonical=TINYEXON_ALLOW_NON_CANONICAL_DONOR,
            non_canonical_min_pssm_score=TINYEXON_NON_CANONICAL_MIN_DONOR_PSSM_SCORE)
    orfX.scan_orf_for_pssm_splice_sites(
            splicetype="acceptor",
            min_pssm_score=TINYEXON_MIN_ACCEPTOR_PSSM_SCORE,
            allow_non_canonical=TINYEXON_ALLOW_NON_CANONICAL_ACCEPTOR,
            non_canonical_min_pssm_score=TINYEXON_NON_CANONICAL_MIN_ACCEPTOR_PSSM_SCORE)

    # return list with exons
    exons = []


    # do some input data processing on preceding_donor
    if preceding_donor == None:
        # preceding donor MUST be set!
        message = "preceding_donor is not a `SpliceDonorGT` or `SpliceDonor` object"
        raise InproperlyAppliedArgument, message
    elif preceding_donor.__class__.__name__ in ['SpliceDonorGT','SpliceDonor']:
        pass
    else:
        message = "preceding_donor is not a `SpliceDonorGT` or `SpliceDonor` object, but a `%s`" % preceding_donor.__class__.__name__
        raise InproperlyAppliedArgument, message

    # do some input data processing on subsequent_acceptor
    if subsequent_acceptor == None:
        # subsequent acceptor MUST be set
        message = "subsequent_acceptor is not a `SpliceAcceptorAG` or `SpliceAcceptor` object"
        raise InproperlyAppliedArgument, message
    elif subsequent_acceptor.__class__.__name__ in ['SpliceAcceptorAG','SpliceAcceptor']:
        pass
    else:
        message = "subsequent_acceptor is not a `SpliceAcceptorAG` or `SpliceAcceptor` object, but a `%s`" % subsequent_acceptor.__class__.__name__
        raise InproperlyAppliedArgument, message

    # check phases of acceptor and donor
    if subsequent_acceptor.phase not in [0,1,2]:
        raise UnexpectedSpliceSitePhase
    if preceding_donor.phase not in [0,1,2]:
        raise UnexpectedSpliceSitePhase

    # some further integrity check on integer arguments
    for variable in ( max_tinyexon_nt_length, min_tinyexon_nt_length,
    max_tinyexon_intron_nt_length, min_tinyexon_intron_nt_length):
        try:
            variable = int(variable)
            if variable <= 0:
                raise 
        except:
            message = "a variable is NOT a positive integer as expected"
            raise InproperlyAppliedArgument, message


    # most quickest scan possible: are there donors & acceptors?
    if orfX._donor_sites == [] or orfX._acceptor_sites == []:
        # no exons possible because splice sites are missing
        return exons

    # make a list of compatible_acceptor_sites
    compatible_acceptor_sites = []
    for acceptor in orfX._acceptor_sites:
        # TODO: check! do we need a combi of donor and acceptor or acceptor and acceptor?
        if acceptor.phase != preceding_donor.phase:
            continue
        if acceptor.pssm_score < min_acceptor_pssm_score:
            continue
        if acceptor.pos - preceding_donor.pos < min_tinyexon_intron_nt_length:
            # intron to short
            continue
        if acceptor.pos - preceding_donor.pos > max_tinyexon_intron_nt_length:
            # intron to long
            continue
        # if we reach this point, compatible site!
        compatible_acceptor_sites.append( acceptor )

    # make a list of compatible_donor_sites
    compatible_donor_sites = []
    for donor in orfX._donor_sites:
        # TODO: check! do we need a combi of donor and acceptor or donor and donor?
        if donor.phase != subsequent_acceptor.phase:
            continue
        if donor.pssm_score < min_donor_pssm_score:
            continue
        if subsequent_acceptor.pos - donor.pos > max_tinyexon_intron_nt_length:
            # intron to long
            continue
        if subsequent_acceptor.pos - donor.pos < min_tinyexon_intron_nt_length:
            # intron to short
            continue
        # if we reach this point, compatible site!
        compatible_donor_sites.append( donor )

    # and combine sites to exons!
    for acceptor in compatible_acceptor_sites:
        for donor in compatible_donor_sites:
            # length of exon
            exon_length = donor.pos - acceptor.pos
            # continue if exon to short
            if exon_length < min_tinyexon_nt_length: continue
            # continue if exon to long
            if exon_length > max_tinyexon_nt_length: continue

            # check sum of donor and acceptor pssm score
            if (min_total_pssm_score or min_total_pssm_score==0.0) and\
            donor.pssm_score + acceptor.pssm_score < min_total_pssm_score:
                continue

            # make a Exon object
            exon = ExonOnOrf(acceptor,donor,orfX)
            exons.append(exon)

    # return ordered exon list
    return _order_intron_list(exons,order_by=order_by)

# end of function find_tiny_exon_on_orf


def bridge_two_pacbporfs_by_tinyexon(preceding_orf,subsequent_orf,
    preceding_donor_sites=[],
    subsequent_acceptor_sites=[],
    orflist=[],
    max_tinyexon_nt_length=TINYEXON_MAX_NT_LENGTH,
    min_tinyexon_nt_length=TINYEXON_MIN_NT_LENGTH,
    max_tinyexon_intron_nt_length=TINYEXON_MAX_INTRON_NT_LENGTH,
    min_tinyexon_intron_nt_length=TINYEXON_MIN_INTRON_NT_LENGTH,
    min_donor_pssm_score=TINYEXON_MIN_DONOR_PSSM_SCORE,
    min_acceptor_pssm_score=TINYEXON_MIN_ACCEPTOR_PSSM_SCORE,
    min_total_pssm_score=TINYEXON_MIN_TOTAL_PSSM_SCORE):
    """
    Bridge two `neighbouring` Orfs by a tinyexon by applying preceding donors and subsequent acceptors

    @type  preceding_orf: Orf object
	@param preceding_orf: Orf object that contains preceding_donor_site(s)

    @type  subsequent_orf: Orf object
	@param subsequent_orf: Orf object that contains subsequent_acceptor_site(s)

    @type  preceding_donor_sites: list
	@param preceding_donor_sites: list with SpliceDonorGT and/or SpliceDonor objects

    @type  subsequent_acceptor_sites: list
	@param subsequent_acceptor_sites: list with SpliceAcceptorAG and/or SpliceAcceptor objects

    @type  orflist: list
	@param orflist: list with Orf objects

    @type  max_tinyexon_nt_length: integer
	@param max_tinyexon_nt_length: positive integer, largest length of tinyexon in nt

    @type  min_tinyexon_nt_length: integer
	@param min_tinyexon_nt_length: positive integer, smallest length of tinyexon in nt

    @type  max_tinyexon_intron_nt_length: integer
    @param max_tinyexon_intron_nt_length: positive integer, largest length of intron around tinyexon in nt

    @type  min_tinyexon_intron_nt_length: integer
    @param min_tinyexon_intron_nt_length: positive integer, smallest length of intron around tinyexon in nt

    @type  min_total_pssm_score: float or None
	@param min_total_pssm_score: minimal sum of donor - acceptor pssm score pair of tinyexon

    @type  min_donor_pssm_score: float or None
	@param min_donor_pssm_score: minimal donor pssm score of tinyexon

    @type  min_acceptor_pssm_score: float or None
	@param min_acceptor_pssm_score: minimal acceptor pssm score of tinyexon

    @rtype:  list
	@return: list of tuples ( preceding_intron, tinyexon, subsequent_intron )

    @attention: Global vars that have to be set upon usage:
        MIN_DONOR_PSSM_SCORE
        MIN_ACCEPTOR_PSSM_SCORE
        # and all TINYEXON variable named
        TINYEXON_MAX_NT_LENGTH                          
        TINYEXON_MIN_NT_LENGTH                          
        TINYEXON_MAX_INTRON_NT_LENGTH                   
        TINYEXON_MIN_INTRON_NT_LENGTH                   
        TINYEXON_MIN_PSSM_SCORE                         
        TINYEXON_MIN_DONOR_PSSM_SCORE                   
        TINYEXON_MIN_ACCEPTOR_PSSM_SCORE                
        TINYEXON_ALLOW_NON_CANONICAL_DONOR              
        TINYEXON_ALLOW_NON_CANONICAL_ACCEPTOR           
        TINYEXON_NON_CANONICAL_MIN_PSSM_SCORE           
        TINYEXON_NON_CANONICAL_MIN_DONOR_PSSM_SCORE     
        TINYEXON_NON_CANONICAL_MIN_ACCEPTOR_PSSM_SCORE  

    """
    if not preceding_donor_sites:
        return []
    if not subsequent_acceptor_sites:
        return []
    if not orflist:
        return []

    # return dictionary with exon coordinates (keys) and exons/introns
    returnexons = {}
    min_preceding_donor_sites_pos     = min([ d.pos for d in preceding_donor_sites ])
    max_subsequent_acceptor_sites_pos = max([ a.pos for a in subsequent_acceptor_sites ]) 
    for orfX in orflist:
        # check if orf is correctly positions towards the splice sites' extremes
        if orfX.endPY   <= min_preceding_donor_sites_pos: continue
        if orfX.startPY >= max_subsequent_acceptor_sites_pos: continue

        # if here, we can try to make a bridge by a tinyexon
        for donor in preceding_donor_sites:
            # orf not correctly positions towards the donor site
            if orfX.endPY <= donor.pos: continue

            # check pssm_score of donor site
            # TODO: this is in fact the donor on the normal, large orf
            # TODO: do we want to check this pssm score?
            if donor.pssm_score < min_donor_pssm_score: continue

            for acceptor in subsequent_acceptor_sites:
                if orfX.startPY >= acceptor.pos: continue

                # check pssm_score of acceptor site
                # TODO: this is in fact the acceptor on the normal, large orf
                # TODO: do we want to check this pssm score?
                if acceptor.pssm_score < min_acceptor_pssm_score: continue

                # okay, now try to bridge it!
                exons = find_tiny_exon_on_orf(orfX,order_by='total_pssm',
                        max_tinyexon_nt_length=max_tinyexon_nt_length,
                        min_tinyexon_nt_length=min_tinyexon_nt_length,
                        max_tinyexon_intron_nt_length=max_tinyexon_intron_nt_length,
                        min_tinyexon_intron_nt_length=min_tinyexon_intron_nt_length,
                        min_donor_pssm_score=min_donor_pssm_score,
                        min_acceptor_pssm_score=min_acceptor_pssm_score,
                        min_total_pssm_score=min_total_pssm_score,
                        preceding_donor=donor,
                        subsequent_acceptor=acceptor
                        )
                # and append to returnexons
                for tinyexon in exons:

                    # make preceding intron
                    shared_nts_A = "TODO"
                    preceding_intron = IntronConnectingOrfs(
                        donor,tinyexon.acceptor,
                        shared_nts_A,preceding_orf,tinyexon.orf )

                    # make subsequent intron
                    shared_nts_B = "TODO"
                    subsequent_intron = IntronConnectingOrfs(
                        tinyexon.donor, acceptor,
                        shared_nts_B,tinyexon.orf,subsequent_orf )

                    # and append to exons
                    key = ( tinyexon.acceptor.pos, tinyexon.donor.pos )
                    #returnexons.append( ( preceding_intron, tinyexon, subsequent_intron ) )
                    if key not in returnexons.keys():
                        returnexons[key] = tinyexon

    # and return the list of intron/exon/intron
    return _order_intron_list( returnexons.values() )

# end of function bridge_two_pacbporfs_by_tinyexon





def bridge_two_pacbporfs_by_two_tinyexons(preceding_orf,subsequent_orf,
    preceding_donor_sites=[],
    subsequent_acceptor_sites=[],
    orflist=[],
    max_tinyexon_nt_length=TINYEXON_MAX_NT_LENGTH,
    min_tinyexon_nt_length=TINYEXON_MIN_NT_LENGTH,
    max_tinyexon_intron_nt_length=TINYEXON_MAX_INTRON_NT_LENGTH,
    min_tinyexon_intron_nt_length=TINYEXON_MIN_INTRON_NT_LENGTH,
    min_donor_pssm_score=TINYEXON_MIN_DONOR_PSSM_SCORE,
    min_acceptor_pssm_score=TINYEXON_MIN_ACCEPTOR_PSSM_SCORE,
    min_total_pssm_score=TINYEXON_MIN_TOTAL_PSSM_SCORE):
    """
    Bridge two `neighbouring` Orfs by TWO tinyexons by applying preceding donors and subsequent acceptors

    @type  preceding_orf: Orf object
	@param preceding_orf: Orf object that contains preceding_donor_site(s)

    @type  subsequent_orf: Orf object
	@param subsequent_orf: Orf object that contains subsequent_acceptor_site(s)

    @type  preceding_donor_sites: list
	@param preceding_donor_sites: list with SpliceDonorGT and/or SpliceDonor objects

    @type  subsequent_acceptor_sites: list
	@param subsequent_acceptor_sites: list with SpliceAcceptorAG and/or SpliceAcceptor objects

    @type  orflist: list
	@param orflist: list with Orf objects

    @type  max_tinyexon_nt_length: integer
	@param max_tinyexon_nt_length: positive integer, largest length of tinyexon in nt

    @type  min_tinyexon_nt_length: integer
	@param min_tinyexon_nt_length: positive integer, smallest length of tinyexon in nt

    @type  max_tinyexon_intron_nt_length: integer
    @param max_tinyexon_intron_nt_length: positive integer, largest length of intron around tinyexon in nt

    @type  min_tinyexon_intron_nt_length: integer
    @param min_tinyexon_intron_nt_length: positive integer, smallest length of intron around tinyexon in nt

    @type  min_total_pssm_score: float or None
	@param min_total_pssm_score: minimal sum of donor - acceptor pssm score pair of tinyexon

    @type  min_donor_pssm_score: float or None
	@param min_donor_pssm_score: minimal donor pssm score of tinyexon

    @type  min_acceptor_pssm_score: float or None
	@param min_acceptor_pssm_score: minimal acceptor pssm score of tinyexon

    @rtype:  list
	@return: list of tuples ( preceding_intron, tinyexon, subsequent_intron )

    @attention: Global vars that have to be set upon usage:
        MIN_DONOR_PSSM_SCORE
        MIN_ACCEPTOR_PSSM_SCORE
        # and all TINYEXON variable named
        TINYEXON_MAX_NT_LENGTH                          
        TINYEXON_MIN_NT_LENGTH                          
        TINYEXON_MAX_INTRON_NT_LENGTH                   
        TINYEXON_MIN_INTRON_NT_LENGTH                   
        TINYEXON_MIN_PSSM_SCORE                         
        TINYEXON_MIN_DONOR_PSSM_SCORE                   
        TINYEXON_MIN_ACCEPTOR_PSSM_SCORE                
        TINYEXON_ALLOW_NON_CANONICAL_DONOR              
        TINYEXON_ALLOW_NON_CANONICAL_ACCEPTOR           
        TINYEXON_NON_CANONICAL_MIN_PSSM_SCORE           
        TINYEXON_NON_CANONICAL_MIN_DONOR_PSSM_SCORE     
        TINYEXON_NON_CANONICAL_MIN_ACCEPTOR_PSSM_SCORE  

    """
    if not preceding_donor_sites:
        return []
    if not subsequent_acceptor_sites:
        return []
    if not orflist:
        return []

    # make the first tinyexon.
    tinyexonlist = bridge_two_pacbporfs_by_tinyexon(preceding_orf,subsequent_orf,
            preceding_donor_sites=preceding_donor_sites,
            subsequent_acceptor_sites=subsequent_acceptor_sites,
            orflist=orflist,
            max_tinyexon_nt_length=TINYEXON_MAX_NT_LENGTH,
            min_tinyexon_nt_length=TINYEXON_MIN_NT_LENGTH,
            max_tinyexon_intron_nt_length=TINYEXON_MAX_INTRON_NT_LENGTH,
            min_tinyexon_intron_nt_length=TINYEXON_MIN_INTRON_NT_LENGTH,
            min_donor_pssm_score=TINYEXON_MIN_DONOR_PSSM_SCORE,
            min_acceptor_pssm_score=TINYEXON_MIN_ACCEPTOR_PSSM_SCORE,
            min_total_pssm_score=TINYEXON_MIN_TOTAL_PSSM_SCORE)


    total_list = []
    for tinyexon in tinyexonlist:

        # try to bridge another one RIGTH of this tinyexon
        right_tinyexonlist = bridge_two_pacbporfs_by_tinyexon(tinyexon.orf,subsequent_orf,
                preceding_donor_sites=[ tinyexon.donor ],
                subsequent_acceptor_sites=subsequent_acceptor_sites,
                orflist=orflist,
                max_tinyexon_nt_length=TINYEXON_MAX_NT_LENGTH,
                min_tinyexon_nt_length=TINYEXON_MIN_NT_LENGTH,
                max_tinyexon_intron_nt_length=TINYEXON_MAX_INTRON_NT_LENGTH,
                min_tinyexon_intron_nt_length=TINYEXON_MIN_INTRON_NT_LENGTH,
                min_donor_pssm_score=TINYEXON_MIN_DONOR_PSSM_SCORE,
                min_acceptor_pssm_score=TINYEXON_MIN_ACCEPTOR_PSSM_SCORE,
                min_total_pssm_score=TINYEXON_MIN_TOTAL_PSSM_SCORE)

        # try to bridge another one LEFT of this tinyexon
        left_tinyexonlist = bridge_two_pacbporfs_by_tinyexon(preceding_orf,tinyexon.orf,
                preceding_donor_sites=preceding_donor_sites,
                subsequent_acceptor_sites=[tinyexon.acceptor],
                orflist=orflist,
                max_tinyexon_nt_length=TINYEXON_MAX_NT_LENGTH,
                min_tinyexon_nt_length=TINYEXON_MIN_NT_LENGTH,
                max_tinyexon_intron_nt_length=TINYEXON_MAX_INTRON_NT_LENGTH,
                min_tinyexon_intron_nt_length=TINYEXON_MIN_INTRON_NT_LENGTH,
                min_donor_pssm_score=TINYEXON_MIN_DONOR_PSSM_SCORE,
                min_acceptor_pssm_score=TINYEXON_MIN_ACCEPTOR_PSSM_SCORE,
                min_total_pssm_score=TINYEXON_MIN_TOTAL_PSSM_SCORE)
 
        # append to total list
        total_list.extend([ (tinyexon,te) for te in right_tinyexonlist])
        total_list.extend([ (te,tinyexon) for te in left_tinyexonlist])

    # make total_list unique and order it
    tmp = {}
    for e1,e2 in total_list:
        coordkey = ( e1.acceptor.pos, e1.donor.pos, e2.acceptor.pos, e2.donor.pos )
        if not tmp.has_key(coordkey):
            sumscore = - ( e1.acceptor.pssm_score + e1.donor.pssm_score + e2.acceptor.pssm_score + e2.donor.pssm_score )
            tmp[coordkey] = ( sumscore, e1, e2 )
    elist = tmp.values()
    elist.sort()
    # and return the list
    return [ (e1,e2) for (score,e1,e2) in elist ]


# end of function bridge_two_pacbporfs_by_two_tinyexons



def _best_splicesite_by_phase(phase,sites):
    """
    Helper function in tinyexon prediction; returns highest pssm site in a range of sites of a specific phase

    @type  phase: int
	@param phase: requested phase; any of (0,1,2)

    @type  sites: list
	@param sites: list of SpliceAcceptors or SpliceDonors; objects must have attributes ``pssm_score`` and ``pos``

    @rtype:  SpliceAcceptor or SpliceDonor object
	@return: highest pssm site in a range of sites of a specific phase
    """
    correct_phase = []
    for site in sites:
        if site.phase == phase:
            correct_phase.append( (-site.pssm_score,site) )
    correct_phase.sort()
    return correct_phase[0][1]

# end of _best_splicesite_by_phase


def find_intermediary_codingblockgraph_with_tinyexon(graphL,graphR,input={},similaritymatrix=None,min_bitscore_ratio=0.3):
    """
    """
    tinyexon_crossdata = {}
    tinyexons_seen = 0

    for org in graphL.organism_set():
        theOrfL = graphL.get_orfs_of_graph(organism=org)[0]
        theOrfR = graphR.get_orfs_of_graph(organism=org)[0]
        # continue if identical orfs
        # TODO: maybe check as well for spanning ranges?
        # TODO: in theory, a tinyorf can exist on this orf as well...
        if theOrfL.id == theOrfR.id: continue
        msrL = graphL.minimal_spanning_range(organism=org)
        msrR = graphR.minimal_spanning_range(organism=org)

        # check for get eligable donors on orfL and acceptors on orfR
        if org in graphL._splicedonorgraph.organism_set() and\
        org in graphR._spliceacceptorgraph.organism_set():
            eligable_donors    = graphL._splicedonorgraph.get_organism_objects(org)
            eligable_acceptors = graphR._spliceacceptorgraph.get_organism_objects(org)
            orflist            = input[org]['orfs'].orfs

            # search for tinyexons
            tinyexonlist = bridge_two_pacbporfs_by_tinyexon(theOrfL,theOrfR,
                    preceding_donor_sites= eligable_donors,
                    subsequent_acceptor_sites= eligable_acceptors,
                    orflist=orflist
                    )

            doubletinyexons = bridge_two_pacbporfs_by_two_tinyexons(theOrfL,theOrfR,
                    preceding_donor_sites= eligable_donors,
                    subsequent_acceptor_sites= eligable_acceptors,
                    orflist=orflist
                    )

        else:
            # not donors and acceptors on both orfs!
            return []


        # Order the tinyexons with respect to which orf they are located on.
        # for now, IGNORE tinyexons on the both left and right orf it self!!
        orf2tinyexons = {}
        for tinyexon in tinyexonlist:
            if tinyexon.orf.id in [ theOrfL.id, theOrfR.id ]:
                continue
            if orf2tinyexons.has_key(tinyexon.orf.id):
                orf2tinyexons[tinyexon.orf.id].append(tinyexon)
            else:
                orf2tinyexons[tinyexon.orf.id] = [ tinyexon ]

        # loop over the unique orfids on which tinyexons are predicted
        for orfid, telist in orf2tinyexons.iteritems():
            # loop over all other organisms (except the organism itself)
            for otherorg in graphL.organism_set():
                if otherorg == org: continue
                orgkey = [org,otherorg]
                orgkey.sort()
                _orgkey_reversed = False
                if orgkey != [org,otherorg]: _orgkey_reversed = True
                orgkey = tuple(orgkey)
                if not tinyexon_crossdata.has_key(orgkey):
                    tinyexon_crossdata[orgkey] = {'accepted_pacbs': {} }
                orfL = graphL.get_orfs_of_graph(organism=otherorg)[0]
                orfR = graphR.get_orfs_of_graph(organism=otherorg)[0]

                # main list for all similarities on this orfid
                similaritiesL = []
                similaritiesR = []
                for tinyexon in telist:
                    # get protein query sequence from tinyorf
                    query_dna    = tinyexon.orf.inputgenomicsequence[tinyexon.acceptor.pos:tinyexon.donor.pos]
                    query        = dna2proteinbyframe(query_dna, (3 - tinyexon.acceptor.phase) % 3 )
                    query_aa_pos = tinyexon.acceptor.pos / 3

                    _similaritiesL = similaritymatrix.scansbjct(query,orfL.protein_sequence,min_bitscore_ratio=min_bitscore_ratio)
                    if orfL.id == orfR.id:
                        _similaritiesR = []
                    else:
                        _similaritiesR = similaritymatrix.scansbjct(query,orfR.protein_sequence,min_bitscore_ratio=min_bitscore_ratio)
                    # Append to all similarities on this orfid; append the tinyexon itself too
                    # in order to place the similarity back to a specific tinyexon.
                    # This is needed because there can be >1 tinyexon on the same orf...
                    _similaritiesL = [ (_data,tinyexon) for _data in _similaritiesL ]
                    _similaritiesR = [ (_data,tinyexon) for _data in _similaritiesR ]
                    similaritiesL.extend(_similaritiesL)
                    similaritiesR.extend(_similaritiesR)

                # re-order the similarities because they can contain data from 2 tinyexons (on the same orf)
                # ordering is performed on ``ratio * bitscore``
                # this - kind of - evalue calculation enables a preferation for longer matches
                similaritiesL = _order_similarities(similaritiesL)
                similaritiesR = _order_similarities(similaritiesR)

                # Now make pacbporfs of only the BEST tinyexon and its
                # similarity on another organism
                TAKE_BEST_SIMILARITIES = 2
                for ( ( ratio, sbjct_pos, q_seq, match, s_seq, bitscore), tinyexon ) in similaritiesL[0:TAKE_BEST_SIMILARITIES]:
                    sbjct_aa_pos = sbjct_pos+orfL.protein_startPY
                    query_aa_pos = tinyexon.acceptor.pos / 3
                    if _orgkey_reversed:
                        ###print s_seq, "'%s'" % match, ratio, orfL.id
                        pacbpkey = (bitscore, len(query), orfL.id, tinyexon.orf.id )
                        pacbp    = pacb.PacbP(input=(s_seq,q_seq,sbjct_aa_pos,query_aa_pos))
                        pacbporf = pacb.conversion.pacbp2pacbporf(pacbp,orfL,tinyexon.orf)
                    else:
                        ###print q_seq, "'%s'" % match, ratio, orfL.id
                        pacbpkey = (bitscore, len(query), tinyexon.orf.id, orfL.id )
                        pacbp    = pacb.PacbP(input=(q_seq,s_seq,query_aa_pos,sbjct_aa_pos))
                        pacbporf = pacb.conversion.pacbp2pacbporf(pacbp,tinyexon.orf,orfL)

                    tinyexons_seen+=1
                    pacbporf.extend_pacbporf_after_stops()
                    tinyexon_crossdata[orgkey]['accepted_pacbs'][pacbpkey] = pacbporf

                for ( ( ratio, sbjct_pos, q_seq, match, s_seq, bitscore), tinyexon ) in similaritiesR[0:TAKE_BEST_SIMILARITIES]:
                    sbjct_aa_pos = sbjct_pos+orfR.protein_startPY
                    query_aa_pos = tinyexon.acceptor.pos / 3
                    if _orgkey_reversed:
                        pacbpkey = (bitscore, len(query), orfR.id, tinyexon.orf.id )
                        pacbp    = pacb.PacbP(input=(s_seq,q_seq,sbjct_aa_pos,query_aa_pos))
                        pacbporf = pacb.conversion.pacbp2pacbporf(pacbp,orfR,tinyexon.orf)
                    else:
                        pacbpkey = (bitscore, len(query), tinyexon.orf.id, orfR.id )
                        pacbp    = pacb.PacbP(input=(q_seq,s_seq,query_aa_pos,sbjct_aa_pos))
                        pacbporf = pacb.conversion.pacbp2pacbporf(pacbp,tinyexon.orf,orfR)

                    tinyexons_seen+=1
                    pacbporf.extend_pacbporf_after_stops()
                    tinyexon_crossdata[orgkey]['accepted_pacbs'][pacbpkey] = pacbporf


    if not tinyexons_seen:
        return []
    else:
        # add the nodes/edges from the input graphs as well
        for ( (a,b,c,d),n1,n2 ) in graphL.pacbps.keys():
            orgkey   = ( n1[0], n2[0] )
            pacbpdna = graphL.pacbps[( (a,b,c,d),n1,n2 )]
            if not tinyexon_crossdata.has_key(orgkey):
                tinyexon_crossdata[orgkey] = {'accepted_pacbs': {} }
            tinyexon_crossdata[orgkey]['accepted_pacbs'][(a,b,c,d)] = pacbpdna

        for ( (a,b,c,d),n1,n2 ) in graphR.pacbps.keys():
            orgkey   = ( n1[0], n2[0] )
            pacbpdna = graphR.pacbps[( (a,b,c,d),n1,n2 )]
            if not tinyexon_crossdata.has_key(orgkey):
                tinyexon_crossdata[orgkey] = {'accepted_pacbs': {} }
            tinyexon_crossdata[orgkey]['accepted_pacbs'][(a,b,c,d)] = pacbpdna
    

    # make graph, remove to low connected nodes and split in complete graphs
    tinyexonsg = create_pacbpcollectiongraph_from_crossdata(tinyexon_crossdata)
    tinyexonsg.remove_low_connectivity_nodes(min_connectivity=2)
    splitted_tinyexongraphs = tinyexonsg.find_fully_connected_subgraphs(
                edges=4,
                max_missing_edges=0 )

    # now remove the graphs that are graphL and graphR ;-)
    graphLnodes = graphL.get_nodes()
    graphLnodes.sort()
    graphRnodes = graphR.get_nodes()
    graphRnodes.sort()
    for pos in range(0,len(splitted_tinyexongraphs)):
        tegNodes = splitted_tinyexongraphs[pos].get_nodes()
        tegNodes.sort()
        if tegNodes == graphLnodes:
            splitted_tinyexongraphs.pop(pos)
            break
    for pos in range(0,len(splitted_tinyexongraphs)):
        tegNodes = splitted_tinyexongraphs[pos].get_nodes()
        tegNodes.sort()
        if tegNodes == graphRnodes:
            splitted_tinyexongraphs.pop(pos)
            break

    # make ListOfCodingBlockGraphs
    cbgList = ListOfCodingBlockGraphs(splitted_tinyexongraphs,
            input=input,
            crossdata=tinyexon_crossdata
            )

    # do all what is needed to create K(s) CBGs of these
    cbgList.harvest_pacbps_from_crossdata()
    cbgList.split_codingblock_on_alternatives_in_pacbps_dict(
            filter_for_msr=True,
            filter_for_omsr=True,
            )
    # remove non-compatible CBGs
    cbgList.remove_incompatible_cbgs(
            minimal_node_count=len(input),
            minimal_edge_count=len(tinyexon_crossdata),
            filter_for_msr=True,
            filter_for_omsr=True
            )

    # get list of accepted TinyExonCbgs 
    accepted_tegs = cbgList.codingblockgraphs

    # and update weights by minimal spanning region
    for teg in accepted_tegs: teg.update_edge_weights_by_minimal_spanning_range()

    # and check if they can be placed IN BETWEEN graphL and graphR
    # TODO some prints
    final_graphs_with_tinyexons = []
    for teg in accepted_tegs:

        test_codingblock_order, rejected_graphs = make_consensus_genestructure_from_compatible_pacb_graphs(
                [graphL,graphR,teg],None)

        print "checking hypo TEG:", teg.get_ordered_nodes(), "of", len(accepted_tegs), "len of join", len(test_codingblock_order)


        #empty_input = {}
        #for org in teg.organism_set(): empty_input[org] = None
        #tmpGSG = GenestructureOfCodingBlockGraphs(empty_input)
        #tmpGSG.add_codingblocks([graphL,graphR,teg])
        #print "tinyexon tmp check:", len(test_codingblock_order), len(tmpGSG), teg.get_ordered_nodes()


        wt_after  = teg.total_weight()
        if len(test_codingblock_order) == 3:
            teg_nodes = teg.get_nodes()
            teg_nodes.sort()
            middle = test_codingblock_order[1]
            middle_nodes = middle.get_nodes()
            middle_nodes.sort()
            if middle_nodes == teg_nodes:
                # yahoo, this one is 100% okay!
                final_graphs_with_tinyexons.append( teg )

    if len(final_graphs_with_tinyexons)==1:
        print final_graphs_with_tinyexons[0].get_nodes()
        return final_graphs_with_tinyexons
    elif len(final_graphs_with_tinyexons)>1:
        print "### WARNING!!!! more than 1 tinyexon graph is found."
        print "### WARNING!!!! however, only single one is returned."
        print "### WARNING!!!! returning >1 can cause errors..."
        return [ final_graphs_with_tinyexons[0] ]
    else:
        return []

# end of function find_intermediary_codingblockgraph_with_tinyexon


def harvest_intermediary_tinyexon_subgraphs(
        genestructure_graphs,input,
        similaritymatrix=None):
    """
    """
    finding_new = True
    startat = 0
    while finding_new:
        for pos in range(startat,len(genestructure_graphs)):
            sg = genestructure_graphs[pos]
            prev = None
            next = None
            if pos >= 1: prev = genestructure_graphs[pos-1]
            if pos < len(genestructure_graphs)-1: next = genestructure_graphs[pos+1]

            # make DonorSiteCollectionGraph and AcceptorSiteCollectionGraph, NO projected sites included yet!!!
            sg.harvest_elegiable_donor_sites(projected_donors={},next=next)
            sg.harvest_elegiable_acceptor_sites(projected_acceptors={},prev=prev)

            if prev:
                tinyexonsg = find_intermediary_codingblockgraph_with_tinyexon(
                        prev,sg,input=input,
                        similaritymatrix=similaritymatrix)
                if tinyexonsg:
                    tinysg = tinyexonsg[0]
                    # make DonorSiteCollectionGraph and AcceptorSiteCollectionGraph, NO projected sites included yet!!!
                    tinysg.harvest_elegiable_donor_sites(projected_donors={},next=sg)
                    tinysg.harvest_elegiable_acceptor_sites(projected_acceptors={},prev=prev)

                    # important! do not forget to rescan ``prev`` and ``sg``
                    # for splice site signale; it has been performed
                    # on hindsight woth the wrong next sg!
                    prev.harvest_elegiable_donor_sites(projected_donors={},next=tinysg)
                    sg.harvest_elegiable_acceptor_sites(projected_acceptors={},prev=tinysg)

                    # insert into genestructure graphs at the correct position
                    genestructure_graphs.insert( pos, tinysg )

                    # set startat pointer to current point in the list genestructure_graphs
                    # all previous sg's are scanned for tinyexons
                    startat = pos+1
                    break
        else:
            finding_new = False

    # return the (updated) list with subgraphs)
    return genestructure_graphs

# end of function harvest_intermediary_tinyexon_subgraphs



def _order_similarities(sims):
    """
    """
    data = []
    for ( ( ratio, sbjct_pos, q_seq, match, s_seq, bitscore), tinyexon ) in sims:
        score = float(ratio) * float(bitscore)
        data.append( ( -score, ( ratio, sbjct_pos, q_seq, match, s_seq, bitscore), tinyexon ) )
    data.sort()
    return [ ( b, c ) for (a,b,c) in data ]

# end of function _order_similarities

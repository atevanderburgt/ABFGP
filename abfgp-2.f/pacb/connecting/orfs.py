"""
Functions for connecting PacPORFs into a genemodel 
"""
# Module metadata
__authors__ = "Ate van der Burgt"
__license__ = "MIT"


# Imports
from pacb.validators import IsOrf

# Gene structure imports
from gene.intron import IntronConnectingOrfs, _order_intron_list
from gene.splicesite import get_shared_nucleotides_at_splicesite
from gene.exon import ExonOnOrf, FirstExonOnOrf
from gene.gene_exceptions import UnexpectedSpliceSitePhase, InproperlyAppliedArgument
from pacb.connecting.functions import _update_kwargs

# Global variable Imports
from settings.splicesites import (
    MIN_DONOR_PSSM_SCORE,
    ALLOW_NON_CANONICAL_DONOR,
    NON_CANONICAL_MIN_DONOR_PSSM_SCORE,
    MIN_ACCEPTOR_PSSM_SCORE,
    ALLOW_NON_CANONICAL_ACCEPTOR,
    NON_CANONICAL_MIN_ACCEPTOR_PSSM_SCORE,
    )
from settings.genestructure import (
    MIN_INTRON_NT_LENGTH,
    MAX_INTRON_NT_LENGTH,
    OPTIMAL_BRACNHPOINT_TO_ACCEPTOR_DISTANCE,
    MAXIMAL_OPTIMAL_BRACNHPOINT_TO_ACCEPTOR_DISTANCE,
    )


# import all TINYEXON stuff...
from settings.genestructure import *
from settings.splicesites import (
    KWARGS_PROJECTED_TINYEXON,
    KWARGS_STOPLESS_3N_INTRONS,
    )


# import all SignalP first exon stuff
SIGNALP_FIRSTEXON_MIN_NT_LENGTH                      = 35
SIGNALP_FIRSTEXON_MAX_NT_LENGTH                      = 200
SIGNALP_FIRSTEXON_MIN_TSS_PSSM_SCORE                 = 5.5
SIGNALP_FIRSTEXON_MIN_DONOR_PSSM_SCORE               = 4.5
SIGNALP_FIRSTEXON_ALLOW_NON_CANONICAL_DONOR          = False
SIGNALP_FIRSTEXON_NON_CANONICAL_MIN_DONOR_PSSM_SCORE = 0.0
SIGNALP_FIRSTEXON_TCODE_IS_CODING                    = None
SIGNALP_FIRSTEXON_TCODE_IS_NONCODING                 = False
SIGNALP_FIRSTEXON_TCODE_MIN_SCORE                    = False


def merge_orfs_with_intron(orfD,orfA,
    max_intron_nt_length         = MAX_INTRON_NT_LENGTH,
    min_intron_nt_length         = MIN_INTRON_NT_LENGTH,
    min_donor_pssm_score         = MIN_DONOR_PSSM_SCORE,
    min_acceptor_pssm_score      = MIN_ACCEPTOR_PSSM_SCORE,
    allow_non_canonical_donor    = ALLOW_NON_CANONICAL_DONOR,
    allow_non_canonical_acceptor = ALLOW_NON_CANONICAL_ACCEPTOR,
    non_canonical_min_donor_pssm_score    = NON_CANONICAL_MIN_DONOR_PSSM_SCORE,
    non_canonical_min_acceptor_pssm_score = NON_CANONICAL_MIN_ACCEPTOR_PSSM_SCORE,
    min_donor_pos=None,
    max_donor_pos=None,
    min_acceptor_pos=None,
    max_acceptor_pos=None,
    order_by = 'length',**kwargs):
    """
    Merge 2 Orf objects by introns

    @attention: **kwargs can contain other (here) unnecessarily arguments

    @type  orfD: Orf object
    @param orfD: Orf object that has to deliver a PSSM donor object

    @type  orfA: Orf object
    @param orfA: Orf object that has to deliver a PSSM acceptor object

    @type  max_intron_nt_length: integer
    @param max_intron_nt_length: maximal length (nt) of the intron
    
    @type  min_intron_nt_length: integer
    @param min_intron_nt_length: minimal length (nt) of the intron

    @type  min_donor_pssm_score: float
    @param min_donor_pssm_score: minimal pssm score of donor splice site

    @type  min_acceptor_pssm_score: float
    @param min_acceptor_pssm_score: minimal pssm score of acceptor splice site

    @type  allow_non_canonical_donor: Boolean
    @param allow_non_canonical_donor: search for non-canonical donor sites too

    @type  allow_non_canonical_acceptor: Boolean
    @param allow_non_canonical_acceptor: search for non-canonical acceptor splice sites too

    @type  non_canonical_min_donor_pssm_score: float
    @param non_canonical_min_donor_pssm_score: minimal pssm score of non-canonical donor

    @type  non_canonical_min_acceptor_pssm_score: float
    @param non_canonical_min_acceptor_pssm_score: minimal pssm score of non-canonical acceptor 

    @rtype:  list
    @return: list with introns
    """
    # input validation
    IsOrf(orfD)
    IsOrf(orfA)

    # scan for splice sites (if not already done -> is checked in function)
    orfD.scan_orf_for_pssm_splice_sites(
            splicetype="donor",
            min_pssm_score=min_donor_pssm_score,
            allow_non_canonical=allow_non_canonical_donor,
            non_canonical_min_pssm_score=non_canonical_min_donor_pssm_score)
    orfA.scan_orf_for_pssm_splice_sites(
            splicetype="acceptor",
            min_pssm_score=min_acceptor_pssm_score,
            allow_non_canonical=allow_non_canonical_acceptor,
            non_canonical_min_pssm_score=non_canonical_min_acceptor_pssm_score)

    # return list with introns
    introns = []

    # most quickest scan possible: are there donors & acceptors?
    if orfD._donor_sites == [] or orfA._acceptor_sites == []:
        # no introns possible because splice sites are missing
        return introns

    # very quick scan: are exons not to far from each other?
    if max_intron_nt_length and\
    (orfA._acceptor_sites[0].pos - orfD._donor_sites[0].pos) > max_intron_nt_length:
        # no introns possible that can bridge this gap
        return introns

    for donor in orfD._donor_sites:
        if not allow_non_canonical_donor and not donor.is_canonical():
            continue
        elif donor.is_canonical() and donor.pssm_score < min_donor_pssm_score:
            continue
        elif not donor.is_canonical() and donor.pssm_score < non_canonical_min_donor_pssm_score:
            continue
        elif (min_donor_pos or min_donor_pos==0) and donor.pos < min_donor_pos:
            continue
        elif (max_donor_pos or max_donor_pos==0) and donor.pos > max_donor_pos:
            continue
        else:
            # donor site accepted
            pass 

        for acceptor in orfA._acceptor_sites:
            if not allow_non_canonical_acceptor and not acceptor.is_canonical():
                continue
            elif acceptor.is_canonical() and acceptor.pssm_score < min_acceptor_pssm_score:
                continue
            elif not acceptor.is_canonical() and acceptor.pssm_score < non_canonical_min_acceptor_pssm_score:
                continue
            elif (min_acceptor_pos or min_acceptor_pos==0) and acceptor.pos < min_acceptor_pos:
                continue
            elif (max_acceptor_pos or max_acceptor_pos==0) and acceptor.pos > max_acceptor_pos:
                continue
            else:
                # acceptor site accepted
                pass 

            # generate intron length and phase variable
            intron_length = acceptor.pos - donor.pos
            intron_phase  = intron_length % 3

            # check phase compatibilty (1) of splice sites
            if donor.phase != acceptor.phase: continue
            # check phase compatibilty (2) of splice sites
            if ( intron_phase + orfD.frame ) % 3 != orfA.frame % 3: continue

            # check if intron length is in between the boundaries
            if max_intron_nt_length and intron_length > max_intron_nt_length: continue
            if min_intron_nt_length and intron_length < min_intron_nt_length: continue

            # okay, if we reach this point, we have a valid intron
            shared_nts = get_shared_nucleotides_at_splicesite(
                    orfA,orfD,acceptor,donor
                    )

            # make a IntronConnectingOrfs object
            intron = IntronConnectingOrfs(donor,acceptor,shared_nts,orfD,orfA)
            introns.append(intron)

    # return ordered intron list
    return _order_intron_list(introns,order_by=order_by)

# end of function merge_orfs_with_intron


def merge_orfs_with_tinyexon(preceding_orf,subsequent_orf,
    preceding_donor_sites=[],
    subsequent_acceptor_sites=[],
    orflist=[],
    max_tinyexon_nt_length=TINYEXON_MAX_NT_LENGTH,
    min_tinyexon_nt_length=TINYEXON_MIN_NT_LENGTH,
    max_tinyexon_intron_nt_length=TINYEXON_MAX_INTRON_NT_LENGTH,
    min_tinyexon_intron_nt_length=TINYEXON_MIN_INTRON_NT_LENGTH,
    min_donor_pssm_score=TINYEXON_MIN_DONOR_PSSM_SCORE,
    min_acceptor_pssm_score=TINYEXON_MIN_ACCEPTOR_PSSM_SCORE,
    min_total_pssm_score=TINYEXON_MIN_TOTAL_PSSM_SCORE,
    **kwargs):
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

    # return list with (intron,tinyexon,intron) tuples
    returnexons = []
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
                    shared_nts_A = get_shared_nucleotides_at_splicesite(
                            tinyexon.orf,preceding_orf,
                            tinyexon.acceptor,donor
                            )
                    preceding_intron = IntronConnectingOrfs(
                        donor,tinyexon.acceptor,
                        shared_nts_A,preceding_orf,tinyexon.orf )

                    # make subsequent intron
                    shared_nts_B = get_shared_nucleotides_at_splicesite(
                            subsequent_orf,tinyexon.orf,
                            acceptor,tinyexon.donor
                            )

                    subsequent_intron = IntronConnectingOrfs(
                        tinyexon.donor, acceptor,
                        shared_nts_B,tinyexon.orf,subsequent_orf )

                    # and append to exons
                    returnexons.append( ( preceding_intron, tinyexon, subsequent_intron ) )

    # and return the list of intron/exon/intron
    return returnexons

# end of function merge_orfs_with_tinyexon


def find_tiny_exon_on_orf(orfX,order_by='total_pssm',
    max_tinyexon_nt_length       =TINYEXON_MAX_NT_LENGTH,
    min_tinyexon_nt_length       =TINYEXON_MIN_NT_LENGTH,
    max_tinyexon_intron_nt_length=TINYEXON_MAX_INTRON_NT_LENGTH,
    min_tinyexon_intron_nt_length=TINYEXON_MIN_INTRON_NT_LENGTH,
    min_donor_pssm_score         =TINYEXON_MIN_DONOR_PSSM_SCORE,
    min_acceptor_pssm_score      =TINYEXON_MIN_ACCEPTOR_PSSM_SCORE,
    min_total_pssm_score         =TINYEXON_MIN_TOTAL_PSSM_SCORE,
    allow_non_canonical_donor    = TINYEXON_ALLOW_NON_CANONICAL_DONOR,
    allow_non_canonical_acceptor = TINYEXON_ALLOW_NON_CANONICAL_ACCEPTOR,
    non_canonical_min_donor_pssm_score    = TINYEXON_NON_CANONICAL_MIN_DONOR_PSSM_SCORE,
    non_canonical_min_acceptor_pssm_score = TINYEXON_NON_CANONICAL_MIN_ACCEPTOR_PSSM_SCORE,
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

    @type  allow_non_canonical_donor: Boolean
    @param allow_non_canonical_donor: search for non-canonical donor sites too

    @type  allow_non_canonical_acceptor: Boolean
    @param allow_non_canonical_acceptor: search for non-canonical acceptor splice sites too

    @type  non_canonical_min_donor_pssm_score: float
    @param non_canonical_min_donor_pssm_score: minimal pssm score of non-canonical donor

    @type  non_canonical_min_acceptor_pssm_score: float
    @param non_canonical_min_acceptor_pssm_score: minimal pssm score of non-canonical acceptor 

    @type  order_by: TODO
	@param order_by: TODO
    """
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
                raise "WUF... WUF..."
        except:
            message = "a variable is NOT a positive integer as expected"
            raise InproperlyAppliedArgument, message

    # scan for splice sites on this (tiny) orf
    orfX.scan_orf_for_pssm_splice_sites(
            splicetype="donor",
            min_pssm_score=min_donor_pssm_score,
            allow_non_canonical=allow_non_canonical_donor,
            non_canonical_min_pssm_score=non_canonical_min_donor_pssm_score)
    orfX.scan_orf_for_pssm_splice_sites(
            splicetype="acceptor",
            min_pssm_score=min_acceptor_pssm_score,
            allow_non_canonical=allow_non_canonical_acceptor,
            non_canonical_min_pssm_score=non_canonical_min_acceptor_pssm_score)

    # return list with exons
    exons = []

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


def get_potention_first_exons_on_orf(orfX,order_by='length',
    max_tinyexon_nt_length      = TINYEXON_MAX_NT_LENGTH,
    min_tinyexon_nt_length      = 2,
    min_tss_pssm_score          = 4.0,
    min_donor_pssm_score        = TINYEXON_MIN_DONOR_PSSM_SCORE,
    allow_non_canonical_donor   = TINYEXON_ALLOW_NON_CANONICAL_DONOR,
    non_canonical_min_donor_pssm_score = TINYEXON_NON_CANONICAL_MIN_DONOR_PSSM_SCORE,
    min_total_pssm_score=TINYEXON_MIN_TOTAL_PSSM_SCORE,**kwargs):
    """
    Predict all possible tiny leading exons on this Orf

    @attention: **kwargs can contain other (here) unnecessarily arguments
    """
    # scan for donor sites on this (tiny) orf
    orfX.scan_orf_for_pssm_splice_sites(
            splicetype="donor",
            min_pssm_score=min_donor_pssm_score,
            allow_non_canonical=allow_non_canonical_donor,
            non_canonical_min_pssm_score=non_canonical_min_donor_pssm_score)

    # scan for TSS sites on this (tiny) orf
    orfX.scan_orf_for_pssm_tss(min_pssm_score=min_tss_pssm_score)

    # return list with leading exons
    leadingexons = []

    # loop over the TSS sites
    for tss in orfX._tss_sites:
        if (min_tss_pssm_score or min_tss_pssm_score==0.0) and\
        tss.pssm_score < min_tss_pssm_score:
            continue
        # loop over the donor sites
        for donor in orfX._donor_sites:
            if donor.pos < tss.pos:
                continue
            if (min_donor_pssm_score or min_donor_pssm_score == 0.0) and\
            donor.pssm_score < min_donor_pssm_score:
                # splicesite doesn't have a high enough pssm_score
                continue
            if (donor.pos - tss.pos) < min_tinyexon_nt_length:
                # short atg-exon to short
                continue
            if (donor.pos - tss.pos) > max_tinyexon_nt_length:
                # short atg-exon to long
                continue

            # make LeadingExonOnOrf object
            leadingexons.append( FirstExonOnOrf(tss,donor,orfX) )

    # order by summed pssm score
    tmp = [ (exon.acceptor.pssm_score+exon.donor.pssm_score,exon) for exon in leadingexons ]
    tmp.sort()
    tmp.reverse()
    leadingexons = [ exon for score,exon in tmp ]

    # return the found exons
    return leadingexons

# end of function get_potention_first_exons_on_orf


def get_signalp_first_exon_on_orf(orfX,order_by='length',
    max_tinyexon_nt_length      = SIGNALP_FIRSTEXON_MAX_NT_LENGTH,
    min_tinyexon_nt_length      = SIGNALP_FIRSTEXON_MIN_NT_LENGTH,
    min_tss_pssm_score          = SIGNALP_FIRSTEXON_MIN_TSS_PSSM_SCORE,
    min_donor_pssm_score        = SIGNALP_FIRSTEXON_MIN_DONOR_PSSM_SCORE,
    tcode_is_coding             = SIGNALP_FIRSTEXON_TCODE_IS_CODING,
    tcode_is_noncoding          = SIGNALP_FIRSTEXON_TCODE_IS_NONCODING,
    tcode_min_score             = SIGNALP_FIRSTEXON_TCODE_MIN_SCORE, 
    allow_non_canonical_donor   = SIGNALP_FIRSTEXON_ALLOW_NON_CANONICAL_DONOR,
    non_canonical_min_donor_pssm_score =\
    SIGNALP_FIRSTEXON_NON_CANONICAL_MIN_DONOR_PSSM_SCORE,**kwargs):
    """
    Predict the single most likely first exon with SignalPeptide on this Orf

    @attention: **kwargs can contain other (here) unnecessarily arguments
    """
    # scan for donor sites on this (tiny) orf
    orfX.scan_orf_for_pssm_splice_sites(
            splicetype="donor",
            min_pssm_score=min_donor_pssm_score,
            allow_non_canonical=allow_non_canonical_donor,
            non_canonical_min_pssm_score=non_canonical_min_donor_pssm_score)

    # scan for TSS sites on this (tiny) orf
    orfX.scan_orf_for_pssm_tss(min_pssm_score=min_tss_pssm_score)

    # scan Orf for SignalPeptide(s); this has been done in the ABGP genelocusdir.
    # so, only check if signalp's are present or not on this orf
    if orfX._has_signalp_sites_predicted:
        if orfX._signalp_sites:
            pass
        else:
            return None # no SignalPetides present!
    else:
        # here, make code that does SignalP predictions. For now, omit
        return None

    # return variable (temp. a list) with SignalPexon
    signalpexons = []

    # loop over the TSS sites
    for tss in orfX._tss_sites:
        if tss.pos not in [ sP.tss.pos for sP in orfX._signalp_sites]:
            continue
        if (min_tss_pssm_score or min_tss_pssm_score==0.0) and\
        tss.pssm_score < min_tss_pssm_score:
            continue
        # loop over the donor sites
        for donor in orfX._donor_sites:
            if donor.pos < tss.pos:
                continue
            if (min_donor_pssm_score or min_donor_pssm_score == 0.0) and\
            donor.pssm_score < min_donor_pssm_score:
                # splicesite doesn't have a high enough pssm_score
                continue
            if (donor.pos - tss.pos) < min_tinyexon_nt_length:
                # short atg-exon to short
                continue
            if (donor.pos - tss.pos) > max_tinyexon_nt_length:
                # short atg-exon to long
                continue

            # make LeadingExonOnOrf object
            firstExonObj = FirstExonOnOrf(tss,donor,orfX)

            # filter on Tcode properties
            if tcode_is_coding and not firstExonObj.is_tcode_coding():
                continue
            if tcode_is_noncoding == False and firstExonObj.is_tcode_noncoding():
                continue
            if tcode_min_score and firstExonObj.tcode() < tcode_min_score:
                continue

            # if here, append to signalpexon list
            signalpexons.append( firstExonObj )

    if not signalpexons:
        return None
    elif len(signalpexons) == 1:
        return signalpexons[0]
    else:
        # order by tss position -> take most frontal
        tmp = [ (exon.acceptor.pos,exon) for exon in signalpexons ]
        tmp.sort()
        # remove all signalp exons that do not start with most frontal TSS
        while tmp[-1][0] != tmp[0][0]: tmp.pop()
        # second, order by donor position -> take longest one
        tmp = [ (exon.donor.pos,exon) for (pos,exon) in tmp ]
        tmp.sort()
        tmp.reverse()
        return tmp[0][1]

# end of function get_signalp_first_exon_on_orf


def get_potential_tiny_exons_on_orf(orfX,order_by='length',
    max_tinyexon_nt_length      =TINYEXON_MAX_NT_LENGTH,
    min_tinyexon_nt_length      =TINYEXON_MIN_NT_LENGTH,
    min_donor_pssm_score        =TINYEXON_MIN_DONOR_PSSM_SCORE,
    min_acceptor_pssm_score     =TINYEXON_MIN_ACCEPTOR_PSSM_SCORE,
    allow_non_canonical_donor   = TINYEXON_ALLOW_NON_CANONICAL_DONOR,
    allow_non_canonical_acceptor= TINYEXON_ALLOW_NON_CANONICAL_ACCEPTOR,
    non_canonical_min_donor_pssm_score    = TINYEXON_NON_CANONICAL_MIN_DONOR_PSSM_SCORE,
    non_canonical_min_acceptor_pssm_score = TINYEXON_NON_CANONICAL_MIN_ACCEPTOR_PSSM_SCORE,
    min_total_pssm_score=TINYEXON_MIN_TOTAL_PSSM_SCORE,**kwargs):
    """
    Predict all possible tiny exons on this Orf

    @attention: **kwargs can contain other (here) unnecessarily arguments

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

    @type  allow_non_canonical_donor: Boolean
    @param allow_non_canonical_donor: search for non-canonical donor sites too

    @type  allow_non_canonical_acceptor: Boolean
    @param allow_non_canonical_acceptor: search for non-canonical acceptor splice sites too

    @type  non_canonical_min_donor_pssm_score: float
    @param non_canonical_min_donor_pssm_score: minimal pssm score of non-canonical donor

    @type  non_canonical_min_acceptor_pssm_score: float
    @param non_canonical_min_acceptor_pssm_score: minimal pssm score of non-canonical acceptor 

    @type  order_by: TODO
	@param order_by: TODO

    @rtype:  list
    @return: list with tinyexons
    """

    # scan for splice sites on this (tiny) orf
    forced = orfX._donor_sites == []
    orfX.scan_orf_for_pssm_splice_sites(
            splicetype="donor",forced=forced,
            min_pssm_score=min_donor_pssm_score,
            allow_non_canonical=allow_non_canonical_donor,
            non_canonical_min_pssm_score=non_canonical_min_donor_pssm_score)
    forced = orfX._acceptor_sites == []
    orfX.scan_orf_for_pssm_splice_sites(
            splicetype="acceptor",forced=forced,
            min_pssm_score=min_acceptor_pssm_score,
            allow_non_canonical=allow_non_canonical_acceptor,
            non_canonical_min_pssm_score=non_canonical_min_acceptor_pssm_score)

    # return list with exons
    tinyexons = []

    # most quickest scan possible: are there donors & acceptors?
    if orfX._donor_sites == [] or orfX._acceptor_sites == []:
        # no exons possible because splice sites are missing
        return tinyexons

    # and combine sites to exons!
    for acceptor in orfX._acceptor_sites:
        for donor in orfX._donor_sites:
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

            # (re) check individual PSSM scores; in case this Orf had
            # already pre-assigned slice sites, potentially stricter
            # parameters are not applied!
            if donor.is_canonical() and\
            (min_donor_pssm_score or min_donor_pssm_score == 0.0) and\
            donor.pssm_score < min_donor_pssm_score:
                continue
            if not donor.is_canonical() and\
            (non_canonical_min_donor_pssm_score or non_canonical_min_donor_pssm_score == 0.0) and\
            donor.pssm_score < non_canonical_min_donor_pssm_score:
                continue
            if (min_acceptor_pssm_score or min_acceptor_pssm_score == 0.0) and\
            acceptor.pssm_score < min_acceptor_pssm_score:
                continue

            # make a Exon object
            exon = ExonOnOrf(acceptor,donor,orfX)
            tinyexons.append(exon)

    # return ordered exon list
    return _order_intron_list(tinyexons,order_by=order_by)

# end of function get_potential_tiny_exons_on_orf


def merge_orfs_with_two_tinyexons(preceding_orf,subsequent_orf,
    preceding_donor_sites=[],
    subsequent_acceptor_sites=[],
    orflist=[],**kwargs):
    """
    Bridge two `neighbouring` Orfs by TWO tinyexon by applying preceding donors and subsequent acceptors

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

    @attention: see get_potential_tiny_exons_on_orf for additional **kwargs

    @rtype:  list
	@return: list of tuples ( preceding_intron, tinyexon1, central_intron, tinyexon2, subsequent_intron )

    """
    if not preceding_donor_sites:
        return []
    if not subsequent_acceptor_sites:
        return []
    if not orflist:
        return []

    # edit **kwargs dictionary for some forced attributes
    _update_kwargs(kwargs,KWARGS_PROJECTED_TINYEXON)

    # return list with (intron,tinyexon,intron) tuples
    returntinyexons = []
    tinyexoncollection = []
    tinyexoncombis = []
    min_preceding_donor_sites_pos     = min([ d.pos for d in preceding_donor_sites ])
    max_subsequent_acceptor_sites_pos = max([ a.pos for a in subsequent_acceptor_sites ]) 

    for orfX in orflist:
        # check if orf is correctly positions towards the splice sites' extremes
        min_pos = min_preceding_donor_sites_pos + kwargs['min_tinyexon_intron_nt_length']
        max_pos = max_subsequent_acceptor_sites_pos - kwargs['min_tinyexon_intron_nt_length']
        # if so, do not check this Orf
        if orfX.endPY   <= min_pos: continue
        if orfX.startPY >= max_pos: continue
        # extend the tinyexoncollection
        tinyexoncollection.extend( get_potential_tiny_exons_on_orf(orfX,**kwargs) )

    # make tinyexoncollection ordered on start pos
    tinyexoncollection = _order_intron_list(tinyexoncollection,order_by='donor_pos')
    # donor_pos makes REVERSE ordering; restore this by reversing
    tinyexoncollection.reverse()

    # make 2-elemented tuples of tinyexons which can co-occur together
    for tinyexon1 in tinyexoncollection:
        for pos in range(len(tinyexoncollection)-1,-1,-1):
            tinyexon2 = tinyexoncollection[pos]
            if tinyexon2.donor.pos < tinyexon1.donor.pos: break
            intron_length = tinyexon2.acceptor.pos - tinyexon1.donor.pos
            if intron_length < kwargs['min_tinyexon_intron_nt_length']: continue
            if intron_length > kwargs['max_tinyexon_intron_nt_length']: continue
            if tinyexon1.donor.phase != tinyexon2.acceptor.phase: continue
            # if here, elegiable combi!
            intron = IntronConnectingOrfs(
                    tinyexon1.donor,tinyexon2.acceptor,
                    get_shared_nucleotides_at_splicesite(
                            subsequent_orf,preceding_orf,
                            tinyexon2.acceptor,tinyexon1.donor
                            ),
                    preceding_orf,subsequent_orf)
            totlen = tinyexon1.length+tinyexon2.length
            combi = ( totlen, tinyexon1, intron, tinyexon2 )
            tinyexoncombis.append( combi )

    # return an ordered list based on length
    tinyexoncombis.sort()
    return [ (exon1,intron,exon2) for l,exon1,intron,exon2 in tinyexoncombis ]

# end of function merge_orfs_with_two_tinyexons


def find_stopless3n_introns_on_orf(orfObj,
    has_branchpoint = False,
    has_polypyrimidine = False,
    order_by = 'length',**kwargs):
    """
    Find potential stopless3n introns on this orf

    @attention: **kwargs can contain other (here) unnecessarily arguments
    @attention: **kwargs are required in the merge_orfs_with_intron() function

    @type  orfObj: Orf object
    @param orfObj: Orf object which is scanned for stopless3n introns

    @rtype:  list
    @return: list with introns
    """
    # input validation
    IsOrf(orfObj)

    # edit **kwargs dictionary for some forced attributes
    _update_kwargs(kwargs,KWARGS_STOPLESS_3N_INTRONS)

    # find stopless3nintrons
    stopless3nintrons = merge_orfs_with_intron(orfObj,orfObj,**kwargs)

    # filter for presence of branchpoint / polypyrimidine tracks
    if has_branchpoint or has_polypyrimidine:
        filtered = []
        for intron in stopless3nintrons:
            intron.assign_bp_and_ppts()
            if has_branchpoint and not intron.branchpoint:
                continue
            intron_bp_dist = intron.get_branchpoint_nt_distance()
            if has_branchpoint and intron_bp_dist == None:
                continue
            intron_bp_optimality = min([ abs(offset-intron_bp_dist) for offset in OPTIMAL_BRACNHPOINT_TO_ACCEPTOR_DISTANCE ])
            if has_branchpoint and intron_bp_optimality > MAXIMAL_OPTIMAL_BRACNHPOINT_TO_ACCEPTOR_DISTANCE:
                continue 
            if has_polypyrimidine and not (intron.ppt5p or intron.ppt3p):
                continue
            # if here, accepted!
            filtered.append( intron )
    else:
        filtered = stopless3nintrons

    # return ordered intron list
    return _order_intron_list(filtered,order_by=order_by)

# end of function find_stopless3n_introns_on_orf

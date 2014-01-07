"""
PacbPORF connection by mapping introns on a more complex way
"""

# Module metadata
__authors__ = "Ate van der Burgt"
__license__ = "MIT"

# Imports
from pacb.validators import IsPacbPORF
from pacb.ordering import order_list_by_attribute as olba
from pacb.connecting.orfs import merge_orfs_with_intron, merge_orfs_with_tinyexon
from pacb.conversion import pacbp2pacbporf
from pacb.exceptions import CoordinateOutOfRange
from pacb.connecting.functions import (
    _tinyexon_list_2_dict,
    _update_kwargs,
    _filter_for_entropy,
    _filter_for_alignable_splice_sites,
    _score_introns_obtained_by_mapping,
    set_apps_intron_query,
    set_apps_intron_sbjct,
    )

# Other Imports

# Python Imports
from sets import Set

# Global variable Imports
from settings.splicesites import (
    KWARGS_MAPPED_INTRON,
    KWARGS_MAPPED_INTRON_STRACC,
    KWARGS_MAPPED_INTRON_STRDON,
)


def merge_pacbporfs_by_query_tinyexon_and_sbjct_intron(pacbporfD,pacbporfA,
    orfSetObjQ,verbose=False,**kwargs):
    """ """
    # input validation
    IsPacbPORF(pacbporfD)
    IsPacbPORF(pacbporfA)

    orfSetObjS = []

    # edit **kwargs dictionary for some forced attributes
    _update_kwargs(kwargs,KWARGS_MAPPED_INTRON)
    if not kwargs.has_key('aligned_site_max_triplet_distance'):
        kwargs['aligned_site_max_triplet_distance'] = kwargs['max_aa_offset']

    # settings for minimal alignment entropy score
    # TODO TODO -> THIS MUST BE FIXED TO A NICE THRESHOLD VALUE!!!
    min_donor_site_alignment_entropy = 0.1
    min_acceptor_site_alignment_entropy = 0.1

    resultlistQ = merge_orfs_with_tinyexon(
            pacbporfD.orfQ,pacbporfA.orfQ,
            preceding_donor_sites=pacbporfD.orfQ._donor_sites,
            subsequent_acceptor_sites=pacbporfA.orfQ._acceptor_sites,
            orflist=orfSetObjQ.orfs,**kwargs)
    resultlistS = merge_orfs_with_intron(
            pacbporfD.orfS,pacbporfA.orfS,**kwargs)

    # translate resultlists to dict: key == exon, value = [ {intronsD},{intronsS} ]
    resultdictQ,key2exonQ = _tinyexon_list_2_dict(resultlistQ)

    # get unique list of donors & acceptors
    donorQ = olba( list(Set([inD.donor for inD,te,inA in resultlistQ ])), order_by='pos', reversed=True)
    accepQ = olba( list(Set([inA.acceptor for inD,te,inA in resultlistQ ])), order_by='pos')

    donorS = olba( list(Set([intron.donor for intron in resultlistS ])), order_by='pos', reversed=True)
    accepS = olba( list(Set([intron.acceptor for intron in resultlistS ])), order_by='pos')

    ## filter for alignable donor & acceptor sites
    kwargs['allow_non_canonical']               = True # True
    kwargs['aligned_site_max_triplet_distance'] = 0     # 2
    algdonors = _filter_for_alignable_splice_sites(donorQ,donorS,pacbporfD,**kwargs)
    algacceps = _filter_for_alignable_splice_sites(accepQ,accepS,pacbporfA,**kwargs)

    # remove sites at the wrong interface of the pabcporfs
    # required to avoid CoordinateOutOfrange exceptions lateron
    algdonors = _filter_for_entropy(algdonors,pacbporfD,'donor',
                min_alignment_entropy=0.0)
    algacceps = _filter_for_entropy(algacceps,pacbporfA,'acceptor',
                min_alignment_entropy=0.0)


    # remove sites with to low alignment entropy
    algdonors = _filter_for_entropy(algdonors,pacbporfD,'donor',
                min_alignment_entropy=min_donor_site_alignment_entropy)
    algacceps = _filter_for_entropy(algacceps,pacbporfA,'acceptor',
                min_alignment_entropy=min_acceptor_site_alignment_entropy)


    results = []

    ## search for positive cases with aligned donors
    ##  ========.....=..=========
    ##  ========.....==========
    for algDq,algDs in algdonors:
        for (intronDQ,tinyexon,intronAQ) in resultlistQ:
            if intronDQ.donor.pos != algDq.pos: continue
            for intronS in resultlistS:
                if intronS.donor.pos != algDs.pos:
                    # not current aligned site
                    continue
                if (intronS.acceptor.pos / 3) < pacbporfA.sbjct_start:
                    # mapping 5' outside of extended PacbP -> continue
                    continue
                if (intronS.acceptor.pos + tinyexon.length) > ((pacbporfA.sbjct_end*3)+pacbporfA.orfS.frame):
                    # mapping 3' outside of extended PacbP -> continue
                    continue

                # calculate distance between the acceptor sites;
                # this distance should match tinyexon.length
                try:
                    distAnt = pacbporfA.get_distance_aligned_nucleotide_positions(
                                query = intronAQ.acceptor.pos,
                                sbjct = intronS.acceptor.pos
                                )
                except CoordinateOutOfRange:
                    # not succesfull in obtaining distance -> no proper mapping
                    continue

                if distAnt == tinyexon.length:
                    # get sbjct sequence for comparison
                    # OLD code; gave pacb.exceptions.CoordinateOutOfRange
                    # errors when (by chance) the sbjct Orf was not long enough
                    # extended in the PacbPORF. Code is now replaced by
                    # literal coordinate retrieval on the sbjct Orf object
                    ###startPos,_phase = pacbporfA.dnaposition_sbjct(intronS.acceptor.pos)
                    ###stopPos,_phase = pacbporfA.dnaposition_sbjct(intronS.acceptor.pos + tinyexon.length)
                    ###start = pacbporfA._positions[startPos].sbjct_pos
                    ###stop  = pacbporfA._positions[stopPos].sbjct_pos
                    start = intronS.acceptor.pos / 3
                    stop  = (intronS.acceptor.pos + tinyexon.length - pacbporfA.orfS.frame) / 3

                    if intronS.phase != 0:
                        # correct for first AA which is broken by the splicesite
                        start+=1
                    if stop <= start:
                        # tinyexon is so tiny that is does not have a single
                        # full aligned AA -> discard here
                        continue

                    # actually get the sbjct sequence
                    # TODO: get sbjctseq directly from ORF !?
                    queryseq,match,sbjctseq,coords = pacbporfA.alignmentpart_by_sbjct(start,stop)
                    sbjctseq = sbjctseq.replace("-","")
                    sbjct_start = coords[2]

                    ############################################################
                    if verbose:
                        print "PACBPseqs:", tinyexon.proteinsequence(), sbjctseq,
                        print (tinyexon.acceptor.phase, tinyexon.donor.phase, tinyexon.length),
                        print (intronS.phase, start, stop)
                    ############################################################

                    if len(tinyexon.proteinsequence()) != len(sbjctseq):
                        # IMPORTANT!! this should not be possible here (but it occurs)
                        # TODO: solve this. For now, continue by breaking out
                        continue

                    # initialize extended tinyexon PacbPORF
                    from pacb import PacbP
                    pacbp = PacbP(input=(
                            tinyexon.proteinsequence(),
                            sbjctseq,
                            (tinyexon.acceptor.pos / 3),
                            sbjct_start,
                            ) )

                    # continue if bitscore is <0
                    if pacbp.bitscore < 0: continue
                    pacbp.strip_unmatched_ends()
                    # continue if no fraction could be aligned
                    if len(pacbp) == 0: continue
                    tinypacbporf = pacbp2pacbporf(pacbp,tinyexon.orf,pacbporfA.orfS)
                    tinypacbporf.extend_pacbporf_after_stops()

                    ############################################################
                    if verbose:
                        print "#### ALIGNED DONORS",
                        print ( intronS.acceptor.pos, intronS.acceptor.pos + tinyexon.length )
                        print "####", pacbporfD
                        print "####", tinyexon, tinyexon.proteinsequence(),sbjctseq
                        print "####", pacbporfA
                        print "####", tinypacbporf
                        tinypacbporf.print_protein_and_dna()
                        print intronDQ
                        print intronAQ
                        print intronS
                        print ""
                    ############################################################

                    ############################################################
                    # set some meta-data properties to the intron objects
                    ############################################################
                    _score_introns_obtained_by_mapping(
                            intronDQ,intronS,pacbporfD,
                            tinypacbporf,source='ABGPmappingTE')
                    # set distance -> 0nt
                    intronAQ._distance = 0
                    # add Alignment Positional Periphery Score into objects
                    succes = set_apps_intron_query(intronAQ,tinypacbporf,pacbporfA)
                    # set GFF fsource attribute for recognition of intron sources
                    intronAQ._gff['fsource'] = "ABGPprojectingTE"

                    # create _linked_to_xxx attributes
                    intronDQ._linked_to_pacbporfs = [ tinypacbporf ]
                    intronAQ._linked_to_pacbporfs = [ tinypacbporf ]
                    intronS._linked_to_pacbporfs = [ tinypacbporf ]
                    intronDQ._linked_to_introns   = [ intronAQ ]
                    intronAQ._linked_to_introns   = [ intronDQ ]

                    # append to results
                    results.append((intronDQ,intronS,tinypacbporf,intronAQ,None))


    ## search for positive cases with aligned acceptors
    ##  ========..=.....==========
    ##    =========.....==========
    for algAq,algAs in algacceps:
        for (intronDQ,tinyexon,intronAQ) in resultlistQ:
            if intronAQ.acceptor.pos != algAq.pos: continue
            for intronS in resultlistS:
                if intronS.acceptor.pos != algAs.pos:
                    # not current aligned site
                    continue
                if (intronS.donor.pos / 3) > pacbporfD.sbjct_end:
                    # mapping 3' outside of extended PacbP -> continue
                    continue
                if (intronS.donor.pos - tinyexon.length) <= ((pacbporfD.sbjct_start*3)+pacbporfD.orfS.frame):
                    # mapping 5' outside of extended PacbP -> continue
                    continue
                if intronDQ.donor.pos < pacbporfD._get_original_alignment_pos_start().query_dna_start:
                    # mapping 5' outside of PacbP -> continue
                    continue

                # calculate distance between the donor sites;
                # this distance should match tinyexon.length
                try:
                    distDnt = pacbporfD.get_distance_aligned_nucleotide_positions(
                                query = intronDQ.donor.pos,
                                sbjct = intronS.donor.pos
                                )
                except CoordinateOutOfRange:
                    # not succesfull in obtaining distance -> no proper mapping
                    continue

                if distDnt == tinyexon.length:
                    # get sbjct sequence for comparison
                    # OLD code; gave pacb.exceptions.CoordinateOutOfRange
                    # errors when (by chance) the sbjct Orf was not long enough
                    # extended in the PacbPORF. Code is now replaced by
                    # literal coordinate retrieval on the sbjct Orf object
                    ###startPos,_phase = pacbporfD.dnaposition_sbjct(intronS.donor.pos - tinyexon.length, forced_return=True)
                    ###stopPos,_phase = pacbporfD.dnaposition_sbjct(intronS.donor.pos)
                    ###start = pacbporfD._positions[startPos].sbjct_pos
                    ###stop  = pacbporfD._positions[stopPos].sbjct_pos
                    ###print start, stop, intronS.phase, pacbporfD.orfS.frame, tinyexon.acceptor.phase
                    start = ( intronS.donor.pos - tinyexon.length - pacbporfD.orfS.frame) / 3
                    stop  = intronS.donor.pos / 3

                    if tinyexon.acceptor.phase != 0:
                        # correct for first AA which is broken by the splicesite
                        start+=1
                    if stop <= start:
                        # tinyexon is so tiny that is does not have a single
                        # full aligned AA -> discard here
                        continue

                    # actually get the sbjct sequence
                    # TODO: get sbjctseq directly from ORF !?
                    queryseq,match,sbjctseq,coords = pacbporfD.alignmentpart_by_sbjct(start,stop)
                    sbjctseq = sbjctseq.replace("-","")
                    sbjct_start = coords[2]

                    ############################################################
                    if verbose or len(tinyexon.proteinsequence()) != len(sbjctseq):
                        print "PACBPseqs:", tinyexon.proteinsequence(), sbjctseq,
                        print (tinyexon.acceptor.phase, tinyexon.donor.phase, tinyexon.length),
                        print (intronS.phase, start, stop)
                    ############################################################

                    if len(tinyexon.proteinsequence()) != len(sbjctseq):
                        # IMPORTANT!! this should not be possible here (but it occurs)
                        # TODO: solve this. For now, continue by breaking out
                        continue

                    # initialize extended tinyexon PacbPORF
                    from pacb import PacbP
                    pacbp = PacbP(input=(
                            tinyexon.proteinsequence(),
                            sbjctseq,
                            (tinyexon.acceptor.pos / 3),
                            sbjct_start,
                            ) )

                    # continue if bitscore is <0
                    if pacbp.bitscore < 0: continue
                    pacbp.strip_unmatched_ends()
                    # continue if no fraction could be aligned
                    if len(pacbp) == 0: continue
                    tinypacbporf = pacbp2pacbporf(pacbp,tinyexon.orf,pacbporfD.orfS)
                    tinypacbporf.extend_pacbporf_after_stops()

                    ############################################################
                    if verbose:
                        print "#### ALIGNED ACCEPTORS"
                        print "####", pacbporfD
                        print "####", tinyexon, tinyexon.proteinsequence(),sbjctseq
                        print "####", pacbporfA
                        print "####", tinypacbporf
                        tinypacbporf.print_protein_and_dna()
                        print intronDQ
                        print intronAQ
                        print intronS
                        print ""
                    ############################################################

                    ############################################################
                    # set some meta-data properties to the intron objects
                    ############################################################
                    _score_introns_obtained_by_mapping(
                            intronAQ,intronS,tinypacbporf,pacbporfA,
                            source='ABGPmappingTE')
                    # set distance -> 0nt
                    intronDQ._distance = 0
                    # add Alignment Positional Periphery Score into objects
                    succes = set_apps_intron_query(intronDQ,pacbporfD,tinypacbporf)
                    # set GFF fsource attribute for recognition of intron sources
                    intronDQ._gff['fsource'] = "ABGPprojectingTE"

                    # create _linked_to_xxx attributes
                    intronDQ._linked_to_pacbporfs = [ tinypacbporf ]
                    intronAQ._linked_to_pacbporfs = [ tinypacbporf ]
                    intronS._linked_to_pacbporfs = [ tinypacbporf ]
                    intronDQ._linked_to_introns   = [ intronAQ ]
                    intronAQ._linked_to_introns   = [ intronDQ ]

                    # append to results
                    results.append((intronDQ,None,tinypacbporf,intronAQ,intronS))

    # return the results
    return results

# end of function merge_pacbporfs_by_query_tinyexon_and_sbjct_intron


def merge_pacbporfs_by_sbjct_tinyexon_and_query_intron(pacbporfD,pacbporfA,
    orfSetObjS,verbose=False,**kwargs):
    """ """
    # input validation
    IsPacbPORF(pacbporfD)
    IsPacbPORF(pacbporfA)

    orfSetObjQ = []

    # edit **kwargs dictionary for some forced attributes
    _update_kwargs(kwargs,KWARGS_MAPPED_INTRON)
    if not kwargs.has_key('aligned_site_max_triplet_distance'):
        kwargs['aligned_site_max_triplet_distance'] = kwargs['max_aa_offset']

    # settings for minimal alignment entropy score
    # TODO TODO -> THIS MUST BE FIXED TO A NICE THRESHOLD VALUE!!!
    min_donor_site_alignment_entropy = 0.1
    min_acceptor_site_alignment_entropy = 0.1

    resultlistQ = merge_orfs_with_intron(
            pacbporfD.orfQ,pacbporfA.orfQ,**kwargs)
    resultlistS = merge_orfs_with_tinyexon(
            pacbporfD.orfS,pacbporfA.orfS,
            preceding_donor_sites=pacbporfD.orfS._donor_sites,
            subsequent_acceptor_sites=pacbporfA.orfS._acceptor_sites,
            orflist=orfSetObjS.orfs,**kwargs)

    # translate resultlists to dict: key == exon, value = [ {intronsD},{intronsS} ]
    resultdictS,key2exonS = _tinyexon_list_2_dict(resultlistS)

    # get unique list of donors & acceptors
    donorS = olba( list(Set([inD.donor for inD,te,inA in resultlistS ])), order_by='pos', reversed=True)
    accepS = olba( list(Set([inA.acceptor for inD,te,inA in resultlistS ])), order_by='pos')

    donorQ = olba( list(Set([intron.donor for intron in resultlistQ ])), order_by='pos', reversed=True)
    accepQ = olba( list(Set([intron.acceptor for intron in resultlistQ ])), order_by='pos')

    ## filter for alignable donor & acceptor sites
    kwargs['allow_non_canonical']               = True # True
    kwargs['aligned_site_max_triplet_distance'] = 0     # 2
    algdonors = _filter_for_alignable_splice_sites(donorQ,donorS,pacbporfD,**kwargs)
    algacceps = _filter_for_alignable_splice_sites(accepQ,accepS,pacbporfA,**kwargs)

    # remove sites at the wrong interface of the pabcporfs
    # required to avoid CoordinateOutOfrange exceptions lateron
    algdonors = _filter_for_entropy(algdonors,pacbporfD,'donor',
                min_alignment_entropy=0.0)
    algacceps = _filter_for_entropy(algacceps,pacbporfA,'acceptor',
                min_alignment_entropy=0.0)


    # remove sites with to low alignment entropy
    algdonors = _filter_for_entropy(algdonors,pacbporfD,'donor',
                min_alignment_entropy=min_donor_site_alignment_entropy)
    algacceps = _filter_for_entropy(algacceps,pacbporfA,'acceptor',
                min_alignment_entropy=min_acceptor_site_alignment_entropy)


    results = []

    ## search for positive cases with aligned donors
    ##  ========.....==========
    ##  ========.....=..=========
    for algDq,algDs in algdonors:
        for (intronDS,tinyexon,intronAS) in resultlistS:
            if intronDS.donor.pos != algDs.pos: continue
            for intronQ in resultlistQ:
                if intronQ.donor.pos != algDq.pos:
                    # not current aligned site
                    continue
                if (intronQ.acceptor.pos / 3) < pacbporfA.query_start:
                    # mapping 5' outside of extended PacbP -> continue
                    continue
                if (intronQ.acceptor.pos + tinyexon.length) > ((pacbporfA.query_end*3)+pacbporfA.orfQ.frame):
                    # mapping 3' outside of extended PacbP -> continue
                    continue
                if intronAS.acceptor.pos > pacbporfA._get_original_alignment_pos_end().sbjct_dna_end:
                    # mapping 3' outside of PacbP -> continue
                    continue

                # calculate distance between the acceptor sites;
                # this distance should match tinyexon.length
                try:
                    distAnt = pacbporfA.get_distance_aligned_nucleotide_positions(
                                sbjct = intronAS.acceptor.pos,
                                query = intronQ.acceptor.pos
                                )
                except CoordinateOutOfRange:
                    # not succesfull in obtaining distance -> no proper mapping
                    continue

                if distAnt == tinyexon.length:
                    ## get query sequence for comparison
                    #startPos,_phase = pacbporfA.dnaposition_query(intronQ.acceptor.pos)
                    #stopPos,_phase = pacbporfA.dnaposition_query(intronQ.acceptor.pos + tinyexon.length)
                    #start = pacbporfA._positions[startPos].query_pos
                    #stop  = pacbporfA._positions[stopPos].query_pos

                    # get sbjct sequence for comparison
                    # OLD code; gave pacb.exceptions.CoordinateOutOfRange
                    # errors when (by chance) the sbjct Orf was not long enough
                    # extended in the PacbPORF. Code is now replaced by
                    # literal coordinate retrieval on the sbjct Orf object
                    ###startPos,_phase = pacbporfA.dnaposition_sbjct(intronS.acceptor.pos)
                    ###stopPos,_phase = pacbporfA.dnaposition_sbjct(intronS.acceptor.pos + tinyexon.length)
                    ###start = pacbporfA._positions[startPos].sbjct_pos
                    ###stop  = pacbporfA._positions[stopPos].sbjct_pos
#                    start = intronS.acceptor.pos / 3
#                    stop  = (intronS.acceptor.pos + tinyexon.length - pacbporfA.orfS.frame) / 3

                    ## get query sequence for comparison
                    start = intronQ.acceptor.pos / 3
                    stop  = (intronQ.acceptor.pos + tinyexon.length - pacbporfA.orfQ.frame) / 3





                    if intronQ.phase != 0:
                        # correct for first AA which is broken by the splicesite
                        start+=1
                    if stop <= start:
                        # tinyexon is so tiny that is does not have a single
                        # full aligned AA -> discard here
                        continue

                    # actually get the query sequence
                    # TODO: get queryseq directly from ORF !?
                    queryseq,match,sbjctseq,coords = pacbporfA.alignmentpart_by_query(start,stop)
                    queryseq = queryseq.replace("-","")
                    query_start = coords[0]

                    ############################################################
                    if verbose:
                        print "S1 PACBPseqs:", tinyexon.proteinsequence(), queryseq,
                        print (tinyexon.acceptor.phase, tinyexon.donor.phase, tinyexon.length),
                        print (intronQ.phase, start, stop),
                        print ("Q",pacbporfD.orfQ.id,pacbporfA.orfQ.id),
                        print ("S",pacbporfD.orfS.id,pacbporfA.orfS.id)
                    ############################################################

                    if len(tinyexon.proteinsequence()) != len(queryseq):
                        # IMPORTANT!! this should not be possible here (but it occurs)
                        # TODO: solve this. For now, continue by breaking out
                        continue

                    # initialize extended tinyexon PacbPORF
                    from pacb import PacbP
                    pacbp = PacbP(input=(
                            queryseq,
                            tinyexon.proteinsequence(),
                            query_start,
                            (tinyexon.acceptor.pos / 3),
                            ) )

                    # continue if bitscore is <0
                    if pacbp.bitscore < 0: continue
                    pacbp.strip_unmatched_ends()
                    # continue if no fraction could be aligned
                    if len(pacbp) == 0: continue
                    tinypacbporf = pacbp2pacbporf(pacbp,pacbporfA.orfQ,tinyexon.orf)
                    tinypacbporf.extend_pacbporf_after_stops()

                    ############################################################
                    if verbose:
                        print "#### ALIGNED DONORS",
                        print ( intronQ.acceptor.pos, intronQ.acceptor.pos + tinyexon.length )
                        print "####", pacbporfD
                        print "####", tinyexon, tinyexon.proteinsequence(),queryseq
                        print "####", pacbporfA
                        print "####", tinypacbporf
                        tinypacbporf.print_protein_and_dna()
                        print intronDS
                        print intronAS
                        print intronQ
                        print ""
                    ############################################################

                    ############################################################
                    # set some meta-data properties to the intron objects
                    ############################################################
                    _score_introns_obtained_by_mapping(
                            intronDS,intronQ,pacbporfD,
                            tinypacbporf,source='ABGPmappingTE')
                    # set distance -> 0nt
                    intronAS._distance = 0
                    # add Alignment Positional Periphery Score into objects
                    succes = set_apps_intron_sbjct(intronAS,tinypacbporf,pacbporfA)
                    # set GFF fsource attribute for recognition of intron sources
                    intronAS._gff['fsource'] = "ABGPprojectingTE"

                    # create _linked_to_xxx attributes
                    intronDS._linked_to_pacbporfs = [ tinypacbporf ]
                    intronAS._linked_to_pacbporfs = [ tinypacbporf ]
                    intronQ._linked_to_pacbporfs  = [ tinypacbporf ]
                    intronDS._linked_to_introns   = [ intronAS ]
                    intronAS._linked_to_introns   = [ intronDS ]

                    # append to results
                    results.append((intronQ,intronDS,tinypacbporf,None,intronAS))

    ## search for positive cases with aligned acceptors
    ##    =========.....==========
    ##  ========..=.....==========
    for algAq,algAs in algacceps:
        for (intronDS,tinyexon,intronAS) in resultlistS:
            if intronAS.acceptor.pos != algAs.pos: continue
            for intronQ in resultlistQ:
                if intronQ.acceptor.pos != algAq.pos:
                    # not current aligned site
                    continue
                if (intronQ.donor.pos / 3) > pacbporfD.query_end:
                    # mapping 3' outside of extended PacbP -> continue
                    continue
                if (intronQ.donor.pos - tinyexon.length) <= ((pacbporfD.query_start*3)+pacbporfD.orfQ.frame):
                    # mapping 5' outside of extended PacbP -> continue
                    continue
                if intronDS.donor.pos < pacbporfD._get_original_alignment_pos_start().sbjct_dna_start:
                    # mapping 5' outside of PacbP -> continue
                    continue

                # calculate distance between the donor sites;
                # this distance should match tinyexon.length
                try:
                    distDnt = pacbporfD.get_distance_aligned_nucleotide_positions(
                                sbjct = intronDS.donor.pos,
                                query = intronQ.donor.pos
                                )
                except CoordinateOutOfRange:
                    # not succesfull in obtaining distance -> no proper mapping
                    continue

                if distDnt == tinyexon.length:
                    # get sbjct sequence for comparison
                    # OLD code; gave pacb.exceptions.CoordinateOutOfRange
                    # errors when (by chance) the sbjct Orf was not long enough
                    # extended in the PacbPORF. Code is now replaced by
                    # literal coordinate retrieval on the sbjct Orf object
                    ###startPos,_phase = pacbporfD.dnaposition_sbjct(intronS.donor.pos - tinyexon.length, forced_return=True)
                    ###stopPos,_phase = pacbporfD.dnaposition_sbjct(intronS.donor.pos)
                    ###start = pacbporfD._positions[startPos].sbjct_pos
                    ###stop  = pacbporfD._positions[stopPos].sbjct_pos
                    ###print start, stop, intronS.phase, pacbporfD.orfS.frame, tinyexon.acceptor.phase
                    start = ( intronQ.donor.pos - tinyexon.length - pacbporfD.orfQ.frame) / 3
                    stop  = intronQ.donor.pos / 3

                    if tinyexon.acceptor.phase != 0:
                        # correct for first AA which is broken by the splicesite
                        start+=1
                    if stop <= start:
                        # tinyexon is so tiny that is does not have a single
                        # full aligned AA -> discard here
                        continue

                    # actually get the sbjct sequence
                    # TODO: get queryseq directly from ORF !?
                    queryseq,match,sbjctseq,coords = pacbporfD.alignmentpart_by_query(start,stop)
                    queryseq = queryseq.replace("-","")
                    query_start = coords[0]

                    ############################################################
                    if verbose:
                        print "S2 PACBPseqs:", tinyexon.proteinsequence(), queryseq,
                        print (tinyexon.acceptor.phase, tinyexon.donor.phase, tinyexon.length),
                        print (intronQ.phase, start, stop),
                        print ("Q",pacbporfD.orfQ.id,pacbporfA.orfQ.id),
                        print ("S",pacbporfD.orfS.id,pacbporfA.orfS.id)
                    ############################################################

                    if len(tinyexon.proteinsequence()) != len(queryseq):
                        # IMPORTANT!! this should not be possible here (but it occurs)
                        # TODO: solve this. For now, continue by breaking out
                        continue

                    if query_start < 0:
                        # a match at the far front of an Orf (overlapping with the
                        # stop codon. Ignore this hit here
                        continue

                    # initialize extended tinyexon PacbPORF
                    from pacb import PacbP
                    pacbp = PacbP(input=(
                            queryseq,
                            tinyexon.proteinsequence(),
                            query_start,
                            (tinyexon.acceptor.pos / 3),
                            ) )

                    # continue if bitscore is <0
                    if pacbp.bitscore < 0: continue
                    pacbp.strip_unmatched_ends()
                    # continue if no fraction could be aligned
                    if len(pacbp) == 0: continue
                    tinypacbporf = pacbp2pacbporf(pacbp,pacbporfD.orfQ,tinyexon.orf)
                    tinypacbporf.extend_pacbporf_after_stops()

                    ############################################################
                    if verbose:
                        print "#### ALIGNED ACCEPTORS"
                        print "####", pacbporfD
                        print "####", tinyexon, tinyexon.proteinsequence(),queryseq
                        print "####", pacbporfA
                        print "####", tinypacbporf
                        tinypacbporf.print_protein_and_dna()
                        print intronDS
                        print intronAS
                        print intronQ
                        print ""
                    ############################################################

                    ############################################################
                    # set some meta-data properties to the intron objects
                    ############################################################
                    _score_introns_obtained_by_mapping(
                            intronAS,intronQ,tinypacbporf,pacbporfA,
                            source='ABGPmappingTE')
                    # set distance -> 0nt
                    intronDS._distance = 0
                    # add Alignment Positional Periphery Score into objects
                    succes = set_apps_intron_sbjct(intronDS,pacbporfD,tinypacbporf)
                    # set GFF fsource attribute for recognition of intron sources
                    intronDS._gff['fsource'] = "ABGPprojectingTE"

                    # create _linked_to_xxx attributes
                    intronDS._linked_to_pacbporfs = [ tinypacbporf ]
                    intronAS._linked_to_pacbporfs = [ tinypacbporf ]
                    intronQ._linked_to_pacbporfs  = [ tinypacbporf ]
                    intronDS._linked_to_introns   = [ intronAS ]
                    intronAS._linked_to_introns   = [ intronDS ]

                    # append to results
                    results.append((None,intronDS,tinypacbporf,intronQ,intronAS))

    # return the results
    return results

# end of function merge_pacbporfs_by_sbjct_tinyexon_and_query_intron




def merge_pacbporfs_by_sbjct_equal_length_exon_and_query_intron(pacbporfD,pacbporfA,
    orfSetObjS,verbose=False,**kwargs):
    """ """
    # input validation
    IsPacbPORF(pacbporfD)
    IsPacbPORF(pacbporfA)

    orfSetObjQ = []

    # edit **kwargs dictionary for some forced attributes
    _update_kwargs(kwargs,KWARGS_MAPPED_INTRON)
    if not kwargs.has_key('aligned_site_max_triplet_distance'):
        kwargs['aligned_site_max_triplet_distance'] = kwargs['max_aa_offset']


    kwargs['max_tinyexon_nt_length'] = 90
    kwargs['max_aa_offset'] = 0
    kwargs['min_donor_pssm_score'] = 1.0
    kwargs['min_acceptor_pssm_score'] = 1.0

    # settings for minimal alignment entropy score
    # TODO TODO -> THIS MUST BE FIXED TO A NICE THRESHOLD VALUE!!!
    min_donor_site_alignment_entropy = 0.1
    min_acceptor_site_alignment_entropy = 0.1

    resultlistQ = merge_orfs_with_intron(
            pacbporfD.orfQ,pacbporfA.orfQ,**kwargs)
    resultlistS = merge_orfs_with_tinyexon(
            pacbporfD.orfS,pacbporfA.orfS,
            preceding_donor_sites=pacbporfD.orfS._donor_sites,
            subsequent_acceptor_sites=pacbporfA.orfS._acceptor_sites,
            orflist=orfSetObjS.orfs,**kwargs)

    results = []

    for (intronDS,tinyexon,intronAS) in resultlistS:
        startSPos,_phase = pacbporfD.dnaposition_sbjct(intronDS.donor.pos,forced_return=True)
        stopSPos, _phase = pacbporfA.dnaposition_sbjct(intronAS.acceptor.pos,forced_return=True)
        # avoid IndexError (CoordOutOfRange)
        if startSPos >= len(pacbporfD._positions): continue
        if stopSPos  >= len(pacbporfA._positions): continue
        if startSPos < 0: continue
        if stopSPos  < 0: continue
        s2queryNTstart = pacbporfD._positions[startSPos].query_dna_start + intronDS.donor.phase
        s2queryNTstop  = pacbporfA._positions[stopSPos].query_dna_start + intronAS.acceptor.phase
        for intronQ in resultlistQ:
            if intronQ.donor.pos <= s2queryNTstart or intronQ.acceptor.pos >= s2queryNTstop:
                # not a scaffold of introns/exons
                continue

            # now check the equal length criterion
            cig5Pexon_length = intronQ.donor.pos - s2queryNTstart
            cig3Pexon_length = s2queryNTstop - intronQ.acceptor.pos
            length_dif = (cig3Pexon_length+cig5Pexon_length) - tinyexon.length
            if abs(length_dif) > kwargs['max_aa_offset']*3:
                continue
            if cig5Pexon_length <= 4:
                # this can be anything => non-alignable sites etc
                continue
            if cig3Pexon_length <= 4:
                # this can be anything => non-alignable sites etc
                continue

            # get sequences and check if PacbP(ORF)s can be created
            cig5Pquery_start = (s2queryNTstart - pacbporfD.orfQ.frame) / 3
            cig5Pquery_stop  = (intronQ.donor.pos- pacbporfD.orfQ.frame) / 3

            cig5Pquery_start = pacbporfD.orfQ.dnapos2aapos(s2queryNTstart)
            cig5Pquery_stop  = pacbporfD.orfQ.dnapos2aapos(intronQ.donor.pos)

            cig5Pquery_seq = pacbporfD.orfQ.getaas(
                    abs_pos_start = cig5Pquery_start,
                    abs_pos_end =  cig5Pquery_stop )
            if intronQ.phase != 0:
                cig5Pquery_seq = cig5Pquery_seq[1:]
                cig5Pquery_start +=1

            cig3Pquery_start = (intronQ.acceptor.pos- pacbporfA.orfQ.frame)/3
            cig3Pquery_stop  = (s2queryNTstop- pacbporfA.orfQ.frame)/3

            cig3Pquery_start = pacbporfD.orfQ.dnapos2aapos(intronQ.acceptor.pos)
            cig3Pquery_stop  = pacbporfD.orfQ.dnapos2aapos(s2queryNTstop)

            cig3Pquery_seq = pacbporfA.orfQ.getaas(
                    abs_pos_start = cig3Pquery_start ,
                    abs_pos_end = cig3Pquery_stop )

            # skip if no sequence left
            if not cig5Pquery_seq: continue
            if not cig3Pquery_seq: continue

            ####################################################################
            if verbose:
                print cig5Pexon_length, cig3Pexon_length, tinyexon.length,
                print cig3Pexon_length+cig5Pexon_length == tinyexon.length
                print (intronQ.donor.pos, intronQ.acceptor.pos), intronQ
                print (intronDS.donor.pos, intronDS.acceptor.pos), intronDS
                print tinyexon, tinyexon.proteinsequence(),
                print cig5Pquery_seq, cig3Pquery_seq,
                print len(tinyexon.proteinsequence()),
                print len(cig5Pquery_seq),
                print len(cig3Pquery_seq)
                print (intronAS.donor.pos, intronAS.acceptor.pos), intronAS
            ####################################################################

            # initialize extended tinyexon PacbPORF
            from pacb import PacbP
            cig5Ppacbp = PacbP(input=(
                    cig5Pquery_seq,
                    tinyexon.proteinsequence()[0:len(cig5Pquery_seq)],
                    cig5Pquery_start,
                    (tinyexon.acceptor.pos / 3),
                    ) )

            # continue if bitscore is <0
            if cig5Ppacbp.bitscore < 0: continue
            cig5Ppacbp.strip_unmatched_ends()
            # continue if no fraction could be aligned
            if len(cig5Ppacbp) <= 1: continue
            cig5Ptinypacbporf = pacbp2pacbporf(cig5Ppacbp,pacbporfD.orfQ,tinyexon.orf)
            cig5Ptinypacbporf.extend_pacbporf_after_stops()


            # initialize extended tinyexon PacbPORF
            cig3Ppacbp = PacbP(input=(
                    cig3Pquery_seq,
                    tinyexon.proteinsequence()[-len(cig3Pquery_seq):],
                    cig3Pquery_start,
                    ((tinyexon.acceptor.pos+cig5Pexon_length) / 3),
                    ) )

            # continue if bitscore is <0
            if cig3Ppacbp.bitscore < 0: continue
            cig3Ppacbp.strip_unmatched_ends()
            # continue if no fraction could be aligned
            if len(cig3Ppacbp) <= 1: continue
            cig3Ptinypacbporf = pacbp2pacbporf(cig3Ppacbp,pacbporfA.orfQ,tinyexon.orf)
            cig3Ptinypacbporf.extend_pacbporf_after_stops()

            ####################################################################
            if verbose:
                print cig5Ptinypacbporf
                cig5Ptinypacbporf.print_protein_and_dna()
                print cig3Ptinypacbporf
                cig3Ptinypacbporf.print_protein_and_dna()
                print ""
            ####################################################################

            # set distance -> 0nt
            intronQ._distance  = 0
            intronDS._distance = 0
            intronAS._distance = 0

            # add Alignment Positional Periphery Score into objects
            succes = set_apps_intron_sbjct(intronDS,pacbporfD,cig5Ptinypacbporf)
            succes = set_apps_intron_sbjct(intronAS,cig3Ptinypacbporf,pacbporfA)
            succes = set_apps_intron_query(intronQ,cig5Ptinypacbporf,cig3Ptinypacbporf)

            # set GFF fsource attribute for recognition of intron sources
            intronDS._gff['fsource'] = "ABGPprojectingTE"
            intronDS._gff['fsource'] = "ABGPprojectingTE"
            intronQ._gff['fsource']  = "ABGPprojectingTE"

            # create _linked_to_xxx attributes
            intronDS._linked_to_pacbporfs = [ cig5Ptinypacbporf,cig3Ptinypacbporf ]
            intronAS._linked_to_pacbporfs = [ cig5Ptinypacbporf,cig3Ptinypacbporf ]
            intronQ._linked_to_pacbporfs  = [ cig5Ptinypacbporf,cig3Ptinypacbporf ]
            intronDS._linked_to_introns   = [ intronAS ]
            intronAS._linked_to_introns   = [ intronDS ]

            # store to results!
            results.append((intronDS,cig5Ptinypacbporf,intronQ,cig3Ptinypacbporf,intronAS))

    # no results? done here!
    if not results: return results

    # finally, return only the highest scoring here
    bitscores = [ b.bitscore+d.bitscore for (a,b,c,d,e) in results ]
    max_score = max(bitscores)
    return [ results[bitscores.index(max_score)] ]

# end of function merge_pacbporfs_by_sbjct_equal_length_exon_and_query_intron



def merge_pacbporfs_by_query_equal_length_exon_and_sbjct_intron(pacbporfD,pacbporfA,
    orfSetObjQ,verbose=False,**kwargs):
    """ """
    # input validation
    IsPacbPORF(pacbporfD)
    IsPacbPORF(pacbporfA)

    orfSetObjS = []

    # edit **kwargs dictionary for some forced attributes
    _update_kwargs(kwargs,KWARGS_MAPPED_INTRON)
    if not kwargs.has_key('aligned_site_max_triplet_distance'):
        kwargs['aligned_site_max_triplet_distance'] = kwargs['max_aa_offset']


    kwargs['max_tinyexon_nt_length'] = 90
    kwargs['max_aa_offset'] = 0
    kwargs['min_donor_pssm_score'] = 1.0
    kwargs['min_acceptor_pssm_score'] = 1.0

    # settings for minimal alignment entropy score
    # TODO TODO -> THIS MUST BE FIXED TO A NICE THRESHOLD VALUE!!!
    min_donor_site_alignment_entropy = 0.1
    min_acceptor_site_alignment_entropy = 0.1

    resultlistS = merge_orfs_with_intron(
            pacbporfD.orfS,pacbporfA.orfS,**kwargs)
    resultlistQ = merge_orfs_with_tinyexon(
            pacbporfD.orfQ,pacbporfA.orfQ,
            preceding_donor_sites=pacbporfD.orfQ._donor_sites,
            subsequent_acceptor_sites=pacbporfA.orfQ._acceptor_sites,
            orflist=orfSetObjQ.orfs,**kwargs)

    results = []

    for (intronDQ,tinyexon,intronAQ) in resultlistQ:
        startQPos,_phase = pacbporfD.dnaposition_query(intronDQ.donor.pos,forced_return=True)
        stopQPos, _phase = pacbporfA.dnaposition_query(intronAQ.acceptor.pos,forced_return=True)
        # avoid IndexError (CoordOutOfRange)
        if startQPos >= len(pacbporfD._positions): continue
        if stopQPos  >= len(pacbporfA._positions): continue
        if startQPos < 0: continue
        if stopQPos  < 0: continue
        q2sbjctNTstart = pacbporfD._positions[startQPos].sbjct_dna_start + intronDQ.donor.phase
        q2sbjctNTstop  = pacbporfA._positions[stopQPos].sbjct_dna_start + intronAQ.acceptor.phase
        for intronS in resultlistS:
            if intronS.donor.pos <= q2sbjctNTstart or intronS.acceptor.pos >= q2sbjctNTstop:
                # not a scaffold of introns/exons
                continue

            # now check the equal length criterion
            cig5Pexon_length = intronS.donor.pos - q2sbjctNTstart
            cig3Pexon_length = q2sbjctNTstop - intronS.acceptor.pos
            length_dif = (cig3Pexon_length+cig5Pexon_length) - tinyexon.length
            if abs(length_dif) > kwargs['max_aa_offset']*3:
                continue
            if cig5Pexon_length <= 4:
                # this can be anything => non-alignable sites etc
                continue
            if cig3Pexon_length <= 4:
                # this can be anything => non-alignable sites etc
                continue

            # get sequences and check if PacbP(ORF)s can be created
            cig5Psbjct_start = pacbporfD.orfS.dnapos2aapos(q2sbjctNTstart)
            cig5Psbjct_stop  = pacbporfD.orfS.dnapos2aapos(intronS.donor.pos)

            cig5Psbjct_seq = pacbporfD.orfS.getaas(
                    abs_pos_start = cig5Psbjct_start,
                    abs_pos_end =  cig5Psbjct_stop )
            if intronS.phase != 0:
                cig5Psbjct_seq = cig5Psbjct_seq[1:]
                cig5Psbjct_start +=1

            cig3Psbjct_start = pacbporfD.orfS.dnapos2aapos(intronS.acceptor.pos)
            cig3Psbjct_stop  = pacbporfD.orfS.dnapos2aapos(q2sbjctNTstop)

            cig3Psbjct_seq = pacbporfA.orfS.getaas(
                    abs_pos_start = cig3Psbjct_start ,
                    abs_pos_end = cig3Psbjct_stop )

            # skip if no sequence left
            if not cig5Psbjct_seq: continue
            if not cig3Psbjct_seq: continue

            ####################################################################
            if verbose:
                print cig5Pexon_length, cig3Pexon_length, tinyexon.length,
                print cig3Pexon_length+cig5Pexon_length == tinyexon.length
                print (intronS.donor.pos, intronS.acceptor.pos), intronS
                print (intronDQ.donor.pos, intronDQ.acceptor.pos), intronDQ
                print tinyexon, tinyexon.proteinsequence(),
                print cig5Psbjct_seq, cig3Psbjct_seq,
                print len(tinyexon.proteinsequence()),
                print len(cig5Psbjct_seq),
                print len(cig3Psbjct_seq)
                print (intronAQ.donor.pos, intronAQ.acceptor.pos), intronAQ
            ####################################################################

            # initialize extended tinyexon PacbPORF
            from pacb import PacbP
            cig5Ppacbp = PacbP(input=(
                    tinyexon.proteinsequence()[0:len(cig5Psbjct_seq)],
                    cig5Psbjct_seq,
                    (tinyexon.acceptor.pos / 3),
                    cig5Psbjct_start,
                    ) )

            # continue if bitscore is <0
            if cig5Ppacbp.bitscore < 0: continue
            cig5Ppacbp.strip_unmatched_ends()
            # continue if no fraction could be aligned
            if len(cig5Ppacbp) <= 1: continue
            cig5Ptinypacbporf = pacbp2pacbporf(cig5Ppacbp,tinyexon.orf,pacbporfD.orfS)
            cig5Ptinypacbporf.extend_pacbporf_after_stops()

            # initialize extended tinyexon PacbPORF
            cig3Ppacbp = PacbP(input=(
                    tinyexon.proteinsequence()[-len(cig3Psbjct_seq):],
                    cig3Psbjct_seq,
                    ((tinyexon.acceptor.pos+cig5Pexon_length) / 3),
                    cig3Psbjct_start,
                    ) )

            # continue if bitscore is <0
            if cig3Ppacbp.bitscore < 0: continue
            cig3Ppacbp.strip_unmatched_ends()
            # continue if no fraction could be aligned
            if len(cig3Ppacbp) <= 1: continue
            cig3Ptinypacbporf = pacbp2pacbporf(cig3Ppacbp,tinyexon.orf,pacbporfA.orfS)
            cig3Ptinypacbporf.extend_pacbporf_after_stops()

            ####################################################################
            if verbose:
                print cig5Ptinypacbporf
                cig5Ptinypacbporf.print_protein_and_dna()
                print cig3Ptinypacbporf
                cig3Ptinypacbporf.print_protein_and_dna()
                print ""
            ####################################################################

            # set distance -> 0nt
            intronS._distance  = 0
            intronDQ._distance = 0
            intronAQ._distance = 0

            # add Alignment Positional Periphery Score into objects
            succes = set_apps_intron_query(intronDQ,pacbporfD,cig5Ptinypacbporf)
            succes = set_apps_intron_query(intronAQ,cig3Ptinypacbporf,pacbporfA)
            succes = set_apps_intron_sbjct(intronS,cig5Ptinypacbporf,cig3Ptinypacbporf)

            # set GFF fsource attribute for recognition of intron sources
            intronDQ._gff['fsource'] = "ABGPprojectingTE"
            intronAQ._gff['fsource'] = "ABGPprojectingTE"
            intronS._gff['fsource']  = "ABGPprojectingTE"

            # create _linked_to_xxx attributes
            intronDQ._linked_to_pacbporfs = [ cig5Ptinypacbporf,cig3Ptinypacbporf ]
            intronAQ._linked_to_pacbporfs = [ cig5Ptinypacbporf,cig3Ptinypacbporf ]
            intronS._linked_to_pacbporfs  = [ cig5Ptinypacbporf,cig3Ptinypacbporf ]
            intronDQ._linked_to_introns   = [ intronAQ ]
            intronAQ._linked_to_introns   = [ intronDQ ]

            # store to results!
            results.append((intronDQ,cig5Ptinypacbporf,intronS,cig3Ptinypacbporf,intronAQ))

    # no results? done here!
    if not results: return results

    # finally, return only the highest scoring here
    bitscores = [ b.bitscore+d.bitscore for (a,b,c,d,e) in results ]
    max_score = max(bitscores)
    return [ results[bitscores.index(max_score)] ]

# end of function merge_pacbporfs_by_query_equal_length_exon_and_sbjct_intron



def merge_pacbporfs_with_conserved_acceptor_introns(pacbporfD,pacbporfA,verbose=False,**kwargs):
    """
    Merge 2 PacbPORF objects by introns by enforcing conserved acceptors

    @attention: see orfs.merge_orfs_with_intron for **kwargs
    @attention: see functions._filter_for_alignable_splice_sites for **kwargs
    @attention: see functions._filter_for_entropy for **kwargs

    @type  pacbporfD: PacbPORF object
    @param pacbporfD: PacbPORF object that has to deliver PSSM donor objects

    @type  pacbporfA: PacbPORF object
    @param pacbporfA: PacbPORF object that has to deliver PSSM acceptor objects

    @type  verbose: Boolean
    @param verbose: print status/debugging messages to STDOUT

    @rtype:  list
    @return: list with ( intron, intron ), in query and sbjct
    """
    # input validation
    IsPacbPORF(pacbporfD)
    IsPacbPORF(pacbporfA)

    # edit **kwargs dictionary for some forced attributes
    _update_kwargs(kwargs,KWARGS_MAPPED_INTRON_STRACC)
    if not kwargs.has_key('aligned_site_max_triplet_distance'):
        kwargs['aligned_site_max_triplet_distance'] = kwargs['max_aa_offset']

    # settings for minimal alignment entropy score
    min_donor_site_alignment_entropy = 0.05     # KEEP little bit >0.0
    min_acceptor_site_alignment_entropy = 0.05  # KEEP little bit >0.0


    # get list of introns
    intronsQ = merge_orfs_with_intron(pacbporfD.orfQ,pacbporfA.orfQ,**kwargs)
    intronsS = merge_orfs_with_intron(pacbporfD.orfS,pacbporfA.orfS,**kwargs)

    # get unique list of donors & acceptors
    donorQ = olba( list(Set([inQ.donor for inQ in intronsQ ])), order_by='pos')
    donorS = olba( list(Set([inS.donor for inS in intronsS ])), order_by='pos')
    accepQ = olba( list(Set([inQ.acceptor for inQ in intronsQ ])), order_by='pos')
    accepS = olba( list(Set([inS.acceptor for inS in intronsS ])), order_by='pos')

    ############################################################################
    if verbose:
        print "dQ1", [ d.pos for d in donorQ ], "aQ1", [ a.pos for a in accepQ ],
        print kwargs['aligned_site_max_triplet_distance']
        print "dS1", [ d.pos for d in donorS ], "aS1", [ a.pos for a in accepS ],
        print kwargs['aligned_site_max_triplet_distance']
    ############################################################################

    # filter for alignable donor & acceptor sites
    kwargs['allow_non_canonical'] = kwargs['allow_non_canonical_donor']
    algdonors = _filter_for_alignable_splice_sites(donorQ,donorS,pacbporfD,**kwargs)
    kwargs['allow_non_canonical'] = kwargs['allow_non_canonical_acceptor']
    # set to 100% positionaly conserved acceptor
    bckp_triplet_dist = kwargs['aligned_site_max_triplet_distance']
    kwargs['aligned_site_max_triplet_distance'] = 0
    algacceps = _filter_for_alignable_splice_sites(accepQ,accepS,pacbporfA,**kwargs)
    kwargs['aligned_site_max_triplet_distance'] = bckp_triplet_dist

    ############################################################################
    if verbose:
        print "dQ2", [ _dq.pos for (_dq,_ds) in algdonors ],
        print "aQ2", [ _aq.pos for (_aq,_as) in algacceps ]
        print "dS2", [ _ds.pos for (_dq,_ds) in algdonors ],
        print "aS2", [ _as.pos for (_aq,_as) in algacceps ]
    ############################################################################

    # remove sites with to low alignment entropy
    algdonors = _filter_for_entropy(algdonors,pacbporfD,'donor',
                min_alignment_entropy=min_donor_site_alignment_entropy)
    algacceps = _filter_for_entropy(algacceps,pacbporfA,'acceptor',
                min_alignment_entropy=min_acceptor_site_alignment_entropy)

    ############################################################################
    if verbose:
        print "dQ3", [ _dq.pos for (_dq,_ds) in algdonors ],
        print "aQ3", [ _aq.pos for (_aq,_as) in algacceps ]
        print "dS3", [ _ds.pos for (_dq,_ds) in algdonors ],
        print "aS3", [ _as.pos for (_aq,_as) in algacceps ]
    ############################################################################


    # make unique position lists for quick lookup in intron lists
    dQpl = Set([ dQ.pos for dQ,dS in algdonors ])
    dSpl = Set([ dS.pos for dQ,dS in algdonors ])
    aQpl = Set([ aQ.pos for aQ,aS in algacceps ])
    aSpl = Set([ aS.pos for aQ,aS in algacceps ])

    # check exterior boundaries of PacbPORFs
    sposD = pacbporfD._get_original_alignment_pos_start()
    eposD = pacbporfD._get_original_alignment_pos_end()
    sposA = pacbporfA._get_original_alignment_pos_start()
    eposA = pacbporfA._get_original_alignment_pos_end()

    # now make list of alignable introns
    algintrons = []
    for intQ in intronsQ:
        # check if intron falls within the PacbPORF aligned area
        if intQ.donor.pos <= sposD.query_dna_start: continue
        if intQ.acceptor.pos >= eposA.query_dna_end: continue
        if intQ.donor.pos in dQpl and intQ.acceptor.pos in aQpl:
            # Query intron occurs in list of alignable splice sites!
            for intS in intronsS:
                # check if intron falls within the PacbPORF aligned area
                if intS.donor.pos <= sposD.sbjct_dna_start: continue
                if intS.acceptor.pos >= eposA.sbjct_dna_end: continue
                if intS.donor.pos in dSpl and intS.acceptor.pos in aSpl:
                    # Sbjct intron occurs as well in alignable splice sites!
                    if (intQ.donor,intS.donor) in algdonors and\
                    (intQ.acceptor,intS.acceptor) in algacceps:
                        # Sbjct & Query Donor & Acceptor are alignable!
                        algintrons.append( ( intQ, intS ) )

    ############################################################################
    # set some meta-data properties to the intron objects
    ############################################################################
    for intQ,intS in algintrons:
        distDnt = pacbporfD.get_distance_aligned_nucleotide_positions(
                        query = intQ.donor.pos, sbjct = intS.donor.pos
                        )
        distAnt = pacbporfA.get_distance_aligned_nucleotide_positions(
                        query = intQ.acceptor.pos, sbjct = intS.acceptor.pos
                        )

        # final distance check. kwargs['aligned_site_max_triplet_distance']
        # is applied on donor and acceptor site. This distance measured on the
        # protein sequence can be DOUBLED in case distDnt / distAnt are
        # opposite (+ and -). Check here if the protein sequence gap is
        # as well <= kwargs['aligned_site_max_triplet_distance'].
        if abs(distAnt - distDnt) > kwargs['aligned_site_max_triplet_distance']*3:
            continue

        # add distance score to introns
        intQ._distance = abs(distDnt) + abs(distAnt)
        intS._distance = abs(distDnt) + abs(distAnt)

        # add Alignment Positional Periphery Score into objects
        succes = set_apps_intron_query(intQ,pacbporfD,pacbporfA)
        succes = set_apps_intron_sbjct(intS,pacbporfD,pacbporfA)

        # set GFF fsource attribute for recognition of intron sources
        intQ._gff['fsource'] = "ABGPmapping"
        intS._gff['fsource'] = "ABGPmapping"

        ########################################################################
        if verbose:
            # some printing....
            print "Aligned introns:", ( intQ.donor.pos, intQ.acceptor.pos ) ,
            print ( intS.donor.pos, intS.acceptor.pos ),
            print "DIST:", distDnt, distAnt,
            print "[%s]" % kwargs['aligned_site_max_triplet_distance'],
            print "ENTROPY: %1.2f %1.2f" % (intQ._apps_donor, intQ._apps_accep),
            print "PSSM: (%1.2f %1.2f) (%1.2f %1.2f)" % (
                intQ.donor.pssm_score, intS.donor.pssm_score,
                intQ.acceptor.pssm_score, intS.acceptor.pssm_score,
                )
        ########################################################################

    # return lists of aligned introns
    return algintrons

# end of function merge_pacbporfs_with_conserved_acceptor_introns



def merge_pacbporfs_with_conserved_donor_introns(pacbporfD,pacbporfA,**kwargs):
    """
    @attention: see _merge_pacbporfs_with_single_conserved_splicesite function
                for argument documentation
    """
    return _merge_pacbporfs_with_single_conserved_splicesite(
                pacbporfD,pacbporfA,"donor",**kwargs)

# end of function merge_pacbporfs_with_conserved_donor_introns


def _merge_pacbporfs_with_single_conserved_splicesite(pacbporfD,pacbporfA,
    donororacceptor,verbose=False,**kwargs):
    """
    Merge 2 PacbPORF objects by introns by enforcing conserved donors

    @attention: see orfs.merge_orfs_with_intron for **kwargs
    @attention: see functions._filter_for_alignable_splice_sites for **kwargs
    @attention: see functions._filter_for_entropy for **kwargs

    @type  pacbporfD: PacbPORF object
    @param pacbporfD: PacbPORF object that has to deliver PSSM donor objects

    @type  pacbporfA: PacbPORF object
    @param pacbporfA: PacbPORF object that has to deliver PSSM acceptor objects

    @type  donororacceptor: string
    @param donororacceptor: literal 'donor' or 'acceptor'

    @type  verbose: Boolean
    @param verbose: print status/debugging messages to STDOUT

    @rtype:  list
    @return: list with ( intron, intron ), in query and sbjct
    """
    # input validation
    IsPacbPORF(pacbporfD)
    IsPacbPORF(pacbporfA)

    # edit **kwargs dictionary for some forced attributes
    _update_kwargs(kwargs,KWARGS_MAPPED_INTRON_STRDON)
    if not kwargs.has_key('aligned_site_max_triplet_distance'):
        kwargs['aligned_site_max_triplet_distance'] = kwargs['max_aa_offset']

    # settings for minimal alignment entropy score
    min_donor_site_alignment_entropy = 0.05     # KEEP little bit >0.0
    min_acceptor_site_alignment_entropy = 0.05  # KEEP little bit >0.0

    # get list of introns
    intronsQ = merge_orfs_with_intron(pacbporfD.orfQ,pacbporfA.orfQ,**kwargs)
    intronsS = merge_orfs_with_intron(pacbporfD.orfS,pacbporfA.orfS,**kwargs)

    # get unique list of donors & acceptors
    donorQ = olba( list(Set([inQ.donor for inQ in intronsQ ])), order_by='pos')
    donorS = olba( list(Set([inS.donor for inS in intronsS ])), order_by='pos')
    accepQ = olba( list(Set([inQ.acceptor for inQ in intronsQ ])), order_by='pos')
    accepS = olba( list(Set([inS.acceptor for inS in intronsS ])), order_by='pos')

    ############################################################################
    if verbose:
        print "dQ1", [ d.pos for d in donorQ ], "aQ1", [ a.pos for a in accepQ ],
        print kwargs['aligned_site_max_triplet_distance']
        print "dS1", [ d.pos for d in donorS ], "aS1", [ a.pos for a in accepS ],
        print kwargs['aligned_site_max_triplet_distance']
    ############################################################################


    if donororacceptor == 'acceptor':
        # filter for alignable donor sites with large offset
        kwargs['allow_non_canonical'] = kwargs['allow_non_canonical_donor']
        algdonors = _filter_for_alignable_splice_sites(donorQ,donorS,pacbporfD,**kwargs)

        # set to 100% positionaly conserved acceptor
        kwargs['allow_non_canonical'] = kwargs['allow_non_canonical_acceptor']
        bckp_triplet_dist = kwargs['aligned_site_max_triplet_distance']
        kwargs['aligned_site_max_triplet_distance'] = 0
        algacceps = _filter_for_alignable_splice_sites(accepQ,accepS,pacbporfA,**kwargs)
        kwargs['aligned_site_max_triplet_distance'] = bckp_triplet_dist

    elif donororacceptor == 'donor':
        # filter for alignable acceptor sites with large offset
        kwargs['allow_non_canonical'] = kwargs['allow_non_canonical_acceptor']
        algacceps = _filter_for_alignable_splice_sites(accepQ,accepS,pacbporfA,**kwargs)

        # set to 100% positionaly conserved donor
        kwargs['allow_non_canonical'] = kwargs['allow_non_canonical_donor']
        bckp_triplet_dist = kwargs['aligned_site_max_triplet_distance']
        kwargs['aligned_site_max_triplet_distance'] = 0
        algdonors = _filter_for_alignable_splice_sites(donorQ,donorS,pacbporfD,**kwargs)
        kwargs['aligned_site_max_triplet_distance'] = bckp_triplet_dist

    else:
        # problem. variable donororacceptor should be donor or acceptor
        raise "WHAHAHAHAHHAHAHAHAHA"


    ############################################################################
    if verbose:
        print "dQ2", [ _dq.pos for (_dq,_ds) in algdonors ],
        print "aQ2", [ _aq.pos for (_aq,_as) in algacceps ]
        print "dS2", [ _ds.pos for (_dq,_ds) in algdonors ],
        print "aS2", [ _as.pos for (_aq,_as) in algacceps ]
    ############################################################################

    # remove sites with to low alignment entropy
    algdonors = _filter_for_entropy(algdonors,pacbporfD,'donor',
                min_alignment_entropy=min_donor_site_alignment_entropy)
    algacceps = _filter_for_entropy(algacceps,pacbporfA,'acceptor',
                min_alignment_entropy=min_acceptor_site_alignment_entropy)

    ############################################################################
    if verbose:
        print "dQ3", [ _dq.pos for (_dq,_ds) in algdonors ],
        print "aQ3", [ _aq.pos for (_aq,_as) in algacceps ]
        print "dS3", [ _ds.pos for (_dq,_ds) in algdonors ],
        print "aS3", [ _as.pos for (_aq,_as) in algacceps ]
    ############################################################################


    # make unique position lists for quick lookup in intron lists
    dQpl = Set([ dQ.pos for dQ,dS in algdonors ])
    dSpl = Set([ dS.pos for dQ,dS in algdonors ])
    aQpl = Set([ aQ.pos for aQ,aS in algacceps ])
    aSpl = Set([ aS.pos for aQ,aS in algacceps ])

    # check exterior boundaries of PacbPORFs
    sposD = pacbporfD._get_original_alignment_pos_start()
    eposD = pacbporfD._get_original_alignment_pos_end()
    sposA = pacbporfA._get_original_alignment_pos_start()
    eposA = pacbporfA._get_original_alignment_pos_end()

    # now make list of alignable introns
    algintrons = []
    for intQ in intronsQ:
        # check if intron falls within the PacbPORF aligned area
        if intQ.donor.pos <= sposD.query_dna_start: continue
        if intQ.acceptor.pos >= eposA.query_dna_end: continue
        if intQ.donor.pos in dQpl and intQ.acceptor.pos in aQpl:
            # Query intron occurs in list of alignable splice sites!
            for intS in intronsS:
                # check if intron falls within the PacbPORF aligned area
                if intS.donor.pos <= sposD.sbjct_dna_start: continue
                if intS.acceptor.pos >= eposA.sbjct_dna_end: continue
                if intS.donor.pos in dSpl and intS.acceptor.pos in aSpl:
                    # Sbjct intron occurs as well in alignable splice sites!
                    if (intQ.donor,intS.donor) in algdonors and\
                    (intQ.acceptor,intS.acceptor) in algacceps:
                        # Sbjct & Query Donor & Acceptor are alignable!
                        algintrons.append( ( intQ, intS ) )

    ############################################################################
    # set some meta-data properties to the intron objects
    ############################################################################
    for intQ,intS in algintrons:
        distDnt = pacbporfD.get_distance_aligned_nucleotide_positions(
                        query = intQ.donor.pos, sbjct = intS.donor.pos
                        )
        distAnt = pacbporfA.get_distance_aligned_nucleotide_positions(
                        query = intQ.acceptor.pos, sbjct = intS.acceptor.pos
                        )

        # final distance check. kwargs['aligned_site_max_triplet_distance']
        # is applied on donor and acceptor site. This distance measured on the
        # protein sequence can be DOUBLED in case distDnt / distAnt are
        # opposite (+ and -). Check here if the protein sequence gap is
        # as well <= kwargs['aligned_site_max_triplet_distance'].
        if abs(distAnt - distDnt) > kwargs['aligned_site_max_triplet_distance']*3:
            continue

        # add distance score to introns
        intQ._distance = abs(distDnt) + abs(distAnt)
        intS._distance = abs(distDnt) + abs(distAnt)

        # add Alignment Positional Periphery Score into objects
        succes = set_apps_intron_query(intQ,pacbporfD,pacbporfA)
        succes = set_apps_intron_sbjct(intS,pacbporfD,pacbporfA)

        # set GFF fsource attribute for recognition of intron sources
        intQ._gff['fsource'] = "ABGPmapping"
        intS._gff['fsource'] = "ABGPmapping"

        ########################################################################
        if verbose:
            # some printing....
            print "Aligned introns:", ( intQ.donor.pos, intQ.acceptor.pos ) ,
            print ( intS.donor.pos, intS.acceptor.pos ),
            print "DIST:", distDnt, distAnt,
            print "[%s]" % kwargs['aligned_site_max_triplet_distance'],
            print "ENTROPY: %1.2f %1.2f" % (intQ._apps_donor, intQ._apps_accep),
            print "PSSM: (%1.2f %1.2f) (%1.2f %1.2f)" % (
                intQ.donor.pssm_score, intS.donor.pssm_score,
                intQ.acceptor.pssm_score, intS.acceptor.pssm_score,
                )
        ########################################################################

    # return lists of aligned introns
    return algintrons

# end of function _merge_pacbporfs_with_single_conserved_splicesite

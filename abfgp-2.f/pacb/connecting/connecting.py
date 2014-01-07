"""
PacbPORF connecting into a coding gene structure by the *best* connection
"""

# Module metadata
__authors__ = "Ate van der Burgt"
__license__ = "MIT"

# Imports
from pacb.validators import IsPacbPORF
from pacb.ordering import order_list_by_attribute as olba
from pacb.connecting.functions import (
    _update_kwargs,
    _filter_stopless_3n_introns,
    )

from pacb.connecting.sequenceerror import (
    merge_pacbporf_with_sequenceerror_in_query,
    merge_pacbporf_with_sequenceerror_in_sbjct,
    )

from pacb.connecting.mapping import (
    merge_pacbporfs_with_introns,
    merge_pacbporfs_by_tinyexons,
    merge_pacbporfs_with_closeby_independant_introns,
    merge_pacbporfs_with_phase_shift_introns,
    _filter_aligned_introns_on_pssm_entropy_combination,
    _filter_aligned_stopless_3n_introns,
    )
from pacb.connecting.projecting import (
    merge_pacbporfs_by_intron_tinyexon_intron_in_query,
    merge_pacbporfs_by_intron_tinyexon_intron_in_sbjct,
    merge_pacbporfs_by_intron_in_query,
    merge_pacbporfs_by_intron_in_sbjct,
    merge_pacbporfs_by_inframe_intron_in_query,
    merge_pacbporfs_by_inframe_intron_in_sbjct,
    merge_pacbporfs_by_two_tinyexons_in_query,
    merge_pacbporfs_by_two_tinyexons_in_sbjct,
    )
from pacb.connecting.mapping_complex import (
    merge_pacbporfs_by_query_tinyexon_and_sbjct_intron,
    merge_pacbporfs_by_sbjct_tinyexon_and_query_intron,
    merge_pacbporfs_by_sbjct_equal_length_exon_and_query_intron,
    merge_pacbporfs_by_query_equal_length_exon_and_sbjct_intron,
    merge_pacbporfs_with_conserved_acceptor_introns,
    merge_pacbporfs_with_conserved_donor_introns,
    )
from pacb.connecting.bridgeing import (
    merge_pacbporfs_with_query_intron_bridgeing
    )


# Python Imports
from sets import Set

# Global variable Imports of splice site settings
from settings.splicesites import KWARGS_SPLICESITES
from settings.pacbp import PACBPORF_HIGH_GAP_RATIO_THRESHOLD

def merge_pacbporfs(pacbporfD,pacbporfA,queryOrfSetObj,sbjctOrfSetObj,
    allow_query_projecting=True,
    allow_sbjct_projecting=True,
    allow_query_mapping=True,
    allow_sbjct_mapping=True,
    allow_projecting=True,
    allow_mapping=True,
    verbose=False):
    """
    Merge 2 PacbPORF objects with an interface into a gene structure

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

    # edit/create **kwargs dictionary for some forced attributes
    kwargs = {}
    _update_kwargs(kwargs,KWARGS_SPLICESITES)
   
    # deal with allow_xxx attributes
    if not allow_projecting:
        allow_query_projecting=False
        allow_sbjct_projecting=False
    if not allow_mapping:
        allow_query_mapping=False
        allow_sbjct_mapping=False
 
    # check if Orf objects of PacbPORFS are identical
    queryOrfsIdentical = pacbporfD.orfQ.id == pacbporfA.orfQ.id
    sbjctOrfsIdentical = pacbporfD.orfS.id == pacbporfA.orfS.id

    # return data structure of introns
    introns = { 'query': [], 'sbjct': [] }

    # Scan Orfs for splice sites.
    # This has probably been performed before, but when not done,
    # cached donor & acceptor sites lists seems to be empty -> no introns
    pacbporfD.orfQ.scan_orf_for_pssm_splice_sites(
            splicetype="donor",
            min_pssm_score=kwargs['min_donor_pssm_score'],
            allow_non_canonical=kwargs['allow_non_canonical_donor'],
            non_canonical_min_pssm_score=kwargs['non_canonical_min_donor_pssm_score'])
    pacbporfD.orfS.scan_orf_for_pssm_splice_sites(
            splicetype="donor",
            min_pssm_score=kwargs['min_donor_pssm_score'],
            allow_non_canonical=kwargs['allow_non_canonical_donor'],
            non_canonical_min_pssm_score=kwargs['non_canonical_min_donor_pssm_score'])
    pacbporfA.orfQ.scan_orf_for_pssm_splice_sites(
            splicetype="acceptor",
            min_pssm_score=kwargs['min_acceptor_pssm_score'],
            allow_non_canonical=kwargs['allow_non_canonical_acceptor'],
            non_canonical_min_pssm_score=kwargs['non_canonical_min_acceptor_pssm_score'])
    pacbporfA.orfS.scan_orf_for_pssm_splice_sites(
            splicetype="acceptor",
            min_pssm_score=kwargs['min_acceptor_pssm_score'],
            allow_non_canonical=kwargs['allow_non_canonical_acceptor'],
            non_canonical_min_pssm_score=kwargs['non_canonical_min_acceptor_pssm_score'])

    if not queryOrfsIdentical and not sbjctOrfsIdentical:

        introns1 = merge_pacbporfs_with_introns(
                        pacbporfD,pacbporfA)
        # filter for **best** candidates based on PSSM/entropy combination
        introns1 = _filter_aligned_introns_on_pssm_entropy_combination(introns1)


        if pacbporfD.gap_ratio_score() < PACBPORF_HIGH_GAP_RATIO_THRESHOLD and\
        pacbporfA.gap_ratio_score() < PACBPORF_HIGH_GAP_RATIO_THRESHOLD:
            introns2 = merge_pacbporfs_with_closeby_independant_introns(
                            pacbporfD,pacbporfA)
            introns3 = merge_pacbporfs_with_phase_shift_introns(
                            pacbporfD,pacbporfA)
            introns4 = merge_pacbporfs_by_tinyexons(
                            pacbporfD,pacbporfA,
                            queryOrfSetObj,sbjctOrfSetObj)
    
            introns5 = merge_pacbporfs_by_query_tinyexon_and_sbjct_intron(
                            pacbporfD,pacbporfA,queryOrfSetObj)
    
            introns6 = merge_pacbporfs_by_sbjct_tinyexon_and_query_intron(
                            pacbporfD,pacbporfA,sbjctOrfSetObj)
    
            introns7 = merge_pacbporfs_by_sbjct_equal_length_exon_and_query_intron(
                            pacbporfD,pacbporfA,sbjctOrfSetObj)
    
            introns8 = merge_pacbporfs_by_query_equal_length_exon_and_sbjct_intron(
                            pacbporfD,pacbporfA,queryOrfSetObj)
        else:
            # do not allow more complex intron merging
            introns2 = {}
            introns3 = {}
            introns4 = {}
            introns5 = {}
            introns6 = {}
            introns7 = {}
            introns8 = {}


        introns9 = merge_pacbporfs_with_conserved_acceptor_introns(
                        pacbporfD,pacbporfA)
        # filter for **best** candidates based on PSSM/entropy combination
        introns9 = _filter_aligned_introns_on_pssm_entropy_combination(introns9)

        introns10 = merge_pacbporfs_with_conserved_donor_introns(
                        pacbporfD,pacbporfA)
        # filter for **best** candidates based on PSSM/entropy combination
        introns10 = _filter_aligned_introns_on_pssm_entropy_combination(introns10)


        # store introns obtained by most simplest case projecting/mapping
        introns['query'].extend( Set([ intrQ for (intrQ,intrS) in introns1 ]) )
        introns['sbjct'].extend( Set([ intrS for (intrQ,intrS) in introns1 ]) )

        # only store introns from intron2 that are NOT encountered already in introns1
        keysQ = [ (intron.donor.pos,intron.acceptor.pos) for intron in introns['query'] ]
        keysS = [ (intron.donor.pos,intron.acceptor.pos) for intron in introns['sbjct'] ]
        for (intrQ,intrS,cigpacbp) in introns2:
            k1 = (intrQ.donor.pos,intrQ.acceptor.pos)
            k2 = (intrS.donor.pos,intrS.acceptor.pos)
            if k1 not in keysQ and k2 not in keysS:
                introns['query'].append(intrQ)
                introns['sbjct'].append(intrS)

        # only store introns from intron3 that are NOT encountered already in introns1
        keysQ = [ (intron.donor.pos,intron.acceptor.pos) for intron in introns['query'] ]
        keysS = [ (intron.donor.pos,intron.acceptor.pos) for intron in introns['sbjct'] ]
        for (intrQ,intrS) in introns3:
            k1 = (intrQ.donor.pos,intrQ.acceptor.pos)
            k2 = (intrS.donor.pos,intrS.acceptor.pos)
            if k1 not in keysQ and k2 not in keysS:
                introns['query'].append(intrQ)
                introns['sbjct'].append(intrS)

        # only store introns from intron4 that are NOT encountered already in introns1
        keysQ = [ (intron.donor.pos,intron.acceptor.pos) for intron in introns['query'] ]
        keysS = [ (intron.donor.pos,intron.acceptor.pos) for intron in introns['sbjct'] ]
        for (intrQ,intrS,pacbporf,intrQ2,intrS2) in introns4:
            k1 = (intrQ.donor.pos,intrQ.acceptor.pos)
            k2 = (intrS.donor.pos,intrS.acceptor.pos)
            k3 = (intrQ2.donor.pos,intrQ2.acceptor.pos)
            k4 = (intrS2.donor.pos,intrS2.acceptor.pos)
            if k1 not in keysQ and k2 not in keysS and k3 not in keysQ and k4 not in keysS:
                introns['query'].append(intrQ)
                introns['sbjct'].append(intrS)
                introns['query'].append(intrQ2)
                introns['sbjct'].append(intrS2)

        # only store introns from intron5 that are NOT encountered already in introns1
        keysQ = [ (intron.donor.pos,intron.acceptor.pos) for intron in introns['query'] ]
        keysS = [ (intron.donor.pos,intron.acceptor.pos) for intron in introns['sbjct'] ]
        for (intrQ,intrS,pacbporf,intrQ2,intrS2) in introns4:
            if intrQ:   k1 = (intrQ.donor.pos,intrQ.acceptor.pos)
            else:       k1 = None
            if intrS:   k2 = (intrS.donor.pos,intrS.acceptor.pos)
            else:       k2 = None
            if intrQ2:  k3 = (intrQ2.donor.pos,intrQ2.acceptor.pos)
            else:       k3 = None
            if intrS2:  k4 = (intrS2.donor.pos,intrS2.acceptor.pos)
            else:       k4 = None
            if k1 not in keysQ and k2 not in keysS and k3 not in keysQ and k4 not in keysS:
                introns['query'].append(intrQ)
                introns['sbjct'].append(intrS)
                introns['query'].append(intrQ2)
                introns['sbjct'].append(intrS2)

        # only store introns from intron6 that are NOT encountered already in introns1
        keysQ = [ (intron.donor.pos,intron.acceptor.pos) for intron in introns['query'] ]
        keysS = [ (intron.donor.pos,intron.acceptor.pos) for intron in introns['sbjct'] ]
        for (intrQ,intrS,pacbporf,intrQ2,intrS2) in introns6:
            if intrQ:   k1 = (intrQ.donor.pos,intrQ.acceptor.pos)
            else:       k1 = None
            if intrS:   k2 = (intrS.donor.pos,intrS.acceptor.pos)
            else:       k2 = None
            if intrQ2:  k3 = (intrQ2.donor.pos,intrQ2.acceptor.pos)
            else:       k3 = None
            if intrS2:  k4 = (intrS2.donor.pos,intrS2.acceptor.pos)
            else:       k4 = None
            if k1 not in keysQ and k2 not in keysS and k3 not in keysQ and k4 not in keysS:
                introns['query'].append(intrQ)
                introns['sbjct'].append(intrS)
                introns['query'].append(intrQ2)
                introns['sbjct'].append(intrS2)

        # remove the 'None' in introns['sbjct'] due to latest addition
        while None in introns['query']: introns['query'].remove(None)
        while None in introns['sbjct']: introns['sbjct'].remove(None)

        # only store introns from intron7 that are NOT encountered already in introns1
        keysQ = [ (intron.donor.pos,intron.acceptor.pos) for intron in introns['query'] ]
        keysS = [ (intron.donor.pos,intron.acceptor.pos) for intron in introns['sbjct'] ]
        for (intrS,pacbporf1,intrQ,pacbporf2,intrS2) in introns7:
            k1 = (intrQ.donor.pos,intrQ.acceptor.pos)
            k2 = (intrS.donor.pos,intrS.acceptor.pos)
            k3 = (intrS2.donor.pos,intrS2.acceptor.pos)
            if k1 not in keysQ and k2 not in keysS and k3 not in keysS:
                introns['query'].append(intrQ)
                introns['sbjct'].append(intrS)
                introns['sbjct'].append(intrS2)

        # only store introns from intron8 that are NOT encountered already in introns1
        keysQ = [ (intron.donor.pos,intron.acceptor.pos) for intron in introns['query'] ]
        keysS = [ (intron.donor.pos,intron.acceptor.pos) for intron in introns['sbjct'] ]
        for (intrQ,pacbporf1,intrS,pacbporf2,intrQ2) in introns8:
            k1 = (intrQ.donor.pos,intrQ.acceptor.pos)
            k2 = (intrS.donor.pos,intrS.acceptor.pos)
            k3 = (intrQ2.donor.pos,intrQ2.acceptor.pos)
            if k1 not in keysQ and k2 not in keysS and k3 not in keysQ:
                introns['query'].append(intrQ)
                introns['query'].append(intrQ2)
                introns['sbjct'].append(intrS)


        # only store introns from introns9 that are NOT encountered already in introns1
        keysQ = [ (intron.donor.pos,intron.acceptor.pos) for intron in introns['query'] ]
        keysS = [ (intron.donor.pos,intron.acceptor.pos) for intron in introns['sbjct'] ]
        for (intrQ,intrS) in introns9:
            k1 = (intrQ.donor.pos,intrQ.acceptor.pos)
            k2 = (intrS.donor.pos,intrS.acceptor.pos)
            if k1 == (2163, 2283):
                print "STRACC", k1, intrQ, k1 not in keysQ
                print "STRACC", k1, intrS, k2 not in keysS
            # do NOT check if any of the introns is present yet;
            # allow addition of each of these
            if k1 not in keysQ:
                introns['query'].append(intrQ)
            if k2 not in keysS:
                introns['sbjct'].append(intrS)

        # only store introns from introns10 that are NOT encountered already in introns1
        keysQ = [ (intron.donor.pos,intron.acceptor.pos) for intron in introns['query'] ]
        keysS = [ (intron.donor.pos,intron.acceptor.pos) for intron in introns['sbjct'] ]
        for (intrQ,intrS) in introns10:
            k1 = (intrQ.donor.pos,intrQ.acceptor.pos)
            k2 = (intrS.donor.pos,intrS.acceptor.pos)
            if k1 == (1642, 1858):
                print "STRDON", k1, intrQ, k1 not in keysQ
                print "STRDON", k1, intrS, k2 not in keysS
            # do NOT check if any of the introns is present yet;
            # allow addition of each of these
            if k1 not in keysQ:
                introns['query'].append(intrQ)
            if k2 not in keysS:
                introns['sbjct'].append(intrS)



        # finally, do the bridging thingy
        introns0 = merge_pacbporfs_with_query_intron_bridgeing(pacbporfD,pacbporfA)

        # only store introns from introns0 that are NOT encountered already in introns1
        keysQ = [ (intron.donor.pos,intron.acceptor.pos) for intron in introns['query'] ]
        for intrQ in introns0:
            if intrQ.coords() not in keysQ:
                introns['query'].append(intrQ)


        #introns['query'].extend([ intrQ for (intrQ,intrS) in introns1 ] )
        #introns['query'].extend([ intrQ for (intrQ,intrS,cigpacbp) in introns2 ] )
        #introns['query'].extend([ intrQ for (intrQ,intrS) in introns3 ] )
        #introns['query'].extend([ intrQ for (intrQ,a,b,c,d) in introns4 ] )
        #introns['query'].extend([ intrQ for (a,b,c,intrQ,d) in introns4 ] )
        #introns['query'].extend([ intrQ for (intrQ,a,b,c,d) in introns5 ] )
        #introns['query'].extend([ intrQ for (a,b,c,intrQ,d) in introns5 ] )
        #introns['sbjct'].extend([ intrS for (intrQ,intrS) in introns1 ] )
        #introns['sbjct'].extend([ intrS for (intrQ,intrS,cigpacbp) in introns2 ] )
        #introns['sbjct'].extend([ intrS for (intrQ,intrS) in introns3 ] )
        #introns['sbjct'].extend([ intrS for (a,intrS,b,c,d) in introns4 ] )
        #introns['sbjct'].extend([ intrS for (a,b,c,d,intrS) in introns4 ] )
        #introns['sbjct'].extend([ intrS for (a,intrS,b,c,d) in introns5 ] )
        #introns['sbjct'].extend([ intrS for (a,b,c,d,intrS) in introns5 ] )

        # remove the 'None' in introns['sbjct'] due to latest addition
        while None in introns['query']: introns['query'].remove(None)
        while None in introns['sbjct']: introns['sbjct'].remove(None)


    elif not queryOrfsIdentical:
        seqerror = merge_pacbporf_with_sequenceerror_in_query(
                        pacbporfD,pacbporfA)
        introns1 = merge_pacbporfs_by_intron_in_query(
                        pacbporfD,pacbporfA)


        if pacbporfD.gap_ratio_score() < PACBPORF_HIGH_GAP_RATIO_THRESHOLD and\
        pacbporfA.gap_ratio_score() < PACBPORF_HIGH_GAP_RATIO_THRESHOLD:
            introns2 = merge_pacbporfs_by_intron_tinyexon_intron_in_query(
                            pacbporfD,pacbporfA,queryOrfSetObj)
            introns3 = merge_pacbporfs_by_two_tinyexons_in_query(
                            pacbporfD,pacbporfA,queryOrfSetObj)
        else:
            # do not allow more complex intron merging
            introns2 = {}
            introns3 = {}

        # store sequencerror if it exists
        if seqerror: introns['query'].append( seqerror )

        # store introns obtained by most simplest case projecting/mapping
        introns['query'].extend([ prj.projected_introns[0] for prj in introns1 ] )

        # only store introns from intron2 that are NOT encountered already in introns1
        keys = [ (intron.donor.pos,intron.acceptor.pos) for intron in introns['query'] ]
        for (intr1,intr2,exon) in introns2:
            k1 = (intr1.donor.pos,intr1.acceptor.pos)
            k2 = (intr2.donor.pos,intr2.acceptor.pos)
            if k1 not in keys and k2 not in keys:
                introns['query'].append(intr1)
                introns['query'].append(intr2)

        # only store introns from intron2 that are NOT encountered already in introns1
        keys = [ (intron.donor.pos,intron.acceptor.pos) for intron in introns['query'] ]
        for (intr1,intr2,intr3,exon1,exon2) in introns3:
            k1 = (intr1.donor.pos,intr1.acceptor.pos)
            k2 = (intr2.donor.pos,intr2.acceptor.pos)
            k3 = (intr3.donor.pos,intr3.acceptor.pos)
            if k1 not in keys and k2 not in keys and k3 not in keys:
                introns['query'].append(intr1)
                introns['query'].append(intr2)
                introns['query'].append(intr3)


        if not introns['query'] and allow_sbjct_mapping and allow_query_mapping:
            # just bridge Orfs by **best** intron(s).
            introns0 = merge_pacbporfs_with_query_intron_bridgeing(pacbporfD,pacbporfA)

            # potential stopless 3n intron in SBJCT
            introns1 = merge_pacbporfs_with_introns(
                            pacbporfD,pacbporfA)
            # filter for **best** candidates based on PSSM/entropy combination
            introns1 = _filter_aligned_introns_on_pssm_entropy_combination(introns1)
            # apply stopless3n intron filtering
            introns1 = _filter_aligned_stopless_3n_introns(introns1)

            introns2 = merge_pacbporfs_with_closeby_independant_introns(
                            pacbporfD,pacbporfA)

            if pacbporfD.gap_ratio_score() < PACBPORF_HIGH_GAP_RATIO_THRESHOLD and\
            pacbporfA.gap_ratio_score() < PACBPORF_HIGH_GAP_RATIO_THRESHOLD:
                introns3 = merge_pacbporfs_with_phase_shift_introns(
                                pacbporfD,pacbporfA)
                # filter for **best** candidates based on PSSM/entropy combination
                introns3 = _filter_aligned_introns_on_pssm_entropy_combination(introns3)
                # apply stopless3n intron filtering
                introns3 = _filter_aligned_stopless_3n_introns(introns3)

            else:
                # do not allow more complex intron merging
                introns3 = {}


            # only store introns from that are NOT encountered already
            keys = [ intron.coords() for intron in introns['query'] ]
            for intrQ,intrS in introns1:
                if intrQ.coords() not in keys:
                    introns['query'].append(intrQ)
                    keys = [ intron.coords() for intron in introns['query'] ]
            for (intrQ,intrS,cigpacbp) in introns2:
                if intrQ.coords() not in keys:
                    introns['query'].append(intrQ)
                    keys = [ intron.coords() for intron in introns['query'] ]
            for intrQ,intrS in introns3:
                if intrQ.coords() not in keys:
                    introns['query'].append(intrQ)
                    keys = [ intron.coords() for intron in introns['query'] ]
            for intron in introns0:
                if intron.coords() not in keys:
                    introns['query'].append(intron)
                    keys = [ intron.coords() for intron in introns['query'] ]


            keys = [ intron.coords() for intron in introns['sbjct'] ]
            for intrQ,intrS in introns1:
                if intrS.coords() not in keys:
                    introns['query'].append(intrS)
                    keys = [ intron.coords() for intron in introns['sbjct'] ]
            for (intrQ,intrS,cigpacbp) in introns2:
                if intrS.coords() not in keys:
                    introns['query'].append(intrS)
                    keys = [ intron.coords() for intron in introns['sbjct'] ]
            for intrQ,intrS in introns3:
                if intrS.coords() not in keys:
                    introns['query'].append(intrS)
                    keys = [ intron.coords() for intron in introns['sbjct'] ]

        elif not introns['query']:

            # just bridge Orfs by **best** intron(s).
            introns0 = merge_pacbporfs_with_query_intron_bridgeing(pacbporfD,pacbporfA)
            # only store introns from that are NOT encountered already
            keys = [ intron.coords() for intron in introns['query'] ]
            for intron in introns0:
                if intron.coords() not in keys:
                    introns['query'].append(intron)
        else:
            # projecting introns yielded results; do not try mapping
            pass


    elif not sbjctOrfsIdentical:
        introns1 = merge_pacbporfs_by_intron_in_sbjct(
                        pacbporfD,pacbporfA)

        if pacbporfD.gap_ratio_score() < PACBPORF_HIGH_GAP_RATIO_THRESHOLD and\
        pacbporfA.gap_ratio_score() < PACBPORF_HIGH_GAP_RATIO_THRESHOLD:
            introns2 = merge_pacbporfs_by_intron_tinyexon_intron_in_sbjct(
                            pacbporfD,pacbporfA,sbjctOrfSetObj)
            introns3 = merge_pacbporfs_by_two_tinyexons_in_sbjct(
                            pacbporfD,pacbporfA,sbjctOrfSetObj)
        else:
            # do not allow more complex intron merging
            introns2 = {}
            introns3 = {}

        # store introns obtained by most simplest case projecting/mapping
        introns['sbjct'].extend([ prj.projected_introns[0] for prj in introns1 ] )

        # only store introns from intron2 that are NOT encountered already in introns1
        keys = [ (intron.donor.pos,intron.acceptor.pos) for intron in introns['sbjct'] ]
        for (intr1,intr2,exon) in introns2:
            k1 = (intr1.donor.pos,intr1.acceptor.pos)
            k2 = (intr2.donor.pos,intr2.acceptor.pos)
            if k1 not in keys and k2 not in keys:
                introns['sbjct'].append(intr1)
                introns['sbjct'].append(intr2)

        # only store introns from intron2 that are NOT encountered already in introns1
        keys = [ (intron.donor.pos,intron.acceptor.pos) for intron in introns['sbjct'] ]
        for (intr1,intr2,intr3,exon1,exon2) in introns3:
            k1 = (intr1.donor.pos,intr1.acceptor.pos)
            k2 = (intr2.donor.pos,intr2.acceptor.pos)
            k3 = (intr3.donor.pos,intr3.acceptor.pos)
            if k1 not in keys and k2 not in keys and k3 not in keys:
                introns['sbjct'].append(intr1)
                introns['sbjct'].append(intr2)
                introns['sbjct'].append(intr3)


        if not introns['sbjct'] and allow_sbjct_mapping and allow_query_mapping:
            # potential stopless 3n intron in QUERY
            introns1 = merge_pacbporfs_with_introns(
                            pacbporfD,pacbporfA)
            # filter for **best** candidates based on PSSM/entropy combination
            introns1 = _filter_aligned_introns_on_pssm_entropy_combination(introns1)
            # apply stopless3n intron filtering
            introns1 = _filter_aligned_stopless_3n_introns(introns1)

            introns2 = merge_pacbporfs_with_closeby_independant_introns(
                            pacbporfD,pacbporfA)


            if pacbporfD.gap_ratio_score() < PACBPORF_HIGH_GAP_RATIO_THRESHOLD and\
            pacbporfA.gap_ratio_score() < PACBPORF_HIGH_GAP_RATIO_THRESHOLD:
                introns3 = merge_pacbporfs_with_phase_shift_introns(
                                pacbporfD,pacbporfA)
                # filter for **best** candidates based on PSSM/entropy combination
                introns3 = _filter_aligned_introns_on_pssm_entropy_combination(introns3)
                # apply stopless3n intron filtering
                introns3 = _filter_aligned_stopless_3n_introns(introns3)
            else:
                # do not allow more complex intron merging
                introns3 = {}

            # store introns
            introns['query'].extend( Set([ intrQ for (intrQ,intrS) in introns1 ]) )
            introns['sbjct'].extend( Set([ intrS for (intrQ,intrS) in introns1 ]) )
            introns['query'].extend([ intrQ for (intrQ,intrS,cigpacbp) in introns2 ] )
            introns['query'].extend([ intrQ for (intrQ,intrS) in introns3 ] )
            introns['sbjct'].extend([ intrS for (intrQ,intrS,cigpacbp) in introns2 ] )
            introns['sbjct'].extend([ intrS for (intrQ,intrS) in introns3 ] )
        else:
            # projecting introns yielded results; do not try mapping
            pass

    elif queryOrfsIdentical and sbjctOrfsIdentical:
        if allow_query_mapping:
            introns1 = merge_pacbporfs_by_inframe_intron_in_query(
                            pacbporfD,pacbporfA)
        else:
            # no mapping (unigene or continious alignment provided)
            introns1 = []

        if allow_sbjct_mapping:
            introns2 = merge_pacbporfs_by_inframe_intron_in_sbjct(
                            pacbporfD,pacbporfA)
        else:
            # no mapping (unigene or continious alignment provided)
            introns2 = []

        if allow_sbjct_mapping and allow_query_mapping:
            introns3 = merge_pacbporfs_with_introns(
                            pacbporfD,pacbporfA)
            # filter for **best** candidates based on PSSM/entropy combination
            introns3 = _filter_aligned_introns_on_pssm_entropy_combination(introns3)
            # apply stopless3n intron filtering
            introns3 = _filter_aligned_stopless_3n_introns(introns3)

        else:
            # no mapping (unigene or continious alignment provided)
            introns3 = []

        #introns4 = merge_pacbporfs_with_closeby_independant_introns(
        #                pacbporfD,pacbporfA)
        #introns5 = merge_pacbporfs_with_phase_shift_introns(
        #                pacbporfD,pacbporfA)

        introns['query'].extend([ prj.projected_introns[0] for prj in introns1 ] )
        introns['sbjct'].extend([ prj.projected_introns[0] for prj in introns2 ] )
        introns['query'].extend([ intrQ for (intrQ,intrS) in introns3 ] )
        introns['sbjct'].extend([ intrS for (intrQ,intrS) in introns3 ] )

    else:
        # none of these cases; allow_projecting or allow_mapping == False!
        pass


    # Filter for stopless3n introns
    introns['query'] = _filter_stopless_3n_introns(introns['query'])
    introns['sbjct'] = _filter_stopless_3n_introns(introns['sbjct'])

    # return list of introns
    return introns

# end of function merge_pacbporfs



def _update_P_introns(allintrondict,intronlist):
    """ store introns obtained by mapping """
    pass

def _update_Q_introns(allintrondict,intronlist):
    """ store introns obtained by query projecting """
    pass

def _update_S_introns(allintrondict,intronlist):
    """ store introns obtained by sbjct projecting """
    pass

def _update_QS_introns(allintrondict,intronlist):
    """ store Closeby Independant Gained introns """
    pass

def _update_PP_introns(allintrondict,intronlist):
    """ store tinyexon introns obtained by mapping """
    pass

def _update_QQ_introns(allintrondict,intronlist):
    """ store query tinyexon introns obtained by projecting """
    pass

def _update_SS_introns(allintrondict,intronlist):
    """ store sbjct tinyexon introns obtained by projecting """
    pass


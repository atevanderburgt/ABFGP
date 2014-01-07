"""
PacbPORF connection by mapping introns
"""

# Module metadata
__authors__ = "Ate van der Burgt"
__license__ = "MIT"

# Imports
from pacb.validators import IsPacbPORF
from pacb.ordering import order_list_by_attribute as olba
from pacb.connecting.orfs import merge_orfs_with_intron, merge_orfs_with_tinyexon
from pacb.conversion import pacbp2pacbporf, pacbp_from_clustalw
from pacb.connecting.functions import (
    _update_kwargs,
    _filter_for_entropy,
    _filter_for_alignable_splice_sites,
    _score_introns_obtained_by_mapping,
    set_apps_intron_query,
    set_apps_intron_sbjct,
    )

# Other Imports
from lib_clustalw import clustalw

# Python Imports
from sets import Set

# Global variable Imports
from settings.splicesites import (
    KWARGS_MAPPED_INTRON,
    KWARGS_CLOSEBY_INDEPENDANT_INTRON_GAIN,
    KWARGS_PHASE_SHIFT_INTRON,
)


def merge_pacbporfs_with_introns(pacbporfD,pacbporfA,verbose=False,**kwargs):
    """
    Merge 2 PacbPORF objects by introns

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
    _update_kwargs(kwargs,KWARGS_MAPPED_INTRON)
    if not kwargs.has_key('aligned_site_max_triplet_distance'):
        kwargs['aligned_site_max_triplet_distance'] = kwargs['max_aa_offset']

    # settings for minimal alignment entropy score
    min_donor_site_alignment_entropy = 0.0
    min_acceptor_site_alignment_entropy = 0.0

    # calculate maximal/minimal donor/acceptor site position based on alignment
    ELEGIABLE_SPLICE_SITE_AA_RANGE = 75

    qdr = pacbporfD.alignment_dna_range_query()
    qar = pacbporfA.alignment_dna_range_query()
    min_donor_query_pos = max([ min(qdr), max(qdr)-(ELEGIABLE_SPLICE_SITE_AA_RANGE*3) ])
    max_accep_query_pos = min([ max(qar), min(qar)+(ELEGIABLE_SPLICE_SITE_AA_RANGE*3) ])

    sdr = pacbporfD.alignment_dna_range_sbjct()
    sar = pacbporfA.alignment_dna_range_sbjct()
    min_donor_sbjct_pos = max([ min(sdr), max(sdr)-(ELEGIABLE_SPLICE_SITE_AA_RANGE*3) ])
    max_accep_sbjct_pos = min([ max(sar), min(sar)+(ELEGIABLE_SPLICE_SITE_AA_RANGE*3) ])

    # get list of introns
    #intronsQ = merge_orfs_with_intron(pacbporfD.orfQ,pacbporfA.orfQ,
    #        min_donor_pos   =min_donor_query_pos,
    #        max_acceptor_pos=max_accep_query_pos,**kwargs)
    #intronsS = merge_orfs_with_intron(pacbporfD.orfS,pacbporfA.orfS,
    #        min_donor_pos   =min_donor_sbjct_pos,
    #        max_acceptor_pos=max_accep_sbjct_pos,**kwargs)

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
        print "dQ1", [ d.pos for d in donorQ ], "aQ1", [ a.pos for a in accepQ ]
        print "dS1", [ d.pos for d in donorS ], "aS1", [ a.pos for a in accepS ]
    ############################################################################

    # filter for alignable donor & acceptor sites
    kwargs['allow_non_canonical'] = kwargs['allow_non_canonical_donor']
    algdonors = _filter_for_alignable_splice_sites(donorQ,donorS,pacbporfD,**kwargs)
    kwargs['allow_non_canonical'] = kwargs['allow_non_canonical_acceptor']
    algacceps = _filter_for_alignable_splice_sites(accepQ,accepS,pacbporfA,**kwargs)

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

    # now make list of aligable introns
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

# end of function merge_pacbporfs_with_introns


def _filter_aligned_stopless_3n_introns(algintrons,verbose=False):
    """ """
    # The PacbPORFs around stopless introns are likely to be somewhat
    # extended on their tails (more sequence to be possibly alignable ;-)
    # As a result, alignment entropy positioning is useless (0.0)
    # This function filters for shorter stopless 3n introns at the cost
    # of longer 3n introns.

    # return empty list if not aligned introns were found
    if not algintrons: return algintrons

    # calculate summed PSSM score of both introns
    for pos in range(0,len(algintrons)):
        intQ,intS = algintrons[pos]
        summed_pssm = intQ.donor.pssm_score + intQ.acceptor.pssm_score +\
                      intS.donor.pssm_score + intS.acceptor.pssm_score
        algintrons[pos] = (summed_pssm,intQ,intS)
    # order by summed_pssm
    algintrons.sort()
    algintrons.reverse()
    best_pssm_score,bestQ,bestS = algintrons[0]
    best_donor_entropy = max([bestQ._apps_donor,bestS._apps_donor])
    best_accep_entropy = max([bestQ._apps_accep,bestS._apps_accep])
    if bestQ.is_stopless_3n_intron() and bestS.is_stopless_3n_intron():
        best_3n_length = max([bestQ.length,bestS.length])
        stopless3n = "QS"
    elif bestQ.is_stopless_3n_intron():
        best_3n_length = bestQ.length
        stopless3n = "Q"
    elif bestS.is_stopless_3n_intron():
        best_3n_length = bestS.length
        stopless3n = "S"
    else:
        # no stopless introns! -> function was called accidentially !?
        # translate list back to list of introns only (remove summed_pssm)
        algintrons = [ (intQ,intS) for (summed_pssm,intQ,intS) in algintrons ]
        # and return the list of (unfiltered) algintrons
        return algintrons

    # select aligned introns to be removed from the list of candidates
    # based on lower summed PSSM AND/OR entropy <= 0.0
    remove_list_index_ids = []
    for pos in range(1,len(algintrons)):
        summed_pssm,intQ,intS = algintrons[pos]
        if summed_pssm < best_pssm_score:
            if stopless3n == "Q" and best_3n_length < intQ.length:
                # longer stopless 3n query intron -> remove
                remove_list_index_ids.append(pos)
            elif stopless3n == "S" and best_3n_length < intS.length:
                # longer stopless 3n sbjct intron -> remove
                remove_list_index_ids.append(pos)
            elif stopless3n == "QS" and best_3n_length < max([intQ.length,intS.length]):
                # longer stopless 3n query/sbjct intron -> remove
                remove_list_index_ids.append(pos)
            else:
                pass

    # remove aligned introns that are removed based on stopless 3n criterion
    remove_list_index_ids.reverse()
    for pos in remove_list_index_ids: algintrons.pop(pos)
    # translate list back to list of introns only (remove summed_pssm)
    algintrons = [ (intQ,intS) for (summed_pssm,intQ,intS) in algintrons ]

    # return lists of aligned introns
    return algintrons

# end of function _filter_aligned_stopless_3n_introns


def _filter_aligned_introns_on_pssm_entropy_combination(algintrons,verbose=False):
    """
    """
    # return empty list if not aligned introns were found
    if not algintrons: return algintrons

    # calculate summed PSSM score of both introns
    for pos in range(0,len(algintrons)):
        intQ,intS = algintrons[pos]
        summed_pssm = intQ.donor.pssm_score + intQ.acceptor.pssm_score +\
                      intS.donor.pssm_score + intS.acceptor.pssm_score
        algintrons[pos] = (summed_pssm,intQ,intS)
    # order by summed_pssm
    algintrons.sort()
    algintrons.reverse()
    best_pssm_score,bestQ,bestS = algintrons[0]
    best_donor_entropy = max([bestQ._apps_donor,bestS._apps_donor])
    best_accep_entropy = max([bestQ._apps_accep,bestS._apps_accep])


    # select aligned introns to be removed from the list of candidates
    # based on lower summed PSSM AND/OR entropy <= 0.0
    remove_list_index_ids = []
    for pos in range(1,len(algintrons)):
        summed_pssm,intQ,intS = algintrons[pos]
        donor_entropy = max([intQ._apps_donor,intS._apps_donor])
        accep_entropy = max([intQ._apps_accep,intS._apps_accep])
        if summed_pssm < (best_pssm_score * 0.8):
            if intQ._apps_donor <= 0.0 or intQ._apps_accep <= 0.0 or\
            intS._apps_donor <= 0.0 or intS._apps_accep <= 0.0:
                # introns have both poor PSSM and poor entropy
                remove_list_index_ids.append(pos)
        elif summed_pssm < best_pssm_score:
            if best_donor_entropy > 0.0 and donor_entropy <= 0.0:
                # introns have okay PSSM but (very) poor donor entropy
                remove_list_index_ids.append(pos)
            elif best_accep_entropy > 0.0 and accep_entropy <= 0.0:
                # introns have okay PSSM but (very) poor acceptor entropy
                remove_list_index_ids.append(pos)
            else:
                pass
        else:
            pass
    # remove aligned introns
    remove_list_index_ids.reverse()
    for pos in remove_list_index_ids: algintrons.pop(pos)
    # translate list back to list of introns only (remove summed_pssm)
    algintrons = [ (intQ,intS) for (summed_pssm,intQ,intS) in algintrons ]

    ########################################################################
    if verbose:
        for (intQ,intS) in algintrons:
            print "Filtered introns:", ( intQ.donor.pos, intQ.acceptor.pos ) ,
            print ( intS.donor.pos, intS.acceptor.pos ),
            print "ENTROPY: %1.2f %1.2f" % (intQ._apps_donor, intQ._apps_accep),
            print "PSSM: (%1.2f %1.2f) (%1.2f %1.2f)" % (
                intQ.donor.pssm_score, intS.donor.pssm_score,
                intQ.acceptor.pssm_score, intS.acceptor.pssm_score,
                )
    ########################################################################

    # return lists of aligned introns
    return algintrons

# end of function _filter_aligned_introns_on_pssm_entropy_combination


def _tinyexon_list_2_dict(resultlist):
    """ """
    resultdict = {}
    key2exon = {}
    for (intronD,tinyexon,intronA) in resultlist:
        exon_key = (tinyexon.orf.id, tinyexon.donor.pos, tinyexon.acceptor.pos)
        if not resultdict.has_key(exon_key):
            resultdict[exon_key] = [ {}, {} ]
            key2exon[exon_key] = tinyexon
        intronDkey = (intronD.orfDonor.id,intronD.donor.pos)
        intronAkey = (intronA.orfAcceptor.id,intronA.acceptor.pos)
        if not resultdict[exon_key][0].has_key(intronDkey):
            resultdict[exon_key][0][intronDkey] = intronD
        if not resultdict[exon_key][1].has_key(intronAkey):
            resultdict[exon_key][1][intronAkey] = intronA
    return resultdict, key2exon
# end of function _tinyexon_list_2_dict


def merge_pacbporfs_by_tinyexons(pacbporfD,pacbporfA,
    orfSetObjQ,orfSetObjS,verbose=False,**kwargs):
    """ """
    # input validation
    IsPacbPORF(pacbporfD)
    IsPacbPORF(pacbporfA)

    # edit **kwargs dictionary for some forced attributes
    _update_kwargs(kwargs,KWARGS_MAPPED_INTRON)
    if not kwargs.has_key('aligned_site_max_triplet_distance'):
        kwargs['aligned_site_max_triplet_distance'] = kwargs['max_aa_offset']

    # settings for minimal alignment entropy score
    min_donor_site_alignment_entropy = 0.0
    min_acceptor_site_alignment_entropy = 0.0

    resultlistQ = merge_orfs_with_tinyexon(
            pacbporfD.orfQ,pacbporfA.orfQ,
            preceding_donor_sites=pacbporfD.orfQ._donor_sites,
            subsequent_acceptor_sites=pacbporfA.orfQ._acceptor_sites,
            orflist=orfSetObjQ.orfs,**kwargs)
    resultlistS = merge_orfs_with_tinyexon(
            pacbporfD.orfS,pacbporfA.orfS,
            preceding_donor_sites=pacbporfD.orfS._donor_sites,
            subsequent_acceptor_sites=pacbporfA.orfS._acceptor_sites,
            orflist=orfSetObjS.orfs,**kwargs)

    # translate resultlists to dict: key == exon, value = [ {intronsD},{intronsS} ]
    resultdictQ,key2exonQ = _tinyexon_list_2_dict(resultlistQ)
    resultdictS,key2exonS = _tinyexon_list_2_dict(resultlistS)

    # get unique list of donors & acceptors
    donorQ = olba( list(Set([inD.donor for inD,te,inA in resultlistQ ])), order_by='pos')
    donorS = olba( list(Set([inD.donor for inD,te,inA in resultlistS ])), order_by='pos')
    accepQ = olba( list(Set([inA.acceptor for inD,te,inA in resultlistQ ])), order_by='pos')
    accepS = olba( list(Set([inA.acceptor for inD,te,inA in resultlistS ])), order_by='pos')

    ## filter for alignable donor & acceptor sites
    kwargs['allow_non_canonical']               = True # True
    kwargs['aligned_site_max_triplet_distance'] = 0     # 2
    algdonors = _filter_for_alignable_splice_sites(donorQ,donorS,pacbporfD,**kwargs)
    algacceps = _filter_for_alignable_splice_sites(accepQ,accepS,pacbporfA,**kwargs)

    # settings for minimal alignment entropy score
    # TODO TODO -> THIS MUST BE FIXED TO A NICE THRESHOLD VALUE!!!
    min_donor_site_alignment_entropy = 0.1
    min_acceptor_site_alignment_entropy = 0.1


    # remove sites with to low alignment entropy
    algdonors = _filter_for_entropy(algdonors,pacbporfD,'donor',
                min_alignment_entropy=min_donor_site_alignment_entropy)
    algacceps = _filter_for_entropy(algacceps,pacbporfA,'acceptor',
                min_alignment_entropy=min_acceptor_site_alignment_entropy)

    # return list: intronQD,intronSD,tinyexon,intronAQ,intronAS
    return_list = []

    ############################################################################
    if verbose:
        print "bridges constructed: ORFS:",
        print (pacbporfD.orfQ.id,pacbporfA.orfQ.id),
        print (pacbporfD.orfS.id,pacbporfA.orfS.id),
        print len(resultdictQ), len(resultdictS),
        print ( len(resultlistQ), len(donorQ), len(accepQ) ),
        print ( len(resultlistS), len(donorS), len(accepS) ),
        print ( len(algdonors), len(algacceps) )
    ############################################################################

    for keyQ,tinyexonQ in key2exonQ.iteritems():
        for keyS,tinyexonS in key2exonS.iteritems():
            if tinyexonQ.donor.phase != tinyexonS.donor.phase:
                continue
            if tinyexonQ.acceptor.phase != tinyexonS.acceptor.phase:
                continue
            if tinyexonQ.length != tinyexonS.length:
                continue
            # if here, then tinyexons of identical structure


            ####################################################################
            if verbose:
                print tinyexonQ.length, tinyexonQ.donor.phase,
                print ( len(resultdictQ[keyQ][0]), len(resultdictQ[keyQ][1]) ),
                print ( len(resultdictS[keyS][0]), len(resultdictS[keyS][1]) ),
                print tinyexonQ,
                print tinyexonQ.proteinsequence(), tinyexonS.proteinsequence(),
                print tinyexonS.acceptor.pssm_score + tinyexonS.donor.pssm_score
            ####################################################################

            donor_introns = []
            acceptor_introns = []
            for intronDQkey, intronDQ in resultdictQ[keyQ][0].iteritems():
                if intronDQ.donor.pos not in [ dQ.pos for dQ,dS in algdonors ]:
                    continue
                for intronDSkey, intronDS in resultdictS[keyS][0].iteritems():
                    if intronDS.donor.pos not in [ dS.pos for dQ,dS in algdonors ]:
                        continue
                    # check if they exists as aligned sites
                    alignedkey = ( intronDQ.donor.pos, intronDS.donor.pos )
                    if alignedkey not in [ (dQ.pos, dS.pos) for dQ,dS in algdonors ]:
                        continue
                    # if here, we have a set of introns 5' of the tinyexon
                    # which are perfectly alignable!
                    donor_introns.append((intronDQ,intronDS))

            for intronAQkey, intronAQ in resultdictQ[keyQ][1].iteritems():
                if intronAQ.acceptor.pos not in [ aQ.pos for aQ,aS in algacceps ]:
                    continue
                for intronASkey, intronAS in resultdictS[keyS][1].iteritems():
                    if intronAS.acceptor.pos not in [ aS.pos for aQ,aS in algacceps ]:
                        continue
                    # check if they exists as aligned sites
                    alignedkey = ( intronAQ.acceptor.pos, intronAS.acceptor.pos )
                    if alignedkey not in [ (aQ.pos, aS.pos) for aQ,aS in algacceps ]:
                        continue
                    # if here, we have a set of introns 3' of the tinyexon
                    # which are perfectly alignable!
                    acceptor_introns.append((intronAQ,intronAS))

            if not len(donor_introns) or not len(acceptor_introns):
                # no aligned 5' && aligned 3' introns
                continue

            # initialize extended tinyexon PacbPORF
            from pacb import PacbP
            pacbp = PacbP(input=( 
                    tinyexonQ.proteinsequence(),
                    tinyexonS.proteinsequence(),
                    tinyexonQ.protein_start(),
                    tinyexonS.protein_start(),
                    ) )
            pacbp.strip_unmatched_ends()
            # continue if no fraction could be aligned
            if len(pacbp) == 0: continue
            tinypacbporf = pacbp2pacbporf(pacbp,tinyexonQ.orf,tinyexonS.orf)
            tinypacbporf.extend_pacbporf_after_stops()

            ####################################################################
            if verbose:
                print tinypacbporf
                tinypacbporf.print_protein_and_dna()
                print len(donor_introns), len(acceptor_introns),
                print max([ dQ.donor.pssm_score+dS.donor.pssm_score for dQ,dS in donor_introns]),
                print max([ aQ.acceptor.pssm_score+aS.acceptor.pssm_score for aQ,aS in acceptor_introns])
            ####################################################################


            # if here, we have accepted tinyexon bridges!
            # gather them and store to return_list
            for intronDQkey, intronDQ in resultdictQ[keyQ][0].iteritems():
                if intronDQ.donor.pos not in [ dQ.pos for dQ,dS in algdonors ]:
                    continue
                for intronDSkey, intronDS in resultdictS[keyS][0].iteritems():
                    if intronDS.donor.pos not in [ dS.pos for dQ,dS in algdonors ]:
                        continue
                    for intronAQkey, intronAQ in resultdictQ[keyQ][1].iteritems():
                        if intronAQ.acceptor.pos not in [ aQ.pos for aQ,aS in algacceps ]:
                            continue
                        for intronASkey, intronAS in resultdictS[keyS][1].iteritems():
                            if intronAS.acceptor.pos not in [ aS.pos for aQ,aS in algacceps ]:
                                continue
                            ####################################################
                            # set some meta-data properties to the intron objects
                            ####################################################
                            _score_introns_obtained_by_mapping(
                                    intronDQ,intronDS,pacbporfD,
                                    tinypacbporf,source='ABGPmappingTE')
                            _score_introns_obtained_by_mapping(
                                    intronAQ,intronAS,tinypacbporf,
                                    pacbporfA,source='ABGPmappingTE')
                            # create _linked_to_xxx attributes
                            intronDQ._linked_to_pacbporfs = [ tinypacbporf ]
                            intronAQ._linked_to_pacbporfs = [ tinypacbporf ]
                            intronDS._linked_to_pacbporfs = [ tinypacbporf ]
                            intronAS._linked_to_pacbporfs = [ tinypacbporf ]
                            intronDQ._linked_to_introns   = [ intronAQ ]
                            intronAQ._linked_to_introns   = [ intronDQ ]
                            intronDS._linked_to_introns   = [ intronAS ]
                            intronAS._linked_to_introns   = [ intronDS ]
                            # append to tmp result list
                            return_list.append(
                                (intronDQ,intronDS,tinypacbporf,intronAQ,intronAS)
                                )

    # check if there are >1 candidate tiny exons
    # currently, we choose only to return the **best** mapped tinyexon 
    if len(return_list) == 0:
        pass
    elif len(return_list) == 1:
        pass
    else:
        # only take the highest scoring candidate here 
        min_distance = min([ (a._distance+d._distance) for a,b,c,d,e in return_list ])
        pos2score = []
        for (intronDQ,intronDS,tinypacbporf,intronAQ,intronAS) in return_list:
            if (intronDQ._distance + intronAQ._distance) > min_distance:
                pos2score.append( 0.0 )
            else:
                # calculate overall pssm score
                total_pssm = 0.0
                total_pssm += intronDQ.donor.pssm_score
                total_pssm += intronDQ.acceptor.pssm_score
                total_pssm += intronDS.donor.pssm_score
                total_pssm += intronDS.acceptor.pssm_score
                total_pssm += intronAQ.donor.pssm_score
                total_pssm += intronAQ.acceptor.pssm_score
                total_pssm += intronAS.donor.pssm_score
                total_pssm += intronAS.acceptor.pssm_score
                pos2score.append( total_pssm )
        # get highest score and linked tinyexon
        max_score = max(pos2score)
        return_list = [ return_list[pos2score.index(max_score)] ]

    ############################################################################
    # some printing in verbose mode
    if verbose and return_list:
        (intronDQ,intronDS,tinypacbporf,intronAQ,intronAS) = return_list[0]
        print "BEST MAPPED TINYEXON:"
        print tinypacbporf
        print tinypacbporf.query, intronDQ._distance, intronAQ._distance,
        print ( intronDQ.donor.pos, intronDQ.acceptor.pos ),
        print ( intronDS.donor.pos, intronDS.acceptor.pos ),
        print ( intronAQ.donor.pos, intronAQ.acceptor.pos ),
        print ( intronAS.donor.pos, intronAS.acceptor.pos )
    ############################################################################

    # return the result list
    return return_list

# end of function merge_pacbporfs_by_tinyexons


def merge_pacbporfs_with_closeby_independant_introns(pacbporfD,pacbporfA,
    verbose=False,**kwargs):
    """
    Merge 2 PacbPORF objects by closeby independant gained introns

    @attention: see pacb.connecting.merge_orfs_with_intron for **kwargs)

    @type  pacbporfD: PacbPORF object
    @param pacbporfD: PacbPORF object that has to deliver PSSM donor objects

    @type  pacbporfA: PacbPORF object
    @param pacbporfA: PacbPORF object that has to deliver PSSM acceptor objects

    @type  verbose: Boolean
    @param verbose: print status/debugging messages to STDOUT

    @rtype:  list
    @return: list with ( intronQ, intronS, CIGexonPacbPORF )
    """
    # input validation
    IsPacbPORF(pacbporfD)
    IsPacbPORF(pacbporfA)

    # edit **kwargs dictionary for some forced attributes
    kwargs['allow_phase_shift'] = True
    _update_kwargs(kwargs,KWARGS_CLOSEBY_INDEPENDANT_INTRON_GAIN)
    if not kwargs.has_key('aligned_site_max_triplet_distance'):
        kwargs['aligned_site_max_triplet_distance'] = kwargs['cig_max_aa_length']

    # run regular merge_pacbporfs_with_introns function
    alg_introns = merge_pacbporfs_with_introns(pacbporfD,pacbporfA,verbose=verbose,**kwargs)
    cig_introns = []

    if verbose:
        print "introns::", len(alg_introns), "cig_max_aa_length:", kwargs['cig_max_aa_length'], kwargs['aligned_site_max_triplet_distance']

    # check if there is length congruence between the cig_introns
    for intQ,intS in alg_introns:
        dQpos, dQphase = pacbporfD.dnaposition_query(intQ.donor.pos,forced_return=True)
        dSpos, dSphase = pacbporfD.dnaposition_sbjct(intS.donor.pos,forced_return=True)
        aQpos, aQphase = pacbporfA.dnaposition_query(intQ.acceptor.pos,forced_return=True)
        aSpos, aSphase = pacbporfA.dnaposition_sbjct(intS.acceptor.pos,forced_return=True)
        distDnt = (dQpos*3 + dQphase) - (dSpos*3 + dSphase)
        distAnt = (aQpos*3 + aQphase) - (aSpos*3 + aSphase)
        ########################################################################
        if verbose:
            print (intQ.donor.pos, intQ.acceptor.pos),
            print (intS.donor.pos, intS.acceptor.pos),
            print distDnt, distAnt, kwargs['max_nt_offset']
        ########################################################################
        if abs(distDnt-distAnt) > kwargs['max_nt_offset']:
            # intermediate ciigPacbPORF has query vs sbjct length discrepancy
            # *3 for AA2nt coordinate conversion, +2 to allow different phases
            # e.g. phase difference can give 1AA+2nt difference
            continue
        if intQ.donor.phase == intS.donor.phase and\
        (distDnt/3) <= kwargs['aligned_site_max_triplet_distance']:
            # a regularly merged intron combination
            continue
        if intQ.acceptor.phase == intS.acceptor.phase and\
        (distAnt/3) <= kwargs['aligned_site_max_triplet_distance']:
            # a regularly merged intron combination
            continue
        if abs(distDnt) <= 5 or abs(distDnt) <= 5:
            # most likely a splice site phase shift, not a c.i.g.
            continue

        if abs(distDnt/3) >= kwargs['cig_min_aa_length'] and\
        abs(distAnt/3) >= kwargs['cig_min_aa_length'] and\
        abs(distDnt/3) <= kwargs['cig_max_aa_length'] and\
        abs(distAnt/3) <= kwargs['cig_max_aa_length']:
            # putatively a closeby independant (intron) gain
            cig_introns.append( ( intQ, intS ) )

    ############################################################################
    if verbose:
        for intQ,intS in cig_introns:
            print "cig?:", (intQ.donor.pos, intQ.acceptor.pos),
            print (intS.donor.pos, intS.acceptor.pos)
    ############################################################################


    # return variable to store found positive cases of CIG into
    found_cig_list = []

    # check if there is some sequence similarity
    for intQ,intS in cig_introns:
        # get alignment positions around query & sbjcts splice sites
        dQpos, dQphase = pacbporfD.dnaposition_query(intQ.donor.pos,forced_return=True)
        dSpos, dSphase = pacbporfD.dnaposition_sbjct(intS.donor.pos,forced_return=True)
        aQpos, aQphase = pacbporfA.dnaposition_query(intQ.acceptor.pos,forced_return=True)
        aSpos, aSphase = pacbporfA.dnaposition_sbjct(intS.acceptor.pos,forced_return=True)
        distD = dQpos - dSpos
        distA = aQpos - aSpos
        distDnt = (dQpos*3 + dQphase) - (dSpos*3 + dSphase)
        distAnt = (aQpos*3 + aQphase) - (aSpos*3 + aSphase)

        if distDnt > 0:   # then, distAnt is as well > 0
            # QUERY is extended on the donor side
            #mode   = "SQ"
            #qStart = pacbporfD._positions[dSpos].query_pos
            #qEnd   = qStart + distD
            #sStart = pacbporfA._positions[aSpos].sbjct_pos
            #sEnd   = sStart + distD
            #qSeq = pacbporfD.orfQ.getaas(abs_pos_start=qStart,abs_pos_end=qEnd)
            #sSeq = pacbporfA.orfS.getaas(abs_pos_start=sStart,abs_pos_end=sEnd)
            mode  = "SQ"
            qEnd  = pacbporfD.orfQ.dnapos2aapos(intQ.donor.pos)
            qStart= qEnd - max([distA,distD])
            sStart= pacbporfA.orfS.dnapos2aapos(intS.acceptor.pos)
            sEnd  = sStart + max([distA,distD])
            qSeq  = pacbporfD.orfQ.getaas(abs_pos_start=qStart,abs_pos_end=qEnd)
            sSeq  = pacbporfA.orfS.getaas(abs_pos_start=sStart,abs_pos_end=sEnd)

        else: # distDnt and distAnt are < 0
            ## SBJCT is extended on the donor site
            #mode   = "QS"
            #qStart = pacbporfA._positions[aQpos].query_pos
            #qEnd   = qStart - distA
            #sStart = pacbporfD._positions[dQpos].sbjct_pos
            #sEnd   = sStart - distA
            #qSeq = pacbporfA.orfQ.getaas(abs_pos_start=qStart, abs_pos_end=qEnd)
            #sSeq = pacbporfD.orfS.getaas(abs_pos_start=sStart, abs_pos_end=sEnd)
            mode  = "QS"
            qStart= pacbporfA.orfQ.dnapos2aapos(intQ.acceptor.pos)
            qEnd  = qStart - min([distA,distD])
            sEnd  = pacbporfD.orfS.dnapos2aapos(intS.donor.pos)
            sStart= sEnd + min([distA,distD])
            qSeq  = pacbporfA.orfQ.getaas(abs_pos_start=qStart,abs_pos_end=qEnd)
            sSeq  = pacbporfD.orfS.getaas(abs_pos_start=sStart,abs_pos_end=sEnd)


        headerQ = "query_%s_%s_%s" % (qStart,qEnd,qSeq)
        headerS = "sbjct_%s_%s_%s" % (sStart,sEnd,sSeq)
        headerQ = headerQ[0:20] # truncate to prevent error
        headerS = headerS[0:20] # truncate to prevent error
        if verbose:
            print mode, (distD,distA), qSeq, sSeq, headerQ, headerS, distDnt, distAnt,
            print dQpos, aQpos, dSpos, aSpos
        if not qSeq: continue # superfluous check-doublecheck for sequence
        if not sSeq: continue # superfluous check-doublecheck for sequence

        ####################################################
        # make PacbPORF with ClustalW
        ####################################################
        # align the sequences with clustalw
        seqs = { headerQ: qSeq, headerS: sSeq }
        (alignedseqs,alignment) = clustalw(seqs=seqs)

        # make pacbp from clustalw alignment
        pacbp = pacbp_from_clustalw(
                    alignment=(
                            alignedseqs[headerQ],
                            alignment,
                            alignedseqs[headerS]
                            ),
                    coords=(qStart,qEnd,sStart,sEnd)
                    )

        if not pacbp: continue

        # strip unaligned fraction of this pacbp object, then check length
        pacbp.strip_unmatched_ends()

        if len(pacbp) < kwargs['cig_min_aa_length']:
            continue
        if len(pacbp) > kwargs['cig_max_aa_length']:
            continue

        if pacbp:
            # initialize extended tiny PacbPORF caused by c.i.g.
            if distDnt > 0:
                cig_pacbporf = pacbp2pacbporf(pacbp,pacbporfD.orfQ,pacbporfA.orfS)
            else:
                cig_pacbporf = pacbp2pacbporf(pacbp,pacbporfA.orfQ,pacbporfD.orfS)
            cig_pacbporf.extend_pacbporf_after_stops()
            ####################################################################
            if verbose:
                print pacbp, len(pacbp)
                print cig_pacbporf
                print "CIG:", intQ
                print "CIG:", intS
                print distD, distA, distDnt, distAnt
                cig_pacbporf.print_protein_and_dna()
            ####################################################################

            ####################################################################
            # set some meta-data properties to the intron objects
            ####################################################################


            # add distance score to introns
            # The distance set in merge_pacbporfs_with_introns is large;
            # it is the actual distance between the splice sites. In CIG,
            # the measure for distance is the length difference between
            # the offset between query and sbjct measured on the cig_pacbporf
            intQ._distance = abs(distDnt-distAnt)
            intS._distance = abs(distDnt-distAnt)
    
            if distDnt > 0:   # then, distAnt is as well > 0
                # QUERY is extended on the donor side
                # add Alignment Positional Periphery Score into objects
                succes = set_apps_intron_query(intQ,cig_pacbporf,pacbporfA)
                succes = set_apps_intron_sbjct(intS,pacbporfD,cig_pacbporf)
            else:
                # SBJCT is extended on the donor side
                # add Alignment Positional Periphery Score into objects
                succes = set_apps_intron_query(intQ,pacbporfD,cig_pacbporf)
                succes = set_apps_intron_sbjct(intS,cig_pacbporf,pacbporfA)

            # set GFF fsource attribute for recognition of intron sources
            intQ._gff['fsource'] = "ABGPcig"
            intS._gff['fsource'] = "ABGPcig"

            # create _linked_to_xxx attributes
            intQ._linked_to_pacbporfs = [ cig_pacbporf ]
            intS._linked_to_pacbporfs = [ cig_pacbporf ]


            # append to found_cig_list
            found_cig_list.append( ( intQ, intS, cig_pacbporf ) )

        else:
            # no alignment possible -> try next
            continue
    
    # return lists of closeby_independant_introns
    return found_cig_list

# end of function merge_pacbporfs_with_closeby_independant_introns


def merge_pacbporfs_with_phase_shift_introns(pacbporfD,pacbporfA,
    verbose=False,**kwargs):
    """
    Merge 2 PacbPORF objects by introns of which one underwent a phase shift

    @attention: see pacb.connecting.merge_orfs_with_intron for **kwargs)

    @type  pacbporfD: PacbPORF object
    @param pacbporfD: PacbPORF object that has to deliver PSSM donor objects

    @type  pacbporfA: PacbPORF object
    @param pacbporfA: PacbPORF object that has to deliver PSSM acceptor objects

    @type  verbose: Boolean
    @param verbose: print status/debugging messages to STDOUT

    @rtype:  list
    @return: list with ( intronQ, intronS, CIGexonPacbPORF )
    """
    # input validation
    IsPacbPORF(pacbporfD)
    IsPacbPORF(pacbporfA)

    # edit **kwargs dictionary
    kwargs['allow_phase_shift'] = True
    _update_kwargs(kwargs,KWARGS_PHASE_SHIFT_INTRON)
    if not kwargs.has_key('aligned_site_max_triplet_distance'):
        kwargs['aligned_site_max_triplet_distance'] = kwargs['max_aa_distance']

    # run regular merge_pacbporfs_with_introns function
    alg_introns = merge_pacbporfs_with_introns(pacbporfD,pacbporfA,**kwargs)
    psh_introns = []

    # check if there is length congruence between the cig_introns
    for intQ,intS in alg_introns:
        # check phase equilibrium -> if equal, no phase shift
        if intQ.donor.phase == intS.donor.phase:
            continue

        ########################################################################
        # set some meta-data properties to the intron objects
        # attribute _distance is already set in merge_pacbporfs_with_introns
        # attribute(s) ~APPS are already set in merge_pacbporfs_with_introns
        ########################################################################

        # set GFF fsource attribute for recognition of intron sources
        intQ._gff['fsource'] = "ABGPphs"
        intS._gff['fsource'] = "ABGPphs"

        # putatively a phase shifted intron pair
        psh_introns.append( ( intQ, intS ) )
    
    # return lists of phase shifted introns
    return psh_introns

# end of function merge_pacbporfs_with_phase_shift_introns


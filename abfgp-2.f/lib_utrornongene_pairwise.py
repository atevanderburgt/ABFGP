"""
"""

# Module metadata
__authors__ = "Ate van der Burgt"
__license__ = "MIT"

# ABGP Imports
from graphAbgp import (
    PacbpCollectionGraph,
    InwardsPointingCodingBlockGraph,
    )
from graphAbgp.graph_pacbpcollection import _delete_pacbp


from lib_synteny_pairwise import get_first_and_final_inwpcbg_pos
from lib_pcg2blocks import (
    PCG2inwpCBGS,
    print_inwpcbgstructure,
    pcg2gtg_by_bitscore,
    pcg2gtg_by_identity,
    )

# Python Imports
from sets import Set

# Global Variable Imports
NONGENE_GTG_MAX_DIFFERENCE      = -0.475
NONGENE_GTG_MAX_MSR_RATIO       = 0.33
NONGENE_GTG_MAX_GTG_RATIO       = 0.85
NONCODINGNONGENE_5p_INWPCBG_MAX_IDENTITYRATIO = 0.75
NONCODINGNONGENE_3p_INWPCBG_MAX_IDENTITYRATIO = 0.75

from settings.genestructure import (
    MIN_INTERGENIC_NT_LENGTH,
    MAX_INTERGENIC_MIN_NT_LENGTH,
    AVERAGE_INTERGENIC_MIN_NT_LENGTH,
    )
from settings.executables import TCODE_MAX_NONCODING


from settings.translationalstartsites import TSS_MIN_PSSM_SCORE



def detect_and_remove_single_nonfinal_inwpcbg(inwpcbgs,PCG,GENE_IDENTIFIER_SET,
    verbose=False):
    """
    Allow deletion of a very shitty, single inwpCBG from the end of the list
    """
    # we need at least 2 inwpCBGs in order to remove one of them
    if len(inwpcbgs) <= 1: return False

    lastInwpCBG = inwpcbgs[-1]
    prevInwpCBG = inwpcbgs[-2]

    lastNodeList = [ lastInwpCBG.get_organism_nodes(org)[0] for org in\
                lastInwpCBG.organism_set().intersection(GENE_IDENTIFIER_SET) ]
    prevNodeList = [ prevInwpCBG.get_organism_nodes(org)[0] for org in\
                prevInwpCBG.organism_set().intersection(GENE_IDENTIFIER_SET) ]

    # identical nodes -> do not delete. Only go for very obvious things
    if Set(lastNodeList).intersection(prevNodeList): return False

    ntdistdict    = prevInwpCBG.nt_spacing_between_codingblocks([lastInwpCBG])
    tcodedistdict = prevInwpCBG.tcode_spacing_between_codingblocks([lastInwpCBG])

    check1 = prevInwpCBG.count_orfs_labeled_as_annotated_exon() >\
             lastInwpCBG.count_orfs_labeled_as_annotated_exon() 
    check2 = prevInwpCBG.get_bitscore() > lastInwpCBG.get_bitscore() 
    check3 = len(prevNodeList) > len(lastNodeList)
    check4 = float(lastInwpCBG.count_orfs_labeled_as_annotated_exon()) /\
             float(len(GENE_IDENTIFIER_SET)) <= 0.33
    if ntdistdict:
        check5 = sum(ntdistdict.values())/float(len(ntdistdict)) >\
             MIN_INTERGENIC_NT_LENGTH
    else:
        check5 = False
    if tcodedistdict:
        check6 = sum(tcodedistdict.values())/float(len(tcodedistdict)) <\
             TCODE_MAX_NONCODING
    else:
        check6 = False
    check7 = prevInwpCBG.get_projected_tailing_stop_aa_difference() <\
             lastInwpCBG.get_projected_tailing_stop_aa_difference()
    check8 = prevInwpCBG.get_projected_tailing_stop_nonaligned_aa_difference()<\
             lastInwpCBG.get_projected_tailing_stop_nonaligned_aa_difference()

    checklist = [check1,check2,check3,check4,check5,check6,check7,check8]

    ############################################################################
    if verbose: print "NonFinal inwpCBG check:", checklist
    ############################################################################
   
    if checklist.count(False) == 0:
        nonfinalPCG = PacbpCollectionGraph(crossdata={},blastmatrix=PCG._blastmatrix)
        # place all PacbPORFs in the nonfinalPCG
        for (pacbpkey,nodeQ,nodeS), pacbporf in lastInwpCBG.pacbps.iteritems():
            # add to noncodingnongenePCG
            nonfinalPCG.add_node(nodeQ)
            nonfinalPCG.add_node(nodeS)
            nonfinalPCG.add_edge(nodeQ,nodeS,wt=pacbporf.bitscore)
            nonfinalPCG.pacbps[(pacbpkey,nodeQ,nodeS)] = pacbporf
            # remove from main PCG
            _delete_pacbp(PCG,(pacbpkey,nodeQ,nodeS))
        # return nonfinalPCG
        return nonfinalPCG
    else:
        return False

# end of function detect_and_remove_single_nonfinal_inwpcbg


def detect_and_remove_single_nonfirst_inwpcbg(inwpcbgs,PCG,GENE_IDENTIFIER_SET,
    verbose=False):
    """
    Allow deletion of a very shitty, single inwpCBG from the start of the list
    """
    # we need at least 2 inwpCBGs in order to remove one of them
    if len(inwpcbgs) <= 1: return False

    firstInwpCBG = inwpcbgs[0]
    nextInwpCBG = inwpcbgs[1]

    firstNodeList = [ firstInwpCBG.get_organism_nodes(org)[0] for org in\
                firstInwpCBG.organism_set().intersection(GENE_IDENTIFIER_SET) ]
    nextNodeList = [ nextInwpCBG.get_organism_nodes(org)[0] for org in\
                nextInwpCBG.organism_set().intersection(GENE_IDENTIFIER_SET) ]

    # identical nodes -> do not delete. Only go for very obvious things
    if Set(firstNodeList).intersection(nextNodeList): return False

    ntdistdict    = firstInwpCBG.nt_spacing_between_codingblocks([nextInwpCBG])
    tcodedistdict = firstInwpCBG.tcode_spacing_between_codingblocks([nextInwpCBG])

    # make a long list of checks which should be True in case
    # firstInwpCBG is *NOT* the first exon of this gene structure
    check1 = nextInwpCBG.count_orfs_labeled_as_annotated_exon() >\
             firstInwpCBG.count_orfs_labeled_as_annotated_exon() 
    check2 = nextInwpCBG.get_bitscore() > firstInwpCBG.get_bitscore() 
    check3 = len(nextNodeList) > len(firstNodeList)
    check4 = float(firstInwpCBG.count_orfs_labeled_as_annotated_exon()) /\
             float(len(GENE_IDENTIFIER_SET)) <= 0.33
    if ntdistdict:
        check5 = sum(ntdistdict.values())/float(len(ntdistdict)) >\
             MIN_INTERGENIC_NT_LENGTH
    else:
        check5 = False
    if tcodedistdict:
        check6 = sum(tcodedistdict.values())/float(len(tcodedistdict)) <\
             TCODE_MAX_NONCODING
    else:
        check6 = False
    check7 = nextInwpCBG.count_orfs_labeled_as_first_exon() >=\
             firstInwpCBG.count_orfs_labeled_as_first_exon()
    check8 = firstInwpCBG.count_orfs_labeled_as_annotated_exon() == 0 
    check9 = nextInwpCBG.get_average_upstream_methionine_pssm_score() >\
             firstInwpCBG.get_average_upstream_methionine_pssm_score() 

    checklist = [check1,check2,check3,check4,check5,check6,check7,check8,check9]

    ############################################################################
    if verbose or True: print "NonFirst inwpCBG check:", checklist
    ############################################################################
   
    if checklist.count(False) <= 1:
        nonfirstPCG = PacbpCollectionGraph(crossdata={},blastmatrix=PCG._blastmatrix)
        # place all PacbPORFs in the nonfirstPCG
        for (pacbpkey,nodeQ,nodeS), pacbporf in firstInwpCBG.pacbps.iteritems():
            # add to noncodingnongenePCG
            nonfirstPCG.add_node(nodeQ)
            nonfirstPCG.add_node(nodeS)
            nonfirstPCG.add_edge(nodeQ,nodeS,wt=pacbporf.bitscore)
            nonfirstPCG.pacbps[(pacbpkey,nodeQ,nodeS)] = pacbporf
            # remove from main PCG
            _delete_pacbp(PCG,(pacbpkey,nodeQ,nodeS))
        # return nonfirstPCG
        return nonfirstPCG
    else:
        return False

# end of function detect_and_remove_single_nonfirst_inwpcbg


def detect_and_remove_gtgdiscrepancy(inwpcbgs,PCG,GENE_IDENTIFIER_SET,verbose=True):
    """ """

    # if empty list or empty PCG provided: return False
    if not inwpcbgs or not PCG or PCG.node_count() == 0: return False

    # get target organism identifier
    target = inwpcbgs[0]._get_target_organism()

    # Make *the* GTG of the strongest X informant species
    # X depends on the maximum number of gene informants (GENE_IDENTIFIER_SET);
    # unigene informants are not taken into account here.
    # X is defined here by:
    # -- at least 3 informants (for very small number of informants)
    # -- optimally half of the total numers of informants
    # -- at most 8 informants
    min_gtg_node_count = 3 + 1
    max_gtg_node_count = 8 + 1
    gtg_size = min([(len(GENE_IDENTIFIER_SET)-1)/2, max_gtg_node_count])
    gtg_size = max([min_gtg_node_count,gtg_size])

    btGTG = pcg2gtg_by_bitscore(PCG,target,identifier_list=GENE_IDENTIFIER_SET)
    ntGTG = pcg2gtg_by_identity(PCG,target,identifier_list=GENE_IDENTIFIER_SET)

    # TEMP solution because OrganismGraph != OrganismStarGraph
    # make bitscore ordered list of nodes
    bitscore_ordered_nodes = []
    for (tNode,iNode),wt in btGTG.weights.iteritems():
        if tNode==target: bitscore_ordered_nodes.append( ( wt, iNode ) )
    bitscore_ordered_nodes.sort() 
    #if verbose: print "btGTG::", bitscore_ordered_nodes

    while ntGTG.node_count() > gtg_size:
        # next line causes errors because OrganismGraph != OrganismStarGraph
        # this causes the target node in rare cases to be assigned as the weakest node
        # informant = btGTG.weakest_connected_node()
        (wt,informant) = bitscore_ordered_nodes.pop(0)
        btGTG.del_node(informant)
        ntGTG.del_node(informant)
        if verbose: print "btGGT.weakest_connected_node() ==", informant, btGTG.get_ordered_nodes()

    ############################################################################
    if verbose:
        print "ntGTG:", ntGTG.get_ordered_nodes(), 
        for node in ntGTG.get_ordered_nodes():
            if node == target: continue
            print "%1.2f" % ntGTG.weights[(target,node)],
        print ""
    ############################################################################

    # detect inwpCBGs which are probably the result of intron alignments
    gtgdiscrepancy_internal_inwpcbg_list = assign_internal_nongene_alignments(inwpcbgs,ntGTG)


    # detect inwpCBGs with strong discrepancy to this GTG
    gtgdiscrepancy_inwpcbg_list = assign_gtgdiscrepancy_inwpcbgs(inwpcbgs,ntGTG)

    # merge both lists
    if gtgdiscrepancy_internal_inwpcbg_list:
        if not gtgdiscrepancy_inwpcbg_list:
            gtgdiscrepancy_inwpcbg_list.extend(gtgdiscrepancy_internal_inwpcbg_list)
        else:
            for inwpcbg in gtgdiscrepancy_internal_inwpcbg_list:
                check_str = str(inwpcbg)
                if check_str not in [ str(gtgdiscrCBG) for gtgdiscrCBG in gtgdiscrepancy_inwpcbg_list ]:
                    gtgdiscrepancy_inwpcbg_list.append( inwpcbg )

    if not gtgdiscrepancy_inwpcbg_list:
        return False

    # get list of inwpCBGs that have NO discrepancy
    correct_inwpcbg_list = []
    check_str_list = []
    for discrinwpCBG in gtgdiscrepancy_inwpcbg_list:
        check_str_list.append( str(discrinwpCBG) )
    for inwpcbg in inwpcbgs:
        if str(inwpcbg) not in check_str_list:
            correct_inwpcbg_list.append( inwpcbg )

    # get all pacbp keys belonging to gtgdiscrepancy inwpcbgs ONLY
    gtgdiscrepancy_pacbpkeys = []
    for discrinwpCBG in gtgdiscrepancy_inwpcbg_list:
        for pacbpkey in discrinwpCBG.pacbps.keys():
            # check if this pacbpkey is occuring in a non-removed inwpCBG
            is_occurring_in_correct_inwpcbg = False
            for inwp in correct_inwpcbg_list:
                if pacbpkey in inwp.pacbps.keys():
                    is_occurring_in_correct_inwpcbg = True
                    break
            # if is_occurring_in_correct_inwpcbg, continue and do not delete
            if is_occurring_in_correct_inwpcbg:
                continue
            # store to gtgdiscrepancy_pacbpkeys when not stored already
            if pacbpkey not in gtgdiscrepancy_pacbpkeys:
                gtgdiscrepancy_pacbpkeys.append(pacbpkey)


    # place all gtgdiscrepancy_pacbpkeys and PacbPORFs in the gtgdiscrepancyPCG
    # and, at the same time, remove from the main PCG
    gtgdiscrepancyPCG = PacbpCollectionGraph(crossdata={},blastmatrix=PCG._blastmatrix)
    for key in gtgdiscrepancy_pacbpkeys:
        if key not in PCG.pacbps.keys():
            # !?!? TODO why not present in the PCG !?!?!
            # anyway, continue here to avoid KeyError
            # This PacbPORF was to be deleted rigth here,
            # so it is not an extreme disaster. But... scary ;-)
            continue
        (pacbpkey,nodeQ,nodeS) = key
        pacbporf = PCG.pacbps[key]
        # add to gtgdiscrepancyPCG
        gtgdiscrepancyPCG.add_node(nodeQ)
        gtgdiscrepancyPCG.add_node(nodeS)
        gtgdiscrepancyPCG.add_edge(nodeQ,nodeS,wt=pacbporf.bitscore)
        gtgdiscrepancyPCG.pacbps[(pacbpkey,nodeQ,nodeS)] = pacbporf

        # remove from main PCG
        _delete_pacbp(PCG,key)


    # return gtgdiscrepancyPCG
    return gtgdiscrepancyPCG

# end of function detect_and_remove_gtgdiscrepancy


def detect_and_remove_utrornonegene_inwpcbgs(inwpcbgs,PCG,verbose=True):
    """ """

    # if empty list or empty PCG provided: return False
    if not inwpcbgs or not PCG or PCG.node_count() == 0: return False

    # MAKE SHURE ALL Orfs HAVE PREDICTED TSS SITES!!
    for inwpCBG in inwpcbgs: inwpCBG.scan_orfs_for_pssm_tss(min_pssm_score=TSS_MIN_PSSM_SCORE)

    # get target organism identifier
    target = inwpcbgs[0]._get_target_organism()

    # detect inwpCBGs which are most likely 5' and 3' non coding or non gene
    ncng_5p_list = assign_utrornongene5p_inwpcbgs(inwpcbgs)
    ncng_3p_list = assign_utrornongene3p_inwpcbgs(inwpcbgs)
    ncng_list = ncng_5p_list
    ncng_list.extend(ncng_3p_list)

    # return False in no inwpcbgs are assigned
    if not ncng_list: return False

    # get list of inwpCBGs that are NON ncng
    correct_inwpcbg_list = []
    check_str_list = []
    for discrinwpCBG in ncng_list:
        check_str_list.append( str(discrinwpCBG) )
    for inwpcbg in inwpcbgs:
        if str(inwpcbg) not in check_str_list:
            correct_inwpcbg_list.append( inwpcbg )

    # get all pacbp keys belonging to noncoding / nongene inwpcbgs ONLY
    ncng_pacbpkeys = []
    for ncnginwpCBG in ncng_list:
        for pacbpkey in ncnginwpCBG.pacbps.keys():
            # check if this pacbpkey is occuring in a non-removed inwpCBG
            is_occurring_in_correct_inwpcbg = False
            for inwp in correct_inwpcbg_list:
                if pacbpkey in inwp.pacbps.keys():
                    is_occurring_in_correct_inwpcbg = True
                    break
            # if is_occurring_in_correct_inwpcbg, continue and do not delete
            if is_occurring_in_correct_inwpcbg:
                continue
            # store to gtgdiscrepancy_pacbpkeys when not stored already
            if pacbpkey not in ncng_pacbpkeys:
                ncng_pacbpkeys.append(pacbpkey)


    # place all ncng_pacbpkeys and PacbPORFs in the noncodingnongenePCG
    # and, at the same time, remove from the main PCG
    noncodingnongenePCG = PacbpCollectionGraph(crossdata={},blastmatrix=PCG._blastmatrix)
    for key in ncng_pacbpkeys:
        (pacbpkey,nodeQ,nodeS) = key
        pacbporf = PCG.pacbps[key]
        # add to noncodingnongenePCG
        noncodingnongenePCG.add_node(nodeQ)
        noncodingnongenePCG.add_node(nodeS)
        noncodingnongenePCG.add_edge(nodeQ,nodeS,wt=pacbporf.bitscore)
        noncodingnongenePCG.pacbps[(pacbpkey,nodeQ,nodeS)] = pacbporf
        # remove from main PCG
        _delete_pacbp(PCG,key)

    # return noncodingnongenePCG
    return noncodingnongenePCG

# end of function detect_and_remove_utrornonegene_inwpcbgs



def assign_internal_nongene_alignments(inwpcbgs,GTG,exclude_annotated=False,verbose=True):
    """
    TODO TODO: this function must be moved to another location.
    TODO TODO: better place in inwpCBGs/blocks filtering
    """ 
    # return empty list when no inwpcbgs applied
    if not inwpcbgs: return []

    # get target organism identifier
    target = inwpcbgs[0]._get_target_organism()

    # return list with inwpcbgs
    gtgdiscrepancy_inwpcbg_list = []

    # get most likely first & final inwpCBG pointer in inwpcbgs list
    posFirst,posFinal = get_first_and_final_inwpcbg_pos(inwpcbgs)

    # check if posFirst,posFinal+1 isa non-empty range
    if not range(posFirst,posFinal+1): return [] 

    # get info on the *best* covered inwpCBG
    best_nt_identity = max([ inwpcbgs[pos].get_nt_identity() for pos in range(posFirst,posFinal+1) ])
    best_bitscore    = max([ inwpcbgs[pos].get_bitscore() for pos in range(posFirst,posFinal+1) ])
    best_annot_cnt   = max([ inwpcbgs[pos].count_orfs_labeled_as_annotated_exon() for pos in range(posFirst,posFinal+1) ])
    best_bits_per_aa = max([ float( inwpcbgs[pos].get_bitscore() ) / float( sum([pf.get_unextended_length() for pf in inwpcbgs[pos].pacbps.values() ]) ) for pos in range(posFirst,posFinal+1) ])

    for pos in range(posFirst,posFinal+1):
        # get this inwpCBG and 
        thisInwpCBG = inwpcbgs[pos]
        if pos > 0: prevInwpCBG = inwpcbgs[pos-1]
        else:       prevInwpCBG = None
        if pos < len(inwpcbgs)-1: nextInwpCBG = inwpcbgs[pos+1]
        else:                     nextInwpCBG = None

        if prevInwpCBG and Set(prevInwpCBG.get_nodes()).intersection(thisInwpCBG.get_nodes()):
            continue
        if nextInwpCBG and Set(nextInwpCBG.get_nodes()).intersection(thisInwpCBG.get_nodes()):
            continue

        tot_length  = sum([pf.get_unextended_length() for pf in thisInwpCBG.pacbps.values() ])
        bits        = thisInwpCBG.get_bitscore()
        bits_per_aa = float(bits)/float(tot_length)


        if bits_per_aa/best_bits_per_aa < 0.55 and\
        float(thisInwpCBG.count_orfs_labeled_as_annotated_exon()) /\
        float(best_annot_cnt) <= 0.50:
            minsr = thisInwpCBG.minimal_spanning_range(organism=target)
            print " __XX__", pos, thisInwpCBG, thisInwpCBG.get_nt_identity(), bits, float(bits)/float(tot_length), thisInwpCBG.count_orfs_labeled_as_annotated_exon()
            if prevInwpCBG and len(prevInwpCBG.maximal_spanning_range(organism=target).intersection(minsr)) == len(minsr): 
                ################################################################
                if verbose:
                    print "PREV::", len(prevInwpCBG.maximal_spanning_range(organism=target).intersection(minsr)), len(minsr)
                ################################################################
                gtgdiscrepancy_inwpcbg_list.append( thisInwpCBG )
                continue
            if nextInwpCBG and len(nextInwpCBG.maximal_spanning_range(organism=target).intersection(minsr)) == len(minsr):
                ################################################################
                if verbose:
                    print "NEXT::", len(nextInwpCBG.maximal_spanning_range(organism=target).intersection(minsr)), len(minsr)
                ################################################################
                gtgdiscrepancy_inwpcbg_list.append( thisInwpCBG )
                continue

    # return list with conflicts
    return gtgdiscrepancy_inwpcbg_list

# end of function assign_internal_nongene_alignments


def assign_gtgdiscrepancy_inwpcbgs(inwpcbgs,GTG,exclude_annotated=True,verbose=True):
    """ """
    # return empty list when no inwpcbgs applied
    if not inwpcbgs: return []

    # get target organism identifier
    target = inwpcbgs[0]._get_target_organism()

    # return list with inwpcbgs
    gtgdiscrepancy_inwpcbg_list = []

    if exclude_annotated:
        # get most likely first & final inwpCBG pointer in inwpcbgs list
        posFirst,posFinal = get_first_and_final_inwpcbg_pos(inwpcbgs)
        range_5p_test = range(0,posFirst)
        range_3p_test = range(posFinal+1,len(inwpcbgs))
        protected_target_orfid_list = []
        for inwpCBG in inwpcbgs[posFirst:posFinal+1]:
            if inwpCBG.count_orfs_labeled_as_annotated_exon() > 0:
                protected_target_orfid_list.append( inwpCBG.get_orfs_of_graph(organism=target)[0].id )
    else:
        range_5p_test = []
        range_3p_test = []
        protected_target_orfid_list = []

    ############################################################################
    if verbose and exclude_annotated:
        print "NOT-excluded:", range_5p_test, range_3p_test
    ############################################################################

    # detect UTR or nongene / noncoding inwpCBGS
    for pos in range(0,len(inwpcbgs)):
        if exclude_annotated and pos in range_5p_test:
            pass
        elif exclude_annotated and pos in range_3p_test:
            pass
        elif exclude_annotated and inwpcbgs[pos].count_orfs_labeled_as_annotated_exon() == 0:
            # in the middle of the annotated geen structure, but not a single
            # Orf annotated as an exon. Asses for gtg difference too!
            pass
        elif exclude_annotated:
            continue
        else:
            pass


        # get this inwpCBG and 
        thisInwpCBG = inwpcbgs[pos]

        # ignore if the target's Orf is belonging to a `protected` Orf
        if protected_target_orfid_list and\
        thisInwpCBG.get_orfs_of_graph(organism=target)[0].id in\
        protected_target_orfid_list:
            continue

        # ignore inwpCBGs which are very likely (poor quality) SignalP alignments
        cntSP = float(thisInwpCBG.count_orfs_with_signalpeptides())
        if cntSP/(thisInwpCBG.count_genomic_informants()+1) > 0.66 and\
        thisInwpCBG.get_signalp_score() > 0.75:
            continue

        # create its GeneTreeGraph
        gtg = pcg2gtg_by_identity(thisInwpCBG,target)

        # step 1. Do the gtg/GTG difference check
        difference = _relative_gtg_difference(gtg,GTG,target)

        if difference < NONGENE_GTG_MAX_DIFFERENCE:
            # step 2. Do the CEXPANDER check
            if thisInwpCBG.node_count() <= 2:
                gtgdiscrepancy_inwpcbg_list.append(thisInwpCBG)
                ################################################################
                if verbose:
                    print pos, "thisInwpCBG", "gtg2GTGdiff:: %1.3f < %1.3f" % (
                        difference,NONGENE_GTG_MAX_DIFFERENCE),
                    print thisInwpCBG.get_organism_nodes(target)[0]
                ################################################################
            elif thisInwpCBG.get_cexpander_uniformly_aligned_count() == 0:
                gtgdiscrepancy_inwpcbg_list.append(thisInwpCBG)
                ################################################################
                if verbose:
                    print pos, "thisInwpCBG", "gtg2GTGdiff:: %1.3f < %1.3f" % (
                        difference,NONGENE_GTG_MAX_DIFFERENCE),
                    print thisInwpCBG.get_organism_nodes(target)[0]
                ################################################################
            else:
                # cexpander check is succesfull, GTGdifference claims
                # the aligment is bogus. Do a more elaborate check on
                # some other variables of thisInwpCBG

                # calculate the difference between minsr & maxsr lengths
                node      = thisInwpCBG.get_organism_nodes(target)[0]
                minsr     = thisInwpCBG.minimal_spanning_range_sizes()[node]
                maxsr     = thisInwpCBG.maximal_spanning_range_sizes()[node]
                msr_ratio = float(minsr)/float(maxsr)

                # calculate the ratio between average weights of gtg and GTG
                average_wt_gtg = _pairwise_gtg_average_weight(gtg,target)
                average_wt_GTG = _pairwise_gtg_average_weight(GTG,target)
                gtg_ratio = average_wt_gtg / average_wt_GTG

                if msr_ratio < NONGENE_GTG_MAX_MSR_RATIO and\
                gtg_ratio < NONGENE_GTG_MAX_GTG_RATIO:
                    gtgdiscrepancy_inwpcbg_list.append(thisInwpCBG)
                    ################################################################
                    if verbose:
                        print pos, "thisInwpCBG", "gtg2GTGdiff:: %1.3f < %1.3f" % (
                            difference,NONGENE_GTG_MAX_DIFFERENCE),
                        print thisInwpCBG.get_organism_nodes(target)[0]
                    ################################################################
                else:
                    pass
        else:
            pass


    # return the gtgdiscrepancy_inwpcbg_list
    return gtgdiscrepancy_inwpcbg_list

# end of function assign_gtgdiscrepancy_inwpcbgs



def _relative_gtg_difference(gtg,GTG,target):
    """ Obtain relative graph difference between a gtg and *the* GTG """
    differences = []
    total_GTG_weights = []
    for informant in GTG.organism_set().intersection(gtg.organism_set()):
        if informant == target: continue
        wtGTG = GTG.weights[(target,informant)]
        wtgtg = gtg.weights[(target,informant)]
        differences.append( wtgtg - wtGTG )
        total_GTG_weights.append( wtGTG )
    for informant in GTG.organism_set().difference(gtg.organism_set()):
        wtGTG = GTG.weights[(target,informant)]
        differences.append( - wtGTG )
        total_GTG_weights.append( wtGTG )
    # return (float) sum of relative difference
    return sum(differences)/sum(total_GTG_weights)

# end of function _relative_gtg_difference


def _pairwise_gtg_total_weight(gtg,target):
    """ """
    total_weights = []
    for informant in gtg.organism_set():
        if informant == target: continue
        total_weights.append( gtg.weights[(target,informant)] )
    return sum(total_weights)

# end of function _pairwise_gtg_total_weight


def _pairwise_gtg_average_weight(gtg,target):
    """ """
    total_weights = []
    for informant in gtg.organism_set():
        if informant == target: continue
        total_weights.append( gtg.weights[(target,informant)] )
    return sum(total_weights) / (gtg.node_count()-1)

# end of function _pairwise_gtg_average_weight


# end of function _pairwise_gtg_total_weight

def assign_utrornongene5p_inwpcbgs(inwpcbgs,verbose=True):
    """ """
    # return empty list when no inwpcbgs applied
    if not inwpcbgs: return []

    # get most likely first & final inwpCBG pointer in inwpcbgs list
    posFirst,posFinal = get_first_and_final_inwpcbg_pos(inwpcbgs)

    # return variable list
    noncoding_inwpcbg_list = []

    # get data of most likely first inwpCBG
    max_cntAnnot        = max([inwp.count_orfs_labeled_as_annotated_exon() for inwp in inwpcbgs ])
    firstInwpCBG        = inwpcbgs[posFirst]
    first_cnt_is_first  = firstInwpCBG.count_orfs_labeled_as_first_exon()
    first_identityscore = firstInwpCBG.get_identityscore()
    first_upstrTSScnt   = [ pf.has_upstream_tss() for pf in firstInwpCBG.pacbps.values() ].count(True)
    if (max_cntAnnot-1) == 0:
        # avoid ZeroDivisionError
        first_upstrTSSratio = 0.0
    else:
        first_upstrTSSratio = float(first_upstrTSScnt) / (max_cntAnnot-1)

    # range of inwpCBGs which are checked for deletion
    range_5p_test = range(0,posFirst)

    # detect UTR or nongene / noncoding inwpCBGS
    for pos in range(0,len(inwpcbgs)):
        if pos not in range_5p_test: continue

        # get this inwpCBG and get statistics
        inwpCBG = inwpcbgs[pos]
        cntFirst = inwpCBG.count_orfs_labeled_as_first_exon()
        # calculated differently as for *the* firstCBG
        cntAnnot = float(inwpCBG.organism_set_size())

        # break when to putatively first is reached
        if cntFirst == first_cnt_is_first: break

        # do is_coding() test
        iscoding = inwpCBG.is_coding()
        # calculate cnt/ratio of upstrTSS sites
        this_upstrTSScnt   = [ pf.has_upstream_tss() for pf in inwpCBG.pacbps.values() ].count(True)
        if (max_cntAnnot-1) == 0:
            this_upstrTSSratio = 0.0
        else:
            this_upstrTSSratio = float(this_upstrTSScnt) / (max_cntAnnot-1.0)

        ########################################################################
        if verbose:
            print pos, range_5p_test, "5'UTR analyses:", inwpCBG, iscoding,
            print "coverage:",cntAnnot/max_cntAnnot,
            print "upstrTSS cnt - ratio: %s - %1.2f" % (
                    this_upstrTSScnt, this_upstrTSSratio)
        ########################################################################

        if not iscoding:
            # inwpCBGs most likely not coding alignments -> remove
            noncoding_inwpcbg_list.append(inwpCBG)
            continue

        # check relative position towards the current firstInwpCBG
        # position is measured in actual nt distance and tcode 'distance':
        # the lowest scoring TCODE window in between these inwpCBGs
        tcodedistdict = inwpCBG.tcode_spacing_between_codingblocks([firstInwpCBG])

        if len(tcodedistdict) >= 3:
            _tmp = [ (v,k) for k,v in tcodedistdict.iteritems() ]
            _tmp.sort()
            del( tcodedistdict[_tmp[0][1]] )
            del( tcodedistdict[_tmp[-1][1]] )


        ########################################################################
        if verbose and len(tcodedistdict) >= 1:
            print pos, sum(tcodedistdict.values())/len(tcodedistdict),
            print TCODE_MAX_NONCODING, firstInwpCBG
        ########################################################################

        # continue when coverage is to high
        if cntAnnot/max_cntAnnot >= 0.40: continue

        if len(tcodedistdict)>0 and\
        sum(tcodedistdict.values())/len(tcodedistdict) <= TCODE_MAX_NONCODING:
            noncoding_inwpcbg_list.append(inwpCBG)
            continue


        if this_upstrTSScnt == 0:
            # no upstream TSS sites at all -> remove!
            noncoding_inwpcbg_list.append(inwpCBG)
            continue
    
        # do a furter check for unlikely first inwpCBG blocks
        if inwpCBG.get_identityscore() / first_identityscore <=\
        NONCODINGNONGENE_5p_INWPCBG_MAX_IDENTITYRATIO and\
        inwpCBG.organism_set_size() < first_cnt_is_first and\
        first_upstrTSSratio==0.0:
            noncoding_inwpcbg_list.append(inwpCBG)
            continue


        # do a furter check for unlikely first inwpCBG blocks
        if inwpCBG.get_identityscore() / first_identityscore <=\
        NONCODINGNONGENE_5p_INWPCBG_MAX_IDENTITYRATIO and\
        inwpCBG.organism_set_size() < first_cnt_is_first and\
        (first_upstrTSSratio!=0.0 and (this_upstrTSSratio / first_upstrTSSratio) < 0.6):
            noncoding_inwpcbg_list.append(inwpCBG)
            continue

        # do a final check for unlikely first inwpCBG blocks
        # all parameters must be (slightly) poorer
        if inwpCBG.get_identityscore() < first_identityscore and\
        inwpCBG.organism_set_size() < first_cnt_is_first and\
        inwpCBG.get_average_upstream_methionine_pssm_score() <\
        firstInwpCBG.get_average_upstream_methionine_pssm_score():
            noncoding_inwpcbg_list.append(inwpCBG)
            continue

    # return the noncoding_inwpcbg_list
    return noncoding_inwpcbg_list

# end of function assign_utrornongene5p_inwpcbgs


def assign_utrornongene3p_inwpcbgs(inwpcbgs,verbose=True):
    """ """
    # return empty list when no inwpcbgs applied
    if not inwpcbgs: return []

    # get most likely first & final inwpCBG pointer in inwpcbgs list
    posFirst,posFinal = get_first_and_final_inwpcbg_pos(inwpcbgs)

    # return variable list
    noncoding_inwpcbg_list = []

    # get data of most likely first inwpCBG
    max_cntAnnot        = max([inwp.count_orfs_labeled_as_annotated_exon() for inwp in inwpcbgs ])
    finalInwpCBG        = inwpcbgs[posFinal]
    final_cnt_is_final  = finalInwpCBG.count_orfs_labeled_as_final_exon()
    final_identityscore = finalInwpCBG.get_identityscore()
    final_prjtls_aadif  = finalInwpCBG.get_projected_tailing_stop_aa_difference()
    final_prjtls_nonad  = finalInwpCBG.get_projected_tailing_stop_nonaligned_aa_difference()

    # range of inwpCBGs which are checked for deletion
    range_3p_test = range(posFinal+1,len(inwpcbgs))

    # detect UTR or nongene / noncoding inwpCBGS
    for pos in range(0,len(inwpcbgs)):
        if pos not in range_3p_test: continue

        # get this inwpCBG and get statistics
        inwpCBG = inwpcbgs[pos]
        cntFinal = inwpCBG.count_orfs_labeled_as_first_exon()
        # calculated differntly as for *the* firstCBG
        cntAnnot = float(inwpCBG.organism_set_size())

        # break when to putatively first is reached
        if cntFinal == final_cnt_is_final: break

        # remove poorly covered inwpCBGs with low identityscore and not
        # having a likely stop codon
        if cntAnnot/max_cntAnnot < 0.80 and\
        inwpCBG.get_identityscore() / final_identityscore <=\
        NONCODINGNONGENE_3p_INWPCBG_MAX_IDENTITYRATIO and\
        inwpCBG.organism_set_size() < final_cnt_is_final and\
        final_prjtls_aadif < inwpCBG.get_projected_tailing_stop_aa_difference() and\
        final_prjtls_nonad < inwpCBG.get_projected_tailing_stop_nonaligned_aa_difference():
            # not contribution to the gene structure at all....
            noncoding_inwpcbg_list.append(inwpCBG)
            continue

        # check relative position towards the current finalInwpCBG
        # position is measured in actual nt distance and tcode 'distance':
        # the lowest scoring TCODE window in between these inwpCBGs
        ntdistdict    = finalInwpCBG.nt_spacing_between_codingblocks([inwpCBG])
        tcodedistdict = finalInwpCBG.tcode_spacing_between_codingblocks([inwpCBG])

        # remove highest & lowest distance and then do stats on remaining dists
        if len(ntdistdict) >= 3:
            _tmp = [ (v,k) for k,v in ntdistdict.iteritems() ]
            _tmp.sort()
            del( ntdistdict[_tmp[0][1]] )
            del( ntdistdict[_tmp[-1][1]] )
        if len(tcodedistdict) >= 3:
            _tmp = [ (v,k) for k,v in tcodedistdict.iteritems() ]
            _tmp.sort()
            del( tcodedistdict[_tmp[0][1]] )
            del( tcodedistdict[_tmp[-1][1]] )

        # do 3 checks.
        # 1) are average and maximum intergenic distances are bridged?
        # 2) is stop codon projection a deterioration?
        # 4) does tcodedistance suggest bridging of a non-coding stretch?
        check_1 = sum(ntdistdict.values())/len(ntdistdict) >=\
                  AVERAGE_INTERGENIC_MIN_NT_LENGTH and\
                  max(ntdistdict.values()) >= MAX_INTERGENIC_MIN_NT_LENGTH
        check_2 = final_prjtls_aadif <\
                  inwpCBG.get_projected_tailing_stop_aa_difference() and\
                  final_prjtls_nonad <\
                  inwpCBG.get_projected_tailing_stop_nonaligned_aa_difference()
        check_3 = len(tcodedistdict)>0 and\
                  sum(tcodedistdict.values())/len(tcodedistdict) <=\
                  TCODE_MAX_NONCODING

        if [ check_1, check_2, check_3 ].count(True) >= 2:
            # not contribution to the gene structure at all....
            noncoding_inwpcbg_list.append(inwpCBG)
            continue


        # do is_coding() test
        iscoding = inwpCBG.is_coding()

        ########################################################################
        if verbose:
            print pos, "3'UTR analyses:", inwpCBG, iscoding,
            print cntAnnot/max_cntAnnot
        ########################################################################

        if not iscoding:
            # probably non-coding inwp CBG alignment block
            noncoding_inwpcbg_list.append(inwpCBG)
            continue


    # return the noncoding_inwpcbg_list
    return noncoding_inwpcbg_list


# end of function assign_utrornongene3p_inwpcbgs

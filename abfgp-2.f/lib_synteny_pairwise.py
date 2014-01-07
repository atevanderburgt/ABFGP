"""
"""

# Module metadata
__authors__ = "Ate van der Burgt"
__license__ = "MIT"

# ABGP Imports
from graphAbgp import (
    GeneTreeGraph,
    PacbpCollectionGraph,
    InwardsPointingCodingBlockGraph,
    )
from graphAbgp.graph_pacbpcollection import _delete_pacbp



from lib_pcg2blocks import (
    PCG2inwpCBGS,
    print_inwpcbgstructure,
    )

# Python Imports


# Global Variable Imports
from settings.genestructure import (
    MIN_INTERGENIC_NT_LENGTH,
    MAX_INTERGENIC_MIN_NT_LENGTH,
    AVERAGE_INTERGENIC_MIN_NT_LENGTH,
    )


    
def inwpcbgs2knownorftransitiondata(inwpcbglist):
    """
    """
    _data = []
    for inwpCBG in inwpcbglist:
        _data.append( set(inwpCBG.get_orfs_labeled_as_annotated_exon().keys()) )

    cnts = []
    for offset in range(0,len(inwpcbglist)):
        cnt5p = set()
        cnt3p = set()
        for pos in range(0,offset+1): cnt5p.update(_data[pos])
        for pos in range(offset+1,len(inwpcbglist)): cnt3p.update(_data[pos])
        cnts.append( ( len(cnt5p), len(cnt3p) ) )

    return cnts

# end of function inwpcbgs2knownorftransitiondata
        
    
    
def detect_and_remove_synteny(inwpcbgs,PCG,GENE_IDENTIFIER_SET,verbose=True):
    """ """
    MIN_OBSERVED_VS_EXPECTED_RATIO = 0.20

    observed_organism_subcombis = []
    syntenic_subinwpcbgs = []

    # detect syntenic genes in MAIN inwpCBGs,
    # without taking strongest informants by GTG analyses
    syntenic_inwpcbgs = assign_syntenic_inwpcbgs(inwpcbgs)

    for syntinwpcbg in syntenic_inwpcbgs:
        syntenic_subinwpcbgs.append(syntinwpcbg)

    for inwpCBG in inwpcbgs:
        # omit inwpCBGs with annotated exons/orfs
        if inwpCBG.count_orfs_labeled_as_annotated_exon() >= 2: continue
        target = inwpCBG._get_target_organism()

        # make a (artificially fully connected) GeneTreeGraph
        gtg = GeneTreeGraph()
        gtg.add_node(target)
        for (pacbpkey,nodeQ,nodeS),pacbporf in inwpCBG.pacbps.iteritems():
            orgS = inwpCBG.organism_by_node(nodeS)
            if orgS not in GENE_IDENTIFIER_SET: continue
            gtg.add_node(orgS)
        for (pacbpkey,nodeQ,nodeS),pacbporf in inwpCBG.pacbps.iteritems():
            orgQ = inwpCBG.organism_by_node(nodeQ)
            orgS = inwpCBG.organism_by_node(nodeS)
            if orgS not in GENE_IDENTIFIER_SET: continue
            gtg.add_edge( orgQ, orgS, wt = pacbporf.bitscore )
    
            # make artificially missed edges between the informants
            for org in inwpCBG.organism_set():
                if org not in [orgQ,orgS] and org in GENE_IDENTIFIER_SET:
                    if gtg.has_edge( orgS, org ) and\
                    gtg.weights[(orgS, org)] > pacbporf.bitscore:
                        gtg.set_edge_weight(orgS,org,wt = pacbporf.bitscore)
                    else:
                        gtg.add_edge( orgS, org, wt = pacbporf.bitscore )
    
        # omit (nearly) empty genetreegraphs
        if gtg.node_count() <= 1: continue

        # remove (much) weaker connected nodes as expected from the gtg
        while gtg.get_nodes() and MIN_OBSERVED_VS_EXPECTED_RATIO >\
        min( [ gtg.get_node_weighted_connectivity_observed_vs_expected(node) for node in gtg.get_nodes() ]):
            node = gtg.weakest_connected_node()
            gtg.del_node(node)
    
        # check if already tested before; present in observed_organism_subcombis
        if gtg.get_ordered_nodes() in observed_organism_subcombis: continue
    
        # store to already tested organism subcombinations
        observed_organism_subcombis.append( gtg.get_ordered_nodes() )
    
        # create a subPCG of these organisms
        subPCG = PacbpCollectionGraph(crossdata={},
                    blastmatrix=PCG._blastmatrix)
        for (pacbpkey,nodeQ,nodeS), pacbporf in PCG.pacbps.iteritems():
            (orgQ,orfQid),(orgS,orfSid) = nodeQ,nodeS
            if orgQ not in gtg.get_nodes(): continue
            if orgS not in gtg.get_nodes(): continue
            subPCG.add_node(nodeQ)
            subPCG.add_node(nodeS)
            subPCG.add_edge(nodeQ,nodeS,wt=pacbporf.bitscore)
            subPCG.pacbps[(pacbpkey,nodeQ,nodeS)] = pacbporf
    
        # make inwpCBGs of this subPCG
        subinwpcbgs = PCG2inwpCBGS(subPCG)

        # check if there are subinwpcbgs
        if not subinwpcbgs: continue

        ########################################################################
        #if verbose:
        #    print "subPCG organism set:", gtg.get_ordered_nodes()
        #    print_inwpcbgstructure(subinwpcbgs,gtg.get_ordered_nodes())
        ########################################################################
    
        # create a subInwardsPointingCodingBlockGraph of these organisms
        #subinwpCBG = InwardsPointingCodingBlockGraph()
        #for (pacbpkey,nodeQ,nodeS), pacbporf in inwpCBG.pacbps.iteritems():
        #    (orgQ,orfQid),(orgS,orfSid) = nodeQ,nodeS
        #    if orgQ not in gtg.get_nodes(): continue
        #    if orgS not in gtg.get_nodes(): continue
        #    subinwpCBG.add_node(nodeQ)
        #    subinwpCBG.add_node(nodeS)
        #    subinwpCBG.add_edge(nodeQ,nodeS,wt=pacbporf.bitscore)
        #    subinwpCBG.pacbps[(pacbpkey,nodeQ,nodeS)] = pacbporf

        # detect syntenic genes in this subinwpcbgs
        syntenic_inwpcbgs = assign_syntenic_inwpcbgs(subinwpcbgs)

        for syntinwpcbg in syntenic_inwpcbgs:
            syntenic_subinwpcbgs.append(syntinwpcbg)
            ####################################################################
            if verbose:
                print "SYNTENIC!!", syntinwpcbg, syntinwpcbg.get_ordered_nodes()
                for subCBG in subinwpcbgs:
                    print "syntenic in:", subCBG, subCBG.get_ordered_nodes()
            ####################################################################

    if not syntenic_subinwpcbgs:
        return False

    # cleanup all inwpCBGs from the syntenic subInwpCBGs
    syntenic_pacbpkeys = []
    for syntinwpcbg in syntenic_subinwpcbgs:
        node_set = syntinwpcbg.node_set()
        for inwpCBG in inwpcbgs:
            if not node_set.difference(inwpCBG.node_set()):
                for pacbpkey in inwpCBG.pacbps.keys():
                    if pacbpkey not in syntenic_pacbpkeys:
                        syntenic_pacbpkeys.append(pacbpkey)

    # place all syntenic_pacbpkeys and PacbPORFs in the syntenicPCG
    # and, at the same time, remove from the main PCG
    syntenicPCG = PacbpCollectionGraph(crossdata={},blastmatrix=PCG._blastmatrix)
    for key in syntenic_pacbpkeys:
        (pacbpkey,nodeQ,nodeS) = key
        pacbporf = PCG.pacbps[key]
        # add to syntenicPCG
        syntenicPCG.add_node(nodeQ)
        syntenicPCG.add_node(nodeS)
        syntenicPCG.add_edge(nodeQ,nodeS,wt=pacbporf.bitscore)
        syntenicPCG.pacbps[(pacbpkey,nodeQ,nodeS)] = pacbporf
        # remove from main PCG
        _delete_pacbp(PCG,key)

    # return syntenicPCG
    return syntenicPCG

# end of function detect_and_remove_synteny



def get_first_and_final_inwpcbg_pos(inwpcbgs):
    """ """
    # assign (most likely) first and final inwpCBG in the list of inwpcbgs
    # based on applied gene annotations (majority voting)
    max_is_first_list = []
    max_is_final_list = []
    for inwpCBG in inwpcbgs:
        cntFirst = inwpCBG.count_orfs_labeled_as_first_exon()
        cntFinal = inwpCBG.count_orfs_labeled_as_final_exon()
        cntAnnot = inwpCBG.count_orfs_labeled_as_annotated_exon()
        cntInfrm = inwpCBG.organism_set_size()
        pos = len(max_is_first_list)
        max_is_first_list.append((cntFirst,cntAnnot,cntInfrm,pos))
        pos = len(max_is_final_list)
        max_is_final_list.append((cntFinal,cntAnnot,cntInfrm,pos))

    max_is_first_list.sort()
    max_is_first_list.reverse()
    max_is_final_list.sort()
    max_is_final_list.reverse()

    # get position of most likely first & final inwpCBG
    first_inwpcbg_pos = max_is_first_list[0][-1]
    final_inwpcbg_pos = max_is_final_list[0][-1]

    # check to which Orf object the first position belongs; potentially
    # the next position makes more sense
    if first_inwpcbg_pos > 0 and first_inwpcbg_pos < len(inwpcbgs) -1:
        prev, this, next = inwpcbgs[first_inwpcbg_pos-1:first_inwpcbg_pos+2]
        if prev._get_target_node() == this._get_target_node() and\
        next._get_target_node() != this._get_target_node() and\
        prev.count_orfs_labeled_as_annotated_exon() <\
        this.count_orfs_labeled_as_annotated_exon() and\
        next.count_orfs_labeled_as_annotated_exon() >\
        this.count_orfs_labeled_as_annotated_exon() and\
        next.organism_set_size() >\
        this.organism_set_size():
            # increase first_inwpcbg_pos
            first_inwpcbg_pos+=1
            
    # check to which Orf object the final position belongs; potentially
    # the next position makes more sense
    while final_inwpcbg_pos < len(inwpcbgs) -1:
        this, next = inwpcbgs[final_inwpcbg_pos:final_inwpcbg_pos+2]
        if next._get_target_node() == this._get_target_node():
            # increase final_inwpcbg_pos
            final_inwpcbg_pos+=1
        else:
            break
            
    # return first and final integer list position
    return (first_inwpcbg_pos,final_inwpcbg_pos)

# end of function get_first_and_final_inwpcbg_pos


def assign_syntenic_inwpcbgs(inwpcbgs):
    """ """

    # get most likely first & final inwpCBG pointer in inwpcbgs list
    posFirst,posFinal = get_first_and_final_inwpcbg_pos(inwpcbgs)
    # return variable list
    syntenic_inwpcbg_list = []

    # get target organism from inwpCBG(s)
    target = inwpcbgs[0]._get_target_organism()

    range_5p_test = range(0,posFirst)
    range_3p_test = range(posFinal+1,len(inwpcbgs))

    # obtain knownorftransitiondata to detect synteny
    knownorftransitiondata = inwpcbgs2knownorftransitiondata(inwpcbgs)

    print "assign_syntenic_inwpcbgs::", range_5p_test, "[%s]" % len(inwpcbgs), range_3p_test
    
    # detect exons of 5' syntenic genes
    for pos in range(0,len(inwpcbgs)-1):
        if pos in range_5p_test:
            pass
        else:
           continue

        thisInwpCBG = inwpcbgs[pos]
        nextInwpCBG = inwpcbgs[pos+1]
        
        # get KnownOrfTransition ratios
        thisKOTD5p,thisKOTD3p = knownorftransitiondata[pos]
        nextKOTD5p,nextKOTD3p = knownorftransitiondata[pos+1]
        thisKOTDratio = float(thisKOTD3p) / float(max([thisKOTD5p,1]))
        nextKOTDratio = float(nextKOTD3p) / float(max([nextKOTD5p,1]))
        try:
            ratioKOTD = thisKOTDratio / nextKOTDratio
        except ZeroDivisionError:
            ratioKOTD = 0.0
            
        if not thisInwpCBG.node_set().difference(nextInwpCBG.node_set()):
            # identical node sets shared -> this cannot be assigned to synteny
            continue

        if thisInwpCBG.node_set().intersection(nextInwpCBG.node_set()):
            # any nodes shared -> no syntheny
            continue

        # get distance info
        ntdistdict = thisInwpCBG.nt_spacing_between_codingblocks([nextInwpCBG])
    
        if min(ntdistdict.values()) >= MIN_INTERGENIC_NT_LENGTH and\
        max(ntdistdict.values()) >= MAX_INTERGENIC_MIN_NT_LENGTH:
            syntenic_inwpcbg_list.append( thisInwpCBG )
        elif min(ntdistdict.values()) >= MIN_INTERGENIC_NT_LENGTH and\
        sum(ntdistdict.values())/len(ntdistdict) >= AVERAGE_INTERGENIC_MIN_NT_LENGTH:
            syntenic_inwpcbg_list.append( thisInwpCBG )
        # thisKOTDratio >= 4.0 means 75% of the Orfs explained -> <25% unexplained
        # ratioKOTD >= 4.0 means much relatively 4x more Orfs explained if split here
        # LARGEST nt distance larger than MIN_INTERGENIC_NT_LENGTH
        elif thisKOTDratio >= 4.0 and ratioKOTD >= 4.0 and\
        max(ntdistdict.values()) >= MIN_INTERGENIC_NT_LENGTH:
            syntenic_inwpcbg_list.append( thisInwpCBG )
        else:
            pass

    # detect exons of 3' syntenic genes
    for pos in range(1,len(inwpcbgs)):
        if pos in range_3p_test:
            pass
        else:
            continue

        thisInwpCBG = inwpcbgs[pos]
        prevInwpCBG = inwpcbgs[pos-1]
        
        if not thisInwpCBG.node_set().difference(prevInwpCBG.node_set()):
            # identical node sets shared -> this cannot be assigned to synteny
            continue

        if thisInwpCBG.node_set().intersection(prevInwpCBG.node_set()):
            # any nodes shared -> no syntheny
            continue
            
        # get distance info
        ntdistdict = prevInwpCBG.nt_spacing_between_codingblocks([thisInwpCBG])
    
        if min(ntdistdict.values()) >= MIN_INTERGENIC_NT_LENGTH and\
        max(ntdistdict.values()) >= MAX_INTERGENIC_MIN_NT_LENGTH:
            syntenic_inwpcbg_list.append( thisInwpCBG )
        elif min(ntdistdict.values()) >= MIN_INTERGENIC_NT_LENGTH and\
        sum(ntdistdict.values())/len(ntdistdict) >= AVERAGE_INTERGENIC_MIN_NT_LENGTH:
            syntenic_inwpcbg_list.append( thisInwpCBG )
        else:
            pass

    # return the syntenic_inwpcbg_list
    return syntenic_inwpcbg_list

# end of function assign_syntenic_inwpcbgs

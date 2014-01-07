################################################################################
### Graph conversion functions                                              ####
################################################################################

# graphAbgp Imports
from graph_genetree import GeneTreeGraph
from graph_alignedsites import *
import ordering
from exceptions import *

# Import Pacb class
import pacb

# Import Python functions
from copy import deepcopy

def SpliceSiteCollectionGraph2AlignedSpliceSiteGraph(gra,max_node_count=None):
    """
    """
    if gra.node_count() == 1:
        # a Collection of 1 node can never become an Aligned graph!
        return gra
    elif gra.phase() in [0,1,2] and gra.connectivitysaturation() == 1.0:
        # Aligned Donor, Acceptor or Splice Site Collection Graph
        if gra.__class__.__name__ == 'DonorSiteCollectionGraph':
            newgra = AlignedDonorSiteGraph(
                    max_node_count=max_node_count,
                    aligned_site_aa_offset=gra.ALIGNED_SITE_AA_OFFSET,
                    min_pssm_score=gra.MIN_PSSM_SCORE 
                    )
        elif gra.__class__.__name__ == 'AcceptorSiteCollectionGraph':
            newgra = AlignedAcceptorSiteGraph(
                    max_node_count=max_node_count,
                    aligned_site_aa_offset=gra.ALIGNED_SITE_AA_OFFSET,
                    min_pssm_score=gra.MIN_PSSM_SCORE 
                    )
        elif gra.__class__.__name__ == 'SpliceSiteCollectionGraph':
            newgra = AlignedSpliceSiteGraph(
                    max_node_count=max_node_count,
                    aligned_site_aa_offset=gra.ALIGNED_SITE_AA_OFFSET,
                    min_pssm_score=gra.MIN_PSSM_SCORE 
                    )
        else:
            message = "wrong input object encountered: %s" % gra.__class__.__name__
            raise InproperlyAppliedArgument, message

    elif gra.connectivitysaturation() == 1.0 and gra.phase() == None:
        # Aligned CodingBlockStart or CodingBlockEnd CollectionGraph
        if gra.__class__.__name__ == 'DonorSiteCollectionGraph':
            newgra = AlignedCbg3pBoundarySiteGraph(
                    max_node_count=max_node_count,
                    aligned_site_aa_offset=gra.ALIGNED_SITE_AA_OFFSET,
                    min_pssm_score=gra.MIN_PSSM_SCORE 
                    )
        elif gra.__class__.__name__ == 'AcceptorSiteCollectionGraph':
            newgra = AlignedCbg5pBoundarySiteGraph(
                    max_node_count=max_node_count,
                    aligned_site_aa_offset=gra.ALIGNED_SITE_AA_OFFSET,
                    min_pssm_score=gra.MIN_PSSM_SCORE 
                    )
        elif gra.__class__.__name__ == 'SpliceSiteCollectionGraph':
            newgra = AlignedCbgBoundarySiteGraph(
                    max_node_count=max_node_count,
                    aligned_site_aa_offset=gra.ALIGNED_SITE_AA_OFFSET,
                    min_pssm_score=gra.MIN_PSSM_SCORE 
                    )
        else:
            message = "wrong input object encountered: %s" % gra.__class__.__name__
            raise InproperlyAppliedArgument, message

    elif gra.connectivitysaturation() == 1.0 and None in gra.phase() and len(gra.phase()) == 2:
        # Mix of CodingBlockStart/Acceptors od CodingBlockEnd/Donors CollectionGraph
        if gra.__class__.__name__ == 'DonorSiteCollectionGraph':
            newgra = AlignedCbg3pMixedSiteGraph(
                    max_node_count=max_node_count,
                    aligned_site_aa_offset=gra.ALIGNED_SITE_AA_OFFSET,
                    min_pssm_score=gra.MIN_PSSM_SCORE 
                    )
        elif gra.__class__.__name__ == 'AcceptorSiteCollectionGraph':
            newgra = AlignedCbg5pMixedSiteGraph(
                    max_node_count=max_node_count,
                    aligned_site_aa_offset=gra.ALIGNED_SITE_AA_OFFSET,
                    min_pssm_score=gra.MIN_PSSM_SCORE 
                    )
        elif gra.__class__.__name__ == 'SpliceSiteCollectionGraph':
            newgra = AlignedCbgMixedSiteGraph(
                    max_node_count=max_node_count,
                    aligned_site_aa_offset=gra.ALIGNED_SITE_AA_OFFSET,
                    min_pssm_score=gra.MIN_PSSM_SCORE 
                    )
        else:
            message = "wrong input object encountered: %s" % gra.__class__.__name__
            raise InproperlyAppliedArgument, message

    elif gra.connectivitysaturation() == 1.0 and None not in gra.phase() and len(gra.phase()) == 2:
        # Mixed phase Aligned site (allowing splice site phase shift)
        if gra.__class__.__name__ == 'DonorSiteCollectionGraph':
            newgra = AlignedDonorSiteWithPhaseShiftGraph(
                    max_node_count=max_node_count,
                    aligned_site_aa_offset=gra.ALIGNED_SITE_AA_OFFSET,
                    min_pssm_score=gra.MIN_PSSM_SCORE
                    )
        elif gra.__class__.__name__ == 'AcceptorSiteCollectionGraph':
            newgra = AlignedAcceptorSiteWithPhaseShiftGraph(
                    max_node_count=max_node_count,
                    aligned_site_aa_offset=gra.ALIGNED_SITE_AA_OFFSET,
                    min_pssm_score=gra.MIN_PSSM_SCORE
                    )
        elif gra.__class__.__name__ == 'SpliceSiteCollectionGraph':
            newgra = AlignedSpliceSiteWithPhaseShiftGraph(
                    max_node_count=max_node_count,
                    aligned_site_aa_offset=gra.ALIGNED_SITE_AA_OFFSET,
                    min_pssm_score=gra.MIN_PSSM_SCORE
                    )
        else:
            message = "wrong input object encountered: %s" % gra.__class__.__name__
            raise InproperlyAppliedArgument, message

    else:
        # not an aligned site, but still a collection
        return gra


    # complete the new graph
    newgra._type        = gra._type
    newgra._node_pssm   = gra._node_pssm
    newgra._node_object = gra._node_object
    newgra.nodes        = gra.nodes
    newgra.weights      = gra.weights
    newgra._edge_binary_entropies = gra._edge_binary_entropies
    return newgra

# end of function SpliceSiteCollectionGraph2AlignedSpliceSiteGraph


def TranslationalStartSiteCollectionGraph2AlignedTranslationalStartSiteGraph(gra,max_node_count=None):
    """
    """
    if gra.node_count() == 1 and gra.edge_count() == 0:
        return gra
    elif gra.connectivitysaturation() == 1.0:
        # make the new graph
        newgra = AlignedTranslationalStartSiteGraph(
            max_node_count=max_node_count,
            aligned_site_aa_offset=gra.ALIGNED_SITE_AA_OFFSET,
            min_pssm_score=gra.MIN_PSSM_SCORE 
            )
        newgra._node_pssm   = gra._node_pssm
        newgra._node_object = gra._node_object
        newgra.nodes        = gra.nodes
        newgra.weights      = gra.weights
        newgra._edge_binary_entropies = gra._edge_binary_entropies
        newgra._tcode5pscore = gra._tcode5pscore
        newgra._tcode3pscore = gra._tcode3pscore
        newgra._TCODE_5P_WINDOWSIZE = gra._TCODE_5P_WINDOWSIZE
        newgra._TCODE_3P_WINDOWSIZE = gra._TCODE_3P_WINDOWSIZE
        return newgra
    else:
        return gra

# end of function TranslationalStartSiteCollectionGraph2AlignedTranslationalStartSiteGraph


def PacbpCollectionGraph2CodingBlockGraph(pcg):
    """
    Convert PacbpCollectionGraph 2 CodingBlockGraph

    @attention: function just converts, error check is not performed here!

    @type  pcg: PacbpCollectionGraph
    @param pcg: PacbpCollectionGraph instance

    @rtype:  CodingBlockGraph
    @return: CodingBlockGraph instance
    """
    from graph_codingblock import CodingBlockGraph
    cbg = CodingBlockGraph()
    cbg.nodes   = pcg.nodes
    cbg.weights = pcg.weights
    cbg.pacbps  = pcg.pacbps
    return cbg

# end of function PacbpCollectionGraph2CodingBlockGraph


def CodingBlockGraph2PacbpCollectionGraph(cbg):
    """
    Convert CodingBlockGraph 2 PacbpCollectionGraph -> a backwards conversion!

    @attention: function just converts, error check is not performed here!

    @type  cbg: CodingBlockGraph 
    @param cbg: CodingBlockGraph instance

    @rtype:  PacbpCollectionGraph 
    @return: PacbpCollectionGraph instance
    """
    from graph_pacbpcollection import PacbpCollectionGraph 
    pcg = PacbpCollectionGraph()
    pcg.nodes   = deepcopy(cbg.nodes)
    pcg.weights = deepcopy(cbg.weights)
    pcg.pacbps  = deepcopy(cbg.pacbps)
    return pcg 

# end of function CodingBlockGraph2PacbpCollectionGraph


def CodingBlockGraph2GeneTreeGraph(cbg):
    """
    Convert CodingBlockGraph 2 GeneTree

    @attention: function just converts, error check is not performed here!

    @type  cbg: CodingBlockGraph
    @param cbg: CodingBlockGraph instance

    @rtype:  GeneTreeGraph
    @return: GeneTreeGraph instance
    """
    gtg = GeneTreeGraph()
    cbgnode2orgnode = {}
    for node in cbg.get_nodes():
        org = cbg._organism_from_node(node)
        gtg.add_node(org)
        # add node/org combi to mapping dict
        cbgnode2orgnode[ node ] = org
    # now add all the edges
    omsr = cbg.overall_minimal_spanning_range()
    for (n1,n2) in cbg.pairwisecrosscombinations_node():
        if cbg.has_edge(n1,n2):
            # get pacbp(orf) object
            thepacbp = cbg.get_pacbps_by_nodes(node1=n1,node2=n2)[0]
            # get relative coordinates of the OMSR part of the alignment
            omsrQs = thepacbp.alignmentposition_by_query_pos( min( omsr[n1] ) )
            omsrQe = thepacbp.alignmentposition_by_query_pos( max( omsr[n1] ) )

            # CHECK these coordinates; pacb.exceptions.CoordinateOutOfRange can occur
            # in freaky cases. They shouldn't, but do without discovered reason.
            # However, in the majority of cases, it is just a 1/few aa offset, which
            # can be easily corrected here.
            if str(omsrQs) == str(pacb.exceptions.CoordinateOutOfRange):
                if thepacbp.__class__.__name__ == 'PacbP':
                    # solve by taking thepacbp.query_start 
                    omsrQs = thepacbp.alignmentposition_by_query_pos( thepacbp.query_start )
                else:
                    # thepacbp.__class__.__name__ in ['PacbPDNA','PacbPORF']:
                    # solve by taking orginal alignment position start
                    omsrQs = thepacbp.alignmentposition_by_query_pos(
                        thepacbp._get_original_alignment_pos_start().query_pos
                        )

                ###########################################################################
                ## print warning message(s)
                #print "WARNING: pacb.exceptions.CoordinateOutOfRange (omsrQs, ", 
                #print "node %s in CodingBlockGraph2GeneTreeGraph" % ( str(n1) )
                #print "WARNING: min(omsr(", min( omsr[n1] ), ")", min(omsr[n1]),
                #print max(omsr[n1]), " taken ->", thepacbp.query_start, omsrQs
                #print "WARNING: ", thepacbp
                ###########################################################################

            if str(omsrQe) == str(pacb.exceptions.CoordinateOutOfRange):
                if thepacbp.__class__.__name__ == 'PacbP':
                    # solve by taking thepacbp.query_end
                    omsrQe = thepacbp.alignmentposition_by_query_pos( thepacbp.query_end )
                else:
                    # thepacbp.__class__.__name__ in ['PacbPDNA','PacbPORF']:
                    # solve by taking orginal alignment position end 
                    omsrQe = thepacbp.alignmentposition_by_query_pos(
                        thepacbp._get_original_alignment_pos_end().query_pos
                        ) + 1  # add +1 to create a python list range coordinate

                ###########################################################################
                ## print warning message(s)
                #print "WARNING: pacb.exceptions.CoordinateOutOfRange (omsrQe, ",
                #print node %s in CodingBlockGraph2GeneTreeGraph" % ( str(n1) )
                #print "WARNING: max(omsr(", max( omsr[n1] ), ")", min(omsr[n1]),
                #print max(omsr[n1]), " taken ->", thepacbp.query_end, omsrQe
                #print "WARNING: ", thepacbp
                ###########################################################################

            else:
                # omsrQe was nicely an integer; add +1 because max(OMSR) is not a range coord
                omsrQe += 1

            # calculate identityscore
            identityscore = pacb.calculate_identityscore( thepacbp.alignment[omsrQs:omsrQe] )
        else:
            # this edge is absent in the CBG!
            # TODO -> this will cause a crash a few lines later
            # by definition, a CBG MUST HAVE ALL EDGES at this stage!
            print "about to crash!!!!"
            print cbg
            print cbg.node_count(), cbg.edge_count(), "missing:", (n1,n2) 
            identityscore = 0.0
        # get organism identifyers from node and add edge
        o1,o2 = cbgnode2orgnode[ n1 ], cbgnode2orgnode[ n2 ]

        # Wt used is identityscore == Identity + 0.5* Similarity
        gtg.add_edge( o1, o2, wt=identityscore )

        # add additional statistics to gtg object. Wt used is
        # identitypercentage is TRUE aa indentity %
        identityperc = pacb.calculate_identity( thepacbp.alignment[omsrQs:omsrQe] )
        gtg._aa_identity_percentages[(o1,o2)] = identityperc
        gtg._aa_identity_percentages[(o2,o1)] = identityperc

        # bitscoreratio is ratio of bits / max bits
        bitscoreratio = pacb.calculate_bitscoreratio(
                thepacbp.query[omsrQs:omsrQe],
                thepacbp.sbjct[omsrQs:omsrQe],
                matrix = thepacbp.MATRIX
                )
        gtg._bitscore_ratios[(o1,o2)] = bitscoreratio
        gtg._bitscore_ratios[(o2,o1)] = bitscoreratio

      
        # ntidentity is obviously nt identity%
        dnaQseq, dnaSseq = thepacbp.get_unextended_aligned_dna_sequences()
        ntidentity = sequence_identity_ratio(dnaQseq,dnaSseq)
        gtg._nt_identity_percentages[(o1,o2)] = ntidentity
        gtg._nt_identity_percentages[(o2,o1)] = ntidentity

    # check if the graph is saturated (complete)
    # if not (organism/node/orf missing), add this as a zero-wt edge
    gtg.makecompletegraph(wt=0.0)
    # and return this new genetree graph
    return gtg

# end of function CodingBlockGraph2GeneTreeGraph


def LowSimilarityRegionCodingBlockGraph2GeneTreeGraph(cbg):
    """
    Convert LowSimilarityRegion 2 GeneTree

    @attention: function just converts, error check is not performed here!

    @type  cbg: LowSimilarityRegion
    @param cbg: LowSimilarityRegion instance

    @rtype:  GeneTreeGraph
    @return: GeneTreeGraph instance
    """
    gtg = GeneTreeGraph()
    cbgnode2orgnode = {}
    for node in cbg.get_nodes():
        org = cbg._organism_from_node(node)
        gtg.add_node(org)
        # add node/org combi to mapping dict
        cbgnode2orgnode[ node ] = org
    # now add all the edges
    omsr = cbg.overall_minimal_spanning_range()
    for (n1,n2) in cbg.pairwisecrosscombinations_node():
        if cbg.has_edge(n1,n2):
            # get pacbp(orf) object
            pacbps = cbg.get_pacbps_by_nodes(node1=n1,node2=n2)
            if pacbps:
                identityscore = pacbps[0].identityscore
            else:
                # this edge has no pacbp in the lsrCBG -> happens often
                identityscore = 0.0
        else:
            # this edge is absent in the lsrCBG!
            identityscore = 0.0
        # get organism identifyers from node and add edge
        o1,o2 = cbgnode2orgnode[ n1 ], cbgnode2orgnode[ n2 ]
        gtg.add_edge( o1, o2, wt=identityscore )
    # check if the graph is saturated (complete)
    # if not (organism/node/orf missing), add this as a zero-wt edge
    gtg.makecompletegraph(wt=0.0)
    # and return this new genetree graph
    return gtg

# end of function LowSimilarityRegionCodingBlockGraph2GeneTreeGraph



def pacbpCollection2AcceptedCodingBlockGraphs(pacbpCollection,gtg=None,prev=None,next=None,
    max_cbg_gtg_topo_dif=None,
    max_cbg_gtg_abs_dif=None,
    min_cbg_gtg_id_ratio=None):
    """
    """
    # make splitted subgraphs from PacbpCollection
    # cbgs must have collection.organism_set_size()-1 nodes
    # and no missing edges are aloued
    # to the number of nodes missing in the input splittedCBG)
    exact_cbg_node_count = pacbpCollection.organism_set_size()
    exact_cbg_edge_count = exact_cbg_node_count - 1
    dpcPacbpCollection = deepcopy(pacbpCollection)
    splitted_subgraphs = pacbpCollection.find_fully_connected_subgraphs(
                edges=exact_cbg_edge_count,
                max_missing_edges=0
                )


    # get pacbps for the splitted subgraphs and update edge weights
    completed_subgraphs = []
    for spl in splitted_subgraphs:
        # only deal with complete CBGs, not incomplete or collections
        if spl.node_count() != exact_cbg_node_count: continue
        if spl.__class__.__name__ == 'PacbpCollectionGraph': continue
        if spl.connectivitysaturation() < 1.0: continue

        # harvest pacbps from the deepcopied PacbpCollection
        spl.harvest_pacbps_from_pacbpcollection(dpcPacbpCollection)
        if not spl.has_overall_minimal_spanning_range(): continue
        if not spl.has_all_pacbps(): continue
        spl.update_edge_weights_by_minimal_spanning_range()
        completed_subgraphs.append(spl)

    # order graphs by total weight
    completed_subgraphs = ordering.order_graphlist_by_total_weight(completed_subgraphs)
    # and re-order on node occurrence: if a neighboring node is incorporated -> more likely!
    completed_subgraphs = ordering.reorder_cbgs_on_node_occurrence(completed_subgraphs,prev=prev,next=next)

    accepted_cbgs = []
    for spl in completed_subgraphs:
        if gtg:
            if max_cbg_gtg_topo_dif:
                topo_dif = gtg.graphalignmentdifference( spl.genetree() )
                if topo_dif > max_cbg_gtg_topo_dif:
                    continue
            if max_cbg_gtg_abs_dif:
                abs_dif  = gtg.absolutegraphalignmentdifference( spl.genetree() )
                if abs_dif > max_cbg_gtg_abs_dif:
                    continue
            if min_cbg_gtg_id_ratio:
                identity_ratio = spl.genetree().identity() / gtg.identity()
                if identity_ratio < min_cbg_gtg_id_ratio:
                    continue
        # if this point is reached: splitted cbg is accepted!
        accepted_cbgs.append(spl)

    # return the accepted_cbgs
    return accepted_cbgs

# end of function pacbpCollection2AcceptedCodingBlockGraphs 


def sequence_identity_ratio(seq1,seq2):
    """
    """
    cnt_ident = 0.0
    len1, len2 = len(seq1),len(seq2)
    for i in range(0,min([len1, len2])):
        if seq1[i].lower() == seq2[i].lower():
            cnt_ident += 1.0
    return cnt_ident / max([len1, len2]) 

# end of function sequence_identity_ratio

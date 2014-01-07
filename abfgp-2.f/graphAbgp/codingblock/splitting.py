"""
Functions for splitting of CodingBlockGraphs
"""

# Module metadata
__authors__ = "Ate van der Burgt"
__license__ = "MIT"

# General graphAbgp Imports
from graphAbgp.lowsimilaritycodingblock import LowSimilarityRegionCodingBlockGraph
from exceptions import InproperlyAppliedArgument
from recombination import pairwise, cross 

# CodingBlockGraph functionality Imports
from spanningrangedifference import (
    spanningrange_difference,
    has_spanningrange_difference, 
    has_left_spanningrange_difference,
    has_rigth_spanningrange_difference,
    )

# Abgp imports
from abgp_warnings import UnexpectedEventWarning
from gene.intron import _filter_intron_list
from lib_stopwatch import StopWatch
import pacb

# Python Imports
from sets import Set
from copy import deepcopy

# Global variables
from settings.genestructure import MIN_INTRON_NT_LENGTH
from settings.codingblockgraph import (
    CBG_MIN_AA_LENGTH,
    CBG_SPRDIF_MIN_NODE_COUNT,
    # inframe intron splitting settings
    INFRAME_INTRON_IN_CBG_GAP_SIZE,
    INFRAME_INTRON_IN_CBG_GAP_SIZE_ALONE,
    INFRAME_INTRON_IN_CBG_LENGTH_DISCREPANCY,
    INFRAME_INTRON_IN_CBG_WINDOW_AA_SIZE,
    INFRAME_INTRON_IN_CBG_WINDOW_MIN_SIMILARITY_SCORE,
    INFRAME_INTRON_IN_CBG_OBSERVED_VS_EXPECTED,
    INFRAME_INTRON_IN_CBG_MIN_TOTAL_PSSM,
    )


def split_codingblock_on_spanningrange_difference(cbg,
    verbose=False,side=None,**kwargs):
    """
    Split this CodingBlockGraph on its spanningrange difference

    @attention: see the spanningrange_difference() for arguments documentation

    @type  cbg: CodingBlockGraph
    @param cbg: CodingBlockGraph instance

    @rtype:  list
    @return: [ cbgUNSPLITTED ], [ cbgL, cbgR ] or [ cbgL, lsrCBG, cbgR ]
    """
    # obtain the spanningrange_difference
    coords = spanningrange_difference(cbg,side=side,**kwargs)

    # if no SPRDIF of the requested properties -> no split possibke
    if not coords: return [ cbg ]

    # if here, start splitting the codingblockgraph
    bckp_cbg = deepcopy(cbg)
    cbg.clear_cache()

    ####################################################################
    if verbose: print "sprdifCBG splitting", side, "coords::", coords
    ####################################################################

    if side in ['left','5p']:
        # and make a split of on the left side of the OMSR
        (cbgL,cbgR,lowsimilarityedges) = _split_on_left_sprdif(cbg,coords)
        if cbgL and cbgL.connectivitysaturation() < 1.0:
            ####################################################################
            #if verbose:
            #    # Old debugging information printed to STDOUT
            #    # Maybe usefull in future debugging ;-)
            #    print "splitL,cbgL.consat < 1.0" #, cbgL
            #    print cbgL.get_ordered_nodes(), cbgL.edge_count(),
            #    print cbgL.weights.keys(), cbgL.connectivitysaturation(),
            #    print cbgL.has_overall_minimal_spanning_range()
            #    cbgL.printmultiplealignment()
            #    print "-"*60
            #    for k,pacbporf in cbgL.pacbps.iteritems():
            #        pacbp = pacb.conversion.pacbporf2pacbp(pacbporf)
            #        print k
            #        print pacbp
            #        pacbp.print_protein(_linesize=100)
            #    print "-"*60
            #    print bckp_cbg
            #    for k,pacbporf in bckp_cbg.pacbps.iteritems():
            #        pacbp = pacb.conversion.pacbporf2pacbp(pacbporf)
            #        print k
            #        print pacbp
            #        pacbp.print_protein(_linesize=100)
            #    print "-"*60
            #    print ""
            ####################################################################
            # Create edges and remove the weakest connected node
            # This situation occurs when only a few of the PACBPS
            # support the sprdif area
            cbgL.make_pacbps_for_missing_edges(use_maxsr=True)
            weakestnode = cbgL.weakest_connected_node()
            cbgL.del_node(weakestnode)
            cbgL._update_after_changes()

    elif side in ['rigth','3p']:
        # and make a split of on the rigth side of the OMSR
        (cbgL,cbgR,lowsimilarityedges) = _split_on_rigth_sprdif(cbg,coords)
        if cbgR and cbgR.connectivitysaturation() < 1.0:
            ####################################################################
            #if verbose:
            #    # Old debugging information printed to STDOUT
            #    # Maybe usefull in future debugging ;-)
            #    print "splitR,cbgR.consat < 1.0" #, cbgR
            #    print cbgR.get_ordered_nodes(), cbgR.edge_count(),
            #    print cbgR.weights.keys(), cbgR.connectivitysaturation(),
            #    print cbgR.has_overall_minimal_spanning_range()
            #    cbgR.printmultiplealignment()
            #    print "-"*60
            #    for k,pacbporf in cbgR.pacbps.iteritems():
            #        pacbp = pacb.conversion.pacbporf2pacbp(pacbporf)
            #        print k
            #        print pacbp
            #        pacbp.print_protein(_linesize=100)
            #    print "-"*60
            #    print bckp_cbg 
            #    for k,pacbporf in bckp_cbg.pacbps.iteritems():
            #        pacbp = pacb.conversion.pacbporf2pacbp(pacbporf)
            #        print k
            #        print pacbp
            #        pacbp.print_protein(_linesize=100)
            #    print "-"*60
            #    print ""
            ####################################################################
            # Create edges and remove the weakest connected node
            # This situation occurs when only a few of the PACBPS
            # support the sprdif area
            cbgR.make_pacbps_for_missing_edges(use_maxsr=True)
            weakestnode = cbgR.weakest_connected_node()
            cbgR.del_node(weakestnode)
            cbgR._update_after_changes()

    else:
        # incorrectly applied side argument; exception already
        # caught in spanningrange_difference() function
        pass

    # check splitting status; if cbgL or cbgR == None,
    # then a pacb.exceptions.CoordinateOutOfRange occurred

    if not _is_compatible_sprdif_split(cbgL,cbgR,side=side,**kwargs):
        # just return the original cbg
        return [ bckp_cbg ]

    else:
        # yes, succesfull split!
        cbgL._update_after_changes()
        cbgR._update_after_changes()
        cbgL.update_edge_weights_by_minimal_spanning_range()
        cbgR.update_edge_weights_by_minimal_spanning_range()
        # update the IS_SPLITTED status of the split
        cbgL.IS_LAST        = False
        cbgL.IS_SPLITTED    = True
        cbgL.IS_3P_SPLITTED = True
        cbgR.IS_FIRST       = False
        cbgR.IS_SPLITTED    = True
        cbgR.IS_5P_SPLITTED = True

        # check for lowsimilarityedges
        if lowsimilarityedges:
            # create an intermediate LowSimilarityRegion CBG
            lsrCBG = create_intermediate_lowsimilarity_region(cbgL,cbgR)
            # return the result: 3 splitted graphs!
            return [ cbgL, lsrCBG, cbgR ]
        else:
            return [ cbgL, cbgR ]

# end of function split_codingblock_on_spanningrange_difference


def iteratively_split_codingblock_on_spanningrange_difference(cbg,**kwargs):
    """
    Split this CodingBlockGraph iteratively on its spanningrange difference

    @attention: see the spanningrange_difference() for arguments documentation

    @type  cbg: CodingBlockGraph
    @param cbg: CodingBlockGraph instance

    @rtype:  list
    @return: list of one (no split) or more (succesfull split) CBGs
    """
    results = [ cbg ]
    split_was_succesfull = True
    split_count = 0

    while split_was_succesfull:
        # split the CBG on SPRDIF
        splits = split_codingblock_on_spanningrange_difference(cbg,**kwargs)

        if len(splits) == 1:
            # HARD-reset splits[0] to results; Python is call-by-reference
            # and `cbg` can be changed in split_codingblock_on_etc..
            if side in ['left','5p']:
                results[0] = splits[0] # set back to FIRST -> left splitting
            else:
                results[-1] = splits[0] # set back to LAST -> rigth splitting
            split_was_succesfull = False
        else:
            if side in ['left','5p']:
                # split on the left sprdif; splits is either:
                # a 2-sized list of 2 CBGs
                # a 3-sized list of a CBG, a lsrCBG and a CBG
                # the possible lsrCBG is dealth with lateron
                mostleft = splits[0]
                left     = splits[-1]
                if left.IS_FIRST:
                    mostleft.IS_FIRST   = True
                mostleft.IS_LAST        = False
                mostleft.IS_SPLITTED    = True
                mostleft.IS_3P_SPLITTED = True
                left.IS_FIRST           = False
                left.IS_SPLITTED        = True
                left.IS_5P_SPLITTED     = True
                if split_count >= 1:
                    left.IS_3P_SPLITTED = True
                # update result list
                results[0] = left
                # check if there is a lsrCBG in between the 2 CBGs in `splits`
                if len(splits) == 3:
                    lsrCBG = splits[1]
                    results.insert(0,lsrCBG)
                results.insert(0,mostleft)
                # increase split count
                split_count+=1
                # and set mostleft to the new to-be-splitted cbg
                cbg = mostleft
            else:
                # split on the rigth sprdif; splits is either:
                # a 2-sized list of 2 CBGs
                # a 3-sized list of a CBG, a lsrCBG and a CBG
                # the possible lsrCBG is dealth with lateron
                rigth     = splits[0]
                mostrigth = splits[-1]
                if rigth.IS_LAST:
                    mostrigth.IS_LAST   = True
                mostrigth.IS_FIRST      = False
                mostrigth.IS_SPLITTED   = True
                mostrigth.IS_5P_SPLITTED= True
                rigth.IS_LAST           = False
                rigth.IS_SPLITTED       = True
                rigth.IS_3P_SPLITTED    = True
                if split_count >= 1:
                    rigth.IS_5P_SPLITTED = True
                # update result list
                results[-1] = rigth
                # check if there is a lsrCBG in between the 2 CBGs in `splits`
                if len(splits) == 3:
                    lsrCBG = splits[1]
                    results.append(lsrCBG)
                results.append(mostrigth)
                # increase split count
                split_count+=1
                # and set mostrigth to the new to-be-splitted cbg
                cbg = mostrigth

    # return the results; a list of >=1 (splitted) CBG(s)
    return results
                
# end of function iteratively_split_codingblock_on_spanningrange_difference


def _is_compatible_sprdif_split(cbgL,cbgR,side=None,
    correct_sprdif_all_nodes=True,
    sprdif_min_node_count=CBG_SPRDIF_MIN_NODE_COUNT,
    sprdif_min_aa_length=CBG_MIN_AA_LENGTH):
    """
    Helper function for split_codingblock_on_spanningrange_difference()

    @attention: see the spanningrange_difference() for arguments documentation

    @type  cbgL: CodingBlockGraph
    @param cbgL: left CodingBlockGraph after the SPRDIF split

    @type  cbgR: CodingBlockGraph
    @param cbgR: rigth CodingBlockGraph after the SPRDIF split

    @rtype:  Boolean
    @return: True or False
    """
    # check splitting status; if cbgL or cbgR == None,
    # then a pacb.exceptions.CoordinateOutOfRange occurred
    if not cbgL or not cbgR:
        return False

    if not cbgL.connectivitysaturation() == 1.0:
        return False

    if not cbgR.connectivitysaturation() == 1.0:
        return False

    # check if the split was indeed succesfull
    checkL1 = cbgL.node_count() >= sprdif_min_node_count
    checkL2 = cbgL.is_complete()
    checkL3 = cbgL.has_overall_minimal_spanning_range()
    checkL4 = min( cbgL.overall_minimal_spanning_range_sizes().values() ) >=\
            sprdif_min_aa_length 
    checkR1 = cbgR.node_count() >= sprdif_min_node_count
    checkR2 = cbgR.is_complete()
    checkR3 = cbgR.has_overall_minimal_spanning_range()
    checkR4 = min( cbgR.overall_minimal_spanning_range_sizes().values() ) >=\
            sprdif_min_aa_length 

    # do not perform these checks on the unsplitted part
    if side in ['left','5p']:
        checkR2, checkR4 = True, True
        # confirm if the original cbg still contains the same number of nodes
        checkX = cbgR.node_count() == bckp_cbg.node_count()
    else:
        checkL2, checkL4 = True, True
        # confirm if the original cbg still contains the same number of nodes
        checkX = cbgL.node_count() == bckp_cbg.node_count()

    # Judgement day ;-) ...
    if False in [ checkL1, checkL2, checkL3, checkL4, checkX,
    checkR1, checkR2, checkR3, checkR4 ]:
        # any of the checks turned out to be False -> incompatible split!
        return False
    else:
        # all checks okay -> compatible split
        return True

# end of function _is_compatible_sprdif_split


def _split_on_left_sprdif(cbg,coords,
    minimal_lowsimilarity_aa_size=0,verbose=False):
    """
    Helper function for ``split_codingblock_on_spanningrange_difference``

    @type  cbg: CodingBlockGraph
    @param cbg: CodingBlockGraph instance

    @type  coords: dictionary
    @param coords: spanningrange difference; Organism identifier as keys,
                    lists of spanningrange coordinates as values

    @type  minimal_lowsimilarity_aa_size: integer
    @param minimal_lowsimilarity_aa_size: smallest length of lsrCBG stretch
                    in a single Organism/Gene

    @type  verbose: Boolean
    @param verbose: print debugging intermediates to STDOUT

    @rtype:  tuple
    @return: ( cbgLeftSPRDIF, cbgOMSR, lowsimilarityedges )
    """
    # create a new/empty CodingBlockGraph
    from graphAbgp import CodingBlockGraph
    newcbg = CodingBlockGraph()
    newcbg.add_nodes(coords.keys())

    # check variable for CoordinateOutOfRange exception and lsr edges dict
    COORDINATE_OUTOFRANGE_EXCEPTION_OCCURRED = False
    lowsimilarityedges = {}

    for (node1,node2) in newcbg.pairwisecrosscombinations_node():
        if cbg.has_edge(node1,node2):
            # get pacbp of this edge and make a deepcopy for the other side
            rigth_pacbp = cbg.get_pacbps_by_nodes(node1=node1,node2=node2)[0]
            left_pacbp  = deepcopy(rigth_pacbp)

            # get the split coords and make a deepcopy for the other side
            splitcoordL = rigth_pacbp.alignmentposition_by_query_pos(
                    max(coords[node1])+1 )
            # make a deepcopy for coordR because the value of L and R
            # can be changed when in a lowsimilarity region
            splitcoordR = deepcopy(splitcoordL)
            # store the original splitcoord for lowsimilarityedges checking
            oriSplitCoord = deepcopy(splitcoordL)

            if splitcoordL == pacb.exceptions.CoordinateOutOfRange:
                # Hmm... unexpected CoordinateOutOfRange execption
                # -> hard-print to STDOUT for debugging
                message = "pacb.exceptions.CoordinateOutOfRange"
                outofrangeWarning = UnexpectedEventWarning(message)
                print outofrangeWarning
                COORDINATE_OUTOFRANGE_EXCEPTION_OCCURRED = True
                break

            # correct split coords when placed in a region with gaps
            while left_pacbp._positions[ splitcoordL - 1 ].isa_gap:
                splitcoordL -= 1
            while rigth_pacbp._positions[ splitcoordR ].isa_gap:
                splitcoordR += 1

            # check if an intermediate low-similarity CBG must be created
            if splitcoordL != splitcoordR:
                if splitcoordR - splitcoordL > minimal_lowsimilarity_aa_size:
                    # a large enough lowsimilarity region
                    lowsimilarityedges[(node1,node2)] = (splitcoordL,splitcoordR)
                else:
                    # no, lowsimilarityregion is to small; reset the correction
                    splitcoordL = oriSplitCoord
                    splitcoordR = oriSplitCoord

            # reset _original_alignment_pos for left and rigth pacbp
            left_pacbp._original_alignment_pos_end = splitcoordL
            rigth_pacbp._original_alignment_pos_start = splitcoordR

            # check if the left Pacbp still has a Non-zero length
            # a zero-length pacbp can occur when the left sprdif is
            # not fully connected!
            if left_pacbp._original_alignment_pos_end <=\
            left_pacbp._original_alignment_pos_start:
                pass
            else:
                # add edge and pacbporf to the new graph
                left_pacbp._score_alignment()
                left_key  = left_pacbp.construct_unique_key(node1,node2)
                newcbg.add_edge(node1,node2,wt=left_pacbp.bits)
                newcbg.pacbps[(left_key,node1,node2)] = left_pacbp

            # replace old pacbp by new, rigth pacbp in the main graph
            rigth_pacbp._score_alignment()
            rigth_key = rigth_pacbp.construct_unique_key(node1,node2)
            cbg.del_edge(node1,node2)
            cbg.add_edge(node1,node2,wt=rigth_pacbp.bits)
            cbg.pacbps[(rigth_key,node1,node2)] = rigth_pacbp

    if COORDINATE_OUTOFRANGE_EXCEPTION_OCCURRED:
        # return None in place of newcbg
        return (None, cbg, {})

    # check if all nodes have edges
    for node in coords.keys():
        if not newcbg.nodes[node]:
            newcbg.del_node(node)

    # check if nodes remaining and if newcbg has OMSR
    if newcbg.node_count() < 2:
        newcbg = None
        lowsimilarityedges = {}
    elif not newcbg.has_overall_minimal_spanning_range():
        newcbg = None
        lowsimilarityedges = {}
    else:
        pass

    ############################################################################
    if verbose:
        print "_lspl coords:",
        print [ (k,(len(v),min(v),"-",max(v))) for k,v in coords.iteritems() ]
        print "_lspl newcbg:", newcbg
        print "_lspl lsredges:", lowsimilarityedges
        print "_lspl inputcbg:", cbg
    ############################################################################

    # return the new left, the existing rigth and lowsimilarityedges
    return (newcbg,cbg,lowsimilarityedges)

# end of function _split_on_left_sprdif


def _split_on_rigth_sprdif(cbg,coords,
    minimal_lowsimilarity_aa_size=0,verbose=False):
    """
    Helper function for ``split_codingblock_on_spanningrange_difference``

    @type  cbg: CodingBlockGraph
    @param cbg: CodingBlockGraph instance

    @type  coords: dictionary
    @param coords: spanningrange difference; Organism identifier as keys,
                    lists of spanningrange coordinates as values

    @type  minimal_lowsimilarity_aa_size: integer
    @param minimal_lowsimilarity_aa_size: smallest length of lsrCBG stretch
                    in a single Organism/Gene

    @type  verbose: Boolean
    @param verbose: print debugging intermediates to STDOUT

    @rtype:  tuple
    @return: ( cbgOMSR, cbgRigthSPRDIF, lowsimilarityedges )
    """
    # create a new/empty CodingBlockGraph
    from graphAbgp import CodingBlockGraph
    newcbg = CodingBlockGraph()
    newcbg.add_nodes(coords.keys())

    # check variable for CoordinateOutOfRange exception and lsr edges dict
    COORDINATE_OUTOFRANGE_EXCEPTION_OCCURRED = False
    lowsimilarityedges = {}

    for (node1,node2) in newcbg.pairwisecrosscombinations_node():
        if cbg.has_edge(node1,node2):
            # get pacbp of this edge and make a deepcopy for the other side
            rigth_pacbp = cbg.get_pacbps_by_nodes(node1=node1,node2=node2)[0]
            left_pacbp  = deepcopy(rigth_pacbp)

            # get the split coords and make a deepcopy for the other side
            splitcoordL = rigth_pacbp.alignmentposition_by_query_pos(
                    min( coords[node1]) )
            # make a deepcopy for coordR because the value of L and R
            # can be changed when in a lowsimilarity region
            splitcoordR = deepcopy(splitcoordL)
            # store the original splitcoord for lowsimilarityedges checking
            oriSplitCoord = deepcopy(splitcoordL)

            if splitcoordL == pacb.exceptions.CoordinateOutOfRange:
                # Hmm... unexpected CoordinateOutOfRange execption
                # -> hard-print to STDOUT for debugging
                message = "pacb.exceptions.CoordinateOutOfRange"
                outofrangeWarning = UnexpectedEventWarning(message)
                print outofrangeWarning
                COORDINATE_OUTOFRANGE_EXCEPTION_OCCURRED = True
                break

            # correct split coords when placed in a region with gaps
            while left_pacbp._positions[ splitcoordL - 1 ].isa_gap:
                splitcoordL -= 1
            while rigth_pacbp._positions[ splitcoordR ].isa_gap:
                splitcoordR += 1

            # check if an intermediate low-similarity CBG must be created
            if splitcoordL != splitcoordR:
                if splitcoordR - splitcoordL > minimal_lowsimilarity_aa_size:
                    # a large enough lowsimilarity region
                    lowsimilarityedges[(node1,node2)] = (splitcoordL,splitcoordR)
                else:
                    # no, lowsimilarityregion is to small; reset the correction
                    splitcoordL = oriSplitCoord
                    splitcoordR = oriSplitCoord

            # reset _original_alignment_pos for left and rigth pacbp
            left_pacbp._original_alignment_pos_end = splitcoordL
            rigth_pacbp._original_alignment_pos_start = splitcoordR

            # check if the rigth Pacbp still has a Non-zero length
            # a zero-length pacbp can occur when the rigth sprdif is
            # not fully connected!
            if rigth_pacbp._original_alignment_pos_end <=\
            rigth_pacbp._original_alignment_pos_start:
                pass
            else:
                # add edge and pacbporf to the new graph
                rigth_pacbp._score_alignment()
                rigth_key = rigth_pacbp.construct_unique_key(node1,node2)
                newcbg.add_edge(node1,node2,wt=rigth_pacbp.bits)
                newcbg.pacbps[(rigth_key,node1,node2)] = rigth_pacbp

            # replace old pacbp by new, left pacbp in the main graph
            left_pacbp._score_alignment()
            left_key  = left_pacbp.construct_unique_key(node1,node2)
            cbg.del_edge(node1,node2)
            cbg.add_edge(node1,node2,wt=left_pacbp.bits)
            cbg.pacbps[(left_key,node1,node2)] = left_pacbp


    if COORDINATE_OUTOFRANGE_EXCEPTION_OCCURRED:
        # return None in place of newcbg
        return (cbg, None, {})

    # check if all nodes have edges
    for node in coords.keys():
        if not newcbg.nodes[node]:
            newcbg.del_node(node)

    # check if nodes remaining and if newcbg has OMSR
    if newcbg.node_count() < 2:
        newcbg = None
        lowsimilarityedges = {}
    elif not newcbg.has_overall_minimal_spanning_range():
        newcbg = None
        lowsimilarityedges = {}
    else:
        pass

    ############################################################################
    if verbose:
        print "_rspl coords:",
        print [ (k,(len(v),min(v),"-",max(v))) for k,v in coords.iteritems() ]
        print "_rspl newcbg:", newcbg
        print "_rspl lsredges:", lowsimilarityedges
        print "_rspl inputcbg:", cbg
    ############################################################################

    # return the existing left, the new rigth and lowsimilarityedges
    return (cbg,newcbg,lowsimilarityedges)

# end of function _split_on_rigth_sprdif


def split_codingblock_on_alternatives_in_pacbps_dict(gra,
    filter_for_msr=True,
    filter_for_omsr=True,
    use_linearisation_filter=True,
    verbose=False):
    """
    Split this CodingBlockGraph on alternatives in its listed PacbPs
    
    @attention: only functional if >=1 edge is supported by >=1 PacbPORF object

    @type  gra: CodingBlockGraph
    @param gra: CodingBlockGraph instance

    @type  filter_for_msr: Boolean
    @param filter_for_msr: if True, return only CBGs with a MSR

    @type  filter_for_omsr: Boolean
    @param filter_for_omsr: if True, return only CBGs with an OMSR

    @type  use_linearisation_filter: Boolean
    @param use_linearisation_filter: True is default and recommended;
                filter pacbps by linearisation (not needed e.g.
                in a hmmsearch completed CBG)

    @type  verbose: Boolean
    @param verbose: print debugging intermediates to STDOUT

    @rtype:  list
    @return: list with 0 or more CodingBlockGraph instances
    """

    ############################################################################
    if verbose:
        stw = StopWatch(name="altPacbpsInCBG")
        stw.start()
    ############################################################################

    # check if PacbPs are present
    if len(gra.pacbps) == 0: return [ gra ]

    # check if edge count equals pacbp count in this CBG
    if gra.edge_count() == len(gra.pacbps):
        # nothing to split here!
        if not filter_for_omsr or\
        (filter_for_omsr and gra.has_overall_minimal_spanning_range()):
            return [ gra ]
        else:
            return [ ]

    # gather edges with >1 pacbp
    edge2pacbpkey = {}
    for (node1,node2) in gra.pairwisecrosscombinations_node():
        # check if these nodes are indeed present as an edge
        if not gra.has_edge(node1,node2): continue

        # get pacbp list and linearize
        pacbps = gra.get_pacbps_by_nodes(node1=node1,node2=node2)
        if use_linearisation_filter:
            linear, rejected = pacb.linearization.linearise_pacbp_list(
                    pacbps,
                    ACCEPTANCE_MARGIN=150,
                    start_with_bests=1,
                    is_weigthed=True,
                    order_by='length'
                    )

            # Remove rejected ones from gra.pacps!
            if rejected:
                for pacbp in rejected:
                    gra.remove_pacbp(pacbp,node1,node2)
        else:
            # no linearisation requested for!
            linear = pacbps

        # Store linear ones as keys in the edge2pacbpkey dict-of-lists
        # in case there are still at least 2 alternatives
        if len(linear) > 1:
            # make a new entry in the new_edge2pacbpkey dict
            edge2pacbpkey[(node1,node2)] = []

            if len(linear) < 3:
                # store the pacbp keys into this new entry
                for pacbp in linear:
                    edge2pacbpkey[(node1,node2)].append( pacbp )
            else:
                # try to remove additional pacbps by checking for repeats!
                withoutrepeats =\
                        pacb.linearization.remove_repeats_from_pacbp_list(
                                linear,overlap_ratio=0.90)
                # and store the pacbp keys into this new entry
                for pacbp in withoutrepeats:
                    edge2pacbpkey[(node1,node2)].append( pacbp )

    # if no duplicated edges found or remaining, return
    if not edge2pacbpkey: return [ gra ]

    ############################################################################
    # if here, then the Pacbp recombination starts that can take very long
    if verbose: stw.lap()
    ############################################################################

    incompatible = []
    orderedkeys = edge2pacbpkey.keys()
    orderedkeys.sort() 
    for (node1A,node2A) in orderedkeys:
        alternativesA = edge2pacbpkey[(node1A,node2A)]
        for (node1B,node2B) in orderedkeys:
            alternativesB = edge2pacbpkey[(node1B,node2B)]
            # ignore identical edges / node combis
            if (node1A,node2A) == (node1B,node2B): continue
            # ignore edges that do not share a single node
            if len(Set([node1A,node2A,node1B,node2B])) == 4: continue
            # check which pacbps share some overlap at their common node
            if   node1A == node1B:  combi = ('Q','Q')
            elif node2A == node1B:  combi = ('S','Q')
            elif node2A == node2B:  combi = ('S','S')
            else:                   combi = ('Q','S')
            for altposA in range(0,len(alternativesA)):
                pacbpA = alternativesA[altposA]
                if combi[0] == 'Q':
                    coordsA = pacbpA.alignment_protein_range_query()
                else:
                    coordsA = pacbpA.alignment_protein_range_sbjct()
                for altposB in range(0,len(alternativesB)):
                    pacbpB = alternativesB[altposB]
                    if combi[1] == 'Q':
                        coordsB = pacbpB.alignment_protein_range_query()
                    else:
                        coordsB = pacbpB.alignment_protein_range_sbjct()
                    # Now compare the coordinate Set ranges
                    # If there is an overlap, these pacbps can possibly/likely
                    # contribute to a CBG with OMSR
                    pacbpcombi = [( (node1A,node2A), altposA ),
                                  ( (node1B,node2B), altposB ) ]
                    pacbpcombi.sort() 
                    pacbpcombi = tuple(pacbpcombi)
                    if len( Set(coordsA).intersection(coordsB) ) <\
                    CBG_MIN_AA_LENGTH:
                        pacbpcombi = [( (node1A,node2A), altposA ),
                                      ( (node1B,node2B), altposB ) ]
                        pacbpcombi.sort()
                        pacbpcombi = tuple(pacbpcombi)
                        incompatible.append( pacbpcombi )
                    else:
                        # compatible so far
                        pass


    # Make a list with for each (ordered) edge a list with the number of PacbP
    # alternatives. 0 means there is only a single PacbP (accesible by
    # zero-based list ID!). Explained by example:
    # K(s) CBG with s=4 means 6 edges. Assume (ordered) edge 2 and 4 have 3 and
    # 2 PacbP alternatives, respectively. After the next lines of code,
    # pacbpsperorderededges will be literaly this:
    # [ [0], [0], [0,1,2], [0], [0,1], [0] ]
    # Note the PacbP alternative list IDs on element 2 [0,1,2] and 4 [0,1]
    pacbpsperorderededges = []
    for (node1,node2) in orderedkeys:
        pacbpsperorderededges.append(
            range(0,len(edge2pacbpkey[(node1,node2)]))
            )

    # Now make all unique recombinations of PacbPs by using the cross() function
    # Explained by example as introduced above: After cross(),
    # allcombis will be literaly this:
    # [
    #   [ 0, 0, 0, 0, 0, 0 ], [ 0, 0, 1, 0, 0, 0 ], [ 0, 0, 2, 0, 0, 0 ],
    #   [ 0, 0, 0, 0, 1, 0 ], [ 0, 0, 1, 0, 1, 0 ], [ 0, 0, 2, 0, 1, 0 ],
    # ]
    allcombis = cross(pacbpsperorderededges)

    ############################################################################
    if verbose:
        print stw.lap(), "incompatible", len(incompatible),
        print "allcombis:", allcombis
    ############################################################################

    # compare all combinations of pacbp alternatives, listed in `allcombis`,
    # with all pairwise pacbp combinations that are incompatible,
    # in `incompatible`
    for ((nq1,ns1),pos1),((nq2,ns2),pos2) in incompatible:
        # find position in element list of `allcombis`
        # example element in allcombis: [0,1,1,0,2,4,1,0]
        # where the index of the element means the node-combi
        # index in `orderedkeys` and the element it self means the
        # pacbp number (in the list) in `edge2pacbpkey[(node1,node2)]`
        posinlist1 = orderedkeys.index((nq1,ns1))
        posinlist2 = orderedkeys.index((nq2,ns2))
        # loop over all combinations in allcombis and check for the
        # occurrence of a pattern-like combi [any,1,any,any,any,4,any,any]
        # when this pattern is observed -> remove this element because it
        # will result in an incompatible combination!
        for combipos in range(len(allcombis)-1,-1,-1):
            if allcombis[combipos][posinlist1] == pos1 and\
            allcombis[combipos][posinlist2] == pos2:
                allcombis.pop(combipos)
        # quit if no combinations are left!
        if len(allcombis) == 0: break

    ############################################################################
    if verbose: print stw.lap(), "remaining allcombis:", allcombis
    ############################################################################

    # compare all combinations of pacbp alternatives, listed in `allcombis`,
    # with all pairwise pacbp combinations that are incompatible,
    # in `incompatible`
    compatible_combis = []
    for combi in allcombis:
        pacbpcombi = []
        for orderedpos in range(0,len(combi)):
            pacbpcombi.append( ( orderedkeys[orderedpos], combi[orderedpos] ) ) 
        # get pairwise pairs of PacbPs
        pairs = pairwise(pacbpcombi)
        for pair in pairs:
            # in any of the pairs is listed as an incompatible pair ->
            # then ignore this combi
            if pair in incompatible:
                break
        else:
            # all PacbP pairs are compatible, so the overall combi is compatible
            compatible_combis.append(combi)

    # create a basal deepcopy of the unambigious edges (no pacbp alternatives)
    basegraph = gra.deepcopy()
    for (node1,node2) in orderedkeys:
        basegraph.del_edge(node1,node2)
        for pacbp in edge2pacbpkey[(node1,node2)]:
            basegraph.remove_pacbp(pacbp,node1,node2)

    ############################################################################
    if verbose:
        print stw.lap(), "incompatible nodes removed B",
        print ( gra.node_count(), gra.edge_count() ),
        print ( basegraph.node_count(), basegraph.edge_count() )
    ############################################################################

    separated_graphs = []
    for compatible in compatible_combis:
        # create a deepcopy of the base CBG
        newgra = basegraph.deepcopy()
        for posinorderedkeys in range(0,len(compatible)):
            (node1,node2) = orderedkeys[posinorderedkeys]
            pacbp = edge2pacbpkey[(node1,node2)][compatible[posinorderedkeys]] 
            newgra.add_edge(node1,node2,wt=pacbp.bitscore)
            pacbpkey = pacbp.construct_unique_key(node1,node2)
            newgra.pacbps[pacbpkey,node1,node2] = pacbp
        if filter_for_msr and not newgra.has_minimal_spanning_range():
            continue 
        if filter_for_omsr and not newgra.has_overall_minimal_spanning_range(): 
            continue
        # if here, then the new separated graph is okay 
        separated_graphs.append( newgra )

    ############################################################################
    if verbose: print stw.lap(), "done!"
    ############################################################################

    # and return the new separated CBGs 
    return separated_graphs 

# end of function split_codingblock_on_alternatives_in_pacbps_dict



def split_hmmsearchcbg_on_lowsimilarity(inputcbg,GTG=None):
    """
    """
    cbg = deepcopy(inputcbg)
    cbg.clear_cache()
    cbg.IS_SPLITTED = False
    cbg.IS_5P_SPLITTED = False
    cbg.IS_3P_SPLITTED = False
    orfs = cbg.get_orfs_of_graph()
    cbg.pacbporfs2pacbps()
    cbgpacbpkeys = cbg.pacbps.keys()
    # split on gaps/lowsimililarity: very high chance! e.g. 2th exon in `abfgp_mgg0180`
    for (key,n1,n2) in cbgpacbpkeys:
        pacbp = cbg.pacbps[(key,n1,n2)]
        # define splitters based on GeneTreeGraph
        TODO = True
        # splitter is a tuple of: (windowsize,identity,similarity)
        # see pacb.splitting._has_splitter for documentation
        splitters = [
            (  8, 0, 0 ),
            ( 10, 1, 0 ),
            ( 12, 0, 2 ),
        ]
        splittedpacbps, is_splitted = pacb.splitting.split_pacb_on_lowsimilarities(
                pacbp,splitters=splitters,may_contain_methione=False)
        if is_splitted:
            del( cbg.pacbps[(key,n1,n2)] )
            for splittedpacb in splittedpacbps:
                splittedkey = splittedpacb.construct_unique_key(n1,n2)
                cbg.pacbps[(splittedkey,n1,n2)] = splittedpacb
    # split on lowsimililarity
    alternatives = cbg.split_codingblock_on_alternatives_in_pacbps_dict(
        filter_for_msr=True,
        filter_for_omsr=False,
        use_linearisation_filter=False # HARD-SET linearisation to False!
        )
    for alt in alternatives:
        # function extend_pacbporfs is useless because we have no ``input`` object here 
        for (k,n1,n2), pacbp in alt.pacbps.iteritems():
            orgQ = alt._organism_from_node(n1)
            orgS = alt._organism_from_node(n2)
            orfQ = orfs[orgQ][0]
            orfS = orfs[orgS][0]
            pacbporf = pacb.conversion.pacbp2pacbporf(pacbp,orfQ,orfS)
            pacbporf.extend_pacbporf_after_stops()
            alt.pacbps[(k,n1,n2)] = pacbporf
            
    # filter for OMSR and order alternatives
    # import GenestructureOfCodingBlockGraphs here; 
    # in header of script creates bogus circular references!
    from graph_genestructure import GenestructureOfCodingBlockGraphs
    partGSG = GenestructureOfCodingBlockGraphs({})
    for alt in alternatives:
        if alt.has_overall_minimal_spanning_range():
            partGSG.add_codingblock(alt)

    # add LowSimilarityRegions in between the splitted ones
    if len(partGSG) > 1:
        # order the alternatives and assign lsrCBGs in between them!
        partGSG = GenestructureOfCodingBlockGraphs({})
        partGSG.add_codingblocks(alternatives)
        for ii in range(len(partGSG)-2,-1,-1):
            lsrCBG = create_intermediate_lowsimilarity_region(
                    partGSG.codingblockgraphs[ii],
                    partGSG.codingblockgraphs[ii+1]
                    )
            partGSG.codingblockgraphs[ii].IS_SPLITTED = True
            partGSG.codingblockgraphs[ii].IS_3P_SPLITTED = True
            partGSG.codingblockgraphs[ii+1].IS_SPLITTED = True
            partGSG.codingblockgraphs[ii+1].IS_5P_SPLITTED = True
            partGSG.codingblockgraphs.insert(ii+1,lsrCBG)

        # HARD-SET no splitted 3P for last, 5P for first
        partGSG.codingblockgraphs[0].IS_5P_SPLITTED = False
        partGSG.codingblockgraphs[-1].IS_3P_SPLITTED = False

        # set the ordered & lsrCBG added list back to alternatives
        alternatives = partGSG.codingblockgraphs

    elif len(partGSG) == 1:
        alternatives = partGSG.codingblockgraphs
    else:
        alternatives = []

    # return the splitted cbgs
    return alternatives

# end of function split_hmmsearchcbg_on_lowsimilarity 


def split_codingblock_by_inframe_intron(gra,
    gap_size=INFRAME_INTRON_IN_CBG_GAP_SIZE,
    length_discrepancy=INFRAME_INTRON_IN_CBG_LENGTH_DISCREPANCY,
    gap_size_alone=INFRAME_INTRON_IN_CBG_GAP_SIZE_ALONE):
    """
    Split this CodingBlockGraph into two seperate CBGs in case an inframe intron
    is detected in one of the PacbPORFs

    @type  gra: CodingBlockGraph
    @param gra: CodingBlockGraph instance

    @type  gap_size: integer
    @param gap_size: positive integer defining the length of a continious gap
                     that triggers inframe intron splitting

    @type  gap_size_alone: integer
    @param gap_size_alone: positive integer defining the length of a continious
                           gap that directly triggers inframe intron splitting

    @type  length_discrepancy: integer
    @param length_discrepancy: positive integer defining the difference in
                               length of sbjct and query of the PacbpOrf that
                               triggers inframe intron splitting

    @rtype:  list
    @return: list wiht CodingBlockGraph instances

    @attention: gap_size and length_discrepancy trigger only if BOTH are True!
    @attention: gap_size_alone triggers on its own if True!
    @attention: at least a single CBG (the one that served as input) is returned
    """

    # make a deepcopy for later use and clear object's cache
    inputgra = deepcopy(gra)
    inputgra.clear_cache()
    inputgra.create_cache()
    gra.clear_cache()

    splitted = {}
    for (key,(org1,orf1),(org2,orf2)),pacbporf in gra.pacbps.iteritems():
        fullkey = (key,(org1,orf1),(org2,orf2))
        if pacbporf.potentially_contains_aligned_intron(
        gap_size=gap_size,length_discrepancy=length_discrepancy,
        gap_size_alone=gap_size_alone):
            splitted_pacbps, is_splitted = pacb.splitting.split_pacb_on_gaps(pacbporf,gapsize=gap_size)
            # what now; proceed if is_splitted or when there are >2 parts??
            if is_splitted and len(splitted_pacbps) >= 2:
                # check for each part if the junction is in the minimal_spanning_range
                for part in splitted_pacbps:
                    if (part.query_end in gra.overall_minimal_spanning_range(organism=org1) and\
                    part.sbjct_end in gra.overall_minimal_spanning_range(organism=org2)) or\
                    (part.query_start in gra.overall_minimal_spanning_range(organism=org1) and\
                    part.sbjct_start in gra.overall_minimal_spanning_range(organism=org2)):
                        # yes, the split took place in the minimal_spanning_range!
                        if splitted.has_key(fullkey):
                            splitted[fullkey].append( part )
                        else:
                            splitted[fullkey] = [ part ]

                # now see what we have found
                if splitted.has_key(fullkey) and len(splitted[fullkey]) >= 2:
                    # yes, a fully compatible split!
                    pass
                elif splitted.has_key(fullkey) and len(splitted[fullkey]) == 1:
                    # nope, only one part in the splitted dict
                    del( splitted[fullkey] )
                else:
                    # nope, all splits failed
                    pass

    # if there is nothing splitted, no inframe intron here
    if not splitted:
        print "INPUT GRAPH will be returned (case0)!!"
        return [ inputgra ]

    #### make a deepcopy for later use and clear object's cache
    ###inputgra = deepcopy(gra)
    ###gra.clear_cache()

    # delete pacbp that is splitted, replace by freshly splitted parts
    for (key,node1,node2), splittedpacbps in splitted.iteritems():
        del( gra.pacbps[(key,node1,node2)] )
        for pacbp in splittedpacbps:
            newkey = pacbp.construct_unique_key(node1,node2)
            gra.pacbps[(newkey,node1,node2)] = pacbp

    # split on alternatives!
    splits = split_codingblock_on_alternatives_in_pacbps_dict(gra,
            filter_for_msr=True,filter_for_omsr=True
            )
    
    # check if split was succesfull, >1 parts are required!
    if len(splits) < 2:
        # no splits or splits got fragmented; assume this is an incompatible split!
        print "INPUT GRAPH will be returned (case1)!!"
        return [ inputgra ]
    else:
        # check if the graphs are nicely orderable in a (part of a ) genestructure.
        emptyinput = {}
        # import GenestructureOfCodingBlockGraphs here;
        # in header of script creates bogus circular references!
        from graph_genestructure import GenestructureOfCodingBlockGraphs
        testGSG = GenestructureOfCodingBlockGraphs(emptyinput)
        testGSG.add_codingblocks(splits)

        # return the splittedgraphs if they are nicely orderable
        if len(testGSG) == 1:
            # not okay! there should be at least 2
            print "INPUT GRAPH will be returned (case2)!!"
            return [ inputgra ]
        else:
            # positive confirmation of a putative inframe intron!
            return testGSG.codingblockgraphs

# end of function split_codingblock_by_inframe_intron


def potentially_contains_aligned_intron(cbg,
    gap_size=INFRAME_INTRON_IN_CBG_GAP_SIZE,
    length_discrepancy=INFRAME_INTRON_IN_CBG_LENGTH_DISCREPANCY,
    intron_window_min_similarity_score=\
            INFRAME_INTRON_IN_CBG_WINDOW_MIN_SIMILARITY_SCORE,
    dif_potentially_observed_vs_expected=\
            INFRAME_INTRON_IN_CBG_OBSERVED_VS_EXPECTED,
    MIN_TOTAL_PSSM_INFRAME_INTRON=INFRAME_INTRON_IN_CBG_MIN_TOTAL_PSSM,
    verbose=False):
    """
    Check if this CodingBlockGraph potentially contains an inframe intron

    @type  cbg: CodingBlockGraph
    @param cbg: CodingBlockGraph instance

    @type  gap_size: integer
    @param gap_size: positive integer defining the length of a continious gap
                     that triggers inframe intron splitting

    @type  length_discrepancy: integer
    @param length_discrepancy: positive integer defining the difference in
                               length of sbjct and query of the PacbpOrf that
                               triggers inframe intron splitting

    @type  intron_window_min_similarity_score: float
    @param intron_window_min_similarity_score: (positive) float, the minimal score
            of the ``positional_alignment_similarity_count`` scores of the range
            of the hypothetical intron

    @type  dif_potentially_observed_vs_expected: integer
    @param dif_potentially_observed_vs_expected: TODO...

    @type  MIN_TOTAL_PSSM_INFRAME_INTRON: float
    @param MIN_TOTAL_PSSM_INFRAME_INTRON: positive float; the minimal sum of the
            donor and acceptor pssm score of the hypothetical inframe intron

    @rtype:  Boolean
    @return: True or False
    """

    # Only graphs taht have an OVERALL minimal spanning range (~=are complete)
    # are taken into account here
    if not cbg.has_overall_minimal_spanning_range():
        return False

    # check the individual pacbporfs for potential inframe introns
    checklist = []
    encountered_nodes = []
    potential_nodes = []
    omsr = cbg.overall_minimal_spanning_range()
    for (key,node1,node2), pacbporf in cbg.pacbps.iteritems():
        if pacbporf.potentially_contains_aligned_intron(
                gap_size=gap_size,length_discrepancy=length_discrepancy):
            # check if it is the sbjct or the query that is having the gap
            # and verify that the gap is indeed in the OMSR.
            # That latter check id not done yet!
            sbjctCoords = ( min(omsr[node2]),max(omsr[node2]) )
            queryCoords = ( min(omsr[node1]),max(omsr[node1]) )
            if pacbporf.query_has_gaps(gap_size=gap_size,abs_coords=queryCoords):
                encountered_nodes.append( node2 )
            if pacbporf.sbjct_has_gaps(gap_size=gap_size,abs_coords=sbjctCoords):
                encountered_nodes.append( node1 )

    if not encountered_nodes:
        # not a single pacbp with potential inframe introns!
        return False

    # check how often a node is encountered in encountered_nodes
    # a True inframe intron **should** have its node reported often
    MIN_NODE_ENCOUNTERED = cbg.organism_set_size() - 1
    MIN_NODE_ENCOUNTERED -= dif_potentially_observed_vs_expected
    ###for node in Set(encountered_nodes):
    ###    print "NODES:", node, encountered_nodes.count(node)
    intron_seen = False
    for node in Set(encountered_nodes):
        ###print node, encountered_nodes.count(node)
        if encountered_nodes.count(node) >= MIN_NODE_ENCOUNTERED:
            # Calculate the overshoot in length of the aligned sequence of
            # this organism compared to the others. Based on this,
            # spatiously large min and max intron sizes are set
            org = cbg._organism_from_node(node)
            omsr = cbg.overall_minimal_spanning_range(organism=org)
            pacbp_length_discrepancy = len(omsr) - min(cbg.overall_minimal_spanning_range_sizes().values())
            est_inframe_intron_length = pacbp_length_discrepancy * 3
            min_inframe_intron_length = max([ int( est_inframe_intron_length * 0.4 ) , MIN_INTRON_NT_LENGTH ])
            max_inframe_intron_length = int( est_inframe_intron_length * 1.6 )

            # search for potential introns that can be responsible for this event
            theorf = cbg.get_orfs_of_graph(organism=org)[0]
            introns = pacb.connecting.merge_orfs_with_intron(theorf,theorf)
            # filter introns for all outside the OMSR, to short, to long and filtered on total pssm_score
            introns = _filter_intron_list(introns,filter_by='length',criterion=max_inframe_intron_length,operator='<=')
            introns = _filter_intron_list(introns,filter_by='length',criterion=min_inframe_intron_length,operator='>=')
            introns = _filter_intron_list(introns,filter_by='donor_pos',criterion=min(omsr)*3,operator='>')
            introns = _filter_intron_list(introns,filter_by='acceptor_pos',criterion=max(omsr)*3,operator='<')
            introns = _filter_intron_list(introns,filter_by='total_pssm',criterion=MIN_TOTAL_PSSM_INFRAME_INTRON,operator='>=')

            # get the Multiple Alignment Track similarity track for this organism
            mat = cbg.get_positional_alignment_similarity_count(organism=org)

            for intron in introns:
                # calculate absolute AA coordinates of this intron
                intron_aa_start = intron.donor.pos / 3
                intron_aa_end   = intron.acceptor.pos / 3
                # calculate the 
                score = sum( [ mat[pos] for pos in range(intron_aa_start,intron_aa_end) ] )
                if score < intron_window_min_similarity_score:
                    if verbose:
                        print "INFRAME INTRON!!", node, intron, score, min_inframe_intron_length, est_inframe_intron_length, max_inframe_intron_length
                    # yep, there definately can be/is an inframe intron in thsi CBG
                    # TODO: right now, it is just delivered upon further cleavage
                    # by another function (split_codingblock_by_inframe_intron)
                    # however, the information of exactly THIS intron should be usefull.
                    # More TODO: now, we break out at the first discovery.
                    # Maybe there are even better (==lower scoring) introns available??
                    #return True
                    intron_seen = True
    else:
        return intron_seen 
    return intron_seen

# end of function potentially_contains_aligned_intron


def create_intermediate_lowsimilarity_region(cbgL,cbgR):
    """
    Create an intermediate LowSimilarityRegionCodingBlockGraph between 2 CBGs

    @type  cbgL: CodingBlockGraph
    @param cbgL: CodingBlockGraph (left / 5' of 2 ordered CBGs)

    @type  cbgR: CodingBlockGraph
    @param cbgR: CodingBlockGraph (rigth / 3' of 2 ordered CBGs)

    @rtype  lsrCBG: LowSimilarityRegionCodingBlockGraph
    @return lsrCBG: intermediate LowSimilarityRegionCodingBlockGraph
    """
    mutual_nodes = Set( cbgL.get_nodes() ).intersection(cbgR.get_nodes())

    lsrCBG = LowSimilarityRegionCodingBlockGraph()
    omsrL = cbgL.overall_minimal_spanning_range()
    omsrR = cbgR.overall_minimal_spanning_range()
    for node in mutual_nodes:
        lsomsr = range( max(omsrL[node])+1, min(omsrR[node]) )
        if lsomsr:
            # no AA coordinates for this Node/Organism;
            # both CBGs glue perfectly together.
            org = cbgL._organism_from_node(node)
            lsrCBG.add_node_and_object( node,
                    cbgL.get_orfs_of_graph(organism=org)[0] )
            lsrCBG.set_node_omsr(node,Set(lsomsr)) 

    # set IS_SPLITTED attributes correctly
    cbgL.IS_SPLITTED    = True
    cbgL.IS_3P_SPLITTED = True
    cbgR.IS_5P_SPLITTED = True
    cbgR.IS_SPLITTED    = True

    # return the lowsimilarity region CBG
    return lsrCBG

# end of function create_intermediate_lowsimilarity_region


def potentially_contains_lowsimilarity_region():
    pass

# end of function potentially_contains_lowsimilarity_region

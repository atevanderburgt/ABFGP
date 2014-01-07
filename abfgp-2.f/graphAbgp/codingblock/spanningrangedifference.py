"""
Functions concerning the SpanningRangeDifference (SPRDIF) of CodingBlockGraphs
"""

# Module metadata
__authors__ = "Ate van der Burgt"
__license__ = "MIT"


# graphAbgp import
from graphAbgp.exceptions import InproperlyAppliedArgument

# Abgp imports
from abgp_warnings import UnexpectedEventWarning

# Global Variable Import
from settings.codingblockgraph import (
    CBG_MIN_AA_LENGTH,
    CBG_SPRDIF_MIN_NODE_COUNT
    )


def has_left_spanningrange_difference(cbg,**kwargs):
    """
    Does CodingBlockGraph has a spanningrange difference on the left/5p side?

    @attention: see the spanningrange_difference() for arguments documentation

    @type  cbg: CodingBlockGraph
    @param cbg: CodingBlockGraph instance

    @rtype:  Boolean
    @return: True or False
    """
    return has_spanningrange_difference(cbg,side='left',**kwargs)

# end of function has_left_spanningrange_difference


def has_rigth_spanningrange_difference(cbg,**kwargs):
    """
    Does CodingBlockGraph has a spanningrange difference on the rigth/3p side?

    @attention: see the spanningrange_difference() for arguments documentation

    @type  cbg: CodingBlockGraph
    @param cbg: CodingBlockGraph instance

    @rtype:  Boolean
    @return: True or False
    """
    return has_spanningrange_difference(cbg,side='rigth',**kwargs)

# end of function has_rigth_spanningrange_difference


def has_spanningrange_difference(cbg,**kwargs):
    """
    Does CodingBlockGraph has a spanningrange difference on the requested side?

    @attention: see the spanningrange_difference() for arguments documentation

    @type  cbg: CodingBlockGraph
    @param cbg: CodingBlockGraph instance

    @rtype:  Boolean
    @return: True or False
    """
    # check for an existing spanningrange_difference
    if spanningrange_difference(cbg,**kwargs):
        return True
    else:
        return False

# end of function has_spanningrange_difference


def spanningrange_difference(cbg,side=None,
    correct_sprdif_all_nodes=True,
    sprdif_min_node_count=CBG_SPRDIF_MIN_NODE_COUNT,
    sprdif_min_aa_length=CBG_MIN_AA_LENGTH):
    """
    Return the spanningrange difference of this CBG on the requested side

    @type  cbg: CodingBlockGraph
    @param cbg: CodingBlockGraph instance

    @type  side: string
    @param side: one of ['left','rigth','5p','3p']

    @type  sprdif_min_node_count: positive integer 
    @param sprdif_min_node_count: minimum number of nodes that must have
                                  a sprdif reported

    @type  sprdif_min_aa_length: positive integer
    @param sprdif_min_aa_length: minimum AA length of a sprdif to be reported

    @type  correct_sprdif_all_nodes: Boolean
    @param correct_sprdif_all_nodes: correct the sprdif by removing the weakest
                                     node if all nodes are in the sprdif

    @rtype:  dictionary
    @return: dictionary with Organism identifier as keys, lists of
             spanningrange coordinates as values
    """
    # check if side argument is applied correctly
    sides = ['left','rigth','5p','3p']
    if side not in sides:
        message = "argument `side` (%s) not in sides: %s" % (side,sides)
        raise InproperlyAppliedArgument, message

    # get OMSR and MAXSR to define the SPRDIF
    omsr  = cbg.overall_minimal_spanning_range()
    maxsr = cbg.maximal_spanning_range()
    coords = {}
    for node in omsr.keys():
        # Get the difference between the MaximumSpanningRange
        # and the OverallMinimalSpanningRange,
        if side in ['left','5p']:
            # take only the coordinates on the LEFT/5p side.
            if len(omsr[node]) == 0 or len(maxsr[node]) == 0:
                # unexpected, but for some reason OMSR of this node is empty
                # do not raise but call a STDOUT warning for later debugging
                message = "empty OMSR (%s) for node %s in CBG %s" % (
                        len(omsr[node]), node, cbg )
                omsrwarning = UnexpectedEventWarning(message)
                print omsrwarning
                pass
            else:
                leftrange = list( set(maxsr[node]).difference(omsr[node]).difference(range(max(omsr[node]),max(maxsr[node])+1)) )
                if len(leftrange) < sprdif_min_aa_length:
                    # SPRDIF is smaller than sprdif_min_aa_length
                    pass
                else:
                    # SPRDIF of proper size
                    coords[node] = leftrange

        elif side in ['rigth','3p']:
            # take only the coordinates on the RIGHT/3p side.
            if len(omsr[node]) == 0 or len(maxsr[node]) == 0:
                # unexpected, but for some reason OMSR of this node is empty
                # do not raise but call a STDOUT warning for later debugging
                message = "empty OMSR (%s) for node %s in CBG %s" % (
                        len(omsr[node]), node, cbg )
                omsrwarning = UnexpectedEventWarning(message)
                print omsrwarning
                pass
            else:
                rightrange = list( set(maxsr[node]).difference(omsr[node]).difference(range(min(maxsr[node]),min(omsr[node])+1)) )
                if len(rightrange) < sprdif_min_aa_length:
                    # SPRDIF is smaller than sprdif_min_aa_length
                    pass
                else:
                    # SPRDIF of proper size
                    coords[node] = rightrange
        else:
            raise "wrong side!!!"

    # SPRDIF node count is smaller than sprdif_min_node_count
    if len(coords) < sprdif_min_node_count: return {}

    # if ALL nodes are reported to have a SPRDIF and
    # correct_sprdif_all_nodes==True -> remove the shortest reported length
    succesfull = True
    if correct_sprdif_all_nodes and len(coords) == cbg.node_count() and\
    len(coords) >= sprdif_min_node_count+1:
        minlength = min( [ len(coords[k]) for k in coords.keys() ] )
        shortest = []
        for k in coords.keys():
            if len(coords[k]) == minlength:
                shortest.append( k )
        if len(shortest) == 1:
            del( coords[shortest[0]] )
        elif len(shortest) == 2:
            if cbg.weakest_connected_node() in shortest:
                del( coords[ cbg.weakest_connected_node() ] )
            else:
                succesfull = False
        else:
            succesfull = False

    # check if consensus MaximumSpanningRange after detailed analyses
    if not succesfull:
        return {}
    else:
        return coords

# end of function spanningrange_difference

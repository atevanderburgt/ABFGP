"""
Functions for Extension of CodingBlockGraphs (CBGs)
"""

# Module metadata
__authors__ = "Ate van der Burgt"
__license__ = "MIT"

# graphAbgp Imports
from exceptions import InproperlyAppliedArgument
import spanningrangedifference

# Python Imports
from sets import Set

# Global variables
from settings.codingblockgraph import CBG_MIN_AA_LENGTH 


def extend_on_spanningrange_difference(cbg,side=None,
    sprdif_min_aa_length=CBG_MIN_AA_LENGTH,verbose=False):
    """
    Extend the PacbPORFs in this CBG on an existing Spanning Range Difference

    @type  cbg: CodingBlockGraph
    @param cbg: CodingBlockGraph instance

    @type  side: string
    @param side: one of ['left','rigth','5p','3p']

    @type  sprdif_min_aa_length: positive integer
    @param sprdif_min_aa_length: minimum AA length of a sprdif to be reported

    @type  verbose: Boolean
    @param verbose: print intermediate info to STDOUT for debugging purposes

    @rtype:  Boolean
    @return: True when some PacbPORFs (and as a result the CBG) are extended
    """
    # Get Nodes & coords that are represented in the spanningrange difference
    sprdif = spanningrangedifference.spanningrange_difference(cbg,side,
            sprdif_min_aa_length=sprdif_min_aa_length,
            sprdif_min_node_count=cbg.node_count(),
            correct_sprdif_all_nodes=False)

    if not sprdif:
        # nothing to be extended -> done
        return False

    # if this point is reached, this CBG is extendable
    # make Sets of the lists in sprdif values
    for node,vlist in sprdif.iteritems(): sprdif[node] = Set(vlist)

    # if here, start extending the codingblockgraph
    # IMPORTANT!!! DO NOT!!! clear_cache() of codingblock!!
    toberemoved = {}

    # loop over the pacbporfs and gather those that have
    # no/little overlap with the observed sprdif on that side
    for (k,n1,n2),pacbporf in cbg.pacbps.iteritems():
        qRange = pacbporf.alignment_protein_range_query()
        sRange = pacbporf.alignment_protein_range_sbjct()
        overlap = len(sprdif[n1].intersection(qRange)) +\
                  len(sprdif[n2].intersection(sRange))
        if overlap < (sprdif_min_aa_length/5)*2:
            # add this key and pacbporf to dict toberemoved
            toberemoved[(k,n1,n2)] = pacbporf
            ####################################################################
            if verbose:
                print pacbporf, "%s-%s" % ( min(sprdif[n1]), max(sprdif[n1]) ),
                print "%s-%s" % ( min(sprdif[n2]), max(sprdif[n2]) ),
                print len(sprdif[n1].intersection(qRange)),
                print len(sprdif[n2].intersection(sRange))
            ####################################################################

    # if no pacbps are registered for removal -> no extention
    if not toberemoved: return False

    # remove the pacbporfs that lack overlap with the sprdif
    for (k,n1,n2) in toberemoved.keys():
        del( cbg.pacbps[(k,n1,n2)] )
        cbg.del_edge(n1,n2)

    # make new pacbps for the missing edges and update edge weights
    cbg.make_pacbps_for_missing_edges(use_maxsr=True)

    IS_EXTENDED = False 
    for (k,n1,n2) in toberemoved.keys():
        if not cbg.has_edge(n1,n2):
            # clustalw creation of pacbp failed aparently
            # set back this one and continue
            pacbporf = toberemoved[(k,n1,n2)]
            cbg.add_edge(n1,n2,wt=pacbporf.bitscore) 
            cbg.pacbps[(k,n1,n2)] = pacbporf 
            if verbose: print "CLUSTALW FAILED!?"
            # now to the next one
            continue

        # get pacbporf of this cbg 
        pacbporf = cbg.get_pacbps_by_nodes(node1=n1,node2=n2)[0]
        if len(pacbporf) > len(toberemoved[(k,n1,n2)]):
            # set pacbporf.source to clustalw-EXTENDED
            pacbporf.source = 'clustalw-EXTENDED' 
            IS_EXTENDED = True
        else:
            # hmmm... it is even truncated! reset to the original one
            cbg.remove_pacbp(pacbporf,n1,n2)
            cbg.pacbps[(k,n1,n2)] = toberemoved[(k,n1,n2)]
            ####################################################################
            if verbose:
                print "FAILED:  ", pacbporf
                print "RESTORED:", toberemoved[(k,n1,n2)] 
            ####################################################################

    # recreate the cache of the CBG
    cbg.clear_cache()
    cbg.create_cache()
    cbg.update_edge_weights_by_minimal_spanning_range()

    # return status of extension (True)
    return IS_EXTENDED

# end of function extend_on_spanningrange_difference


def extend_on_left_spanningrange_difference(cbg,
    sprdif_min_aa_length=CBG_MIN_AA_LENGTH):
    """
    Extend CBG on left spanningrange difference

    @attention: see extend_on_spanningrange_difference() for documentation
    """
    return extend_on_spanningrange_difference(cbg,side='left',
            sprdif_min_aa_length=sprdif_min_aa_length)

# end of function extend_on_left_spanningrange_difference


def extend_on_rigth_spanningrange_difference(cbg,
    sprdif_min_aa_length=CBG_MIN_AA_LENGTH):
    """
    Extend CBG on rigth spanningrange difference

    @attention: see extend_on_spanningrange_difference() for documentation
    """
    return extend_on_spanningrange_difference(cbg,side='rigth',
            sprdif_min_aa_length=sprdif_min_aa_length)

# end of function extend_on_rigth_spanningrange_difference

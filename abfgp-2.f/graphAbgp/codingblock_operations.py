"""
Varia of operations on CodingBlockGraphs
"""

# Module metadata
__authors__ = "Ate van der Burgt"
__license__ = "MIT"

# graphAbgp Imports

# Abgp imports

# Python Imports

# Global Variable Import


def deepcopy_with_removed_nodes(cbg,removenodes):
    """
    Make a (deep)copy of this a CBG, possibly by omitting some nodes

    @attention: PacbP(ORF)s are not deepcopied!

    @type  cbg: CodingBlockGraph
    @param cbg: CodingBlockGraph that serves as source to (deep)copy from

    @type  removenodes: []
    @param removenodes: list with Node Identifiers to omit in the copied CBG

    @rtype  newcbg: CodingBlockGraph
    @return newcbg: copied CodingBlockGraph
    """
    # import CodingBlockGraph here to prevent circular import
    from graphAbgp import CodingBlockGraph
    newcbg = CodingBlockGraph(
        minimal_overal_spanning_range_size=\
                cbg.MINIMAL_OVERAL_SPANNING_RANGE_SIZE,
        alternative_alignment_overlap_ratio=\
                cbg.ALTERNATIVE_ALIGNMENT_OVERLAP_RATIO
        )
    if removenodes in cbg.get_nodes(): removenodes = [ removenodes ]

    # copy nodes & edges
    for node in cbg.get_nodes():
        if node not in removenodes: newcbg.add_node(node)
    for (node1,node2) in cbg.pairwisecrosscombinations_node():
        if node1 not in removenodes and node2 not in removenodes:
            newcbg.add_edge(node1,node2,cbg.weights[(node1,node2)])

    # get pacbp objects from parental CBG
    newcbg.harvest_pacbps_from_pacbpcollection(cbg)
    newcbg.create_cache()
    # done!
    return newcbg

# end of function deepcopy_with_removed_nodes


def _update_cbg_with_pacbporf_replacements(cbg,replacements,verbose=False):
    """
    Replace the PacbPORFs in the CBG with those in an applied pacbpsdict 

    @type  cbg: CodingBlockGraph
    @param cbg: CodingBlockGraph to replace PacbPORFs in

    @type  replacements: dict
    @param replacements: pacbpsdict <dict data structure> with new PacbPORFs 

    @type  verbose: Boolean
    @param verbose: print intermediate states to STDOUT for debugging purposes
    
    @attention: see lib_crossdata for description of the pacbpsdict
    
    @rtype:  NoneBoolean
    @return: None (no repleacements to do), True (updated) or False (no OMSR!)
    """
    # do the replacements of (corrected) pacbporfs
    for (currentkey,nodeQ,nodeS),newpacbporf in replacements.iteritems():
        newkey = newpacbporf.construct_unique_key(nodeQ,nodeS)
        cbg.set_edge_weight(nodeQ,nodeS,wt=newpacbporf.bitscore)
        ####################################################################
        if verbose:
            print "OLD:", cbg.pacbps[(currentkey,nodeQ,nodeS)], currentkey
            print "NEW:", newpacbporf, newkey, nodeQ, nodeS
        ####################################################################
        del( cbg.pacbps[(currentkey,nodeQ,nodeS)] )
        cbg.pacbps[(newkey,nodeQ,nodeS)] = newpacbporf

    # Check if, after the replacement, still an OMSR is available.
    # In malafide CBGs with very poor identity% and many gaps,
    # this function likely breaks the OMSR criterion, thereby
    # providing us with a clear signal that this CBG must be removed!

    if replacements:
        cbg.clear_cache()
        hasomsr = cbg.has_overall_minimal_spanning_range()
        ####################################################################
        if verbose: print "after replacements OMSR?:", hasomsr
        ####################################################################
        if not hasomsr:
            # no OMSR left; return status False. In the parental
            # function, this is a signal for a NoOverallMinimalSpanningRange
            # Exception. In the parental function even one level higher up,
            # this Exception is the signal for removal of this cbg.
            return False

        # If here, then there is still an OMSR in the CBG
        # The CBG is succesfully optimized!
        cbg.create_cache()
        cbg.update_edge_weights_by_minimal_spanning_range()
        ########################################################################
        if verbose:
            omsr = cbg.overall_minimal_spanning_range()
            for node in omsr.keys(): print node,min(omsr[node]),max(omsr[node])
        ########################################################################
        return True
    else:
        return None

# end of function _update_cbg_with_pacbporf_replacements

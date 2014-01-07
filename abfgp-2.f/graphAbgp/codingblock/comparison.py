"""
Functions for comparing (the contents of) CodingBlockGraphs
"""

# Module metadata
__authors__ = "Ate van der Burgt"
__license__ = "MIT"

# graphPlus Imports
import graphPlus.comparison

# Python Imports

def cbgs_identical_pacbp_analysis(first,second):
    """
    Obtain data about identity between PacbPs in 2 CodingBlockGraphs

    @type  first: CodingBlockGraph
    @param first: CodingBlockGraph object

    @type  second: CodingBlockGraph
    @param second: CodingBlockGraph object

    @rtype:  {}, {}
    @return: nodes_seen, max_nodes_seen

    @attention: nodes_seen counts how often nodes' PacbPs are identical
    @attention: max_nodes_seen counts how often nodes' edges exist in both CBGs
    """
    # gather list of nodes present in both CBGs
    mutual_nodes = comparison.mutual_nodes(first,second)
    # set counter dicts for:
    # nodes_seen    -> nodes with identical PacbPs in both CBGs
    # max_node_seen -> nodes that have edges between them in both CBGs
    nodes_seen    = {}
    max_node_seen = {}
    for node in mutual_nodes:
        nodes_seen[node] = 0
        max_node_seen[node] = 0
    # loop over all node combinations / edges
    for nodeQ,nodeS in first.pairwisecrosscombinations_node():
        # only deal with combinations of mutual nodes
        if nodeQ not in mutual_nodes or nodeS not in mutual_nodes:
            continue
        # get PacbP(ORF) of first (donor) and second (acceptor) CBG
        donorPacbp, accepPacbp = None, None
        if first.has_edge(nodeQ,nodeS):
            donorPacbp = first.get_pacbps_by_nodes(node1=nodeQ,node2=nodeS)[0]
        if second.has_edge(nodeQ,nodeS):
            accepPacbp = second.get_pacbps_by_nodes(node1=nodeQ,node2=nodeS)[0]
        # check if both edges exist -> increase max_node_seen
        if accepPacbp and donorPacbp:
            max_node_seen[nodeQ]+=1
            max_node_seen[nodeS]+=1
            # check if Pacbps are identical -> increase nodes_seen
            if donorPacbp.barcode() == accepPacbp.barcode():
                nodes_seen[nodeQ]+=1
                nodes_seen[nodeS]+=1

    # return both dictionaries
    return nodes_seen, max_node_seen

# end of function cbgs_identical_pacbp_analysis
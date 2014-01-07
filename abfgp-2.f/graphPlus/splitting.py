"""
Splitting algorithms for python-graph.
"""

# Module metadata
__authors__ = "Ate van der Burgt"
__license__ = "MIT"


# Imports
from copy import deepcopy
from sets import Set

# Split into non-connected sub-graphs

def get_subgraph_count(graph):
    """
    """
    return len( _nodes_per_subgraph(graph) )

# end of function get_subgraph_count


def split_into_subgraphs(graph):
    """
	Split Graph in non-connected subgraphs.

	@attention: Graph splitting is meaningful only for XXX graphs.

	@type  graph: graph
	@param graph: Graph.

   	@rtype:  list
	@return: List with non-connected subgraphs. 
    """
    subgraph_nodes = _nodes_per_subgraph(graph)

    if len(subgraph_nodes) == 1: return [ graph ]

    # if we reach this point, there are subgraphs
    subgraphs = []
    for nodes in subgraph_nodes:
        sub = deepcopy(graph)
        for each in sub.get_nodes():
            if (not each in nodes):
                sub.del_node(each)

        # append to list of subgraphs
        subgraphs.append( sub )

    # and return list of subgraphs
    return subgraphs

# end of function split_into_subgraphs


def _nodes_per_subgraph(graph):
    """
	@type  graph: graph
	@param graph: Graph.

	@rtype:  list
	@return: List of lists of connected nodes, that are not interconnected (subgraphs)
    """
    accessib = graph.accessibility().values()	# Accessibility matrix values
    if not accessib:
        return []
    else:
        subgraphs = [ Set(accessib[0]) ]
        for acces in accessib[1:]:
            for subgraph in subgraphs:
                if subgraph.intersection(acces):
                    subgraph.update(acces)
                    break
            else:
                subgraphs.append( Set(acces) )
        return [ list(sg) for sg in subgraphs ]

# end of function _nodes_per_subgraph

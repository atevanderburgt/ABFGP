"""
Merging algorithms for python-graph.
"""

# Module metadata
__authors__ = "Ate van der Burgt"
__license__ = "MIT"


# Imports
from copy import deepcopy

def connect_graphs(gra1,gra2,gra1node,gra2node):
    """
    """
    new = deepcopy(gra1)
    cnt = gra1.node_count()
    n2id = {}
    for n in gra2.get_nodes():
        id = new.add_node()
        n2id[n] = id
        new.set_node_label( id, gra2.get_node_label(n) )
    # now copy the weights
    for (n1,n2), weight in self.weights.iteritems():
        id1, id2 = n2id[n1], n2id[n2]
        new.weights[(id1, id2)] = weight
    # and now actually merge the graphs on the nodes
    return merge_nodes( new, gra1node, n2id[gra2node] )

# end of function connect_graphs

def merge_nodes(gra,u,v):
    """
    Merge 2 existing nodes of a graph.
    
    @type  u: number of the node to be merged
	@param u: One node.

    @type  v: number of the node to be merged
	@param v: One node.

    @rtype:   graph
	@return:  Graph with merged node.
    """
    # first check if both nodes exist
    if (u in gra.get_nodes()) and (v in gra.get_nodes()):
        new = deepcopy(gra)
        # create new edges and copy weights
        for n in new.get_node(v):
            if n == u: continue
            new.add_edge(u,n)
            new.set_edge_weight(u, n, new.weights[(n, v)] )
        # delete the edge between these nodes (if there)
    	new.del_edge(u, v)
        # delete the node that is now merged
        new.del_node(v)
        # and return
        return new
    else:
        # not both nodes present; just return the graph!
        return gra

# end of function merge_nodes

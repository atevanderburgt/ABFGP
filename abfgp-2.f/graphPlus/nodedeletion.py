"""
Node deletion algorithms for python-graph.
"""

# Module metadata
__authors__ = "Ate van der Burgt"
__license__ = "MIT"


def del_node(gra, node):
    """
    Remove a node from the graph.

    @attention: this function removes all concerning edges, arrows and weigths.

    @type  node: *
    @param node: One Node Identifier
    """
    if node in gra.get_nodes():
        # delete all edges, weights and outward
        # pointing arrows from this node
        while gra.get_node(node):
            u = gra.get_node(node)[0]
            if gra.has_edge(node,u):
                gra.del_edge(node,u)
            if gra.has_arrow(node,u):
                gra.del_arrow(node,u)
        # scan all nodes to see if there is an
        # inwards pointing arrow
        for each in gra.get_nodes():
            if node in gra.get_node(each):
                gra.del_arrow(each,node)
        # delete the node in the gra
        del( gra.nodes[node] )
        # scan all gra.weights for still occurring references
        edges = gra.weights.keys()
        for u,v in edges:
            if u == node or v == node:
                del( gra.weights[(u,v)] ) 

# end of function del_node


def minimal_node_degree(gra,degree):
    """
    Remove all nodes with a degree less than requested

    @type  gra: graph
    @param gra: Graph object

    @type  degree: integer
    @param degree: integer value (gte 3)
    """
    # some input sanity checks
    if type(degree) != type(int()):
        raise InproperlyAppliedArgument, "variable 'degree' not integer but '%s'" % type(degree)

    # gather nodes with a to low degree
    delnodes = []
    for node in gra.get_nodes():
        if gra.degree(node) < degree:
            delnodes.append(node)

    # delete all the nodes with a to low degree
    for node in delnodes: gra.del_node(node)

# end of function minimal_node_degree


#def remove_low_connectivity_nodes(gra,min_connectivity=2):
#    """
#    Remove all nodes below a certain connectivity threshold
#
#    @type  min_connectivity: integer
#    @param min_connectivity: minimal connectivity to maintain a node
#    """
#    removing = True
#    while removing:
#        # now remove all with to low connectivity
#        for node in gra.get_nodes():
#            if len(gra.nodes[node]) < min_connectivity:
#                gra.del_node(node)
#                # and break out of the node for loop
#                break
#        else:
#            removing = False
#    
## end of remove_low_connectivity_nodes


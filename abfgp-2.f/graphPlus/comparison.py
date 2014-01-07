"""
Comparison algorithms for python-graph.
"""

# Module metadata
__authors__ = "Ate van der Burgt"
__license__ = "MIT"


# Imports
from sets import Set
from copy import deepcopy

def has_identical_topology(gra1,gra2):
    """
    Do both graphs have exactly the same nodes and edges?

	@type  gra1: graph
	@param gra1: Graph.

	@type  gra2: graph
	@param gra2: Graph.

   	@rtype:  Boolean
	@return: True or False
    """
    if has_identical_nodes(gra1,gra2) and gra1.edge_count() == gra2.edge_count():
        # identical nodes and same number of edges
        # now it is worth doing a more elaborate comparison
        if is_subgraph_of(gra1,gra2):
            return True
        else:
            return False

    else:
        return False

# end of function has_identical_topology


def has_identical_nodes(gra1,gra2):
    """
    Do both graphs have exactly the same nodes?

	@type  gra1: graph
	@param gra1: Graph.

	@type  gra2: graph
	@param gra2: Graph.

   	@rtype:  Boolean
	@return: True or False
    """
    if gra1.node_count() != gra2.node_count():
        return False
    else:
        if len(common_nodes(gra1,gra2)) == gra1.node_count():
            return True
        else:
            return False

# end of function has_identical_nodes


def common_nodes(gra1,gra2):
    """
	@type  gra1: graph
	@param gra1: Graph.

	@type  gra2: graph
	@param gra2: Graph.

   	@rtype:  list
	@return: List with nodes present in both graphs 
    """
    return list( Set(gra1.get_nodes()).intersection(gra2.get_nodes()) )

# end of function common_nodes


def mutual_nodes(gra1,gra2):
    """
    Return a list of nodes shared by both graphs
    """
    return common_nodes(gra1,gra2)

# end of function mutual_nodes


def is_subgraph_of(gra1,gra2):
    """
    Is gra1 fully covered in gra2 (means all nodes and all edges present?

	@type  gra1: graph
	@param gra1: Graph.

	@type  gra2: graph
	@param gra2: Graph.

   	@rtype:  Boolean
	@return: True or False
    """
    if not Set( common_nodes(gra1,gra2) ) == Set( gra1.get_nodes() ):
        # not all nodes of gra1 are present in the
        # potential `parental` graph gra2
        return False
    for (n1,n2) in gra1.weights.keys():
        if not gra2.has_edge(n1,n2):
            # edge of gra1 is not present in the
            # potential `parental` graph gra2
            return False
    else:
        # all edges of gra1 are present in the
        # potential `parental` graph gra2
        return True

# end of function is_subgraph_of


def identical_graph_weight_difference(gra1,gra2,ignore_zero_edges=False):
    """
    Relative weight difference of topologically identical graphs.

	@type  gra1: graph
	@param gra1: Graph.

	@type  gra2: graph
	@param gra2: Graph.

   	@rtype:  Positive float or False
	@return: Positive float or False if graphs are not topologically identical

    @attention: Function is only meaningfull for graphs without arrows.
    @attention: Only applicable when both graphs have identical topology.
    """
    # This function kind of 'aligns' the weights of two (near) identical graphs
    # and returns a positive float of 0.0 to 1.0 for their difference, with 0.0
    # meaning 100% identical and 1.0 meaninging completely different.
    #
    # An example.
    # Two identical and complete graphs g1,g2 have nodes A,B,C
    # So, both graphs have (symmetrical) edges ab,bc,ca.
    # Graph g1 has edge weights 8,11,1, graph g2 has weights 60,22,18.
    # The relative weights of graph g1 are 0.40, 0.55, 0.05
    # The relative weights of graph g2 are 0.60, 0.22, 0.18
    # The identical_graph_weight_difference score measures:
    # abs( 0.40 - 0.60 ) + abs( 0.55 - 0.22 ) + abs( 0.05 - 0.18 )
    # 0.20 + 0.33 + 0.13 = 0.68
    # Divided by two     = 0.34
    # The final division is to limit the outcome in the range 0.0-1.0
    # Because, imagine the same graphs with weights 1,0,0 and 0,70,30.
    # This would result in absolute relative difference of 1.0, 0.7 and 0.3,
    # summing up to a total of 2.0. 

    # First, check if both graphs have identical topology.
    if not has_identical_topology(gra1,gra2):
        return False
    # compare the dictionaries of weights
    gra1tw = float( gra1.total_weight() )
    gra2tw = float( gra2.total_weight() )
    diffs = []
    seen_edges = []
    for edge, weigth in gra1.weights.iteritems():
        if ignore_zero_edges and weigth == 0.0: continue
        if edge in seen_edges: continue
        (node1,node2) = edge
        seen_edges.append(edge)
        seen_edges.append((node2,node1))
        if ignore_zero_edges and gra2.weights[edge] == 0.0: continue
        # calculate relative weights; check for tw==0.0
        if gra1tw == 0.0: ratio1 = 0.0
        else:             ratio1 = float(weigth)/gra1tw
        if gra2tw == 0.0: ratio2 = 0.0
        else:             ratio2 = float(gra2.weights[edge])/gra2tw 
        # calculate difference between relative weights
        relative_difference = ratio1 - ratio2
        # append to list of absolute differences
        diffs.append( abs( relative_difference ) )
    # and return the sum of all differences, divided by two
    return sum(diffs) / 2.0

# end of function identical_graph_weight_difference


def identical_graph_absolute_weight_difference(gra1,gra2,ignore_zero_edges=False):
    """
    Absolute weight difference of topologically identical graphs.

    @type  gra1: graph
    @param gra1: Graph.

    @type  gra2: graph
    @param gra2: Graph.

    @rtype:  Positive float or False
    @return: Positive float or False if graphs are not topologically identical

    @attention: Function is only meaningfull for graphs without arrows.
    @attention: Only applicable when both graphs have identical topology.
    @attention: for detailed explanation, see `identical_graph_weight_difference`
    """
    # First, check if both graphs have identical topology.
    if not has_identical_topology(gra1,gra2):
        return False
    # compare the dictionaries of weights
    diffs = []
    seen_edges = []
    for edge, weigth in gra1.weights.iteritems():
        if ignore_zero_edges and weigth == 0.0: continue
        if edge in seen_edges: continue
        (node1,node2) = edge
        seen_edges.append(edge)
        seen_edges.append((node2,node1))
        if ignore_zero_edges and gra2.weights[edge] == 0.0: continue
        # calculate difference between weights
        diffs.append( float(weigth) - float(gra2.weights[edge]) )
    # and return the sum of all differences, divided by two
    return sum(diffs) / 2.0

# end of function identical_graph_absolute_weight_difference


def deepcopywithrelativeweights(gra):
    """
    Return as a graph with relative weights, means sum(weights)==1.0 

	@type  gra: graph
	@param gra: Graph.

	@rtype:  Graph
	@return: Graph.
    """
    gratw = float( gra.total_weight() )
    _gra = deepcopy(gra)
    for edge, weigth in _gra.weights.iteritems():
        _gra.weights[edge] = float(weigth)/gratw
    return _gra

# end of function deepcopywithrelativeweights


def phylogeneticaltree(gra):
    """
    Return a phylogenetical tree in bracket-string-represention of an organism graph

	@type  gra: graph
	@param gra: Graph.

   	@rtype:  String
	@return: phylogenetical tree string
    """
    _gra = deepcopy(gra)
    while _gra.node_count() > 2:
        (n1,n2) = _gra.strongest_edge()
        #_gra.del_edge(n1,n2)
        combi_node = (n1,n2)
        # add new combination node
        _gra.add_node( combi_node )
        for node in _gra.get_nodes():
            if node == (n1,n2): continue
            if node not in (n1,n2):
                combi_wt = float( _gra.get_edge_weight(node,n1) + _gra.get_edge_weight(node,n2) ) / 2.0
                _gra.add_edge(node,combi_node,wt=combi_wt)
                _gra.del_edge(node,n1)
                _gra.del_edge(node,n2)
        # delete the old nodes
        _gra.del_node(n1)
        _gra.del_node(n2)

    # return the phylogetetical string representation
    return tuple( _gra.get_nodes() )

# end of function phylogeneticaltree

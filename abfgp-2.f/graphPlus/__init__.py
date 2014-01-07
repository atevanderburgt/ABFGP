# Copyright (c) 2008 Ate van der Burgt <ate.vanderburgt@wur.nl>
#
# Permission is hereby granted, free of charge, to any person
# obtaining a copy of this software and associated documentation
# files (the "Software"), to deal in the Software without
# restriction, including without limitation the rights to use,
# copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the
# Software is furnished to do so, subject to the following
# conditions:

# The above copyright notice and this permission notice shall be
# included in all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
# OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
# HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
# WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
# OTHER DEALINGS IN THE SOFTWARE.


"""
python-graph-Plus

A library for working with graphs in Python.
A Plus-version of the python-graph package of Pedro Matiello
"""


# Module metadata
__authors__ = "Ate van der Burgt"
__license__ = "MIT"


# Imports
import graph
import splitting
import comparison
import merging
import completegraph
import nodedeletion
from copy import deepcopy


class graphPlus(graph.graph):
    """
    GraphPlus class.

    Most Basic operations are defined in the Graph class
    Additional Basic operations are defined in this file.
    Algorithms should refer to external files.
    """

    def __init__(self):
        """
        Initialize a graph.
        """
        graph.graph.__init__(self)

    # end of function __init__


    def get_ordered_nodes(self):
        """
        Return ordered iterative list of nodes

       	@rtype:  list
    	@return: ordered list of Node identifiers
        """
        nodes = self.get_nodes()
        nodes.sort()
        return nodes

    # end of function


    def has_node(self,node):
        """
        Is this node identifier occurring in this graph?

        @type  node: *
        @param node: One Node Identifier

        @rtype:  Boolean 
        @return: True or False 
        """
        return node in self.get_nodes()

    # end of function has_node 


    def deepcopyasunweighted(self):
        """
        Returns an unweighted (weight=1) graph
    
       	@rtype:  graph
    	@return: deepcopied, unweighted version of the input graph
        """
        _gra = deepcopy(self)
        for edge,wt in _gra.weights.iteritems():
            _gra.weights[edge] = 1.0
        return _gra
    
    # end of function deepcopyasunweighted


    def get_node_weights(self,node):
        """
        Returns list of weights of edges connected to a node

        @type  node: *
        @param node: One Node Identifier

        @rtype:  list
        @return: List of weights of outgoing edges of this node
        """
        return [ self.weights[(u,node)] for u in self.get_node(node) ]

    # end of function get_node_weights


    def get_node_weighted_connectivity(self,node):
        """
        Returns relative contribution of this node to all the edges' weight in the graph

        @type  node: *
        @param node: One Node Identifier

        @rtype:  float
        @return: relative contribution of this node to the graph 
        """
        if not self.weights: return 0.0
        # The final division by 2.0 enables the relative contribution of this node
        # The sum of the outcome of this function for all nodes in the graph will be 1.0
        if self.total_weight() == 0.0:
            # avoid the ZeroDivisionError
            return 0.0
        else:
            return ( float(sum(self.get_node_weights(node)))/self.total_weight() ) / 2.0

    # end of function get_node_weighted_connectivity


    def get_node_weighted_connectivity_observed_vs_expected(self,node):
        """
        Returns ratio of relative contribution of this node to all the edges' weight in the graph

        @type  node: *
        @param node: One Node Identifier

        @rtype:  float
        @return: relative contribution of this node to the graph, as a ratio of expected contribution
        """
        if not self.weights: return 0.0
        # The final division by 2.0 enables the relative contribution of this node
        # The sum of the outcome of this function for all nodes in the graph will be 1.0
        return self.get_node_weighted_connectivity(node) * self.node_count() 

    # end of function get_node_weighted_connectivity_observed_vs_expected


    def lowest_node_degree(self):
        """
        Return lowest node degree

        @rtype:  integer
        @return: degree of node with least outgoing edges
        """
        return min([ self.degree(node) for node in self.get_nodes() ])

    # end of function lowest_node_degree


    def highest_node_degree(self):
        """
        Return highest node degree

        @rtype:  integer
        @return: degree of node with most outgoing edges
        """
        return max([ self.degree(node) for node in self.get_nodes() ])

    # end of function highest_node_degree


    def weakest_connected_node(self):
        """
        Return (relatively) weakest connected node

        @attention: only usefull if weights != 1 are assigned

        @rtype:  *
        @return: Node Identifier
        """
        if not self.get_nodes(): return None
        tmp = []
        for node in self.get_nodes():
            tmp.append( ( self.get_node_weighted_connectivity(node), node ) )
        tmp.sort()
        return tmp[0][1]

    # end of function weakest_connected_node


    def strongest_connected_node(self):
        """
        Return (relatively) strongest connected node

        @attention: only usefull if weights != 1 are assigned

        @rtype:  *
        @return: Node Identifier
        """
        if not self.get_nodes(): return None
        tmp = []
        for node in self.get_nodes():
            tmp.append( ( self.get_node_weighted_connectivity(node), node ) )
        tmp.sort()
        tmp.reverse()
        return tmp[0][1]

    # end of function strongest_connected_node


    def weakest_edge(self):
        """
        Return weakest edge

        @attention: only usefull if weights are assigned

        @rtype:  tuple
        @return: An edge ( NodeIdentifierA, NodeIdentifierB )
        """

        if self.weights:
            lowest = min(self.weights.values())
            for edge,wt in self.weights.iteritems():
                if wt == lowest:
                    return edge
        else:
            return None

    # end of function weakest_edge


    def strongest_edge(self):
        """
        Return strongest edge

        @attention: only usefull if weights are assigned

        @rtype:  tuple
        @return: An edge ( NodeIdentifierA, NodeIdentifierB )
        """

        if self.weights:
            highest = max(self.weights.values())
            for edge,wt in self.weights.iteritems():
                if wt == highest:
                    return edge
        else:
            return None

    # end of function strongest_edge


    def total_weight(self):
        """
        Return total weight represented by all edges

        @attention: only 100% correct when no arrows are present

        @rtype:  number
        @return: Total weight
        """
        return sum(self.weights.values())/2

    # end of function total_weight


    def average_weight(self):
        """
        Return total weight devided by number of edges

        @attention: only 100% correct when no arrows are present

        @rtype:  number
        @return: Average weight per edge
        """
        if self.edge_count() == 0:
            return 0.0
        else:
            return float( self.total_weight() ) / float( self.edge_count() )

    # end of function average_weight


    def node_connectivity_counts(self):
        """
        @rtype:  list
        @return: list with counts how many nodes have how many edges (0-based)
        """
        ll = [0]*len(self)
        for node in self.get_nodes():
            cnt = len(self.get_node(node))
            ll[cnt]+=1
        return ll

    # end of function node_connectivity_counts


    def add_node(self, node):
        """
        Create a single node.

        @type  node: *
        @param node: One Node identifier
        """
        self.nodes[node] = []

    # end of function add_node


    def remove_edgeless_nodes(self):
        """
        Remove all nodes without edges

        @rtype:  integer
        @return: number of removed nodes 
        """

        del_nodes = []
        for node in self.nodes.keys():
            if not self.nodes[node]:
                del_nodes.append(node)
        for node in del_nodes: self.del_node(node)
        return len(del_nodes)

    # end of function remove_edgeless_nodes


    def remove_low_connectivity_nodes(self,min_connectivity=2):
        """
        Remove all nodes below a certain connectivity threshold

        @type  min_connectivity: integer
        @param min_connectivity: minimal connectivity to maintain a node
        """
        removing = True
        while removing:
            # now remove all with to low connectivity
            for node in self.get_nodes():
                if len(self.nodes[node]) < min_connectivity:
                    self.del_node(node)
                    # and break out of the node for loop
                    break
            else:
                removing = False

    # end of remove_low_connectivity_nodes


    ########################################################
    #### Functions for node deletion                    ####
    ########################################################

    def del_node(self, node):
        """
        Remove a node from the graph.

        @attention: for documentation see nodedeletion.del_node()
        """
        nodedeletion.del_node(self,node)

    # end of function del_node

    def minimal_node_degree(self,degree):
        """
        Remove all nodes with a degree less than requested

        @attention: for documentation see nodedeletion.minimal_node_degree()
        """
        nodedeletion.minimal_node_degree(self,degree)

    # end of function minimal_node_degree 

    def remove_low_connectivity_nodes(self,min_connectivity=2):
        """
        Remove all nodes below a certain connectivity threshold

        @attention: old alias for nodedeletion.minimal_node_degree()
        @attention: for documentation see nodedeletion.minimal_node_degree()
        """
        nodedeletion.minimal_node_degree(self,min_connectivity)

    # end of function remove_low_connectivity_nodes

    ########################################################
    #### Functions for complete graphs                  ####
    ########################################################

    def is_complete(self):
        """
        Is this graph a complete K(s) graph ( fully connected )? 

        @attention: for documentation see completegraph.is_complete()
        """
        completegraph.is_complete(self)

    # end of function is_complete 

    def is_fully_connected(self):
        """
        Is this graph a complete K(s) graph ( fully connected )?

        @attention: alias function for completegraphs.is_complete()
        @attention: for documentation see completegraph.is_complete()
        """
        completegraph.is_complete(self)

    # end of function is_fully_connected 

    def makecompletegraph(self, wt=0.0):
        """
        Make the graph complete, i.e. set all non-existing edges to wt(=0.0) 

        @attention: only 100% correct when no arrows are present
        @attention: for documentation see completegraph.makecompletegraph()
        """
        completegraph.makecompletegraph(self,wt=wt)

    # end of function makecompletegraph 

    def remove_incomplete_components(self, s=None):
        """
        Remove all nodes and edges that are not part of a Complete K(s) graph 

        @attention: for documentation see completegraph.remove_incomplete_components()
        @attention: function is overwritten in graphAbgp.OrganismGraph class (a coloured graph)
        """
        completegraph.remove_incomplete_components(self,s=s)

    # end of function remove_incomplete_components

    ########################################################
    #### Functions for ....                             ####
    ########################################################


    def has_edge(self, u, v):
        """
        Does edge (u,v) connecting nodes u and v exist?

        @type  u: *
        @param u: One node Identifier.

        @type  v: *
        @param v: Other node Identifier.

        @rtype:  True or False
        @return: True or False
        """
        if u in self.get_nodes() and v in self.get_nodes():
            if (v in self.nodes[u]) and (u in self.nodes[v]):
                if self.weights.has_key((u, v)) and self.weights.has_key((v, u)):
                    return True
                else:
                    return False
            else:
                return False
        else:
            return False

    # end of function has_edge


    def has_arrow(self, u, v):
        """
        Does arrow (u,v) pointing from node u to node v exist?

        @type  u: *
        @param u: One node Identifier

        @type  v: *
        @param v: Other node Identifier.

        @rtype:  True or False
        @return: True or False
        """
        if (v in self.nodes[u]) and self.weights.has_key((u, v)):
            return True
        else:
            return False


    def set_edge_weight(self, u, v, wt):
        """
        Set the weight of an edge.

        @type  u: *
        @param u: One node Identifier.

        @type  v: *
        @param v: Other node Identifier.

        @type  wt: number (integer or float)
        @param wt: Weight of this edge.
        """
        self.weights[(u, v)] = wt
        self.weights[(v, u)] = wt

    # end of function set_edge_weight


    def get_subgraph_count(self):
        """
        Get number of non-connected subgraphs.
        """
        return splitting.get_subgraph_count(self)

    # end of function get_subgraph_count


    def node_count(self):
        """ """
        return len(self)

    # end of function node_count


    def edge_count(self):
        """
        Get number of edges in the graph.

        @attention: Function is only meaningfull for graphs without arrows.

        @rtype:  number
        @return: Number of edges in the graph
        """
        return len(self.weights)/2

    # end of function edge_count


    def max_edge_count(self):
        """
        Get maximal number of edges possible.

        @rtype:  number
        @return: Maximum number of possible edges
        """
        return sum(range(0,len(self)))

    # end of function max_edge_count


    def connectivitysaturation(self):
        """
        Get relative connectivity saturation.

        @rtype:  float
        @return: ratio of maximum number of possible edges
        """
        if not self.max_edge_count(): return 1.0
        # divide count by maximum number
        return float( self.edge_count() ) / float( self.max_edge_count() )

    # end of function connectivitysaturation


    def degree(self,node):
        """
        Return the degree of a node.

        @type  node: *
        @param node: One node Identifier.

        @rtype:  number
        @return: Number of edges connected to this node.
        """
        return len(self.get_node(node))

    # end of function degree

# end of class graphPlus

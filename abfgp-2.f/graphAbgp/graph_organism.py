################################################################################
### OrganismGraph class                                                     ####
################################################################################

# graphAbgp imports
import graphPlus
import pacb.recombination
from exceptions import *

# Python imports
from sets import Set

class OrganismGraph(graphPlus.graphPlus):
    """
    OrganismGraph (OSG) class, inheriting from graphPlus class.
    """

    def __init__(self):
        """
		Initialize a OrganismGraph
        """
        graphPlus.graphPlus.__init__(self)
        # as self.weights, _edge_binary_entropies is a dictionary
        # that stores a score, in this case the binary entropy of a node.
        self._edge_binary_entropies = {}

    # end of function __init__


    def __str__(self):
        """
        Basal string representation of OrganismGraph
        """
        return "<%s wt=%s [%s]>" % (
            self.__class__.__name__,
            self.total_weight(),
            " ".join([ str(node) for node in self.get_ordered_nodes() ]),
            )

    # end of function __str__


    def _organism_from_node(self,node):
        """
        Get the organism identifier from the node identifier.

        @attention: MUST BE OVERRIDDEN IN THE SUBCLASSES, if the returned value is not correct

        @attention: `_organism_from_node` and `organism_by_node` are aliasses

        @type  node: *
        @param node: One Node

        @rtype:  *
        @return: organism identifier (presumably a string)
        """
        return node

    # end of function _organism_from_node 


    def organism_by_node(self,node):
        """
        Get the organism identifier from the node identifier.

        @attention: MUST BE OVERRIDDEN IN THE SUBCLASSES, if the returned value is not correct

        @type  node: *
        @param node: One Node

        @rtype:  *
        @return: organism identifier (presumably a string)
        """
        return node

    # end of function organism_by_node


    def organism_set(self):
        """
        Return of Set of unique organisms/genes present among the nodes

        @rtype:  sets.Set
        @return: Set with unique organisms/genes present among the nodes
        """
        return Set([ self._organism_from_node(node) for node in self.get_nodes() ])

    # end of function organism_set


    def get_organism_edges(self,organism1,organism2):
        """
        Get all edges between these two organisms

        @type  organism1: *
        @param organism1: organism identifier (presumably a string)

        @type  organism2: *
        @param organism2: organism identifier (presumably a string)

        @rtype:  list 
        @return: list with all edges between these two organism identifiers
        """
        selected_edges = []
        for node1 in self.get_organism_nodes(organism1):
            for node2 in self.get_node(node1):
                if self.organism_by_node(node2) == organism2:
                    selected_edges.append( ( node1, node2 ) )
        return selected_edges
        
    # end of function get_organism_edges



    def mutual_nodes(self,other):
        """ """
        return self.node_set().intersection(other.node_set())
    # end of function mutual_nodes


    def different_nodes(self,other):
        """ """
        return self.node_set().difference(other.node_set())
    # end of function different_nodes


    def mutual_organisms(self,other):
        """ """
        return self.organism_set().intersection(other.organism_set())
    # end of function mutual_organisms


    def different_organisms(self,other):
        """ """
        return self.organism_set().difference(other.organism_set())
    # end of function different_organisms


    def organism_set_size(self):
        """
        Return of size of Set of unique organisms/genes present among the nodes

		@rtype:  integer
		@return: number of unique organisms/genes present among the nodes
        """
        return len(self.organism_set())

    # end of function organism_set_size


    def node_set(self):
        """
        Return the nodes of the graph as a Set i.s.o. a list

        @rtype:  sets.Set
        @return: Set with nodes occuring in the graph
        """
        return Set( self.get_nodes() )

    # end of function node_set


    def node_set_size(self):
        """
        Return the number of nodes in graph

        @rtype:  integer
        @return: number of nodes

        @attention: alias for node_count(), added for uniform naming with organism__set_(size)
        """
        return self.node_count()

    # end of function node_set_size


    def organism_edge_set(self):
        """
        Return a Set of ORDERED tuples of oganism combinations of existing edges
        """
        orgedgeset = Set()
        for node1,node2 in self.weights.keys():
            combi = [node1,node2]
            combi.sort()
            if combi != [node1,node2]: continue
            org1 = self.organism_by_node(node1)
            org2 = self.organism_by_node(node2)
            orgedgeset.add( ( org1, org2 ) )
        # return the created  orgedgeset
        return orgedgeset

    # end of function organism_edge_set


    def crosscombinationcount(self):
        """
        Get the number of unique cross combinations in the graph

        @rtype:  integer
        @return: number of unique combinations (of organisms)
        """
        return sum(range(0,self.organism_set_size()))

    # end of function crosscombinationcount



    def pairwisecrosscombinations_organism(self):
        """
        Get unique cross combinations of organisms in the graph

        @rtype:  list
        @return: list of all unique combinations of two organisms

        @attention: returns ORDERED list of ORDERED tuples of organism combinations
        """
        return pacb.recombination.pairwise( list(self.organism_set()) )

    # end of function pairwisecrosscombinations_organism


    def pairwisecrosscombinations_node(self):
        """
        Get unique cross combinations of nodes in the graph

        @rtype:  list
        @return: list of all unique combinations of two nodes

        @attention: returns ORDERED list of ORDERED tuples of node combinations
        """
        return pacb.recombination.pairwise( self.get_nodes() )

    # end of function pairwisecrosscombinations_node


    def organism_counts(self):
        """
        Count the number of nodes (values) per organism (keys), returned in a dictionary

        @rtype:  dictionairy
        @return: number of nodes (values) per organism (keys) in a dictionary
        """
        d = {}
        for org in self.organism_set():
            d[org] = [ self._organism_from_node(node) for node  in self.get_nodes() ].count(org)
        # and return
        return d

    # end of function organism_set


    def remove_edgeless_nodes(self,organism=None):
        """
        Remove all nodes without edges

        @type  organism: *
        @param organism: organism identifier (presumably a string)

        @rtype:  integer
        @return: number of removed nodes 
        """
        if not organism:
            return graphPlus.graphPlus.remove_edgeless_nodes(self) 
        else:
            del_nodes = []
            for node in self.get_organism_nodes(organism):
                if not self.nodes[node]:
                    del_nodes.append(node)
            for node in del_nodes: self.del_node(node)
            return len(del_nodes)

    # end of function remove_edgeless_nodes


    def remove_organism_nodes(self,organism):
        """
        Remove all the nodes of this organism from this graph

        @type  organism: *
        @param organism: organism identifier (presumably a string)
        """
        if organism not in self.organism_set():
            pass
        else:
            for node in self.get_organism_nodes(organism):
                self.del_node(node)

    # end of function remove_organism_nodes


    def remove_incomplete_components(self,s=None):
        """
        Remove all nodes and edges that are not part of a Complete K(s) graph 
    
        @type  self: graph
        @param self: graph object
    
        @type  s: integer
        @param s: integer value (gte 3)
    
        @attention: OVERWRITES graphPlus.completegraph.remove_incomplete_components (OrganismGraph isa coloured graph) 
        """
        # by default, s is equal to the number of organisms/genes
        if not s: s = self.organism_set_size()
        
        # first, perform the basal remove_incomplete_components 
        graphPlus.graphPlus.remove_incomplete_components(self,s=s)

        # second, remove nodes that do not have links to all other organisms
        delnodes = []
        for node in self.get_nodes():
            if len(Set([ self._organism_from_node(n) for n in self.get_node(node)])) < s-1:
                delnodes.append(node)
        # remove all the failing nodes
        for node in delnodes: self.del_node(node)

    # end of function remove_incomplete_components


    def total_binary_entropy(self):
        """
        Total binary entropy of all nodes in this graph.
        For explanation of binary entropy, see ... TODO

        @rtype:  float
        @return: total binary entropy of all the nodes in this graph
        """
        tot = 0.0
        for v1,v2 in self._edge_binary_entropies.values():
            tot = tot + v1 + v2
        # devide by 4; each edge is present twice and has 2 values!
        return tot / 4.0

    # end of function total_binary_entropy


    def average_binary_entropy(self):
        """
        Total binary entropy averaged by number of edges

        @rtype:  float
        @return: average binary entropy of all edges in the graph (0.0..1.0)
        """
        if self.node_count() == 0:
            return 0.0
        elif self.edge_count() == 0:
            return 0.0
        else:
            return self.total_binary_entropy() / float(self.edge_count())

    # end of function average_binary_entropy


    def _update_edge_binary_entropies(self):
        """
        Remove edges from self._edge_binary_entropies that are no longer
        present in self.weights (after deletion of nodes/edges)

        @attention: USE THIS FUNCTION AFTER DELETION OF NODES FROM THE GRAPH
        """
        keys_weights = self.weights.keys()
        keys_binary  = self._edge_binary_entropies.keys()
        # check if binary key is present in weights
        for kb in keys_binary:
            if kb not in keys_weights:
                del( self._edge_binary_entropies[kb] )

    # end of function _update_edge_binary_entropies


    def remove_non_uniformly_connected_nodes(self):
        """
        Removes all nodes that are not connected to all organisms,
        except of course the organism from this node it self.
        """
        del_nodes = []
        for node in self.get_nodes():
            if not self.is_node_connected_to_all_organisms(node):
                del_nodes.append(node)
        for node in del_nodes: self.del_node(node)

    # end of function remove_non_uniformly_connected_nodes


    def is_node_connected_to_all_organisms(self,node):
        """
        Is this node connected to at least one node of every organism?

        @type  node: *
        @param node: One Node

        @rtype:  Boolean
        @return: True or False
        """
        connectedto = len(Set([ self._organism_from_node(n) for n in self.get_node(node) ]))
        if connectedto == self.organism_set_size() - 1:
            return True
        else:
            return False

    # end of function is_node_connected_to_all_organisms


    def get_organism_nodes(self,organism):
        """
        Get all the nodes of a specific organism from the graph

        @type  organism: *
        @param organism: organism identifier (presumably a string)

        @rtype:  list
        @return: list of all nodes from this organism in this graph
        """
        if organism not in self.organism_set():
            raise OrganismNotPresentInGraph, organism
        tmp = []
        for node in self.get_nodes():
            if self._organism_from_node(node) != organism:
                continue
            else:
                tmp.append( node )
        # order the nodes and return
        tmp.sort()
        return tmp

    # end of function get_organism_nodes


    def strongest_organism_node(self,organism):
        """
        Get the strongest connected node from a specific organism.

        @type  organism: *
        @param organism: organism identifier (presumably a string)

        @rtype:  *
        @return: One Node identifier
        """
        if organism not in self.organism_set():
            raise OrganismNotPresentInGraph
        tmp = [ ( self.get_node_weighted_connectivity(node), node ) for node in self.get_organism_nodes(organism) ]
        tmp.sort()
        return tmp[-1][1]

    # end of function strongest_organism_node


    def minimal_organism_node_degree(self,organism,degree):
        """
        Remove all nodes of this Organism with a degree less than requested
        
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
        for node in self.get_organism_nodes(organism):
            if self.degree(node) < degree:
                #print "MOND::", node 
                delnodes.append(node)
    
        # delete all the nodes with a to low degree
        for node in delnodes: self.del_node(node)
    
    # end of function minimal_organism_node_degree


    def togff(self,organism=None,gff={},**kwargs):
        """
        Return gff tuple(s) of this graph.

        @attention: MUST BE OVERRIDDEN IN SUBCLASSES
        """
        pass

    # end of function togff


    def _column9data2string(self,data):
        """
        """
        column9string = ""
        # take care for the formatting of column9data
        tmp = []
        if data:
            for key,value in data.iteritems():
                tmp.append( "; %s '%s'" % (key,value) )
            column9string = "".join(tmp)
        # and return the string
        return column9string

    # end of function _column9data2string

# end of class OrganismGraph

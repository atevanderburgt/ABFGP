"""
Pssm-type object functions joined in a class that is inherited in graph classes for Aligment Based Gene Predictions
"""

# Module metadata
__authors__ = "Ate van der Burgt"
__license__ = "MIT"


# Imports
from sets import Set

# Abgp imports
import ordering

class BasalPSSMObjectGraphFunctions:
    """
    """
    def _update_after_changes(self):
        """
        @attention: MUST BE OVERRIDDEN IN THE SUBCLASSES, if the returned value is not correct
        """
        self._update_edge_binary_entropies()
        self._update_node_object()
        self._update_node_pssm()

    # end of function _update_after_changes


    def add_node_and_object(self,node,object):
        """
        Create a single node with a accompagnying object

        @type  node: *
        @param node: One Node identifier

        @type  object: *
        @param object: The object represented by the Node identifier
        """
        if node not in self.get_nodes():
            self.nodes[node] = []
            self._node_object[node] = object
            self._node_pssm[node] = object.pssm_score

    # end of function add_node_and_object


    def _update_node_object(self):
        """
        """
        keys_nodes  = self.get_nodes()
        keys_object = self._node_object.keys()
        # check if object key is present in weights
        for kp in keys_object:
            if kp not in keys_nodes:
                del( self._node_object[kp] )

    # end of function _update_node_pssm


    def _update_node_pssm(self):
        """
        """
        keys_nodes   = self.get_nodes()
        keys_pssm    = self._node_pssm.keys()
        # check if pssm key is present in weights
        for kp in keys_pssm:
            if kp not in keys_nodes:
                del( self._node_pssm[kp] )

    # end of function _update_node_pssm


    def _organism_from_node(self,node):
        """
        Get the organism identifier from the node identifier.

        @attention: MUST BE OVERRIDDEN IN THE SUBCLASSES, if the returned value is not correct

        @type  node: *
        @param node: One Node

        @rtype:  *
        @return: organism identifier (presumably a string)
        """
        return node[0]

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
        return node[0]

    # end of function organism_by_node


    def get_organism_objects(self,organism,order_by=''):
        """
        Get all the OBJECTS of a specific organism from the graph

        @type  organism: *
        @param organism: organism identifier (presumably a string)

        @rtype:  list
        @return: list of all OBJECTS from this organism in this graph
        """
        return ordering.order_list_by_attribute(
                [ self._node_object[node] for node in self.get_organism_nodes(organism) ],
                order_by=order_by
                )

    # end of function get_organism_objects


    def get_node_object(self,node):
        """
        Get the ONLY object of a specific node from the graph

        @type  node: *
        @param node: One Node identifier

        @rtype:  *
        @return: Object that represents this Node
        """
        return self._node_object[node]

    # end of function get_node_object


    def total_pssm(self):
        """
        Get the sum of the individual PSSM scores of the sites in tis graph

        @rtype:  float
        @return: total PSSM score of all (splicesite) nodes in the graph
        """
        return sum( self._node_pssm.values() )

    # end of function total_pssm


    def average_pssm(self):
        """
        Get the average of the PSSM scores of the sites in tis graph

        @rtype:  float
        @return: average PSSM score of the (splicesite) nodes in the graph
        """
        if self.node_count() == 0:
            return 0.0
        else:
            return self.total_pssm() / float(len(self._node_pssm))

    # end of function average_pssm


    def relative_binary_entropy(self):
        """
        TODO
        """
        pass

    # end of function relative_binary_entropy


    def relative_pssm_score(self):
        """
        """
        return self.average_pssm() / self.max_pssm_score()

    # end of function relative_pssm_score


    def relative_distance_weight(self):
        """
        """
        return self.average_weight()

    # end of function relative_distance_weight


    ########################################################################
    ### other Functions                                                  ###
    ########################################################################

    def lowestpssmscoreingraph(self):
        """ """
        node = self.lowestpssmscorenodeingraph()
        if not node: return None
        else:        return self._node_object[node].pssm_score

    # end of function lowestpssmscoreingraph


    def lowestpssmscoreobjectingraph(self):
        """ """
        node = self.lowestpssmscorenodeingraph()
        if not node: return None
        else:        return self._node_object[node]

    # end of function lowestpssmscoreobjectingraph


    def lowestpssmscorenodeingraph(self):
        """ """
        lowestnode = None
        for node, obj in self._node_object.iteritems():
            if not lowestnode:
                lowestnode = node
            elif obj.pssm_score < self._node_object[lowestnode].pssm_score:
                lowestnode = node
            else:
                pass
        return lowestnode

    # end of function lowestpssmscorenodeingraph


    def highestpssmscoreingraph(self):
        """ """
        node = self.highestpssmscorenodeingraph()
        if not node: return None        
        else:        return self._node_object[node].pssm_score

    # end of function highestpssmscoreingraph


    def highestpssmscoreobjectingraph(self):
        """ """
        node = self.highestpssmscorenodeingraph()
        if not node: return None
        else:        return self._node_object[node]

    # end of function highestpssmscoreobjectingraph


    def highestpssmscorenodeingraph(self):
        """ """
        highestnode = None
        for node, obj in self._node_object.iteritems():
            if not highestnode:
                highestnode = node
            elif obj.pssm_score > self._node_object[highestnode].pssm_score:
                highestnode = node
            else:
                pass
        return highestnode

    # end of function highestpssmscorenodeingraph

# end of class BasalPSSMObjectGraphFunctions

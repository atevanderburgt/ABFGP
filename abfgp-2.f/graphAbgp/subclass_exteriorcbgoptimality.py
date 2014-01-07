"""
Functions to asses the optimality of the AlignedStopCodonGraph and the
AlignedTranslationalStartSiteGraph.
"""
# graphAbgp Imports
from exceptions import *

# Python Imports

class ExteriorCbgOptimalityAnalyses:

    def is_optimal_tcode_ratio(self,organism=None,node=None):
        """
        Is the AlignedStopCodonGraph optimal for its tcode ratio (for the given organism or node) ?

        @type  organism: *
        @param organism: One Organism Identifier (or None)

        @type  node: *
        @param node: One Node Identifier (or None)

        @rtype:  NoneBoolean
        @return: True, None or False
        """
        # get organism from node
        if node: organism = self.organism_by_node(node)

        # return NoneBoolean depending on value of total_weight
        ratio = self.average_tcode_entropy()
        if ratio >= self._optimal_max_tcode:        return True
        elif ratio >= self._optimal_min_tcode:
            if organism and not organism in self.organism_set():
                return False
            elif organism:
                theNode = self.node_by_organism(organism)
                ratio = self.average_tcode_entropy(
                        tcode5p=self._tcode5pscore[theNode],
                        tcode3p=self._tcode3pscore[theNode] )
                print "tcode", organism, ratio
                if ratio>=self._optimal_max_tcode:  return True
                else:                               return None
            else:                                   return None
        else:
            if organism and not organism in self.organism_set():
                return False
            elif organism:
                theNode = self.node_by_organism(organism)
                ratio = self.average_tcode_entropy(
                        tcode5p=self._tcode5pscore[theNode],
                        tcode3p=self._tcode3pscore[theNode] )
                if ratio>=self._optimal_max_tcode:  return True
                elif ratio>=self._optimal_min_tcode:return None
                else:                               return False
            else:                                   return False

    # end of function is_optimal_tcode_ratio


    def has_organism_at_least_single_minimal_weight(self,organism,wt):
        """
        """
        if max(self.weights.values()) < wt: return False
        wts = []
        node = self.node_by_organism(organism)
        for node1,node2 in self.weights.keys():
            if node1==node:
                wts.append( self.weights[(node1,node2)] )
        if max(wts) >= wt: return True
        else:              return False

    # end of function has_organism_at_least_single_minimal_weight


    def is_optimal_weight(self,organism=None,node=None):
        """
        Is the AlignedStopCodonGraph optimal for its total weight (for the given organism or node) ?

        @type  organism: *
        @param organism: One Organism Identifier (or None)

        @type  node: *
        @param node: One Node Identifier (or None)

        @rtype:  NoneBoolean
        @return: True, None or False
        """
        # get organism from node
        if node: organism = self.organism_by_node(node)

        # return NoneBoolean depending on value of total_weight
        wt = self.total_weight() / float(self.edge_count())
        if wt >= self._optimal_max_weight:
            return True
        elif wt >= self._optimal_min_weight:
            if organism and not organism in self.organism_set():
                return False
            elif organism:
                if self.has_organism_at_least_single_minimal_weight(organism,1.0):
                    return True
                else:
                    return None
            else:   return None
        else:
            if organism and not organism in self.organism_set():
                return False
            if organism:
                if self.has_organism_at_least_single_minimal_weight(organism,1.0):
                    return True
                else:
                    return False
            else:   return False

    # end of function is_optimal_weight


    def is_optimal_gtgweakestnode(self,organism=None,node=None):
        """
        Is the AlignedStopCodonGraph optimal for the weakestnoderatio of its GeneTreeGraph (for the given organism or node) ?

        @type  organism: *
        @param organism: One Organism Identifier (or None)

        @type  node: *
        @param node: One Node Identifier (or None)

        @rtype:  NoneBoolean
        @return: True, None or False
        """
        # get organism from node
        if node: organism = self.organism_by_node(node)

        # check if _codingblockgraph attribute is applied
        if not self._codingblockgraph:
            message = "attribute '%s' NoneValue in class '%s' object '%s'" % (
                    '_codingblockgraph', self.__class__.__name__, str(self) )
            raise MissingAttributeValue, message

        # return NoneBoolean depending on value of relative
        # connectivity of the weakest node in the GTG
        # representation of the CBG that belongs to this node
        gtg = self._codingblockgraph.genetree()
        gtgweakestnode = gtg.weakest_connected_node()
        gtgweakestorg  = self.organism_by_node(gtgweakestnode)
        ratio = gtg.get_node_weighted_connectivity_observed_vs_expected(gtgweakestnode)
        if ratio >= self._optimal_max_gtgweakest:     return True
        elif ratio >= self._optimal_min_gtgweakest:
            if organism and not organism in self.organism_set():
                                                      return False
            elif organism and gtgweakestorg!=organism:return True
            else:                                     return None
        else:
            if organism and not organism in self.organism_set():
                                                      return False
            elif organism and gtgweakestorg!=organism:return None
            else:                                     return False

    # end of function is_optimal_gtgweakestnode


    def nonebooleanintegermapper(self,noneboolean):
        """
        """
        return nonebooleanintegermapper(noneboolean)
    # end of function nonebooleanintegermapper

# end of function nonebooleanintegermapper



# end of class ExteriorCbgOptimalityAnalyses

########################################################################
### Helper Functions                                                 ###
########################################################################

def nonebooleanintegermapper(noneboolean):
    """
    """
    mapper = {True: 1, None: 0, False: -1}
    try: return mapper[noneboolean]
    except KeyError: raise "nonebooleanintegermapper::KeyError('%s')" % noneboolean

# end of function nonebooleanintegermapper



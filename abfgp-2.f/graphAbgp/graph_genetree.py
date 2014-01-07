################################################################################
### GeneTreeGraph class                                                     ####
################################################################################

# graphAbgp imports
from graph_organism import OrganismGraph
from abgp_analyses import printgtganalyses
from exceptions import *
import graphPlus

# Python imports

# Global variables


class GeneTreeGraph(OrganismGraph):
    """
    GeneTreeGraph (GTG) class, inheriting from OrganismGraph class.
    """
    def __init__(self):
        """
        Initialize a GeneTreeGraph 
        """
        # Initialize as an OrganismGraph
        OrganismGraph.__init__(self)
        # self.weights contains identityscores
        # add similar dicts for aa_identity, bitscore_ratios and nt_identity
        self._aa_identity_percentages = {}
        self._bitscore_ratios         = {} 
        self._nt_identity_percentages = {}

    # end of function __init__


    def node_by_organism(self,organism):
        """
        Get the node identifier belonging to the organism identifier.

        @attention: possible in CodingBlockGraphs, NOT in PacbpCollectionGraphs

        @type  organism: *
        @param organism: Organism identifier

        @rtype:  *
        @return: Node identifier
        """
        for node in self.get_nodes():
            if organism == self._organism_from_node(node):
                return node
        else:
            raise OrganismNotPresentInGraph

    # end of function node_by_organism


    def pairwisecrosscombinations_organism(self,order_by=None):
        """
        Get unique cross combinations of organisms in the graph, ordered on request
    
        @rtype:  list
        @return: list of all unique combinations of two organisms
    
        @attention: returns ORDERED list of ORDERED tuples of organism combinations
        @attention: orderability overwrites OrganismGraph.pairwisecrosscombinations_organism
        @attention: use order_by='identity' to obtain identity ordering except alphabetical
        """
        if order_by == 'identity':
            combinations = []
            for a in self.organism_set():
                for b in self.organism_set():
                    if a == b: continue
                    combi = [a,b]
                    combi.sort()
                    if combi == [a,b]:
                        wt = self.weights[(a,b)]
                        combinations.append( ( wt, tuple(combi) ) )
            # return the ordered (by wt) list of unique combinations
            combinations.sort()
            combinations.reverse()
            return [ combi for (wt,combi) in combinations ]
        else:
            # return basal (alphabetical) ordered pairs of organisms
            return OrganismGraph.pairwisecrosscombinations_organism(self)

    # end of function pairwisecrosscombinations_organism


    def identity(self,organism=None,node=None,edge=None):
        """
        Get the identity of this GeneTreeGraph

        @type  organism: *
        @param organism: organism identifier (or None)

        @type  node: *
        @param node: node identifier (or None)

        @type  edge: tuple 
        @param edge: tuple of 2 Node identifiers (or None)

        @attention: when no organism specified, this function is an alias for `average_weight()`
        """
        if organism:
            if organism not in self.organism_set():
                raise OrganismNotPresentInGraph
            else:
                node = self.node_by_organism(organism)
                return self.get_node_weighted_connectivity(node) * self.total_weight()
        elif node:
            if node not in self.get_nodes():
                raise NodeNotPresentInGraph
            else:
                return self.get_node_weighted_connectivity(node) * self.total_weight()
        elif edge:
            if not self.weights.has_key(edge):
                raise EdgeNotPresentInGraph
            else:
                return self.weights[edge]
        else:
            return self.average_weight()

    # end of function identity


    def aaidentity(self,organism=None,node=None,edge=None):
        """
        Get the AAidentity of this GeneTreeGraph

        @type  organism: *
        @param organism: organism identifier (or None)

        @type  node: *
        @param node: node identifier (or None)

        @type  edge: tuple
        @param edge: tuple of 2 Node identifiers (or None)

        @attention: when no Organism, Node or Edge specified, average aaidentity% is returned
        """
        if not self._aa_identity_percentages: return 0.0
        if organism:
            if organism not in self.organism_set():
                raise OrganismNotPresentInGraph
            else:
                node = self.node_by_organism(organism)
        elif node:
            if node not in self.get_nodes():
                raise NodeNotPresentInGraph
        elif edge:
            if not self.weights.has_key(edge):
                raise EdgeNotPresentInGraph
            else:
                return self._aa_identity_percentages[edge]
        else:
            pass
        # check if a node was specified
        if node:
            ll = []
            for (n1,n2), aaident in self._aa_identity_percentages.iteritems():
                if n1 == node: ll.append(aaident)
            return sum(ll)/len(ll)
        else:
            return sum(self._aa_identity_percentages.values())/len(self._aa_identity_percentages) 

    # end of function aaidentity


    def ntidentity(self,organism=None,node=None,edge=None):
        """
        Get the DNA identity% of this GeneTreeGraph

        @type  organism: *
        @param organism: organism identifier (or None)

        @type  node: *
        @param node: node identifier (or None)

        @type  edge: tuple
        @param edge: tuple of 2 Node identifiers (or None)

        @attention: when no Organism, Node or Edge specified, average ntidentity% is returned
        """
        if not self._nt_identity_percentages: return 0.0
        if organism:
            if organism not in self.organism_set():
                raise OrganismNotPresentInGraph
            else:
                node = self.node_by_organism(organism)
        elif node:
            if node not in self.get_nodes():
                raise NodeNotPresentInGraph
        elif edge:
            if not self.weights.has_key(edge):
                raise EdgeNotPresentInGraph
            else:
                return self._nt_identity_percentages[edge]
        else:
            pass
        # check if a node was specified
        if node:
            ll = []
            for (n1,n2), ntident in self._nt_identity_percentages.iteritems():
                if n1 == node: ll.append(ntident)
            return sum(ll)/len(ll)
        else:
            return sum(self._nt_identity_percentages.values())/len(self._nt_identity_percentages)

    # end of function ntidentity


    def bitscoreratio(self,organism=None,node=None,edge=None):
        """
        Get the bitscore ratio of this GeneTreeGraph
        
        @type  organism: *
        @param organism: organism identifier (or None)
        
        @type  node: *
        @param node: node identifier (or None)

        @type  edge: tuple
        @param edge: tuple of 2 Node identifiers (or None)

        @attention: when no Organism, Node or Edge specified, average bitscoreratio is returned
        """ 
        if not self._bitscore_ratios: return 0.0
        if organism:
            if organism not in self.organism_set():
                raise OrganismNotPresentInGraph
            else:
                node = self.node_by_organism(organism)
        elif node:
            if node not in self.get_nodes():
                raise NodeNotPresentInGraph
        elif edge:
            if not self.weights.has_key(edge):
                raise EdgeNotPresentInGraph
            else:
                return self._bitscore_ratios[edge]
        else:
            pass
        # check if a node was specified
        if node:
            ll = []
            for (n1,n2), bsr in self._bitscore_ratios.iteritems():
                if n1 == node: ll.append(bsr)
            return sum(ll)/len(ll)
        else:
            return sum(self._bitscore_ratios.values())/len(self._bitscore_ratios)

    # end of function bitscoreratio 


    def bitscoreratios(self):
        """
        Return reversed ordered list of bitscore-integers (0.56->56)
        """
        bsrs = [ int(self._bitscore_ratios[(n1,n2)]*100) for (n1,n2) in self.pairwisecrosscombinations_node() ]
        bsrs.sort()
        bsrs.reverse()
        return bsrs

    # end of function bitscoreratios


    def aaidentities(self):
        """ 
        Return reversed ordered list of aaidentity-integers (0.56->56)
        """
        aaids = [ int(self._aa_identity_percentages[(n1,n2)]*100) for (n1,n2) in self.pairwisecrosscombinations_node() ]
        aaids.sort()
        aaids.reverse()
        return aaids 

    # end of function aaidentities 


    def ntidentities(self):
        """
        Return reversed ordered list of DNA identity-integers (0.56->56)
        """
        ntids = [ int(self._nt_identity_percentages[(n1,n2)]*100) for (n1,n2) in self.pairwisecrosscombinations_node() ]
        ntids.sort()
        ntids.reverse()
        return ntids

    # end of function ntidentities


    def __str__(self):
        """
        """
        # order the nodes per weighted connectedness
        tmp = []
        for node in self.get_nodes():
            tmp.append( ( self.get_node_weighted_connectivity(node), node ) )
        tmp.sort()
        tmp.reverse()
        nodes = [ self._organism_from_node(node) for (score,node) in tmp ]
        # order the weights as well
        wts = []
        for (node1,node2) in self.pairwisecrosscombinations_node():
            wts.append( int( self.weights[(node1,node2)]*100 ) )
        wts.sort()
        wts.reverse()
        return "<%s identity=%1.2f [%s] (%s)>" % (
            self.__class__.__name__,
            self.average_weight(),
            " ".join([ str(node) for node in nodes ]),
            "-".join([str(w) for w in wts]),
            )

    # end of function __str__


    def phylogeneticaltree(self):
        """
        """
        return graphPlus.comparison.phylogeneticaltree(self)

    # end of function phylogeneticaltree


    def graphalignmentdifference(self,gtg,**kwargs):
        """
        """
        if gtg.__class__.__name__ == 'GeneTreeGraph':
            pass
        elif gtg.__class__.__name__ == 'CodingBlockGraph':
            gtg = CodingBlockGraph2GeneTreeGraph(gtg)
        else:
            message = "argument `gtg` is not a GeneTreeGraph or CodingBlockGraph object"
            raise InproperlyAppliedArgument, message 
        return graphPlus.comparison.identical_graph_weight_difference(self,gtg,**kwargs)

    # end of function graphalignmentdifference


    def absolutegraphalignmentdifference(self,gtg,**kwargs):
        """
        """
        if gtg.__class__.__name__ == 'GeneTreeGraph':
            pass
        elif gtg.__class__.__name__ == 'CodingBlockGraph':
            gtg = CodingBlockGraph2GeneTreeGraph(gtg)
        else:
            message = "argument `gtg` is not a GeneTreeGraph or CodingBlockGraph object"
            raise InproperlyAppliedArgument, message
        absdif = graphPlus.comparison.identical_graph_absolute_weight_difference(self,gtg,**kwargs)
        if gtg.node_count() == 0:
            return 0.0
        else:
            return absdif / float( self.edge_count() )

    # end of function graphalignmentdifference


    def printgtganalyses(self,prefix=None):
        """
        Print tab delimited line of GeneTreeGraph characteristics

        @attention: see abgp_analyses.printcbganalyses()
        """
        printgtganalyses(self,prefix=prefix)

    # end of function printgtganalyses

# end of class GeneTreeGraph

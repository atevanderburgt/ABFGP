################################################################################
### AlignedStopCodonGraph class                                             ####
################################################################################

# graphAbgp imports
from graph_organism import OrganismGraph
from subclass_tcodedata import GraphTcodeDataAccesFunctions
from subclass_exteriorcbgoptimality import ExteriorCbgOptimalityAnalyses
import codingblock_collectionharvesting
from exceptions import *

# Abgp Imports
from settings.alignedstopcodongraph import *
from settings.executables import (
    TCODE_MAX_NONCODING,
    TCODE_MIN_CODING, 
)

# Python imports


class AlignedStopCodonGraph(OrganismGraph,GraphTcodeDataAccesFunctions,ExteriorCbgOptimalityAnalyses):
    """
    AlignedStopCodonGraph (ASCG) class, inheriting from graphPlus class.
    """

    def __init__(self,cbg,tcode_5p_windowsize=201,tcode_3p_windowsize=201):
        """
        Initialize a AlignedStopCodonGraph
        """
        OrganismGraph.__init__(self)
        # attribute for storing the CBG itself
        self._codingblockgraph = cbg 
        # attributes for TCODE data
        self._tcode5pscore = {}
        self._tcode3pscore = {}
        self._TCODE_5P_WINDOWSIZE = tcode_5p_windowsize
        self._TCODE_3P_WINDOWSIZE = tcode_3p_windowsize

        # is_optimal_xxx thresholds
        self._optimal_min_tcode      = ALIGNEDSTOPCODONGRAPH_OPTIMALITY_MIN_TCODE
        self._optimal_max_tcode      = ALIGNEDSTOPCODONGRAPH_OPTIMALITY_MAX_TCODE
        self._optimal_min_weight     = ALIGNEDSTOPCODONGRAPH_OPTIMALITY_MIN_WEIGHT
        self._optimal_max_weight     = ALIGNEDSTOPCODONGRAPH_OPTIMALITY_MAX_WEIGHT 
        self._optimal_min_gtgweakest = ALIGNEDSTOPCODONGRAPH_OPTIMALITY_MIN_GTGWEAKEST
        self._optimal_max_gtgweakest = ALIGNEDSTOPCODONGRAPH_OPTIMALITY_MAX_GTGWEAKEST

        # run function codingblock_collectionharvesting.align_stop_codons
        # to fully initialize the object
        self = codingblock_collectionharvesting.align_stop_codons(cbg,self)

    # end of function __init__

    ########################################################################
    ### Build-in functions                                               ###
    ########################################################################

    def __str__(self):
        """ Nicely formatted oneliner reprenting this object """
        nodes = self.get_ordered_nodes()
        nodes = [ '%s(%s):%s' % (node[0],node[1],node[3]) for node in nodes ]
        nodes =  " ".join([ n for n in nodes])
        # and return a nicely formatted string
        return "<%s N%s (tw:%2.1f betw:%2.1f tcode:%1.2f-(%1.2f)-%1.2f) distOMSR:%s [%s]>" % (
                self.__class__.__name__,
                self.node_count(),
                self.total_weight(),
                self.total_binary_entropy(),
                self.average_5p_tcode_score(),
                self.average_tcode_entropy(),
                self.average_3p_tcode_score(),
                self.stopcodon2omsrdistance(),
                nodes,
                )

    # end of function __str__

    ########################################################################
    ### Functions mapped from the CodingBlockGraphs                      ###
    ########################################################################

    def _organism_from_node(self,node):
        """ @attention: see OrganismGraph._organism_from_node """
        return node[0]
    # end of function _organism_from_node 

    def organism_by_node(self,node):
        """ @attention: see OrganismGraph.organism_by_node """
        return node[0]
    # end of function organism_by_node

    def node_by_organism(self,organism):
        """ @attention: see OrganismGraph.node_by_organism """
        for node in self.get_nodes():
            if node[0] == organism:
                return node
        else:
            raise OrganismNotPresentInGraph
    # end of function organism_by_node

    def _update_after_changes(self):
        """
        @attention: MUST BE OVERRIDDEN IN THE SUBCLASSES, if the returned value is not correct
        """
        self._update_edge_binary_entropies()

    # end of function _update_after_changes


    def average_tcode_entropy(self,tcode5p=None,tcode3p=None):
        """
        Calculate ratio between Coding (left) and Non-Coding (rigth) side of Stop-Codon
        """
        if tcode5p and tcode3p:
            return average_tcode_entropy_stopcodon(tcode5p=tcode5p,tcode3p=tcode3p)
        elif self.get_nodes():
            return average_tcode_entropy_stopcodon(
                    tcode5p = sum(self._tcode5pscore.values()) / float(self.node_count()),
                    tcode3p = sum(self._tcode3pscore.values()) / float(self.node_count()),
                    )
        else:
            return average_tcode_entropy_stopcodon()

    # end of function average_tcode_entropy
    

    def togff(self,gff={},organism=None):
        """
        Create gff tuple for aligned stop codon of a specific organism

        @type  gff: dictionary
        @param gff: overwrite default gff data, keys: ('fstrand','fhase','fref',etc...)

        @type  organism: * (presumably string)
        @param organism: Organism identifier to make the gff for

        @rtype:  tuple
        @return: gff tuple with 9 elements
        """
        # get node of this organism
        orgnode = self.get_organism_nodes(organism)[0]
        # default gff data
        _gff = {
            'fstrand'   : '+',
            'fphase'    : '.',
            'fref'      : 'none',
            'fsource'   : 'ABFGP',
            'fmethod'   : 'aligned_stop_codon',
            'gclass'    : 'AlignedStopCodon',
            'gname'     : "-".join([ str(elem) for elem in orgnode ]),
            'column9data' : {}
        }
        # overwrite defaults with gff dictionary passed to this function
        _gff.update(gff)
        # and return data
        return (
            _gff['fref'],
            _gff['fsource'],
            _gff['fmethod'],
            orgnode[3]+1,
            orgnode[3]+3,
            "%1.2f" % self.total_weight(),
            _gff['fstrand'],
            _gff['fphase'],
            "%s %s%s" % ( _gff['gclass'], _gff['gname'], self._column9data2string(_gff['column9data']) )
        )

    # end of function togff

    ########################################################################
    ### Function to asses the optimality of the AlignedStopCodonGraph    ###
    ########################################################################

    def is_optimal(self,organism=None,node=None):
        """
        Is the AlignedStopCodonGraph optimal (for the given organism or node) ?

        @type  organism: *
        @param organism: One Organism Identifier (or None)

        @type  node: *
        @param node: One Node Identifier (or None)

        @rtype:  NoneBoolean
        @return: True, None or False
        """
        # get organism from node
        if node: organism = self.organism_by_node(node)

        # first check: how much does the StopCodons in average stick out to the OMSR?
        omsrdist   = self.stopcodon2omsrdistance()
        # get nonebooleans for tcode, weight and gtgweakest
        tcode      = self.nonebooleanintegermapper( self.is_optimal_tcode_ratio() )
        weight     = self.nonebooleanintegermapper( self.is_optimal_weight() ) 
        gtgweakest = self.nonebooleanintegermapper( self.is_optimal_gtgweakestnode() )
        totalsum   = sum([tcode,weight,gtgweakest]) 

        # now translate into NoneBoolean outcome; weight MUST be >=0 too!
        if omsrdist > MAX_OMSR_DIST_FUNCTION(self) and organism:    pass
        elif omsrdist > MAX_OMSR_DIST_FUNCTION(self):               return False
        elif totalsum > 0 and weight >= 0:                          return True
        elif totalsum == 0 and min([tcode,weight,gtgweakest]) == 0: return True
        elif not organism:                                          return False
        else:                                                       pass 

        # if this point in the function is reached, do a check for
        # the organism that is asked for
        if omsrdist > MAX_OMSR_DIST_FUNCTION(self):
            # check the MAXSR for this specific organism
            maxsrdist = self.stopcodon2maxsrdistance(organism)
            if maxsrdist > ( MAX_OMSR_DIST_FUNCTION(self) / 2 ):
                return False
            else:
                pass

        tcode      = self.nonebooleanintegermapper( self.is_optimal_tcode_ratio(organism=organism) )
        weight     = self.nonebooleanintegermapper( self.is_optimal_weight(organism=organism) )
        gtgweakest = self.nonebooleanintegermapper( self.is_optimal_gtgweakestnode(organism=organism) )
        totalsum   = sum([tcode,weight,gtgweakest])

        # now translate into NoneBoolean outcome
        if totalsum >= 1:                                           return True
        elif totalsum == 0:                                         return None 
        else:                                                       return False 

    # end of function is_optimal


    def stopcodon2omsrdistance(self,organism=None,node=None):
        """
        Obtain the (AA) distance between the StopCodon and the OMSR of the CBG

        @type  organism: *
        @param organism: One Organism Identifier (or None)

        @type  node: *
        @param node: One Node Identifier (or None)

        @rtype:  Integer (Positive) 
        @return: AA distance between CBGs OMSR and the StopCodon of the Orf(s) of this CBG
        """
        # get organism from node
        if node: organism = self.organism_by_node(node)
        # get the OMSR of this CBG
        omsr = self._codingblockgraph.overall_minimal_spanning_range()
        dists = [ ]
        for org in self.organism_set():
            if organism and org != organism: continue
            stopnode = self.node_by_organism(org)
            cbgnode  = ( stopnode[0], stopnode[1] )
            aapos    = stopnode[2]
            # calculate distance; correct omsr+1 (omsr isa python list range)
            # and correct dist by -1 for stop codon that is not a AA it self
            dist     = aapos - ( max(omsr[cbgnode])+1 ) -1
            dists.append( dist )
        # return average distance as an integer
        return int( float(sum(dists)) / len(dists) )

    # end of function stopcodon2omsrdistance 


    def stopcodon2maxsrdistance(self,organism):
        """
        Obtain the (AA) distance between the StopCodon and the MAXSR of the CBG for this organism

        @type  organism: *
        @param organism: One Organism Identifier (or None)

        @rtype:  Integer (Positive)
        @return: AA distance between CBGs MAXSR and the StopCodon of the Orf(s) of this CBG
        """
        # get the MAXSR of this CBG
        maxsr    = self._codingblockgraph.maximal_spanning_range(organism=organism)
        stopnode = self.node_by_organism(organism)
        aapos    = stopnode[2]

        # calculate distance; correct omsr+1 (omsr isa python list range)
        # and correct dist by -1 for stop codon that is not a AA it self
        dist     = aapos - ( max(maxsr)+1 ) -1
        # return distance as an integer
        return dist

    # end of function stopcodon2maxsrdistance

# end of class AlignedStopCodonGraph


########################################################################
### Helper Functions for AlignedStopCodonGraph                       ###
########################################################################

def average_tcode_entropy_stopcodon(tcode5p=0.0,tcode3p=0.0):
    """
    Ratio between Coding (left) and Non-Coding (rigth) side of Stop-Codon

    @type  tcode5p: positive float (0.450 - 1.320)
    @param tcode5p: average tcode score 5'/left/upstream of specific site

    @type  tcode3p: positive float (0.450 - 1.320)
    @param tcode3p: average tcode score 3'/right/downstream of specific site

    @rtype:  float (-0.450 - 2.160)
    @return: score attached to the likelihood that site is a stopcodon; >1.0 is likely

    """
    return 1.0 + ( TCODE_MAX_NONCODING - tcode3p ) * 2.0 + ( tcode5p -  TCODE_MAX_NONCODING )

# end of function average_tcode_entropy_stopcodon


def nonebooleanintegermapper(noneboolean):
    """
    """
    mapper = {True: 1, None: 0, False: -1}
    try: return mapper[noneboolean]
    except KeyError: raise "nonebooleanintegermapper::KeyError('%s')" % noneboolean

# end of function nonebooleanintegermapper


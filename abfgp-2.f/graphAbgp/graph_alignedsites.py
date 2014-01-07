################################################################################
### AlignedPssmObjectGraph class and inheriting classes                     ####
################################################################################

# graphAbgp imports
from graph_organism import OrganismGraph
from subclass_pssmobjects import BasalPSSMObjectGraphFunctions
from subclass_tcodedata import GraphTcodeDataAccesFunctions
from graph_pssmcollections import _format_splicesite_nodes_to_txt
import codingblock_collectionharvesting
from exceptions import *


from graph_alignedstopcodongraph import AlignedStopCodonGraph
from graph_alignedtranslationalstartsitegraph import AlignedTranslationalStartSiteGraph

# Abgp Imports
from settings.alignedstopcodongraph import *


# Python imports
from sets import Set


class AlignedPssmObjectGraph(OrganismGraph,BasalPSSMObjectGraphFunctions):
    """
    AlignedPssmObjectGraph (APOG) class, inheriting from OrganismGraph and BasalPSSMObjectGraphFunctions class.
    """
    def __init__(self,max_node_count=None,min_pssm_score=None,aligned_site_aa_offset=None):
        """
        Initialize a AlignedPssmObjectGraph 
        """
        OrganismGraph.__init__(self)
        self._node_pssm = {}
        self._node_object = {}
        self.MIN_PSSM_SCORE = min_pssm_score
        self.MAX_NODE_COUNT = max_node_count
        self.ALIGNED_SITE_AA_OFFSET = aligned_site_aa_offset

    # end of function __init__

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

    def add_edge(self, u, v, wt=1, entropy=1.0):
        """
        """
        OrganismGraph.add_edge(self,u,v,wt=wt)
        self._edge_binary_entropies[(u,v)] = (entropy,entropy)
        self._edge_binary_entropies[(v,u)] = (entropy,entropy)

    # end of function add_edge

    def makecompletegraph(self):
        """
        """
        OrganismGraph.makecompletegraph(self)
        for u,v in self.weights.keys():
            if not self._edge_binary_entropies.has_key((u,v)):
                self._edge_binary_entropies[(u,v)] = (1.0,1.0)

    # end of function makecompletegraph

    ########################################################################
    ### other Functions                                                  ###
    ########################################################################

    def is_aligned(self):
        """
        """
        if len(self.nodes) == 0:
            return False
        else:
            if self.connectivitysaturation() == 1.0:
                return True
            else:
                return False

    # end of function is_aligned


    def cumulative_score(self,weight_pssm_part=1.5,weight_entropy_part=30.0,weight_weight_part=10.0,weight_edges_part=1.0):
        """
        Cumulative score of a group of aligned PSSM objects

        @type  weight_pssm_part:    float
        @param weight_pssm_part:    relative contribution of pssm score [ ~-3.0 .. ~10.0 ]

        @type  weight_entropy_part: float
        @param weight_entropy_part: relative contribution of binary entropy [ 0.0 .. 1.0 ]

        @type  weight_weight_part:  float
        @param weight_weight_part:  relative contribution of graph weight [ 0.0 .. 1.0 ]

        @type  weight_edges_part:  float
        @param weight_edges_part:  relative contribution of edge count [ 1,3,6,10,15 for n=2,3,4,5,6 ]


        @rtype:  float
        @return: cumulative score of this aligned set of sites

        @attention: CAN BE OVERRIDDEN IN THE SUBCLASSES!!
        """
        if self.node_count() == 0:
            return 0.0
        else:
            node_ratio = float(self.node_count()) / float(self.MAX_NODE_COUNT)
            pssm_correction = self.total_pssm() / abs( self.total_pssm() )
            return ( (
                ( float(self.edge_count())      * weight_edges_part ) +\
                ( self.average_binary_entropy() * weight_entropy_part ) +\
                ( self.average_pssm()           * weight_pssm_part ) +\
                ( self.average_weight()         * weight_weight_part ) ) *\
                node_ratio ) * pssm_correction

    # end of function cumulative_score


    def __str__(self):
        """ """
        nodes = self.get_nodes()
        nodes.sort()
        # try some string shortening here
        try:    nodes = [ '%s_%s' % (node[0],node[-1]) for node in nodes ]
        except: pass

        # hack for AlignedSpliceSiteGraphs
        if self.__class__.__name__ in [
                'AlignedSpliceSiteGraph',
                'AlignedDonorSiteGraph',
                'AlignedAcceptorSiteGraph',
                'AlignedSpliceSiteWithPhaseShiftGraph',
                'AlignedDonorSiteWithPhaseShiftGraph',
                'AlignedAcceptorSiteWithPhaseShiftGraph',
                'AlignedCbg3pMixedSiteGraph',
                'AlignedCbg5pMixedSiteGraph',
                'AlignedCbg5pBoundarySiteGraph',
                'AlignedCbg3pBoundarySiteGraph',
                ]:
            nodes = _format_splicesite_nodes_to_txt(self)

        # try phase grabbing; not all pssm objects have self.phase()!
        try:    phase = self.phase()
        except: phase = ''
        # try float formatting of cumulative_score
        try:    cumulative = "%2.2f" % self.cumulative_score()
        except: cumulative = str(self.cumulative_score())
        return "<%s N%s-O%s-E%s [o=%s d=%s p=%1.1f], %s score=%s (csat: %1.1f tw: %1.1f betw: %1.1f pssmtw: %1.1f) %s >" % (
                self.__class__.__name__,
                self.node_count(),
                self.organism_set_size(),
                self.edge_count(),
                self.MAX_NODE_COUNT,
                self.ALIGNED_SITE_AA_OFFSET,
                self.MIN_PSSM_SCORE,
                phase,
                cumulative,
                self.connectivitysaturation(),
                self.total_weight(),
                self.total_binary_entropy(),
                self.total_pssm(),
                nodes
                )

    # end of function __str__


    def togff(self,gff={},organism=None):
        """
        Create gff tuple for aligned sites of a specific organism

        @type  gff: dictionary
        @param gff: overwrite default gff data, keys: ('fstrand','fphase','fref',etc...)

        @type  organism: * (presumably string)
        @param organism: Organism identifier to make the gff for

        @rtype:  tuple
        @return: gff tuple with 9 elements
        """
        if organism in self.organism_set():
            # get node of this organism
            orgnode = self.get_organism_nodes(organism)[0]
            # get (single and only) object of this organism
            item = self.get_organism_objects(organism)[0]
            gff['gname'] = "-".join([ str(elem) for elem in orgnode ])
            # add some metadata
            if not gff.has_key('column9data'): gff['column9data'] = {}
            gff['column9data'].update( {
                    'BinaryEntropy':    "%1.2f" % self.total_binary_entropy(),
                    'PSSMScore':        "%1.2f" % self.total_pssm(),
                    'TotalWeight':      "%1.2f" % self.total_weight(),
                    } )
            # get list of the gff upon hard-updating the fscore
            gff_tup = list( item.togff( gff=gff ) )
            # IMPORTANT!! update the fscore of this object
            gff_tup[5] = self.cumulative_score()
            # and return gff tuple of this item
            return tuple( gff_tup )
        else:
            # return empty tuple; this will be filtered out lateron
            return tuple()

    # end of function togff

# end of class AlignedPssmObjectGraph


class AlignedStartCodonGraph(OrganismGraph,GraphTcodeDataAccesFunctions):
    """
    AlignedStartCodonGraph (ASCG) class, inheriting from graphPlus class.
    """

    def __init__(self,tcode_5p_windowsize=201,tcode_3p_windowsize=201):
        """
		Initialize a AlignedStartCodonGraph
        """
        OrganismGraph.__init__(self)
        self._tcode5pscore = {}
        self._tcode3pscore = {}
        self._TCODE_5P_WINDOWSIZE = tcode_5p_windowsize
        self._TCODE_3P_WINDOWSIZE = tcode_3p_windowsize

    # end of function __init__


    def _update_after_changes(self):
        """
        @attention: MUST BE OVERRIDDEN IN THE SUBCLASSES, if the returned value is not correct
        """
        self._update_edge_binary_entropies()

    # end of function _update_after_changes


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


    def average_tcode_entropy(self):
        """
        Calculate ratio between Coding (rigth) and Non-Coding (left) side of Start-Codon
        """
        if self.get_nodes:
            return average_tcode_entropy_startcodon(
                    tcode5p = sum(self._tcode5pscore.values()) / float(self.node_count()),
                    tcode3p = sum(self._tcode3pscore.values()) / float(self.node_count()),
                    )
        else:
            return average_tcode_entropy_startcodon()

    # end of function average_tcode_entropy


# end of class AlignedStartCodonGraph



class AlignedCbgBoundarySiteGraph(AlignedPssmObjectGraph):
    """
    AlignedCbgBoundarySiteGraph (ABSG) class, inheriting from AlignedPssmObjectGraph
    """
    def __init__(self,**kwargs):
        """
		Initialize a AlignedCbgBoundarySiteGraph
        """
        AlignedPssmObjectGraph.__init__(self,**kwargs)
        self._node_phase = {}
        # which type? -> not needed in this type of graph, for compatibilty purposes
        self._type = None

    # end of function __init__


    def is_aligned(self):
        """
        @attention: overrides AlignedPssmObjectGraph.is_aligned()
        """
        if len(self.nodes) == 0:    return False
        else:                       return True

    # end of function is_aligned


    def phase(self):
        """
        Phase function of AlignedCbgBoundarySiteGraph; return None (this type of sites have a None phase)
        """
        return None

    # end of function phase


    def cumulative_score(self):
        """
        Cumulative score of a group of aligned AlignedSpliceSites

        @attention: See AlignedPssmObjectGraph.cumulative_score, CAN BE OVERRIDDEN HERE IF DESIRED!!
        """
        return AlignedPssmObjectGraph.cumulative_score(self)

    # end of function cumulative_score


    def togff(self,gff={},organism=None):
        """
        Create gff tuple for aligned splice site of a specific organism

        @attention: See AlignedPssmObjectGraph.togff, CAN BE OVERRIDDEN HERE IF DESIRED!!
        """
        return AlignedPssmObjectGraph.togff(self,organism=organism,gff=gff)

    # end of function togff


# end of class AlignedCbgBoundarySiteGraph


class AlignedCbgMixedSiteGraph(AlignedCbgBoundarySiteGraph):
    """
    AlignedCbgMixedSiteGraph (AMSG) class, inheriting from AlignedCbgBoundarySiteGraph
    """
    def phase(self):
        """
        Phase function of AlignedCbgMixedSiteGraph; return the unique phase of the splice sites in the mixed collection
        """
        if len(self._node_object) == 0:
            # no objects in graph (yet)
            return None
        else:
            phases = list( Set( [ o.phase for o in self._node_object.values() ] ) )
            if len(phases) == 1:
                # hmmm not expected, but not a (very) big disaster
                return phases[0]
            elif len(phases) == 2 and None in phases:
                # a list of [None, x] with x in [0,1,2]
                return phases
            else:
                raise IncompatibleSpliceSitePhases

    # end of function phase

# end of class AlignedCbgMixedSiteGraph



class AlignedCbg5pBoundarySiteGraph(AlignedCbgBoundarySiteGraph):
    pass

# end of class AlignedCbg5pBoundarySiteGraph


class AlignedCbg3pBoundarySiteGraph(AlignedCbgBoundarySiteGraph):
    pass

# end of class AlignedCbg3pBoundarySiteGraph


class AlignedCbg5pMixedSiteGraph(AlignedCbgMixedSiteGraph):
    pass

# end of class AlignedCbg5pMixedSiteGraph


class AlignedCbg3pMixedSiteGraph(AlignedCbgMixedSiteGraph):
    pass

# end of class AlignedCbg3pMixedSiteGraph


class AlignedSpliceSiteGraph(AlignedPssmObjectGraph):
    """
    AlignedSpliceSiteGraph (ASSG) class, inheriting from AlignedPssmObjectGraph
    """

    def __init__(self,**kwargs):
        """
		Initialize a AlignedSpliceSiteGraph
        """
        AlignedPssmObjectGraph.__init__(self,**kwargs)
        self._node_phase = {}
        # which type? donors or acceptors?
        self._type = None    # to be overwritten in the subclasses

    # end of function __init__


    def is_aligned(self):
        """
        @attention: overrides AlignedPssmObjectGraph.is_aligned()
        """
        if len(self.nodes) == 0:
            return False
        else:
            try:
                self.phase()
                return True
            except IncompatibleSpliceSitePhases:
                return False

    # end of function is_aligned


    def phase(self):
        """
        Phase function of AlignedSpliceSiteGraph; return single unique phase

        """
        if len( Set( [ o.phase for o in self._node_object.values() ] ) ) != 1:
            raise IncompatibleSpliceSitePhases
        elif len(self._node_object) == 0:
            # no objects in graph (yet)
            return None
        else:
            return self._node_object.values()[0].phase

    # end of function phase


    def cumulative_score(self):
        """
        Cumulative score of a group of aligned AlignedSpliceSites

        @attention: See AlignedPssmObjectGraph.cumulative_score, CAN BE OVERRIDDEN HERE IF DESIRED!!
        """
        return AlignedPssmObjectGraph.cumulative_score(self)

    # end of function cumulative_score


    def togff(self,gff={},organism=None):
        """
        Create gff tuple for aligned splice site of a specific organism

        @attention: See AlignedPssmObjectGraph.togff, CAN BE OVERRIDDEN HERE IF DESIRED!!
        """
        return AlignedPssmObjectGraph.togff(self,organism=organism,gff=gff)

    # end of function togff


    def replace_splice_site_object(self,oldnode,oldobj,newobj,verbose=False):
        """
        HARD-replace a certain object with another splice site object

        @attention: use with care!
        """
        # a splice site node is represented as ( orgID, orfID, pos )
        newnode = ( oldnode[0], oldnode[1], newobj.pos )

        ########################################################################
        if verbose:
            print self
            print "oldnode:", oldnode
            print "newnode:", newnode
        ########################################################################

        # first replace all the edge attributes
        edges = self.weights.keys()
        current_average_weight = self.average_weight()
        for (n1,n2) in edges:
            if n1==oldnode:
                # create new dict entries
                if current_average_weight == 1.0:
                    # here, we can accurately define a new
                    # weight value based on distance
                    dist = float( abs(oldobj.pos - newobj.pos) / 3 )
                    self.weights[(newnode,n2)] = 1.0 / dist
                else:
                    self.weights[(newnode,n2)] =\
                        self.weights[(n1,n2)]
                # just replace edge binary entropy values;
                # this is not 100% correct
                self._edge_binary_entropies[(newnode,n2)] =\
                    self._edge_binary_entropies[(n1,n2)]
                # remove old dict entries
                del( self.weights[(n1,n2)] )
                del( self._edge_binary_entropies[(n1,n2)] )
            elif n2==oldnode:
                # create new dict entries
                if current_average_weight == 1.0:
                    # here, we can accurately define a new
                    # weight value based on distance
                    dist = float( abs(oldobj.pos - newobj.pos) / 3 )
                    self.weights[(n1,newnode)] = 1.0 / dist
                else:
                    self.weights[(n1,newnode)] =\
                        self.weights[(n1,n2)]
                self._edge_binary_entropies[(n1,newnode)] =\
                    self._edge_binary_entropies[(n1,n2)]
                # remove old dict entries
                del( self.weights[(n1,n2)] )
                del( self._edge_binary_entropies[(n1,n2)] )
            else:
                pass

        # remove all the oldnode attributes
        del( self.nodes[oldnode] )
        del( self._node_pssm[oldnode] )
        del( self._node_object[oldnode] )

        # add the newnode attributes
        self.nodes[newnode] = self.get_nodes()
        self._node_pssm[newnode] = newobj.pssm_score
        self._node_object[newnode] = newobj

        ########################################################################
        if verbose: print self
        ########################################################################

    # end of function replace_splice_site_object

# end of class AlignedSpliceSiteGraph


class AlignedDonorSiteGraph(AlignedSpliceSiteGraph):
    """
    AlignedDonorSiteGraph (ADSG) class, inheriting from AlignedSpliceSiteGraph class.
    """
    def __init__(self,**kwargs):
        """
        Initialize as a AlignedSpliceSiteGraph
        """
        AlignedSpliceSiteGraph.__init__(self,**kwargs)
        # which type donors or acceptors?
        self._type = "Donor"

    # end of function __init__

# end of class AlignedDonorSiteGraph


class AlignedAcceptorSiteGraph(AlignedSpliceSiteGraph):
    """
    AlignedAcceptorSiteGraph (AASG) class, inheriting from AlignedSpliceSiteGraph class.
    """
    def __init__(self,**kwargs):
        """
        Initialize as a AlignedSpliceSiteGraph
        """
        AlignedSpliceSiteGraph.__init__(self,**kwargs)
        # which type donors or acceptors?
        self._type = "Acceptor"

    # end of function __init__

# end of class AlignedAcceptorSiteGraph


class AlignedSpliceSiteWithPhaseShiftGraph(AlignedSpliceSiteGraph):
    """
    AlignedSpliceSiteWithPhaseShiftGraph (ASSPSG) class, inheriting from AlignedSpliceSiteGraph class.
    """
    def __init__(self,**kwargs):
        """
        Initialize as a AlignedSpliceSiteGraph
        """
        AlignedSpliceSiteGraph.__init__(self,**kwargs)

    # end of function __init__


    def phase(self):
        """
        Phase function of AlignedSpliceSiteWithPhaseShiftGraph; return ordered list of two phases 
        """
        phases = list( Set( [ o.phase for o in self._node_object.values() ] ) )
        phases.sort() 
        if len(phases)== 1:
            raise "UNEXPECTED uniform phase in %s" % self.__class__.__name__ 
        elif len(phases)==2 and None not in phases:
            return phases
        elif len(self._node_object) == 0:
            # no objects in graph (yet)
            return None
        else:
           raise "UNEXPECTED phases in %s: %s" % (self.__class__.__name__, phases)

    # end of function phase

# end of class AlignedSpliceSiteWithPhaseShiftGraph


class AlignedDonorSiteWithPhaseShiftGraph(AlignedSpliceSiteWithPhaseShiftGraph):
    """
    AlignedDonorSiteGraph (ADSG) class, inheriting from AlignedSpliceSiteWithPhaseShiftGraph class.
    """
    def __init__(self,**kwargs):
        """
        Initialize as a AlignedSpliceSiteWithPhaseShiftGraph 
        """
        AlignedSpliceSiteWithPhaseShiftGraph.__init__(self,**kwargs)
        # which type donors or acceptors?
        self._type = "Donor"

    # end of function __init__

# end of class AlignedDonorSiteWithPhaseShiftGraph 


class AlignedAcceptorSiteWithPhaseShiftGraph(AlignedSpliceSiteWithPhaseShiftGraph):
    """
    AlignedAcceptorSiteWithPhaseShiftGraph (AASG) class, inheriting from AlignedSpliceSiteWithPhaseShiftGraph class.
    """
    def __init__(self,**kwargs):
        """
        Initialize as a AlignedSpliceSiteWithPhaseShiftGraph 
        """
        AlignedSpliceSiteWithPhaseShiftGraph.__init__(self,**kwargs)
        # which type donors or acceptors?
        self._type = "Acceptor"

    # end of function __init__

# end of class AlignedAcceptorSiteWithPhaseShiftGraph 


################################################################################
### Helper function for AlignedPssmObjectGraph class and inheriting classes ####
################################################################################


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
    TCODE_MAX_NONCODING = 0.740
    TCODE_MIN_CODING    = 0.950
    return 1.0 + ( TCODE_MAX_NONCODING - tcode3p ) * 2.0 + ( tcode5p -  TCODE_MAX_NONCODING )

# end of function average_tcode_entropy_stopcodon


def average_tcode_entropy_startcodon(tcode5p=0.0,tcode3p=0.0):
    """
    Ratio between Coding (rigth) and Non-Coding (left) side of Start-Codon

    @type  tcode5p: positive float (0.450 - 1.320)
    @param tcode5p: average tcode score 5'/left/upstream of specific site

    @type  tcode3p: positive float (0.450 - 1.320)
    @param tcode3p: average tcode score 3'/right/downstream of specific site

    @rtype:  float (-0.450 - 2.160)
    @return: score attached to the likelihood that site is a startcodon; >1.0 is likely

    """
    TCODE_MAX_NONCODING = 0.740
    TCODE_MIN_CODING    = 0.950
    return 1.0 + ( TCODE_MAX_NONCODING - tcode5p ) * 2.0 + ( tcode3p -  TCODE_MAX_NONCODING )

# end of function average_tcode_entropy_startcodon


def nonebooleanintegermapper(noneboolean):
    """
    """
    mapper = {True: 1, None: 0, False: -1}
    try: return mapper[noneboolean]
    except KeyError: raise "nonebooleanintegermapper::KeyError('%s')" % noneboolean

# end of function nonebooleanintegermapper


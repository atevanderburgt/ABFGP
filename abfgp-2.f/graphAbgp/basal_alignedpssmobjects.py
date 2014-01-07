################################################################################
### BasalAlignedPssmObjectGraph class                                       ####
################################################################################

# graphAbgp imports
from graph_organism import OrganismGraph
from subclass_pssmobjects import BasalPSSMObjectGraphFunctions
from exceptions import *

# Python imports


class BasalAlignedPssmObjectGraph(OrganismGraph,BasalPSSMObjectGraphFunctions):
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
                'AlignedSpliceSiteoWithPhaseShiftGraph',
                'AlignedDonorSiteWithPhaseShiftGraph',
                'AlignedAcceptorSiteWithPhaseShiftGraph',
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

# end of class BasalAlignedPssmObjectGraph


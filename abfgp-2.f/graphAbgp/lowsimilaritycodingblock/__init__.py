################################################################################
### LowSimilarityRegionCodingBlockGraph class                               ####
################################################################################

# graphAbgp imports
from graphAbgp.graph_pacbpcollection import PacbpCollectionGraph
from graphAbgp.codingblock.genetree import GeneTreeOfCodingBlockFunctions
from graphAbgp.codingblock.sequenceretrieval import CodingBlockGraphSequenceRetievalFunctions
from graphAbgp.exceptions import *
import graphAbgp.ordering as ordering
import graphAbgp.conversion as conversion

# Abfgp Imports
from lib_cexpander import cexpanderanalyses

# Python imports
from sets import Set

# Import Global variables and settings
from settings.pacbp import ALTERNATIVE_ALIGNMENT_OVERLAP_RATIO


class LowSimilarityRegionCodingBlockGraph(PacbpCollectionGraph,GeneTreeOfCodingBlockFunctions,CodingBlockGraphSequenceRetievalFunctions):
    """
    """
    def __init__(self,alternative_alignment_overlap_ratio=ALTERNATIVE_ALIGNMENT_OVERLAP_RATIO):
        """
		Initialize a LowSimilarityRegionCodingBlockGraph
        """
        # Initialize from PacbpCollectionGraph, NOT from CodingBlockGraph
        PacbpCollectionGraph.__init__(self)
        self._node_object = {}          # in use for orf objects
        self.pacbps       = {}          # optional, not required
                                        # N.B. key of pacbps is different as
                                        # in CodingBlockGraph.
                                        # NOT: (pacb_identifier, node1, node2)
                                        # BUT: (node1, node2)

        # short name tag for printing
        self._short_name = "lsrCBG"

        # set special attributes needed
        self._splicedonorgraph    = None
        self._spliceacceptorgraph = None
        self._stopcodongraph      = None
        self._startcodongraph     = None

        # dicts for data of VISTA like tracks
        self._paoc_per_organism = {}    # attribute stays EMPTY!!
        self._pasc_per_organism = {}    # attribute stays EMPTY!!

        # attribute needed to set the GeneTreeGraph into
        self._GENETREE      = None

        # attribute for cexpanderanalyses data
        self._cexpander     = None

        # attributes that are needed once added into a genestructure
        self.IS_FIRST       = False
        self.IS_LAST        = False
        self.IS_IGNORED     = False
        self.IS_SPLITTED    = True      # this class is always IS_SPLITTED
        self.IS_5P_SPLITTED = True      # this class is always IS_SPLITTED
        self.IS_3P_SPLITTED = True      # this class is always IS_SPLITTED

        # set LowSimilarityRegion specific attributes
        self._omsr = {}
        self._forced_3p_ends = {}
        self._forced_5p_ends = {}

    # end of function __init__


    def __str__(self):
        """ """
        try:
            return "<%s wt=%s length:%s..%s [%s]>" % (
                self._short_name,
                self.total_weight(),
                min(self.overall_minimal_spanning_range_sizes().values()),
                max(self.overall_minimal_spanning_range_sizes().values()),
                " ".join( [ "%s(%s):%s..%s" % (node[0],node[1],min(self._omsr[node]),max(self._omsr[node])) for node in self.get_ordered_nodes() ] )
                )
        except:
            return "<%s wt=%s NO OMSR !? [%s]>" % (
                self._short_name,
                self.total_weight(),
                self.get_ordered_nodes() )

    # end of function __str__


    def _update_after_changes(self):
        """
        @attention: MUST BE OVERRIDDEN IN THE SUBCLASSES, if the returned value is not correct
        """
        self._update_node_object()

    # end of function _update_after_changes


    def node_by_organism(self,organism):
        """
        Get the node identifier belonging to the organism identifier.

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


    def set_node_omsr(self,node,omsr_range):
        """
        Set an OMSR for this node

        @type  node: *
        @param node: One Node identifier

        @type  omsr_range: Set()
        @param omsr_range: continious range (as Set) representing the OMSR
        """
        if node not in self.get_nodes():
            raise NodeNotPresentInGraph
        # set the omsr_range (actually a Set) to self._omsr
        self._omsr[node] = omsr_range

    # end of function set_node_omsr


    def set_organism_omsr(self,organism,omsr_range):
        """
        Set an OMSR for this organism

        @type  organism: *
        @param organism: organism identifier (presumably a string)

        @type  omsr_range: Set()
        @param omsr_range: continious range (as Set) representing the OMSR
        """
        if organism not in self.organism_set():
            raise OrganismNotPresentInGraph
        node = self.node_by_organism(organism)
        self.set_node_omsr(node,omsr_range)

    # end of function set_organism_omsr


    def add_node_and_object(self,node,object):
        """
        Create a single node with a accompagnying object. Object ISA Orf

        @type  node: *
        @param node: One Node identifier

        @type  object: *
        @param object: The Orf object represented by the Node identifier
        """
        if node not in self.get_nodes():
            self.nodes[node] = []
            self._node_object[node] = object

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

    # end of function _update_node_object


    def get_sequence_of_organism(self,organism):
        """ 
        Get sequence of this organism of this lsrCBG

        @type  organism: *
        @param organism: organism identifier (presumably a string)

		@rtype:  string
		@return: AA sequence string
        """
        if organism not in self.organism_set():
            raise OrganismNotPresentInGraph
        orf = self.get_organism_object(organism)
        omsr = self.overall_minimal_spanning_range(organism=organism)
        return orf.getaas(abs_pos_start=min(omsr),abs_pos_end=max(omsr)+1)

    # end of function get_sequence_of_organism


    def get_sequence_of_node(self,node):
        """ 
        Get sequence of this node of this lsrCBG

        @type  node: *
        @param node: One Node identifier

		@rtype:  string
		@return: AA sequence string
        """
        if node not in self.get_nodes():
            raise NodeNotPresentInGraph
        orf = self.get_node_object(node)
        omsr = self.overall_minimal_spanning_range(node=node)
        return orf.getaas(abs_pos_start=min(omsr),abs_pos_end=max(omsr)+1)

    # end of function get_sequence_of_node


    def have_all_starts(self):
        """
        Do all orfs in this lsrCBG have start codons?

        @rtype:  Boolean
        @return: True or False
        """
        for orf in self._node_object.values():
            if orf.has_start() == False:
                return False
        else:
            return True

    # end of function have_all_starts


    def have_all_starts_upstream_of_omsr(self,omsr_offset=10):
        """
        Are there (putative) startcodons (Methionines) upstream of the OMSR in all orfs?

        @type  omsr_offset: integer
        @param omsr_offset: offset into 3' side of OMSR to take into account as well

        @rtype:  Boolean
        @return: True or False
        """
        # first, do the quickest possible check: are there Methionines?
        if not self.have_all_starts(): return False

        omsr = self.overall_minimal_spanning_range()
        for node,orf in self._node_object.iteritems():
            offset = min(omsr[node]) + omsr_offset
            if not orf.has_start_upstream_of(offset):
                return False
        else:
            return True

    # end of function have_all_starts_upstream_of_omsr


    ####################################################################
    #### Functions for accessing the (Orf) objects                  ####
    ####################################################################

    def get_organism_object(self,organism):
        """
        Get the ONLY object of a specific organism from the graph

        @type  organism: *
        @param organism: organism identifier (presumably a string)

        @rtype:  *
        @return: Object that represents this organism
        """
        return self._node_object[self.node_by_organism(organism)]

    # end of function get_organism_object


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


    def get_organism_objects(self,organism,order_by=''):
        """
        Get all the OBJECTS of a specific organism from the graph

        @type  organism: *
        @param organism: organism identifier (presumably a string)

        @rtype:  list
        @return: list of all OBJECTS from this organism in this graph

        @attention: not recommended to use in this class, added for compatibily with other classes
        """
        return ordering.order_list_by_attribute(
                [ self._node_object[node] for node in self.get_organism_nodes(organism) ],
                order_by=order_by
                )

    # end of function get_organism_objects


    def get_orfs_of_graph(self,organism=None,node=None):
        """
        Get all the Orf objects of this graph

        @type  organism: *
        @param organism: Organism identifier (or None)

        @type  node: *
        @param node: Node identifier (or None)

		@rtype:  dictionary (or list if organism or node is specified)
		@return: dictionary with organisms (keys) and list of Orf objects (values),
                 or only list of Orf objects if an organism identifier was specified
        """
        orfs = {}
        for n in self.get_nodes():
            org = self._organism_from_node(n)
            orf = self.get_node_object(n)
            orfs[org] = [ orf ]
        # now check what to return
        if organism and organism not in self.organism_set():
            raise OrganismNotPresentInGraph
        elif node and node not in self.get_nodes():
            raise NodeNotPresentInGraph
        elif organism:
            return orfs[organism]
        elif node:
            return orfs[self._organism_from_node(node)]
        else:
            return orfs

    # end of function get_orfs_of_graph


    def set_genetree(self):
        """
        Convert LowSimilarityRegionCodingBlockGraph 2 GeneTree

        @attention: OVERWRITES codingblock_genetree.set_genetree

        @rtype:  GeneTreeGraph
        @return: GeneTreeGraph instance
        """
        self._GENETREE = conversion.LowSimilarityRegionCodingBlockGraph2GeneTreeGraph(self)
        return self._GENETREE

    # end of function genetree


    def distance_between_codingblocks(self,other,organism=None,node=None):
        """
        Distance in AA between two CodingBlockGraphs

        @type  organism: *
        @param organism: Organism identifier (or None)

        @type  node: *
        @param node: Node identifier (or None)

        @rtype:  dictionary (or integer if organism or node is specified)
        @return: dictionary with organisms (keys) and AA-distance between CBGs (values),
                 or only a distance if an Organism or Node identifier was specified

        @attention: other (CBG) is supposed to be 3p/rigth of 'self' (current CBG)
        """

        # get overall minimal spanning ranges
        omsrSelf  = self.overall_minimal_spanning_range()
        omsrOther = other.overall_minimal_spanning_range()
        distances = {}
        for org in self.organism_set().intersection(other.organism_set()):
            if organism and org != organism: continue
            nodeA = self.node_by_organism(org)
            nodeB = other.node_by_organism(org)
            distA = min( omsrSelf[nodeA] ) - max( omsrOther[nodeB] )
            distB = min( omsrOther[nodeB] ) - max( omsrSelf[nodeA] )
            distances[org] = max( [ distA, distB ] ) - 1

        # return distance only of a specific organism is requested for
        if organism:
            return distances[organism]
        else:
            return distances

    # end of function distance_between_codingblocks


    def potentially_contains_inframe_intron(self,min_inframe_intron_aa_size=10,max_inframe_introns=3,verbose=False):
        """
        Is this LowSimilarityRegion in fact an inframe intron in at least a single species?

        @type  min_inframe_intron_aa_size: integer
        @param min_inframe_intron_aa_size: minimal AA length of (inframe) intron

        @type  max_inframe_introns: integer
        @param max_inframe_introns: maximum number of inframe introns; the more, the unlikely!!

        @type  verbose: Boolean
        @param verbose: print status messages to STDOUT

        @rtype:  list
        @return: List of Organism indentifiers in which a potential inframe intron exists
        """
        omsr_sizes = self.overall_minimal_spanning_range_sizes()
        sizes = omsr_sizes.values()
        sizes.sort()

        if len(sizes) == 1:
            dif_sizes = max(sizes)
        else:
            dif_sizes = max(sizes) - min(sizes)

        # break out if there is no size difference -> nothing to improve!
        if dif_sizes == 0: return []

        # return if there is only a single one
        if len(sizes) == 1:
            if dif_sizes >= min_inframe_intron_aa_size:
                # a single organism remains that potentially is an inframe intron
                return list( self.organism_set() ) 
            else:
                # nope, not long enough to represent an inframe intron
                return []

        if verbose:
            print self
            print "0 inframe introns", sizes, "non-introns:", sizes, "introns:", [], "intron lengths:", []

        inframe_intron_count = 0
        for num_introns in range(1,max_inframe_introns+1):
            # check if list sizes is long enough
            if num_introns > len(sizes): break

            # split sizes in 2 lists of non-introns and (potential) inframe introns
            non_introns     = sizes[0:-num_introns]
            inframe_introns = sizes[-num_introns:]

            # calculate an average size for the non-introns
            if non_introns:
                average_size_non_introns = sum(non_introns) / len(non_introns) 
            else:
                # all orgs/orfs/genes in the lsrCBG are listed as introns, no non-introns left!
                average_size_non_introns = 0
            # correct the (potential) inframe-introns by this length
            inframe_introns_corrected_length = [ ifint - average_size_non_introns for ifint in inframe_introns ]
            
            print num_introns, "inframe introns", sizes, "non-introns:", non_introns, "introns:", inframe_introns, "intron lengths:", inframe_introns_corrected_length

            if not False in [ ifint >= min_inframe_intron_aa_size for ifint in inframe_introns_corrected_length ]:
                inframe_intron_count = num_introns
            else:
                # break out -> this is not a likely inframe intron anymore
                break

        if inframe_intron_count:
            # Gather the Organism identifiers in which an inframe intron
            # inproves this lsrCBG; take the `inframe_intron_count` longest
            # nodes/organisms in this lsrCBG
            sizes = omsr_sizes.values()
            sizes.sort()
            organisms_with_potential_inframe_intron = Set()
            for pos in range(-inframe_intron_count,0):
                size = sizes[pos]
                for node, omsr_size in omsr_sizes.iteritems():
                    if omsr_size == size:
                        organisms_with_potential_inframe_intron.add(
                                self._organism_from_node(node)
                                )
            if verbose: print list(organisms_with_potential_inframe_intron)
            return list(organisms_with_potential_inframe_intron)
        else:
            return []

    # end of function potentially_contains_inframe_intron


    ####################################################################
    #### Functions that simulate the outcome as expected in         ####
    #### a regular CodingBlockGraph for MSR, OMSR and MAXSR         ####
    ####################################################################


    def omsr_starts(self):
        """
        Return a dictionary of OMSR start coodinates (keys==nodes)

        @rtype:  dictionary 
        @return: dictionary with nodes (*, presumably tuple) as keys, OMSR start coords (integer) as values
        """
        retdict = {}
        for node,omsr in self.overall_minimal_spanning_range().iteritems():
            retdict[node] = min(omsr)
        return retdict

    # end of function omsr_starts


    def omsr_ends(self):
        """
        Return a dictionary of OMSR end coodinates (keys==nodes)

        @rtype:  dictionary
        @return: dictionary with nodes (*, presumably tuple) as keys, OMSR end coords (integer) as values
        """
        retdict = {}
        for node,omsr in self.overall_minimal_spanning_range().iteritems():
            retdict[node] = max(omsr)
        return retdict

    # end of function omsr_ends


    def has_overall_minimal_spanning_range(self):
        """
        Has this LowSimilarityRegionCodingBlockGraph an OMSR

		@rtype:  Boolean
		@return: True or False
        """
        if self._omsr:  return True
        else:           return False

    # end of function has_overall_minimal_spanning_range


    def overall_minimal_spanning_range(self,organism=None,node=None):
        """
        Get the OMSR of all nodes, or only of the node of the specified organism or node

        @type  organism: *
        @param organism: organism identifier (or None)

        @type  node: *
        @param node: Node identifier (or None)

		@rtype:  dictionary (or Set if organism is specified)
		@return: dictionary with nodes (keys) and spanning range Sets (values),
                 or only the spanning range Set if an organism identifier was specified
        """
        if organism and organism not in self.organism_set():
            raise OrganismNotPresentInGraph
        elif node and node not in self.get_nodes():
            raise NodeNotPresentInGraph
        elif node and self._omsr.has_key(node):
            return self._omsr[node]
        elif organism and self._omsr.has_key( self.node_by_organism(organism) ) :
            return self._omsr[self.node_by_organism(organism)]
        elif organism:
            return Set()
        elif node:
            return Set()
        else:
            return self._omsr

    # end of function overall_minimal_spanning_range


    def overall_minimal_spanning_range_sizes(self):
        """
        Get the sizes of the OMSR of all nodes

        @attention: only functional if pacbp's are present in the graph
                    (perform the function `harvest_pacbps_from_crossdata`)

        @attention: only meaningfull for (nearly) Fully Connected Graphs

		@rtype:  dictionary
		@return: dictionary with nodes (keys) and size of spanning range Sets (values)

        """
        d = {}
        for node, r in self.overall_minimal_spanning_range().iteritems():
            d[node] = len(r)
        return d

    # end of function overall_minimal_spanning_range_sizes


    def minimal_spanning_range(self,organism=None,node=None):
        """
        @attention: see function overall_minimal_spanning_range
        """
        return self.overall_minimal_spanning_range(organism=organism,node=node)

    # end of function minimal_spanning_range


    def maximal_spanning_range(self,organism=None,node=None):
        """
        @attention: see function overall_minimal_spanning_range
        """
        return self.overall_minimal_spanning_range(organism=organism,node=node)

    # end of function maximal_spanning_range


    def has_minimal_spanning_range(self):
        """
        @attention: see function has_overall_minimal_spanning_range
        """
        return self.has_overall_minimal_spanning_range()

    # end of function has_minimal_spanning_range


    def minimal_spanning_range_sizes(self):
        """
        @attention: see function overall_minimal_spanning_range_sizes
        """
        return self.overall_minimal_spanning_range_sizes()

    # end of function minimal_spanning_range_sizes


    def maximal_spanning_range_sizes(self):
        """
        @attention: see function overall_minimal_spanning_range_sizes
        """
        return self.overall_minimal_spanning_range_sizes()

    # end of function maximal_spanning_range_sizes


    def superinposed_left_maximal_spanning_range(self,organism=None):
        """
        lsrCBG have no/empty superinposed_maximal_spanning_range
        """
        smaxsr = {}
        for node in self._omsr.keys():
            smaxsr[node] = Set()
        if organism:
            return smaxsr[ self.node_by_organism(organism) ]
        else:
            return smaxsr

    # end of function superinposed_left_maximal_spanning_range


    def superinposed_rigth_maximal_spanning_range(self,organism=None):
        """
        @attention: see function superinposed_left_maximal_spanning_range
        """
        return self.superinposed_left_maximal_spanning_range()

    # end of function superinposed_rigth_maximal_spanning_range


    def has_superinposed_left_maximal_spanning_range(self,organism=None):
        """
        lsrCBG have no/empty superinposed_maximal_spanning_range
        """
        return False

    # end of function superinposed_left_maximal_spanning_range


    def has_superinposed_rigth_maximal_spanning_range(self,organism=None):
        """
        @attention: see function superinposed_left_maximal_spanning_range
        """
        return False

    # end of function has_superinposed_rigth_maximal_spanning_range



    def sources(self):
        """
        Get a dict with counts of sources of PacbP(ORF)s -> lsrPACBP 

        @rtype:  {}
        @return: dict with sources (keys) and occurence counts (values) of pacbp.source of PacbP(ORF)s

        @attention: see pacbp/__init__.py for recognized sources (blastp,unknown,hmmsearch,clustalw,clustalw-EXTENDED,lsrPACBP...)
        """
        return { 'lsrPACBP': self.edge_count() }

    # end of function sources


    def mutual_nodes(self,other):
        """ """
        return self.node_set().intersection(other.node_set())
    # end of function mutual_nodes

    def different_nodes(self,other):
        """ """
        return self.node_set().difference(other.node_set())
    # end of function different_nodes

    def organisms_with_different_orfs(self,other):
        """ """
        return [ self.organism_by_node(node) for node in self.different_nodes(other) ]
    # end of function organisms_with_different_orfs


    def cexpanderanalyses(self,**kwargs):
        """                 
        @attention: see lib_cexpander.cexpanderanalyses for function documentation
        """         
        self._cexpander = cexpanderanalyses(self,**kwargs)

    # end of function cexpanderanalyses

    ####################################################################
    #### Dummy Functions for compatibility with CodingBlockGraphs   ####
    ####################################################################

    def align_stop_codons(self):
        pass

    def extend_pacbporfs(self,*args):
        pass

    def update_edge_weights_by_minimal_spanning_range(self):
        pass

    def multiplealignment(self,**kwargs):
        pass

    def printmultiplealignment(self,**kwargs):
        pass

    def create_cache(self,**kwargs):
        pass

    def improvealignment(self,**kwargs):
        return False

    def correct_pacbpgaps_nearby_omsr(self,**kwargs):
        return False

# end of class LowSimilarityRegionCodingBlockGraph

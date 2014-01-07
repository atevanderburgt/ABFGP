################################################################################
#### GenestructureOfCodingBlockGraphs class                                 ####
################################################################################

# graphAbgp Imports
import graphAbgp 
from graph_organism import OrganismGraph
from genestructure_printing import GenestructureOfCodingBlockGraphsPrintingFunctions
from genestructure_cbginterface import CodingBlockGraphInterfaceFunctions
from genestructure_firstcbg import FirstCodingBlockGraphFunctions
from genestructure_finalcbg import FinalCodingBlockGraphFunctions
from genestructure_cbgoptimization import CbgInGeneStructureOptimizationFunctions
from genestructure_intermediatecbg import (
    IntermediateCodingBlockGraphFunctions,
    intermediateCBG_node_comparison,
    )
from genestructure_cbgremoval import CodingBlockGraphRemovalFunctions
from genestructure_intergenecity import IntergenecityFunctions
from genestructure_lsrcbgs import LowSimilarityRegionCodingBlockGraphFunctions
from subclass_genestructureorfmodel import GeneStructureOrfModelFunctions
from genestructure_cbgaddition import AddCodingBlockGraphFunctions, findmostlikelyCBG2GSGinsert 
from graph_lowsimilaritycodingblock import LowSimilarityRegionCodingBlockGraph
from graph_exoncollection import *
from lib_orfset import OrfSet
from exceptions import *
import ordering
import recombination
import conversion
import codingblock_splitting
import codingblock_operations

# Gene Imports
from gene.stop import StopCodon
from gene.codingblock import CodingBlockStart, CodingBlockEnd
from gene.exon import FirstExonOnOrf

# Abgp Imports
from intronprojection import cbgs_identical_pacbp_analysis
from codingblockgraphinterface import CodingBlockGraphInterface
from lib_shortleadingexononorf import find_leading_exon_on_orf
from lib_shorttailingexononorf import find_tailing_exon_on_orf
from lib_tinyexononorf import find_intermediary_codingblockgraph_with_tinyexon, find_tiny_exon_on_orf, bridge_two_pacbporfs_by_tinyexon
from lib_hmm import *
from lib_stopwatch import StopWatch
from lib_clustalw import (
    clustalw,
    strip_alignment_for_exterior_gaps,
    sprdif2clustalw2cbg
    )
import lib_cexpander

# Pacb class Import
from pacb import swap_query_and_sbjct
from pacb.conversion import pacbporf2pacbp, pacbp2pacbporf, pacbp_from_clustalw

# Python imports
from sets import Set
from copy import deepcopy
from time import time

# Global variables
from settings.codingblockgraph import *
from settings.genestructure import *
from settings.inframeintron import *
from settings.translationalstartsites import *

class GenestructureOfCodingBlockGraphs(
        OrganismGraph,
        GeneStructureOrfModelFunctions,
        AddCodingBlockGraphFunctions,
        CodingBlockGraphRemovalFunctions,
        GenestructureOfCodingBlockGraphsPrintingFunctions,
        CbgInGeneStructureOptimizationFunctions,
        CodingBlockGraphInterfaceFunctions,
        FirstCodingBlockGraphFunctions,
        FinalCodingBlockGraphFunctions,
        IntergenecityFunctions,
        IntermediateCodingBlockGraphFunctions,
        LowSimilarityRegionCodingBlockGraphFunctions
        ):
    """
    """
    def __init__(self,input):
        """
        """
        self.codingblockgraphs = []
        self.input = input
        # attributes for projected introns and status of it
        self._projected_donor_sites    = []
        self._projected_acceptor_sites = []
        self._HAS_INTRONS_PROJECTED    = False
        # attribute to store the codingblock gff data into
        self._gff_current_codingblocks = {}
        for organism in self.input.keys():
            self._gff_current_codingblocks[organism] = []
        # how many nodes are expected?
        self.EXACT_SG_NODE_COUNT = len(self.input)
        self.EXACT_SG_EDGE_COUNT = sum(range(0,self.EXACT_SG_NODE_COUNT))

        # attribute for storing current GeneTreeGraph
        self._GENETREE = None

        # thresholds for CBG acceptance into GSG: regular
        self.MAX_CBG_GTG_TOPO_DIF = MAX_CBG_GTG_TOPO_DIF
        self.MAX_CBG_GTG_ABS_DIF  = MAX_CBG_GTG_ABS_DIF
        self.MIN_CBG_GTG_ID_RATIO = MIN_CBG_GTG_ID_RATIO

        # thresholds for CBG aaceptance into GSG: first CBG sprdif
        self.MAX_CBG_FIRST_SPRDIF_GTG_TOPO_DIF = MAX_CBG_FIRST_SPRDIF_GTG_TOPO_DIF
        self.MAX_CBG_FIRST_SPRDIF_GTG_ABS_DIF  = MAX_CBG_FIRST_SPRDIF_GTG_ABS_DIF
        self.MIN_CBG_FIRST_SPRDIF_GTG_ID_RATIO = MIN_CBG_FIRST_SPRDIF_GTG_ID_RATIO

        # thresholds for CBG aaceptance into GSG: final CBG sprdif
        self.MAX_CBG_FINAL_SPRDIF_GTG_TOPO_DIF = MAX_CBG_FINAL_SPRDIF_GTG_TOPO_DIF
        self.MAX_CBG_FINAL_SPRDIF_GTG_ABS_DIF  = MAX_CBG_FINAL_SPRDIF_GTG_ABS_DIF
        self.MIN_CBG_FINAL_SPRDIF_GTG_ID_RATIO = MIN_CBG_FINAL_SPRDIF_GTG_ID_RATIO

        # thresholds for CBG acceptance into GSG: HMMsearch CBGs
        self.MAX_CBG_HMM_COMPLETION_GTG_TOPO_DIF = MAX_CBG_HMM_COMPLETION_GTG_TOPO_DIF 
        self.MAX_CBG_HMM_COMPLETION_GTG_ABS_DIF  = MAX_CBG_HMM_COMPLETION_GTG_ABS_DIF
        self.MIN_CBG_HMM_COMPLETION_GTG_ID_RATIO = MIN_CBG_HMM_COMPLETION_GTG_ID_RATIO

        # thresholds for CBG acceptance into GSG: large sprdif splitting 
        self.MAX_CBG_LARGE_SPRDIF_GTG_TOPO_DIF   = MAX_CBG_LARGE_SPRDIF_GTG_TOPO_DIF 
        self.MAX_CBG_LARGE_SPRDIF_GTG_ABS_DIF    = MAX_CBG_LARGE_SPRDIF_GTG_ABS_DIF
        self.MIN_CBG_LARGE_SPRDIF_GTG_ID_RATIO   = MIN_CBG_LARGE_SPRDIF_GTG_ID_RATIO 


    # end of function __init__

    def __str__(self):
        """ """
        return "<GSG: %s CBGs [AAs=%s tw=%s]>" % (
                len(self),
                self.total_length(),
                self.total_weight()
                )

    # end of function __str__

    ########################################################################
    ### Build-in functions to enshure behavious of GSG asa list of CBGs  ###
    ### at the moment, Read-only (no __setitem__ or __delitem__ )        ###
    ########################################################################

    def __len__(self):
        """ Return the length of list self.codingblockgraphs """
        return len(self.codingblockgraphs)
    # end of function __len__

    def __iter__(self):
        """ Return an iterable of self.codingblockgraphs """
        return iter(self.codingblockgraphs)
    # end of function __iter__

    def __getitem__(self,x):
        """ Return an xth element of self.codingblockgraphs """
        return self.codingblockgraphs[x]
    # end of function __getitem__

    ########################################################################
    ### Functions mapped from the CodingBlockGraphs                      ###
    ########################################################################

    def organism_by_node(self,node):
        """
        Get the organism identifier from the node identifier.

        @attention: overwrites OrganismGraph.organism_by_node()

        @type  node: *
        @param node: One Node

        @rtype:  *
        @return: organism identifier (presumably a string)
        """
        return node[0]

    # end of function organism_by_node


    def organism_set(self):
        return Set(self.input.keys())
    # end of function organism_set

    def organism_set_size(self):
        return len(self.input)
    # end of function organism_set_size

    ########################################################################
    ### Functions for GeneTreeGraph accessibility                        ###
    ########################################################################

    def genetree(self):
        """
        Return the GeneTreeGraph of this GenestructureGraph object

        @rtype:  GeneTreeGraph
        @return: GeneTreeGraph instance
        """
        if len(self):
            #if self._GENETREE:
            #    return self._GENETREE
            #else:
            scoreordered = ordering.order_graphlist_by_total_weight(self.codingblockgraphs)
            for cbg in scoreordered:
                if cbg.node_count() == self.EXACT_SG_NODE_COUNT:
                    return cbg.genetree()
            else:
                # not a single cbg in the GSG is a K(s) CBG, all K(s-x) !?
                return None
        else:
            return None

    # end of function genetree


    def get_genetree(self):
        """
        Return the GeneTreeGraph of this GenestructureGraph object

        @attention: alias for function genetree()

        @rtype:  GeneTreeGraph
        @return: GeneTreeGraph instance
        """
        return self.genetree()

    # end of function get_genetree


    def set_genetree(self):
        """
        Recreate the GeneTreeGraph of this GenestructureGraph object

        @rtype:  GeneTreeGraph
        @return: GeneTreeGraph instance
        """
        if len(self):
            scoreordered = ordering.order_graphlist_by_total_weight(self.codingblockgraphs)
            self._GENETREE = scoreordered[0].genetree()
            return self._GENETREE
        else:
            return None

    # end of function genetree


    def cbgpositerator(self,reversed=False):
        """
        Return an iterable list of CBG positions in the GSG
        """
        if reversed:
            return range(len(self)-1,-1,-1)

        else:
            return range(0,len(self))

    # end of function cbgpositerator


    def initialize_first_added_cbg(self):
        """
        TODO: TO BE DEPREACTED! nothing happens currently in this function...
        """
        return True


        # skip all this.....

        gtg = self.genetree()

        # attributes for threshold values for topological difference
        topo_dif_attributes = [
                'MAX_CBG_GTG_TOPO_DIF',
                'MAX_CBG_FINAL_SPRDIF_GTG_TOPO_DIF',
                ]

        # attributes for threshold values for absolute topological difference
        abs_dif_attributes = [
                'MAX_CBG_GTG_ABS_DIF',
                'MAX_CBG_FINAL_SPRDIF_GTG_ABS_DIF',
                ]

        # attributes for threshold values for identity ratio
        id_ratio_attributes = [
                'MIN_CBG_GTG_ID_RATIO',
                'MIN_CBG_FINAL_SPRDIF_GTG_ID_RATIO',
                ]

        for attr in topo_dif_attributes:
            newval = MAX_GTG_TOPO_DIF_FUNCTION(getattr(self,attr),gtg)
            setattr(self,attr,newval)


        for attr in abs_dif_attributes:
            newval = MAX_GTG_ABS_DIF_FUNCTION(getattr(self,attr),gtg)
            setattr(self,attr,newval)

        for attr in id_ratio_attributes:
            newval = MIN_GTG_ID_RATIO_FUNCTION(getattr(self,attr),gtg)
            setattr(self,attr,newval)

        #self.MAX_CBG_GTG_TOPO_DIF = MAX_GTG_TOPO_DIF_FUNCTION(self.MAX_CBG_GTG_TOPO_DIF,gtg)
        #self.MAX_CBG_GTG_ABS_DIF  = MAX_GTG_ABS_DIF_FUNCTION(self.MAX_CBG_GTG_ABS_DIF,gtg)
        #self.MIN_CBG_GTG_ID_RATIO = MIN_GTG_ID_RATIO_FUNCTION(self.MIN_CBG_GTG_ID_RATIO,gtg)


    # end of function initialize_first_added_cbg 

    ########################################################################
    ### Functions that summate CBG properties of the GeneStructure       ###
    ########################################################################

    def cbg_count(self):
        """
        Return the number of CBGs in the GSG

        @rtype:  integer
        @return: number of CBGs in the GSG
        """
        if len(self) == 0:
            return 0
        else:
            cnt = 0
            for cbg in self:
                if cbg.__class__.__name__ ==\
                'CodingBlockGraph':
                    cnt+=1
            return cnt

    # end of function cbg_count


    def lsrcbg_count(self):
        """
        Return the number of lsrCBGs in the GSG

        @rtype:  integer
        @return: number of lsrCBGs in the GSG 
        """
        if len(self) == 0:
            return 0
        else:
            cnt = 0
            for cbg in self:
                if cbg.__class__.__name__ ==\
                'LowSimilarityRegionCodingBlockGraph':
                    cnt+=1
            return cnt

    # end of function lsrcbg_count


    def total_weight(self):
        """
        Return the sum of the total_weight() of all CBGs in the GSG

        @rtype:  integer
        @return: summated total weight (==OMSR bitscore of edges) of all CBGs
        """
        return sum([cbg.total_weight() for cbg in self])

    # end of function total_weight


    def total_length(self):
        """
        Return the overall AA length of all CBGs in the GSG

        @rtype:  integer
        @return: summated AA length of all CBGs
        """
        if len(self)==0:
            return 0
        else:
            return min(self.overall_minimal_spanning_range_sizes().values())

    # end of function total_length

    ########################################################################
    ### Functions used in optimal informant selection                    ###
    ########################################################################

    def nearest_n_track(self,organism):
        """
        """
        if organism not in self.organism_set():
            raise OrganismNotPresentInGraph
        if self.input[organism]['genomeseq'].lower().find("n") == -1:
            return False
        # if here, then there are N-symbols in the sequence
        omsr = self.overall_minimal_spanning_range(organism)
        min_omsr = min(omsr)*3
        max_omsr = max(omsr)*3+3 
        dist_3p = self.input[organism]['genomeseq'].lower().find("n",max_omsr)
        if dist_3p > -1: dist_3p = dist_3p - max_omsr
        else:            dist_3p = False
        dist_5p = self.input[organism]['genomeseq'].lower()[0:min_omsr].rfind("n")
        if dist_5p > -1: dist_5p = min_omsr - dist_5p
        else:            dist_5p = False
        dist_intern = self.input[organism]['genomeseq'].lower()[min_omsr:max_omsr].find("n")
        if dist_intern > -1: dist_intern = 0
        else:                dist_intern = False

        distances = [ dist_3p, dist_5p, dist_intern ]
        if distances == [ False, False, False ]:
            return False
        else:
            while False in distances: distances.remove(False)
            return min(distances)

    # end of function nearest_n_track


    def orfidstructure(self,organism):
        """
        """
        if organism not in self.organism_set():
            raise OrganismNotPresentInGraph
        # loop over the CBGs in the GSG and grab orfids
        orfidlist = []
        for cbg in self:
            if organism in cbg.organism_set():
                orf = cbg.get_orfs_of_graph(organism=organism)[0]
                if not orfidlist:
                    orfidlist.append( orf.id )
                elif orfidlist[-1] != orf.id:
                    orfidlist.append( orf.id )
                else:
                    # lsrCBG or CBG with idenitcal Orf as previous one
                    pass
            else:
                # lsrCBG without this Organism identifier
                pass

        # return list of orf ids
        return orfidlist

    # end of function orfidstructure


    def is_organism_orf_cbginterface_optimal(self,*args):
        """
        Is the cbgIF between these 2 Orfs of this Organism identifier optimal?

        @attention: see get_organism_orf_cbginterface for argument documentation

        @rtype:  Boolean
        @return: True or False
        """
        cbgIF = self.get_organism_orf_cbginterface(*args)
        if not cbgIF:
            return False
        elif cbgIF.is_optimal():
            return True
        else:
            return False

    # end of function is_organism_orf_cbginterface_optimal


    def get_organism_orf_cbginterface_score(self,*args):
        """
        Get the score of the cbgIF between these 2 Orfs of this Organism identifier.

        @attention: see get_organism_orf_cbginterface for argument documentation

        @rtype:  integer
        @return: element of (2,1,0,-1)
                 2 optimal,
                 1 compatible and optimal donor or acceptor
                 0 compatible
                -1 incompatible
        """
        cbgIF = self.get_organism_orf_cbginterface(*args)
        if not cbgIF:
            return -1
        elif cbgIF.is_optimal():
            return 2
        elif cbgIF.is_compatible():
            if cbgIF.is_optimal_donor():
                return 1
            elif cbgIF.is_optimal_acceptor():
                return 1
            else:
                return 0
        else:
            return -1
        
    # end of function get_organism_orf_cbginterface_score


    def get_organism_orf_cbginterface(self,organism,orfAid,orfBid):
        """
        Get the cbgIF between these 2 Orfs of this Organism identifier.

        @type  organism: * (string)
        @param organism: Organism identifier

        @type  orfAid: integer
        @param orfAid: Orf integer identifier (>=0)

        @type  orfBid: integer
        @param orfBid: Orf integer identifier (>=0)

        @attention: apply Orf identifiers as the order in the GSG!

        @rtype:  CodingBlockGraphInterface
        @return: CodingBlockGraphInterface or None
        """
        if organism not in self.organism_set():
            raise OrganismNotPresentInGraph, organism

        # get orfIDs in the GSG for this Organism identifier
        orfidstruct = self.orfidstructure(organism)

        # check if Orfs /nodes are present in the GSG: (organism,orfid) == Node!
        if not orfAid in orfidstruct:
            raise NodeNotPresentInGraph, ( organism, orfAid )
        if not orfBid in orfidstruct:
            raise NodeNotPresentInGraph, ( organism, orfBid )

        # check if this orf interface exists (orfAid followed by orfBid)
        if len(self) <= 1:
            return None
        if orfidstruct.index(orfAid) != orfidstruct.index(orfBid) -1:
            return None

        # get the cbgIF that belongs to this - existing - Orf/Node interface
        cbgIF = None
        for b in range(1,len(self)):
            # check if this CBG has the requested orfBid 
            cbgB = self.codingblockgraphs[b]
            if organism not in cbgB.organism_set(): continue
            orfB = cbgB.get_orfs_of_graph(organism=organism)[0]
            if orfB.id != orfBid: continue

            # check if the previous CBG has the requested orfAid 
            cbgA = self.codingblockgraphs[b-1]            
            if organism not in cbgA.organism_set(): continue
            orfA = cbgA.get_orfs_of_graph(organism=organism)[0]
            if orfA.id != orfAid: continue

            # if here, we found the cbgIF
            cbgIF = cbgA._CBGinterface3p
            break

        # return the cbgIF
        return cbgIF

    # end of function get_organism_orf_cbginterface


    def is_organism_orf_cbginterface_continious_in_other_organism(self,organism,orfAid,orfBid,otherorganism):
        """
        Is the cbgIF between these 2 Orfs of this Organism continious in another Organism?

        @attention: see get_organism_orf_cbginterface for argument documentation

        @rtype:  Boolean
        @return: True or False
        """
        cbgIF = self.get_organism_orf_cbginterface(organism,orfAid,orfBid)
        if not cbgIF:
            return False
        if otherorganism not in cbgIF.donorCBG.organism_set():
            raise OrganismNotPresentInGraph, otherorganism
        if otherorganism not in cbgIF.acceptorCBG.organism_set():
            raise OrganismNotPresentInGraph, otherorganism
        orfA = cbgIF.donorCBG.get_orfs_of_graph(organism=otherorganism)[0]
        orfB = cbgIF.acceptorCBG.get_orfs_of_graph(organism=otherorganism)[0]
        if orfA.id == orfB.id:
            return True
        else:
            return False

    # end of function is_organism_orf_cbginterface_continious_in_other_organism


    def overall_aa_identityscore(self,**kwargs):
        """
        Get the overall AA identity score % of all the CBGs in the GSG

        @attention: see GeneTreeGraph.identity() for argument documentation

        @rtype:  float
        @return: aa identityscore ratio (0.0-1.0)
        """
        cbgs_aaident = []
        cbgs_lengths = []
        for cbg in self:
            if cbg.__class__.__name__ != 'CodingBlockGraph':
                continue
            cbgs_lengths.append( cbg.omsrlength() )
            cbgs_aaident.append( cbg.get_genetree().identity(**kwargs) )
            cbgs_aaident[-1] = cbgs_aaident[-1] * cbgs_lengths[-1]
        if not cbgs_aaident:
            return 0.0
        else:
            return sum(cbgs_aaident) / sum(cbgs_lengths)

    # end of function overall_aa_identity


    def overall_aa_identity(self,**kwargs):
        """
        Get the overall AA identity % of all the CBGs in the GSG

        @attention: see GeneTreeGraph.aaidentity() for argument documentation

        @rtype:  float
        @return: aa identity ratio (0.0-1.0)
        """
        cbgs_aaident = []
        cbgs_lengths = []
        for cbg in self:
            if cbg.__class__.__name__ != 'CodingBlockGraph':
                continue
            cbgs_lengths.append( cbg.omsrlength() )
            cbgs_aaident.append( cbg.get_genetree().aaidentity(**kwargs) )
            cbgs_aaident[-1] = cbgs_aaident[-1] * cbgs_lengths[-1]
        if not cbgs_aaident:
            return 0.0
        else:
            return sum(cbgs_aaident) / sum(cbgs_lengths)

    # end of function overall_aa_identity


    def overall_nt_identity(self,**kwargs):
        """
        Get the overall NT identity score % of all the CBGs in the GSG

        @attention: see GeneTreeGraph.ntidentity() for argument documentation

        @rtype:  float
        @return: nt identityscore ratio (0.0-1.0)
        """
        cbgs_ntident = []
        cbgs_lengths = []
        for cbg in self:
            if cbg.__class__.__name__ != 'CodingBlockGraph':
                continue
            cbgs_lengths.append( cbg.omsrlength() )
            cbgs_ntident.append( cbg.get_genetree().ntidentity(**kwargs) )
            cbgs_ntident[-1] = cbgs_ntident[-1] * cbgs_lengths[-1]
        if not cbgs_ntident:
            return 0.0
        else:
            return sum(cbgs_ntident) / sum(cbgs_lengths)

    # end of function overall_nt_identity

    ########################################################################
    ### Functions for MSR / OMSR / MAXSR of the GeneStructure            ###
    ########################################################################

    def minimal_spanning_range(self,organism):
        """
        Get the minimal spanning range in AA coordinates for a specific organism

        @type  organism: *
        @param organism: organism identifier (or None)

        @rtype:  Set
        @return: Set of all protein AA coordinates of all the codingblockgraphs

        @attention: alias for overall_minimal_spanning_range()
        """
        return self.overall_minimal_spanning_range(organism)

    # end of function minimal_spanning_range


    def overall_minimal_spanning_range(self,organism=None):
        """
        Get the overall minimal spanning range in AA coordinates for a specific organism

        @type  organism: *
        @param organism: organism identifier (or None)

        @rtype:  Set or dict
        @return: Set of all protein AA coordinates of all CBGs of an Organism or a dict
                 with Sets for each Organism identifier 
        """
        omsr = Set()
        allomsrs = {}
        for cbg in self:
            if organism:
                if organism in cbg.organism_set():
                    omsr.update( cbg.overall_minimal_spanning_range(organism=organism) )
            else:
                for node,thisomsr in cbg.overall_minimal_spanning_range().iteritems():
                    org = self.organism_by_node(node)
                    if allomsrs.has_key(org):
                        allomsrs[org].update(thisomsr)
                    else:
                        # create a new node-key in allomsrs dict.
                        # IMPORTANT! use deepcopy on `thisomsr`, because
                        # Python is call-by-reference, and otherwise the OMSR
                        # of the CBG that can be a cached dictionary, will be updated!
                        allomsrs[org] = deepcopy(thisomsr)
        if organism:
            return omsr
        else:
            return allomsrs

    # end of function overall_minimal_spanning_range


    def overall_minimal_spanning_range_sizes(self,organism=None):
        """
        Get the overall minimal spanning range size for a specific organism

        @type  organism: *
        @param organism: organism identifier (or None)

        @rtype:  int or dict
        @return: integer AA length of cumulative OMSR of the CBGs of an Organism or
                 a dict with integer values for each Organism identifier
        """
        if organism:
            return len( self.overall_minimal_spanning_range(organism=organism) )
        else:
            allomsrs = self.overall_minimal_spanning_range()
            for org in allomsrs.keys():
                allomsrs[org] = len( allomsrs[org] )
            return allomsrs

    # end of function overall_minimal_spanning_range_sizes



    def omsr2mask(self,organism,border_aa_offset=5,minimal_aa_length=10):
        """
        """
        coords = []
        for sg in self.codingblockgraphs:
            if sg.IS_IGNORED: continue
            if organism in sg.organism_set():
                omsr = sg.overall_minimal_spanning_range(organism=organism)
                if len(omsr) > border_aa_offset*2 + minimal_aa_length:
                    coordslice = ( min(omsr)+border_aa_offset, max(omsr)+1-border_aa_offset )
                    if coordslice not in coords:
                        coords.append( coordslice )

        # return the list of masking coordinates tuples
        return coords

    # end of function overall_minimal_spanning_range


    def complete_or_remove_small_graphs(self,verbose=False,
        minimal_ks_ksminx_omsr_ratio = 0.33,
        force_creation_of_first=True,  # not used yet!
        omit_cbgif_check=False,
        omit_first=False,
        omit_last=False,
        omit_exterior=False):
        """
        Complete all K(s-x) graphs with hmm searches or remove them when it failed
        """

        # counters for number of completed/deleted K(s-x) graphs
        cbgs_completed_cnt, cbgs_deleted_cnt = 0, 0

        ########################################################################
        if verbose: print [ cbg.node_count() for cbg in self ]
        ########################################################################
        minimal_observed_ksminxsize = min([ cbg.node_count() for cbg in self ])

        for ksminxsize in range(self.EXACT_SG_NODE_COUNT-1,minimal_observed_ksminxsize-1,-1):
            ksminx_cbgs = []
            has_eq_ksminx = Set([cbg.node_count() == ksminxsize for cbg in self])
            has_lt_ksminx = Set([cbg.node_count() < ksminxsize for cbg in self])
            if not True in has_eq_ksminx:
                continue
            if not True in has_lt_ksminx:
                pass
            else:
                # remove all the K(s-x) CBGs lt ksminxsize
                for i in self.cbgpositerator(reversed=True):
                    # if requested for, omit exterior CBGs
                    if (omit_exterior or omit_first) and i == 0:
                        continue
                    if (omit_exterior or omit_last) and i == len(self)-1:
                        continue 
                    if self.codingblockgraphs[i]._short_name == 'lsrCBG':
                        continue
                    if self.codingblockgraphs[i].node_count() < ksminxsize:
                        # temporarily remove this K(s-x) lt ksminxsize CBG
                        ksminx_cbgs.append( self.codingblockgraphs.pop(i) )
                    elif self.codingblockgraphs[i].node_count() == ksminxsize:
                        # do not remove, but place in ksminx_cbgs
                        ksminx_cbgs.append( self.codingblockgraphs[i] )
                    else:
                        pass  # CBG is K(s)CBG

            # now start completing the K(s-x) CBGs in the GSG
            ####################################################################
            if verbose:
                print "STARTING K(s-x) CBGs",
                print ksminxsize, self, len(ksminx_cbgs)
            ####################################################################
            ksminxgraphs_in_gsg = True
            while ksminxgraphs_in_gsg:
                # `pos`  is the variable name of current K(s-x) CBG to complete
                pos = None
                # find elegiable K(s-x) CBG pos to complete
                for testpos in self.cbgpositerator(reversed=True):
                    # if requested for, omit exterior CBGs
                    if (omit_exterior or omit_first) and testpos == 0:
                        continue
                    if (omit_exterior or omit_last)  and testpos == len(self)-1:
                        continue 
                    # get this CBG from the ordered list
                    this = self.codingblockgraphs[testpos]
                    # ignore lsrCBGs and CBGs that are already Ks graphs
                    if this.node_count() == self.EXACT_SG_NODE_COUNT: continue
                    if this.__class__.__name__ == 'LowSimilarityRegionCodingBlockGraph':
                        continue

                    # get the neighbouring CBGS
                    next,prev = None,None
                    if testpos > 0:
                        prev = self.codingblockgraphs[testpos-1]
                    if testpos < len(self)-1:
                        next = self.codingblockgraphs[testpos+1]
                    # if prev isa K(s) CBG, then start optimizing this K(s-x) CBG
                    if prev and prev.node_count() == self.EXACT_SG_NODE_COUNT:
                        # start optimizing this K(s-x) CBG!
                        pos = testpos
                        break
                    # if next isa K(s) CBG, then start optimizing this K(s-x) CBG
                    if next and next.node_count() == self.EXACT_SG_NODE_COUNT:
                        # start optimizing this K(s-x) CBG!
                        pos = testpos
                        break
    
                else:
                    # if end of for loop is reached: no K(s-x) graphs
                    # are remaining -> break while loop!
                    ksminxgraphs_in_gsg = False
                    break
    
                ################################################################
                # if this point is reached, we start to complete this K(s-x) CBG
                ################################################################
                # get this, prev and next cbg
                this = self.codingblockgraphs[pos]
                if pos > 0:
                    prev = self.codingblockgraphs[pos-1]
                if pos < len(self)-1:
                    next = self.codingblockgraphs[pos+1]
    
                ####################################################################
                if verbose:
                    print "\n", pos, "K(s-x) CBG completion", self
                    print prev
                    print this
                    print next
                ####################################################################
                    
                # complete with cbghmmsearch2pacbpcollection
                pacbpCollection = cbghmmsearch2pacbpcollection(this,self.input,
                        next=next,prev=prev,
                        hmmsearch_num_hits=3,     # normal == 3
                        max_intron_nt_length=250, # IMPORTANT!!! keep this number rather low.
                                                  # when a more biological realistic maximal length is choosen
                                                  # e.g. >=500nt, this HMMcompleted cbg is likely to get abandoned
                                                  # when by chance a wrong orf is linked to the HMMcbg
                        )
    
                #################################################################### 
                if verbose:
                    node_set = pacbpCollection.node_set()
                    print "elegiable HMM nodes:", node_set.difference(this.node_set())
                ####################################################################
    
                # make splitted subgraphs from PacbpCollection
                # cbgs must have self.EXACT_SG_NODE_COUNT nodes
                # and no missing edges are allowed
                # to the number of nodes missing in the input splittedCBG)
                dpcPacbpCollection = deepcopy(pacbpCollection)
                splitted_subgraphs = pacbpCollection.find_fully_connected_subgraphs(
                            edges=self.EXACT_SG_NODE_COUNT-1 ,
                            max_missing_edges=0
                            )
    
                # get pacbps for the splitted subgraphs and update edge weights
                completed_subgraphs = []
                for spl in splitted_subgraphs:
                    if spl.node_count() != self.EXACT_SG_NODE_COUNT: continue
                    # harvest pacbps from the deepcopied PacbpCollection
                    spl.harvest_pacbps_from_pacbpcollection(dpcPacbpCollection)
                    if not spl.has_overall_minimal_spanning_range(): continue
                    spl.update_edge_weights_by_minimal_spanning_range()
                    # create cexpander data in CBG object
                    spl.cexpanderanalyses()
                    # check the cexpander binarystring for only zeros
                    if spl._cexpander.binarystring.count("1") == 0:
                        # no uniformly matched AAs in this CBG -> ignore by continue
                        continue
                    # and optimize this novel cbg on the cexpander binarystring
                    try:
                        status = lib_cexpander.cexpander_checkCBG4omsrbordergaps(spl)
                    except lib_cexpander.ZeroUniformlyAlignedPositions:
                        # multiple alignment crashed after 5' optimization step
                        # no significant alignment remains -> ignore by continue
                        continue
                    except NoOverallMinimalSpanningRange:
                        # cexpander showed large non-uniform AA regions and after
                        # correction no OMSR left -> ignore by continue 
                        continue 
                    # check if the K(s) and K(s-x) CBG omsrlength do not differ to much
                    if float(spl.omsrlength()) / float(this.omsrlength()) < minimal_ks_ksminx_omsr_ratio:
                        #############################################################
                        if verbose:
                            print "OMITTED ks/ksminx OMSR ratio",
                            print "%1.3f" % ( float(spl.omsrlength()) / float(this.omsrlength()) ),
                            print "< %1.3f" %  minimal_ks_ksminx_omsr_ratio
                        #############################################################
                        continue
                    # if here, then add to completed_subgraphs list 
                    completed_subgraphs.append(spl)
    
                # check if there are hmm completed solutions
                if not completed_subgraphs:
                    # this CBG cannot be completed with HMM at all.
                    # remove this K(s-x) codingblockgraph
                    this = self.codingblockgraphs.pop(pos)
                    cbgs_deleted_cnt += 1
                    # whipe out cbgInterface objects here
                    if pos < len(self):
                        self.codingblockgraphs[pos]._CBGinterface5p = None
                    if pos > 0:
                        self.codingblockgraphs[pos-1]._CBGinterface3p = None
                    # continue to the postion in the GSG 
                    continue
    
                ####################################################################
                if verbose or True:
                    print "elegiable HMM completed K(s-x) graphs (%s):" % len(completed_subgraphs)
                    for completedcbg in completed_subgraphs: print completedcbg
                ####################################################################

                if prev and next:
                    partGSG = GenestructureOfCodingBlockGraphs(self.input)
                    partGSG.codingblockgraphs = [prev,next]
                elif prev:
                    partGSG = GenestructureOfCodingBlockGraphs(self.input)
                    partGSG.codingblockgraphs = [prev]
                elif next:
                    partGSG = GenestructureOfCodingBlockGraphs(self.input)
                    partGSG.codingblockgraphs = [next]
                else:
                    pass # cannot happen
   
                if omit_cbgif_check:
                    # just place the highest scoring (top listed) completed CBG
                    # into the GSG. THIS IS VERY DANGEROUS becasue FP CBGs are
                    # easily created here. Use the omit_cbgif_check parameter with care...
                    if prev and next:
                        partGSG.codingblockgraphs.insert(1,completed_subgraphs[0])
                    elif prev:
                        partGSG.codingblockgraphs.append(completed_subgraphs[0])
                    elif next:
                        partGSG.codingblockgraphs.insert(0,completed_subgraphs[0])
                    else:
                        pass
                    # create the CBGifs
                    partGSG.create_cbginterfaces()
                    # update the slice into the main GSG
                    stapos, endpos = pos+0, pos+1
                    if prev: stapos-=1
                    if next: endpos+=1
                    self.codingblockgraphs.__setslice__(stapos,endpos,partGSG.codingblockgraphs)
                    cbgs_completed_cnt += 1
                    ####################################################################
                    if verbose:
                        print "K(s-x) CBG completed and added, omit_cbgif_check:"
                        print completed_subgraphs[0]
                    ####################################################################
                    # continue here!
                    continue
 
                # copy genetree object from main GSG
                partGSG._GENETREE = self._GENETREE
                curpartgsglen = len(partGSG)
                partGSG = findmostlikelyCBG2GSGinsert(partGSG,completed_subgraphs,verbose=verbose)
                if verbose:
                    print partGSG, "graphs remaining:", len(completed_subgraphs)
                if len(partGSG) > curpartgsglen:
                    # update the slice into the main GSG
                    stapos, endpos = pos+0, pos+1
                    if prev: stapos-=1
                    if next: endpos+=1
                    self.codingblockgraphs.__setslice__(stapos,endpos,partGSG.codingblockgraphs)
                    cbgs_completed_cnt += 1
                else:
                    # remove this K(s-x) codingblockgraph
                    this = self.codingblockgraphs.pop(pos)
                    cbgs_deleted_cnt += 1
                    # whipe out cbgInterface objects here
                    if pos < len(self):
                        self.codingblockgraphs[pos]._CBGinterface5p = None
                    if pos > 0:
                        self.codingblockgraphs[pos-1]._CBGinterface3p = None
                    # continue to the postion in the GSG 
                    continue

            # Here, the `while ksminxgraphs_in_gsg:` loop is ended
            ####################################################################
            if verbose:
                print "ENDING K(s-x) CBGs",
                print ksminxsize, self, len(ksminx_cbgs)
            ####################################################################
            # replace all ksminx_cbgs in the GSG and go to the next
            if ksminxsize > 3:
                # place back all the ksminx_cbgs in the GSG;
                # omit_conditional_addition because these K(s-x) came from
                # the original GSG!
                for ksmincbg in ksminx_cbgs:
                    self.add_codingblock(ksmincbg,omit_conditional_addition=True)
                ################################################################
                if verbose:
                    print "replaced previously removed K(s-x) CBGs",
                    print ksminxsize, self, len(ksminx_cbgs)
                ################################################################
            else:
                # ksminxsize == 3 is the smallest K(s-x) CBG that can exist
                # do not place back the previously abandoned K(s-x) CBGs
                # stored in list ksminx_cbgs
                pass


        # reset FIRST & LAST cbgs; things could have been messed up ..
        self.finalize_genestructure()

        # return number of how much K(s-x) graphs completed/deleted
        return (cbgs_completed_cnt, cbgs_deleted_cnt)

    # end of function complete_or_remove_small_graphs 


    def remove_reversecomplement_cbgs(self,overlap_ratio = None,
        ignore_ksminx_cbgs=False, ignore_ks_cbgs=False,
        verbose=False):
        """
        Remove all CBGs that are ReversecomplementCodingBlockGraph

        @attention: function removes as well all non-CodingBlockGraph objects

        @rtype:  ( list, list )
        @return: [ deleted CBGs ], [ ReversecomplementCodingBlockGraphs ]
        """
        deleted = [] 
        reversedcbgs = []
        for i in range(len(self)-1,-1,-1):
            cbg = self[i]
            if ignore_ks_cbgs and cbg.node_count() == self.EXACT_SG_NODE_COUNT:
                continue
            if ignore_ksminx_cbgs and cbg.node_count() < self.EXACT_SG_NODE_COUNT:
                continue
            if cbg.__class__.__name__ == 'LowSimilarityRegionCodingBlockGraph':
                continue
            IS_DELETED = False 
            for _revcbg in reversedcbgs:
                if cbg.is_reversecomplement(revcbg=_revcbg):
                    deleted.append( self.codingblockgraphs.pop(i) )
                    #############################################################
                    if verbose: print "d", cbg, "\n", _revcbg
                    #############################################################
                    IS_DELETED = True
                    break
            # continue if cbg is deleted
            if IS_DELETED: continue
            # obtain reverse complement CBG
            revcbg = cbg.get_reversecomplement()
            if cbg.is_reversecomplement(revcbg=revcbg):
                reversedcbgs.append( revcbg )
                deleted.append( self.codingblockgraphs.pop(i) )
                #############################################################
                if verbose: print "D", cbg, "\n", revcbg
                #############################################################
            else:
                #############################################################
                if verbose: print "A", cbg, "\n", revcbg
                #############################################################
                pass

        # return list with (unique) reversedcbgs & list of deleted CBGs
        return deleted, reversedcbgs
        
    # end of function remove_reversecomplement_cbgs


    def replace_frameshifted_cbgs(self,
        cbg_max_aa_length = 20,
        verbose=False):
        """ """
        frameshifted = 0
        removed      = 0 
        for pos in self.cbgpositerator(reversed=True):
            cbg = self.codingblockgraphs[pos]
            if cbg.__class__.__name__ == 'LowSimilarityRegionCodingBlockGraph':
                continue
            if cbg.omsrlength() > cbg_max_aa_length:
                continue
            prevCBG, nextCBG = None, None
            if pos > 0:
                prevCBG = self.codingblockgraphs[pos-1]
                if prevCBG.__class__.__name__ ==\
                'LowSimilarityRegionCodingBlockGraph':
                    continue
            if pos < len(self)-1:
                nextCBG = self.codingblockgraphs[pos+1]
                if nextCBG.__class__.__name__ ==\
                'LowSimilarityRegionCodingBlockGraph':
                    continue
 
            ####################################################################
            if verbose: print "FScheck:", cbg
            ####################################################################
            fscbg = cbg.get_frameshifted_cbg(self.input,verbose=verbose)
            if fscbg:
                if prevCBG and not fscbg.node_set().symmetric_difference(prevCBG.node_set()):
                    # fscbg ISA cbg that is already in the GSG!
                    removedCBG = self.remove_cbg_by_pos(pos)
                    removed+=1
                    continue
                if nextCBG and not fscbg.node_set().symmetric_difference(nextCBG.node_set()):
                    # fscbg ISA cbg that is already in the GSG!
                    removedCBG = self.remove_cbg_by_pos(pos)
                    removed+=1
                    continue

                # if here, then replace cbg by fscbg
                fscbg.IS_FIRST = cbg.IS_FIRST
                fscbg.IS_LAST  = cbg.IS_LAST
                self.codingblockgraphs[pos] = fscbg
                if pos > 0:
                    self.codingblockgraphs[pos-1]._CBGinterface3p = None
                if pos < len(self)-1: 
                    self.codingblockgraphs[pos+1]._CBGinterface5p = None
                frameshifted += 1

        # correct for incorrect lsrCBGs after CBH changes
        if frameshifted or removed: self.remove_incorrect_lsrcbgs()

        # return counter how much CBGs are frameshifted & removed
        return ( frameshifted, removed )

    # end of function replace_frameshifted_cbgs


    def removeweakconnectednode(self,weakorgnode):
        """
        Remove a (weakly) connected organism/gene from all the CBGs in the GSG 

        @type  weakorgnode: * 
        @param weakorgnode: organism identifier (in practice a string) 
        """
        # try if node still present in self.input
        # python is call-by-reference; it should/can be removed already!
        if self.input.has_key(weakorgnode): del( self.input[weakorgnode] )
        # update numbers of desired nodes and edges
        self.EXACT_SG_NODE_COUNT = len(self.input)
        self.EXACT_SG_EDGE_COUNT = sum(range(0,self.EXACT_SG_NODE_COUNT))

        # remove nodes from CBGs in GSG (self)
        for cbg in self:
           if weakorgnode in cbg.organism_set():
               node = cbg.node_by_organism(weakorgnode)
               cbg.del_node(node)
               if cbg.__class__.__name__ == 'LowSimilarityRegionCodingBlockGraph': continue
               # clear and recreate cache!
               cbg.clear_cache()
               cbg.update_edge_weights_by_minimal_spanning_range()
               cbg.create_cache()
        # check if all CBGs have at least 2 nodes; when an incomplete graph of only 2 nodes
        # had on of its nodes removed, it is no longer a CBG!
        for pos in range(len(self)-1,-1,-1):
            if self[pos].node_count() == 1:
                # remove it!
                self.pop(pos)
        # check for doublets in CBGs -> occurs because intron present in removed organism
        for pos in range(len(self)-1,0,-1):
            this = self[pos]
            prev = self[pos-1]
            if not this.node_set().symmetric_difference(prev.get_nodes()):
                if str(this) == str(prev):
                    popped = self.codingblockgraphs.pop(pos)
        
        # finally, try to extend some of the CBGs of removal of the node
        self.extend_cbgs()
        # and update the genetree
        self.set_genetree()

    # end of function removeweakconnectednode

    ########################################################################
    ### Functions for splitting CBGs in GSG on spanningrange difference  ###
    ########################################################################

    def split_cbgs_on_lowconnected_spanningrange_difference(self,side=None,
        sprdif_min_aa_length=CBG_LARGE_RIGTH_SPRDIF_MIN_AA_LENGTH,
        cbg_min_node_count=None,
        verbose=False):
        """
        Split CBGs on lowconnected spanningrangedifference covered by all nodes and try to assign a new CBG in the splitted area

        @type  side: string
        @param side: 'left' or 'rigth', meaning a split on left/5p or rigth/3p side

        @type  sprdif_min_aa_length: integer
        @param sprdif_min_aa_length: minimal length of the sprdif in aa's

        @type  cbg_min_node_count: integer
        @param cbg_min_node_count: minimal number of nodes in a CBG to be elegiable for trying a split
        """
        # input integrity check
        if side not in ('left','rigth'):
            message = "`side` must be 'left' or 'rigth', not '%s'" % side
            raise InproperlyAppliedArgument, message

        # loop backwards over the genestructure in case of an insert
        for pos in range(len(self)-1,-1,-1):
            cbg = self[pos]
            # ignore lsrCBGs
            if cbg.__class__.__name__ == 'LowSimilarityRegionCodingBlockGraph': continue
            # ignore when less that a specific node count
            if cbg_min_node_count and cbg.node_count() < cbg_min_node_count: continue
            # ignore CBGs with a 3p/5p split (in case the side is rigth/left)
            if side == 'left'  and cbg.IS_5P_SPLITTED: continue
            if side == 'rigth' and cbg.IS_3P_SPLITTED: continue

            # obtain the spanningrange difference
            sprdif = codingblock_splitting.spanningrange_difference(cbg,side,
                    correct_sprdif_all_nodes=False,
                    sprdif_min_aa_length=sprdif_min_aa_length,
                    sprdif_min_node_count=cbg.organism_set_size(),
                    )
            # if no sprdif of the requested properties -> continue
            if not sprdif: continue
         
            # gather actual sequences of the sprdif area and store the coords
            # of these sequences for later use
            seqs = {}
            node2coords = {}
            for node,vlist in sprdif.iteritems():
                (org,orfid) = node
                orf = cbg.get_orfs_of_graph(organism=org)[0]
                seq = orf.getaas(abs_pos_start=min(vlist),abs_pos_end=max(vlist)+1)
                hdr = cbg.stringrepr(node)
                seqs[hdr] = seq
                node2coords[node] = ( min(vlist) , max(vlist)+1 )

            # some printing in verbose mode
            if verbose:
                print "CBG lowconnected sprdif", side, cbg
                for k,vlist in sprdif.iteritems(): print k,len(vlist),min(vlist),max(vlist),
                print ""
                for hdr,seq in seqs.iteritems():
                    print hdr,"\t",seq

            # make a new CBG by doing pairwise clustalw2pacbp formation
            newCBG = graphAbgp.CodingBlockGraph()
            newCBG.add_nodes(cbg.get_nodes())
            for hdr1,hdr2 in codingblock_splitting.pairwise(seqs.keys()):
                # get node from stringrepr and get Orfs by nodes
                node1 = cbg.stringrepr2node(hdr1)
                node2 = cbg.stringrepr2node(hdr2)
                orfQ  = cbg.get_orfs_of_graph(organism=node1[0])[0]
                orfS  = cbg.get_orfs_of_graph(organism=node2[0])[0]
                # make pacbp with clustalw
                clwseqs = { hdr1: seqs[hdr1], hdr2: seqs[hdr2] }
                aligned_seqs, _alignment = clustalw(seqs=clwseqs)
                pacbp = pacbp_from_clustalw(
                        alignment=(aligned_seqs[hdr1],_alignment,aligned_seqs[hdr2]),
                        coords=(node2coords[node1][0],node2coords[node1][1],node2coords[node2][0],node2coords[node2][1])
                        )
                if pacbp:
                    # make extended pacbporf
                    pacbporf = pacbp2pacbporf(pacbp,orfQ,orfS)
                    pacbporf.extend_pacbporf_after_stops()
                    # construct pacbp unique key and get wt as edge weight
                    pkey = pacbporf.construct_unique_key(node1,node2)
                    wt = pkey[0]
                    newCBG.add_edge(node1,node2,wt=wt)
                    newCBG.pacbps[(pkey,node1,node2)] = pacbporf
                    if verbose: print node1, node2, pacbp

            # check if we have a fully connected CBG (case 1) with omsr (case 2) of proper length (case 3)
            if newCBG.connectivitysaturation() != 1.0 or\
            not newCBG.has_overall_minimal_spanning_range() or\
            newCBG.omsrlength() < sprdif_min_aa_length:
                ################################################################
                if verbose:
                    print  "IGNORED: complete",  newCBG.connectivitysaturation(),
                    print "omsr", newCBG.has_overall_minimal_spanning_range(),
                    print "omsrlength", newCBG.omsrlength()
                ################################################################
                # no, not a valid newCBG -> continue with the next one
                continue

            # if here, start adding this newCBG to the GSG
            newCBG.update_edge_weights_by_minimal_spanning_range()

            ################################################################
            if verbose:
                print "NEWCBG:", newCBG
                newCBG.printmultiplealignment()
            ################################################################

            topo_dif = self.genetree().graphalignmentdifference( newCBG.genetree() )
            abs_dif  = self.genetree().absolutegraphalignmentdifference( newCBG.genetree() )
            id_ratio = newCBG.genetree().identity() / self.genetree().identity()
            tcode    = newCBG.msr_tcode_score()

            ################################################################
            if verbose:
                print "omsrlength:", newCBG.omsrlength(),
                print "TCODE: %1.2f ratio: %1.3f topodif: %1.3f absdif %1.3f" % (
                    tcode,id_ratio,topo_dif,abs_dif)
            ################################################################


            status = self.add_codingblock(newCBG,only_try_adding=True,
                    max_cbg_gtg_topo_dif=MAX_CBG_LOWCONNECTED_SPRDIF_GTG_TOPO_DIF,
                    max_cbg_gtg_abs_dif=MAX_CBG_LOWCONNECTED_SPRDIF_GTG_ABS_DIF,
                    min_cbg_gtg_id_ratio=MIN_CBG_LOWCONNECTED_SPRDIF_GTG_ID_RATIO,
                    min_tcode_omsr=MIN_CBG_LOWCONNECTED_TCODE_OMSR,
                    )

            if status:
                # newCBG can be placed in the GSG and forfills the requirements!
                if side == 'left':
                    lsrCBG = codingblock_splitting.create_intermediate_lowsimilarity_region( newCBG, cbg )
                    if not lsrCBG.node_count(): lsrCBG = None
                    cbg.IS_SPLITTED = True
                    cbg.IS_5P_SPLITTED = True
                    if cbg.IS_FIRST: newCBG.IS_FIRST = True
                    cbg.IS_FIRST = False
                    newCBG.IS_SPLITTED = True
                    newCBG.IS_3P_SPLITTED = True
                else:
                    lsrCBG = codingblock_splitting.create_intermediate_lowsimilarity_region( cbg, newCBG )
                    if not lsrCBG.node_count(): lsrCBG = None
                    cbg.IS_SPLITTED = True
                    cbg.IS_3P_SPLITTED = True
                    if cbg.IS_LAST: newCBG.IS_LAST = True
                    cbg.IS_LAST = False
                    newCBG.IS_SPLITTED = True
                    newCBG.IS_5P_SPLITTED = True

                # actually add the newCBG and possibly a lsrCBG
                status = self.add_codingblock(newCBG,
                        max_cbg_gtg_topo_dif=MAX_CBG_LOWCONNECTED_SPRDIF_GTG_TOPO_DIF,
                        max_cbg_gtg_abs_dif=MAX_CBG_LOWCONNECTED_SPRDIF_GTG_ABS_DIF,
                        min_cbg_gtg_id_ratio=MIN_CBG_LOWCONNECTED_SPRDIF_GTG_ID_RATIO,
                        min_tcode_omsr=MIN_CBG_LOWCONNECTED_TCODE_OMSR,
                        )
                if status and lsrCBG:
                    status = self.add_codingblock(lsrCBG,log_debug=False)
                    ################################################################
                    if verbose:
                        print "ADDED", newCBG
                        print "ADDED:", status, lsrCBG
                    ################################################################


    # end of function split_cbgs_on_lowconnected_spanningrange_difference


    def split_cbgs_on_spanningrange_difference(self,side=None,
        sprdif_min_aa_length=CBG_MIN_AA_LENGTH,
        sprdif_min_node_count=CBG_SPRDIF_MIN_NODE_COUNT,
        sprdif_min_identity=0.0,     # 0.0 means by default no check on this property!
        cbg_min_identity=0.0,        # 0.0 means by default no check on this property!
        cbg_min_node_count=None,verbose=False,
        do_cexander_allzeros_check = True,
        ignore_first_cbg=False,
        ignore_final_cbg=False,
        perform_cbgif_optimization=False,              # new 01/12/2009 !!
        remove_lower_scoring_overlapping_cbgs=False,
        sprdif_vs_omsr_min_identityratio=CBG_SPRDIF_VS_OMSR_MIN_IDENTITYSCORE_RATIO):
        """             
        Split CBGs on spanningrangedifference and try to assign a new CBG in the splitted area
            
        @type  side: string 
        @param side: 'left' or 'rigth', meaning a split on left/5p or rigth/3p side
                
        @type  sprdif_min_aa_length: integer
        @param sprdif_min_aa_length: minimal length of the sprdif in aa's
            
        @type  sprdif_min_node_count: integer
        @param sprdif_min_node_count: minimal number of nodes that support the sprdif

        @type  sprdif_min_identity: float
        @param sprdif_min_identity: minimal identity (0.0-1.0) of the splitted CBG to be accepted/tried

        @type  cbg_min_identity: float 
        @param cbg_min_identity: minimal identity (0.0-1.0) of the CBG to allow a split on sprdif

        @type  cbg_min_node_count: integer
        @param cbg_min_node_count: minimal number of nodes in a CBG to be elegiable for trying a split

        @type  do_cexander_allzeros_check: Boolean
        @param do_cexander_allzeros_check: omit novel CBGs that have not a single uniformly alignable AA position

        @type  ignore_first_cbg: Boolean
        @param ignore_first cbg: do or ignore the first CBG in the list / the CBG labelled as IS_FIRST

        @type  ignore_final_cbg: Boolean
        @param ignore_final cbg: do or ignore the final CBG in the list / the CBG labelled as IS_LAST
        """ 
        # input integrity check
        if side not in ('left','rigth'):
            message = "`side` must be 'left' or 'rigth', not '%s'" % side
            raise InproperlyAppliedArgument, message
        if sprdif_min_node_count < 2:
            message = "`sprdif_min_node_count` must be >=2, not '%s'" % sprdif_min_node_count
            raise InproperlyAppliedArgument, message

        stw = StopWatch(name='StwTestSprdif(MIN_AA_LENGT:%s)' % sprdif_min_aa_length)
        stw.start()

        # counters for how much splits have been tried/done
        sprdif_splits_done_cnt  = 0
        sprdif_splits_tried_cnt = 0

        # loop backwards in case a splitted CBG is inserted!
        for pos in self.cbgpositerator(reversed=True):
            cbg = self.codingblockgraphs[pos]
            # skip IGNORED and LowSimilarityRegion codingblocks
            if cbg.IS_IGNORED: continue
            if cbg.__class__.__name__ == 'LowSimilarityRegionCodingBlockGraph': continue
            if cbg_min_node_count and cbg.node_count() < cbg_min_node_count: continue
            # skip CBGs that are optimal already in terms of predicted splice sites
            if side == 'left' and cbg._CBGinterface5p and cbg._CBGinterface5p.is_optimal(): continue
            if side == 'rigth' and cbg._CBGinterface3p and cbg._CBGinterface3p.is_optimal(): continue
            ###if side == 'rigth' and self.cbginterface_is_optimal_donor(cbg): continue
            ###if side == 'left' and self.cbginterface_is_optimal_acceptor(cbg): continue

            # skip CBGs that have a to low identity when this is asked for
            if cbg_min_identity > 0.0 and cbg_min_identity > cbg.genetree().identity(): continue
            # skip first/final CBGs if requested for
            if ignore_first_cbg and cbg.IS_FIRST:       continue
            if ignore_first_cbg and pos == 0:           continue
            if ignore_final_cbg and cbg.IS_LAST:        continue
            if ignore_final_cbg and pos == len(self)-1: continue

            prev, next = None, None
            if side == 'left':
                # split on the LEFT side; next==currentCBG
                if pos >= 1: prev = self.codingblockgraphs[pos-1]
                else:        prev = None
                # set variables to confirm the splitted CBG agains
                compareCBG = prev
                pos_in_splits = 0
            else:
                # split on the rigth side; prev==currentCBG
                if pos < len(self)-1: next = self.codingblockgraphs[pos+1]
                else:                 next = None
                # set variables to confirm the splitted CBG agains
                compareCBG = next
                pos_in_splits = -1

            # do not check for sprdifference when the next/prev (compareCBG)
            # isa LowSimilarityRegionCodingBlockGraph
            if compareCBG and compareCBG.__class__.__name__ == 'LowSimilarityRegionCodingBlockGraph':
                continue

            # the nodes that are represented in the spanningrange difference
            sprdif = codingblock_splitting.spanningrange_difference(cbg,side,
                        sprdif_min_aa_length=sprdif_min_aa_length,
                        sprdif_min_node_count=sprdif_min_node_count)
            sprdifnodes = sprdif.keys()

            if not sprdifnodes:
                if verbose:
                    print "no sprdif nodes", pos, (sprdif_min_aa_length, sprdif_min_node_count)
                continue

            # only (try to) perform a split in there is a difference in nodes between
            # sprdifnodes and the nodes in compareCBG
            if compareCBG and not Set(sprdifnodes).difference(compareCBG.get_nodes()):
                continue

            # gather information of the distance between the CBGs
            if prev:
                # gather distance per organism between codingblocks
                aadistomsr  = prev.distance_between_codingblocks(cbg)
                aadistmaxsr = prev.maxsr_distance_between_codingblocks(cbg)
                # gather mutual nodes & organisms
                mutual_nodes = Set(prev.get_nodes()).intersection(cbg.get_nodes())
                mutual_orgs  = [ cbg.organism_by_node(node) for node in mutual_nodes ]
            elif next:
                # gather distance per organism between codingblocks
                aadistomsr  = cbg.distance_between_codingblocks(next)
                aadistmaxsr = cbg.maxsr_distance_between_codingblocks(next)
                # gather mutual nodes & organisms
                mutual_nodes = Set(cbg.get_nodes()).intersection(next.get_nodes())
                mutual_orgs  = [ cbg.organism_by_node(node) for node in mutual_nodes ]
            else:
                mutual_nodes = []
                mutual_orgs  = []

            # gather status_list for potential spanningrange difference
            # each organism has a single vote: True, None, False
            sprdif_status_list = []
            for org in mutual_orgs:
                if aadistomsr[org] < sprdif_min_aa_length: 
                    sprdif_status_list.append(False)
                else:
                    node = cbg.node_by_organism(org)
                    if node in sprdifnodes:
                        sprdif_status_list.append(True)
                    else:
                        sprdif_status_list.append(None)
            for org in cbg.organism_set().difference(mutual_orgs):
                node = cbg.node_by_organism(org)
                if node in sprdifnodes:
                    sprdif_status_list.append(True)
                else:
                    sprdif_status_list.append(None)

            ############################################################################
            if verbose:
                print stw.lap(), "STATUSsprdif:", side, "side", sprdif_status_list, " must have >1 Trues and 0 Falses!"
                print cbg
                print "sprdifnodes:", sprdifnodes
                if prev or next:
                    mutual_nodes = list(mutual_nodes)
                    mutual_nodes.sort()
                    for node in cbg.get_ordered_nodes():
                        if node not in mutual_nodes:
                            org = cbg.organism_by_node(node)
                            if not aadistomsr.has_key(org): continue
                            omsr = ( min(cbg.overall_minimal_spanning_range(node=node)),
                                     max(cbg.overall_minimal_spanning_range(node=node))
                                   )
                            sprdiftxt = ""
                            if node in sprdifnodes:
                                sprdiftxt = min(sprdif[node]), max(sprdif[node]), len(sprdif[node])
                            print " ", org, "\tomsr:", omsr, "\taadist: ", aadistomsr[org],
                            print "\tdistMAXSR:", aadistmaxsr[org], sprdiftxt
                        else:
                            org = cbg.organism_by_node(node)
                            if not aadistomsr.has_key(org): continue
                            omsr = ( min(cbg.overall_minimal_spanning_range(node=node)),
                                     max(cbg.overall_minimal_spanning_range(node=node))
                                   )
                            sprdiftxt = ""
                            if node in sprdifnodes:
                                sprdiftxt = min(sprdif[node]), max(sprdif[node]), len(sprdif[node])
                            print "M", org, "\tomsr:", omsr, "\taadist: ", aadistomsr[org],
                            print "\tdistMAXSR:", aadistmaxsr[org], sprdiftxt

                else:
                    for node in cbg.get_nodes():
                        org = cbg.organism_by_node(node)
                        omsr = ( min(cbg.overall_minimal_spanning_range(node=node)),
                                 max(cbg.overall_minimal_spanning_range(node=node))
                               )
                        sprdiftxt = ""
                        if node in sprdifnodes:
                            sprdiftxt = min(sprdif[node]), max(sprdif[node]), len(sprdif[node])
                        print " ", org, "\tomsr:", omsr, "\taadist: ",
                        print "n.a.", "\tdistMAXSR:", "n.a.", sprdiftxt
            ############################################################################

            # check the status list; only continue when NO False and at least a single True
            if not (True in sprdif_status_list and not False in sprdif_status_list): continue
 

            # reset sprdif_status_list to empty; scoring starts all over
            sprdif_status_list = []
            # Check the pacbporfs for identityscore in their sprdif range
            # When the ratio between sprdif / omsr identityscore ratio is to low
            # this is likely the event of a non-relevant run-through alignment,
            # either in less-similar intronic region (gaps!) or complete nonsense alignment
            # Although it seems a lot of work/calculation time, this step detects many
            # sprdifs that will never yield a succesfull new CBG. Detecting this at this
            # point in facts saves a lot of time (starting with the cbg deepcopy step!)
            for nodeA,nodeB in cbg.pairwisecrosscombinations_node():
                if not (nodeA in sprdifnodes and nodeB in sprdifnodes): continue
                omsr = ( min(cbg.overall_minimal_spanning_range(node=nodeA)), max(cbg.overall_minimal_spanning_range(node=nodeA)) )
                pacbporf = cbg.get_pacbps_by_nodes(nodeA,nodeB)[0]
                omsr_identityscore = pacbporf.identityscore_slice_by_abs_protein_query( omsr[0], omsr[1]+1 )
                if len(sprdif[nodeA]) <= len(sprdif[nodeB]):
                    bitscore = pacbporf.bitscore_slice_by_abs_protein_query( min(sprdif[nodeA]), max(sprdif[nodeA])+1 )
                    sprdif_identityscore = pacbporf.identityscore_slice_by_abs_protein_query( min(sprdif[nodeA]), max(sprdif[nodeA])+1 )
                    ###if verbose: print "takenQ", nodeA, min(sprdif[nodeA]), max(sprdif[nodeA])+1
                else:
                    bitscore = pacbporf.bitscore_slice_by_abs_protein_sbjct( min(sprdif[nodeB]), max(sprdif[nodeB])+1 )
                    sprdif_identityscore = pacbporf.identityscore_slice_by_abs_protein_sbjct( min(sprdif[nodeB]), max(sprdif[nodeB])+1 )
                    ###if verbose: print "takenS", nodeA, min(sprdif[nodeB]), max(sprdif[nodeB])+1

                # calculate ratio between the identityscore in the sprdif area
                # vc the identityscore in the omsr area. When this ratio is
                # to low, most likely a non-significant alignment!
                try:
                    ratio = sprdif_identityscore / omsr_identityscore
                    if ratio >= sprdif_vs_omsr_min_identityratio:
                        sprdif_status_list.append(True)
                    else:
                        sprdif_status_list.append(False)
                except ZeroDivisionError:
                    # hmm.. omsr_identityscore == 0.0 !?
                    # just assume okay, but htis is hard to believe ;-) 
                    ratio = 1.0
                    sprdif_status_list.append(True)

                ########################################################################
                if verbose:
                    print nodeA, nodeB, "length:", pacbporf.length, "bits (all, omsr, split):", pacbporf.bitscore,
                    print pacbporf.bitscore_slice_by_abs_protein_query( omsr[0], omsr[1]+1 ),
                    print bitscore,
                    print "identscore (all, omsr,split):",
                    print "%1.3f, %1.3f, %1.3f" % (pacbporf.identityscore, omsr_identityscore,sprdif_identityscore),
                    print "RATIO: %1.3f" % ratio
                ########################################################################

            # check the sprdif_status_list; now, at least a single True is a pass!
            if verbose: print stw.lap(), "sprdif_ratio_list:", sprdif_status_list
            if not True in sprdif_status_list: continue

            # here we start trying a split on sprdif -> increase counter
            sprdif_splits_tried_cnt+=1

            # split with novel ClustalW code.
            newsplittedcbg = sprdif2clustalw2cbg(cbg,sprdif,verbose=False)

            ####################################################################
            if verbose:
                print stw.lap(), pos, "splitted", "PREV .. split .. NEXT"
                print "sprdif:", sprdif
                print prev
                print newsplittedcbg
                print next
            ####################################################################

            # for backwards compatibility with previous code store -- single --
            # CBG in a list and create variable (name) bckp_cbg
            if newsplittedcbg:
                splits = [ newsplittedcbg ]
                bckp_cbg = cbg
                if side == 'left':  next = bckp_cbg
                else:               prev = bckp_cbg

                if newsplittedcbg.node_count() > 2:
                    # create as well the `minimal` sprdif by limiting
                    # to the 2 strongest nodes
                    minimalsplittedcbg = newsplittedcbg.deepcopy()
                    while minimalsplittedcbg.node_count() > 2:
                        weakestnode = minimalsplittedcbg.weakest_connected_node()                
                        minimalsplittedcbg.del_node(weakestnode)
                        minimalsplittedcbg.update_edge_weights_by_minimal_spanning_range()
                        minimalsplittedcbg.create_cache()
                    splits.append( minimalsplittedcbg )
                    ############################################################
                    if verbose:
                        print "WEAKEST NODE REMOVED:"
                        print minimalsplittedcbg
                    ############################################################

            else:
                continue


            prevCBG, nextCBG = None, None  # new variables for prev and next CBG
            prev_cbg_is_omitted = False    # Boolean if prev_cbg is `removed` from partGSG
            next_cbg_is_omitted = False    # Boolean if next_cbg is `removed` from partGSG
            partgsgCBGpos_start = None     # slice offset of the partGSG (based on its size!)
            partgsgCBGpos_end   = None     # slice offset of the partGSG (based on its size!)
            # create a partialGSG to insert the potential novel CBGs into
            if side == 'left' and prev:
                # default modus. Make a partialGSG of prev/bckp
                partGSG = mutualorgsofcbgs2partialGSG(input=self.input,cbgs=[prev,bckp_cbg])
                nextCBG = bckp_cbg 
                prevCBG = prev
                partgsgCBGpos_start = pos-1 # based on variable pos in most outern for-loop
                partgsgCBGpos_end   = pos+1 # based on variable pos in most outern for-loop

                # check if the sprdif splitted CBG overlaps with prev
                bestCBG = splits[0]
                if remove_lower_scoring_overlapping_cbgs and\
                bestCBG.total_weight() > prev.total_weight() * 2:
                    # best splitted CBG is higher scoring than currently prev CBG
                    if not partGSG.add_codingblock(bestCBG,only_try_adding=True,
                    omit_conditional_addition=True):
                        # although higher scoring, the splitted CBG can not be
                        # placed in the GSG! Remove `prev` from the partGSG
                        partGSG = mutualorgsofcbgs2partialGSG(input=self.input,cbgs=[bckp_cbg])
                        prevCBG = None
                        prev_cbg_is_omitted = True
                        partgsgCBGpos_start = pos   # based on variable pos in most outern for-loop
                        partgsgCBGpos_end   = pos+1 # based on variable pos in most outern for-loop
                    
            elif side == 'rigth' and next:
                # default modus. Make a partialGSG of bckp/next CBG
                partGSG = mutualorgsofcbgs2partialGSG(input=self.input,cbgs=[bckp_cbg,next])
                prevCBG = bckp_cbg
                nextCBG = next
                partgsgCBGpos_start = pos   # based on variable pos in most outern for-loop
                partgsgCBGpos_end   = pos+2 # based on variable pos in most outern for-loop

                # check if the sprdif splitted CBG overlaps with next 
                bestCBG = splits[0]
                if remove_lower_scoring_overlapping_cbgs and\
                bestCBG.total_weight() > next.total_weight() * 2:
                    # best splitted CBG is higher scoring than currently next CBG
                    if not partGSG.add_codingblock(bestCBG,only_try_adding=True,
                    omit_conditional_addition=True):
                        # although higher scoring, the splitted CBG can not be
                        # placed in the GSG! Remove `next` from the partGSG
                        partGSG = mutualorgsofcbgs2partialGSG(input=self.input,cbgs=[bckp_cbg])
                        nextCBG = None
                        next_cbg_is_omitted = True
                        partgsgCBGpos_start = pos   # based on variable pos in most outern for-loop
                        partgsgCBGpos_end   = pos+1 # based on variable pos in most outern for-loop

            else:
                # no prev or next CBG available in the GSG on the splitted side
                partGSG = mutualorgsofcbgs2partialGSG(input=self.input,cbgs=[bckp_cbg])
                partgsgCBGpos_start = pos   # based on variable pos in most outern for-loop
                partgsgCBGpos_end   = pos+1 # based on variable pos in most outern for-loop

            # copy genetree object from main GSG to partGSG
            partGSG._GENETREE = self._GENETREE
            curpartgsglen = len(partGSG)

            ####################################################################
            if verbose:
                print partGSG, side, "remove_lower_scoring_overlapping_cbgs:",
                print remove_lower_scoring_overlapping_cbgs
                print "nextCBG:", nextCBG
                print "prevCBG:", prevCBG
                print "OMITTED??, nextCBG", next_cbg_is_omitted, 
                print "prevCBG", prev_cbg_is_omitted
            ####################################################################

            # list with all accepted CBGs to be placed into the GSG
            all_accepted_cbgs = []

            # loop over the splits
            for splittedCBG in splits:
                if splittedCBG.__class__.__name__ == 'LowSimilarityRegionCodingBlockGraph':
                    continue
    
                # try to add into the genestructure; when it does not fit
                # into the partGSG -> no sense in continueing here
                if not partGSG.add_codingblock(splittedCBG,only_try_adding=True,
                    omit_conditional_addition=True):
                    ####################################################################
                    if verbose: print "not placeable!:", splittedCBG
                    ####################################################################
                    continue

                ########################################################################
                if verbose:
                    print "->split:", splittedCBG
                    print "->split:", splittedCBG.node_count(),
                    print splittedCBG.edge_count(),
                    print splittedCBG.genetree().identity()
                    splittedCBG.printmultiplealignment()
                    splittedCBG.cexpanderanalyses()
                    print splittedCBG._cexpander.binarystring,
                    print "projected_on:", splittedCBG._cexpander.projected_on,
                    print stw.lap()
                ########################################################################

                if sprdif_min_identity > 0.0 and sprdif_min_identity > splittedCBG.genetree().identity():
                    ########################################################################
                    if verbose:
                        print "sprdif_min_identity (%1.2f) >" % sprdif_min_identity,
                        print "splittedCBG.genetree().identity() -> ignore"
                    ########################################################################
                    continue


                if not perform_cbgif_optimization:
                    # Default modus; complete with cbghmmsearch2pacbpcollection
                    pacbpCollection = cbghmmsearch2pacbpcollection(splittedCBG,self.input,
                        next=nextCBG,prev=prevCBG,
                        pacbp_min_length=sprdif_min_aa_length,
                        hmmsearch_num_hits=3,verbose=verbose,
                        )
                else:
                    # special case; perform_cbgif_optimization make artificial input dict
                    dummyinput = {}
                    for org in self.input.keys():
                        partialorflist = OrfSet()
                        if nextCBG and org in nextCBG.organism_set():
                            orf = nextCBG.get_orfs_of_graph(organism=org)[0]
                            partialorflist.orfs.append(orf)
                        if prevCBG and org in prevCBG.organism_set():
                            orf = prevCBG.get_orfs_of_graph(organism=org)[0]
                            if orf.id not in [ o.id for o in partialorflist.orfs]:
                                partialorflist.orfs.append(orf)
                        # create key in dummyinput dict
                        dummyinput[org] = {
                                'orfs'      : partialorflist,
                                'genomeseq' : self.input[org]['genomeseq'],
                        }
                    ############################################################
                    if verbose:
                        print stw.lap(), "perform_cbgif_optimization input"
                    ############################################################
                    # Now do cbghmmsearch2pacbpcollection
                    pacbpCollection = cbghmmsearch2pacbpcollection(splittedCBG,
                        dummyinput,next=nextCBG,prev=prevCBG,
                        pacbp_min_length=sprdif_min_aa_length,
                        hmmsearch_num_hits=2,verbose=verbose
                        )


                ################################################################
                if verbose:
                    print "pacbpCollection(%s,%s) created, %s" % (
                        pacbpCollection.node_count(),   
                        pacbpCollection.edge_count(),
                        stw.lap()
                        )
                    for org in pacbpCollection.organism_set():
                        print org, [ node[1] for node in pacbpCollection.get_organism_nodes(org) ]
                ################################################################

                # get list of accepted CBGs
                accepted =  conversion.pacbpCollection2AcceptedCodingBlockGraphs(
                        pacbpCollection,prev=prevCBG,next=nextCBG)

                ################################################################
                if verbose:
                    print len(accepted), "CBGs converted from pacbpCollection"
                ################################################################

                for acc_cbg in accepted:
                    # update the acc_cbg on edge_weight and create cache
                    acc_cbg.update_edge_weights_by_minimal_spanning_range()
                    acc_cbg.create_cache()

                    # check if all organisms/genes are covered;
                    # a missing node/org/gene in pacbpCollection
                    # results automatically in a - non detected - missing piece
                    if acc_cbg.node_count() != self.EXACT_SG_NODE_COUNT:
                        ########################################################
                        if verbose:
                           print "DISCARDED != self.EXACT_SG_NODE_COUNT",
                           print acc_cbg.get_ordered_nodes()
                        ########################################################
                        continue

                    # check if this CBG does not share all nodes (orfs)
                    # with its parent
                    if cbg.get_ordered_nodes() == acc_cbg.get_ordered_nodes():
                        ########################################################
                        if verbose:
                           print "DISCARDED identical nodes:",
                           print acc_cbg.get_ordered_nodes()
                        ########################################################
                        continue

                    # check cexpander binarystring for all zeros.
                    if do_cexander_allzeros_check and\
                    acc_cbg._cexpander.binarystring.find("1") == -1:
                        ########################################################
                        if verbose:
                           print "DISCARDED cexpander 0000",
                           print acc_cbg.get_ordered_nodes()
                        ########################################################
                        continue


                    # remove non-uniformly aligned fraction in the cbgs
                    try:
                        status = lib_cexpander.cexpander_checkCBG4omsrbordergaps(acc_cbg)
                    except lib_cexpander.ZeroUniformlyAlignedPositions:
                        ########################################################
                        if verbose:
                           print "DISCARDED cexpander.checkCBG4omsrbordergaps 0000",
                           print acc_cbg.get_ordered_nodes()
                        ########################################################
                        continue
                    except NoOverallMinimalSpanningRange:
                        ########################################################
                        if verbose:
                           print "DISCARDED cexpander.checkCBG4omsrbordergaps ",
                           print "NO OMSR", acc_cbg.get_ordered_nodes()
                        ########################################################
                        continue

                    # remove CBGs with spurious node sets compared to its neighbours
                    nodestat = intermediateCBG_node_comparison(prevCBG,acc_cbg,nextCBG)
                    if nodestat==False:
                        ########################################################
                        if verbose:
                           print "DISCARDED intermediateCBG_node_comparison ",
                           print acc_cbg
                        ########################################################
                        continue
 
                    # if here, this intermediate CBG is accepted as a possible insert
                    # Replace proper pacbporfs from the parental CBG
                    replacements = acc_cbg._recrute_pacbporfs_from_parental_cbg(bckp_cbg,verbose=verbose)
                    if not acc_cbg.has_overall_minimal_spanning_range():
                        ########################################################
                        if verbose:
                           print "DISCARDED recrute_pacbporfs_from_parent CBGs ",
                           print "NO OMSR", acc_cbg
                        ########################################################
                        continue

                    # place in list of suitable CBGs to check for acceptance
                    all_accepted_cbgs.append( acc_cbg)

                    ############################################################
                    if verbose:
                        print acc_cbg
                        acc_cbg.printmultiplealignment()
                        for trf in acc_cbg._cexpander._transferblocks:
                            print trf.binarystring, trf.projected_on
                    ############################################################


            # all_accepted_cbgs isa list with K(s) CBGs that are potential
            # enrichments for the GSG. Try adding them to the GSG by
            # finding most likely CBG2GSGinsert
            partGSG = findmostlikelyCBG2GSGinsert(partGSG,all_accepted_cbgs,verbose=verbose)
            ####################################################################
            if verbose:
                print partGSG, "graphs remaining:",
                print len(all_accepted_cbgs), stw.lap()
            ####################################################################
            if len(partGSG) > curpartgsglen:
                # update the slice into the main GSG
                self.codingblockgraphs.__setslice__(
                        partgsgCBGpos_start,
                        partgsgCBGpos_end,
                        partGSG.codingblockgraphs )
                # update counter that novel CBGs are created
                sprdif_splits_done_cnt += len(partGSG) - curpartgsglen

                if next_cbg_is_omitted or prev_cbg_is_omitted:
                    # There was a total_weight discrepancy observed.
                    # Next/prev CBG was omitted from the partGSG.
                    # Chance is very high now that a overlap in the GSG is created.
                    # Resolve it right here & now
                    print "WARNING!!!:, next/prev cbg_is_omitted", next_cbg_is_omitted, prev_cbg_is_omitted
                    print self
                    for thecbg in self: print thecbg
                    removed_cnt = self.remove_overlapping_cbgs(verbose=True)
                    print self, removed_cnt
                    for thecbg in self: print thecbg

                ################################################################
                if verbose:
                    print "CBG(s) added in GSG after sprdif", side
                    for cbg in partGSG: print cbg 
                ################################################################
                continue
            else:
                # no succesfull split performed.
                # DO NOT FORGET TO SET BACK THE ORIGINAL INTERFACE!
                # python is call-by-reference, so this updates the
                # cbgIFs in the main GSG
                if len(partGSG) >= 2:
                    partGSG.clear_central_cbginterfaces()
                    partGSG.create_cbginterfaces()
                else:
                    partGSG._CBGinterface5p = None
                    partGSG._CBGinterface3p = None
                # continue with the next try-to-split
                continue

        # correct possible erroneous setted FIRST/LAST CBGs
        if sprdif_splits_done_cnt: self.finalize_genestructure()

        # return how much new cbgs are added
        return sprdif_splits_done_cnt

    # end of function split_cbgs_on_spanningrange_difference


    def WORKING_split_cbgs_on_spanningrange_difference(self,side=None,
        sprdif_min_aa_length=CBG_MIN_AA_LENGTH,
        sprdif_min_node_count=CBG_SPRDIF_MIN_NODE_COUNT,
        sprdif_min_identity=0.0,     # 0.0 means by default no check on this property!
        cbg_min_identity=0.0,        # 0.0 means by default no check on this property!
        cbg_min_node_count=None,verbose=False,
        do_cexander_allzeros_check = True,
        ignore_first_cbg=False,
        ignore_final_cbg=False,
        perform_cbgif_optimization=False, # new 01/12/2009 !!
        remove_lower_scoring_overlapping_cbgs=False,
        sprdif_vs_omsr_min_identityratio=CBG_SPRDIF_VS_OMSR_MIN_IDENTITYSCORE_RATIO):
        """             
        Split CBGs on spanningrangedifference and try to assign a new CBG in the splitted area
            
        @type  side: string 
        @param side: 'left' or 'rigth', meaning a split on left/5p or rigth/3p side
                
        @type  sprdif_min_aa_length: integer
        @param sprdif_min_aa_length: minimal length of the sprdif in aa's
            
        @type  sprdif_min_node_count: integer
        @param sprdif_min_node_count: minimal number of nodes that support the sprdif

        @type  sprdif_min_identity: float
        @param sprdif_min_identity: minimal identity (0.0-1.0) of the splitted CBG to be accepted/tried

        @type  cbg_min_identity: float 
        @param cbg_min_identity: minimal identity (0.0-1.0) of the CBG to allow a split on sprdif

        @type  cbg_min_node_count: integer
        @param cbg_min_node_count: minimal number of nodes in a CBG to be elegiable for trying a split

        @type  do_cexander_allzeros_check: Boolean
        @param do_cexander_allzeros_check: omit novel CBGs that have not a single uniformly alignable AA position

        @type  ignore_first_cbg: Boolean
        @param ignore_first cbg: do or ignore the first CBG in the list / the CBG labelled as IS_FIRST

        @type  ignore_final_cbg: Boolean
        @param ignore_final cbg: do or ignore the final CBG in the list / the CBG labelled as IS_LAST
        """ 
        # input integrity check
        if side not in ('left','rigth'):
            message = "`side` must be 'left' or 'rigth', not '%s'" % side
            raise InproperlyAppliedArgument, message
        if sprdif_min_node_count < 2:
            message = "`sprdif_min_node_count` must be >=2, not '%s'" % sprdif_min_node_count
            raise InproperlyAppliedArgument, message

        stw = StopWatch(name='StwTestSprdif(MIN_AA_LENGT:%s)' % sprdif_min_aa_length)
        stw.start()

        # counters for how much splits have been tried/done
        sprdif_splits_done_cnt  = 0
        sprdif_splits_tried_cnt = 0

        # loop backwards in case a splitted CBG is inserted!
        for pos in self.cbgpositerator(reversed=True):
            cbg = self.codingblockgraphs[pos]
            # skip IGNORED and LowSimilarityRegion codingblocks
            if cbg.IS_IGNORED: continue
            if cbg.__class__.__name__ == 'LowSimilarityRegionCodingBlockGraph': continue
            if cbg_min_node_count and cbg.node_count() < cbg_min_node_count: continue
            # skip CBGs that are optimal already in terms of predicted splice sites
            if side == 'left' and cbg._CBGinterface5p and cbg._CBGinterface5p.is_optimal(): continue
            if side == 'rigth' and cbg._CBGinterface3p and cbg._CBGinterface3p.is_optimal(): continue
            ###if side == 'rigth' and self.cbginterface_is_optimal_donor(cbg): continue
            ###if side == 'left' and self.cbginterface_is_optimal_acceptor(cbg): continue

            # skip CBGs that have a to low identity when this is asked for
            if cbg_min_identity > 0.0 and cbg_min_identity > cbg.genetree().identity(): continue
            # skip first/final CBGs if requested for
            if ignore_first_cbg and cbg.IS_FIRST:       continue
            if ignore_first_cbg and pos == 0:           continue
            if ignore_final_cbg and cbg.IS_LAST:        continue
            if ignore_final_cbg and pos == len(self)-1: continue

            prev, next = None, None
            if side == 'left':
                # split on the LEFT side; next==currentCBG
                if pos >= 1: prev = self.codingblockgraphs[pos-1]
                else:        prev = None
                # set variables to confirm the splitted CBG agains
                compareCBG = prev
                pos_in_splits = 0
            else:
                # split on the rigth side; prev==currentCBG
                if pos < len(self)-1: next = self.codingblockgraphs[pos+1]
                else:                 next = None
                # set variables to confirm the splitted CBG agains
                compareCBG = next
                pos_in_splits = -1

            # do not check for sprdifference when the next/prev (compareCBG)
            # isa LowSimilarityRegionCodingBlockGraph
            if compareCBG and compareCBG.__class__.__name__ == 'LowSimilarityRegionCodingBlockGraph':
                continue

            # the nodes that are represented in the spanningrange difference
            sprdif = codingblock_splitting.spanningrange_difference(cbg,side,
                        sprdif_min_aa_length=sprdif_min_aa_length,
                        sprdif_min_node_count=sprdif_min_node_count)
            sprdifnodes = sprdif.keys()

            if not sprdifnodes:
                if verbose:
                    print "no sprdif nodes", pos, (sprdif_min_aa_length, sprdif_min_node_count)
                continue

            # only (try to) perform a split in there is a difference in nodes between
            # sprdifnodes and the nodes in compareCBG
            if compareCBG and not Set(sprdifnodes).difference(compareCBG.get_nodes()):
                continue

            # gather information of the distance between the CBGs
            if prev:
                # gather distance per organism between codingblocks
                aadistomsr  = prev.distance_between_codingblocks(cbg)
                aadistmaxsr = prev.maxsr_distance_between_codingblocks(cbg)
                # gather mutual nodes & organisms
                mutual_nodes = Set(prev.get_nodes()).intersection(cbg.get_nodes())
                mutual_orgs  = [ cbg.organism_by_node(node) for node in mutual_nodes ]
            elif next:
                # gather distance per organism between codingblocks
                aadistomsr  = cbg.distance_between_codingblocks(next)
                aadistmaxsr = cbg.maxsr_distance_between_codingblocks(next)
                # gather mutual nodes & organisms
                mutual_nodes = Set(cbg.get_nodes()).intersection(next.get_nodes())
                mutual_orgs  = [ cbg.organism_by_node(node) for node in mutual_nodes ]
            else:
                mutual_nodes = []
                mutual_orgs  = []

            # gather status_list for potential spanningrange difference
            # each organism has a single vote: True, None, False
            sprdif_status_list = []
            for org in mutual_orgs:
                if aadistomsr[org] < sprdif_min_aa_length: 
                    sprdif_status_list.append(False)
                else:
                    node = cbg.node_by_organism(org)
                    if node in sprdifnodes:
                        sprdif_status_list.append(True)
                    else:
                        sprdif_status_list.append(None)
            for org in cbg.organism_set().difference(mutual_orgs):
                node = cbg.node_by_organism(org)
                if node in sprdifnodes:
                    sprdif_status_list.append(True)
                else:
                    sprdif_status_list.append(None)

            ############################################################################
            if verbose:
                print stw.lap(), "STATUSsprdif:", side, "side", sprdif_status_list, " must have >1 Trues and 0 Falses!"
                print cbg
                print "sprdifnodes:", sprdifnodes
                if prev or next:
                    mutual_nodes = list(mutual_nodes)
                    mutual_nodes.sort()
                    for node in cbg.get_ordered_nodes():
                        if node not in mutual_nodes:
                            org = cbg.organism_by_node(node)
                            if not aadistomsr.has_key(org): continue
                            omsr = ( min(cbg.overall_minimal_spanning_range(node=node)),
                                     max(cbg.overall_minimal_spanning_range(node=node))
                                   )
                            sprdiftxt = ""
                            if node in sprdifnodes:
                                sprdiftxt = min(sprdif[node]), max(sprdif[node]), len(sprdif[node])
                            print " ", org, "\tomsr:", omsr, "\taadist: ", aadistomsr[org],
                            print "\tdistMAXSR:", aadistmaxsr[org], sprdiftxt
                        else:
                            org = cbg.organism_by_node(node)
                            if not aadistomsr.has_key(org): continue
                            omsr = ( min(cbg.overall_minimal_spanning_range(node=node)),
                                     max(cbg.overall_minimal_spanning_range(node=node))
                                   )
                            sprdiftxt = ""
                            if node in sprdifnodes:
                                sprdiftxt = min(sprdif[node]), max(sprdif[node]), len(sprdif[node])
                            print "M", org, "\tomsr:", omsr, "\taadist: ", aadistomsr[org],
                            print "\tdistMAXSR:", aadistmaxsr[org], sprdiftxt

                else:
                    for node in cbg.get_nodes():
                        org = cbg.organism_by_node(node)
                        omsr = ( min(cbg.overall_minimal_spanning_range(node=node)),
                                 max(cbg.overall_minimal_spanning_range(node=node))
                               )
                        sprdiftxt = ""
                        if node in sprdifnodes:
                            sprdiftxt = min(sprdif[node]), max(sprdif[node]), len(sprdif[node])
                        print " ", org, "\tomsr:", omsr, "\taadist: ",
                        print "n.a.", "\tdistMAXSR:", "n.a.", sprdiftxt
            ############################################################################

            # check the status list; only continue when NO False and at least a single True
            if not (True in sprdif_status_list and not False in sprdif_status_list): continue
 

            # reset sprdif_status_list to empty; scoring starts all over
            sprdif_status_list = []
            # Check the pacbporfs for identityscore in their sprdif range
            # When the ratio between sprdif / omsr identityscore ratio is to low
            # this is likely the event of a non-relevant run-through alignment,
            # either in less-similar intronic region (gaps!) or complete nonsense alignment
            # Although it seems a lot of work/calculation time, this step detects many
            # sprdifs that will never yield a succesfull new CBG. Detecting this at this
            # point in facts saves a lot of time (starting with the cbg deepcopy step!)
            for nodeA,nodeB in cbg.pairwisecrosscombinations_node():
                if not (nodeA in sprdifnodes and nodeB in sprdifnodes): continue
                omsr = ( min(cbg.overall_minimal_spanning_range(node=nodeA)), max(cbg.overall_minimal_spanning_range(node=nodeA)) )
                pacbporf = cbg.get_pacbps_by_nodes(nodeA,nodeB)[0]
                omsr_identityscore = pacbporf.identityscore_slice_by_abs_protein_query( omsr[0], omsr[1]+1 )
                if len(sprdif[nodeA]) <= len(sprdif[nodeB]):
                    bitscore = pacbporf.bitscore_slice_by_abs_protein_query( min(sprdif[nodeA]), max(sprdif[nodeA])+1 )
                    sprdif_identityscore = pacbporf.identityscore_slice_by_abs_protein_query( min(sprdif[nodeA]), max(sprdif[nodeA])+1 )
                    ###if verbose: print "takenQ", nodeA, min(sprdif[nodeA]), max(sprdif[nodeA])+1
                else:
                    bitscore = pacbporf.bitscore_slice_by_abs_protein_sbjct( min(sprdif[nodeB]), max(sprdif[nodeB])+1 )
                    sprdif_identityscore = pacbporf.identityscore_slice_by_abs_protein_sbjct( min(sprdif[nodeB]), max(sprdif[nodeB])+1 )
                    ###if verbose: print "takenS", nodeA, min(sprdif[nodeB]), max(sprdif[nodeB])+1

                # calculate ratio between the identityscore in the sprdif area
                # vc the identityscore in the omsr area. When this ratio is
                # to low, most likely a non-significant alignment!
                ratio = sprdif_identityscore / omsr_identityscore
                if ratio >= sprdif_vs_omsr_min_identityratio:
                    sprdif_status_list.append(True)
                else:
                    sprdif_status_list.append(False)

                ########################################################################
                if verbose:
                    print nodeA, nodeB, "length:", pacbporf.length, "bits (all, omsr, split):", pacbporf.bitscore,
                    print pacbporf.bitscore_slice_by_abs_protein_query( omsr[0], omsr[1]+1 ),
                    print bitscore,
                    print "identscore (all, omsr,split):",
                    print "%1.3f, %1.3f, %1.3f" % (pacbporf.identityscore, omsr_identityscore,sprdif_identityscore),
                    print "RATIO: %1.3f" % ratio
                ########################################################################

            # check the sprdif_status_list; now, at least a single True is a pass!
            if verbose: print stw.lap(), "sprdif_ratio_list:", sprdif_status_list
            if not True in sprdif_status_list: continue

            # here we start trying a split on sprdif -> increase counter
            sprdif_splits_tried_cnt+=1

            ####################################################################
            ### NEW 30/11/2009
            ### replace - slow- code that cat yield >1 splitted CBG with code
            ### that quickly produces a single split
            ### this saves a very time-consuming deepcopy....
            ###
            #### make a deepcopy before splitting
            ###bckp_cbg = deepcopy(cbg)
            ###cbg.clear_cache()
            ###
            ###if side == 'left':
            ###    # split on the LEFT side; next==currentCBG
            ###    next = bckp_cbg
            ###else:
            ###    # split on the rigth side; prev==currentCBG
            ###    prev = bckp_cbg
            ###
            #### split on spanningsrange difference
            ###splits = cbg.iteratively_split_codingblock_on_spanningrange_difference(
            ###            side=side,
            ###            sprdif_min_aa_length=sprdif_min_aa_length,
            ###            sprdif_min_node_count=sprdif_min_node_count,
            ###            )
            #### restore the original CBG (splitting might have caused changes)
            ###self.codingblockgraphs[pos] = bckp_cbg
            ####################################################################
            ###if verbose:
            ###    print stw.lap(), pos, "splitted", "PREV .. splits %s.. NEXT" % (
            ###        len(splits) )
            ###    print prev
            ###    for spl in splits: print spl
            ###    print next
            ####################################################################
            ###
            #### no splits -> no action here!
            ###if len(splits)==1:
            ###    if verbose: print pos, "IGNORED, no splits"
            ###    continue
            ###
            #### now do several checks upon the most outward new splitted CBG (pos_in_splits)
            ###if len(splits) >= 1 and compareCBG and not splits[pos_in_splits].node_set().difference(compareCBG.get_nodes()):
            ###    # new splitted CBG has exactly the same nodes as the next/prev (variable compareCBG)
            ###    # in the genestructure. This is NOT a missed extention -> no action here
            ###    if verbose: print pos, "IGNORED, no node difference"
            ###    continue
            ###
            #### now deal with the splitted; first get rif of the 'original' cbg in splits
            ###if side == 'left': _oricbg = splits.pop()
            ###else:              _oricbg = splits.pop(0) 
            ###
            #### order by total graph weight
            ###splits = ordering.order_graphlist_by_total_weight(splits)
            ###
            ####################################################################

            # split with novel ClustalW code.
            newsplittedcbg = sprdif2clustalw2cbg(cbg,sprdif,verbose=False)

            ####################################################################
            if verbose:
                print stw.lap(), pos, "splitted", "PREV .. split .. NEXT"
                print prev
                print newsplittedcbg
                print next
            ####################################################################

            # for backwards compatibility with previous code store -- single --
            # CBG in a list and create variable (name) bckp_cbg
            if newsplittedcbg:
                splits = [ newsplittedcbg ]
                bckp_cbg = cbg
                if side == 'left':  next = bckp_cbg
                else:               prev = bckp_cbg
            else:
                continue


            prevCBG, nextCBG = None, None  # new variables for prev and next CBG
            prev_cbg_is_omitted = False    # Boolean if prev_cbg is `removed` from partGSG
            next_cbg_is_omitted = False    # Boolean if next_cbg is `removed` from partGSG
            partgsgCBGpos_start = None     # slice offset of the partGSG (based on its size!)
            partgsgCBGpos_end   = None     # slice offset of the partGSG (based on its size!)
            # create a partialGSG to insert the potential novel CBGs into
            if side == 'left' and prev:
                # default modus. Make a partialGSG of prev/bckp
                partGSG = mutualorgsofcbgs2partialGSG(input=self.input,cbgs=[prev,bckp_cbg])
                nextCBG = bckp_cbg 
                prevCBG = prev
                partgsgCBGpos_start = pos-1 # based on variable pos in most outern for-loop
                partgsgCBGpos_end   = pos+1 # based on variable pos in most outern for-loop

                # check if the sprdif splitted CBG overlaps with prev
                bestCBG = splits[0]
                if remove_lower_scoring_overlapping_cbgs and\
                bestCBG.total_weight() > prev.total_weight() * 2:
                    # best splitted CBG is higher scoring than currently prev CBG
                    if not partGSG.add_codingblock(bestCBG,only_try_adding=True,
                    omit_conditional_addition=True):
                        # although higher scoring, the splitted CBG can not be
                        # placed in the GSG! Remove `prev` from the partGSG
                        partGSG = mutualorgsofcbgs2partialGSG(input=self.input,cbgs=[bckp_cbg])
                        prevCBG = None
                        prev_cbg_is_omitted = True
                        partgsgCBGpos_start = pos   # based on variable pos in most outern for-loop
                        partgsgCBGpos_end   = pos+1 # based on variable pos in most outern for-loop
                    
            elif side == 'rigth' and next:
                # default modus. Make a partialGSG of bckp/next CBG
                partGSG = mutualorgsofcbgs2partialGSG(input=self.input,cbgs=[bckp_cbg,next])
                prevCBG = bckp_cbg
                nextCBG = next
                partgsgCBGpos_start = pos   # based on variable pos in most outern for-loop
                partgsgCBGpos_end   = pos+2 # based on variable pos in most outern for-loop

                # check if the sprdif splitted CBG overlaps with next 
                bestCBG = splits[0]
                if remove_lower_scoring_overlapping_cbgs and\
                bestCBG.total_weight() > next.total_weight() * 2:
                    # best splitted CBG is higher scoring than currently next CBG
                    if not partGSG.add_codingblock(bestCBG,only_try_adding=True,
                    omit_conditional_addition=True):
                        # although higher scoring, the splitted CBG can not be
                        # placed in the GSG! Remove `next` from the partGSG
                        partGSG = mutualorgsofcbgs2partialGSG(input=self.input,cbgs=[bckp_cbg])
                        nextCBG = None
                        next_cbg_is_omitted = True
                        partgsgCBGpos_start = pos   # based on variable pos in most outern for-loop
                        partgsgCBGpos_end   = pos+1 # based on variable pos in most outern for-loop

            else:
                # no prev or next CBG available in the GSG on the splitted side
                partGSG = mutualorgsofcbgs2partialGSG(input=self.input,cbgs=[bckp_cbg])
                partgsgCBGpos_start = pos   # based on variable pos in most outern for-loop
                partgsgCBGpos_end   = pos+1 # based on variable pos in most outern for-loop

            # copy genetree object from main GSG to partGSG
            partGSG._GENETREE = self._GENETREE
            curpartgsglen = len(partGSG)

            ####################################################################
            if verbose:
                print partGSG, side, "remove_lower_scoring_overlapping_cbgs:",
                print remove_lower_scoring_overlapping_cbgs
                print "nextCBG:", nextCBG
                print "prevCBG:", prevCBG
                print "OMITTED??, nextCBG", next_cbg_is_omitted, 
                print "prevCBG", prev_cbg_is_omitted
            ####################################################################

            # list with all accepted CBGs to be placed into the GSG
            all_accepted_cbgs = []

            # loop over the splits
            for splittedCBG in splits:
                if splittedCBG.__class__.__name__ == 'LowSimilarityRegionCodingBlockGraph':
                    continue
    
                # try to add into the genestructure; when it does not fit
                # into the partGSG -> no sense in continueing here
                if not partGSG.add_codingblock(splittedCBG,only_try_adding=True,
                    omit_conditional_addition=True):
                    ####################################################################
                    if verbose: print "not placeable!:", splittedCBG
                    ####################################################################
                    continue

                ########################################################################
                if verbose:
                    print "->split:", splittedCBG
                    print "->split:", splittedCBG.node_count(),
                    print splittedCBG.edge_count(),
                    print splittedCBG.genetree().identity()
                    splittedCBG.printmultiplealignment()
                    splittedCBG.cexpanderanalyses()
                    print splittedCBG._cexpander.binarystring,
                    print "projected_on:", splittedCBG._cexpander.projected_on,
                    print stw.lap()
                ########################################################################

                if sprdif_min_identity > 0.0 and sprdif_min_identity > splittedCBG.genetree().identity():
                    ########################################################################
                    if verbose:
                        print "sprdif_min_identity (%1.2f) >" % sprdif_min_identity,
                        print "splittedCBG.genetree().identity() -> ignore"
                    ########################################################################
                    continue


                if not perform_cbgif_optimization:
                    # Default modus; complete with cbghmmsearch2pacbpcollection
                    pacbpCollection = cbghmmsearch2pacbpcollection(splittedCBG,self.input,
                        next=nextCBG,prev=prevCBG,
                        pacbp_min_length=sprdif_min_aa_length,
                        hmmsearch_num_hits=3,verbose=verbose,
                        )
                else:
                    # special case; perform_cbgif_optimization make artificial input dict
                    dummyinput = {}
                    for org in self.input.keys():
                        partialorflist = OrfSet()
                        if nextCBG and org in nextCBG.organism_set():
                            orf = nextCBG.get_orfs_of_graph(organism=org)[0]
                            partialorflist.orfs.append(orf)
                        if prevCBG and org in prevCBG.organism_set():
                            orf = prevCBG.get_orfs_of_graph(organism=org)[0]
                            if orf.id not in [ o.id for o in partialorflist.orfs]:
                                partialorflist.orfs.append(orf)
                        # create key in dummyinput dict
                        dummyinput[org] = {
                                'orfs'      : partialorflist,
                                'genomeseq' : self.input[org]['genomeseq'],
                        }
                    ############################################################
                    if verbose:
                        print stw.lap(), "perform_cbgif_optimization input"
                    ############################################################
                    # Now do cbghmmsearch2pacbpcollection
                    pacbpCollection = cbghmmsearch2pacbpcollection(splittedCBG,
                        dummyinput,next=nextCBG,prev=prevCBG,
                        pacbp_min_length=sprdif_min_aa_length,
                        hmmsearch_num_hits=2,verbose=verbose
                        )


                ################################################################
                if verbose:
                    print "pacbpCollection(%s,%s) created, %s" % (
                        pacbpCollection.node_count(),   
                        pacbpCollection.edge_count(),
                        stw.lap()
                        )
                    for org in pacbpCollection.organism_set():
                        print org, [ node[1] for node in pacbpCollection.get_organism_nodes(org) ]
                ################################################################

                # get list of accepted CBGs
                accepted =  conversion.pacbpCollection2AcceptedCodingBlockGraphs(
                        pacbpCollection,prev=prevCBG,next=nextCBG)

                ################################################################
                if verbose:
                    print len(accepted), "CBGs converted from pacbpCollection"
                    for accCBG in accepted: print accCBG
                ################################################################


                for acc_cbg in accepted:
                    # update the acc_cbg on edge_weight and create cache
                    acc_cbg.update_edge_weights_by_minimal_spanning_range()
                    acc_cbg.create_cache()

                    # check if all organisms/genes are covered;
                    # a missing node/org/gene in pacbpCollection
                    # results automatically in a - non detected - missing piece
                    if acc_cbg.node_count() != self.EXACT_SG_NODE_COUNT: continue

                    # check cexpander binarystring for all zeros.
                    if do_cexander_allzeros_check and\
                    acc_cbg._cexpander.binarystring.find("1") == -1:
                        # ignore all zeros here...
                        continue

                    # replace proper pacbporfs from the parental CBG
                    replacements = acc_cbg._recrute_pacbporfs_from_parental_cbg(bckp_cbg)
                    # place in list of suitable CBGs to check for acceptance
                    all_accepted_cbgs.append( acc_cbg)

            # all_accepted_cbgs isa list with K(s) CBGs that are potential
            # enrichments for the GSG. Try adding them to the GSG by
            # finding most likely CBG2GSGinsert
            partGSG = findmostlikelyCBG2GSGinsert(partGSG,all_accepted_cbgs,verbose=verbose)
            if verbose:
                print partGSG, "graphs remaining:", len(all_accepted_cbgs), stw.lap()

            if len(partGSG) > curpartgsglen:
                # update the slice into the main GSG
                self.codingblockgraphs.__setslice__(
                        partgsgCBGpos_start,
                        partgsgCBGpos_end,
                        partGSG.codingblockgraphs )
                # update counter that novel CBGs are created
                sprdif_splits_done_cnt += len(partGSG) - curpartgsglen

                if next_cbg_is_omitted or prev_cbg_is_omitted:
                    # There was a total_weight discrepancy observed.
                    # Next/prev CBG was omitted from the partGSG.
                    # Chance is very high now that a overlap in the GSG is created.
                    # Resolve it right here & now
                    print "WARNING!!!:, next/prev cbg_is_omitted", next_cbg_is_omitted, prev_cbg_is_omitted
                    print self
                    for thecbg in self: print thecbg
                    removed_cnt = self.remove_overlapping_cbgs(verbose=True)
                    print self, removed_cnt
                    for thecbg in self: print thecbg

                ################################################################
                if verbose:
                    print "CBG(s) added in GSG after sprdif", side
                    for cbg in partGSG: print cbg 
                ################################################################
                continue
            else:
                # no succesfull split performed.
                # DO NOT FORGET TO SET BACK THE ORIGINAL INTERFACE!
                # python is call-by-reference, so this updates the
                # cbgIFs in the main GSG
                if len(partGSG) >= 2:
                    partGSG.clear_central_cbginterfaces()
                    partGSG.create_cbginterfaces()
                else:
                    partGSG._CBGinterface5p = None
                    partGSG._CBGinterface3p = None
                # continue with the next try-to-split
                continue

        # return how much new cbgs are added
        return sprdif_splits_done_cnt

    # end of function split_cbgs_on_spanningrange_difference


    def split_cbgs_on_left_spanningrange_difference(self,**kwargs):
        """
        Split CBGs in genestructure on left spanningrange difference

        @attention: see split_cbgs_on_spanningrange_difference for documentation
        """
        return self.split_cbgs_on_spanningrange_difference(side='left',**kwargs)

    # end of function split_cbgs_on_left_spanningrange_difference


    def split_cbgs_on_rigth_spanningrange_difference(self,**kwargs):
        """
        Split CBGs in genestructure on rigth spanningrange difference

        @attention: see split_cbgs_on_spanningrange_difference for documentation
        """
        return self.split_cbgs_on_spanningrange_difference(side='rigth',**kwargs)

    # end of function split_cbgs_on_rigth_spanningrange_difference



    def create_large_intermediate_cbg_for_scaffold_gap(self,side=None,
        sprdif_min_aa_length=CBG_MIN_AA_LENGTH,
        sprdif_min_node_count=CBG_SPRDIF_MIN_NODE_COUNT,
        sprdif_min_identity=0.0,     # 0.0 means by default no check on this property!
        cbg_min_identity=0.0,        # 0.0 means by default no check on this property!
        cbg_min_node_count=None,
        do_cexander_allzeros_check = True,
        verbose=False,
        ):
        """

        @type  sprdif_min_aa_length: integer
        @param sprdif_min_aa_length: minimal length of the sprdif in aa's

        @type  sprdif_min_node_count: integer
        @param sprdif_min_node_count: minimal number of nodes that support the sprdif

        @type  sprdif_min_identity: float
        @param sprdif_min_identity: minimal identity (0.0-1.0) of the splitted CBG to be accepted/tried

        @type  cbg_min_identity: float
        @param cbg_min_identity: minimal identity (0.0-1.0) of the CBG to allow a split on sprdif

        @type  cbg_min_node_count: integer
        @param cbg_min_node_count: minimal number of nodes in a CBG to be elegiable for trying a split

        @type  do_cexander_allzeros_check: Boolean
        @param do_cexander_allzeros_check: omit novel CBGs that have not a single uniformly alignable AA position
        """
        # input integrity check
        if sprdif_min_node_count < 2:
            message = "`sprdif_min_node_count` must be >=2, not '%s'" % sprdif_min_node_count
            raise InproperlyAppliedArgument, message

        NEW_INTERMEDIATE_CBG_IS_CREATED = False

        # loop backwards in case a splitted CBG is inserted!
        for pos in range(len(self)-1,0,-1):
            next = self.codingblockgraphs[pos]
            prev = self.codingblockgraphs[pos-1]
            # skip IGNORED and LowSimilarityRegion codingblocks
            if next.IS_IGNORED: continue
            if prev.IS_IGNORED: continue
            if next.__class__.__name__ == 'LowSimilarityRegionCodingBlockGraph': continue
            if prev.__class__.__name__ == 'LowSimilarityRegionCodingBlockGraph': continue
            # only continue in case CBGs have a certain number of nodes
            # recommended is to use here Ks graph node count!
            if cbg_min_node_count and next.node_count() < cbg_min_node_count: continue
            if cbg_min_node_count and prev.node_count() < cbg_min_node_count: continue
            # skip CBGs that have a to low identity when this is asked for
            if cbg_min_identity > 0.0 and cbg_min_identity > prev.genetree().identity(): continue
            if cbg_min_identity > 0.0 and cbg_min_identity > next.genetree().identity(): continue

            # the nodes that are represented in the spanningrange difference
            sprdifprev = codingblock_splitting.spanningrange_difference(
                    prev, 'rigth',
                    sprdif_min_aa_length=sprdif_min_aa_length,
                    sprdif_min_node_count=sprdif_min_node_count)
            sprdifprevnodes = sprdifprev.keys()
            sprdifnext = codingblock_splitting.spanningrange_difference(
                    next, 'left',
                    sprdif_min_aa_length=sprdif_min_aa_length,
                    sprdif_min_node_count=sprdif_min_node_count)
            sprdifnextnodes = sprdifnext.keys()
            mutualnodes = []
            sprdifmutual = {}

            for prevNode in sprdifprevnodes:
                if prevNode in sprdifnextnodes:
                    mutualnodes.append(prevNode)
                    sprdifrange = Set(sprdifprev[prevNode])
                    sprdifrange = sprdifrange.intersection(sprdifnext[prevNode])
                    sprdifmutual[prevNode] = sprdifrange

            ########################################################
            if verbose:
                print pos, "prev  :", [ (node,len(sprd)) for node,sprd in sprdifprev.iteritems() ]
                print pos, "next  :", [ (node,len(sprd)) for node,sprd in sprdifnext.iteritems() ]
                print pos, "mutual:", [ (node,len(sprd)) for node,sprd in sprdifmutual.iteritems() ]
            ########################################################

            if len(mutualnodes) < sprdif_min_node_count: continue
            if min([ len(vl) for vl in sprdifmutual.values()]) < sprdif_min_aa_length: continue
            # TODO impleemtn some other checks here!

            # gather sequence concerning the scaffold gap of the mutual nodes
            SCAFFOLD_GAP_OMSR_OFFSET = 0
            seqs, orfs, coords = {}, {}, {}
            for node in mutualnodes:
                org = prev.organism_by_node(node)
                sta = min( sprdifmutual[node] ) - SCAFFOLD_GAP_OMSR_OFFSET
                end = max( sprdifmutual[node] ) + SCAFFOLD_GAP_OMSR_OFFSET
                orf = prev.get_orfs_of_graph(organism=org)[0]
                seq = orf.getaas(abs_pos_start=sta,abs_pos_end=end)
                seqs[org]   = seq
                orfs[org]   = orf
                coords[org] = [sta,end]

            # do clustalw and strip_alignment_for_exterior_gaps
            (_algseqs,_algm) = clustalw(seqs=seqs)
            _algseqs,_algm,coords = strip_alignment_for_exterior_gaps(_algseqs,_algm,coords)

            # translate the clustalw alignment into an artificial CBG
            from graph_codingblock import CodingBlockGraph
            newcbg = CodingBlockGraph()
            newcbg.add_nodes(mutualnodes)
            for nodeA,nodeB in newcbg.pairwisecrosscombinations_node():
                orgA       = prev.organism_by_node(nodeA)
                orgB       = prev.organism_by_node(nodeB)
                # _algseqs keys are organisms, not nodes!
                alignment  = ( _algseqs[orgA], _algm, _algseqs[orgB] )
                paircoords = ( coords[orgA][0], coords[org][1], coords[orgB][0], coords[orgB][1] )
                pacbp = pacb.conversion.pacbp_from_clustalw(alignment=alignment,coords=paircoords)
                pacbporf = pacb.conversion.pacbp2pacbporf(pacbp,orfs[orgA],orfs[orgB])
                wt = pacbporf.bitscore
                pacbpkey = pacbporf.construct_unique_key(nodeA,nodeB)
                newcbg.add_edge(nodeA,nodeB,wt=wt)
                newcbg.pacbps[(pacbpkey,nodeA,nodeB)] = pacbporf
            newcbg.update_edge_weights_by_minimal_spanning_range()


            ########################################################
            if verbose:
                print "INTERMEDIATE HMM INPUT CBG:"
                print newcbg
                newcbg.printmultiplealignment()
            ########################################################

            # complete with cbghmmsearch2pacbpcollection
            pacbpCollection = cbghmmsearch2pacbpcollection(newcbg,self.input,
                    next=next,prev=prev,
                    hmmsearch_num_hits=3, # normal == 3
                    max_intron_nt_length=250, # IMPORTANT!!! keep this number rather low.
                                              # when a more biological realistic maximal length is choosen
                                              # e.g. >=500nt, this HMMcompleted cbg is likely to get abandoned
                                              # when by chance a wrong orf is linked to the HMMcbg
                    )


            # get list of accepted CBGs
            all_accepted_cbgs =  []
            for acc_cbg in conversion.pacbpCollection2AcceptedCodingBlockGraphs(pacbpCollection,prev=prev,next=next):
                acc_cbg.update_edge_weights_by_minimal_spanning_range()
                acc_cbg.create_cache()
                # check cexpander binarystring for all zeros.
                if do_cexander_allzeros_check and\
                acc_cbg._cexpander.binarystring.find("1") == -1:
                    # ignore all zeros here...
                    continue
                # if here, then this split2hmmcompletedCBG is `accepted`
                all_accepted_cbgs.append( acc_cbg)

            # order graphs by total weight and re-order on node occurrence:
            # if a neighboring node is incorporated -> more likely!
            all_accepted_cbgs = ordering.order_cbgraphlist(all_accepted_cbgs)
            all_accepted_cbgs = ordering.reorder_cbgs_on_node_occurrence(all_accepted_cbgs,prev=prev,next=next)

            ########################################################
            if verbose:
                print "PREV", prev
                print "scaf", newcbg
                print "NEXT", next
            ########################################################

            for intermediatecbg in all_accepted_cbgs:
                ratio = intermediatecbg.genetree().identity() / newcbg.genetree().identity()
                if verbose: print intermediatecbg, ratio
                if ratio < 0.70: continue
                status = self.add_codingblock(intermediatecbg)
                if status:
                    NEW_INTERMEDIATE_CBG_IS_CREATED = True

        if NEW_INTERMEDIATE_CBG_IS_CREATED:
            # recreate lsrCBGs in case newCBG is created
            self.create_intermediary_lsrcbgs()

        # return the Boolean weather or not something is created
        return NEW_INTERMEDIATE_CBG_IS_CREATED

    # end of function create_large_intermediate_cbg_for_scaffold_gap


    ########################################################################
    ### Functions for obtaining gff of codingblocks/genestructure        ###
    ########################################################################

    def get_codingblock_gff(self,organism,gff={},detailed=True):
        """ TODO: create functionality here """
        # IMPORTANT! Reset current self._gff_current_codingblocks[organism]
        self._gff_current_codingblocks[organism] = []

        # some variables needed
        gffdata = []
        gff['gclassandname'] = "%s %s" % ( gff['gclass'], gff['fref'] )
        # loop over the codingblocksgraphs
        for sg in self.codingblockgraphs:
            startCoord, endCoord = None, None

            # continue if this sg has status ignored
            if sg.IS_IGNORED:
                gffdata.append( () )
                continue

            if organism not in sg.organism_set():
                # TODO TODO TODO
                # FIX THIS! the fact that there are codingblocks in other organisms,
                # means we have to do **something** for this organism to.
                # What? -> ab initio gene prediction, not alignment based!
                # What? -> DNA alignment based to detect possible frameshifts
                # What? -> HMM searches of protein alignment
                gffdata.append( () )
                continue

            # get start coordinate of codingblock gff track
            try:
                (site,lStatus) = _get_cbg_left_site(sg,organism)
                if site: startCoord = site.pos
            except:
                print sg
                (site,lStatus) = _get_cbg_left_site(sg,organism)

            #if not site and sg.IS_FIRST:
            #    noATG = True            # no eligiable tss found although this is the first codingblock graph
            #elif not site:
            #    noACCEPTOR = True       # no eligiable acceptor site found
            #else:
            #    startCoord = site.pos    # just get the coordinate of the site

            # get end coordinate of codingblock gff track
            try:
                (site,rStatus) = _get_cbg_right_site(sg,organism)
                if site: endCoord = site.pos
            except:
                print sg
                (site,rStatus) = _get_cbg_right_site(sg,organism)

            #if not site and sg.IS_LAST:
            #    noTGA = True            # no eligiable stop codon found although this is the last codingblock graph
            #elif not site:
            #    noDONOR = True          # no eligiable donor site found
            #else:
            #    endCoord = site.pos     # just get the coordinate of the site

            # now see if the coordinates are resolved
            # 'Fragment'  is used when no site is found at all
            # 'Estimated' is used when the boundary is estimated
            # prefix is assigned to 3p/5p/5p3p to show which boundary is a problem
            fmethod = gff['fmethod']
            #fstrand = '+'
            #prefix  = ''
            _fscoremapper  = { False: 0, True: 1 }
            _prefixmapper  = { (True,True): '', (True,False): '3p', (False,True): '5p', (False,False): '5p3p' }
            _fstrandmapper = { (True,True): '+', (True,False): '+', (False,True): '-', (False,False): '+' }
            fscore  = _fscoremapper[ self._codingblock_prediction_status(sg,organism) ]
            prefix  = _prefixmapper[(lStatus,rStatus)]
            fstrand = _fstrandmapper[(lStatus,rStatus)]

            #if (lStatus,rStatus) == (True,True):
            #    prefix = ''
            #elif (lStatus,rStatus) == (True,False):
            #    prefix = '3p'
            #elif (lStatus,rStatus) == (False,True):
            #    prefix = '5p'
            #    fstrand = '-'
            #else:
            #    prefix = '5p3p'

            if not (startCoord or startCoord==0) or not endCoord:
                startCoord = min(sg.minimal_spanning_range(organism=organism)) * 3
                endCoord = max(sg.minimal_spanning_range(organism=organism)) * 3
                fmethod = 'Fragment5p3p'+fmethod
            elif not (startCoord or startCoord==0):
                startCoord = min(sg.minimal_spanning_range(organism=organism)) * 3
                fmethod = 'Fragment5p'+fmethod
                fstrand = '-'
            elif not endCoord:
                endCoord = max(sg.minimal_spanning_range(organism=organism)) * 3
                fmethod = 'Fragment3p'+fmethod
            else:
                if prefix:
                    fmethod = 'Estimated'+prefix+fmethod

            if detailed and sg.__class__.__name__ == 'LowSimilarityRegionCodingBlockGraph':
                fmethod = sg.__class__.__name__
                c9tv_txt = [ ( '; LowSimilarityRegionCodingBlockGraph', str(sg) ) ]
                # append the Orfs of this CodingBlock
                for org in sg.organism_set():
                    c9tv_txt.append( ('Orf', str(sg.get_orfs_of_graph(organism=org)[0]) ) )
                # append the boundary5p and boundar3p object
                boundary5p = _get_cbg_boundary5p_object(sg)
                boundary3p = _get_cbg_boundary3p_object(sg)
                c9tv_txt.append( ('Boundary5p',str(boundary5p) ) )
                c9tv_txt.append( ('Boundary3p',str(boundary3p) ) )
                c9tv_txt = "; ".join( [ "%s '%s'" % (a,b) for (a,b) in c9tv_txt ] )
                # replace '<' '>' symbols for '&lt;' and '&gt;' to enable proper html visualisation
                c9tv_txt = c9tv_txt.replace('<','[ ').replace('>',' ]')
            elif detailed:
                # make nicely layouted c9tv data for this gff line
                c9tv_txt = [ ( '; CodingBlockGraph', str(sg) ) ]
                c9tv_txt.append( ( 'AVERAGE-TCODE-OMSR', "%1.3f" % sg.msr_tcode_score() ) )
                # append the Orfs of this CodingBlock
                for org in sg.organism_set():
                    c9tv_txt.append( ('Orf', str(sg.get_orfs_of_graph(organism=org)[0]) ) )
                # append cbg prediction stati
                for org in sg.organism_set():
                    predstatus = self._codingblock_prediction_status(sg,org)
                    c9tv_txt.append( ('CbgAbfgpOptimalityStatus', "%s-%s" % (org,predstatus)))
                # append the PAcbPORFs of this CodingBlock
                for k,pacbporf in sg.pacbps.iteritems():
                    c9tv_txt.append( ('PacbpOrf',str(pacbporf) ) )
                # append the boundary5p and boundar3p object
                boundary5p = _get_cbg_boundary5p_object(sg)
                boundary3p = _get_cbg_boundary3p_object(sg)
                c9tv_txt.append( ('GeneTreeGraphCBG', sg.genetree() ) )
                c9tv_txt.append( ('GeneTreeGraph', self.genetree() ) )
                c9tv_txt.append( ('Boundary5p',str(boundary5p) ) )
                c9tv_txt.append( ('Boundary3p',str(boundary3p) ) )
                c9tv_txt = "; ".join( [ "%s '%s'" % (a,b) for (a,b) in c9tv_txt ] )
                # replace '<' '>' symbols for '&lt;' and '&gt;' to enable proper html visualisation
                c9tv_txt = c9tv_txt.replace('<','[ ').replace('>',' ]')
            else:
                c9tv_txt = ""

            # create gff line and append to gffdata
            col9 = "%s %s_%s_%s" % ( gff['gclass'], gff['fref'], startCoord+1, endCoord )
            gffline = ( gff['fref'], gff['fsource'], fmethod, startCoord+1, endCoord, fscore, fstrand, '.', "%s%s" % (col9,c9tv_txt)  )
            gffdata.append( gffline )

        # done; set to self._gff_current_codingblocks[organism]
        # and return the gffdata list
        self._gff_current_codingblocks[organism] = gffdata
        return gffdata

    # end of function get_codingblock_gff


    def sources(self):
        """
        """
        retd = {}
        for pacbp in self.pacbps.values():
            if retd.has_key(pacbp.source): retd[pacbp.source] += 1
            else:                          retd[pacbp.source] =  1
        return retd

    # end of function sources


    def extend_cbgs(self,verbose=False):
        """
        extend CBGs on partial aligned regions
        """
        for i in range(0,len(self)):
            cbg = self.codingblockgraphs[i]
            if cbg.IS_IGNORED: continue
            if cbg.__class__.__name__ == 'LowSimilarityRegionCodingBlockGraph': continue
            prev, next = None, None
            if i>0:           prev = self.codingblockgraphs[i-1]
            if i<len(self)-1: next = self.codingblockgraphs[i+1]
            if prev and len(prev.node_set().intersection(cbg.get_nodes())) == cbg.node_count():
                statusL = False
            elif cbg.IS_5P_SPLITTED:
                statusL = False
            else:
                if verbose: print "extCBG L:", cbg
                statusL = cbg.extend_on_left_spanningrange_difference(
                        sprdif_min_aa_length=CBG_EXTEND_LEFT_SPRDIF_MIN_AA_LENGTH,
                        verbose=verbose)
            if next and len(next.node_set().intersection(cbg.get_nodes())) == cbg.node_count():
                statusR = False
            elif cbg.IS_3P_SPLITTED:
                statusR = False
            else:
                if verbose: print "extCBG R:", cbg
                statusR = cbg.extend_on_rigth_spanningrange_difference(
                        sprdif_min_aa_length=CBG_EXTEND_RIGTH_SPRDIF_MIN_AA_LENGTH,
                        verbose=verbose)

    # end of function extend_cbgs


    def correct_lsrcbgs_after_optimization(self,verbose=False):
        """
        Recreate lsrCBGs in places where the OMSR is changed due to OMSR optimization

        @attention: this function is harmless in a sence that it can be called as often as desired

        @rtype:  Boolean 
        @return: True or False weather or not a lsrCBG is recreated
        """
        # return status boolean weather or not a lsrCBG is changed 
        RETURN_STATUS_LSRCBG_IS_CHANGED = False
            
        # loop BACKWARDS over the CBGs; not required (no inserts), but common in other functions too 
        for i in range(len(self)-2,0,-1):
            thiscbg = self.codingblockgraphs[i]
            if thiscbg.__class__.__name__ != 'LowSimilarityRegionCodingBlockGraph': continue 

            # get the surrounding CBGs
            prevcbg = self.codingblockgraphs[i-1]
            nextcbg = self.codingblockgraphs[i+1]

            # check omsr coordinate difference
            lsromsrstarts = thiscbg.omsr_starts()
            lsromsrends   = thiscbg.omsr_ends()
            for k in lsromsrstarts.keys(): lsromsrstarts[k]-=1
            for k in lsromsrends.keys(): lsromsrends[k]+=1

            # check if OSMR ends are identical
            if5Pstatus = prevcbg.omsr_ends() == lsromsrstarts
            if3Pstatus = nextcbg.omsr_starts() == lsromsrends 

            ########################################################################
            if verbose:
                print prevcbg
                print thiscbg
                print nextcbg
                print lsromsrstarts, if5Pstatus 
                print lsromsrends, if3Pstatus
            ########################################################################

            if not if5Pstatus or not if3Pstatus:
                # recreate the lsrCBG -> coords have changed!
                newlsrCBG = codingblock_splitting.create_intermediate_lowsimilarity_region(prevcbg,nextcbg)
                self.codingblockgraphs[i] = newlsrCBG
                # clear (potentially present) cbgIF objects
                prevcbg._forced_3p_ends = {}
                nextcbg._forced_5p_ends = {}
                prevcbg._CBGinterface3p = None
                nextcbg._CBGinterface5p = None
                # set the return boolean to True
                RETURN_STATUS_LSRCBG_IS_CHANGED = True

                ########################################################################
                if verbose:
                    # check omsr coordinate difference
                    lsromsrstarts = newlsrCBG.omsr_starts()
                    lsromsrends   = newlsrCBG.omsr_ends()
                    for k in lsromsrstarts.keys(): lsromsrstarts[k]-=1
                    for k in lsromsrends.keys(): lsromsrends[k]+=1
                    if5Pstatus = prevcbg.omsr_ends() == lsromsrstarts
                    if3Pstatus = nextcbg.omsr_starts() == lsromsrends 
                    print "# CHANGED!!!", newlsrCBG 
                    print "# changed:::", lsromsrstarts, if5Pstatus        
                    print "# changed:::", lsromsrends, if3Pstatus
                ########################################################################

        # return the boolean
        return RETURN_STATUS_LSRCBG_IS_CHANGED

    # end of function correct_lsrcbgs_after_optimization


    def join_false_inframe_introns(self):
        """
        Join (false) inframe introns by creating a lsrCBG in between CBGs with exactly the same nodes
        """
        # return status boolean weather or not a lsrCBG is added
        RETURN_STATUS_LSRCBG_IS_ADDED = False

        # loop BACKWARDS over the CBGs in case of an insert
        for i in range(len(self)-1,0,-1):
            # get combinations of 2 neighbouring CBGs
            (first,second) = self.codingblockgraphs[i-1:i+1]
            poscombi = (i-1,i)
            # ignore if one of them IS_IGNORED
            if True in [ cbg.IS_IGNORED for cbg in (first,second) ]:
                continue 
            # ignore if combination is created by a split
            if [ first.IS_3P_SPLITTED, second.IS_5P_SPLITTED ] == [True,True]:
                continue
            # ignore if one of them is a lsrCBG
            if 'LowSimilarityRegionCodingBlockGraph' in [ first.__class__.__name__, second.__class__.__name__ ]:
                continue
            # ignore if not all mutual nodes
            if first.node_set().symmetric_difference(second.get_nodes()):
                continue

            # If this point is reached, `first` and `second` are CBGs with
            # exactly the same nodes. There are 2 options possible:
            # H0: continious exons in all species, with a region of low similarity in between them
            # H1: non-continious exon in at least a single species

            # get some data on these neighboring CBGs that
            distances = first.distance_between_codingblocks(second)
            (nodes_seen, max_node_seen) = cbgs_identical_pacbp_analysis(first,second)

            # TMP HARD_SET to True!
            # TODO: solve this lateron; we want a special place where
            # hypothetical inframe introns are mapped.
            # at this point, by setting hard to True (and not having written the code to
            # scan the lsrCBGs yet, not a single inframe intron will be found again anymore! 

            lowsimilarity_region = True

            #if 0 in nodes_seen.values():
            #    lowsimilarity_region = False
            #if sum(max_node_seen.values()) - sum(nodes_seen.values()) > len(nodes_seen.keys()):
            #    lowsimilarity_region = False
            if max(distances.values()) < 0:
                lowsimilarity_region = False


            #### tmp printing
            #### tmp printing
            print "JOINED ANALYSES:", (i-1,i)
            print "DISTANCES:",distances
            print nodes_seen
            print max_node_seen
            print "lowsimilarity_region:", lowsimilarity_region
            #### tmp printing
            #### tmp printing

            if lowsimilarity_region:
                ##firstOMSR   = first.overall_minimal_spanning_range()
                ##secondOMSR  = second.overall_minimal_spanning_range()
                ##lsrCBG = LowSimilarityRegionCodingBlockGraph()
                ##for node in nodes_seen.keys():
                ##    lsomsr = range( max(firstOMSR[node])+1, min(secondOMSR[node]) )
                ##    if lsomsr: 
                ##        org = first.organism_by_node(node)
                ##        lsrCBG.add_node_and_object( node, first.get_orfs_of_graph(organism=org)[0] )
                ##        lsrCBG.set_node_omsr(node,Set(lsomsr)) 

                # create intermediate lsrCBG
                lsrCBG = codingblock_splitting.create_intermediate_lowsimilarity_region(first,second)

                # check if this lsrCBG got any nodes (lsromsr not added!)
                if lsrCBG.get_nodes():
                    print lsrCBG
                    print "potential inframe intron:", lsrCBG.potentially_contains_inframe_intron() 
                    # update the status of CBG first and second
                    first.IS_SPLITTED     = True
                    first.IS_3P_SPLITTED  = True
                    second.IS_SPLITTED    = True
                    second.IS_5P_SPLITTED = True
                    # insert the LowSimilarityRegionCodingBlockGraph on the proper position
                    self.codingblockgraphs.insert(i,lsrCBG)
                    RETURN_STATUS_LSRCBG_IS_ADDED = True
                else:
                    # hmmmm .... does not happen a lot, but it can happen
                    # 1) lowconnected sprdif function did not detect this similarity (or had not been performed)
                    # 2) large sprdif function did not detect this similarity either (or had not been performed)
                    # the result however is a seamless (or slighly overlapping) extension of the neighboring CBG
                    # without any gap in between this neighbour and the new detected one
                    print "MaiorWarning: in join_false_inframe_introns function: NO lsrCBG CREATED!!!"
                    print first
                    first.printmultiplealignment()
                    print second
                    second.printmultiplealignment()

        # return the status wetaher or not a lsrCBG is added
        return RETURN_STATUS_LSRCBG_IS_ADDED

    # end of function join_false_inframe_introns


    def get_genestructure_gff(self,organism=None,gff={}):
        """ TODO: create functionality here """
        # if no gff dictionary is given, use dummy default values
        if not gff: gff = { 'gclass': None, 'fref': None, 'fsource': None, 'fmethod': None }

        # check if gff of this codingblock is already precalculated
        if not self._gff_current_codingblocks[organism]:
            codingblocksgff = self.get_codingblock_gff(organism,gff=gff)
        else:
            codingblocksgff = self._gff_current_codingblocks[organism]

        # some variables needed
        gffdata = []
        coords = []
        gff['gclassandname'] = "%s %s" % ( gff['gclass'], gff['fref'] )
        # loop over the gene structure codingblock graphs
        prev_is_resolved = True
        for pos in range(0,len(self)):
            # get cbg and gffdata of codingblock
            cbg = self.codingblockgraphs[pos]
            if cbg.IS_IGNORED: continue  # skip IGNORED codingblocks

            if organism not in cbg.organism_set():
                # TODO TODO TODO
                # FIX THIS! the fact that there are codingblocks in other organisms,
                # means we have to do **something** for this organism to.
                # What? -> ab initio gene prediction, not alignment based!
                # What? -> DNA alignment based to detect possible frameshifts
                # What? -> HMM searches of protein alignment
                # In current implementation (april 2010) this can never happen,
                # but in furure implementations this might be a possibility again!
                continue

            cbggff = codingblocksgff[pos]

            # get the boundaries objects
            (lBound,lStatus) = _get_cbg_left_site(cbg,organism)
            (rBound,rStatus) = _get_cbg_right_site(cbg,organism)

            #print organism, "\t", pos, lBound,lStatus
            #print organism, "\t", pos, rBound,rStatus

            # boundaries to coordinates
            lCoord,rCoord = None,None
            if not lBound:  lCoord = cbggff[3]  # TODO!! erroneous coordinate
            else:           lCoord = lBound.pos + 1 # Python2GFF conversion
            if not rBound:  rCoord = cbggff[4]  # TODO!! erroneous coordinate
            else:           rCoord = rBound.pos

            if pos == 0:
                coords.append( [ lCoord, rCoord ] )
            else:
                if prev_is_resolved:
                    # now see if this is an intron position or not!
                    if lBound.__class__.__name__ == 'NoneType':
                        coords.append( [ lCoord, rCoord ] )
                    elif lBound.__class__.__name__ == 'SpliceAcceptor':
                        coords.append( [ lCoord, rCoord ] )
                    elif lBound.__class__.__name__ == 'ProjectedSpliceAcceptor':
                        coords[-1].extend( [ lCoord, rCoord ] )
                    else:
                        # what else !?!?!?
                        print "WARNING WARNING, unexpected class: ", lBound.__class__.__name__
                        coords.append( [ lCoord, rCoord ] )
                else:
                    # previous site was shit; solve this: TODO!
                    coords.append( [ lCoord, rCoord ] )

            # set the `prev_is_resolved` variable
            if rBound.__class__.__name__ in ['ProjectedSpliceDonor','SpliceDonor']:
                prev_is_resolved = True
            else:
                prev_is_resolved = False


        # get list of ABGP prediction outcome on each of the interfaces, S and E
        # this to make the final classifier (and layout the gene structure tracks
        ABGP_PREDICTION_LIST   = self.abgp_prediction_status(organism,verbose=False)
        if False in ABGP_PREDICTION_LIST:           ABGP_PREDICTION_STATUS = False
        elif ABGP_PREDICTION_LIST.count(None) <= 1: ABGP_PREDICTION_STATUS = True
        elif ABGP_PREDICTION_LIST.count(None) <= 2: ABGP_PREDICTION_STATUS = None
        else:                                       ABGP_PREDICTION_STATUS = False

        if ABGP_PREDICTION_STATUS == True:    fscore = 1
        elif ABGP_PREDICTION_STATUS == False: fscore = -1
        else:                                 fscore = 0

        # loop over the coord items
        for coordlist in coords:
            startCoord = min(coordlist) # already in GFF coordinate system!
            endCoord   = max(coordlist)
            # create gff line and append to gffdata
            gffline = ( gff['fref'], gff['fsource'], gff['fmethod'],
                        startCoord, endCoord, fscore, '+', '.',
                        gff['gclassandname']  )
            gffdata.append( gffline )

        # return the gffdata list
        return gffdata

    # end of function get_genestructure_gff


    def finalize_genestructure(self):
        """ Finalize; (re)set IS_FIRST and IS_LAST if not done properly """

        first_is_seen = False
        for pos in range(0,len(self)):
            # get current CBG
            cbg  = self.codingblockgraphs[pos]
            if cbg.IS_IGNORED:
                cbg.IS_FIRST = False
                cbg.IS_LAST  = False
            elif not first_is_seen:
                cbg.IS_FIRST  = True
                first_is_seen = True
            else:
                cbg.IS_FIRST = False

        last_is_seen  = False
        for pos in range(len(self)-1,-1,-1):
            # get current CBGs
            cbg  = self.codingblockgraphs[pos]
            if cbg.IS_IGNORED:
                cbg.IS_FIRST = False
                cbg.IS_LAST  = False
            elif not last_is_seen:
                cbg.IS_LAST  = True
                last_is_seen = True
            else:
                cbg.IS_LAST  = False

    # end of function finalize_genestructure


    def search_for_intermediary_tinyexon(self,similaritymatrix=None):
        """
        Search CBG junctions in genestructure for intermediary tiny exons
        """
        finding_new = True
        startat = 0
        while finding_new:
            for pos in range(startat,len(self)):
                sg = self.codingblockgraphs[pos]
                prev = None
                if pos >= 1: prev = self.codingblockgraphs[pos-1]
    
                # confirm that organism sets in both CBGs are identical (no missing organisms / partial CBGs)
                if prev and not prev.organism_set().symmetric_difference(sg.organism_set()):

                    # make DonorSiteCollectionGraph and AcceptorSiteCollectionGraph, NO projected sites included yet!!!
                    sg.harvest_elegiable_acceptor_sites(projected_acceptors={},prev=prev)
                    prev.harvest_elegiable_donor_sites(projected_donors={},next=sg)

                    tinyexonsg = find_intermediary_codingblockgraph_with_tinyexon(
                            prev,sg,input=self.input,
                            similaritymatrix=similaritymatrix)
                    if tinyexonsg:
                        # insert into genestructure graphs at the correct position
                        self.codingblockgraphs.insert( pos, tinyexonsg[0] )
    
                        # set startat pointer to current point in the list genestructure_graphs
                        # all previous sg's are scanned for tinyexons
                        startat = pos+1
                        break
            else:
                finding_new = False
    
    # end of search_for_intermediary_tinyexon


    def _codingblock_prediction_status(self,cbg,organism,verbose=False):
        """ """
        if cbg.__class__.__name__ == 'LowSimilarityRegionCodingBlockGraph': return True

        # check if cbgIFs are defined already. If not -> no status to be assigned yet!
        if cbg._CBGinterface5p and cbg._CBGinterface3p:
            pass
        elif cbg._CBGinterface5p and cbg.IS_LAST:
            pass
        elif cbg._CBGinterface3p and cbg.IS_FIRST:
            pass
        elif cbg.IS_FIRST and cbg.IS_LAST:
            pass
        else:
            # return None -> this is likely to result in Raised Erros lateron..
            # CBGs MUST fall in the above 4 categories in this step
            return None

        # get start & end coordinate of codingblock gff track
        try:
            (siteL,lStatus) = _get_cbg_left_site(cbg,organism)
            startCoord      = siteL.pos
        except:
            startCoord, lStatus = None, False
        try:
            (siteR,rStatus) = _get_cbg_right_site(cbg,organism)
            endCoord        = siteR.pos
        except:
            endCoord, rStatus = None, False 

        # TODO: these thresholds should be made omsrlength() dependant
        MAXIMAL_CBG_GFF2OMSR_MATCH_RATIO          = 0.35 # 35% of the alignment
        MAXIMAL_CBG_GFF2OMSR_NT_LENGTH_DIFFERENCE = 24   # 8 amino acids

        cbg_boundary_combinations = [
            ('CodingBlockStart',       'CodingBlockEnd'),
            ('TranslationalStartSite', 'CodingBlockEnd'),
            ('CodingBlockStart',       'StopCodon'),
        ]
        if startCoord == None or endCoord == None:
            # hmm... this is rather quick conclusion.
            # What in case of a (bona fide), very long CBG
            # but with a very poor cbgIF? Then it is penalized here...
            # TODO: make a better estimated judgement in this case
            status = False
        elif lStatus == True and rStatus == True and\
        (siteL.__class__.__name__,siteR.__class__.__name__) in\
        cbg_boundary_combinations:
            # a boundary that is restricted by a 'forced' site;
            # bordering lsrCBGs, splitted CBGs, Start- and StopCodons
            # Here, status is always (?) okay, leave the prediction
            # accurracy estimation to the cbgIFs and aligned starts/stops
            status = True
        else:
            # get omsr range of this cbg
            if  siteR.__class__.__name__ == 'StopCodon':
                omsr   = cbg.omsr2orfend(organism=organism) 
            else:
                omsr   = cbg.overall_minimal_spanning_range(organism=organism)
            ntomsr     = Set(range(min(omsr)*3,(max(omsr)+1)*3))
            gffomsr    = Set(range(startCoord,endCoord+1))
            ntomsrlen  = len(ntomsr)
            gffomsrlen = len(gffomsr)
            # compare sets -> big overlap required!
            difference = ntomsr.symmetric_difference(gffomsr)
            difflen    = len(difference)

            # get ratio's; check for ZeroDivisionError. These errors
            # can occur on very obscure CBGs where the start is located
            # AFTER the stop. As I said: very obscure CBGs....
            try: ratioA = float(difflen) / float(ntomsrlen)
            except ZeroDivisionError: ratioA = 0.0
            try: ratioB  = float(difflen) / float(gffomsrlen)
            except ZeroDivisionError: ratioB = 0.0

            ####################################################################
            if verbose:
                print cbg,organism
                print siteL,lStatus
                print siteR,rStatus
                # prevent possible ValueError: min() max() arg is an empty sequence
                if len(ntomsr)  == 0: ntomsr = [0]
                if len(gffomsr) == 0: gffomsr= [0] 
                print "OMSR:", min(ntomsr), max(ntomsr),
                print "GFF:", min(gffomsr), max(gffomsr),
                print "ratio: %1.3f" % max([ratioA,ratioB]),
                print "difflen:", difflen
            ####################################################################

            if max([ratioA,ratioB]) > MAXIMAL_CBG_GFF2OMSR_MATCH_RATIO and\
            difflen > MAXIMAL_CBG_GFF2OMSR_NT_LENGTH_DIFFERENCE:
                # append status False!
                status = False
            else:
                status = True

        # return status
        return status

    # end of function _codingblock_prediction_status 


    def abgp_prediction_status(self,organism,verbose=False):
        """ """
        checklist = []

        # check the first CBG
        checklist.append( self._codingblock_prediction_status(self.codingblockgraphs[0],organism) )

        # analyses of all the CBGInterface objects
        for cbg in self.codingblockgraphs[1:]:
            # asses the CBGInterface object
            cbgIF = cbg._CBGinterface5p 
            if cbgIF.is_optimal():
                checklist.append( True )
            elif cbgIF.is_forced_for_organism(organism):
                checklist.append( True )
            elif cbgIF.is_majority_of_donors_optimal():
                if cbgIF.is_optimal_donor(organism=organism):
                    checklist.append( True )
                elif cbgIF.is_donor_improveable(organism=organism):
                    checklist.append( None ) # hmm...there is an alternative! 
                else:
                    checklist.append( True ) # no site of this phase higher scoring
            elif cbgIF.is_lacking_organism_node(organism=organism):
                # no site(s) for this organism/gene
                checklist.append( False )
            elif not cbgIF._optimal_aligned_donor or not cbgIF._optimal_aligned_acceptor:
                # no optimal donor or acceptor assigned?
                # things can't get much worse than that ;-)
                checklist.append( False )
            elif cbgIF._optimal_aligned_donor.phase() != cbgIF._optimal_aligned_acceptor.phase():
                # no phase agreement -> the most shittyest cbgIF possible
                # TODO: still, there can be correctly predicted sites for
                # some of the organisms. This must be implemented
                checklist.append( False )
            elif cbgIF._optimal_aligned_donor.phase() == cbgIF._optimal_aligned_acceptor.phase():
                # there are sites for this organism/gene
                # and the phase sof the optimal sites are in agreement
                # Now check if the best phase for this organism is this phase
                phase = cbgIF._optimal_aligned_donor.phase() 
                if phase in [0,1,2] and cbgIF.is_compatible():
                    # identical phases and sites available for all organisms/genes
                    if cbgIF.is_highest_pssm_phase(phase) and cbgIF.is_highest_pssm_phase(phase,organism=organism):
                        checklist.append( True )
                    else:
                        checklist.append( False )
                elif phase in [0,1,2] and not cbgIF.is_compatible():
                    # identical phases but not a complete aligned site
                    # check here if it seems okay
                    donororgset = Set() 
                    acceporgset = Set()
                    donororgset.update(cbgIF._forced_3p_ends.keys())
                    donororgset.update(cbgIF._optimal_aligned_donor.organism_set())
                    acceporgset.update(cbgIF._forced_5p_ends.keys())
                    acceporgset.update(cbgIF._optimal_aligned_acceptor.organism_set())
                    if cbgIF.organism_set_size() == len(acceporgset) and len(acceporgset) == len(donororgset):
                        # forced sites and aligned sites together form an `aligned`  graph 
                        dsite = cbgIF._optimal_aligned_donor.get_organism_objects(organism)[0]
                        asite = cbgIF._optimal_aligned_acceptor.get_organism_objects(organism)[0]
                        if dsite.__class__.__name__ == 'ProjectedSpliceDonor' and asite.__class__.__name__ == 'ProjectedSpliceAcceptor':
                            # both projected sites -> seems a good alignment 
                            checklist.append( True )
                        elif dsite.__class__.__name__ == 'SpliceDonor' and asite.__class__.__name__ == 'SpliceAcceptor':
                            # seems compatible, again. Check if they represent the best pssm sites or not
                            if cbgIF.is_highest_pssm_phase(phase,organism=organism) and\
                            cbgIF.is_highest_pssm_donor(organism=organism) and\
                            cbgIF.is_highest_pssm_acceptor(organism=organism):
                                checklist.append( True )   # the best scoring pssm site
                            else:
                                checklist.append( False )  # not the best scoring pssm site
                        else:
                           # again, a big mess. No compatible sites.
                            checklist.append( False ) 
                    else:
                        # Nope, a big mess -> False for all!
                        checklist.append( False )
                else:
                    # more complicated case -> phase shift. Check identity of phase shift
                    # TODO: solve the pahse shift matching puzzle!
                    #print "PHASE SHOFT PHASE SHIFT", self
                    #for node, obj in cbgIF._optimal_aligned_donor._node_object.iteritems():
                    #    print node, obj
                    #for node, obj in cbgIF._optimal_aligned_acceptor._node_object.iteritems():
                    #    print node, obj
                    checklist.append( False )
            else:
                # ALL OTHER THINGIES NOT TRUE! no phase compatibility! But, that does not mean that all structures are False!
                print "ALL OTHER FAILED", cbfgIF
                checklist.append( False )

            # check the CBG itself
            checklist.append( self._codingblock_prediction_status(cbg,organism) )

        # analyses of the aligned stop codons
        checklist.append( self.get_final_cbg()._stopcodongraph.is_optimal(organism=organism) )

        if verbose: print self.get_final_cbg()._stopcodongraph

        # analyses of the aligned TSS / start codons
        firstcbg = self.get_first_cbg()
        if verbose: print firstcbg._startcodongraph

        if firstcbg._startcodongraph and firstcbg._startcodongraph.alignedsites:
            besttssgraph = firstcbg._startcodongraph.alignedsites[0]
            if besttssgraph.__class__.__name__ == 'TranslationalStartSiteCollectionGraph':
                checklist.insert( 0, False )
            else:
                besttssgraph._codingblockgraph = firstcbg # HARD-SET first cbg as codingblockgraph attribute
                checklist.insert( 0, besttssgraph.is_optimal(organism=organism) )
        else:
            checklist.insert( 0, False )
            print "NO ALIGNED _startcodongraph"

        return checklist

    # end of function abgp_prediction_status


    def obtain_noncbgif_pssm_objects(self,verbose=False):
        """
        Obtain (non-CBGInterface) PSSM objects in the GSG
        """
        # get the firstCBG for site scanning
        firstCBG = self.get_first_cbg()

        # make TranslationalStartSiteCollectionGraph
        firstCBG.harvest_elegiable_tss_sites()
        # do site alignment of TSS
        firstCBG._startcodongraph.collection2alignedsites()

        ########################################################################
        if verbose: print "B", firstCBG._startcodongraph
        ########################################################################

        # set cbg attribute to all aligned sites -> required for is_optimal check
        for atssgra in firstCBG._startcodongraph.alignedsites:
            atssgra._codingblockgraph = firstCBG

        ########################################################################
        if verbose: print "A", firstCBG._startcodongraph
        ########################################################################
        
        # make an acceptorsite collectiongraph for the first CBG
        # no Interface, no prev CBG, no projections, no forced ends!
        firstCBG.harvest_elegiable_acceptor_sites(
            projected_acceptors={},forced_codingblock_ends={},prev=None)
        # do site alignment of acceptors
        firstCBG._spliceacceptorgraph.collection2alignedsites()
        
        # align all the STOPcodon graphs (is hardly time-consuming)
        for cbg in self: cbg.align_stop_codons()

        # get the finalCBG for site scanning
        finalCBG = self.get_final_cbg()
        
        # make a donorsite collectiongraph for the final CBG
        # no Interface, no next CBG, no projections, no forced ends!
        finalCBG.harvest_elegiable_donor_sites(
            projected_donors={},forced_codingblock_ends={},next=None)
        # do site alignment of donors
        finalCBG._splicedonorgraph.collection2alignedsites()

    # end of function obtain_noncbgif_pssm_objects

# end of class GenestructureOfCodingBlockGraphs


################################################################################
#### Helper functions for GenestructureOfCodingblocks class                 ####
################################################################################



def _get_cbg_boundary5p_object(cbg):
    """ """
    if cbg.IS_FIRST:
        if cbg._startcodongraph.alignedsites:
            return cbg._startcodongraph.alignedsites[0]
        else:
            return "NO_ALIGNED_STARTCODONS!"
    else:
        if cbg._CBGinterface5p:
            if cbg._CBGinterface5p._forced_5p_ends:
                # forced_ends is always first choice; lsrCBG and related
                return cbg._CBGinterface5p._forced_5p_ends
            elif cbg._CBGinterface5p._optimal_aligned_acceptor:
                # take optimal acceptor when available
                return cbg._CBGinterface5p._optimal_aligned_acceptor
            elif cbg._CBGinterface5p._spliceacceptorgraph and\
            cbg._CBGinterface5p._spliceacceptorgraph.alignedsites:
                # final options: return the highest scoring aligned site
                return cbg._CBGinterface5p._spliceacceptorgraph.alignedsites[0]
            else:
                # shit... acceptorgraph or no aligned sites
                return "NO_ALIGNED_ACCEPTORS"
        else:
            # no cbgIF at all !?
            raise "AAAAP"

# end of function _get_cbg_boundary5p_object

            
def _get_cbg_boundary3p_object(cbg):
    """ """ 
    if cbg.IS_LAST:
        return cbg._stopcodongraph 
    else:   
        if cbg._CBGinterface3p:
            if cbg._CBGinterface3p._forced_3p_ends:
                # forced_ends is always first choice; lsrCBG and related
                return cbg._CBGinterface3p._forced_3p_ends
            elif cbg._CBGinterface3p._optimal_aligned_donor:
                # take optimal donor when available
                return cbg._CBGinterface3p._optimal_aligned_donor
            elif cbg._CBGinterface3p._splicedonorgraph and\
            cbg._CBGinterface3p._splicedonorgraph.alignedsites:
                # final options: return the highest scoring aligned site
                return cbg._CBGinterface3p._splicedonorgraph.alignedsites[0]
            else:
                # shit... no donorgraph or no aligned sites
                return "NO_ALIGNED_DONORS"
        else:
            # no cbgIF at all !?
            raise "AAAAP"

# end of function _get_cbg_boundary3p_object


#def _get_cbg_boundary5p_object(cbg):
#    """ """
#    if cbg.IS_FIRST:
#        if cbg._startcodongraph.alignedsites:
#            return cbg._startcodongraph.alignedsites[0]
#        else:
#            return "NO_ALIGNED_STARTCODONS!"
#    else:
#        if cbg._forced_5p_ends:
#            return cbg._forced_5p_ends
#        elif cbg._spliceacceptorgraph.alignedsites:
#            return cbg._spliceacceptorgraph.alignedsites[0]
#        else:
#            return "NO_ALIGNED_ACCEPTORS"
#
## end of function _get_cbg_boundary5p_object
#
#
#def _get_cbg_boundary3p_object(cbg):
#    """ """
#    if cbg.IS_LAST:
#        return cbg._stopcodongraph
#    else:
#        if cbg._forced_3p_ends:
#            return cbg._forced_3p_ends
#        elif cbg._splicedonorgraph.alignedsites:
#            return cbg._splicedonorgraph.alignedsites[0]
#        else:
#            return "NO_ALIGNED_DONORS"
#
## end of function _get_cbg_boundary3p_object


def _get_cbg_left_site(cbg,organism):
    """
    Return tuple is of shape (site, Boolean)
        site is a PSSM/GFF site or None
        Boolean is True or False, weather or not the site is obtained from a True alignment or a False estimation
    """
    if cbg.IS_FIRST:
        if not cbg._startcodongraph:
            # no startcodon graph set yet!
            # create an OMSR based CodingBlockStart object
            minntomsr = min(cbg.overall_minimal_spanning_range(organism=organism))*3
            cbgStart = CodingBlockStart(minntomsr)
            return ( cbgStart, False )
        elif cbg._startcodongraph.alignedsites and\
        cbg._startcodongraph.alignedsites[0].__class__.__name__ ==\
        'TranslationalStartSiteCollectionGraph':
            # the best and only `aligned` site is the collection it self
            if organism in cbg._startcodongraph.organism_set():
                # return the 'optimal' site
                tss = cbg._startcodongraph.get_optimal_single_site(organism)
                return ( tss, False )
            else:
                # hmm.. .we can not find an eligiable tss although this is the first codingblock graph
                return ( None, False )
        elif cbg._startcodongraph.alignedsites and organism in cbg._startcodongraph.alignedsites[0].organism_set():
            tss = cbg._startcodongraph.alignedsites[0].get_organism_objects(organism)[0]
            return ( tss, True )
        elif cbg._startcodongraph.alignedsites and organism in cbg._startcodongraph.organism_set():
            # return the 'optimal' site, given an optimal aligned site
            # where this organism is missing from
            tss = cbg._startcodongraph.get_optimal_single_site(organism)
            return ( tss, False )
        elif cbg._startcodongraph and organism in cbg._startcodongraph.organism_set():
            # return the 'optimal' site
            tss = cbg._startcodongraph.get_optimal_single_site(organism)
            return ( tss, False )
        else:
            # hmm.. .we can not find an eligiable tss although this is the first codingblock graph
            return ( None, False )
    elif not cbg._CBGinterface5p:
        # hmm. Probably the first CBG that is not labeled with IS_FIRST yet...
        # create an OMSR based CodingBlockStart object AND status False
        pass 
    else:
        # get the acceptor object from the CBGInterface
        return cbg._CBGinterface5p.get_acceptor_object(organism)

# end of function _get_cbg_left_site


def _get_cbg_right_site(cbg,organism):
    """
    Return tuple is of shape (site, Boolean)
        site is a PSSM/GFF site or None
        Boolean is True or False, weather or not the site is obtained from a True alignment or a False estimation
    """
    if cbg.IS_LAST:
        tgaNode = cbg._stopcodongraph.get_organism_nodes(organism)[0]
        (_org,_orf,stopQpos,stopQdnapos) = tgaNode
        tgaSite = StopCodon( stopQdnapos )
        return ( tgaSite, True )
    else:
        # get the donor object from the CBGInterface
        return cbg._CBGinterface3p.get_donor_object(organism)

# end of function _get_cbg_right_site


def mutualorgsofcbgs2partialGSG(input={},cbgs=[],verbose=False):
    """
    """
    ####################################################
    if verbose:
        stw = StopWatch( name='mutualorgsofcbgs2partialGSG' )
        stw.start()
    ####################################################
    allmutualorgs = cbgs[0].organism_set()
    removenodes = [ [] ]
    for i in range(1,len(cbgs)):
        cbg = cbgs[i]
        allmutualorgs.intersection_update(cbg.organism_set())
        removenodes.append( [] )


    ####################################################
    if verbose:
        print stw.lap(), "allmutualorgs:", allmutualorgs
    ####################################################

    for i in range(0,len(cbgs)):
        cbg = cbgs[i]
        for node in cbg.get_nodes():
            if cbg.organism_by_node(node) not in allmutualorgs:
                removenodes[i].append( node )

    ####################################################
    if verbose:
        print stw.lap() 
        for i in range(0,len(cbgs)):
            cbg = cbgs[i]
            print cbg
            print removenodes[i]
    ####################################################

    gsgcbgslist = []
    for i in range(0,len(cbgs)):
        cbg = cbgs[i]
        removenodelist = removenodes[i]
        if removenodelist:
            newCBG = codingblock_operations.deepcopy_with_removed_nodes(cbg,removenodelist)
        else:
            newCBG = cbg
        gsgcbgslist.append(newCBG)

    ####################################################
    if verbose:
        print stw.lap() 
        for i in range(0,len(cbgs)):
            newCBG = gsgcbgslist[i]
            print newCBG
    ####################################################


    mutualorgsinput = {}
    for org,datastruct in input.iteritems():
        if org in allmutualorgs:
            mutualorgsinput[org] = datastruct


    # now create a partialGSG
    partialGSG = GenestructureOfCodingBlockGraphs(mutualorgsinput)
    partialGSG.codingblockgraphs = gsgcbgslist
    ####################################################
    if verbose: print stw.lap(), partialGSG
    ####################################################
    return partialGSG

# end of function mutualorgsofcbgs2partialGSG


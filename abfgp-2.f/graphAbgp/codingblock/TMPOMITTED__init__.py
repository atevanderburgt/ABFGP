"""
################################################################################
### CodingBlockGraph class                                                  ####
################################################################################

Several size properties of a CodingBlockGraph are introduced here:

+ Minimal Spanning Range (MSR)

+ Overall Minimal Spanning Range (OMSR)

    The AA coordinate range of a CBG that is covered by all pairwise alignment
    objects (PacbPORFS). Each node/organism in a CBG has its own OMSR. The
    OMSR function only supports CBGs with PacbPORF objects, not PacbP objects.
    The OMSR is supposed to correspond to a large extend to the actual
    intron-exon boundary which is present in at least a single node/organism
    in the group of genes. MSR and OMSR can be identical, but mostly the OMSR
    will be a smaller subset of the MSR.

    An example of the applicability of the OMSR is in cases of large
    evolutionary distances in a group of genes, of which a few are more closely
    related. Highly similar genes/proteins have high(er) conservation rates
    in the intronic regions, and as a consequence their PacbP(ORF)s are likely
    to run into the intronic region. In PacbP(ORF)s between species with larger
    evolutionary distances, the PacbP(ORF)s are more likely to stop nearby the
    actual intron-exon boundary. The coordinate range of these more distantly
    related genes will yield an OMSR which indeed corresponds to the exonic
    part of the genes.

+ Maximal Spanning Range (MAXSPR)

    The largest AA coordinate range of a CBG, that is covered by at least a
    single pairwise alignment object (PacbPORFS). Each node/organism in a CBG
    has its own MAXSR. The MAXSPR function supports both PacbPORF and
    PacbP objects. Obviously, MSR, OMSR and SPRDIF are subsets of the MAXSPR.

+ Spanning Range Difference (SPRDIF)

    The difference between the MAXSR and the OMSR, defined on one side of
    the CBG (left/5p side or rigth/3p side). Each node/organism in a CBG has
    its own spanningrange difference. The SPRDIF functions (left/rigth) only
    supports CBGs with PacbPORF objects, not PacbP objects. A SPRDIF does
    typically not exist for all nodes/organisms in the CBG.

    An example of the applicability of the SPRDIF is in cases where at least
    2 genes lack a small/tiny exon which is present in one or more other genes.
    In this case, a SPRDIF is likely to exist for the nodes/organisms which
    lack the small exon. This SPRDIF can be used to construct a HMMER profile
    which can be used to search in the Orfs in the other group of genes for
    presence of a small/tiny exon sequence that corresponds to the SPRDIF.





+ Minimal Spanning Range (MINSPR)


"""

# Module metadata
__authors__ = "Ate van der Burgt"
__license__ = "MIT"


# graphAbgp imports
from graphAbgp.graph_pacbpcollection import PacbpCollectionGraph
import codingblock_printing
import splitting as codingblock_splitting
import extension as codingblock_extension
import codingblock_optimization
import codingblock_collectionharvesting as cbg2colgra
import codingblock_genetree
import codingblock_codingperiodicity
import codingblock_sequenceretrieval
import conversion
from exceptions import *
from graph_alignedsites import AlignedStopCodonGraph


# Pacb imports
import pacb

# Abgp imports
from lib_clustalw import clustalw
from lib_cexpander import cexpanderanalyses
from lib_crossdata import remove_alternatives_from_pacbps_dict

# Python imports
from sets import Set
from copy import deepcopy

# Import Global variables and settings
from settings.codingblockgraph import *
from settings.pacbp import ALTERNATIVE_ALIGNMENT_OVERLAP_RATIO
from settings.sitealignment import (
    ALIGNED_DONOR_MAX_TRIPLET_DISTANCE,
    ALIGNED_ACCEPTOR_MAX_TRIPLET_DISTANCE
    )
from settings.splicesites import (
    DONORPSSMOPTIONS,
    ACCEPTORPSSMOPTIONS,
    )


class CodingBlockGraph(PacbpCollectionGraph,
    codingblock_printing.CodingBlockGraphPrintingFunctions,
    codingblock_genetree.GeneTreeOfCodingBlockFunctions,
    codingblock_sequenceretrieval.CodingBlockGraphSequenceRetievalFunctions):
    """
    CodingBlockGraph (CBG) class, inheriting from graphPlus class.
    """

    def __init__(self,
        minimal_overal_spanning_range_size=\
                CBG_MINIMAL_OVERAL_SPANNING_RANGE_SIZE,
        alternative_alignment_overlap_ratio=\
                ALTERNATIVE_ALIGNMENT_OVERLAP_RATIO):
        """
        Initialize a CodingBlockGraph
        """
        # Initialize from PacbpCollectionGraph
        PacbpCollectionGraph.__init__(self)

        # short name tag for printing
        self._short_name = "CBG"

        # lists for sites
        self._stopcodongraph = None
        self._startcodongraph = None
        self._splicedonorgraph = None
        self._spliceacceptorgraph = None

        # some GLOBAL variables that are required in class functions
        self.ALTERNATIVE_ALIGNMENT_OVERLAP_RATIO =\
                alternative_alignment_overlap_ratio

        # 08/04/09 smallest omsr size must be 3 (means 4 AAs overall aligned)
        # this is defined in settings.codingblockgraph
        self.MINIMAL_OVERAL_SPANNING_RANGE_SIZE =\
                minimal_overal_spanning_range_size 

        # attribute needed to set the GeneTreeGraph into
        self._GENETREE      = None

        # attributes to cache OMSR and related data structures
        self._omsr  = {}    # Overall Minimal Spanning Range
        self._msr   = {}    # Minimal Spanning Range
        self._maxsr = {}    # MAXimal Spanning Range

        # attributes to cache other CBG properties
        self._tcode         = None
        self._cexpander     = None

        # attributes that are needed once added into a genestructure
        self.IS_FIRST       = False
        self.IS_LAST        = False
        self.IS_IGNORED     = False
        self.IS_SPLITTED    = False
        self.IS_5P_SPLITTED = False
        self.IS_3P_SPLITTED = False

        # attributes that are filled when this CBG is neighboring a
        # LowSimilarityRegionCodingBlockGraph
        self._forced_3p_ends = {}
        self._forced_5p_ends = {}

        # attributes for CodingBlockGraphInterface objects
        self._CBGinterface3p = None
        self._CBGinterface5p = None

        # THESE ATTRIBUTES WILL BE DEPRECATED SOON!!!!
        # NOT USED ANYMORE AND VERY TIME CONSUMING.
        # dicts for data of VISTA like tracks
        self._paoc_per_organism = {} # positional_alignment_occurrence_count
        self._pasc_per_organism = {} # positional_alignment_similarity_count

    # end of function __init__


    def __str__(self):
        """ """
        organdorfs = self.get_ordered_nodes()
        tcode    = ""
        identity = " " # single space!
        if self.has_overall_minimal_spanning_range():
            allmsr = self.overall_minimal_spanning_range()

            # TODO: this happens rarely, but e.g. in `example6` full mode
            # TODO: although OMSR is reported, allmsr lacks (at least) 1 key!
            if len(allmsr.keys()) != len(organdorfs):
                msrlen = min( self.minimal_spanning_range_sizes().values() )

                for k,pacbporf in self.pacbps.iteritems():
                    print k, pacbporf

                themsr = [ "??OMSR??" ]
                themsr.extend( organdorfs )
                themsr.extend( [ "".join( ( org,"(",str(orf),")",":",str(min(allmsr[(org,orf)])),"..",str(max(allmsr[(org,orf)])) ) ) for (org,orf) in allmsr.keys() ] )
                return "<FCcbg?? wt=%s len=%s cs=%1.1f [%s]>" % (
                        self.total_weight(),
                        msrlen,
                        self.connectivitysaturation(),
                        " ".join([str(elem) for elem in themsr])
                        )
                # TODO: this piece of code has to be solved!!

            themsr = [ "OMSR" ]
            themsr.extend( [ "".join( ( org,"(",str(orf),")",":",str(min(allmsr[(org,orf)])),"..",str(max(allmsr[(org,orf)])) ) ) for (org,orf) in organdorfs ] )
            msrlen   = min( self.overall_minimal_spanning_range_sizes().values() )
            tcode    = " TC=%1.2f "  % self.msr_tcode_score()
            identity = " id=%1.2f " % self.omsr_identityscore()

        elif self.has_minimal_spanning_range():
            allmsr = self.minimal_spanning_range()
            themsr = [ "MSR" ]
            themsr.extend( [ "".join( ( org,"(",str(orf),")",":",str(min(allmsr[(org,orf)])),"..",str(max(allmsr[(org,orf)])) ) ) for (org,orf) in organdorfs ] )
            msrlen = min( self.minimal_spanning_range_sizes().values() )
        else:
            msrlen = 0
            themsr = ["NO OMSR"]
            themsr.extend( [ "".join( ( org,"(",str(orf),")" ) ) for (org,orf) in organdorfs ] )

        # calculate pacbp overshoot
        pacbp_overshoot = " "
        if len(self.pacbps) - self.edge_count() :
            pacbp_overshoot = " pacbp:+%s " % ( len(self.pacbps) - self.edge_count() )

        # return the string representation of this CBG
        return "<%s wt=%s%slen=%s cs=%1.1f%s%s[%s]>" % (
                self._short_name,
                self.total_weight(),
                identity,
                msrlen,
                self.connectivitysaturation(),
                tcode,
                pacbp_overshoot,
                " ".join(themsr)
                )
    # end of function __str__

    ############################################################################
    #### Overwritten functions from PacbpCollectionGraph / OrganismGraph    ####
    ############################################################################

    def _update_after_changes(self):
        """
        Update hidden/cached data in this CBG after a change in its PacbP(ORF)s
        
        @attention: OVERWRITES PacbpCollectionGraph._update_after_changes()
        """
        self._update_pacbps()
        self.clear_cache()

    # end of function _update_after_changes


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

    ############################################################################
    #### Functions concerning caching / copying the CBG                     ####
    ############################################################################

    def deepcopy(self):
        """
        Create a QUICK deepcopy of this CBG

        @rtype:  CodingBlockGraph
        @return: ``deepcopied`` CodingBlockGraph

        @attention: no actual deepcopy performed on the PacbP(ORF)s !!
        """
        # initialize a novel CBG instance
        dpcpCBG = CodingBlockGraph(
                minimal_overal_spanning_range_size=\
                 self.MINIMAL_OVERAL_SPANNING_RANGE_SIZE,
                alternative_alignment_overlap_ratio=\
                 self.ALTERNATIVE_ALIGNMENT_OVERLAP_RATIO
                )

        # copy the nodes
        dpcpCBG.add_nodes(self.get_nodes())

        # copy the edges
        for (node1,node2) in dpcpCBG.pairwisecrosscombinations_node():
            if self.has_edge(node1,node2):
                # copy edge from input CBG gra
                dpcpCBG.add_edge(node1,node2,wt=self.weights[(node1,node2)])

        # copy the pacbps
        for k,pacbp in self.pacbps.iteritems(): dpcpCBG.pacbps[k] = pacbp

        # return the ``deepcopied`` CBG
        return dpcpCBG

    # end of function deepcopy


    def clear_cache(self):
        """
        Cleanup the CBGs cached attributes (after a CBG change or update)
        """
        self._paoc_per_organism = {}
        self._pasc_per_organism = {}
        self._GENETREE = None
        self._omsr  = {}
        self._msr   = {}
        self._maxsr = {}

    # end of function clear_cache


    def create_cache(self):
        """
        Create cached attributes of this CodingBlockGraph (saves time)
        """
        self.clear_cache()
        self._msr   = self.minimal_spanning_range()
        self._maxsr = self.maximal_spanning_range()
        self._omsr  = self.overall_minimal_spanning_range()
        self.cexpanderanalyses()
        gt = self.set_genetree()
        
    # end of function create_cache

    ############################################################################
    #### Functions concerning the CodingBlockGraphInterfaces of the CBG     ####
    ############################################################################

    def copy_5pcbginterface_from_othercbg(self,othercbg):
        """
        Copy all 5' cbgIF properties from another CBG to this CBG

        @type  othercbg: CodingBlockGraph 
        @param othercbg: CodingBlockGraph 
        """
        self._CBGinterface5p  = othercbg._CBGinterface5p
        self._startcodongraph = othercbg._startcodongraph
        self.IS_5P_SPLITTED   = othercbg.IS_5P_SPLITTED
        self.IS_FIRST         = othercbg.IS_FIRST
        self._forced_5p_ends  = othercbg._forced_5p_ends
        if self.IS_5P_SPLITTED: self.IS_SPLITTED = True

    # end of function copy_5pcbginterface_from_othercbg 


    def copy_3pcbginterface_from_othercbg(self,othercbg):
        """
        Copy all 3' cbgIF properties from another CBG to this CBG

        @type  othercbg: CodingBlockGraph
        @param othercbg: CodingBlockGraph
        """
        self._CBGinterface3p = othercbg._CBGinterface3p
        self._stopcodongraph = othercbg._stopcodongraph
        self.IS_3P_SPLITTED  = othercbg.IS_3P_SPLITTED
        self.IS_LAST         = othercbg.IS_LAST
        self._forced_3p_ends = othercbg._forced_3p_ends
        if self.IS_3P_SPLITTED: self.IS_SPLITTED = True

    # end of function copy_3pcbginterface_from_othercbg


    def clear_donor_cbginterface(self):
        """
        Clear the CodingBlockGraphInterface on the 3' / right / donor side
        """
        self._CBGinterface3p = None
        self._forced_3p_ends = {}
        self._splicedonorgraph = None

    # end of function clear_donor_cbginterface


    def clear_acceptor_cbginterface(self):
        """
        Clear the CodingBlockGraphInterface on the 5' / left / acceptor side
        """
        self._CBGinterface5p = None
        self._forced_5p_ends = {}
        self._spliceacceptorgraph = None

    # end of function clear_acceptor_cbginterface


    def minimal_eligable_donor_site_position(self,organism):
        """
        Get the AA position for this Organism Identifier on this
        CodingBlockGraph from where to take donor positions into account

        @attention: See codingblock_collectionharvesting.
                    minimal_eligable_donor_site_position for documentation
        """
        return cbg2colgra.minimal_eligable_donor_site_position(
                self,organism)

    # end of function minimal_eligable_donor_site_position


    def maximal_eligable_acceptor_site_position(self,organism):
        """
        Get the AA position for this Organism Identifier on this
        CodingBlockGraph up to where to take acceptor positions into account

        @attention: See codingblock_collectionharvesting.
                    maximal_eligable_acceptor_site_position for documentation
        """
        return cbg2colgra.maximal_eligable_acceptor_site_position(
                self,organism)

    # end of function maximal_eligable_acceptor_site_position


    def maximal_eligable_donor_site_position(self,organism,nextcbg=None):
        """
        Get the AA position for this Organism Identifier on this
        CodingBlockGraph untill where to take donor positions into account

        @attention: See codingblock_collectionharvesting.
                    maximal_eligable_donor_site_position for documentation
        """
        return cbg2colgra.maximal_eligable_donor_site_position(
                self,organism,nextcbg=nextcbg)

    # end of function maximal_eligable_donor_site_position


    def minimal_eligable_acceptor_site_position(self,organism,prevcbg=None):
        """
        Get the AA position for this Organism Identifier on this
        CodingBlockGraph from where to take acceptor positions into account

        @attention: See codingblock_collectionharvesting.
                    minimal_eligable_acceptor_site_position for documentation
        """
        return cbg2colgra.minimal_eligable_acceptor_site_position(
                self,organism,prevcbg=prevcbg)

    # end of function minimal_eligable_acceptor_site_position


    ############################################################################
    #### Varia of uncategorized functions                                   ####
    ############################################################################

    def sources(self):
        """
        Get a dict with counts of sources of PacbP(ORF)s

        @rtype:  {} 
        @return: dict with sources (keys) and occurence counts (values)
                 of pacbp.source of PacbP(ORF)s 

        @attention: see pacbp/__init__.py for recognized sources, currently:
                    blastp, unknown, hmmsearch, clustalw,
                    clustalw-EXTENDED, clustalw-OPTIMIZED, lsrPACBP
        """
        ret = {}
        for pacbp in self.pacbps.values():
            if ret.has_key(pacbp.source):
                ret[pacbp.source]+=1
            else:
                ret[pacbp.source] =1
        return ret

    # end of function sources


    def update_edge_weights_by_minimal_spanning_range(self):
        """
        Update the weigths of the edges by the bitscore of the OMSR region only

        @attention: see the explanation of Overall Minimal Spanning Range
        """
        nodes_done = []
        for (org,orf), omsr in self.overall_minimal_spanning_range().iteritems():
            try:
                min_omsr = min(omsr)
                max_omsr = max(omsr) + 1
            except:
                # Exception occurred in min/max of OMSR ->
                # print error to STDOUT and than let the exception be raised
                print "UEWBMSP", self
                print "UEWBMSP", self.overall_minimal_spanning_range_sizes()
                min_omsr = min(omsr)
                max_omsr = max(omsr) + 1
            # loop over the pacbps
            for ( (a,b,c,d),(g1,o1),(g2,o2) ), pacbporf in self.pacbps.iteritems():
                if org == g1:
                    bitscore = pacbporf.bitscore_slice_by_abs_protein_query(
                                    min_omsr,max_omsr)
                elif org == g2:
                    bitscore = pacbporf.bitscore_slice_by_abs_protein_sbjct(
                                    min_omsr,max_omsr)
                else:
                    # what else? That is a case that can not happen...
                    continue
                
                # Now see if this node is done already.
                # If so, set the bitscore only if it is higher as the previous.
                # In practice AxB != BxA in terms of the alignment due to
                # (small) gaps. This can however cause quit big differences
                # in bitscore calculation. In stead of taking only the first
                # of a cross combination, take maximum bitscore of both.
                if not ( (g1,o1),(g2,o2) ) in nodes_done:
                    self.set_edge_weight( (g1,o1), (g2,o2), bitscore )
                    # append to done nodes
                    nodes_done.append( ( (g1,o1),(g2,o2) ) )
                    nodes_done.append( ( (g2,o2),(g1,o1) ) )
                else:
                    if bitscore > self.get_edge_weight( (g1,o1), (g2,o2) ):
                        self.set_edge_weight( (g1,o1), (g2,o2), bitscore )
        
    # end of function update_edge_weights_by_minimal_spanning_range


    def organisms_with_different_orfs(self,other):
        """
        Get a list of Organism identifiers with a distinct Orf in the other CBG

        @type  other: CodingBlockGraph
        @param other: CodingBlockGraph

        @rtype:  []
        @return: list with Organism identifiers which have distinct Orfs
        """
        return [ self.organism_by_node(node) for node in self.different_nodes(other) ]
    # end of function organisms_with_different_orfs 


    def cexpanderanalyses(self,**kwargs):
        """
        @attention: see lib_cexpander.cexpanderanalyses for function documentation
        """
        self._cexpander = cexpanderanalyses(self,**kwargs)

    # end of function cexpanderanalyses

    ############################################################################
    #### Functions fr detection of inframe-introns in this CodingBlockGraph ####
    ############################################################################

    def split_codingblock_by_inframe_intron(self,
        gap_size=INFRAME_INTRON_IN_CBG_GAP_SIZE,
        length_discrepancy=INFRAME_INTRON_IN_CBG_LENGTH_DISCREPANCY,
        gap_size_alone=INFRAME_INTRON_IN_CBG_GAP_SIZE_ALONE):
        """
        Split CodingBlockGraph on a potential inframe intron

        @attention: See codingblock_splitting.
                    split_codingblock_by_inframe_intron() for documentation

        @rtype:  []
        @return: Ordered!! list of 2 (succesfull) or 1
                 (non succesfull, unchanged) CodingblockGraphs
        """
        return codingblock_splitting.split_codingblock_by_inframe_intron(
                self,gap_size=gap_size,length_discrepancy=length_discrepancy,
                gap_size_alone=gap_size_alone )

    # end of function split_codingblock_by_inframe_intron


    def potentially_contains_aligned_intron(self,
        gap_size=INFRAME_INTRON_IN_CBG_GAP_SIZE,
        length_discrepancy=INFRAME_INTRON_IN_CBG_LENGTH_DISCREPANCY,
        window_aa_size=INFRAME_INTRON_IN_CBG_WINDOW_AA_SIZE,
        intron_window_min_similarity_score=\
                INFRAME_INTRON_IN_CBG_WINDOW_MIN_SIMILARITY_SCORE,
        dif_potentially_observed_vs_expected=\
                INFRAME_INTRON_IN_CBG_OBSERVED_VS_EXPECTED,
        MIN_TOTAL_PSSM_INFRAME_INTRON=INFRAME_INTRON_IN_CBG_MIN_TOTAL_PSSM):
        """
        Check if this CodingBlockGraph potentially contains an inframe intron

        @rtype:  Boolean
        @return: True or False

        @attention: window_aa_size is DEPRECATED!!!
        @attention: See codingblock_splitting.
                    potentially_contains_aligned_intron() for documentation
        """
        return codingblock_splitting.potentially_contains_aligned_intron(
                self,gap_size=gap_size,length_discrepancy=length_discrepancy,
                intron_window_min_similarity_score=\
                        intron_window_min_similarity_score,
                dif_potentially_observed_vs_expected=\
                        dif_potentially_observed_vs_expected,
                MIN_TOTAL_PSSM_INFRAME_INTRON=MIN_TOTAL_PSSM_INFRAME_INTRON
                )

    # end of function potentially_contains_aligned_intron

    ############################################################################
    #### Functions concerning the Orfs of this CodingBlockGraph             ####
    ############################################################################

    def orfenddistance(self,side,organism=None,node=None):
        """
        Return the distance between the OMSR and the exterior of the ORFs

        @type  node: *
        @param node: Node identifier (or None)

        @type  organism: *
        @param organism: Organism identifier (or None)

        @rtype:  {} or integer 
        @return: dict with (AA) distances or single integer when organism or node was specified 
        """
        if node: organism = self.organism_by_node(node)
        if organism and organism not in self.organism_set():
            raise OrganismNotPresentInGraph 
        omsr = self.overall_minimal_spanning_range()
        distances = {}
        for thenode in self.get_nodes():
            org = self.organism_by_node(thenode)
            if organism and org!=organism: continue
            theorf = self.get_orfs_of_graph(organism=org)[0]
            if side in ['left','5p']:
                distances[org] = min(omsr[thenode]) - theorf.protein_startPY
            elif side in ['right','3p']:
                distances[org] = theorf.protein_endPY - ( max(omsr[thenode]) + 1 )
            else:
                message = "argument 'side' not in 'left','5p','right','3p' but: '%s'" % side
                raise InproperlyAppliedArgument, message 

        # return single value or dict with distances
        if organism: return distances[organism]
        else:        return distances

    # end of function orfenddistance


    def get_orfs_of_graph(self,organism=None):
        """
        Get all the Orf objects of this graph

        @type  organism: *
        @param organism: organism identifier (or None)

        @rtype:  dictionary (or list if organism is specified)
        @return: dictionary with organisms (keys) and list of Orf objects (values),
                 or only list of Orf objects if an organism identifier was specified
        """
        orfs = {}
        for (key,(org1,orf1),(org2,orf2)),pacbporf in self.pacbps.iteritems():
            if orfs.has_key(org1):
                if orf1 in [ o.id for o in orfs[org1] ]:
                    pass
                else:
                    orfs[org1].append( pacbporf.orfQ )
            else:
                orfs[org1] = [ pacbporf.orfQ ]
            if orfs.has_key(org2):
                if orf2 in [ o.id for o in orfs[org2] ]:
                    pass
                else:
                    orfs[org2].append( pacbporf.orfS )
            else:
                orfs[org2] = [ pacbporf.orfS ]
        # now check what to return
        if organism and organism not in self.organism_set():
            raise OrganismNotPresentInGraph
        elif organism:
            return orfs[organism]
        else:
            return orfs

    # end of function get_orfs_of_graph

    ############################################################################
    #### Functions concerning the PacbP(ORF)s of this CodingBlockGraph      ####
    ############################################################################

    def has_pacbps(self):
        """
        Are there any PacbP(ORF) objects in this CodingBlockGraph?

        @rtype:  Boolean
        @return: Are there any pacbp objects in this CodingBlockGraph? 
        """
        if self.pacbps: return True
        else:           return False

    # end of function has_pacbps


    def has_all_pacbps(self):
        """
        Are all edges upproted by a PacbP(ORF) object in this CodingBlockGraph?

        @rtype:  Boolean
        @return: Are all edges supported by at least 1 PacbP(ORF) objects 
        """
        if not self.pacbps:
            return False 
        elif self.edge_count() > len(self.pacbps):
            return False
        else:
            all_pacbps_present = True 
            for nodeQ,nodeS in self.pairwisecrosscombinations_node():
                if self.has_edge(nodeQ,nodeS):
                    if not self.get_pacbps_by_nodes(node1=nodeQ,node2=nodeS):
                        all_pacbps_present = False
                        break
            return all_pacbps_present

    # end of function has_all_pacbps


    def remove_alternatives_from_pacbps_dict(self):
        """
        Remove PacbP objects that overlap (to much) with each other

        @attention: only applicable when self.edge_count() < len(self.pacbps)
        """
        # only perform this function if graph edge_count
        # not equals self.pacbps size
        if self.edge_count() < len(self.pacbps):
            accepted, rejected = remove_alternatives_from_pacbps_dict(
                self.pacbps,
                overlap_ratio=self.ALTERNATIVE_ALIGNMENT_OVERLAP_RATIO
                )
            if len(accepted) < len(self.pacbps):
                self.pacbps = accepted
        else:
            pass

    # end of function remove_alternatives_from_pacbps_dict


    def keep_only_best_alternative_from_pacbps_dict(self):
        """
        Remove all except the highest scoring PacbP(ORF) when alternatives are present

        @attention: only applicable when self.edge_count() < len(self.pacbps)
        """
        if self.edge_count() < len(self.pacbps):
            accepted, rejected = remove_alternatives_from_pacbps_dict(
                    self.pacbps,
                    keep_only_highest_scoring=True)
            # replace self.pacbps with novel, shortened dict of pacbps
            self.pacbps = accepted 
        else:
            pass

    # end of function keep_only_best_alternative_from_pacbps_dict


    def make_pacbps_for_missing_edges(self,use_maxsr=False,aa_extra_offset=1):
        """
        Try to make PabcpP(ORF) objects with ClustalW when CBG is non-complete

        @type  use_maxsr: Boolean
        @param use_maxsr: when True, use MAXSR in stead of OMSR for novel PacbPs

        @type  aa_extra_offset: integer
        @param aa_extra_offset: extra AA offset around OMSR or MAXSR

        @attention: only applicable when the CBG has an OMSR
        @attention: only applicable when self.connectivitysaturation() < 1.0
        """
        if self.connectivitysaturation() < 1.0:
            make_pacbps_for_missing_edges(self,
                use_maxsr=use_maxsr,
                aa_extra_offset=aa_extra_offset
                )

    # end of function make_pacbps_for_missing_edges


    def harvest_pacbps_from_pacbpcollection(self,pcg):
        """
        Get PacbP(ORF)s from a (parental) PCG or another CBG

        @type  pcg: PacbpCollectionGraph
        @param pcg: PacbpCollectionGraph (or CodingBlockGraph)
        """
        for (node1,node2) in self.pairwisecrosscombinations_node():
            # check if these nodes are indeed present as an edge
            if not self.has_edge(node1,node2): continue
            # get pacbps from pcg
            pacbps = pcg.get_pacbps_by_nodes(node1=node1,node2=node2)
            for pacbp in pacbps:
                key = pacbp.construct_unique_key(node1,node2)
                self.pacbps[(key,node1,node2)] = pacbp

    # end of function harvest_pacbps_from_pacbpcollection


    def split_codingblock_on_alternatives_in_pacbps_dict(self,**kwargs):
        """
        Split CodingBlockGraph on difference between max and min spanning range

        @attention: See codingblock_splitting.
                    split_codingblock_on_alternatives_in_pacbps_dict()
        """
        return codingblock_splitting.split_codingblock_on_alternatives_in_pacbps_dict(self,**kwargs)

    # end of function split_codingblock_on_alternatives_in_pacbps_dict

    ############################################################################
    #### Functions concerning (potential) TSS in this CodingBlockGraph      ####
    ############################################################################

    def have_all_starts(self):
        """
        Do all Orfs in the PacbPs have Methionines?

        @attention: only functional if PacbPORFs are present in the CBG

        @rtype:  Boolean
        @return: True or False
        """
        done = []
        for ((a,b,c,d),n1,n2), pacbporf in self.pacbps.iteritems():
            if n1 not in done:
                if pacbporf.orfQ.has_start() == False:
                    return False
                else:
                    done.append(n1)
            if n2 not in done:
                if pacbporf.orfS.has_start() == False:
                    return False
                else:
                    done.append(n2)
        else:
            return True

    # end of function have_all_starts


    def have_all_starts_upstream_of_omsr(self,omsr_offset=10):
        """
        Are there (putative) TSS-Methionines upstream of the OMSR in all Orfs?

        @type  omsr_offset: integer
        @param omsr_offset: AA offset into 3' side of OMSR to take into account

        @attention: only functional if PacbPORFs are present in the CBG

        @rtype:  Boolean
        @return: True or False
        """
        omsr = self.overall_minimal_spanning_range()
        for org,orflist in self.get_orfs_of_graph().iteritems():
            node = ( org,orflist[0].id )
            offset = min([ min(omsr[node]) + omsr_offset, max(omsr[node]) ])
            if not orflist[0].has_start_upstream_of(offset):
                return False
        else:
            return True

    # end of function have_all_starts_upstream_of_omsr


    def minimal_eligable_tss_position(self,organism):
        """
        Get the position on the Orf of this organism in the CBG from
        where to take TSS positions into account

        @attention: See codingblock_collectionharvesting.
                    minimal_eligable_tss_position for documentation
        """
        return cbg2colgra.minimal_eligable_tss_position(self,organism)

    # end of function minimal_eligable_tss_position


    def maximal_eligable_tss_position(self,organism):
        """
        Get the position on the Orf of this organism in the CBG untill
        where to take TSS positions into account

        @attention: See codingblock_collectionharvesting.
                    maximal_eligable_tss_position for documentation
        """
        return cbg2colgra.maximal_eligable_tss_position(self,organism)
        
    # end of function maximal_eligable_tss_position

    ############################################################################
    #### Helper Functions for MSR, OMSR, MAXSR                              ####
    ############################################################################

    def _get_ranges_by_nodes(self,node1,node2,
        _limit_range_node1=[],
        _limit_range_node2=[]):
        """
        Function in use in minimal_spanning_range() and maximal_spanning_range()

        @type  node1: *
        @param node1: One Node

        @type  node2: *
        @param node2: One Node

        @type  _limit_range_node1: list
        @param _limit_range_node1: integer list, allowed interval for node1

        @type  _limit_range_node2: list
        @param _limit_range_node2: integer list, allowed interval for node2

		@rtype:  ( Set, Set )
		@return: tuple of 2 sets, range for Query and Sbjct

        @attention: node1 MUST be the query node, node2 MUST be the sbjct node!
        @attention: function requires PacbPORF or PacbPDNA objects, not PacbPs

        """
        # set known ranges by pre-limit known ranges
        thisrangeQ = Set()
        thisrangeS = Set()

        # IMPORTANT REMINDER!!! Most of the classes that use this function,
        # are (semi-)complete K(s) Graphs, with a single PacbPobj per edge.
        # In that case the iteration over ``get_pacbps_by_nodes`` will obviously
        # result in a single pacbpobj to be returned.
        # However, this function is called as well from more ``collection``
        # types of graphs: (semi-)complete Graphs with >1 PacbPobj per edge.
        # In an older (<01/01/2009) version of this function, this caused
        # a - for long - undetected flaw in the function's output.
        # All PacbPobject of a single edge must here be treated as a `single`
        # PacbPobject. So, its range is the SUM of all the individual ranges,
        # and not the cross-section!

        for pacbpobj in self.get_pacbps_by_nodes(node1=node1,node2=node2):
            if pacbpobj.__class__.__name__ in ['PacbPORF','PacbPDNA']:
                # PacbPORF or further inhertited object.
                # The alignment might be extended,so
                # so take original alignment coordinates
                spos = pacbpobj._positions[pacbpobj._original_alignment_pos_start]
                epos = pacbpobj._positions[pacbpobj._original_alignment_pos_end-1]
                startQ, endQ = spos.query_pos, epos.query_pos + 1
                startS, endS = spos.sbjct_pos, epos.sbjct_pos + 1

            else:
                # a Pacbp or PacbpDNA object, take start & end
                startQ, endQ = pacbpobj.query_start, pacbpobj.query_end
                startS, endS = pacbpobj.sbjct_start, pacbpobj.sbjct_end

            # if _limit_range_node1 and DNA-type PacbP,
            # then limit by allowed interval _limit_range_node1 
            if pacbpobj.__class__.__name__ in ['PacbPORF','PacbPDNA']:
                if _limit_range_node1:
                    if startQ < min(_limit_range_node1):
                        algQpos = pacbpobj.alignmentposition_by_query_pos(
                                    min(_limit_range_node1) )
                        spos = pacbpobj._positions[ algQpos ]
                        startQ, startS = spos.query_pos, spos.sbjct_pos
                    if endQ-1 > max(_limit_range_node1):
                        algQpos = pacbpobj.alignmentposition_by_query_pos(
                                    max(_limit_range_node1) )
                        epos = pacbpobj._positions[  algQpos ]
                        endQ, endS = epos.query_pos+1, epos.sbjct_pos+1

                if _limit_range_node2:
                    if startS < min(_limit_range_node2):
                        algSpos = pacbpobj.alignmentposition_by_sbjct_pos(
                                    min(_limit_range_node2) )
                        spos = pacbpobj._positions[ algSpos ]
                        startQ, startS = spos.query_pos, spos.sbjct_pos
                    if endS-1 > max(_limit_range_node2):
                        algSpos = pacbpobj.alignmentposition_by_sbjct_pos(
                                    max(_limit_range_node2) ) 
                        epos = pacbpobj._positions[ algSpos ]
                        endQ, endS = epos.query_pos+1, epos.sbjct_pos+1

            else:
                # A Pacbp or PacbpDNA object. Limiting ranges by preset ranges
                # (_limit_range_node1 and _limit_range_node2) is less easy;
                # the sequence has to be checked for occurring gaps etc.
                # For now, leave it here as such. This means that the function
                # overall_minimal_spanning_range() has a different result for
                # PacbP objects than for PacbPORF, PacbPDNA objects
                pass

            # and update the ranges; this can only enlarge the final range
            thisrangeQ.update( range(startQ, endQ) )
            thisrangeS.update( range(startS, endS) )

        # and return the ranges
        return thisrangeQ, thisrangeS

    # end of function _get_ranges_by_nodes

    ############################################################################
    #### Functions concerning the MinimalSpanningRange (MSR)                ####
    ############################################################################

    def minimal_spanning_range(self,organism=None,node=None):
        """
        Get the Minimal Spanning Range (MSR) of this CBG

        @attention: only applicable when PacbP(ORF)s are present in the graph
                    (perform the function harvest_pacbps_from_crossdata())

        @type  organism: *
        @param organism: organism identifier (or None)

        @type  node: *
        @param node: Node identifier (or None)

        @rtype:  dictionary (or Set if organism is specified)
        @return: dictionary with nodes (keys) and spanning range Sets
                 (values), or only the spanning range Set if an organism
                 or node identifier was specified
        """
        if organism and organism not in self.organism_set():
            raise OrganismNotPresentInGraph
        if node and node not in self.get_nodes():
            raise NodeNotPresentInGraph
        # try to get from the objects' cached attributes
        if self._msr:
            if organism:    return self._msr[self.node_by_organism(organism)]
            elif node:      return self._msr[node]
            else:           return self._msr

        initial_ranges = {}
        # first iteration: get minimal spanning range by mutual comparison
        for (node1,node2) in self.pairwisecrosscombinations_node():
            # check if these nodes are indeed present as an edge
            if not self.has_edge(node1,node2): continue
            # get ranges for these nodes
            thisrangeQ, thisrangeS = self._get_ranges_by_nodes(node1,node2)
            # set to ranges
            if not initial_ranges.has_key(node1):
                initial_ranges[node1] = thisrangeQ
            else:
                initial_ranges[node1].intersection_update(thisrangeQ)
            if not initial_ranges.has_key(node2):
                initial_ranges[node2] = thisrangeS
            else:
                initial_ranges[node2].intersection_update(thisrangeS)

        # return ranges or the range of a specific organism or node
        if organism:
            return initial_ranges[self.node_by_organism(organism)]
        elif node:
            return initial_ranges[node]
        else:
            return initial_ranges

    # end of function minimal_spanning_range


    def minimal_spanning_range_sizes(self):
        """
        Get the sizes of the Minimal Spanning Range (MSR) of this CBG

        @attention: only meaningfull for (nearly) K(s) CBGs
        @attention: only applicable when PacbP(ORF)s are present in the graph
                    (perform the function harvest_pacbps_from_crossdata())

        @rtype:  dictionary
        @return: dict with nodes (keys) and size of spanning range Sets (values)

        """
        d = {}
        for node, r in self.minimal_spanning_range().iteritems():
            d[node] = len(r)
        return d

    # end of function minimal_spanning_range_sizes


    def has_minimal_spanning_range(self):
        """
        Has this (nearly) K(s) CBG a minimal spanning range?

        @attention: only functional if pacbp's are present in the graph
                    (perform the function `harvest_pacbps_from_crossdata`)

        @attention: only meaningfull for (nearly) K(s) CBGs 

        @rtype:  Boolean
        @return: True or False
        """
        if self.minimal_spanning_range_sizes().values():
            if min(self.minimal_spanning_range_sizes().values()) == 0:
                return False
            else:
                return True
        else:
            return False

    # end of function minimal_spanning_range_sizes


    def has_msr(self):
        """
        @attention: alias function name of has_minimal_spanning_range
        """
        return self.has_minimal_spanning_range()
        
    # end of function has_msr

    
    ############################################################################
    #### Functions concerning the OverallMinimalSpanningRange (OMSR)        ####
    ############################################################################

    def overall_minimal_spanning_range(self,organism=None,node=None):
        """
        Get the Overall Minimal Spanning Range (OMSR) of this CBG

        @attention: only applicable when PacbP(ORF)s are present in the graph
                    (perform the function harvest_pacbps_from_crossdata())

        @type  organism: *
        @param organism: Organism identifier (or None)

        @type  node: * 
        @param node: Node identifier (or None)

        @rtype:  dictionary (or Set if organism is specified)
        @return: dictionary with nodes (keys) and spanning range Sets
                 (values), or only the spanning range Set if an organism
                 or node identifier was specified
        """
        if organism and organism not in self.organism_set():
            raise OrganismNotPresentInGraph
        if node and node not in self.get_nodes():
            raise NodeNotPresentInGraph
        # try to get from the objects' cached attributes
        if self._omsr:
            if organism:    return self._omsr[self.node_by_organism(organism)]
            elif node:      return self._omsr[node]
            else:           return self._omsr

        # FIRST iteration; get normal minimal spanning range
        initial_ranges = self.minimal_spanning_range()
        final_ranges = deepcopy(initial_ranges)
        # SECOND iteration: get minimal spanning range by overall comparison
        for (node1,node2) in self.pairwisecrosscombinations_node():
            # check if these nodes are indeed present as an edge
            if not self.has_edge(node1,node2): continue
            # get ranges for these nodes
            thisrangeQ, thisrangeS = self._get_ranges_by_nodes(
                    node1,node2,
                    _limit_range_node1=initial_ranges[node1],
                    _limit_range_node2=initial_ranges[node2],
                    )
            # set to ranges
            final_ranges[node1].intersection_update(thisrangeQ)
            final_ranges[node2].intersection_update(thisrangeS)

        # return ranges or the range of a specific organism or node
        if organism:
            return final_ranges[self.node_by_organism(organism)]
        elif node:
            return final_ranges[node]
        else:
            return final_ranges

    # end of function overall_minimal_spanning_range


    def overall_minimal_spanning_range_sizes(self):
        """
        Get the sizes of the Overall Minimal Spanning Range (OMSR) of this CBG

        @attention: only meaningfull for (nearly) K(s) CBGs
        @attention: only applicable when PacbP(ORF)s are present in the graph
                    (perform the function harvest_pacbps_from_crossdata())

        @rtype:  dictionary
        @return: dict with nodes (keys) and size of spanning range Sets (values)
        """
        d = {}
        for node, r in self.overall_minimal_spanning_range().iteritems():
            d[node] = len(r)
        return d

    # end of function overall_minimal_spanning_range_sizes


    def has_overall_minimal_spanning_range(self):
        """
        Has this CBG an Overall Minimal Spanning Range (OMSR)?

        @attention: only meaningfull for (nearly) K(s) CBGs
        @attention: only applicable when PacbP(ORF)s are present in the graph
                    (perform the function harvest_pacbps_from_crossdata())

        @rtype:  Boolean
        @return: True or False
        """
        try:
            if self.overall_minimal_spanning_range_sizes().values():
                # 8/4/09 smallest omsr size must be 3
                if min(self.overall_minimal_spanning_range_sizes().values()) <\
                self.MINIMAL_OVERAL_SPANNING_RANGE_SIZE:
                    return False
                else:
                    return True
            else:
                return False
        except:
            return False

    # end of function has_overall_minimal_spanning_range

    
    def has_omsr(self):
        """
        @attention: alias function name of has_overall_minimal_spanning_range
        """
        return self.has_overall_minimal_spanning_range()

    # end of function has_omsr


    def omsr_starts(self):
        """
        Return a dictionary of OMSR start coodinates (keys==nodes)

        @rtype:  dictionary 
        @return: dictionary with nodes (*, presumably tuple) as keys,
                 OMSR start coords (integer) as values
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
        @return: dictionary with nodes (*, presumably tuple) as keys,
                 OMSR end coords (integer) as values
        """
        retdict = {}
        for node,omsr in self.overall_minimal_spanning_range().iteritems():
            retdict[node] = max(omsr)
        return retdict

    # end of function omsr_ends


    def omsr_length(self):
        """
        Get the (smallest) OMSR length for this CBG

        @attention: no OMSR will return 0

        @rtype:  integer
        @return: smallest OMSR length of this CBG in AA coords
        """
        if self.has_overall_minimal_spanning_range():
            return min(self.overall_minimal_spanning_range_sizes().values())
        else:
            return 0

    # end of function omsr_length


    def omsrlength(self):
        """
        @attention: alias for omsr_length()
        """
        return self.omsr_length()

    # end of function omsrlength


    def omsr_identityscore(self):
        """
        Get the CBGs identityscore in the OMSR region
   
        @attention: a different measure as total_weight()

        @rtype:  float
        @return: indentity score in the OMSR region ( >>0.0 .. <=1.0 )
        """
        omsr = self.overall_minimal_spanning_range()
        wts = []
        for (key,nodeQ,nodeS), pacbporf in self.pacbps.iteritems():
            sta,end = min(omsr[nodeQ]), max(omsr[nodeQ])
            wt = pacbporf.identityscore_slice_by_abs_protein_query(sta,end)
            wts.append(wt)
        return float(sum(wts))/self.edge_count()
    
    # end of function omsr_identityscore

    ############################################################################
    #### Functions concerning the MaximalSpanningRange (MAXSR)              ####
    ############################################################################

    def maximal_spanning_range(self,organism=None,node=None):
        """
        Get the Maximal Spanning Range (MAXSR) of this CBG

        @attention: only meaningfull for (nearly) K(s) CBGs
        @attention: only applicable when PacbP(ORF)s are present in the graph
                    (perform the function harvest_pacbps_from_crossdata())

        @type  organism: *
        @param organism: organism identifier (or None)

        @type  node: *
        @param node: Node identifier (or None)

        @rtype:  dictionary (or Set if organism is specified)
        @return: dictionary with nodes (keys) and spanning range Sets
                 (values), or only the spanning range Set if an organism
                 or node identifier was specified
        """
        if organism and organism not in self.organism_set():
            raise OrganismNotPresentInGraph
        if node and node not in self.get_nodes():
            raise NodeNotPresentInGraph
        # try to get from the objects' cached attributes
        if self._maxsr:
            if organism:    return self._maxsr[self.node_by_organism(organism)]
            elif node:      return self._maxsr[node]
            else:           return self._maxsr

        # create Maximal Spanning Range dictionary
        ranges = {}
        for (node1,node2) in self.pairwisecrosscombinations_node():
            # get ranges for these nodes
            thisrangeQ, thisrangeS = self._get_ranges_by_nodes(node1,node2)
            # set to ranges
            if not ranges.has_key(node1):
                ranges[node1] = thisrangeQ
            else:
                ranges[node1].union_update(thisrangeQ)
            if not ranges.has_key(node2):
                ranges[node2] = thisrangeS
            else:
                ranges[node2].union_update(thisrangeS)

        # return all ranges or the range of a specific organism
        if organism:
            for node, orgrange in ranges.iteritems():
                if self._organism_from_node(node) == organism:
                    return orgrange
        else:
            return ranges

    # end of function maximal_spanning_range


    def maximal_spanning_range_sizes(self):
        """
        Get the sizes of the Maximal Spanning Range (MAXSR) of this CBG

        @attention: only meaningfull for (nearly) K(s) CBGs
        @attention: only applicable when PacbP(ORF)s are present in the graph
                    (perform the function harvest_pacbps_from_crossdata())

        @rtype:  dictionary
        @return: dict with nodes (keys) and size of spanning range Sets (values)
        """
        d = {}
        for node, r in self.maximal_spanning_range().iteritems():
            d[node] = len(r)
        return d

    # end of function maximal_spanning_range_sizes


    def superinposed_left_maximal_spanning_range(self,organism=None):
        """
        Get the LEFT side of the Maximal Spanning Range (MAXSR) of this CBG

        @attention: only meaningfull for (nearly) K(s) CBGs
        @attention: only applicable when PacbP(ORF)s are present in the graph
                    (perform the function harvest_pacbps_from_crossdata())

        @type  organism: *
        @param organism: organism identifier (or None)

        @rtype:  dictionary (or Set if organism is specified)
        @return: dictionary with nodes (keys) and spanning range Sets (values),
                 or only the Set if an organism identifier was specified
        """
        omsr  = self.overall_minimal_spanning_range()
        maxsr = self.maximal_spanning_range()
        smaxsrL = {}
        for node in omsr.keys():
            smaxsrL[node] = min(omsr[node]) - min(maxsr[node])
        max_smaxsrL = max(smaxsrL.values())
        for node in omsr.keys():
            if smaxsrL[node] == max_smaxsrL:
                smaxsrL[node] = Set()
            else:
                coordSta = min(omsr[node])-(max_smaxsrL-smaxsrL[node])
                coordEnd = min(omsr[node])
                smaxsrL[node] = Set(range( coordSta, coordEnd ) )
        if organism:
            return smaxsrL[ self.node_by_organism(organism) ]
        else:
            return smaxsrL

    # end of function superinposed_left_maximal_spanning_range


    def superinposed_rigth_maximal_spanning_range(self,organism=None):
        """
        Get the RIGTH side of the Maximal Spanning Range (MAXSR) of this CBG

        @attention: only meaningfull for (nearly) K(s) CBGs
        @attention: only applicable when PacbP(ORF)s are present in the graph
                    (perform the function harvest_pacbps_from_crossdata())

        @type  organism: *
        @param organism: organism identifier (or None)

        @rtype:  dictionary (or Set if organism is specified)
        @return: dictionary with nodes (keys) and spanning range Sets (values),
                 or only the Set if an organism identifier was specified
        """
        omsr  = self.overall_minimal_spanning_range()
        maxsr = self.maximal_spanning_range()
        smaxsrR = {}
        for node in omsr.keys():
            smaxsrR[node] = max(maxsr[node]) - max(omsr[node])
        max_smaxsrR = max(smaxsrR.values())

        for node in omsr.keys():
            if smaxsrR[node] == max_smaxsrR:
                smaxsrR[node] = Set()
            else:
                coordSta = max(omsr[node])
                coordEnd = max(omsr[node])+(max_smaxsrR-smaxsrR[node])
                smaxsrR[node] = Set( range( coordSta, coordEnd ) )

        if organism:
            return smaxsrR[ self.node_by_organism(organism) ]
        else:
            return smaxsrR

    # end of function superinposed_rigth_maximal_spanning_range


    ####################################################################
    #### Functions for codingperiodicity data                       ####
    ####################################################################

    def codingperiodicitystatistics(self):
        """
        @attention: see codingblock_codingperiodicity.
                    codingperiodicitystatistics for documentation
        """
        return codingblock_codingperiodicity.codingperiodicitystatistics(self)

    # end of function codingperiodicitystatistics


    def printcodingperiodicitystatistics(self):
        """
        @attention: see codingblock_codingperiodicity.
                    printcodingperiodicitystatistics for documentation
        """
        codingblock_codingperiodicity.printcodingperiodicitystatistics(self)

    # end of function printcodingperiodicitystatistics


    ####################################################################
    #### Functions for CodingBlockGraph improvement                 ####
    ####################################################################

    def improvealignment(self,**kwargs):
        """
        @attention: see codingblock_optimization.improvealignment()
        """
        return codingblock_optimization.improvealignment(self,**kwargs)

    # end of function improvealignment


    def correct_pacbpgaps_nearby_omsr(self,**kwargs):
        """
        @attention: see codingblock_optimization.correct_pacbpgaps_nearby_omsr()
                    for documentation
        """
        return codingblock_optimization.correct_pacbpgaps_nearby_omsr(
                self,**kwargs)

    # end of function correct_pacbpgaps_nearby_omsr

    ############################################################################
    #### Functions for TCODE data accessibility in CodingBlockGraphs        ####
    ############################################################################

    def omsr_tcode_score(self,organism=None):
        """ 
        Get the average EMBOSS-TCODE score of the OMSR of this CBG

        @attention: alias for msr_tcode_score
        """
        return self.msr_tcode_score(organism=organism)

    # end of function omsr_tcode_score


    def msr_tcode_score(self,organism=None):
        """
        Get the average EMBOSS-TCODE score of the OMSR of this CBG

        @type  organism: *
        @param organism: Organism identifier (or None)

        @rtype:  float
        @return: average TCODE score, either of the OMSR of all organisms,
                 or of a single Organism identifier when provided as argument
        """
        tcode_scores = {}
        for node, r in self.overall_minimal_spanning_range().iteritems():
            org = self._organism_from_node(node)
            if organism and organism != org: continue
            orf = self.get_orfs_of_graph(organism=org)[0]
            score = orf.tcode_score(start=min(r),end=max(r)+1)
            tcode_scores[org] = score
        if organism and organism not in self.organism_set():
            raise OrganismNotPresentInGraph
        elif organism:
            return tcode_scores[organism]
        else:
            return sum(tcode_scores.values()) / float(len(tcode_scores))

    # end of function msr_tcode_score


    def average_tcode_score_of_position(self,organism=None,aaPos=None):
        """
        TODO: this functionality takes waaaaaaaaaaay to much time (to much iterations)
        """
        # if no organism specified, take strongest connected node
        if not organism:
            organism = self._organism_from_node( self.strongest_connected_node() )
        # get pacbporfs of this organism
        pacbporfs = self.get_pacbps_by_organism(organism)

        scores = []
        is_first = True
        for pacbporf in pacbporfs:
            algPos = pacbporf.alignmentposition_by_query_pos(qposPY=aaPos)
            algPosObj = pacbporf._positions[algPos]
            if is_first:
                # get the EMBOSS-TCODE data about this QUERY position
                if not (algPosObj.isa_gap and algPosObj.query_dna_start==0):
                    scores.append( pacbporfs[0].orfQ.tcode_score(start=aaPos,end=aaPos+1) )
                is_first=False
            # get the EMBOSS-TCODE data about this QUERY position
            if not (algPosObj.isa_gap and algPosObj.sbjct_dna_start==0):
                aaPosSbjct = algPosObj.sbjct_pos
                scores.append( pacbporfs[0].orfQ.tcode_score(start=aaPosSbjct,end=aaPosSbjct+1) )
        if not scores:
            return (0.0, 0)  # apperently something went wrong; for compatibility, return zero        
        else:
            return ( sum(scores)/float(len(scores)), len(scores) )

    # end of function average_tcode_score_of_position


    ############################################################################
    #### Functions for (PSSM) gene structure element collection harvesting  ####
    ############################################################################

    def align_stop_codons(self):
        """
        Create the AlignedStopCodonGraph of this CodingBlockGraph
        """
        self._stopcodongraph = AlignedStopCodonGraph(self)

    # end of function align_stop_codons


    def harvest_elegiable_tss_sites(self,**kwargs):
        """
        Scan the (elegiable) TSSs in this CBG

        @attention: see codingblock_collectionharvesting.
                    harvest_elegiable_tss_sites for documentation
        """
        # perform cbg2colgra.harvest_elegiable_tss_sites function
        tsscolgra = cbg2colgra.harvest_elegiable_tss_sites(self,**kwargs)
        # done; set graph object to object
        self._startcodongraph = tsscolgra

    # end of function harvest_elegiable_tss_sites


    #def harvest_elegiable_donor_sites(self,next=None,
    #    projected_donors={},
    #    forced_codingblock_ends={},
    #    allow_phase_shift=False,
    #    allow_non_canonical=False,
    #    store_all_projected_sites=False,
    #    min_donor_pssm_score = DONORPSSMOPTIONS.MIN_DONOR_PSSM_SCORE,
    #    non_canonical_min_donor_pssm_score =\
    #            DONORPSSMOPTIONS.NON_CANONICAL_MIN_DONOR_PSSM_SCORE,
    #    enlarge_5p_boundary_by=None,
    #    enlarge_3p_boundary_by=None,
    #    aligned_donor_max_triplet_distance=ALIGNED_DONOR_MAX_TRIPLET_DISTANCE):
    #    """
    #    Make the DonorSiteCollectionGraph for this CodingBlockGraph
    #    
    #    @attention: see codingblock_collectionharvesting.
    #                harvest_elegiable_donor_sites() for documentation
    #    """
    #    colgra = cbg2colgra.harvest_elegiable_donor_sites(self,next=next,
    #            MIN_DONOR_PSSM_SCORE=min_donor_pssm_score,
    #            NON_CANONICAL_MIN_DONOR_PSSM_SCORE=\
    #                    non_canonical_min_donor_pssm_score,
    #            ALLOW_NON_CANONICAL_DONOR=allow_non_canonical,
    #            ALIGNED_DONOR_MAX_TRIPLET_DISTANCE=\
    #                    aligned_donor_max_triplet_distance,
    #            enlarge_5p_boundary_by=enlarge_5p_boundary_by, # in AA coords
    #            enlarge_3p_boundary_by=enlarge_3p_boundary_by, # in AA coords
    #            allow_phase_shift=allow_phase_shift,
    #            projected_donors=projected_donors,
    #            store_all_projected_sites=store_all_projected_sites,
    #            forced_codingblock_ends=forced_codingblock_ends
    #            )
    #    # done; set graph object to object
    #    self._splicedonorgraph = colgra
    #
    ## end of function harvest_elegiable_donor_sites


    #def harvest_elegiable_acceptor_sites(self,prev=None,
    #    projected_acceptors={},
    #    forced_codingblock_ends={},
    #    allow_phase_shift=False,
    #    store_all_projected_sites=False,
    #    min_acceptor_pssm_score = ACCEPTORPSSMOPTIONS.MIN_ACCEPTOR_PSSM_SCORE,
    #    allow_non_canonical_acceptor =\
    #            ACCEPTORPSSMOPTIONS.ALLOW_NON_CANONICAL_ACCEPTOR,
    #    non_canonical_min_acceptor_pssm_score =\
    #            ACCEPTORPSSMOPTIONS.NON_CANONICAL_MIN_ACCEPTOR_PSSM_SCORE,
    #    enlarge_5p_boundary_by=None,
    #    enlarge_3p_boundary_by=None,
    #    aligned_acceptor_max_triplet_distance=\
    #            ALIGNED_ACCEPTOR_MAX_TRIPLET_DISTANCE):
    #    """
    #    Make the AcceptorSiteCollectionGraph for this CodingBlockGraph
    #
    #    @attention: see codingblock_collectionharvesting.
    #                harvest_elegiable_acceptor_sites() for documentation
    #    """
    #    colgra = cbg2colgra.harvest_elegiable_acceptor_sites(self,prev=prev,
    #            MIN_ACCEPTOR_PSSM_SCORE = min_acceptor_pssm_score,
    #            ALLOW_NON_CANONICAL_ACCEPTOR = allow_non_canonical_acceptor,
    #            NON_CANONICAL_MIN_ACCEPTOR_PSSM_SCORE =\
    #                    non_canonical_min_acceptor_pssm_score,
    #            ALIGNED_ACCEPTOR_MAX_TRIPLET_DISTANCE=\
    #                    aligned_acceptor_max_triplet_distance,
    #            enlarge_5p_boundary_by=enlarge_5p_boundary_by, # in AA coords
    #            enlarge_3p_boundary_by=enlarge_3p_boundary_by, # in AA coords
    #            allow_phase_shift=allow_phase_shift,
    #            projected_acceptors=projected_acceptors,
    #            store_all_projected_sites=store_all_projected_sites,
    #            forced_codingblock_ends=forced_codingblock_ends,
    #            )
    #    # done; set graph object to object
    #    self._spliceacceptorgraph = colgra
    #
    ## end of function harvest_elegiable_acceptor_sites

    ############################################################################
    #### Functions for MultipleAligmentData PAOC / PASC                     ####
    ############################################################################

    def set_multiple_alignment_data_structures(self):
        """
        Set all the VISTA like tracks
        """
        for organism in self.organism_set():
            # set positional_alignment_occurrence_count (self._paoc_per_organism)
            self.set_positional_alignment_occurrence_count(organism=organism)
            # set positional_alignment_similarity_count (self._pasc_per_organism)
            self.set_positional_alignment_similarity_count(organism=organism)

    # end of function set_multiple_alignment_data_structures


    def get_positional_alignment_occurrence_count(self,organism=None):
        """
        @rtype:  {}
        @return: data structure describing multiple aligned occurrences
        """
        if self._paoc_per_organism.has_key(organism):
            return self._paoc_per_organism[organism]
        else:
            return self.set_positional_alignment_occurrence_count(organism=organism)

    # end of function get_positional_alignment_occurrence_count


    def set_positional_alignment_occurrence_count(self,organism=None):
        """
        @rtype:  {}
        @return: data structure describing multiple aligned occurrences
        """
        if organism not in self.organism_set():
            raise OrganismNotPresentInGraph
        area = {}
        this2others = {}
        for aapos in self.maximal_spanning_range(organism=organism):
            area[aapos] = 0
        for ((a,b,c,d),(org1,orf1),(org2,orf2)), pacbporf in self.pacbps.iteritems():
            if org1 == organism:
                thispacbporf = pacbporf
                othernode = (org2,orf2)
            elif org2 == organism:
                thispacbporf = pacb.swap_query_and_sbjct(pacbporf)
                othernode = (org1,orf1)
            else:
                continue
            # loop over the original alignment positions
            for posid in range(thispacbporf._original_alignment_pos_start,thispacbporf._original_alignment_pos_end):
                posobj = thispacbporf._positions[posid]
                area[posobj.query_pos] += 1
                if this2others.has_key(othernode):
                    this2others[othernode][posobj.sbjct_pos] = posobj.query_pos
                else:
                    this2others[othernode] = { posobj.sbjct_pos: posobj.query_pos }

        for (orgX,orfX), otherrange in this2others.iteritems():
            for ((a,b,c,d),(org1,orf1),(org2,orf2)), pacbporf in self.pacbps.iteritems():
                if (org1,orf1) == (orgX,orfX) and org2 != organism:
                    thispacbporf = pacbporf
                else:
                    continue
                # loop over the original alignment positions
                for posid in range(thispacbporf._original_alignment_pos_start,thispacbporf._original_alignment_pos_end):
                    posobj = thispacbporf._positions[posid]
                    if posobj.query_pos in otherrange.keys():
                        orig_query_pos = otherrange[posobj.query_pos]
                        area[orig_query_pos] += 1


        # devide by the maximum number of counts possible
        combis = self.crosscombinationcount()
        for k in area.keys(): area[k] = min([ float(area[k])/float(combis), 1.0 ])
        # set this data structure to the object
        self._paoc_per_organism[organism] = area
        # return the scores for this CodingBlock
        return area

    # end of function set_positional_alignment_occurrence_count


    def get_gff_of_positional_alignment_occurrence_count(self,organism=None,fref='FREF'):
        """
        @rtype:  []
        @return: list with gff tuples
        """
        # get dict with positional_alignment_occurrence_count data
        occposcounts = self.get_positional_alignment_occurrence_count(organism=organism)
        # get the orf belonging to this organism of this CodingBlockGraph            
        theorf = self.get_orfs_of_graph(organism=organism)[0]
        gffs = []
        for occpos,score in occposcounts.iteritems():
            dnapos = theorf.aapos2dnapos(occpos)
            gff = ( fref, 'abfgp', 'AbfgpMultipleOccurrence', dnapos+1, dnapos+3, score, '+', '.', "AMO %s" % fref )
            gffs.append(gff)
        gffs.sort()
        return gffs

    # end of function get_gff_of_positional_alignment_occurrence_count


    def get_positional_alignment_similarity_count(self,organism=None):
        """
        @rtype:  {}
        @return: data structure describing multiple aligned similarities
        """
        if self._pasc_per_organism.has_key(organism):
            return self._pasc_per_organism[organism]
        else:
            return self.set_positional_alignment_similarity_count(organism=organism)

    # end of function get_positional_alignment_similarity_count


    def set_positional_alignment_similarity_count(self,organism=None):
        """
        @rtype:  {}
        @return: data structure describing multiple aligned similarities
        """
        identity_score   = 1.0
        similarity_score = 0.5
        gap_score        = -1.0
        if organism not in self.organism_set():
            raise OrganismNotPresentInGraph
        area = {}
        this2others = {}
        for aapos in self.maximal_spanning_range(organism=organism):
            area[aapos] = 0
        for ((a,b,c,d),(org1,orf1),(org2,orf2)), pacbporf in self.pacbps.iteritems():
            if org1 == organism:
                thispacbporf = pacbporf
                othernode = (org2,orf2)
            elif org2 == organism:
                thispacbporf = pacb.swap_query_and_sbjct(pacbporf)
                othernode = (org1,orf1)
            else:
                continue
            # loop over the original alignment positions
            for posid in range(thispacbporf._original_alignment_pos_start,thispacbporf._original_alignment_pos_end):
                posobj = thispacbporf._positions[posid]
                # check the status of this aligned position
                # options are identity, similarity, nonthing and gap
                if posobj.isa == 'similarity':
                    area[posobj.query_pos] += similarity_score
                elif posobj.isa == 'identity':
                    area[posobj.query_pos] += identity_score
                elif posobj.isa == 'gap':
                    area[posobj.query_pos] += gap_score
                else:
                    pass
                if this2others.has_key(othernode):
                    this2others[othernode][posobj.sbjct_pos] = posobj.query_pos
                else:
                    this2others[othernode] = { posobj.sbjct_pos: posobj.query_pos }

        for (orgX,orfX), otherrange in this2others.iteritems():
            for ((a,b,c,d),(org1,orf1),(org2,orf2)), pacbporf in self.pacbps.iteritems():
                if (org1,orf1) == (orgX,orfX) and org2 != organism:
                    thispacbporf = pacbporf
                else:
                    continue
                # loop over the original alignment positions
                for posid in range(thispacbporf._original_alignment_pos_start,thispacbporf._original_alignment_pos_end):
                    posobj = thispacbporf._positions[posid]
                    if posobj.query_pos in otherrange.keys():
                        orig_query_pos = otherrange[posobj.query_pos]
                        # increase counter if there is a match on this position
                        if posobj.isa == 'similarity':
                            area[orig_query_pos] += similarity_score
                        elif posobj.isa == 'identity':
                            area[orig_query_pos] += identity_score
                        elif posobj.isa == 'gap':
                            area[orig_query_pos] += gap_score
                        else:
                            pass


        # devide by the maximum number of counts possible
        combis = self.crosscombinationcount()
        for k in area.keys(): area[k] = float(area[k])/float(combis)
        # set this data structure to the object
        self._pasc_per_organism[organism] = area
        # return the scores for this CodingBlock
        return area

    # end of function set_positional_alignment_similarity_count


    def get_gff_of_positional_alignment_similarity_count(self,organism=None,fref='FREF'):
        """
        """
        # get dict with positional_alignment_similarity_count data
        simposcounts = self.get_positional_alignment_similarity_count(organism=organism)
        # get the orf belonging to this organism of this CodingBlockGraph            
        theorf = self.get_orfs_of_graph(organism=organism)[0]
        gffs = []
        for simpos,score in simposcounts.iteritems():
            dnapos = theorf.aapos2dnapos(simpos)
            gff = ( fref, 'abfgp', 'AbfgpMultipleAlignment', dnapos+1, dnapos+3, score, '+', '.', "AMA %s" % fref )
            gffs.append(gff)
        gffs.sort()
        return gffs

    # end of function get_gff_of_positional_alignment_similarity_count

    ############################################################################
    #### Functions for CBG operations on spanningrange difference           ####
    ############################################################################

    def extend_on_left_spanningrange_difference(self,
        sprdif_min_aa_length=CBG_MIN_AA_LENGTH):
        """
        Extend CBG on left spanningrange difference

        @type  sprdif_min_aa_length: positive integer
        @param sprdif_min_aa_length: minimum AA length of sprdif to be extended

        @attention: see codingblock_extension.extend_on_spanningrange_difference
                    for documentation
        """
        return codingblock_extension.extend_on_spanningrange_difference(self,
                side='left',sprdif_min_aa_length=sprdif_min_aa_length)

    # end of function extend_on_left_spanningrange_difference
 

    def extend_on_rigth_spanningrange_difference(self,
        sprdif_min_aa_length=CBG_MIN_AA_LENGTH):
        """
        Extend CBG on rigth spanningrange difference

        @type  sprdif_min_aa_length: positive integer
        @param sprdif_min_aa_length: minimum AA length of sprdif to be extended

        @attention: see codingblock_extension.extend_on_spanningrange_difference
                    for documentation
        """
        return codingblock_extension.extend_on_spanningrange_difference(self,
                side='rigth',sprdif_min_aa_length=sprdif_min_aa_length)

    # end of function extend_on_rigth_spanningrange_difference


    def split_codingblock_on_spanningrange_difference(self,side='left',
        sprdif_min_aa_length=CBG_MIN_AA_LENGTH,
        sprdif_min_node_count=CBG_SPRDIF_MIN_NODE_COUNT):
        """
        Split CBG on difference between maximum and minimum spanning range

        @attention: see codingblock_splitting.
                split_codingblock_on_spanningrange_difference for documentation
        
        @type  side: string
        @param side: one of ['left','rigth','5p','3p']; the side to split the CBG upon

        @type  sprdif_min_node_count: positive integer
        @param sprdif_min_node_count: minimum number of nodes that must have a
                                      sprdif in order to be splitted

        @type  sprdif_min_aa_length: positive integer
        @param sprdif_min_aa_length: minimum AA length of the sprdif in order
                                     to be splitted

        @rtype:  list
        @return: list of one (no split), two (succesfull split) CBGs,
                 or three (a CBG, a lsrCBG and a CBG) CBGs
        """
        return codingblock_splitting.split_codingblock_on_spanningrange_difference(
                self, side=side, sprdif_min_node_count=sprdif_min_node_count,
                sprdif_min_aa_length=sprdif_min_aa_length)

    # end of function split_codingblock_on_spanningrange_difference


    def iteratively_split_codingblock_on_spanningrange_difference(self,
        side='left',
        sprdif_min_aa_length=CBG_MIN_AA_LENGTH,
        sprdif_min_node_count=CBG_SPRDIF_MIN_NODE_COUNT):
        """
        Iteratively split CodingBlockGraph on difference between maximum and minimum spanning range

        @attention: see codingblock_splitting.
                    iteratively_split_codingblock_on_spanningrange_difference
                    for documentation

        @type  side: string
        @param side: one of ['left','rigth','5p','3p']; the side to split the CBG upon

        @type  sprdif_min_node_count: positive integer
        @param sprdif_min_node_count: minimum number of nodes that must have a
                                      sprdif in order to be splitted

        @type  sprdif_min_aa_length: positive integer
        @param sprdif_min_aa_length: minimum AA length of the sprdif in order
                                     to be splitted

        @rtype:  list
        @return: list of one (no split) or more (succesfull split) CBGs
        """
        return codingblock_splitting.iteratively_split_codingblock_on_spanningrange_difference(
                self, side=side, sprdif_min_node_count=sprdif_min_node_count,
                sprdif_min_aa_length=sprdif_min_aa_length)
        
    # end of function split_codingblock_on_spanningrange_difference


    def has_left_spanningrange_difference(self,
        sprdif_min_node_count=CBG_SPRDIF_MIN_NODE_COUNT,
        sprdif_min_aa_length=CBG_MIN_AA_LENGTH):
        """
        Has this CodingBlockGraph a SPRDIF on the left/5p/upstream side?

        @attention: see codingblock_splitting.has_left_spanningrange_difference
        
        @rtype:  Boolean
        @return: True or False
        """
        return codingblock_splitting.has_left_spanningrange_difference(self,
                sprdif_min_node_count=sprdif_min_node_count,
                sprdif_min_aa_length=sprdif_min_aa_length)

    # end of function has_left_spanningrange_difference


    def has_rigth_spanningrange_difference(self,
        sprdif_min_node_count=CBG_SPRDIF_MIN_NODE_COUNT,
        sprdif_min_aa_length=CBG_MIN_AA_LENGTH):
        """
        Has this CodingBlockGraph a SPRDIF on the left/5p/upstream side?

        @attention: see codingblock_splitting.has_left_spanningrange_difference

        @rtype:  Boolean
        @return: True or False
        """
        return codingblock_splitting.has_rigth_spanningrange_difference(self,
                sprdif_min_node_count=sprdif_min_node_count,
                sprdif_min_aa_length=sprdif_min_aa_length)

    # end of function has_rigth_spanningrange_difference


    def left_spanningrange_difference(self,correct_sprdif_all_nodes=True):
        """
        Returns the left spanningrange difference of this CodingBlockGraph

        @attention: see codingblock_splitting.spanningrange_difference

        @rtype:  dict
        @return: dict with nodes as keys and coordinate tuples of
                 deviating spanning range difference as values
        """
        return codingblock_splitting.spanningrange_difference(
                self,'left',correct_sprdif_all_nodes=correct_sprdif_all_nodes)

    # end of function left_spanningrange_difference


    def rigth_spanningrange_difference(self,correct_sprdif_all_nodes=True):
        """
        Returns the rigth spanningrange difference of this CodingBlockGraph

        @attention: see codingblock_splitting.spanningrange_difference

        @rtype:  dict
        @return: dict with nodes as keys and coordinate tuples of
                 deviating spanning range difference as values
        """
        return codingblock_splitting.spanningrange_difference(
                self,'rigth',correct_sprdif_all_nodes=correct_sprdif_all_nodes)

    # end of function rigth_spanningrange_difference


    def get_left_spanningrange_difference(self,**kwargs):
        """ Alias for self.left_spanningrange_difference() """
        return self.left_spanningrange_difference(**kwargs)
    # end of function get_left_spanningrange_difference


    def get_rigth_spanningrange_difference(self,**kwargs):
        """ Alias for self.rigth_spanningrange_difference() """
        return self.rigth_spanningrange_difference(**kwargs)
    # end of function get_rigth_spanningrange_difference

    ############################################################################
    #### Functions measuring the distance between CBGs (in a GSG)           ####
    ############################################################################

    def distance_between_codingblocks(self,other,organism=None):
        """
        Distance in AA between two CodingBlockGraphs

        @type  organism: *
        @param organism: Organism identifier (or None)

        @rtype:  dictionary (or integer if organism or node is specified)
        @return: dictionary with organisms (keys) and AA-distance between CBGs
                 (values), or only a AA-distance if an Organism identifier
                 was specified

        @attention: other CBG is supposed to be 3p/rigth of 'self' (current CBG)
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


    def omsr_distance_between_codingblocks(self,other,organism=None):
        """
        @attention: alias for distance_between_codingblocks() function
        """
        return self.distance_between_codingblocks(other,organism=organism)

    # end of function omsr_distance_between_codingblocks


    def maxsr_distance_between_codingblocks(self,other,organism=None):
        """
        Distance in AA between two CodingBlockGraphs measured on their maxsr

        @type  organism: *
        @param organism: Organism identifier (or None)

        @rtype:  dictionary (or integer if organism or node is specified)
        @return: dictionary with organisms (keys) and AA-distance between CBGs
                 (values), or only a AA-distance if an Organism identifier
                 was specified

        @attention: other CBG is supposed to be 3p/rigth of 'self' (current CBG)
        """

        # get maximal spanning ranges
        maxsrSelf  = self.maximal_spanning_range()
        maxsrOther = other.maximal_spanning_range()
        distances = {}
        for org in self.organism_set().intersection(other.organism_set()):
            if organism and org != organism: continue
            nodeA = self.node_by_organism(org)
            nodeB = other.node_by_organism(org)
            # although other is supposed to be 3p of self, take both
            # orientations into account.
            distA = min( maxsrSelf[nodeA] ) - max( maxsrOther[nodeB] ) - 1
            distB = min( maxsrOther[nodeB] ) - max( maxsrSelf[nodeA] ) - 1
            distances[org] = max( [ 0, distA, distB ] )

        # return distance only of a specific organism is requested for
        if organism:
            return distances[organism]
        else:
            return distances

    # end of function maxsr_distance_between_codingblocks

# end of class CodingBlockGraph


################################################################################
#### Helper functiond for CodingBlockGraph class                            ####
################################################################################

def recover_missing_node(subgraph,crossdata,minimal_required_edges=2):
    """
    Recover a missed node that completes the graph.

    @attention: crossdata dict must contain Pacbps in some rejected keys!
    
    @type  subgraph: PacbpCollectionGraph or CodingBlockGraph
    @param subgraph: PacbpCollectionGraph or CodingBlockGraph

    @type  crossdata: dict
    @param crossdata: crossdata <dict data structure>

    @rtype:  []
    @return: Graph objects with additional edges recruited from crossdata

    @rtype  crossdata: dict
    @return crossdata: crossdata <dict data structure>
    """

    # only perform this function if graph is incomplete BUT fully connected
    if subgraph.connectivitysaturation() != 1.0:
        return ( [ subgraph ], crossdata ) 

    # find which/if organism is missing
    allorganisms = Set()
    missing_org  = None
    for (geneQ,geneS) in crossdata.keys(): allorganisms.update([geneQ,geneS])
    if allorganisms.difference(subgraph.organism_set()):
        missing_org = list(allorganisms.difference(subgraph.organism_set()))[0]
    else:
        return ( [ subgraph ], crossdata )

    elegiable_nodes = {}
    all_potential_nodes = Set()
    potential_pacbps = {}
    for org in subgraph.organism_set():
        elegiable_nodes[org] = Set()
        potential_pacbps[org] = {}

    rejected_keys = [
        'lowscoring_pacbs',
        'rejected_pacbs_wp',    # wrong periodicity
        'rejected_pacbs_ls',    # weaker similarity out-of-frame
        #'rejected_pacbs_lc'    # low connectivity is not filled yet!!
        #'rejected_pacbs_nl',   # non-linear
        #'rejected_pacbs_aa',   # alternative alignment
    ]

    # loop over the organism/gene combinations
    # and gather all nodes that have an edge with
    # any of the known nodes
    for (geneQ,geneS) in crossdata.keys():
        if missing_org not in (geneQ,geneS):
            continue
        for rej_key in rejected_keys:
            for (a,b,orfQ,orfS) in crossdata[(geneQ,geneS)][rej_key].keys():
                if missing_org == geneQ:
                    if (geneS,orfS) in subgraph.get_nodes():
                        elegiable_nodes[geneS].add(orfQ)
                        all_potential_nodes.add(orfQ)
                        potential_pacbps[geneS][(a,b,orfQ,orfS)] =\
                            ( (geneQ,geneS),
                              crossdata[(geneQ,geneS)][rej_key][(a,b,orfQ,orfS)]
                            )
                else:
                    if (geneQ,orfQ) in subgraph.get_nodes():
                        elegiable_nodes[geneQ].add(orfS)
                        all_potential_nodes.add(orfS)
                        potential_pacbps[geneQ][(a,b,orfQ,orfS)] =\
                            ( (geneQ,geneS),
                              crossdata[(geneQ,geneS)][rej_key][(a,b,orfQ,orfS)]
                            )

    # check how often nodes are observed
    retgraphs = []
    for node in all_potential_nodes:
        cnt = sum([node in vset for vset in elegiable_nodes.values() ])
        if cnt >= minimal_required_edges:
            newgra  = deepcopy(subgraph)
            newnode = (missing_org,node)
            newgra.add_node(newnode)
            for org in elegiable_nodes.keys():
                if node in elegiable_nodes[org]:
                    knownnode = newgra.get_organism_nodes(org)[0]
                    for (a,b,orfQ,orfS) in potential_pacbps[org].keys():
                        if (newnode[1],knownnode[1]) in [(orfQ,orfS),(orfS,orfQ)]: 
                            combi,pacbp = potential_pacbps[org][(a,b,orfQ,orfS)]
                            if newgra.has_edge(newnode, knownnode):
                                _wt = newgra.weights[( newnode, knownnode )]
                                if pacbp.bitscore > _wt:
                                    newgra.set_edge_weight(
                                        newnode, knownnode, wt=pacbp.bitscore )
                            else:
                                newgra.add_edge(
                                    newnode, knownnode, wt=pacbp.bitscore )
                            # update the crossdata dictionary
                            crossdata[combi]['accepted_pacbs'][(a,b,orfQ,orfS)] = pacbp
            # harvest the pacbps from crossdata
            newgra.harvest_pacbps_from_crossdata(crossdata)
            # split on alternative pacbps if created
            alternatives =\
                newgra.split_codingblock_on_alternatives_in_pacbps_dict(
                    filter_for_msr=True,
                    filter_for_omsr=True,
                    )
            # and append to the list of graphs to return
            for alternative in alternatives:
                retgraphs.append(alternative)

    # return list AND updated crossdata
    return (retgraphs, crossdata)

# end of function recover_missing_node


def make_pacbps_for_missing_edges(gra,use_maxsr=False,aa_extra_offset=1):
    """
    """
    if gra.has_overall_minimal_spanning_range():
        omsr = gra.overall_minimal_spanning_range()
        has_omsr = True
    else:
        # use maxsr in stead of omsr to create edges
        omsr = gra.maximal_spanning_range()
        has_omsr = False

    if use_maxsr and has_omsr:
        sprL = gra.left_spanningrange_difference(
                correct_sprdif_all_nodes=False)
        sprR = gra.rigth_spanningrange_difference(
                correct_sprdif_all_nodes=False)

    for (node1,node2) in gra.pairwisecrosscombinations_node():
        # check if these are nodes present as an edge
        if gra.has_edge(node1,node2): continue
        # Check if noth nodes are present as omsr keys
        # This can occur when has_omsr=False and maxsr
        # is in fact a mess as well.
        if not omsr.has_key(node1): continue
        if not omsr.has_key(node2): continue

        # if not, then start makeing a Pacbp!
        # get ranges, orgs, orfs, positions and sequences
        org1 = gra._organism_from_node(node1)
        org2 = gra._organism_from_node(node2)
        orf1 = gra.get_orfs_of_graph(organism=org1)[0]
        orf2 = gra.get_orfs_of_graph(organism=org2)[0]

        # get the coords of the to-be-aligned sequence pair from OMSR
        aa1start = min(omsr[node1])-aa_extra_offset
        aa1end   = max(omsr[node1])+1+aa_extra_offset
        aa2start = min(omsr[node2])-aa_extra_offset
        aa2end   = max(omsr[node2])+1+aa_extra_offset

        if use_maxsr and has_omsr:
            # extend the to-be-created alignment to the MAXSR
            if sprL.has_key(node1) and sprL.has_key(node2):
                extensionL = min([len(sprL[node1]),len(sprL[node2])])
                aa1start-=extensionL
                aa2start-=extensionL
            if sprR.has_key(node1) and sprR.has_key(node2):
                extensionR = min([len(sprR[node1]),len(sprR[node2])])
                aa1end+=extensionR
                aa2end+=extensionR

        # check if requested coordinates are not outside of the Orf's range
        # this can happen due to use_maxsr and/or aa_extra_offset
        if aa1start < orf1.protein_startPY: aa1start = orf1.protein_startPY
        if aa1end > orf1.protein_endPY:     aa1end   = orf1.protein_endPY
        if aa2start < orf2.protein_startPY: aa2start = orf2.protein_startPY
        if aa2end > orf2.protein_endPY:     aa2end   = orf2.protein_endPY

        # create headers and fetch sequences from Orf objects
        header1  = "%s_orf_%s_%s_%s" % (org1,orf1.id,aa1start,aa1end)
        header2  = "%s_orf_%s_%s_%s" % (org2,orf2.id,aa2start,aa2end)
        seq1 = orf1.getaas(abs_pos_start=aa1start,abs_pos_end=aa1end)
        seq2 = orf2.getaas(abs_pos_start=aa2start,abs_pos_end=aa2end)

        # check if sequences exist/ at least 1 AA
        if not seq1 and not seq2:
            print "Warning: ZeroProteinSequenceLengthException",
            print "query+", aa1start, aa1end
            print "Warning: ZeroProteinSequenceLengthException",
            print "sbjct+", aa2start, aa2end
            # warning is printed, now continue
            continue
        elif not seq2:
            print "Warning: ZeroProteinSequenceLengthException",
            print "sbjct", aa2start, aa2end
            # warning is printed, now continue
            continue
        elif not seq1:
            print "Warning: ZeroProteinSequenceLengthException",
            print "query", aa1start, aa1end
            # warning is printed, now continue
            continue
        else:
            pass

        # align the sequences with clustalw
        seqs = { header1: seq1, header2: seq2 }
        (alignedseqs,alignment) = clustalw(seqs=seqs)

        # make pacbp from clustalw alignment
        pacbp = pacb.conversion.pacbp_from_clustalw(
                    alignment=(
                            alignedseqs[header1],
                            alignment,
                            alignedseqs[header2]
                            ),
                    coords=(aa1start,aa1end,aa2start,aa2end)
                    )

        if pacbp:
            # make & extend PacbPORF
            pacbporf   = pacb.PacbPORF(pacbp,orf1,orf2)
            pacbporf.extend_pacbporf_after_stops()
            # add edge & pacbporf to CBG
            gra.add_edge(node1,node2,wt=pacbporf.bitscore)
            key = pacbporf.construct_unique_key(node1,node2)
            gra.pacbps[(key,node1,node2)] = pacbporf
        else:
            # pacbp.conversion.pacbp_from_clustalw did
            # not yield any proper alignment
            pass

# end of function make_pacbps_for_missing_edges


def estimate_bitscores_ksminxcbg(ksminxcbg,gtg):
    """
    """
    missing_esges = {}
    diforgnodes = gtg.different_organisms(ksminxcbg)
    omsrseqs = ksminxcbg.getomsrproteinsequences()
    # get aa similarity matrix from (arbitrary pacbp)
    matrix = ksminxcbg.pacbps.values()[0].MATRIX
    for o1,o2 in gtg.pairwisecrosscombinations_organism():
        if o1 in diforgnodes and o2 in diforgnodes:
            # Both organisms not in ksminxcbg. Extrapolation
            # of bitscore is not possible!
            pass
        elif o1 in diforgnodes:
            seq = omsrseqs[o2]
            bitSeq = float(matrix.scorealignment(seq,seq))
            missing_esges[(o1,o2)] = gtg._bitscore_ratios[(o1,o2)] * bitSeq
        elif o2 in diforgnodes:
            seq = omsrseqs[o1]
            bitSeq = float(matrix.scorealignment(seq,seq))
            missing_esges[(o1,o2)] = gtg._bitscore_ratios[(o1,o2)] * bitSeq
        else:
            pass
    # return estimated bitscores for missing edges
    return missing_esges

# end of function estimate_bitscores_ksminxcbg

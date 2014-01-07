################################################################################
### InwardsPointingCodingBlockGraph class                                   ####
################################################################################

# graphAbgp imports
from graph_pacbpcollection import PacbpCollectionGraph
from exceptions import *
from codingblock.sequenceretrieval import CodingBlockGraphSequenceRetievalFunctions
from codingblock.MINSR import MinimalSpanningRange
from codingblock.MAXSR import MaximalSpanningRange
from codingblock.OMSR import OverallMinimalSpanningRange
import codingblock.spanningrangedifference as SPRDIF
import codingblock_splitting

# Pacb imports
from pacb.validators import IsPacbPORF
from pacb.exceptions import NoPacbPORF
from pacb.ordering import relpos2binaryrelpos

# Abgp imports
from lib_clustalw import clustalw, strip_alignment_for_exterior_gaps
from lib_cexpander import runcexpander
from pythonlibs.uniqueness import get_random_string_tag
from pythonlibs.statistics import meanstdv

# Python imports
from sets import Set
from copy import deepcopy

# Import Global variables and settings
from settings.executables import TCODE_MIN_CODING,TCODE_MAX_NONCODING
from settings.codingblockgraph import *
from settings.genestructure import (
        FIRST_ANNOTATED_EXON_LABEL,
        FINAL_ANNOTATED_EXON_LABEL,
        IS_ANNOTATED_EXON_LABEL,
        ORF_IS_UNIGENE_LABEL,
        )

TCODE_SINGLE_MAX_NONCODING   = 0.690
TCODE_SINGLE_MIN_CODING      = TCODE_MIN_CODING    # 0.950
TCODE_COMBI_MAX_NONCODING    = TCODE_MAX_NONCODING # 0.740

NT_RATIO_SINGLE_MAX_NONCODING= 0.85
NT_RATIO_COMBI_MAX_NONCODING = 0.95
NT_RATIO_SINGLE_MIN_CODING   = 1.00



class InwardsPointingCodingBlockGraph(PacbpCollectionGraph,
    MinimalSpanningRange,
    MaximalSpanningRange,
    OverallMinimalSpanningRange,
    CodingBlockGraphSequenceRetievalFunctions):
    """
    InwardsPointingCodingBlockGraph (IWPCBG) class, inheriting from PacbpCollectionGraph class.
    """

    def __init__(self):
        """
        Initialize a CodingBlockGraph
        """
        # Initialize from PacbpCollectionGraph
        PacbpCollectionGraph.__init__(self)

        # short name tag for printing
        self._short_name = "IWPCBG"


        # some GLOBAL variables that are required in class functions
        #self.ALTERNATIVE_ALIGNMENT_OVERLAP_RATIO = alternative_alignment_overlap_ratio

        # 08/04/09 smallest omsr size must be 3 (means 4 AAs overall aligned)
        # this is defined in settings.codingblockgraph
        #self.MINIMAL_OVERAL_SPANNING_RANGE_SIZE = minimal_overal_spanning_range_size 

        # attribute needed to set the GeneTreeGraph into
        #self._GENETREE      = None

        # attributes to cache OMSR and related data structures
        self._omsr  = {}    # Overall Minimal Spanning Range
        self._msr   = {}    # Minimal Spanning Range
        self._maxsr = {}    # MAXimal Spanning Range

    # end of function __init__


    def _get_target_organism(self):
        """ """
        targetNode = self._get_target_node()
        if not targetNode:
            return None
        else:
            return self.organism_by_node(targetNode)

    # end of function _get_target_organism


    def _get_target_node(self):
        if self.pacbps:
            targetNode = self.pacbps.keys()[0][1]
        else:
            targetNode = None
            nodeCnts = dict( [ (node,0) for node in self.get_nodes() ] )
            for nodeQ,nodeS in self.weights.iteritems():
                nodeCnts[nodeQ] +=1
                nodeCnts[nodeS] +=1
            max_cnt = max(nodeCnts.values())
            for node,cnt in nodeCnts.iteritems():
                if cnt == max_cnt:
                    targetNode = node
                    break
        # return targetNode
        return targetNode

    # end of function _get_target_node


    ############################################################################
    #### Basic Function OVERRIDES                                           ####
    ############################################################################

    def pairwisecrosscombinations_node(self):
        """ OVERRIDES PacbpCollectionGraph.pairwisecrosscombinations_node() """
        if not self.weights:
            return []
        # get targetNode -> make pairwise combinations!
        targetNode = self._get_target_node()
        combis = []
        for node in self.get_ordered_nodes():
            if node == targetNode: continue
            combis.append((targetNode,node))
        return combis

    # end of function pairwisecrosscombinations_node


    def has_edge(self,u,v):
        """ OVERRIDES GraphPlus.has_edge() """
        if len(self.nodes) == 2 and ( self.weights.keys() == [(u,v),(v,u)] or self.weights.keys() == [(v,u),(u,v)] ):
            # hmm.... escape. Weights are okay, nut nodes not.
            # Temporarily allow this; it can happen in InwardsPointingCodingblockGraphs
            return True
        else:
            return PacbpCollectionGraph.has_edge(self,u,v)

    # end of function has_edge


    def barcode(self):
        """ Get a semi-unique barcode string of this CBG """
        return "%s_%s_%s" % ( self._short_name, self.total_weight(),
            "_".join([ "%s_%s" % ( node[0],node[1]) for node in self.get_ordered_nodes() ])
            )

    # end of function barcode


    def _update_after_changes(self):
        """
        @attention: MUST BE OVERRIDDEN IN THE SUBCLASSES, if the returned value is not correct
        """
        self._update_pacbps()
        self.clear_cache()

    # end of function _update_after_changes


    def clear_cache(self):
        """
        Cleanup the CBGs cashed attributes; after an update they are False!
        """
        self._GENETREE = None
        self._omsr  = {}
        self._msr   = {}
        self._maxsr = {}

    # end of function clear_cache


    def create_cache(self):
        """
        Create cached attributes of this CodingBlockGraph; saves a lot processing time
        """
        self.clear_cache()
        self._msr   = self.minimal_spanning_range()
        self._maxsr = self.maximal_spanning_range()
        #self._omsr  = self.overall_minimal_spanning_range()
        
    # end of function create_cache


    ############################################################################
    #### OMSR/MAXSR/MINSR Functions OVERRIDES                               ####
    ############################################################################

    def minimal_spanning_range(self,organism=None,node=None):
        """ OVERRIDES MINSR.minimal_spanning_range() """
        if organism and organism not in self.organism_set():
            raise OrganismNotPresentInGraph
        if node and node not in self.get_nodes():
            raise NodeNotPresentInGraph

        # perform basal MINSR.minimal_spanning_range() function
        msr = MinimalSpanningRange.minimal_spanning_range(self)

        # second iteration: update informant MINSR area
        for (nodeQ,nodeS) in self.pairwisecrosscombinations_node():
            # check if these nodes are indeed present as an edge
            if not self.has_edge(nodeQ,nodeS): continue
            minsrQ = msr[nodeQ]
            pacbporf = self.get_pacbps_by_nodes(nodeQ,nodeS)[0]
            minsrSstart = pacbporf.aapos_query2sbjct( min(minsrQ) )
            minsrSend   = pacbporf.aapos_query2sbjct( max(minsrQ) )
            msr[nodeS] = Set(range(minsrSstart,minsrSend+1))

        # return ranges or the range of a specific organism or node
        if organism:
            return msr[self.node_by_organism(organism)]
        elif node:
            return msr[node]
        else:
            return msr

    # end of function minimal_spanning_range


    def overall_minimal_spanning_range(self,**kwargs):
        """ OVERRIDES OMSR.overall_minimal_spanning_range() """
        return self.minimal_spanning_range(**kwargs)

    # end of function overall_minimal_spanning_range


    def has_overall_minimal_spanning_range(self,**kwargs):
        """ OVERRIDES OMSR.has_overall_minimal_spanning_range() """
        return self.has_minimal_spanning_range(**kwargs)

    # end of function has_overall_minimal_spanning_range

    ############################################################################
    #### SPRDIF Functions                                                   ####
    ############################################################################

    def has_left_spanningrange_difference(self,**kwargs):
        """ @attention: codingblock.spanningrangedifference for documentation"""
        return SPRDIF.has_left_spanningrange_difference(self,**kwargs)
    # end of function has_left_spanningrange_difference


    def has_rigth_spanningrange_difference(self,**kwargs):
        """ @attention: codingblock.spanningrangedifference for documentation"""
        return SPRDIF.has_rigth_spanningrange_difference(self,**kwargs)
    # end of function has_rigth_spanningrange_difference


    def left_spanningrange_difference(self,**kwargs):
        """ @attention: codingblock.spanningrangedifference for documentation"""
        return SPRDIF.spanningrange_difference(self,side='left',**kwargs)
    # end of function left_spanningrange_difference


    def rigth_spanningrange_difference(self,**kwargs):
        """ @attention: codingblock.spanningrangedifference for documentation"""
        return SPRDIF.spanningrange_difference(self,side='rigth',**kwargs)
    # end of function rigth_spanningrange_difference


    def sources(self):
        """
        Get a dict with counts of sources of PacbP(ORF)s

        @rtype:  {} 
        @return: dict with sources (keys) and occurence counts (values) of pacbp.source of PacbP(ORF)s 

        @attention: see pacbp/__init__.py for recognized sources (blastp,unknown,hmmsearch,clustalw,clustalw-EXTENDED,clustalw-OPTIMIZED,lsrPACBP...)
        """
        ret = {}
        for pacbp in self.pacbps.values():
            if ret.has_key(pacbp.source):
                ret[pacbp.source]+=1
            else:
                ret[pacbp.source] =1
        return ret

    # end of function sources


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


    def are_all_edges_covered_by_pacbps(self):
        """
        """
        if len(self.pacbps) < self.edge_count():
            return False
        for (node1,node2) in self.pairwisecrosscombinations_node():
            if self.has_edge(node1,node2):
                pacbps = self.get_pacbps_by_nodes(node1=node1,node2=node2)
                if not pacbps:
                    return False
            else:
                pass
        # EOF edges reached; all edges covered by a pacbp
        return True

    # end of function are_all_edges_covered_by_pacbps


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


    def have_all_starts(self):
        """
        Do all orfs in the pacbps have start codons?

        @attention: only functional if pacbp's are present in the graph
                    (perform the function `harvest_pacbps_from_crossdata`)

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
        Are there (putative) startcodons (Methionines) upstream of the OMSR in all orfs?

        @type  omsr_offset: integer
        @param omsr_offset: offset into 3' side of OMSR to take into account as well

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


    ############################################################################
    #### Functions concerning ReversecomplementCodingBlockGraphs            ####
    ############################################################################
    
    def is_reversecomplement(self,revcbg=None):
        """
        """
        return False
    
    # end of function is_reversecomplement
    
    
    def get_reversecomplement(self):
        """
        """
        return None
    
    # end of function is_reversecomplement 


    def get_frameshifted_cbg(self,input,verbose=False):
        """
        """
        return None

    # end of function get_frameshifted_cbg

    ############################################################################
    #### etcetera Functions                                                 ####
    ############################################################################

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
        omsr = self.minimal_spanning_range()
        distances = {}
        for thenode in self.get_nodes():
            org = self.organism_by_node(thenode)
            if organism and org!=organism: continue
            theorf = self.get_orfs_of_graph(organism=org)[0]
            print thenode,min(omsr[thenode]),max(omsr[thenode]), side, theorf
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


    def nt_spacing_between_codingblocks(self,others):
        """
        Spacing (in nt) between two CodingBlockGraphs

        @attention: distance between identical PacbPORFs is set to 0nt

        @rtype:  {} 
        @return: dict with Organism Identifiers (keys) nt distances (values)
        """
        distdict = {}
        maxsr = dict([ (self.organism_by_node(n), r) for n,r in self.maximal_spanning_range().iteritems() ])
        for other in others:
            _maxsr = dict([ (other.organism_by_node(n), r) for n,r in other.maximal_spanning_range().iteritems() ])
            for organism,setrange in _maxsr.iteritems():
                if distdict.has_key(organism): continue
                if not maxsr.has_key(organism): continue
                # calculate distance
                if not setrange.symmetric_difference(maxsr[organism]):
                    distdict[organism] = 0
                else:
                    distdict[organism] = max([0,(min(setrange) - max(maxsr[organism])) * 3])

        return distdict

    # end of function nt_spacing_between_codingblocks


    def tcode_spacing_between_codingblocks(self,others):
        """
        TCODE `Spacing` between two CodingBlockGraphs

        @attention: identical PacbPORFs are omitted

        @rtype:  {} 
        @return: dict with Organism Identifiers (keys) TCODE scores (values)
        """

        if type(others) not in (type([]),type(())): others = [ others ]
        distdict = {}
        target = self._get_target_organism()
        maxsr = dict([ (self.organism_by_node(n), r) for n,r in self.maximal_spanning_range().iteritems() ])
        orfObjs = self.get_orfs_of_graph()

        for other in others:
            _maxsr = dict([ (other.organism_by_node(n), r) for n,r in other.maximal_spanning_range().iteritems() ])
            for organism,setrange in _maxsr.iteritems():
                if distdict.has_key(organism): continue
                if not maxsr.has_key(organism): continue
                # calculate distance
                if not setrange.symmetric_difference(maxsr[organism]):
                    # identical PacbPORFS
                    distdict[organism] = None
                else:
                    prev_dna_end = max(maxsr[organism])*3
                    next_dna_sta = min(setrange)*3
                    if next_dna_sta <= prev_dna_end:
                        # no distance left to calculate TCODE score for
                        distdict[organism] = None
                    else:
                        orfObj = orfObjs[organism][0]
                        tcode = orfObj.find_lowest_scoring_tcode_stretch(
                                orfObj._RAW_TCODE_DATA,
                                prev_dna_end,next_dna_sta)
                        distdict[organism] = tcode


        # remove None values in distdict
        while None in distdict.values():
            for k,v in distdict.iteritems():
                if v == None:
                    del( distdict[k] )
                    break
        # return TCODE distance dictionary
        return distdict

    # end of function tcode_spacing_between_codingblocks


    def distance_between_codingblocks(self,other,**kwargs):
        """
        Distance in AA between two CodingBlockGraphs

        @attention: alias funcion name for omsr_distance_between_codingblocks()
        """
        return self.omsr_distance_between_codingblocks(other,**kwargs)

    # end of function distance_between_codingblocks


    def organisms_with_different_orfs(self,other):
        """ """
        return [ self.organism_by_node(node) for node in self.different_nodes(other) ]
    # end of function organisms_with_different_orfs 


    def multiplealignment(self):
        """ """
        # get sequences & coordinated and rewrite Nodes to Organism identifiers
        seqs,coords = self.get_maxsr_proteinsequences_and_coords()
        coords = dict([ (self.organism_by_node(node),[min(vlist),max(vlist)+1]) for node,vlist in coords.iteritems() ])
        seqs   = dict([ (self.organism_by_node(node),seq) for node,seq in seqs.iteritems() ])

        # align sequences with ClustalW
        (alignedseqs,alignment) = clustalw( seqs= seqs )
        # trim alignment for leading & trailing gaps
        alignedseqs,alignment,coords = strip_alignment_for_exterior_gaps(alignedseqs,alignment,coords)
        # return single string of multilined fasta
        return "\n".join([">%s_orf_%s\n%s" % (k,self.node_by_organism(k)[1],v) for k,v in alignedseqs.iteritems()])
    # end of function multiplealignment


    def multiplealignment_and_coords(self):
        """ """
        # get sequences & coordinated and rewrite Nodes to Organism identifiers
        seqs,coords = self.get_maxsr_proteinsequences_and_coords()
        coords = dict([ (self.organism_by_node(node),[min(vlist),max(vlist)+1]) for node,vlist in coords.iteritems() ])
        seqs   = dict([ (self.organism_by_node(node),seq) for node,seq in seqs.iteritems() ])

        # align sequences with ClustalW
        (alignedseqs,alignment) = clustalw( seqs= seqs )
        # trim alignment for leading & trailing gaps
        alignedseqs,alignment,coords = strip_alignment_for_exterior_gaps(alignedseqs,alignment,coords)
        # return string of multilined fasta AND the coords
        fastastring = "\n".join([">%s_orf_%s\n%s" % (k,self.node_by_organism(k)[1],v) for k,v in alignedseqs.iteritems()])
        return fastastring, coords

    # end of function multiplealignment_and_coords


    def tcode(self):
        """ """
        tcode_scores_query = []
        tcode_scores_sbjct = []
        for pacbporf in self.pacbps.values():
            try:
                IsPacbPORF(pacbporf)
                tcQ = pacbporf.tcode_query()
                tcS = pacbporf.tcode_sbjct()
                if not (tcQ==0.0 or tcS==0.0):
                    # one of the informants is a unigene -> no TCODE data!
                    # TODO: make unigene detectiom here direct, not indirect
                    tcode_scores_query.append( pacbporf.tcode_query() )
                    tcode_scores_sbjct.append( pacbporf.tcode_sbjct() )
            except NoPacbPORF:
                # PacbPORFs are required here!
                tcode_scores_query = []
                tcode_scores_sbjct = []
                break
        if not tcode_scores_query or not tcode_scores_sbjct:
            return None
        else:
            # Append the (maximal) TCODE query score to the sbjct
            # and devide by the total length of this list
            tcode_scores_sbjct.append( max(tcode_scores_query) )
            return sum(tcode_scores_sbjct) / len(tcode_scores_sbjct)

    # end of function tcode


    ############################################################################
    #### Functions 'moved'  from lib_pcg2blocks.print_orfidstruct to here   ####
    ############################################################################

    def has_target_orf_methionine(self):
        """ """
        targetOrg = self._get_target_organism()
        return self.get_orfs_of_graph(organism=targetOrg)[0].has_methionine()

    # end of function has_target_orf_methionine


    def has_informant_orf_methionine(self,informant):
        """ """
        return self.get_orfs_of_graph(organism=informant)[0].has_methionine()

    # end of function has_informant_orf_methionine


    def have_informant_orfs_methionines(self):
        """ """
        targetOrg = self._get_target_organism()
        for org in self.organism_set():
            if org == targetOrg: continue
            if not self.get_orfs_of_graph(organism=org)[0].has_methionine():
                return False
        else:
            return True

    # end of function have_informant_orfs_methionines


    def count_genomic_informants(self):
        """ """
        cnt = 0
        for pacbporf in self.pacbps.values(): 
            if not hasattr(pacbporf.orfS,ORF_IS_UNIGENE_LABEL):
               cnt+=1
        return cnt

    # end of function count_genomic_informants 


    def count_unigene_informants(self):
        """ """
        cnt = 0
        for pacbporf in self.pacbps.values():
            if hasattr(pacbporf.orfS,ORF_IS_UNIGENE_LABEL):
               cnt+=1
        return cnt

    # end of function count_unigene_informants 


    def count_orfs_with_methionines(self):
        """ """
        cnt=0
        for org in self.organism_set():
            if self.get_orfs_of_graph(organism=org)[0].has_methionine():
                cnt+=1
        return cnt

    # end of function count_orfs_with_methionines


    def count_orfs_with_signalpeptides(self):
        """ """
        cnt=0
        for org in self.organism_set():
            theorf = self.get_orfs_of_graph(organism=org)[0]
            if theorf._signalp_sites:
                cnt+=1
        return cnt

    # end of function count_orfs_with_signalpeptides 


    def get_signalp_score(self):
        """ """
        scores = [] 
        for org in self.organism_set():
            theorf = self.get_orfs_of_graph(organism=org)[0]
            if theorf._signalp_sites:
                scores.append( max([ sp.pssm_score for sp in theorf._signalp_sites]) ) 
        if not scores:
            return None
        else:
            return sum(scores)/len(scores) 

    # end of function count_orfs_with_signalpeptides 



    def is_target_orf_labeled_as_first_exon(self):
        """ """
        targetOrg = self._get_target_organism()
        theorf = self.get_orfs_of_graph(organism=targetOrg)[0]
        if hasattr(theorf,FIRST_ANNOTATED_EXON_LABEL):
            return True
        else:
            return False

    # end of function is_target_orf_labeled_as_first_exon


    def is_target_orf_labeled_as_final_exon(self):
        """ """
        targetOrg = self._get_target_organism()
        theorf = self.get_orfs_of_graph(organism=targetOrg)[0]
        if hasattr(theorf,FINAL_ANNOTATED_EXON_LABEL):
            return True
        else:
            return False

    # end of function is_target_orf_labeled_as_final_exon


    def is_target_orf_labeled_as_annotated_exon(self):
        """ """
        targetOrg = self._get_target_organism()
        theorf = self.get_orfs_of_graph(organism=targetOrg)[0]
        if hasattr(theorf,IS_ANNOTATED_EXON_LABEL):
            return True
        else:
            return False

    # end of function is_target_orf_labeled_as_annotated_exon


    def is_informant_orf_labeled_as_first_exon(self,informant):
        """ """
        theorf = self.get_orfs_of_graph(organism=informant)[0]
        if hasattr(theorf,FIRST_ANNOTATED_EXON_LABEL):
            return True
        else:
            return False

    # end of function is_informant_orf_labeled_as_first_exon


    def is_informant_orf_labeled_as_final_exon(self,informant):
        """ """
        theorf = self.get_orfs_of_graph(organism=informant)[0]
        if hasattr(theorf,FINAL_ANNOTATED_EXON_LABEL):
            return True
        else:
            return False

    # end of function is_informant_orf_labeled_as_final_exon


    def is_informant_orf_labeled_as_annotated_exon(self,informant):
        """ """
        theorf = self.get_orfs_of_graph(organism=informant)[0]
        if hasattr(theorf,IS_ANNOTATED_EXON_LABEL):
            return True
        else:
            return False

    # end of function is_informant_orf_labeled_as_annotated_exon


    def count_orfs_labeled_as_first_exon(self):
        """ """
        cnt=0
        for org in self.organism_set():
            theorf = self.get_orfs_of_graph(organism=org)[0]
            if hasattr(theorf,FIRST_ANNOTATED_EXON_LABEL):
                cnt+=1
        return cnt

    # end of function count_orfs_labeled_as_first_exon


    def count_orfs_labeled_as_final_exon(self):
        """ """
        cnt=0
        for org in self.organism_set():
            theorf = self.get_orfs_of_graph(organism=org)[0]
            if hasattr(theorf,FINAL_ANNOTATED_EXON_LABEL):
                cnt+=1
        return cnt

    # end of function count_orfs_labeled_as_final_exon


    def count_orfs_labeled_as_annotated_exon(self):
        """ """
        cnt=0
        for org in self.organism_set():
            theorf = self.get_orfs_of_graph(organism=org)[0]
            if hasattr(theorf,IS_ANNOTATED_EXON_LABEL):
                cnt+=1
        return cnt

    # end of function count_orfs_labeled_as_annotated_exon


    def get_orfs_labeled_as_annotated_exon(self):
        """ """
        orfdict = {}
        for org in self.organism_set():
            theorf = self.get_orfs_of_graph(organism=org)[0]
            if hasattr(theorf,IS_ANNOTATED_EXON_LABEL):
                node = (org,theorf.id)
                orfdict[node] = theorf
        return orfdict

    # end of function get_orfs_labeled_as_annotated_exon


    def print_detail(self):
        """ """
        print "##", self
        try:    tcode = "%1.2f" % self.tcode()
        except: tcode = "n.a."
        pro_nt_id = self.get_nt_identity()
        dna_nt_id = self.get_unguided_nt_identity()
        print "## CODING?: %s TCODE %s ntID %1.2f - %1.2f RATIO: %1.2f" % (
                self.is_coding(),
                tcode,
                pro_nt_id,
                dna_nt_id,
                (pro_nt_id / dna_nt_id),
                )
        for key,pacbporf in self.pacbps.iteritems():
            print "# %1.2f - %1.2f" % (
                    pacbporf.get_nt_identity(),
                    pacbporf.get_unguided_nt_identity() ),
            print pacbporf.is_coding(),
            print "%1.2f [ Q %1.1f S %1.1f ]" % (
                    pacbporf.tcode(),
                    pacbporf.tcode_query(),
                    pacbporf.tcode_sbjct() ),
            print pacbporf, (key[1][0], key[2][0])

    # end of function print_detail


    def get_exon_uniformity_score(self):
        """ """
        if not self.pacbps: return None
        if len(self.pacbps) == 1: return None
        maxsr = self.maximal_spanning_range(organism=self._get_target_organism())
        ratios = []
        for pf in self.pacbps.values():
            ratios.append( float(len(pf.alignment_protein_range_query())) / float(len(maxsr)) )
        mean,std = meanstdv(ratios)
        return std

    # end of function get_exon_uniformity_score


    def get_projected_tailing_stop_aa_difference(self):
        """ """
        distances = []
        for (key,(org1,orf1),(org2,orf2)),pacbporf in self.pacbps.iteritems():
            if not hasattr(pacbporf.orfS,ORF_IS_UNIGENE_LABEL):
                dist = pacbporf.get_projected_tailing_stop_aa_difference()
                distances.append(dist)
        if not distances:
            return None
        else:
            return int( float(sum(distances)) / len(distances) )

    # end of function get_projected_tailing_stop_aa_difference


    def get_projected_tailing_stop_nonaligned_aa_difference(self):
        """ """
        distances = []
        for (key,(org1,orf1),(org2,orf2)),pacbporf in self.pacbps.iteritems():
            if not hasattr(pacbporf.orfS,ORF_IS_UNIGENE_LABEL):
                dist = pacbporf.get_projected_tailing_stop_nonaligned_aa_difference()
                distances.append(dist)
        if not distances:
            return None
        else:
            return int( float(sum(distances)) / len(distances) )

    # end of function get_projected_tailing_stop_nonaligned_aa_difference


    def get_identityscore(self):
        """ """
        ids = [ pf.identityscore for pf in self.pacbps.values() ]
        if ids: return sum(ids)/len(ids)
        else:   return False
    # end of function get_identityscore


    def get_bitscore(self):
        """ """
        if self.pacbps:
            return sum([ pf.bitscore for pf in self.pacbps.values() ])
        else:
            return 0 
    # end of function get_bitscore


    def get_bits_per_aa(self):
        """ """
        if self.pacbps:
            return float(sum([ pf.bitscore for pf in self.pacbps.values() ])) /\
                float(sum([ pf.get_unextended_length() for pf in self.pacbps.values() ]))
        else:
           return 0.0
    # end of function bits_per_aa

    def get_nt_identity(self):
        """ """
        # hmmmmm; do we do this on MINSR or MAXSR?
        if not self.pacbps: return None
        ids = []
        for key,pacbporf in self.pacbps.iteritems():
            ids.append( pacbporf.get_nt_identity() )
        return sum(ids)/len(ids)

    # end of function get_nt_identity


    def get_unguided_nt_identity(self):
        """ """
        # hmmmmm; do we do this on MINSR or MAXSR?
        if not self.pacbps: return None
        ids = []
        for key,pacbporf in self.pacbps.iteritems():
            ids.append( pacbporf.get_unguided_nt_identity() )
        return sum(ids)/len(ids)

    # end of function get_unguided_nt_identity


    def get_cexpander_uniformly_aligned_count(self):
        """ """
        # run cexpander. TODO -> move to one place
        fname = "%s.tmp.cexpander.mfa" % get_random_string_tag()
        fh = open(fname,'w')
        for node,seq in self.getmaxsrproteinsequences().iteritems():
            fh.write( ">%s\n%s\n" % (node,seq))
        fh.close()
        # get cxpdrOutput object; file-cleanup is taken care for
        cxpdrOutput = runcexpander(fname,
                cbalignp_commandline = " -y", output='binary')

        # do cexpander binary string evaluation
        return cxpdrOutput.binarystring.count("1")

    # end of function get_cexpander_uniformly_aligned_count


    def is_coding(self):
        """ """
        # do cexpander binary string evaluation
        cnt = self.get_cexpander_uniformly_aligned_count()
        if cnt>0:   cexpandercheck = True
        else:       cexpandercheck = False

        tcode = self.tcode()
        # do tcode evalution
        if not tcode:                            pass
        elif tcode >= TCODE_SINGLE_MIN_CODING:   return True
        elif tcode < TCODE_SINGLE_MAX_NONCODING:
            if cexpandercheck == False:          return False
            else:                                pass
        else:                                    pass


        # do nt identity evalution
        pro_nt_id = self.get_nt_identity()
        dna_nt_id = self.get_unguided_nt_identity()
        ratio     = pro_nt_id / dna_nt_id

        if ratio < NT_RATIO_SINGLE_MAX_NONCODING:
            if cexpandercheck == False:             return False
            else:                                   pass
        elif ratio >= NT_RATIO_SINGLE_MIN_CODING:   return True
        else:                                       pass

        # now do combinatory check.
        if cexpandercheck == True:
            if ratio >= NT_RATIO_COMBI_MAX_NONCODING:   return True
            elif tcode >= TCODE_COMBI_MAX_NONCODING:    return True
            else:                                       return False
        else:
            if ratio < NT_RATIO_COMBI_MAX_NONCODING:    return False
            elif tcode < TCODE_COMBI_MAX_NONCODING:     return False
            else:                                       return True

    # end of function is_coding


    def get_average_upstream_methionine_pssm_score(self):
        """ """
        tssobjs = [ pf.get_sbjct_upstream_tss() for pf in self.pacbps.values() ]
        tssobjs.append( self.pacbps.values()[0].get_query_upstream_tss() )
        while None in tssobjs: tssobjs.remove(None)
        if not tssobjs:
            return 0.0
        else:
            return sum([ tss.pssm_score for tss in tssobjs ])/len(tssobjs)

    # end of function get_average_upstream_methionine_pssm_score


    def scan_orfs_for_pssm_tss(self,**kwargs):
        """ """
        for pf in self.pacbps.values():
            pf.orfQ.scan_orf_for_pssm_tss(**kwargs)
            pf.orfS.scan_orf_for_pssm_tss(**kwargs)

    # end of function scan_orfs_for_pssm_tss


    def is_positioned_compatibly(self,other):
        """ Are these InwardsPointingCodingBlockGraph positioned compatibly? """
        mutualorgs = self.mutual_organisms(other)
        analyses = []
        for (pacbpkey,nodeQ,nodeS),pacbporf in self.pacbps.iteritems():
            orgQ = self.organism_by_node(nodeQ)
            orgS = self.organism_by_node(nodeS)
            if not orgQ in mutualorgs or not orgS in mutualorgs: continue
            for (_pacbpkey,_nodeQ,_nodeS),_pacbporf in other.pacbps.iteritems():
                _orgQ = other.organism_by_node(_nodeQ)
                _orgS = other.organism_by_node(_nodeS)
                if _orgQ != orgQ or _orgS != orgS: continue
                pos = pacbporf.relatively_positioned_towards(_pacbporf)
                (q1,s1) = relpos2binaryrelpos(pos['Q1'],pos['S1'])
                (q2,s2) = relpos2binaryrelpos(pos['Q2'],pos['S2'])
                # check positioning
                if q1 != s1:   analyses.append( [ False, None ] )
                elif q2 != s2: analyses.append( [ False, None ] )
                else:          analyses.append( [ True, q1[0]==1 ] )

        # Have closer look on the (True,...) tuples.
        # The orientation self->other or other->self is stored in the 2th value
        # These orientations must therefor be identical
        cnt_false = [ a for a,b in analyses ].count(False)
        cnt_ab    = [ b for a,b in analyses ].count(True)
        cnt_ba    = [ b for a,b in analyses ].count(False)
        return_list = []
        for pos in range(0,cnt_false): return_list.append(False)
        if cnt_ab > cnt_ba:
            for pos in range(0,cnt_ba): return_list.append(False)
            for pos in range(0,cnt_ab): return_list.append(True)
        elif cnt_ab < cnt_ba:
            for pos in range(0,cnt_ab): return_list.append(False)
            for pos in range(0,cnt_ba): return_list.append(True)
        else:
            # identical in size
            for pos in range(0,cnt_ab): return_list.append(False)
            for pos in range(0,cnt_ba): return_list.append(True)
        # return the return_list outcome
        return return_list

    # end of function is_positioned_compatibly


    def is_including(self,other):
        """ """
        mutualorgs = self.mutual_organisms(other)
        analyses = []
        for (pacbpkey,nodeQ,nodeS),pacbporf in self.pacbps.iteritems():
            orgQ = self.organism_by_node(nodeQ)
            orgS = self.organism_by_node(nodeS)
            if not orgQ in mutualorgs or not orgS in mutualorgs: continue
            for (_pacbpkey,_nodeQ,_nodeS),_pacbporf in other.pacbps.iteritems():
                _orgQ = other.organism_by_node(_nodeQ)
                _orgS = other.organism_by_node(_nodeS)
                if _orgQ != orgQ or _orgS != orgS: continue
                overlap = pacbporf.overlap(_pacbporf)
                if overlap == 1.0:  analyses.append(True)
                else:               analyses.append(False)

        # return analyses list
        return analyses

    # end of function is_including


# end of class InwardsPointingCodingBlockGraph

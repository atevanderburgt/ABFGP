################################################################################
### CodingBlockGraph class                                                  ####
################################################################################

"""
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


"""

# graphAbgp imports
from graphAbgp.graph_pacbpcollection import PacbpCollectionGraph
from exceptions import *

#from graphAbgp.cbg.genetree import GeneTreeOfCodingBlockFunctions
#from graphAbgp.cbg.sequenceretrieval import CodingBlockGraphSequenceRetievalFunctions
#import graphAbgp.cbg.extension as codingblock_extension
#import graphAbgp.cbg.splitting as codingblock_splitting

from codingblock_genetree import GeneTreeOfCodingBlockFunctions
from codingblock_sequenceretrieval import CodingBlockGraphSequenceRetievalFunctions
import codingblock_extension
import codingblock_splitting

# Pacb imports
import pacb

# Abgp imports
from lib_clustalw import clustalw
from lib_cexpander import cexpanderanalyses

# Python imports
from sets import Set
from copy import deepcopy

# Import Global variables and settings
from settings.codingblockgraph import *
from settings.pacbp import ALTERNATIVE_ALIGNMENT_OVERLAP_RATIO
from settings.sitealignment import ALIGNED_DONOR_MAX_TRIPLET_DISTANCE, ALIGNED_ACCEPTOR_MAX_TRIPLET_DISTANCE
from settings.splicesites import DONORPSSMOPTIONS, ACCEPTORPSSMOPTIONS

class BasicCodingBlockGraph(PacbpCollectionGraph,GeneTreeOfCodingBlockFunctions,CodingBlockGraphSequenceRetievalFunctions):
    """
    BasicCodingBlockGraph (CBG) class, inheriting from PacbpCollectionGraph class.
    """

    def __init__(self,minimal_overal_spanning_range_size=CBG_MINIMAL_OVERAL_SPANNING_RANGE_SIZE,
        alternative_alignment_overlap_ratio=ALTERNATIVE_ALIGNMENT_OVERLAP_RATIO):
        """
        Initialize a CodingBlockGraph
        """
        # Initialize from PacbpCollectionGraph
        PacbpCollectionGraph.__init__(self)

        # short name tag for printing
        self._short_name = "basicCBG"

        # 08/04/09 smallest omsr size must be 3 (means 4 AAs overall aligned)
        # this is defined in settings.codingblockgraph
        self.MINIMAL_OVERAL_SPANNING_RANGE_SIZE = minimal_overal_spanning_range_size 

        # attribute needed to set the GeneTreeGraph into
        self._GENETREE      = None

        # attributes to cache OMSR and related data structures
        self._omsr  = {}    # Overall Minimal Spanning Range
        self._msr   = {}    # Minimal Spanning Range
        self._maxsr = {}    # MAXimal Spanning Range

        # attributes to cache other CBG properties
        self._tcode         = None
        self._cexpander     = None

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
        self._GENETREE  = None
        self._omsr      = {}
        self._msr       = {}
        self._maxsr     = {}
        self._cexpander = None

    # end of function clear_cache


    def create_cache(self):
        """
        Create cached attributes of this CodingBlockGraph; saves a lot processing time
        """
        self.clear_cache()
        self._msr   = self.minimal_spanning_range()
        self._maxsr = self.maximal_spanning_range()
        self._omsr  = self.overall_minimal_spanning_range()
        self.cexpanderanalyses()
        gt = self.set_genetree()
        
    # end of function create_cache


    def omsr_identityscore(self):
        """
        Get the CBG's identityscore in the OMSR region
   
        @attention: a different measure as total_weight()
        """
        omsr = self.overall_minimal_spanning_range()
        wts = []
        for (key,nodeQ,nodeS), pacbporf in self.pacbps.iteritems():
            sta,end = min(omsr[nodeQ]), max(omsr[nodeQ])
            wt = pacbporf.identityscore_slice_by_abs_protein_query(sta,end)
            wts.append(wt)
        return float(sum(wts))/self.edge_count()
    
    # end of function omsr_identityscore


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


    def remove_alternatives_from_pacbps_dict(self):
        """
        """
        # only perform this function if graph edge_count
        # not equals self.pacbps size
        if self.edge_count() < len(self.pacbps):
            accepted, rejected = remove_alternatives_from_pacbps_dict(
                self.pacbps,
                ALTERNATIVE_ALIGNMENT_OVERLAP_RATIO=self.ALTERNATIVE_ALIGNMENT_OVERLAP_RATIO
                )
            if len(accepted) < len(self.pacbps):
                self.pacbps = accepted
        else:
            pass

    # end of function remove_alternatives_from_pacbps_dict


    def keep_only_best_alternative_from_pacbps_dict(self):
        """
        Remove all except the highest scoring PacbP(ORF) when alternatives are present
        """
        if self.edge_count() < len(self.pacbps):
            accepted, rejected = remove_alternatives_from_pacbps_dict(
                    self.pacbps,keep_only_highest_scoring=True)
            # replace self.pacbps with novel, shortened dict of pacbps
            self.pacbps = accepted 
        else:
            pass

    # end of function keep_only_best_alternative_from_pacbps_dict


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


    def _get_ranges_by_nodes(self,node1,node2,_limit_range_node1=[],_limit_range_node2=[]):
        """
        Function in use in minimal_spanning_range() and maximal_spanning_range()

        @type  node1: *
        @param node1: One Node

        @type  node2: *
        @param node2: One Node

        @type  _limit_range_node1: list
        @param _limit_range_node1: list of integers, used as allowed interval for node1

        @type  _limit_range_node2 list
        @param _limit_range_node2: list of integers, used as allowed interval for node2

		@rtype:  ( Set, Set )
		@return: tuple of 2 sets, range for Query and Sbjct

        @attention: node1 MUST be the query node, node2 MUST be the sbjct node!
        @attention: THIS FUNCTION IS ONLY 100% FUNCTIONAL FOR PacbPORF AND PacbPDNA objects!

        """
        # set known ranges by pre-limit known ranges
        thisrangeQ = Set()
        thisrangeS = Set()

        # IMPORTANT REMINDER!!! Most of the classes that use this function,
        # are (semi-)complete Graphs, with a single PacbPobj per edge.
        # In that case the iteration over ``get_pacbps_by_nodes`` will obviously
        # result in a single pacbpobj to be returned.
        # However, this function is called as well from more ``collection`` types
        # of graphs: (semi-)complete Graphs with >1 PacbPobj per edge. In the older
        # version of this function, this caused a - for long - undetected flaw in
        # the function's output. All PacbPobject of a single edge must here
        # be treated as a `single` PacbPobject. So, its range is the SUM of all
        # the individual ranges, and not the cross-section!

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
                        spos = pacbpobj._positions[ pacbpobj.alignmentposition_by_query_pos( min(_limit_range_node1) ) ]
                        startQ, startS = spos.query_pos, spos.sbjct_pos
                    if endQ-1 > max(_limit_range_node1):
                        #print "max(_limit_range_node1)", max(_limit_range_node1)
                        #print "obj.a_b_q_p()", pacbpobj.alignmentposition_by_query_pos( max(_limit_range_node1) )
                        #print node1,node2, _limit_range_node1, _limit_range_node2
                        #print self.get_nodes()
                        #print startQ, endQ
                        #print startS, endS

                        epos = pacbpobj._positions[ pacbpobj.alignmentposition_by_query_pos( max(_limit_range_node1) ) ]
                        endQ, endS = epos.query_pos+1, epos.sbjct_pos+1

                if _limit_range_node2:
                    if startS < min(_limit_range_node2):
                        spos = pacbpobj._positions[ pacbpobj.alignmentposition_by_sbjct_pos( min(_limit_range_node2) ) ]
                        startQ, startS = spos.query_pos, spos.sbjct_pos
                    if endS-1 > max(_limit_range_node2):
                        epos = pacbpobj._positions[ pacbpobj.alignmentposition_by_sbjct_pos( max(_limit_range_node2) ) ]
                        endQ, endS = epos.query_pos+1, epos.sbjct_pos+1

            else:
                # A Pacbp or PacbpDNA object.
                # Limiting ranges by preset ranges (_limit_range_node1 and node2)
                # is less easy; check the sequence for occurring gaps etc.
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


    def minimal_spanning_range(self,organism=None,node=None):
        """
        TODO description of minimal_spanning_range

        @attention: only functional if pacbp's are present in the graph
                    (perform the function `harvest_pacbps_from_crossdata`)

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


    def overall_minimal_spanning_range(self,organism=None,node=None):
        """
        TODO description of overall_minimal_spanning_range

        @attention: only functional if pacbpORFs (not Pacbps) are present in the graph
                    (perform the function `harvest_pacbps_from_crossdata`)

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


    def superinposed_left_maximal_spanning_range(self,organism=None):
        """
        @type  organism: *
        @param organism: organism identifier (or None)

        @rtype:  dictionary (or Set if organism is specified)
        @return: dictionary with nodes (keys) and spanning range Sets (values),
                 or only the spanning range Set if an organism identifier was specified
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
                smaxsrL[node] = Set( range( min(omsr[node])-(max_smaxsrL-smaxsrL[node]), min(omsr[node]) ) )

        if organism:
            return smaxsrL[ self.node_by_organism(organism) ]
        else:
            return smaxsrL

    # end of function superinposed_left_maximal_spanning_range


    def superinposed_rigth_maximal_spanning_range(self,organism=None):
        """
        @type  organism: *
        @param organism: organism identifier (or None)

        @rtype:  dictionary (or Set if organism is specified)
        @return: dictionary with nodes (keys) and spanning range Sets (values),
                 or only the spanning range Set if an organism identifier was specified
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
                smaxsrR[node] = Set( range( max(omsr[node]), max(omsr[node])+(max_smaxsrR-smaxsrR[node]) ) )

        if organism:
            return smaxsrR[ self.node_by_organism(organism) ]
        else:
            return smaxsrR

    # end of function superinposed_rigth_maximal_spanning_range


    def maximal_spanning_range(self,organism=None,node=None):
        """
        TODO description of maximal_spanning_range

        @attention: only functional if pacbp's are present in the graph
                    (perform the function `harvest_pacbps_from_crossdata`)

        @attention: only meaningfull for (nearly) Fully Connected Graphs

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


    def minimal_spanning_range_sizes(self):
        """
        Get the sizes of the minimal spanning range of all pacbp's for each node

        @attention: only functional if pacbp's are present in the graph
                    (perform the function `harvest_pacbps_from_crossdata`)

        @attention: only meaningfull for (nearly) Fully Connected Graphs

		@rtype:  dictionary
		@return: dictionary with nodes (keys) and size of spanning range Sets (values)

        """
        d = {}
        for node, r in self.minimal_spanning_range().iteritems():
            d[node] = len(r)
        return d

    # end of function minimal_spanning_range_sizes


    def overall_minimal_spanning_range_sizes(self):
        """
        Get the sizes of the overall_minimal spanning range of all pacbp's for each node

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


    def maximal_spanning_range_sizes(self):
        """
        Get the sizes of the maximal spanning range of all pacbp's for each node

        @attention: only functional if pacbp's are present in the graph
                    (perform the function `harvest_pacbps_from_crossdata`)

        @attention: only meaningfull for (nearly) Fully Connected Graphs

		@rtype:  dictionary
		@return: dictionary with nodes (keys) and size of spanning range Sets (values)
        """
        d = {}
        for node, r in self.maximal_spanning_range().iteritems():
            d[node] = len(r)
        return d

    # end of function maximal_spanning_range_sizes


    def has_minimal_spanning_range(self):
        """
        Has this (nearly Fully Connected) Graph a minimal spanning range?

        @attention: only functional if pacbp's are present in the graph
                    (perform the function `harvest_pacbps_from_crossdata`)

        @attention: only meaningfull for (nearly) Fully Connected Graphs

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


    def has_overall_minimal_spanning_range(self):
        """
        Has this (nearly Fully Connected) Graph an overall minimal spanning range?

        @attention: only functional if pacbp's are present in the graph
                    (perform the function `harvest_pacbps_from_crossdata`)

        @attention: only meaningfull for (nearly) Fully Connected Graphs

        @rtype:  Boolean
        @return: True or False
        """
        # TODO TODO temporarily set in try/except clause; this must be fixed!!
        try:
            if self.overall_minimal_spanning_range_sizes().values():
                # 8/4/09 smallest omsr size must be 3 (means 4 AAs overall aligned)
                if min(self.overall_minimal_spanning_range_sizes().values()) < self.MINIMAL_OVERAL_SPANNING_RANGE_SIZE:
                    return False
                else:
                    return True
            else:
                return False
        except:
            return False

    # end of function has_overall_minimal_spanning_range


    def has_msr(self):
        """
        @attention: alias function name of has_minimal_spanning_range
        """
        return self.has_minimal_spanning_range()

    # end of function has_msr


    def has_omsr(self):
        """
        @attention: alias function name of has_overall_minimal_spanning_range
        """
        return self.has_overall_minimal_spanning_range()

    # end of function has_omsr


    def omsrlength(self):
        """
        """
        if self.has_overall_minimal_spanning_range():
            return min(self.overall_minimal_spanning_range_sizes().values())
        else:
            return 0

    # end of function omsrlength


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


    def make_pacbps_for_missing_edges(self,use_maxsr=False,aa_extra_offset=1):
        """
        Try to make PabcpP(ORF) objects with ClustalW for absent edges in a non-complete CBG

        @attention: the CBG MUST HAVE AN overall_minimal_spanning_range !!
        """
        if self.connectivitysaturation() < 1.0:
            make_pacbps_for_missing_edges(self,
                use_maxsr=use_maxsr,
                aa_extra_offset=aa_extra_offset
                )

    # end of function make_pacbps_for_missing_edges


    def harvest_pacbps_from_pacbpcollection(self,pcg):
        """
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


    def update_edge_weights_by_minimal_spanning_range(self):
        """
        """
        nodes_done = []
        for (org,orf), omsr in self.overall_minimal_spanning_range().iteritems():
            try:
                min_omsr = min(omsr)
                max_omsr = max(omsr) + 1
            except:
                # Exception occurred in min/max of OMSR -> lof the error and than let the exception be raised
                print "UEWBMSP", self
                print "UEWBMSP", self.overall_minimal_spanning_range_sizes()
                min_omsr = min(omsr)
                max_omsr = max(omsr) + 1
            # loop over the pacbps
            for ( (a,b,c,d),(g1,o1),(g2,o2) ), pacbporf in self.pacbps.iteritems():
                if org == g1:
                    bitscore = pacbporf.bitscore_slice_by_abs_protein_query(min_omsr,max_omsr)
                elif org == g2:
                    bitscore = pacbporf.bitscore_slice_by_abs_protein_sbjct(min_omsr,max_omsr)
                else:
                    # what else? That is a case that can not happen...
                    continue
                # Now see if this node is done already.
                # If so, set the bitscore only if it is higher as the previous one
                # In practice AxB != BxA in terms of the alignment due to (small) gaps.
                # This can however cause quit big differences in bitscore calculation.
                # In stead of taking only the first of a cross combination,
                # take both and take the maximum bitscore.
                if not ( (g1,o1),(g2,o2) ) in nodes_done:
                    self.set_edge_weight( (g1,o1), (g2,o2), bitscore )
                    # append to done nodes
                    nodes_done.append( ( (g1,o1),(g2,o2) ) )
                    nodes_done.append( ( (g2,o2),(g1,o1) ) )
                else:
                    if bitscore > self.get_edge_weight( (g1,o1), (g2,o2) ):
                        self.set_edge_weight( (g1,o1), (g2,o2), bitscore )
        
    # end of function update_edge_weights_by_minimal_spanning_range


    def extend_on_left_spanningrange_difference(self,
        sprdif_min_aa_length=CBG_MIN_AA_LENGTH):
        """
        Extend CBG on left spanningrange difference

        @attention: see codingblock_extension.extend_on_spanningrange_difference for documentation
        """
        return codingblock_extension.extend_on_spanningrange_difference(self,side='left',
                sprdif_min_aa_length=sprdif_min_aa_length)

    # end of function extend_on_left_spanningrange_difference
 

    def extend_on_rigth_spanningrange_difference(self,
        sprdif_min_aa_length=CBG_MIN_AA_LENGTH):
        """
        Extend CBG on rigth spanningrange difference

        @attention: see codingblock_extension.extend_on_spanningrange_difference for documentation
        """
        return codingblock_extension.extend_on_spanningrange_difference(self,side='rigth',
                sprdif_min_aa_length=sprdif_min_aa_length)

    # end of function extend_on_rigth_spanningrange_difference


    def has_left_spanningrange_difference(self,
        sprdif_min_node_count=CBG_SPRDIF_MIN_NODE_COUNT,
        sprdif_min_aa_length=CBG_MIN_AA_LENGTH):
        """
        Has CodingBlockGraph a spanningrange difference on the left/5p/upstream side?

        @attention: see codingblock_splitting.has_left_spanningrange_difference for documentation

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
        Has CodingBlockGraph a spanningrange difference on the left/5p/upstream side?

        @attention: see codingblock_splitting.has_left_spanningrange_difference for documentation

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

        @rtype:  dict
        @return: dict with nodes as keys and coordinate tuples of deviating spanning range difference as values
        """
        return codingblock_splitting.spanningrange_difference(self,'left',correct_sprdif_all_nodes=correct_sprdif_all_nodes)

    # end of function left_spanningrange_difference


    def rigth_spanningrange_difference(self,correct_sprdif_all_nodes=True):
        """
        Returns the rigth spanningrange difference of this CodingBlockGraph

        @rtype:  dict
        @return: dict with nodes as keys and coordinate tuples of deviating spanning range difference as values
        """
        return codingblock_splitting.spanningrange_difference(self,'rigth',correct_sprdif_all_nodes=correct_sprdif_all_nodes)

    # end of function left_spanningrange_difference


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


    def distance_between_codingblocks(self,other,organism=None):
        """
        Distance in AA between two CodingBlockGraphs

        @type  organism: *
        @param organism: Organism identifier (or None)

        @rtype:  dictionary (or integer if organism or node is specified)
        @return: dictionary with organisms (keys) and AA-distance between CBGs (values),
                 or only a distance if an Organism identifier was specified

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
        @return: dictionary with organisms (keys) and AA-distance between CBGs (values),
                 or only a distance if an Organism identifier was specified

        @attention: other (CBG) is supposed to be 3p/rigth of 'self' (current CBG)
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

# end of class BasicCodingBlockGraph

################################################################################
#### Helper functions for CodingBlockGraph class                            ####
################################################################################


def remove_alternatives_from_pacbps_dict(pacbpsdict,keep_only_highest_scoring=False,ALTERNATIVE_ALIGNMENT_OVERLAP_RATIO=0.8):
    """
    Remove PacbP objects that overlap (to much) with each other

    @type:  {}
    @param: pacbpsdict dictionary with predefined structure

    @type  ALTERNATIVE_ALIGNMENT_OVERLAP_RATIO: float
    @param ALTERNATIVE_ALIGNMENT_OVERLAP_RATIO: maximal overlap ratio between PacbPs

    @attention: ALTERNATIVE_ALIGNMENT_OVERLAP_RATIO >> 0.0 and < 1.0

    @rtype  pacbpsdict: dictionary
    @return pacbpsdict: cleaned-up pacbpsdict dictionary

    @rtype  rejected: dictionary
    @return rejected: part of the pacbpsdict dictionary that has been removed
    """
    keys = pacbpsdict.keys()
    keys.sort()
    keys.reverse()
    rejected = {}
    # Capital A/B for arbitrarily first/second pacbp
    # a,b,c,d represents the 4-elem-tuple of a pacbp (bitscore,length,orf1,orf2)
    # g1,o1 represents gene/organism and orf of gene/organisms 1
    # g2,o2 represents gene/organism and orf of gene/organisms 2
    # make a full cross of all vs. all pacbps keys in the dict
    for ( (Aa,Ab,Ac,Ad),(Ag1,Ao1),(Ag2,Ao2) ) in keys:
        for ( (Ba,Bb,Bc,Bd),(Bg1,Bo1),(Bg2,Bo2) ) in keys:
            # ignore identical keys
            if (Aa,Ab,Ac,Ad) == (Ba,Bb,Bc,Bd): continue

            # only do upper triangle of square matrix
            if Aa < Ba: continue

            # ignore alignments between non-identical orf sets
            if ((Ag1,Ao1),(Ag2,Ao2)) != ((Bg1,Bo1),(Bg2,Bo2)): continue

            # yep, identical orf sets; get the PacbP objects
            pacbp1 = pacbpsdict[( (Aa,Ab,Ac,Ad),(Ag1,Ao1),(Ag2,Ao2) )]
            pacbp2 = pacbpsdict[( (Ba,Bb,Bc,Bd),(Bg1,Bo1),(Bg2,Bo2) )]

            # check if argument keep_only_highest_scoring is True
            if keep_only_highest_scoring:
                # just discard the 2th pacbp
                rejected[( (Ba,Bb,Bc,Bd),(Bg1,Bo1),(Bg2,Bo2) )] = pacbp2
                # done here; continue to next comparison
                continue

            # check relative position of PacbP objects
            relpos = pacbp1.relatively_positioned_towards(pacbp2)

            # calculate overlap based on pacbp with lowest bitscore
            overlap_q = float( relpos['Q2'][1] ) / float(sum(relpos['Q2']))
            overlap_s = float( relpos['S2'][1] ) / float(sum(relpos['S2']))

            # check if overlap is sufficiently large to remove
            if max([overlap_q,overlap_s]) >= ALTERNATIVE_ALIGNMENT_OVERLAP_RATIO:
                rejected[( (Ba,Bb,Bc,Bd),(Bg1,Bo1),(Bg2,Bo2) )] = pacbp2
                #print "DELETED:", pacbp2
                #print "KEPT   :", pacbp1
                #print "overlap:", max([overlap_q,overlap_s]), overlap_q, overlap_s

    # now separate the alternative alignments
    # from the original overlaps
    for deletethiskey,value in rejected.iteritems():
        del( pacbpsdict[deletethiskey] )  

    # and return the (new) pacbpsdict and rejected
    return pacbpsdict, rejected

# end of function remove_alternatives_from_pacbps_dict


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
        sprL = gra.left_spanningrange_difference(correct_sprdif_all_nodes=False)
        sprR = gra.rigth_spanningrange_difference(correct_sprdif_all_nodes=False)

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
            #print orf1, node1, orf1.protein_startPY, orf1.protein_endPY, min(omsr[node1]), max(omsr[node1]), aa_extra_offset
            #print orf2, node2, orf2.protein_startPY, orf2.protein_endPY, min(omsr[node2]), max(omsr[node2]), aa_extra_offset
            print "Warning: ZeroProteinSequenceLengthException", "query+", aa1start, aa1end
            print "Warning: ZeroProteinSequenceLengthException", "sbjct+", aa2start, aa2end
            continue
        elif not seq2:
            #print orf2, node2, orf2.protein_startPY, orf2.protein_endPY, min(omsr[node2]), max(omsr[node2]), aa_extra_offset
            print "Warning: ZeroProteinSequenceLengthException", "sbjct", aa2start, aa2end
            continue
        elif not seq1:
            #print orf1, node1, orf1.protein_startPY, orf1.protein_endPY, min(omsr[node1]), max(omsr[node1]), aa_extra_offset
            print "Warning: ZeroProteinSequenceLengthException", "query", aa1start, aa1end
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

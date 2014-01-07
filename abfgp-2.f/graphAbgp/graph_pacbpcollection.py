################################################################################
### PacbpCollectionGraph class                                              ####
################################################################################

# graphAbgp imports
from graph_organism import OrganismGraph
from subclass_pacbpaccessibility import PacbpAccesibilityFunctions
from subclass_sitealignment import find_fully_connected_subgraphs
from exceptions import *
import conversion

# Pacb class Import
import pacb

# BlastP SimilarityMatrix import
from lib_proteinsimilaritymatrix import ProteinSimilarityMatrix
from lib_clustalw import clustalw

# crossdata functions Import
from lib_crossdatafunctions import (
    CROSSDATA_SUBDICT,
    create_genetree_from_crossdata,
    count_pacbps_in_crossdata,
    crossdata2pacbpcollectiongraph,
    remove_nonlinear_pacbs,
    remove_inclusive_pacbps,
    remove_alternative_alignments,
    remove_repetitive_pacbps,
    recover_rejected_pacbps,
    recover_lowscoring_pacbps,
    split_accepted_pacbps_on_gapsize,
    split_lowscoring_pacbps_on_gapsize,
    )

# Python imports
from sets import Set
from copy import deepcopy

# Global variables
from settings.codingblockgraph import (
    MAXIMAL_KSMINX_PACBP_GSGOMSR_OVERLAP
    )
from settings.pacbp import (
    LINEARIZATION_STARTWITHBEST,
    LINEARIZATION_ACCEPTANCE_MARGIN,
    LINEARIZATION_WEIGTHED,
    PACBPS_SPLIT_ON_GAPSIZE,
    ALTERNATIVE_ALIGNMENT_OVERLAP_RATIO,
    REPETITIVE_ALIGNMENT_OVERLAP_RATIO,
    )

from settings.blastp import BLASTP_MATRIX_NAME

# get the blastp matrices
BLASTP_SIMILARITY_MATRIX = ProteinSimilarityMatrix(name=BLASTP_MATRIX_NAME)


class PacbpCollectionGraph(OrganismGraph,PacbpAccesibilityFunctions):
    """
    PacbpCollectionGraph (PCG) class, inheriting from OrganismGraph class.
    """
    def __init__(self,crossdata={},blastmatrix=BLASTP_SIMILARITY_MATRIX):
        """
	Initialize a PacbpCollectionGraph
        """
        # Initialize as an OrganismGraph
        OrganismGraph.__init__(self)
        # attribute to store the PacbP(ORF)s belonging to the edges
        self.pacbps = {}
        # attribute to store the (default) blastmatrix into
        self._blastmatrix = blastmatrix
	# attribute to mirror the crossdata data structure into
	self._crossdata = crossdata
	# translate crossdata dict into nodes, edges & PacbPs in the PCG
	if self._crossdata:
	    self._crossdata2pcg()
	    self.harvest_pacbps_from_crossdata(self._crossdata)

    # end of function __init__


    def deepcopy(self):
        """
        Create a QUICK deepcopy of this PCG 

        @rtype:  PacbpCollectionGraph 
        @return: ``deepcopied`` PacbpCollectionGraph

        @attention: no actual deepcopy performed on the PacbP(ORF)s !!
        @attention: NO COPYING AT ALL IS PERFORMED ON _crossdata ATTRIBUTE!
        """
        # initialize a novel PCG instance
        dpcpPCG = PacbpCollectionGraph(blastmatrix=self._blastmatrix)

        # copy the nodes
        dpcpPCG.add_nodes(self.get_nodes())

        # copy the edges
        for (node1,node2) in self.weights.keys():
            if self.has_edge(node1,node2) and not dpcpPCG.has_edge(node1,node2):
                # copy edge from input PCG 
                dpcpPCG.add_edge(node1,node2,wt=self.weights[(node1,node2)])

        # copy the pacbps
        for k,pacbp in self.pacbps.iteritems(): dpcpCBG.pacbps[k] = pacbp

        # return the ``deepcopied`` PCG 
        return dpcpPCG

    # end of function deepcopy


    def stringrepr(self,node):
        """
        Get the string representation of a (non-string) Node identifier

        @type  node: *
        @param node: One Node

        @rtype:  string 
        @return: string representation of node identifier
        """
        return "%s_%s" % ( node[0], node[1] )

    # end of function stringrepr


    def stringrepr2node(self,string):
        """
        Get the Node of a stringrepres. of a (non-string) Node identifier

        @type:  string
        @param: string representation of node identifier

        @rtype:  *
        @return: One Node
        """
        try:
            parts = string.split("_")
            return ( "_".join(parts[0:-1]), int(parts[-1]) )
        except:
            message = "stringrepr '%s', expected is '<string>_<integer>'" % (
		string )
            raise InValidStringRepresentationOfNode, message

    # end of function stringrepr2node  

    ############################################################################
    #### Overwritten functions from OrganismGraph			    ####
    ############################################################################

    def __str__(self):
        """ """
        try:    blastmatrixname = self._blastmatrix.name
        except: blastmatrixname = None
        pacbp_count = count_pacbps_in_crossdata(
                    self._crossdata,keys=['accepted_pacbps'])
        pacbp_len = None
        pacbp_type = None
        if not pacbp_count:
            pacbp_count = len(self.pacbps.keys())
            pacbp_len = sum([p.length for p in self.pacbps.values()])
            pacbp_type = Set([p.__class__.__name__ for p in self.pacbps.values()])
            if len(pacbp_type) == 1:
                pacbp_type = list(pacbp_type)[0]
            else:
                pacbp_type = "mixed"
        return "<%s N,E,P: %s,%s,%s,%s [%s] [%s] [%s]>" % (
                self.__class__.__name__,
                self.node_count(),
                self.edge_count(),
                pacbp_count,
                pacbp_len,
                pacbp_type,
                blastmatrixname,
                self.organism_set_size()
                )

    # end of function __str__


    def _organism_from_node(self,node):
        """
        Get the organism identifier from the node identifier.

        @attention: MUST BE OVERRIDDEN IN THE SUBCLASSES

        @attention: `_organism_from_node` and `organism_by_node` are aliasses

        @type  node: *
        @param node: One Node

        @rtype:  *
        @return: organism identifier (presumably a string)
        """
        return node[0]

    # end of function _organism_from_node 


    def organism_by_node(self,node):
        """
        Get the organism identifier from the node identifier.

        @attention: MUST BE OVERRIDDEN IN THE SUBCLASSES

        @attention: `_organism_from_node` and `organism_by_node` are aliasses

        @type  node: *
        @param node: One Node

        @rtype:  *
        @return: organism identifier (presumably a string)
        """
        return node[0]

    # end of function organism_by_node 


    def _update_after_changes(self):
        """
        @attention: MUST BE OVERRIDDEN IN THE SUBCLASSES
        """
        self._update_pacbps()

    # end of function _update_after_changes


    def del_edge(self,node1,node2):
        """
        Remove an edge (and its PacbPs) from the graph

        @type  node1: *
        @param node1: Node Identifier

        @type  node2: *
        @param node2: Node Identifier

        @attention: OVERWRITES graph.del_edge()
        """
        OrganismGraph.del_edge(self,node1,node2)
        self._update_pacbps()

    # end of function del_edge

    ############################################################################
    #### Conversion & PCG object integrity handling functions		    ####
    ############################################################################

    def _handle_edges_argument(self,edges):
        """
        TODO: REMOVE FROM HERE -> SHOULD BE PLACED IN sitelaignment classes
        TODO: is placed here temp. for function find_connected_node_sets() 
        """
        # if edges is not applied get by definition form the OrganismGraph
        if not edges: edges = self.organism_set_size() - 1
        # set edges to a minimum of 2
        edges = max([2,edges])
        # return edges argument
        return edges

    # end of function _handle_edges_argument


    def _crossdata2pcg(self):
	"""
	Fill this PCG with nodes & edges from PacbP(ORF)s in crossdata
	"""
	if not self._crossdata: raise NoCrossdataApplied
	self = crossdata2pacbpcollectiongraph(
	    self._crossdata,
	    increment_to_graph=self
	    )
    
    # end of function _crossdata2pcg


    def _get_crossdata(self):
	"""
	Get crossdata attribute from PCG
	
	@rtype  crossdata: dict
	@return crossdata: crossdata <dict data structure>
	"""
	if self._crossdata and self.pacbps:
	    # trust _crossdata in the PCG
	    return self._crossdata
	if self._crossdata:
	    # trust _crossdata in the PCG (no pacbps anyway...)
	    return self._crossdata
	elif self.pacbps:
	    # No filled _crossdata in PCG. Can occur when PCG is empty
	    # or when somebody messed up the object's integrity. Recreate
	    # crossdata and hope for the best...
	    self._crossdata = pacbpcollectiongraph2crossdata(self)
	    return self._crossdata
	else:
	    # Neither pacbps or crossdata applied! Can occur when PCG is empty
	    # or when somebody messed up the object's integrity. Recreate
	    # crossdata and hope for the best...
	    return pacbpcollectiongraph2crossdata(self)

    # end of function _get_crossdata


    def _sync_pcg_after_pacbp_removal(self,crossdata,
	accepted_key="",rejected_key=""):
	"""
	Synchronize edges & pacbps in the PCG after a crossdata dict operation
	in which PacbPs are removed from crossdata.

	@type  crossdata: dict
	@param crossdata: crossdata <dict data structure>

	@attention: accepted_key must be a key of CROSSDATA_SUBDICT
	@attention: rejected_key must be a key of CROSSDATA_SUBDICT

	@type  accepted_key: string
	@param accepted_key: crossdata dict key with accepted pacbps

	@type  rejected_key: string
	@param rejected_key: crossdata dict key with rejected pacbps

	@rtype:  tuple
	@return: ( number of deleted edges, number of deleted PacbP(ORF)s )
	"""
	# update the PCG by removing the edges & pacps that were lost
	deleted_edges  = 0
	deleted_pacbps = 0
        for (geneQ,geneS) in crossdata.keys():
	    accepted_edges = _crossdatasubdict2edgelist(
		    geneQ,geneS,
		    crossdata[(geneQ,geneS)][accepted_key]
		    )
	    rejected_edges = _crossdatasubdict2edgelist(
		    geneQ,geneS,
		    crossdata[(geneQ,geneS)][rejected_key]
		    )
	    for (nodeQ, nodeS) in rejected_edges:
		if (nodeQ, nodeS) not in accepted_edges:
                    if nodeQ in self.get_nodes() and nodeS in self.get_nodes():
                        # remove possibly present pacbps
                        pacbps = self.get_pacbps_by_nodes(node1=nodeQ,node2=nodeS)
                        for pacbp in pacbps:
                            self.remove_pacbp(pacbp,nodeQ,nodeS)
                            deleted_pacbps+=1
		    # remove this edge from the PCG
		    if self.has_edge(nodeQ, nodeS):
			self.del_edge(nodeQ, nodeS)
			deleted_edges+=1
                #elif nodeQ not in self.get_nodes():
                #    # nodeQ apparently deleted in a non-synced operation
                #    pass
                #elif nodeS not in self.get_nodes():
                #    # nodeS apparently deleted in a non-synced operation
                #    pass
		else:
		    # an edge present in both accepted and rejected edges!
		    # check if there are pacbps in PCG
                    if not nodeQ in self.get_nodes(): continue
                    if not nodeS in self.get_nodes(): continue
		    pacbps = self.get_pacbps_by_nodes(node1=nodeQ,node2=nodeS)
		    if not pacbps:
			# no pacbps yet in PCG, so nothing to worry about
			pass
		    elif len(pacbps) == 1:
			# just assume this is okay too; 1 pacbp left for edge
			pass
		    else:
			# Shit. Here we have to check if maybe some of the
			# PacbPs in PCG.pacbps are rejected
			barcodes = [ pacbp.barcode() for pacbp in\
			    crossdata[(geneQ,geneS)]['accepted_pacbs'].values()]
			for pos in range(len(pacbps)-1,-1,-1):
			    pacbp = pacbps[pos]
			    if pacbp.barcode() not in barcodes:
				self.remove_pacbp(pacbp,nodeQ,nodeS)
				deleted_pacbps+=1

	# run as a check-double-check the _update_pacbps function
	pacbp_count = len(self.pacbps)
	self._update_pacbps()
	deleted_edges += ( pacbp_count - len(self.pacbps) )

	# return counters of deleted edges & pacbps
	return ( deleted_edges, deleted_pacbps )

    # end of function _sync_pcg_after_pacbp_removal


    def _sync_pcg_after_pacbp_recovery(self,crossdata,
	accepted_key='accepted_pacbs'):
	"""
	Synchronize edges & pacbps in the PCG after a crossdata dict operation
	in which PacbPs are newly added to crossdata.

	@type  crossdata: dict
	@param crossdata: crossdata <dict data structure>

	@attention: accepted_key must be a key of CROSSDATA_SUBDICT

	@type  accepted_key: string
	@param accepted_key: crossdata dict key with accepted pacbps

	@rtype:  tuple
	@return: ( number of gained edges, number of gained PacbP(ORF)s )
	"""
	# update the PCG by adding new edges & pacps
	gained_edges  = 0
	gained_pacbps = 0
        for (geneQ,geneS) in crossdata.keys():
	    for key, pacbp in crossdata[(geneQ,geneS)][accepted_key].iteritems():
		(bitscore,length,orfQid,orfSid) = key
		nodeQ = (geneQ,orfQid)
		nodeS = (geneS,orfSid)
                if not nodeQ in self.nodes: self.add_node(nodeQ)
                if not nodeS in self.nodes: self.add_node(nodeS)
		if not self.has_edge(nodeQ, nodeS):
		    self.add_edge( nodeQ, nodeS, wt=bitscore )
		    gained_edges+=1
		elif self.weights[(nodeQ, nodeS)] < bitscore:
		    # very unlikely event; thit wt is higher
		    # as currently set wt. Update it here.
		    self.weights[(nodeQ, nodeS)] = bitscore
		    self.weights[(nodeS, nodeQ)] = bitscore
		else:
		    pass

	# return counters of gained edges & pacbps
	return ( gained_edges, gained_pacbps )

    # end of function _sync_pcg_after_pacbp_recovery


    def _update_pacbps(self):
        """
        Remove PacbP (or inheriting) objects that are not in the Graph
	anymore after graph operations (graph splitting, node deletions, etc.)
        """
        todelete = []
        for ((a,b,c,d),n1,n2), pacbp in self.pacbps.iteritems():
            if n1 not in self.get_nodes() or n2 not in self.get_nodes():
                todelete.append( ((a,b,c,d),n1,n2) )
            elif not self.weights.has_key((n1,n2)) or not self.weights.has_key((n2,n1)):
                todelete.append( ((a,b,c,d),n1,n2) )
            #elif not self.has_edge(n1,n2):
            #    print "..delB", n1, n2, self.weights.has_key((n1,n2)) ,self.weights.has_key((n2,n1))
            #    todelete.append( ((a,b,c,d),n1,n2) )
            else:
                pass
        for key in todelete:
            del( self.pacbps[key] )

    # end of function _update_pacbps


    def is_codingblockgraph_compatible(self):
        """
        Can this PacbpCollectionGraph be converted into a CodingBlockGraph ?

	@rtype:  Boolean
	@return: True or False
        """
        compatibility = True
        for org in self.organism_set():
            if len( self.get_organism_nodes(org) ) > 1:
                compatibility = False
                break
        return compatibility

    # end of function is_codingblockgraph_compatible


    def harvest_pacbps_from_crossdata(self,crossdata={}):
        """
        Get all edges from the graph, represented by pacbps, from the
        ``crossdata`` data structure and store them to self.pacbps
	
        @type  crossdata: dict
        @param crossdata: crossdata <data structure dict> with PacbPs
	"""
	if not crossdata:
	    # get crossdata attribute from PCG itself
	    crossdata = self._get_crossdata()

	# harvest the PacbPs from crossdata into the PCG object
        harvest_pacbps_from_crossdata(self,crossdata)

    # end of function harvest_pacbps_from_crossdata

    ############################################################################
    #### Functions that were originally functions on crossdata structure    ####
    ############################################################################

    def get_genetreegraph(self):
	"""
	Get (estimation of the) GeneTreeGraph object of the PacbPs in the graph

	@rtype:  GeneTreeGraph
	@return: GeneTreeGraph instance
	"""
	if not self._crossdata: raise NoCrossdataApplied
	return create_genetree_from_crossdata(self._crossdata)
    
    # end of function get_genetree


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


    def remove_short_pacbps(self,min_aa_length=5,verbose=False):
        """ """
        del_keys = []
        IS_ANY_PACBP_REMOVED = False
        for key,pacbp in self.pacbps.iteritems():
            if pacbp.get_unextended_length() < min_aa_length:
                del_keys.append( key )
                ################################################################
                if verbose:
                    print pacbp, key, pacbp.get_unextended_length(),
                    print "<", min_aa_length
                ################################################################

        # remove all pacbporfs/pacbpdnas listed in del_keys
        for (k,nodeQ,nodeS) in del_keys:
            status = _delete_pacbp(self,(k,nodeQ,nodeS))
            if status == True: IS_ANY_PACBP_REMOVED = True
        # return status of pacbp removal
        return IS_ANY_PACBP_REMOVED

    # end of function remove_short_pacbps


    def remove_noncoding_pacbpdnas(self,verbose=False):
        """ """
        del_keys = []
        IS_ANY_PACBP_REMOVED = False
        for key,pacbpdna in self.pacbps.iteritems():
            pacb.validators.IsPacbPDNA(pacbpdna)
            aa_coding_check = pacbpdna.is_coding()
            # if not (False,False) -> accept as potentially coding
            if aa_coding_check != (False,False): continue
            # check ratio of guided & unguided dna alignment
            dna_id_pro = pacbpdna.get_nt_identity()
            dna_id_dna = pacbpdna.get_unguided_nt_identity()
            if (dna_id_pro/dna_id_dna) < 0.90:
                del_keys.append( key )
                IS_ANY_PACBP_REMOVED = True
                ################################################################
                if verbose:
                    #pacbpdna.unextend_pacbporf()
                    print key
                    print pacbpdna, "TCODE:", pacbpdna.tcode(), 
                    print "ID%", dna_id_dna, dna_id_pro,
                    print "RATIO:", dna_id_pro/dna_id_dna
                    #pacbpdna.print_protein_and_dna()
                    #pacbpdna.extend_pacbporf_after_stops
                    #print ""
                ################################################################

        # remove all pacbporfs/pacbpdnas listed in del_keys
        for (k,nodeQ,nodeS) in del_keys:
            _delete_pacbp(self,(k,nodeQ,nodeS))
        # return status of pacbp removal
        return IS_ANY_PACBP_REMOVED

    # end of function remove_noncoding_pacbpdnas


    def remove_nonlinear_pacbs(self,
	startwithbest=LINEARIZATION_STARTWITHBEST,
	acceptance_margin=LINEARIZATION_ACCEPTANCE_MARGIN,
	is_weigthed=LINEARIZATION_WEIGTHED):
	"""
	Remove PacbP(ORF)s that break a colinear pairwise alignment

	@attention: see lib_crossdatafunctions.remove_nonlinear_pacbs for
		    more detailed documentation

	@type  startwithbest: integer
	@param startwithbest: start with this number of highest scoring Pacbs

	@type  acceptance_margin: integer
	@param acceptance_margin: offset (in AA coordinates) for overlap used
				  during the linearization proces

	@type  is_weigthed: Boolean
	@param is_weigthed: weight linearization by bitscore (default True)
	"""
	# get crossdata attribute
	local_crossdata = self._get_crossdata()

	# linearize the pacps in the dict
	local_crossdata = remove_nonlinear_pacbs(local_crossdata,
		startwithbest=startwithbest,
		acceptance_margin=acceptance_margin,
		is_weigthed=is_weigthed
		)

	# update the PCG by removing the edges & pacps that were lost
	(dedges,dpacbps) = self._sync_pcg_after_pacbp_removal(
		local_crossdata,
		accepted_key='accepted_pacbs',
		rejected_key='rejected_pacbs_nl'
		)
			
    # end of function remove_nonlinear_pacbs


    def remove_inclusive_pacbps(self,**kwargs):
        """
        Remove PacbP(ORF)s that are fully included in another PacbP(ORF)

        @attention: see lib_crossdatafunctions.remove_inclusive_pacbps() for **kwargs
        """
        # get crossdata attribute
        local_crossdata = self._get_crossdata()

        print self

        # remove alternative pacbp alignments in the dict
        local_crossdata = remove_inclusive_pacbps(local_crossdata,**kwargs)

        print self

        # update the PCG by removing the edges & pacps that were lost
        (dedges,dpacbps) = self._sync_pcg_after_pacbp_removal(
                local_crossdata,
                accepted_key='accepted_pacbs',
                rejected_key='rejected_pacbs_aa'
                )

    # end of function remove_inclusive_pacbps


    def remove_alternative_pacbps(self,
	overlap_ratio=ALTERNATIVE_ALIGNMENT_OVERLAP_RATIO):
	"""
	Remove alternative PacbP(ORF) alignments from the PCG

	@type  overlap_ratio: float
	@param overlap_ratio: minimal length overlap of lower scoring PacbP
			      in order to be removed
	"""
	# get crossdata attribute
	local_crossdata = self._get_crossdata()

	# remove alternative pacbp alignments in the dict
	local_crossdata = remove_alternative_alignments(local_crossdata,
		overlap_ratio=overlap_ratio)

	# update the PCG by removing the edges & pacps that were lost
	(dedges,dpacbps) = self._sync_pcg_after_pacbp_removal(
		local_crossdata,
		accepted_key='accepted_pacbs',
		rejected_key='rejected_pacbs_aa'
		)
    
    # end of function remove_alternative_pacbps


    def remove_repetitive_pacbps(self,
        overlap_ratio=REPETITIVE_ALIGNMENT_OVERLAP_RATIO):
        """
        Remove repetitive (overlapping) PacbP(ORF) alignments from the PCG

        @type  overlap_ratio: float
        @param overlap_ratio: minimal length overlap of lower scoring PacbP
                              in order to be removed
        """
        # get crossdata attribute
        local_crossdata = self._get_crossdata()

        # remove alternative pacbp alignments in the dict
        local_crossdata = remove_repetitive_pacbps(local_crossdata,
                overlap_ratio=overlap_ratio)

        # update the PCG by removing the edges & pacps that were lost
        (dedges,dpacbps) = self._sync_pcg_after_pacbp_removal(
                local_crossdata,
                accepted_key='accepted_pacbs',
                rejected_key='rejected_pacbs_aa'
                )

    # end of function remove_repetitive_pacbps



    def recover_lowscoring_pacbps(self,crossdata={},
        genetreegraph=None,blastmatrix=None):
	"""
	Recover small, low-scoring PacbPs based on GTG characteristics

	@type  crossdata: dict
	@param crossdata: crossdata <dict data structure>

	@type  genetreegraph: GeneTreeGraph
	@param genetreegraph: GeneTreeGraph instance

        @type  blastmatrix: BlastMatrix
        @param blastmatrix: BlastMatrix instance

	@attention: crossdata variable overrrides PCG._get_crossdata()
	@attention: genetreegraph variable overrrides PCG.get_genetreegraph()
        @attention: blastmatrix variable overrrides PCG._blastmatrix attribute
	"""
	if not crossdata:
	    # get crossdata attribute from PCG itself
	    crossdata = self._get_crossdata()

	if not genetreegraph:
	    # get genetreegraph attribute from PCG itself
	    genetreegraph = self.get_genetreegraph()

        if not blastmatrix:
            # get blastmatrix object from PCG itself
            blastmatrix = self._blastmatrix

	# recover lowscoring PacbPs
	crossdata = recover_lowscoring_pacbps(
		crossdata,genetreegraph,blastmatrix)

	# update the PCG by adding new edges & pacps
	(aedges,apacbps) = self._sync_pcg_after_pacbp_recovery(crossdata)
    
    # end of function recover_lowscoring_pacbps


    def recover_rejected_pacbps(self,crossdata={},omit_pcg_sync=False):
        """
        Recover absent edges (pacbps) between nodes in the graph that
        are present in the ``crossdata`` data structure.

        @type  crossdata: dict
        @param crossdata: crossdata <data structure dict> with PacbPs

        @type  omit_pcg_sync: Boolean
        @param omit_pcg_sync: default False,set to True in case PCG ~= CBG
        """
	if not crossdata:
	    # get crossdata attribute from PCG itself
	    crossdata = self._get_crossdata()

        if self.connectivitysaturation() < 1.0:
            self, crossdata = recover_rejected_pacbps(self,crossdata)

	# update the PCG by adding new edges & pacps
        if not omit_pcg_sync:
	    (aedges,apacbps) = self._sync_pcg_after_pacbp_recovery(crossdata)

    # end of function recover_rejected_pacbps


    def split_pacbps_on_gapsize(self,gapsize=None):
	"""
	Split PacbPs in the PCG on (GTG obtained) gapsize

	@type  gapsize: integer or None
	@param gapsize: AA length gapsize to split PacbPs on
	
	@attention: when gapsize=None a sensible gapsize will be estimated
		    based on GTG characteristics
	"""
	if gapsize and gapsize >= 1:
	    pass
	if not self._crossdata:
	    gapsize = int(float(PACBPS_SPLIT_ON_GAPSIZE)/0.50)
	else:
	    gtg = self.get_genetreegraph()
	    if gtg.node_count():
		gapsize = int(float(PACBPS_SPLIT_ON_GAPSIZE)/gtg.identity())
	    else:
		# hmm...GTG obtained but it is empty! That means that after
		# crossblasting no K(s) graph can be construted. Either not
		# homologous proteins or a very low similar group op proteins.
		# Therefor, assume a very low id%
		gapsize = int(float(PACBPS_SPLIT_ON_GAPSIZE)/0.50)

	# split Pacbps in crossdata in calculated gapsize
	self._crossdata = split_accepted_pacbps_on_gapsize(
		self._crossdata,
		gapsize=gapsize
		)
    
	# update edges and PacbPs in the PCG itself
	self._sync_pcg_after_pacbp_recovery(self._crossdata)
	self._update_pacbps()
	
    # end of function split_pacbps_on_gapsize


    def split_lowscoring_pacbps_on_gapsize(self,gapsize=None):
	"""
	Split lowscoring PacbPs (in crossdata) and store them in the PCG
	based on (GTG obtained) gapsize

	@type  gapsize: integer or None
	@param gapsize: AA length gapsize to lowscoring split PacbPs on
	
	@attention: when gapsize=None a sensible gapsize will be estimated
		    based on GTG characteristics
	"""
	if gapsize and gapsize >= 1:
	    pass
	if not self._crossdata:
	    gapsize = int(float(PACBPS_SPLIT_ON_GAPSIZE)/0.50)
	else:
	    gtg = self.get_genetreegraph()
	    if gtg.node_count():
		gapsize = int(float(PACBPS_SPLIT_ON_GAPSIZE)/gtg.identity())
	    else:
		# hmm...GTG obtained but it is empty! That means that after
		# crossblasting no K(s) graph can be construted. Either not
		# homologous proteins or a very low similar group op proteins.
		# Therefor, assume a very low id%
		gapsize = int(float(PACBPS_SPLIT_ON_GAPSIZE)/0.50)

	# split Pacbps in crossdata in calculated gapsize
	self._crossdata = split_lowscoring_pacbps_on_gapsize(
		self._crossdata,
		gapsize=gapsize
		)
    
	# update edges and PacbPs in the PCG itself
	self._sync_pcg_after_pacbp_recovery(self._crossdata)
	self._update_pacbps()
	
    # end of function split_lowscoring_pacbps_on_gapsize


    ############################################################################
    #### Functions for subgraph retrieval				    ####
    ############################################################################

    def find_connected_node_sets(self,edges=None,max_missing_edges=None):
        """
        Find (potentially) connected graphs in the input graph.

        @type  edges: number
        @param edges: number of outgoing edges of a node in a FCG

        @type  max_missing_edges: number
        @param max_missing_edges: number of missing edges to allow in a
				  (nearly) K(s) PacbpCollectionGraph

        @rtype:  list
        @return: list with (nearly) K(s) PacbpCollectionGraphs
		 of the requested properties
        """
        # handle edges argument (error check etc.)
        edges = self._handle_edges_argument(edges)

        # list variable to store node sets
        nodesets = []

        # loop over the nodes and find nodegroups
        for node in self.get_ordered_nodes():
            nodes2org = { self.organism_by_node(node): [ node ] }
            for nodeA,nodeB in self.weights.keys():
                if nodeA == node:
                    orgofnode = self.organism_by_node(nodeB)
                    if nodes2org.has_key(orgofnode):
                        nodes2org[orgofnode].append( nodeB )
                    else:
                        nodes2org[orgofnode] = [ nodeB ]
            if len(nodes2org) >= edges+1:
                paths = nodes2org.values()
                paths.sort()
                combis = pacb.recombination.cross(paths)
                for combi in combis:
                    if combi not in nodesets:
                        nodesets.append(combi)

        # TODO TODO TODO TODO
        # TODO make this part more generic -> no HARD-SET PACBPCOL is returned
        # this function is as well usefull in other graph types!
        # do this by making a _get_empty_graph() function that is
        # overwritten in each graph class. Then, a novel graph of each
        # type can be created by a _nodes2XXXgraph() function, that is as
        # well overwritten in each graph class.

        # recombine into sub-PacbpCollectionGraphs
        pacbpcollections = []
        for nodes in nodesets:
            pacbpcol = self._nodes2pacbpsubcollection(nodes)

            ## check edge argument. If pacbpcol is larger, remove node(s)
            #while pacbpcol.node_count() > edges+1:
            #   weakestnode = pacbpcol.weakest_connected_node()
            #   pacbpcol.del_node(weakestnode)
            #if (max_missing_edges or max_missing_edges==0) and\
            #pacbpcol.max_edge_count() - pacbpcol.edge_count() >\
            #max_missing_edges:
            #    pass
            #else:
            #    pacbpcollections.append(pacbpcol)

            if pacbpcol.node_count() > edges+1:
                # obtain connectivity data of the nodes
                connectivity = []
                for node in pacbpcol.get_nodes():
                    connectivity.append(
                        ( self.get_node_weighted_connectivity(node), node ) )
                connectivity.sort()
                # start removing weakest node(s), with using
                # each node once as an anchor point
                node_combinations_seen = []
                for anchornode in pacbpcol.get_nodes():
                    dpcppacbpcol = pacbpcol.deepcopy()
                    pos = 0
                    for pos in range(0,dpcppacbpcol.node_count()):
                        weakestnode = connectivity[pos][1]
                        pos += 1
                        if weakestnode == anchornode:
                            # prevent anchornode from removal!
                            continue
                        if weakestnode in dpcppacbpcol.nodes[anchornode] and\
                        len(dpcppacbpcol.nodes[anchornode]) < edges:
                            # prevent one of the few edges to achornode
                            # to be deleted
                            continue
                        # if here, remove this node
                        dpcppacbpcol.del_node(weakestnode)
                        # check if enough nodes are removed
                        if dpcppacbpcol.node_count() == edges+1:
                            break

                    # check if enough nodes were deletable from dpcppacbpcol
                    if dpcppacbpcol.node_count() > edges+1:
                        # graph topology did not allow node deletion
                        pass
                    elif (max_missing_edges or max_missing_edges==0) and\
                    dpcppacbpcol.max_edge_count() - dpcppacbpcol.edge_count() >\
                    max_missing_edges:
                        # sub-combination has to little edges
                        pass
                    else:
                        # store dpcppacbpcol as a connected node set
                        if dpcppacbpcol.get_ordered_nodes() not in\
                        node_combinations_seen:
                            pacbpcollections.append(dpcppacbpcol)
                            node_combinations_seen.append(
                                    dpcppacbpcol.get_ordered_nodes() )
                        else:
                            # this sub-combination is already seen
                            pass

            else:
                if (max_missing_edges or max_missing_edges==0) and\
                pacbpcol.max_edge_count() - pacbpcol.edge_count() >\
                max_missing_edges:
                    pass
                else:
                    pacbpcollections.append(pacbpcol)


        # return list of pacbpcollections
        return pacbpcollections

    # end of function find_connected_node_sets


    def _nodes2pacbpsubcollection(self,nodes):
        """
	@type  nodes: list
	@param nodes: list with nodes

	@rtype  pacbpcol: PacbpCollectionGraph
	@return pacbpcol: PacbpCollectionGraph
        """
        pacbpcol = PacbpCollectionGraph()
        pacbpcol.add_nodes(nodes)
        for nodeA,nodeB in pacb.recombination.pairwise(nodes):
            if self.has_edge(nodeA,nodeB):
                pacbpcol.add_edge(nodeA,nodeB,wt=self.weights[(nodeA,nodeB)])
        return pacbpcol

    # end of function _nodes2pacbpsubcollection


    def obtain_ksminxgraphs(self,input,crossdata,gsg=None,
        max_ksminx_pacbps_gsgomsr_overlap=MAXIMAL_KSMINX_PACBP_GSGOMSR_OVERLAP,
        method='strict',
        verbose=False):
        """
        Get K(s-x) CBGs from a PacbpCollectionGraph

        @type  input: dict
        @param input: input <data structure dict> that contains lists of orfs

        @type  crossdata: dict
        @param crossdata: crossdata <data structure dict> with PacbPs

        @type  gsg: GenestructureOfCodingBlockGraphs
        @param gsg: GenestructureOfCodingBlockGraphs or None

	@type  max_ksminx_pacbps_gsgomsr_overlap: integer
	@param max_ksminx_pacbps_gsgomsr_overlap: length (in AAs) of alowed
		    overlap between PacbPs in K(s-x) CBGs with the applied gsg

        @type  method: string 
        @param method: 'loose' or 'strict'

        @type  verbose: Boolean
        @param verbose: print messages to STDOUT when True, default False

        @attention: requires global var MAXIMAL_KSMINX_PACBP_GSGOMSR_OVERLAP
	@attention: K(s-x) CBGs that overlap (or coincide) the gsg are discarded

	@rtype:  list
	@return: list with compatible K(s-x) CBGs (with PacbPORFs)
        """
        potential_ksminx_cbgraphs = []
        ksminx_cbgraphs = []
        if gsg:
            # calculate OMSR of the GSG
            gsgOMSR = gsg.overall_minimal_spanning_range()
            EXACT_SG_NODE_COUNT = gsg.EXACT_SG_NODE_COUNT
            # remove edges that can never contribute to GSG enrichments
            self._remove_gsg_edges_upon_obtain_ksminxcbgs(gsg,
                    crossdata=crossdata,verbose=verbose)
            MAX_MISSING_EDGES = 1
        else:
            EXACT_SG_NODE_COUNT = len(input)
            MAX_MISSING_EDGES = 1

        # check if max_ksminx_pacbps_gsgomsr_overlap >= 0 / exists
        if max_ksminx_pacbps_gsgomsr_overlap or\
        max_ksminx_pacbps_gsgomsr_overlap == 0:
            check_pacbp_gsgomsr_overlap = True
        else:
            check_pacbp_gsgomsr_overlap = False

        ########################################################################
        if verbose:
           for org in self.organism_set():
               print "\tPCG nodes:", org, "\t",
               print [ node[1] for node in self.get_organism_nodes(org) ]
        ########################################################################

        # obtain potential_ksminx_cbgraphs from find_connected_node_sets funct.
        if method == 'loose':
            # example K(s)CBG has 5 nodes. Then, here is searched for:
            # 4.3, 3.2, 2.1 edges/missing edges as K(s-x) CBGs
            # In fact 4.3 is not a K(s-x) CBG, but a K(s) CBG that misses
            # many edges
            fcn_settings = [ (i,i-1) for i in range(EXACT_SG_NODE_COUNT-1,1,-1) ]
        else:
            # example K(s)CBG has 5 nodes. Then, here is searched for:
            # 3.1, 2.0 edges/missing edges as K(s-x) CBGs
            fcn_settings = [ (i,i-2) for i in range(EXACT_SG_NODE_COUNT-2,1,-1) ]

        for (edges,max_missing_edges) in fcn_settings:
            potential_ksminx_cbgraphs.extend(
                    self.find_connected_node_sets(
                            edges=edges,
                            max_missing_edges=max_missing_edges
                            ) )
            ####################################################################
            if verbose:
                print "K(s-x) find_connected_node_sets():",
                print len( potential_ksminx_cbgraphs ),
                print "edges:", edges, 
		print "max_missing_edges:", max_missing_edges
            ####################################################################

        if EXACT_SG_NODE_COUNT == 3:
            # only 3 organisms/genes .... 
            potential_ksminx_cbgraphs = self.find_fully_connected_subgraphs(
                    edges=2,max_missing_edges=1)
            ####################################################################
            if verbose:
                print "K(s-x) find_connected_node_sets():",
                print len( potential_ksminx_cbgraphs ),
                print "edges:", 2,
                print "max_missing_edges:", 1
            ####################################################################


        for sg in potential_ksminx_cbgraphs:
            # check if the node_set is not equal to a CBG
	    # that is already in the GSG
            if gsg:
                identical_node_set = False
                for cbg in gsg:
                    if not sg.node_set().difference(cbg.node_set()):
                        identical_node_set = True
                        break
                if identical_node_set: continue
            # check if this node_set is not equal to
	    # a partialCBG that is already accepted
            identical_node_set = False
            for cbg in ksminx_cbgraphs:
                if not sg.node_set().difference(cbg.node_set()):
                    identical_node_set = True
                    break
            if identical_node_set: continue

            # gradually build out sg(CBG) to a fully functional CBG
            sg.recover_rejected_pacbps(crossdata,omit_pcg_sync=True)
            sg.harvest_pacbps_from_crossdata(crossdata)
            ksminxcbg = conversion.PacbpCollectionGraph2CodingBlockGraph(sg)
            ksminxcbg.pacbps2pacbporfs(input)

            if not ksminxcbg.has_overall_minimal_spanning_range():
                continue

            ####################################################################
            if verbose:
                try:    print "cand: %s" % ( str(ksminxcbg) ),
                except: print "cand: %s" % ( ksminxcbg.get_ordered_nodes() ),
                print ksminxcbg.node_count(), ksminxcbg.edge_count(),
                print len(ksminxcbg.pacbps)
            ####################################################################

            # confirm pacbps in gsgOMSR. If not placeable -> ignore/continue
            if gsg and check_pacbp_gsgomsr_overlap:
                is_placeable = True
                cbgomsr = ksminxcbg.overall_minimal_spanning_range()
                for node,omsr in cbgomsr.iteritems():
                    org = ksminxcbg.organism_by_node(node)
                    iso = gsgOMSR[org].intersection(omsr)
                    if len(iso) > MAXIMAL_KSMINX_PACBP_GSGOMSR_OVERLAP:
                        is_placeable=False
                        break
                # not placeable -> continue (and save time!)
                if not is_placeable: continue

            # now build out to a fully functional & compatible (ksminx)CBG
            ksminxcbg.keep_only_best_alternative_from_pacbps_dict()
            ksminxcbg.extend_pacbporfs(input)
            ksminxcbg.make_pacbps_for_missing_edges()
            if ksminxcbg.connectivitysaturation() != 1.0: continue
            if not ksminxcbg.has_overall_minimal_spanning_range(): continue

            # now, we have CBGs with OMSR. Do further checks
            ksminxcbg.update_edge_weights_by_minimal_spanning_range()
            ksminxcbg.cexpanderanalyses()
            if ksminxcbg._cexpander.binarystring.count("1") == 0: continue

            # if this point is reached, this ksminxcbg is accepted as a
            # K(s-x) CBG graph which is a potential enrichment for the GSG
            ksminx_cbgraphs.append( ksminxcbg )
            ####################################################################
            if verbose: print "K(s-x):", ksminxcbg 
            ####################################################################


        ########################################################################
        if verbose:
            print "potential_ksminx_cbgraphs:", len(potential_ksminx_cbgraphs)
            print "ksminx_cbgraphs:", len(ksminx_cbgraphs)
            for org in self.organism_set():
               print "\tPCG nodes:", org, "\t",
               print [ node[1] for node in self.get_organism_nodes(org) ]
        ########################################################################

        # return the list if K(s-x) graphs
        return ksminx_cbgraphs

    # end of function obtain_ksminxgraphs


    def _remove_gsg_edges_upon_obtain_ksminxcbgs(self,GSG,crossdata={},verbose=False):
        """
        """
        if not crossdata: crossdata = self._get_crossdata()
        if not count_pacbps_in_crossdata(crossdata):
            # no pacbps listed in crossdata. If we continue here,
            # ALL PCG edges are deleted. So, break out
            return 0 

        # use the gsgOMSR to identify non-informative edges
        gsgOMSR = GSG.overall_minimal_spanning_range()
        ############################################################
        if verbose:
            print "PCG", self
            print "GSG", GSG
            for org in gsgOMSR.keys():
                print "gsgOMSR", org, "\t", len(gsgOMSR[org])
        ############################################################
        
        # remove edges from the PCG that are not an enrichment
        # for the gsg once in a K(s-x) CBG
        cnt = 0
        if gsgOMSR:
            accepted_edges = []
            for (geneQ,geneS) in crossdata.keys():
                # try to perform the splits on the dict `accepted_pacbps`
                for pacbpkey, pacbp in\
                crossdata[(geneQ,geneS)]['accepted_pacbs'].iteritems():
                    nodeQ = ( geneQ, pacbpkey[2] )
                    nodeS = ( geneS, pacbpkey[3] )
                    coordSetQ = Set(range(pacbp.query_start,pacbp.query_end))
                    coordSetS = Set(range(pacbp.sbjct_start,pacbp.sbjct_end))
                    if not coordSetQ.difference( gsgOMSR[ geneQ ] ):
                        if verbose: print "NO new K(s-x) edge:", nodeQ, nodeS
                        pass
                    elif not coordSetS.difference( gsgOMSR[ geneS ] ):
                        if verbose: print "NO new K(s-x) edge:", nodeQ, nodeS
                        pass
                    else:
                        # this edge is fine; potentially an enrichment
                        accepted_edges.append( ( nodeQ, nodeS) )

            # remove all edges in PCG that are NOT in accepted_edges
            current_edges = self.weights.keys()
            for ( nodeQ, nodeS ) in current_edges:
                edge = [ nodeQ, nodeS ]
                edge.sort()
                if edge != [ nodeQ, nodeS ]: continue
                if ( nodeQ, nodeS ) not in accepted_edges:
                    self.del_edge( nodeQ, nodeS )
                    cnt+=1
                    if verbose: print "edge deleted:", nodeQ, nodeS

        # return number of edges that are deleted
        return cnt

    # end of function _remove_gsg_edges_upon_obtain_ksminxcbgs


    def find_fully_connected_subgraphs(self,edges=None,max_missing_edges=2):
        """
        Find all Fully Connected Graphs (FCG) in the input graph.

        @type  edges: number
        @param edges: number of outgoing edges of a node in a FCG

        @type  max_missing_edges: number
        @param max_missing_edges: number of missing edges to allow in a (nearly)-FCG

        @rtype:  list
        @return: list with CodingBlockGraphs of the requested properties
        """
        # if edges is not applied get by definition from the OrganismGraph
        if not edges: edges = self.organism_set_size() - 1

        retlist = []
        for sg in find_fully_connected_subgraphs(self,edges=edges,max_missing_edges=max_missing_edges):
            if sg.is_codingblockgraph_compatible():
                retlist.append( conversion.PacbpCollectionGraph2CodingBlockGraph(sg) )
            else:
                retlist.append( sg )
        # return the splitted graphs                
        return retlist

    # end of function find_fully_connected_subgraphs


    def recombine_into_codingblockgraphs(self,edges=None,max_missing_edges=2):
        """
        Create all possible CBGs by organism node recombination 

        @type  edges: number
        @param edges: number of outgoing edges of a node in a FCG

        @rtype:  list
        @return: list with CodingBlockGraphs of the requested properties
        """
        from codingblock_splitting import cross

        # if edges is not applied get by definition from the OrganismGraph
        if not edges: edges = self.organism_set_size() - 1

        retlist = []
        # make a cross of all the pacbp positions in the lists of alternatives 
        allcombis = cross([ self.get_organism_nodes(org) for org in self.organism_set() ])
        for combi in allcombis:
            sg = deepcopy(self)
            sg.nodes  = {}
            sg.pacbps = {} 
            sg.add_nodes(combi)
            edges = sg.weights.keys()
            for (k,n1,n2) in self.pacbps.keys():
                (insg1,insg2) =  n1 in sg.get_nodes(), n2 in sg.get_nodes()
                if (insg1,insg2) == (True,True):
                    sg.pacbps[(k,n1,n2)] = self.pacbps[(k,n1,n2)]
                elif True in (insg1,insg2):
                    try:
                        del( sg.weights[(n1,n2)])
                        del( sg.weights[(n2,n1)])
                    except KeyError:
                        continue
            if sg.is_codingblockgraph_compatible():
                cbg = conversion.PacbpCollectionGraph2CodingBlockGraph(sg)
                retlist.append(cbg)
                print cbg, len(cbg.weights), len(cbg.pacbps)

        # and return all the recombinations!
        return order_graphlist_by_total_weight(retlist)

    # end of function recombine_into_codingblockgraphs


    def find_missing_node(self,ORGANISMS,crossdata,minimal_new_connections=2):
        """
        """
        find_missing_node(self,ORGANISMS,crossdata,
		minimal_new_connections=minimal_new_connections)

    # end of function find_missing_node


    def delete_weakest_organism_node(self,organism):
        """
        Delete the weakest connected node of an organism

        @type  organism: *
        @param organism: organism identifier (or None)
        """
        wts = []
        for node in self.get_nodes():
            if self._organism_from_node(node) == organism:
                wts.append((self.get_node_weighted_connectivity(node), node))

        # Quit this function if there are already less than 2 nodes
        # of this organism (0 nodes -> not present, 1 node -> already only one
        if len(wts) < 2: return

        # sort the weights; first one is poorest scoring
        wts.sort()
        del_node = wts[0][1]
        # delete node from the graph
        self.del_node(del_node)
        pacbps_keys = self.pacbps.keys()
        # and cleanup the pacbps dictionairy
        for ( (a,b,c,d),n1,n2 ) in pacbps_keys:
            if del_node in (n1,n2):
                del( self.pacbps[( (a,b,c,d),n1,n2 )] )
            
    # end of function delete_weakest_organism_node


    def remove_lowerasexpected_ntidentity_pacbps(self,GENECOMBIS,input,
        ignore_annotated_exons=True,ratio=0.725,verbose=True):
        """
        Code for pairwise ABFGP mode, january 2011

        @attention: requires PacbPORFs, not PacbPs
        """
        pacbporfs_removed_cnt = 0
        for orgA,orgB in GENECOMBIS:
            if not orgA in self.organism_set(): continue
            if not orgB in self.organism_set(): continue
            # order pacbplist, highest bitscore first
            pacbpsubdict = {}
            for (k,(orgQid,orfQid),(orgSid,orfSid)),pacbporf in self.pacbps.iteritems():
                if orgQid == orgA and orgSid == orgB:
                    pacbpsubdict[(k,(orgQid,orfQid),(orgSid,orfSid))] = pacbporf

            # order pacbporfs on bitscore
            pacbporfs = pacb.ordering.order_list_by_attribute(
                    pacbpsubdict.values(),order_by='bits',reversed=True)

            # continue if no pacbporfs found
            if not pacbporfs: continue

            # get nt-identity of the highest scoring pacbporf alignment
            best_pacbporf_ntidentity = pacbporfs[0].get_nt_identity()

            for key,pacbporf in pacbpsubdict.iteritems():
                if pacbporf.get_nt_identity() / best_pacbporf_ntidentity < ratio:
                    # check if ignore_annotated_exons is set to True
                    if ignore_annotated_exons:
                        if pacbporf.orfQ.id in input[orgA]['orfid-genestructure'] and\
                        pacbporf.orfS.id in input[orgB]['orfid-genestructure']:
                            # do NOT delete this pacbporf
                            continue

                    # remove this PacbPORF from the PCG
                    _delete_pacbp(self,key)
                    pacbporfs_removed_cnt+=1
                    ############################################################
                    if verbose:
                        print "POOR-NT-IDENTITY: %s %1.2f << %1.2f %s" % (
                            orgB,pacbporf.get_nt_identity(),
                            best_pacbporf_ntidentity,
                            pacbporf)
                    ############################################################
                    

        # return number of removed pacbporfs
        return pacbporfs_removed_cnt

    # end of function remove_lowerasexpected_ntidentity_pacbps


    def merge_high_gap_ratio_pacbporfs(self,target,high_gap_ratio_score=0.40,
        enforce_nodes_to_be_merged = [], verbose=False):
        """ """
        # find target nodes that have high_gap_ratio
        high_gap_ratio_nodes = Set()
        for (a,nQ,nS),pacbporf in self.pacbps.iteritems():
            if pacbporf.gap_ratio_score() >= high_gap_ratio_score:
                high_gap_ratio_nodes.add(nQ)

        # check if nodes have been enforced to be merged
        if enforce_nodes_to_be_merged:
            for node in enforce_nodes_to_be_merged:
                if node in self.get_nodes():
                    high_gap_ratio_nodes.add(node)
        
        # return status False if no merges are allowed
        if not high_gap_ratio_nodes: return False

        ########################################################################
        if verbose: print "HIGH GAP RATIO SCORE nodes::", high_gap_ratio_nodes
        ########################################################################

        is_any_pacbporf_merged = False
        for informant in self.organism_set():
            if informant == target: continue
            pacbporfs = self.get_pacbps_by_organisms(target,informant)
            for targetNode in high_gap_ratio_nodes:
                involved_informant_nodes = []
                for pacbporf in pacbporfs:
                    if (target, pacbporf.orfQ.id) == targetNode:
                        involved_informant_nodes.append( ( informant, pacbporf.orfS.id ) )
                # check how many times an informant node is observed
                for informantNode in Set(involved_informant_nodes):
                    if involved_informant_nodes.count(informantNode) == 1:
                        # we are searching here for multiple present node IDs
                        continue
    
                    # get subset of pacbporfs that must be merged
                    pacbporf_subset = self.get_pacbps_by_nodes(targetNode,informantNode)
    
                    # obtain coords that must be merged
                    pfsQsta = min([ min(pf.alignment_protein_range_query()) for pf in pacbporf_subset ])
                    pfsQend = max([ max(pf.alignment_protein_range_query()) for pf in pacbporf_subset ]) + 1
                    pfsSsta = min([ min(pf.alignment_protein_range_sbjct()) for pf in pacbporf_subset ])
                    pfsSend = max([ max(pf.alignment_protein_range_sbjct()) for pf in pacbporf_subset ]) + 1
                    orfQ = pacbporf_subset[0].orfQ
                    orfS = pacbporf_subset[0].orfS
    
                    # obtain coords relative to the Orf objects
                    orfQsta = pfsQsta-orfQ.protein_startPY
                    orfQend = orfQsta + (pfsQend-pfsQsta)
                    orfSsta = pfsSsta-orfS.protein_startPY
                    orfSend = orfSsta + (pfsSend-pfsSsta)
    
                    # get sequence parts of the Orf objects
                    seqpQ = orfQ.protein_sequence[orfQsta:orfQend]
                    seqpS = orfS.protein_sequence[orfSsta:orfSend]

                    # apply an anchor on the to-be-aligned sequence to
                    # enforce its (potentially BLAST-obtained) tails to
                    # be re-aligned as they were
                    anchor_size = 5
                    anchor = "W"*anchor_size
                    seqs = { target: anchor + seqpQ + anchor,
                             informant: anchor + seqpS + anchor }
    
                    # make pacbp from via clustalw alignment
                    (alignedseqs,alignment) = clustalw( seqs=seqs)
                    # remove the anchors from the ClustalW output
                    alignedseqs[target] = alignedseqs[target][anchor_size:-anchor_size]
                    alignedseqs[informant] = alignedseqs[informant][anchor_size:-anchor_size]
                    alignment = alignment[anchor_size:-anchor_size]

                    # define clustalw input for PacbP formation
                    clw2pfinput = ( alignedseqs[target],
                                    alignment,
                                    alignedseqs[informant] )

                    # make pacbp
                    pacbp = pacb.conversion.pacbp_from_clustalw(
                                    alignment=clw2pfinput,
                                    coords= (pfsQsta,pfsQend,pfsSsta,pfsSend) )
    
                    # remove the anchor from the pacbp
                    for iter in range(0,anchor_size):
                        pacbp._strip_one_left_position()
                        pacbp._strip_one_rigth_position()
                    pacbp.strip_unmatched_ends()
    
                    # make extended pacbporf of the pacbp
                    newpacbporf = pacb.conversion.pacbp2pacbporf(pacbp,orfQ,orfS)
                    newpacbporf.extend_pacbporf_after_stops()
                    new_key = newpacbporf.construct_unique_key(targetNode,informantNode)
  
                    ############################################################
                    if verbose: print "MERGED::", informant, newpacbporf 
                    ############################################################
 
                    for pf in pacbporfs:
                        if (newpacbporf.orfQ.id,newpacbporf.orfS.id) ==\
                        (pf.orfQ.id,pf.orfS.id):
                            # delete this pacbp
                            self.remove_pacbp(pf,targetNode,informantNode)
                            ####################################################
                            if verbose: print "merging::", informant, pf
                            ####################################################
 
                    # place the merged pacbporf into the PCG
                    self.pacbps[(new_key,targetNode,informantNode)] = newpacbporf
                    is_any_pacbporf_merged = True

        # return status weather or not pacbporfs have been merged
        return is_any_pacbporf_merged

    # end of function merge_high_gap_ratio_pacbporfs


    def correct_overlaps(self,GENECOMBIS,verbose=False):
        """
        Code for pairwise ABFGP mode, december 2010
        """
        ############################################################################
        ### pacbp overlap correction
        ############################################################################
        self.pacbporfs2pacbps()
        zerosized_pacbp_encountered = False
        for orgA,orgB in GENECOMBIS:
            if not orgA in self.organism_set(): continue
            if not orgB in self.organism_set(): continue
            pacbporfs = pacb.ordering.order_pacbp_list(self.get_pacbps_by_organisms(orgA,orgB))
            for pos in range(1,len(pacbporfs)):
                prevPACBP = pacbporfs[pos-1]
                nextPACBP = pacbporfs[pos]
                prevPACBP,nextPACBP,status = pacb.overlap.correct_overlap_for_sbjct(
                            prevPACBP,nextPACBP,verbose=verbose)
                prevPACBP,nextPACBP,status = pacb.overlap.correct_overlap_for_query(
                            prevPACBP,nextPACBP,verbose=verbose)
                # check if the PacbPORF can still be stripped
                for pb in [ prevPACBP,nextPACBP ]:
                    if pb.length==0: continue
                    if '-' in [ pb.query[0], pb.sbjct[0], pb.query[-1], pb.sbjct[-1] ]:
                        pb.strip_unmatched_ends()
                if prevPACBP.length == 0: zerosized_pacbp_encountered = True
                if nextPACBP.length == 0: zerosized_pacbp_encountered = True


        if zerosized_pacbp_encountered:
            # PacbP(ORF)(s) was/were completely consumed due to overlap correction
            delete_keys = []
            for (pacbpkey,nodeQ,nodeS), pacbpobj in self.pacbps.iteritems():
                if pacbpobj.length == 0:
                    delete_keys.append( (pacbpkey,nodeQ,nodeS) )
                else:
                    # final correction for unmatched ends
                    pacbpobj.strip_unmatched_ends()
            # remove these pacbps
            for (pacbpkey,nodeQ,nodeS) in delete_keys:
                ####################################################################
                if verbose:
                    print "PACBP deleted:", self.pacbps[(pacbpkey,nodeQ,nodeS)],
                    print nodeQ,nodeS
                ####################################################################
                del( self.pacbps[(pacbpkey,nodeQ,nodeS)] )
                # check if still pacbps of this type exist
                if not self.get_pacbps_by_nodes(node1=nodeQ,node2=nodeS):
                    # not the case -> delete the edge!
                    if self.has_edge(nodeQ,nodeS):
                        self.del_edge(nodeQ,nodeS)
    
        # return a meaningless True
        return True
    
    # end of function correct_overlaps



# end of class PacbpCollectionGraph


################################################################################
#### Helper functions for PacbpCollectionGraph class                        ####
################################################################################


def pacbporf2PCG(pacbporf,target,informant,PCG,source=''):
    """ """
    queryNode = (target,pacbporf.orfQ.id)
    sbjctNode = (informant,pacbporf.orfS.id)
    pacbporf.extend_pacbporf_after_stops()
    if source:
        pacbporf.source = source 
        pacbporf._gff['fsource'] = source 
    pacbpkey = pacbporf.construct_unique_key(queryNode,sbjctNode)
    PCG.add_node(queryNode)
    PCG.add_node(sbjctNode)
    PCG.add_edge(queryNode,sbjctNode,wt=pacbpkey[0])
    PCG.pacbps[(pacbpkey,queryNode,sbjctNode)] = pacbporf

# end of function pacbporf2PCG




def _delete_pacbp(PCG,key,pacbporf=None):
    """ """
    # pacbporf is not required in this function
    verbose = False

    # check if not already deleted in a previous block deletion!?
    if not PCG.pacbps.has_key(key): return False
    # check if delete-protected
    if PCG.pacbps[key]._IS_DELETE_PROTECTED: return False
    # if here, start the deletion step(s)
    (pacbpkey,nodeQ,nodeS) = key
    # delete from the PCG
    if verbose: print "_1_", ('afla',145)  in [ c for a,b,c in PCG.pacbps.keys() ], ('aory',144)  in [ c for a,b,c in PCG.pacbps.keys() ]
    del( PCG.pacbps[key] )
    # check if this edge is supported by >1 PacbPORF
    check = [ (n1,n2) == (nodeQ,nodeS) for k,n1,n2 in PCG.pacbps.keys() ]
    if verbose: print "_2_", ('afla',145)  in [ c for a,b,c in PCG.pacbps.keys() ], ('aory',144)  in [ c for a,b,c in PCG.pacbps.keys() ]
    if check.count(True) == 0:
        # single PacbP supports this edge -> delete edge
        if verbose: print "DELETE EDGE:", nodeQ, nodeS, PCG
        if verbose: print "_2a_", ('afla',145)  in [ c for a,b,c in PCG.pacbps.keys() ], ('aory',144)  in [ c for a,b,c in PCG.pacbps.keys() ]
        # if edge still there -> delete it
        #if PCG.has_edge(nodeQ,nodeS): PCG.del_edge(nodeQ,nodeS)
        try:    del( PCG.weights[(nodeQ,nodeS)] )
        except: pass
        try:    del( PCG.weights[(nodeS,nodeQ)] )
        except: pass
        try:    PCG.nodes[nodeQ].remove(nodeS)
        except: pass
        try:    PCG.nodes[nodeS].remove(nodeQ)
        except: pass
        if verbose: print "_2b_", ('afla',145)  in [ c for a,b,c in PCG.pacbps.keys() ], ('aory',144)  in [ c for a,b,c in PCG.pacbps.keys() ]
        # BE CAREFULL!! edges are reported as ordered ones
        # check if nodes can be deleted too
        for node in [nodeQ,nodeS]:
            if node in PCG.get_nodes():
                if not PCG.nodes[node]:
                    if verbose: print "_2c_", node, ('afla',145)  in [ c for a,b,c in PCG.pacbps.keys() ], ('aory',144)  in [ c for a,b,c in PCG.pacbps.keys() ]
                    if node in [ n1 for n1,n2 in PCG.weights.keys() ]:
                        pass
                    elif node in [ n2 for n1,n2 in PCG.weights.keys() ]:
                        pass
                    else:
                        print "___del_node(%s)" % (str(node))
                        PCG.del_node(node)
    if verbose: print "_3_", ('afla',145)  in [ c for a,b,c in PCG.pacbps.keys() ], ('aory',144)  in [ c for a,b,c in PCG.pacbps.keys() ]
    # return status of succesfull deletion
    return True

# end of function _delete_pacbp


def order_graphlist_by_total_weight(gralist):
    """
    Order a list with graphs on their total_weight() function outcome

    @type  gralist: []
    @param gralist: list with graphs objects (with a total_weight() function!)

    @rtype  gralist: []
    @retrun gralist: list with graphs objects ordered on total_weight()
    """
    # sort the subgraphs based on their total_weight
    tmp = []
    for sg in gralist:
        tmp.append( ( sg.total_weight(), sg ) )
    tmp.sort()
    tmp.reverse()
    return [ sg for (tw,sg) in tmp ]

# end of function order_graphlist_by_total_weight


def harvest_pacbps_from_crossdata(cbg,crossdata):
    """
    Harvest PacbP(ORF)s from crossdata to CodingBlockGraph/PacbpCollectionGraph

    @type  cbg: PacbpCollectionGraph or CodingBlockGraph
    @param cbg: object instance

    @type  crossdata: dict
    @param crossdata: crossdata <dict data structure>
    """
    # harvest the correct pacbs from crossdata
    for (g1,o1),(g2,o2) in cbg.weights.keys():
        combi = (g1,g2)
        if crossdata.has_key(combi):
            for (a,b,c,d) in crossdata[combi]['accepted_pacbs'].keys():
                if c==o1 and d==o2:
                    cbg.pacbps[ ( (a,b,c,d),(g1,o1),(g2,o2) ) ] =\
			    crossdata[combi]['accepted_pacbs'][(a,b,c,d)]

        else:
            # weights are given as AxB and BxA,
            # only one out of both combi's is present in crossdata!
            pass

# end of function harvest_pacbps_from_crossdata


def pacbpcollectiongraph2crossdata(pcg,increment_to_dict=None):
    """
    Create (incremental) crossdata dict structure from PacbpCollectionGraph

    @type  pcg: PacbpCollectionGraph
    @param pcg: PacbpCollectionGraph

    @type  increment_to_dict: crossdata dictionary or None
    @param increment_to_dict: crossdata dict to update or None (create novel)

    @rtype  crossdata: dict
    @return crossdata: crossdata <data structure dict> with PacbPs

    @attention: cd1  = pacbpcollectiongraph2crossdata(pcg1)
                pcg2 = crossdata2pacbpcollectiongraph(cd1)
                cd2  = pacbpcollectiongraph2crossdata(pcg2)
                cd1 != cd2
                all PacbP(ORF)s in the PCG are stores into the key
                accepted_pacbs in crossdata. All (extra/hidden) PacbP(ORF)s in
                other keys in crossdata are abandoned in conversions!
    """
    if increment_to_dict:
        # add new values to dict
        crossdata = increment_to_dict
    elif pcg._crossdata:
        # Assume crossdata in PCG to be insynch with PCG and its Pacbps
        # and return it here before any operation is performed on it
	return pcg._crossdata
    else:
        # create a new blank graph
        crossdata = {}
        for orgA, orgB in pcg.pairwisecrosscombinations_organism():
            crossdata[(orgA,orgB)] = deepcopy(CROSSDATA_SUBDICT)

    # Now fill the crossdata dict with PacbPs in the PCG.
    # Trust the PCG.pacbps dict attribute to be insync with nodes & edges in PCG
    for ( (a,b,c,d),(g1,o1),(g2,o2) ), pacbporf in pcg.pacbps.iteritems():
	crossdata[(g1,g2)]['accepted_pacbs'][(a,b,c,d)] = pacbporf
	
    # return crossdata dict
    return crossdata

# end of function pacbpcollectiongraph2crossdata


def _crossdatasubdict2edgelist(geneQ,geneS,subdict):
    """
    Get list of unique (PCG) edges supported by a subdict of crossdata

    @type  geneQ: * (string)
    @param geneQ: Organism identifier of the query Organism/Gene of subdict

    @type  geneS: * (string)
    @param geneS: Organism identifier of the sbjct Organism/Gene of subdict

    @type  subdict: dict
    @param subdict: subdict with PacbPs of crossdata

    @attention: keys of subdict are (bitscore,length,orfQid,orfSid)
    @attention: values of subdict are PacbP(ORF)s
    
    @rtype:  list
    @return: list with (PCG) edges present in a subdict with PacbPs of crossdata
    """
    edge_list = []
    for (a,b,c,d), pacbporf in subdict.iteritems():
	nodeQ = (geneQ,c) # c == query Orf id
	nodeS = (geneS,d) # d == sbjct Orf id
	edge_list.append( (nodeQ, nodeS) )
    return edge_list

# end of function _crossdatasubdict2edgelist




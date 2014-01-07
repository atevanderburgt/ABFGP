################################################################################
### PacbpCollectionGraph class                                              ####
################################################################################

# graphAbgp imports
from graph_organism import OrganismGraph
from graph_pssmcollections import DonorSiteCollectionGraph, AcceptorSiteCollectionGraph, TranslationalStartSiteCollectionGraph
from subclass_pacbpaccessibility import PacbpAccesibilityFunctions
from subclass_pssmobjects import BasalPSSMObjectGraphFunctions
from subclass_sitealignment import find_fully_connected_subgraphs, sort_by_cumulative_score
from exceptions import *
from validators import *
import conversion
import ordering

# Abgp Imports
from lib_intron import *

# Pacb imports
import pacb

# Abgp imports
from lib_clustalw import clustalw

# Python imports
from sets import Set
from copy import deepcopy

# Global variables


class ExonCollectionGraph(OrganismGraph,BasalPSSMObjectGraphFunctions,PacbpAccesibilityFunctions):
    """
    PacbpCollectionGraph (PCG) class, inheriting from OrganismGraph class.
    """
    def __init__(self):
        """
		Initialize a PacbpCollectionGraph
        """
        # Initialize as an OrganismGraph
        OrganismGraph.__init__(self)
        # set extra attributes
        self._node_object = {}
        self._node_pssm = {}

        # needed for backwards compatibilty with PacbpCollectionGraph, CodingBlockGraph
        self.pacbps = {}
        self._omsr = {}

    # end of function __init__

    def __str__(self):
        """ """
        return "<%s wt=%s lengths=(%s) [%s] >" % (
                self.__class__.__name__,
                self.total_weight(),
                ",".join(list(Set([ str(d-c) for (a,b,c,d) in self.get_nodes() ]))),
                " ".join([ "%s:%s(%s-%s)" % (a,b,c,d) for (a,b,c,d) in self.get_ordered_nodes() ])
                )
    # end of function __str__


    def _update_after_changes(self):
        """
        @attention: MUST BE OVERRIDDEN IN THE SUBCLASSES, if the returned value is not correct
        """
        self._update_edge_binary_entropies()
        self._update_node_object()
        self._update_node_pssm()

    # end of function _update_after_changes


    def _organism_from_node(self,node):
        """
        Get the organism identifier from the node identifier.

        @attention: MUST BE OVERRIDDEN IN THE SUBCLASSES, if the returned value is not correct

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

        @attention: MUST BE OVERRIDDEN IN THE SUBCLASSES, if the returned value is not correct

        @type  node: *
        @param node: One Node

        @rtype:  *
        @return: organism identifier (presumably a string)
        """
        return node[0]

    # end of function organism_by_node


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


    def cumulative_score(self):
        """ """
        return self.total_weight()

    # end of function cumulative_score


    def create_edges(self,max_length_nt_difference=15,max_length_ratio=0.5):
        """
        Create edges between nodes in the ExonCollectionGraph
        """
        for orgA,orgB in self.pairwisecrosscombinations_organism():
            for nodeA in self.get_organism_nodes(orgA):
                objA = self._node_object[nodeA]
                for nodeB in self.get_organism_nodes(orgB):
                    objB = self._node_object[nodeB]
                    if None in [ objA.donor.phase, objB.donor.phase ]:
                        # objA.donor and/or objB.donor are CodingBlockEnds
                        pass
                    elif objA.donor.phase == objB.donor.phase:
                        # objA.donor and objB.donor have identical phases
                        pass
                    else:
                        # incompatible phases -> no edge!
                        continue
                    # if here, then create an edge!
                    pssm_wt = objA.pssm_score + objB.pssm_score
                    dist    = abs( objA.length - objB.length )
                    reldist = float(dist)/float(max([objA.length,objB.length]))
                    if dist > max_length_nt_difference and reldist > max_length_ratio:
                        #  length difference to large!
                        continue

                    if None in [ objA.donor.phase, objB.donor.phase ] and dist <= 2:
                        # correct for phase-unequality between CodingBlockEnds and DonorSites
                        dist = 0
                    dist_wt = 1.0 / ( 1.0 + float(dist) )
                    self.add_edge(nodeA,nodeB,wt=pssm_wt*dist_wt)

    # end of function create_edges


    def find_fully_connected_subgraphs(self,edges=None,max_missing_edges=0,iterate=True,verbose=False):
        """
        Find all Fully Connected Graphs (FCG) in the input graph.

        @type  edges: number
        @param edges: number of outgoing edges of a node in a FCG

        @type  max_missing_edges: number
        @param max_missing_edges: number of missing edges to allow in a (nearly)-FCG

        @attention: max_missing_edges is not functional (max_missing_edges==0)

        @rtype:  list
        @return: list with (nearly) FCGs of the requested properties
        """
        subgraphs = self.recombine_into_completegraphs(edges=edges,verbose=verbose)
        if verbose: print "accepted recombined ECGs:", len(subgraphs)
        # if pacbps in parent (==self) ECG -> transfer elegiable pacbps 
        if self.pacbps:
            for sg in subgraphs:
                for node1,node2 in sg.pairwisecrosscombinations_node():
                    pacbps = self.get_pacbps_by_nodes(node1=node1,node2=node2)
                    if pacbps:
                        pacbp = pacbps[0]
                        pacbpkey = pacbp.construct_unique_key(node1,node2)
                        sg.pacbps[(pacbpkey,node1,node2)] = pacbp
        # return the fully connected ECGs
        return subgraphs

    # end of function find_fully_connected_subgraphs


    def recombine_into_completegraphs(self,edges=None,verbose=False):
        """
        Create all possible ExonCollectionGraphs by organism node recombination 

        @type  edges: number
        @param edges: number of outgoing edges of a node in a FCG

        @rtype:  list
        @return: list with ExonCollectionGraph of the requested properties
        """
        from codingblock_splitting import cross

        # if edges is not applied get by definition from the OrganismGraph
        if not edges: edges = self.organism_set_size() - 1

        # currently, this function has only a hard-set max_missing_edges == 0
        max_missing_edges = 0

        retlist = []
        # make a cross of all the pacbp positions in the lists of alternatives 
        allcombis = cross([ self.get_organism_nodes(org) for org in self.organism_set() ])
        if verbose: print "combinations:", len(allcombis)
        # gather a list of missing edges in the ECG
        missing_edges = []
        for (node1,node2) in self.pairwisecrosscombinations_node():
            if self._organism_from_node(node1) == self._organism_from_node(node2): continue
            if not self.has_edge(node1,node2):
                missing_edges.append((node1,node2))
        if verbose: print "missing edges:", len(missing_edges)

        # check for combinations that nodes that are listed as a missing edge
        # these are not relevant because max_missing_edges == 0
        if edges == self.organism_set_size() - 1:
            for pos in range(len(allcombis)-1,-1,-1):
                combi = allcombis[pos]
                for node1,node2 in missing_edges:
                    if node1 in combi and node2 in combi:
                        allcombis.pop(pos)
                        break
        if verbose: print "relevant:", len(allcombis)

        for combi in allcombis:
            sg = ExonCollectionGraph()
            for node in combi:
                # get the exon object from the main ExonCollectionGraph
                exon = self.get_node_object(node)
                # add node & object to the subgraph
                sg.add_node_and_object(node,exon)

            # create the edges in the subgraph
            for (node1,node2) in sg.pairwisecrosscombinations_node():
                if self.has_edge(node1,node2):
                    wt = self.get_edge_weight(node1,node2)
                    sg.add_edge(node1,node2,wt=wt)

            ## now check if is a succesfull recombination
            #if sg.node_count() != self.organism_set_size():
            #    continue
            #if sg.edge_count() < sum(range(0,self.organism_set_size())) - max_missing_edges:
            #    continue

            # remove nodes that have zero edges
            sg.remove_low_connectivity_nodes(min_connectivity=1)

            # do not check nr. of nodes on organism_set_size, but on variable `edges`!
            if sg.node_count() < edges+1:
                continue

            if sg.edge_count() < sum(range(0,edges+1)) - max_missing_edges:
                continue

            # if here -> accepted!
            if edges == self.organism_set_size() - 1:
                retlist.append(sg)
            else:
                # hmm... recombination with allowing missing organisms/nodes
                # that means that there are duplicates in the subgraphs
                # that come to this point. Check if this subgraph is already
                # present in retlist before addition
                for alt in retlist:
                    if len(alt.node_set().difference(sg.node_set())) == 0:
                        break
                else:
                    # not recognized -> add!
                    retlist.append(sg)

        # update the attributes dicts
        for sg in retlist: sg._update_after_changes()
        # and return a ordered/prioritized list
        return sort_by_cumulative_score(retlist)

    # end of function recombine_into_completegraphs


    def donor_phase(self):
        """
        Get the (uniform) phase of the donors in this ECG

        @attention: only use this function when ECG is NOT a FinalExon ECG

        @rtype:  integer
        @return: donor site phase (0,1,2)
        """
        if self.node_count() == 0:
            return None
        else:
            for obj in self._node_object.values():
                if obj.donor.__class__.__name__ != 'CodingBlockEnd':
                    return obj.donor.phase
            else:
                return None
                                
    # end of function donor_phase


    def acceptor_phase(self):
        """
        Get the (uniform) phase of the acceptors in this ECG

        @attention: only use this function when ECG is NOT a FirstExon ECG

		@rtype:  integer
		@return: acceptor site phase (0,1,2)
        """
        if self.node_count() == 0:
            return None
        else:
            for obj in self._node_object.values():
                if obj.acceptor.__class__.__name__ != 'CodingBlockStart':
                    return obj.acceptor.phase
            else:
                # all projected CodingBlockStart objects !?
                # this means a split in an existing Orf,
                # a split in an AA sequence, so phase == 0
                return 0
                                
    # end of function acceptor_phase



    def create_projected_donorsite(self,node):
        """
        @type  node: *
        @param node: One Node
        """
        pass

    # end of function create_projected_donorsite


    def create_projected_acceptorsite(self,organism):
        """ """
        pass
    # end of function create_projected_acceptorsite


    def get_orfs_of_graph(self,organism=None,node=None):
        """
        Get all the Orf objects of this graph

        @type  organism: *
        @param organism: organism identifier (or None)

        @type  node: *
        @param node: node identifier (or None)

        @rtype:  dictionary (or list if organism is specified)
        @return: dictionary with organisms (keys) and list of Orf objects (values),
                 or only list of Orf objects if an organism or node identifier was specified.
        """
        if organism and organism not in self.organism_set():
            raise OrganismNotPresentInGraph
        elif node and node not in self.get_nodes():
            raise NodeNotPresentInGraph
        elif organism:
            orfs = []
            for node, exon in self._node_object.iteritems():
                if organism == self.organism_by_node(node):
                    orfs.append( exon.orf )
            return orfs
        elif node:
            return [ self._node_object[node].orf ]
        else:
            orfs = {}
            for org in self.organism_set(): orfs[org] = []
            for node, exon in self._node_object.iteritems():
                org = self.organism_by_node(node)
                orfs[org].append( exon.orf )
            return orfs

    # end of funtion get_orfs_of_graph


    def overall_minimal_spanning_range(self):
        """ """
        return self._omsr

    # end of function overall_minimal_spanning_range


    def make_pacbps_for_edges(self,aa_extra_offset=1):
        """
        """
        make_pacbps_for_edges(self,aa_extra_offset=aa_extra_offset)

    # end of function make_pacbps_for_edges


    def ecgnode2cbgnode(self,ecgnode):
        """
        convert ECGnode to CBGnode
        """
        return ( ecgnode[0], ecgnode[1] )

    # end of function ecgnode2cbgnode 

# end of class ExonCollectionGraph

################################################################################
#### Helper functions for ExonCollectionGraph class                         ####
################################################################################


def ExonCollectionGraph2TranslationalStartSiteCollectionGraph(gra):
    """
    Convert ECG -> TranslationalStartSiteCollectionGraph

    @attention: only in use when ECG is a FirstExon ECG

    @rtype:  TranslationalStartSiteCollectionGraph
    @return: TranslationalStartSiteCollectionGraph instance to be placed in the CBG
    """
    newgra = TranslationalStartSiteCollectionGraph()
    newgra.ALIGNED_SITE_AA_OFFSET = 10
    newgra.MIN_PSSM_SCORE = -0.0
    for node in gra.get_nodes():
        # exon node is ( org, orf, ntpos), TSS node ( org, orf, aapos, ntpos )
        newnode = ( node[0], node[1], node[2]/3, node[2] )
        tss = gra._node_object[node].acceptor
        newgra.add_node_and_object(newnode,tss)
        newgra._node_pssm[newnode] = tss.pssm_score
    for nodeA,nodeB in newgra.pairwisecrosscombinations_node():
        entropyQ = 1.0
        entropyS = 1.0
        newgra.add_edge(nodeA,nodeB,wt=1.0)
        newgra._edge_binary_entropies[(nodeA,nodeB)] = (entropyQ,entropyS)
        newgra._edge_binary_entropies[(nodeB,nodeA)] = (entropyS,entropyQ)

    # Get tcode data for these start codon nodes
    # Assuming that this is indeed the start-codon,
    for (org,orfid,aaPos,dnaPos) in newgra.get_nodes():
        theorf = gra.get_orfs_of_graph(organism=org)[0]
        # calculate the average TCODE scores for the windows
        ( tcode5p,tcode3p ) = theorf.tcode_entropy_of_pos(
                aaPos,
                window_left=newgra._TCODE_5P_WINDOWSIZE,
                window_right=newgra._TCODE_3P_WINDOWSIZE,
                )
        newgra._tcode5pscore[(org,orfid,aaPos,dnaPos)] = tcode5p
        newgra._tcode3pscore[(org,orfid,aaPos,dnaPos)] = tcode3p
    # return the TranslationalStartSiteCollectionGraph
    return newgra

# end of function ExonCollectionGraph2TranslationalStartSiteCollectionGraph


def ExonCollectionGraph2DonorSiteCollectionGraph(gra):
    """
    Convert ECG -> DonorSiteCollectionGraph

    @attention: only in use when ECG is NOT a FinalExon ECG

    @rtype:  DonorSiteCollectionGraph
	@return: DonorSiteCollectionGraph instance to be placed in the CBG
    """
    newgra = DonorSiteCollectionGraph()
    newgra.ALIGNED_SITE_AA_OFFSET = 10
    newgra.MIN_PSSM_SCORE = -0.0
    for node in gra.get_nodes():
        donor = gra._node_object[node].donor
        if donor.__class__.__name__ == 'CodingBlockEnd':
            phase = gra.donor_phase()
            # return a ProjectedSpliceSite
            projDonor = CodingBlockEnd2ProjectedSpliceDonor(donor,phase=phase)
            newnode = ( node[0], node[1], projDonor.pos )
            newgra.add_node_and_object(newnode,projDonor)
            newgra._node_pssm[newnode] = donor.pssm_score
        else:
            newnode = ( node[0], node[1], donor.pos )
            newgra.add_node_and_object(newnode,donor)
            newgra._node_pssm[newnode] = donor.pssm_score
    for nodeA,nodeB in newgra.pairwisecrosscombinations_node():
        newgra.add_edge(nodeA,nodeB,wt=1.0,entropy=1.0)
    # return the donorsitecollection
    return newgra

# end of function ExonCollectionGraph2DonorSiteCollectionGraph


def ExonCollectionGraph2AcceptorSiteCollectionGraph(gra):
    """
    Convert ECG -> AcceptorSiteCollectionGraph

    @attention: only in use when ECG is NOT a FirstExon ECG

    @rtype:  AcceptorSiteCollectionGraph
    @return: AcceptorSiteCollectionGraph instance to be placed in the CBG
    """
    newgra = AcceptorSiteCollectionGraph()
    newgra.ALIGNED_SITE_AA_OFFSET = 10
    newgra.MIN_PSSM_SCORE = -0.0
    for node in gra.get_nodes():
        accep = gra._node_object[node].acceptor
        if accep.__class__.__name__ == 'CodingBlockStart':
            phase = gra.acceptor_phase()
            # return a ProjectedSpliceSite
            projAccep = CodingBlockStart2ProjectedSpliceAcceptor(accep,phase=phase)
            newnode = ( node[0], node[1], projAccep.pos )
            newgra.add_node_and_object(newnode,projAccep)
            newgra._node_pssm[newnode] = accep.pssm_score
        else:
            newnode = ( node[0], node[1], accep.pos )
            newgra.add_node_and_object(newnode,accep)
            newgra._node_pssm[newnode] = accep.pssm_score
    for nodeA,nodeB in newgra.pairwisecrosscombinations_node():
        entropyQ = 1.0
        entropyS = 1.0
        newgra.add_edge(nodeA,nodeB,wt=1.0)
        newgra._edge_binary_entropies[(nodeA,nodeB)] = (entropyQ,entropyS)
        newgra._edge_binary_entropies[(nodeB,nodeA)] = (entropyS,entropyQ)
    return newgra

# end of function ExonCollectionGraph2AcceptorSiteCollectionGraph 


def ExonCollectionGraph2CodingBlockGraph(gra,is_first=False,is_last=False,firstCBG=None,lastCBG=None):
    """
    Convert ECG -> CodingBlockGraph

    @type  gra: ExonCollectionGraph
    @param gra: ExonCollectionGraph instance

    @type  is_first: Boolean
    @param is_first: True or False (default False)

    @type  is_last: Boolean
    @param is_last: True or False (default False)

    @type  firstCBG: CodingBlockGraph (or None)
    @param firstCBG: ...

    @type  lastCBG: CodingBlockGraph (or None)
    @param lastCBG: ...

    @attention: make shure to specify arguments correctly!

    @rtype:  CodingBlockGraph
    @return: CodingBlockGraph instance or None when failed!
    """
    from graph_codingblock import CodingBlockGraph
    cbg = CodingBlockGraph()
    cbg.MINIMAL_OVERAL_SPANNING_RANGE_SIZE = 1
    for ecgnode in gra.get_nodes():
        cbg.add_node( gra.ecgnode2cbgnode(ecgnode) )


    # make Pacbp objects for the edges if not done yet
    if not gra.pacbps or len(gra.pacbps) != gra.edge_count():
        gra.make_pacbps_for_edges()

    # check if number of pacbps matches number of edges
    if len(gra.pacbps) != gra.edge_count():
        # pacbp creation failed for at least 1 edge -> no CBG!
        return None

    # check if connectivitysaturation == 1.0:
    if gra.connectivitysaturation() != 1.0:
        # no edge listed between some of the nodes!
        return None

    # transfer pacbps and edges to new CBG
    for (key,ecgnode1,ecgnode2), pacbporf in gra.pacbps.iteritems():
        # convert ECGnode to CBGnode
        cbg_node1 = gra.ecgnode2cbgnode(ecgnode1)
        cbg_node2 = gra.ecgnode2cbgnode(ecgnode2)
        bitscore = key[0]
        cbg.add_edge(cbg_node1,cbg_node2,wt=bitscore)
        cbg.pacbps[(key,cbg_node1,cbg_node2)] = pacbporf

    # fix pacbps that are already present in what is currently
    # the first CBG. This is the case for nodes (exons) in the firstExonGraph
    # that end with a CodingBlockEnd.
    if firstCBG:
        replacements = {}
        for (key,node1,node2), pacbporf in cbg.pacbps.iteritems():
            if firstCBG.has_edge(node1,node2):
                startPos = pacbporf._get_original_alignment_pos_start()
                firstpacbporf = firstCBG.get_pacbps_by_nodes(node1=node1,node2=node2)[0]
                firstStartPos = firstpacbporf._get_original_alignment_pos_start()
                if firstStartPos.query_pos <= startPos.query_pos and firstStartPos.sbjct_pos <= startPos.sbjct_pos:
                    replacements[(key,node1,node2)] = firstpacbporf
        if replacements:
            for (key,node1,node2), pacbporf in replacements.iteritems():
                del( cbg.pacbps[(key,node1,node2)] )
                newkey = pacbporf.construct_unique_key(node1,node2)
                cbg.pacbps[(newkey,node1,node2)] = pacbporf

    if lastCBG:
        replacements = {}
        for (key,node1,node2), pacbporf in cbg.pacbps.iteritems():
            if lastCBG.has_edge(node1,node2):
                endPos = pacbporf._get_original_alignment_pos_end()
                lastpacbporf = lastCBG.get_pacbps_by_nodes(node1=node1,node2=node2)[0]
                lastEndPos   = lastpacbporf._get_original_alignment_pos_end()
                if lastEndPos.query_pos >= endPos.query_pos and lastEndPos.sbjct_pos <= endPos.sbjct_pos:
                    replacements[(key,node1,node2)] = lastpacbporf
        if replacements:
            for (key,node1,node2), pacbporf in replacements.iteritems():
                del( cbg.pacbps[(key,node1,node2)] )
                newkey = pacbporf.construct_unique_key(node1,node2)
                cbg.pacbps[(newkey,node1,node2)] = pacbporf

    # Now make shure this CodingBlockGraph is fully compatible in the Genestructure
    # That means, do some site scanning because this newcbg is inserted AFTER
    # all the sitescanning has been done!
    if not cbg.has_overall_minimal_spanning_range():
        ### print "ecg2cbg:", cbg, cbg.edge_count(), len(cbg.pacbps)
        return None 


    # set the footprint of where this CBG came from
    # this can be deleted lateron, but is not required per se
    # this footprint is needed because the surrounding CBGs
    # must get some (additional) forced splice sites that depend
    # on the sites in the ExonCollectionGraph 
    cbg._ExonCollectionGraph = gra

    # DO NOT update edge weights!!
    # this wipes out the alignment evidence outside of the OMSR region
    #cbg.update_edge_weights_by_minimal_spanning_range()

    # make AlignedStopCodonGraph
    cbg.align_stop_codons()


    # create splicedonorgraph & align sites
    cbg._splicedonorgraph = ExonCollectionGraph2DonorSiteCollectionGraph(gra)
    cbg._splicedonorgraph.collection2alignedsites(
            edges=gra.node_count()-1,
            minimal_edges=gra.node_count()-1
            )

    if is_first:
        # create tssgraph & align sites
        dummyTSScolgra = ExonCollectionGraph2TranslationalStartSiteCollectionGraph(gra)
        dummyTSScolgra.collection2alignedsites( edges=gra.node_count()-1 )
        # make the *TRUE* TSScolgra and align them
        cbg.harvest_elegiable_tss_sites()
        cbg._startcodongraph.collection2alignedsites()
        # place the dummyTSScolgra in front of the que
        cbg._startcodongraph.alignedsites.insert( 0, dummyTSScolgra.alignedsites[0] )
        # and assign _codingblockgraph attribute to all alignedsites;
        # this is needed to do an is_optimal() check
        for _algtss in cbg._startcodongraph.alignedsites:
            _algtss._codingblockgraph = cbg

    else:
        # create spliceacceptorgraph & align sites
        cbg._spliceacceptorgraph = ExonCollectionGraph2AcceptorSiteCollectionGraph(gra)
        cbg._spliceacceptorgraph.collection2alignedsites(
                edges=gra.node_count()-1,
                minimal_edges=gra.node_count()-1
                )

    if is_first:
        # make SpliceSiteCollectionGraphs
        ### print "ecg2cbg:(1)", cbg
        cbg.harvest_elegiable_acceptor_sites(projected_acceptors={},forced_codingblock_ends={},prev=None)
        
        # do site alignment of acceptors
        cbg._spliceacceptorgraph.collection2alignedsites(edges=cbg.node_count()-1,minimal_edges=2)

    if is_last:
        # make SpliceSiteCollectionGraphs
        cbg.harvest_elegiable_donor_sites(projected_donors={},forced_codingblock_ends={},next=None)

        # do site alignment of acceptors
        cbg._splicedonorgraph.collection2alignedsites(edges=cbg.node_count()-1,minimal_edges=2)

    # done! return the CBG
    ### print "ecg2cbg:done!", cbg
    return cbg

# end of function ExonCollectionGraph2CodingBlockGraph 


def CodingBlockStart2ProjectedSpliceAcceptor(cbgstart,phase=None):
    """
    Make a ProjectedSpliceAcceptor for the AcceptorGraph of the ExonCollectionGraph -> CBG
    """
    # input validation
    IsProperPhaseValidator(phase)
    IsCodingBlockStartValidator(cbgstart)

    # return a ProjectedSpliceSite
    position = cbgstart.pos - (3-phase)
    dummyAccep = SpliceAcceptorAG(position,phase=phase)
    return ProjectedSpliceAcceptor(position,"XX",
                    dummyAccep,
                    distance=0,
                    entropy=1.0,
                    pssm_score=cbgstart.pssm_score,
                    cases=1)

# end of function CodingBlockStart2ProjectedSpliceAcceptor


def CodingBlockEnd2ProjectedSpliceDonor(cbgend,phase=None):
    """
    Make a ProjectedSpliceDonor for the DonorGraph of the ExonCollectionGraph -> CBG
    """
    IsProperPhaseValidator(phase)
    IsCodingBlockEndValidator(cbgend)

    # return a ProjectedSpliceSite
    position = cbgend.pos - (3-phase)
    dummyDonor = SpliceDonorGT(position,phase=phase)
    return ProjectedSpliceDonor(position,"XX",
                    dummyDonor,
                    distance=0,
                    entropy=1.0,
                    pssm_score=cbgend.pssm_score,
                    cases=1)

# end of function CodingBlockEnd2ProjectedSpliceDonor


def CodingBlockEnd2ProjectedSpliceAcceptor(cbgend,phase=None):
    """
    Make a ProjectedSpliceAceptor for the NEXT graph (which was the original IS_FIRST)
    """
    IsProperPhaseValidator(phase)
    IsCodingBlockEndValidator(cbgend)

    # return a ProjectedSpliceSite
    position = cbgend.pos - (3-phase)
    dummyAccep = SpliceAcceptorAG(position,phase=phase)
    return ProjectedSpliceAcceptor(position,"XX",
                    dummyAccep,
                    distance=0,
                    entropy=1.0,
                    pssm_score=cbgend.pssm_score,
                    cases=1)

# end of function CodingBlockEnd2ProjectedSpliceAcceptor


def CodingBlockStart2ProjectedSpliceDonor(cbgstart,phase=None):
    """
    Make a ProjectedSpliceDonor for the PREVIOUS graph (which was the original IS_LAST)
    """
    IsProperPhaseValidator(phase)
    IsCodingBlockStartValidator(cbgstart)

    # return a ProjectedSpliceSite
    position = cbgstart.pos - (3-phase)
    dummyDonor = SpliceDonorGT(position,phase=phase)
    return ProjectedSpliceDonor(position,"XX",
                    dummyDonor,
                    distance=0,
                    entropy=1.0,
                    pssm_score=cbgstart.pssm_score,
                    cases=1)

# end of function CodingBlockStart2ProjectedSpliceDonor



def make_pacbps_for_edges(gra,aa_extra_offset=1,verbose=False):
    """
    """
    coordsandseqs = {}
    # create dummy omsr attribute!
    # omsr is filled with Exon.acceptor and Exon.donor positions
    # recalculate nt positions to aa positions!
    for node in gra.get_ordered_nodes():
        accep = gra._node_object[node].acceptor
        donor = gra._node_object[node].donor
        aaStart = accep.pos / 3
        if donor.pos - accep.pos % 3 == 0:
            aaEnd = donor.pos / 3
        else:
            aaEnd = (donor.pos / 3) +1
        # get orf, seequence coordinates and sequence itself
        theorg = gra._organism_from_node(node)
        theorf = gra.get_orfs_of_graph(node=node)[0]
        aaStart -= aa_extra_offset
        aaEnd   += aa_extra_offset
        # correct end coordinates when falling outside of Orf
        if aaEnd > theorf.protein_endPY: aaEnd = theorf.protein_endPY
        if aaStart < theorf.protein_startPY: aaStart = theorf.protein_startPY
        theseq = theorf.getaas(abs_pos_start=aaStart,abs_pos_end=aaEnd)
        # store to dict
        coordsandseqs[node] = (theseq,theorg,theorf,aaStart,aaEnd)


    for (node1,node2) in gra.pairwisecrosscombinations_node():
        # check if these are nodes present as an edge
        if not gra.has_edge(node1,node2): continue

        # start makeing a Pacbp from clustalw
        (seq1,org1,orf1,aa1start,aa1end) = coordsandseqs[node1]
        (seq2,org2,orf2,aa2start,aa2end) = coordsandseqs[node2]

        # create headers and fetch sequences from Orf objects
        header1  = "%s_orf_%s_%s_%s" % (org1,orf1.id,aa1start,aa1end)
        header2  = "%s_orf_%s_%s_%s" % (org2,orf2.id,aa2start,aa2end)

        # check if sequences exist/ at least 1 AA
        if not seq1 and not seq2:
            print "Warning: ZeroProteinSequenceLengthException", "S1", aa1start, aa1end, node1, node2, orf1
            print "Warning: ZeroProteinSequenceLengthException", "S2", aa2start, aa2end, node1, node2, orf2
            continue
        elif not seq2:
            print "Warning: ZeroProteinSequenceLengthException", "S2", aa2start, aa2end, node1, node2, orf2
            continue
        elif not seq1:
            print "Warning: ZeroProteinSequenceLengthException", "S1", aa1start, aa1end, node1, node2, orf1
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
            # update edge weight
            #new_wt = pacbporf.bitscore
            # wt was sum(PSSM) * distance ratio
            # now multiply with identityscore (0.0-1.0) float too
            new_wt = pacbporf.identityscore * gra.get_edge_weight(node1,node2)
            gra.set_edge_weight(node1,node2,wt=new_wt)
            # add pacbporf to CBG
            key = pacbporf.construct_unique_key(node1,node2)
            gra.pacbps[(key,node1,node2)] = pacbporf
        else:
            # pacbp.conversion.pacbp_from_clustalw did
            # not yield any proper alignment
            if verbose: print "NO PACBP!!", node1,node2, seq1,seq2
            pass

# end of function make_pacbps_for_edges


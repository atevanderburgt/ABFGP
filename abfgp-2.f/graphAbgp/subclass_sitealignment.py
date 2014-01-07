"""
Sitealignment functions joined in a class that is inherited in graph classes for Aligment Based Gene Predictions
"""

# Module metadata
__authors__ = "Ate van der Burgt"
__license__ = "MIT"

# graphAbgp imports
from exceptions import *
import ordering

# Python Imports
from sets import Set
from copy import deepcopy


class BasalSiteAlignmentFunctions:
    """
    """
    def _handle_edges_argument(self,edges):
        """
        """
        # if edges is not applied get by definition form the OrganismGraph
        if not edges: edges = self.organism_set_size() - 1
        # set edges to a minimum of 2
        edges = max([2,edges])
        # return edges argument
        return edges

    # end of function _handle_edges_argument


    def find_fully_connected_subgraphs(self,edges=None,max_missing_edges=2):
        """
        Find all Fully Connected Graphs (FCG) in the input graph.

        @type  edges: number
        @param edges: number of outgoing edges of a node in a FCG

        @type  max_missing_edges: number
        @param max_missing_edges: number of missing edges to allow in a (nearly)-FCG

        @rtype:  list
        @return: list with (nearly) FCGs of the requested properties
        """
        # handle edges argument (error check etc.)
        edges = self._handle_edges_argument(edges)

        subgraphs = find_fully_connected_subgraphs(self,edges=edges,max_missing_edges=max_missing_edges)

        # TODO: THIS FUNCTIONALLITY IS NOT GENERIC FOR THIS CLASS!
        # TODO: move it to the desired subclasses....
        # update the subgraphs by removing all that is not required anymore!
        for sg in subgraphs:
            sg._edge_binary_entropies = deepcopy(self._edge_binary_entropies)
            sg._update_after_changes()
        # and return a ordered/prioritized list
        return sort_by_cumulative_score(subgraphs)

    # end of function find_fully_connected_subgraphs


    def find_perfectly_conserved_sites(self,edges=None):
        """
        Find all 100% conserved sites in the input graph.
    
        @type  edges: number
        @param edges: number of outgoing edges of a node in a Fully Connected Graph (FCG)
    
        @rtype:  list
        @return: list with (nearly) FCGs of the requested properties
        """
        # handle edges argument (error check etc.)
        edges = self._handle_edges_argument(edges)

        subgraphs = find_perfectly_conserved_sites(self,edges=edges)
        # update the subgraphs by removing all that is not required anymore!
        for sg in subgraphs:
            sg._edge_binary_entropies = deepcopy(self._edge_binary_entropies)
            sg._update_after_changes()
        # and return a ordered/prioritized list
        return sort_by_cumulative_score(subgraphs)

    # end of function find_perfectly_conserved_sites


    def find_conserved_sites(self,edges=None):
        """
        Split the input graph in alignable/conserved sites
    
        @type  edges: number
        @param edges: number of outgoing edges of a node in a Fully Connected Graph (FCG)
    
        @rtype:  list
        @return: list with (nearly) FCGs of the requested properties
        """
        # handle edges argument (error check etc.)
        edges = self._handle_edges_argument(edges)

        subgraphs = find_conserved_sites(self,edges=edges)
        # update the subgraphs by removing all that is not required anymore!
        for sg in subgraphs:
            sg._edge_binary_entropies = deepcopy(self._edge_binary_entropies)
            sg._update_after_changes()
        # and return a ordered/prioritized list
        return sort_by_cumulative_score(subgraphs)

    # end of function find_conserved_sites


    def align_all_sites(self,edges=None,minimal_edges=2):
        """
        Align all sites in a graph iteratively.
    
        @type  edges: number
        @param edges: number of outgoing edges of a node in a Fully Connected Graph (FCG)
    
        @type  minimal_edges: number
        @param minimal_edges: number of minimal outgoing edges, then stop the iterative search
    
        @rtype:  list
        @return: list of subgraphs with aligned sites
        """
        # handle edges argument (error check etc.)
        edges = self._handle_edges_argument(edges)

        subgraphs = align_all_sites(self,edges=edges,minimal_edges=minimal_edges)
        # update the subgraphs by removing all that is not required anymore!
        for sg in subgraphs:
            sg._edge_binary_entropies = deepcopy(self._edge_binary_entropies)
            sg._update_after_changes()
        # and return a ordered/prioritized list
        return sort_by_cumulative_score(subgraphs)

    # end of function align_all_sites

# end of class BasalSiteAlignmentFunctions


###########################################################################
### helper functions for class BasalSiteAlignmentFunctions              ###
###########################################################################


def find_fully_connected_subgraphs(gra,edges=None,max_missing_edges=2,iterate=True):
    """
    Find all Fully Connected Graphs (FCG) in the input graph.

    @type  gra: OrganismGraph or inheriting Graph class
    @param gra: input graph

    @type  edges: number
    @param edges: number of outgoing edges of a node in a FCG

    @type  max_missing_edges: number
    @param max_missing_edges: number of missing edges to allow in a (nearly)-FCG

    @rtype:  list
    @return: list with (nearly) FCGs of the requested properties
    """
    if gra.connectivitysaturation() == 1.0:
        # graph is already fully connected
        return [ gra ]

    if not edges or type(edges) != type(int()) or edges < 0:
        message = "apply argument edges as a positive integer"
        raise InproperlyAppliedArgument, message

    if sum(gra.node_connectivity_counts()[edges:]) < edges+1:
        # nope, not enough nodes with enough
        # edges that can potentially contain
        # a fully connected subgraph
        return [ gra ]

    # total edges in a fully connected graph WITH hub
    edges_with_hub    = sum(range(1,edges+1))

    # (return)list for fully_connected_subgraphs
    connected_sgs = []

    # make a (new) and empty graph
    empty = deepcopy(gra)
    while empty.get_nodes(): empty.del_node( empty.get_nodes()[0] )

    for thenode in gra.get_nodes():
        # check if this node has potentially enough connections
        if len(gra.get_node(thenode)) < edges:
            continue

        # check if all organisms are present in the connections
        if len(Set([ gra._organism_from_node(n) for n in gra.get_node(thenode) ] )) < edges:
            continue

        nodesperorg = { gra._organism_from_node(thenode): [ thenode ] }
        for link in gra.get_node(thenode):
            org = gra._organism_from_node(link)
            if nodesperorg.has_key(org):
                nodesperorg[org].append( link )
            else:
                nodesperorg[org] = [ link ]

        nodecombis = [[]]
        for org,orgnodes in nodesperorg.iteritems():
            if len(orgnodes) == 1:
                for combi in nodecombis: combi.append(orgnodes[0])
            else:
                # more than one node from this organism
                # that means, double,triple,.. the amount of
                # nodecombis that are possible
                currentcombis = deepcopy(nodecombis)
                for combi in nodecombis: combi.append(orgnodes[0])
                for i in range(1,len(orgnodes)):
                    thiscombis = deepcopy(currentcombis)
                    for combi in thiscombis: combi.append(orgnodes[i])
                    nodecombis.extend( thiscombis )

        # now we have unique combination of nodes
        for nodecombi in nodecombis:
            sg = deepcopy(empty)
            sg.add_nodes( nodecombi )
            for n1 in sg.get_nodes():
                for n2 in sg.get_nodes():
                    if gra.has_edge(n1,n2):
                        sg.add_edge(n1,n2,gra.get_edge_weight(n1,n2))

            # now see if this is a subgraph to return
            if sg.max_edge_count() - sg.edge_count() <= max_missing_edges:
                nodecombi = Set(nodecombi)
                check = [ nodecombi.difference(prev.get_nodes()) == Set() for prev in connected_sgs ]
                if not True in check:
                    connected_sgs.append( sg )
            else:
                # To little edges to be accepted.
                # NEWNEW 06/04/2009
                # BUT: when Org(n)=x and graphs of less than x organisms are requested,
                # this is not per se a problem. Just iterate through this function
                # in search for fullfillment of the criterion <= max_missing_edges
                # IMPORTANT: do not forget to set iterate=False; otherwise eternal looping is achieved!
                if iterate and sg.node_count() > edges+1:
                    one_edge_less = sg.node_count()-2
                    for subsg in find_fully_connected_subgraphs(sg,edges=one_edge_less,max_missing_edges=0,iterate=False):
                        if subsg.node_count() >= edges+1 and subsg.max_edge_count() - subsg.edge_count() == 0:
                            nodecombi = Set( subsg.get_nodes() )
                            check = [ nodecombi.difference(seen.get_nodes()) == Set() for seen in connected_sgs ]
                            if not True in check:
                                connected_sgs.append( subsg )

    # NEWNEW 03/03/2009
    # and finally, append the remainer of the input graph
    for sg in connected_sgs:
        for node in sg.get_nodes():
            if node in gra.get_nodes(): gra.del_node( node )
    # check if something is remaining from input gra
    gra.remove_low_connectivity_nodes(min_connectivity=1)
    if gra.get_nodes(): connected_sgs.append(gra)

    # return the fully connected subgraphs
    return connected_sgs

# end of function find_fully_connected_subgraphs


def find_conserved_sites(gra,edges=None):
    """
    Split the input graph in alignable/conserved sites

    @type  gra: StartCodonGraph, StopCodonGraph or AlignedSpliceSiteGraph
    @param gra: input graph

    @type  edges: number
    @param edges: number of outgoing edges of a node in a Fully Connected Graph (FCG)

    @rtype:  list
    @return: list with (nearly) FCGs of the requested properties
    """
    if not edges or type(edges) != type(int()) or edges < 0:
        message = "apply argument edges as a positive integer"
        raise InproperlyAppliedArgument, message

    # first, find all perfectly conserved sites
    subgraphs = gra.find_perfectly_conserved_sites(edges=edges)
    # deepcopy the input graph and delete all nodes present
    # in the perfectly conserved sites in the graphs of `subgraphs`
    new = deepcopy(gra)
    for sg in subgraphs:
        for node in sg.get_nodes():
            new.del_node(node)

    # find fully connected graphs (but not per se perfectly conserved)
    fcgs = find_fully_connected_subgraphs(new,edges=edges,max_missing_edges=0)
    subgraphs.extend( fcgs )

    # remove empty graphs
    for i in range(len(subgraphs)-1,-1,-1):
        if not subgraphs[i].get_nodes():
            emptygra = subgraphs.pop(i)

    # and return 
    return subgraphs

# end of function find_conserved_sites


def find_perfectly_conserved_sites(gra,edges=None):
    """
    Find all 100% conserved sites in the input graph.

    @type  gra: StartCodonGraph, StopCodonGraph or AlignedSpliceSiteGraph
    @param gra: input graph

    @type  edges: positive integer
    @param edges: number of outgoing edges of a node in a Fully Connected Graph (FCG)

    @rtype:  list
    @return: list with (nearly) FCGs of the requested properties
    """
    if not edges or type(edges) != type(int()) or edges < 0:
        mesage = "apply argument edges as a positive integer"
        raise InproperlyAppliedArgument, message

    # deepcopy the input graph, because we are about to delete
    # nodes and edges from it
    new = deepcopy(gra)
    # find and delete all edges that do not have wt==1.0
    # wt=1.0 represent an alignment offset of 0nt (perfectly conserved)
    del_edges = []
    for edge,wt in new.weights.iteritems():
        if wt < 1.0: del_edges.append(edge)
    for (n1,n2) in del_edges:
        if new.has_edge(n1,n2):
            new.del_edge(n1,n2)

    # remove all components that cannot belong to a complete K(s) graph, with s==edges+1
    new.remove_incomplete_components(s=edges+1)
    #print "FPCS:", gra.node_count(), new.node_count(), gra.edge_count(), new.edge_count()
    #for org in gra.organism_set(): print ( org, len(gra.get_organism_nodes(org) ) ),
    #print ""

    subgraphs = []
    # split into Fully Connected Graphs (FCG)
    for sg in new.find_fully_connected_subgraphs(edges=edges,max_missing_edges=0):
        if sg.node_count() == edges+1 and sg.connectivitysaturation() == 1.0:
            # yep, Fully Connected and perfectly conserved
            subgraphs.append(sg)

    # and return
    return subgraphs

# end of function find_perfectly_conserved_sites


def align_all_sites(gra,edges=None,minimal_edges=2):
    """
    Align all sites in a graph iteratively.

    @type  gra: StartCodonGraph, StopCodonGraph or AlignedSpliceSiteGraph
    @param gra: input graph

    @type  edges: number
    @param edges: number of outgoing edges of a node in a Fully Connected Graph (FCG)

    @type  minimal_edges: number
    @param minimal_edges: number of minimal outgoing edges, then stop the iterative search

    @rtype:  list
    @return: list of subgraphs with aligned sites
    """
    if not edges or type(edges) != type(int()) or edges < 0:
        mesage = "apply argument edges as a positive integer"
        raise InproperlyAppliedArgument, message

    # calculate a range of edges
    edgerange = range(minimal_edges,edges+1)
    edgerange.reverse()

    # return list with aligned graphs
    aligned_sites = []

    # first, find all perfectly conserved sites
    # after an interation of X edges, remove all
    # nodes in all the found subgraphs
    new = deepcopy(gra)

    print "AASfunction:", new.node_count(),

    for edge in edgerange:
        pcs = find_perfectly_conserved_sites(new,edges=edge)
        aligned_sites.extend( pcs )
        for sg in pcs:
            for node in sg.get_nodes():
                new.del_node(node)

    print new.node_count(), 

    # second, find all non-perfectly conserved sites
    # after an interation of X edges, remove all
    # nodes in all the found subgraphs
    new = deepcopy(gra)
    for edge in edgerange:
        cs = find_conserved_sites(new,edges=edge)
        aligned_sites.extend( cs )
        for sg in cs:
            for node in sg.get_nodes():
                new.del_node(node)
        print (edge,new.node_count()),
    print ""

    # and return all the aligned sites
    return aligned_sites

# end of function align_all_sites


def sort_by_cumulative_score(sites):
    """
    """
    ret_sites = []
    for site in sites:
        ret_sites.append( ( site.cumulative_score(), site ) )
    ret_sites.sort()
    ret_sites.reverse()
    return [ site for (score,site)in ret_sites ]

# end of sort_by_cumulative_score

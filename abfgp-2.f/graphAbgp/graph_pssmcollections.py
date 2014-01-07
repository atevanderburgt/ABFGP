################################################################################
### PssmObjectCollectionGraph class and inheriting classes                  ####
################################################################################

# Python imports
from sets import Set

# Make sure Numeric is installed somewhere!!!
# For convenience, it is installed in ./requiredmodules
import sys
from os.path import join as osPathJoin
from settings.abgp import MAIN_ABGP_PATH as BASEPATH
sys.path.append(osPathJoin(BASEPATH,"requiredmodules"))
from Numeric import zeros

# graphAbgp imports
from graph_organism import OrganismGraph
from subclass_sitealignment import BasalSiteAlignmentFunctions, sort_by_cumulative_score
from subclass_pssmobjects import BasalPSSMObjectGraphFunctions
import conversion
import ordering
import graphPlus
from exceptions import *


# Global Varibale Imports
from settings.translationalstartsites import TCODE_TSS_5P_WINDOW, TCODE_TSS_3P_WINDOW

class PssmObjectCollectionGraph(OrganismGraph,BasalSiteAlignmentFunctions,BasalPSSMObjectGraphFunctions):
    """
    """
    def __init__(self):
        """
        Initialize a PssmObjectCollectionGraph; objects on orfs on pacbporfs (splice sites, start sites, stop sites)
        """
        OrganismGraph.__init__(self)
        self._node_pssm     = {}
        self._node_object   = {}
        self.MIN_PSSM_SCORE = None
        # list to store aligned sites graphs
        self.alignedsites = []
        # dictionary keeping track of considered range
        self._organism_consideredrange = {}

        # Attributes storing information for thresholds of 
        # how this collection was generated
        self.ALIGNED_SITE_AA_OFFSET = None
        self.MIN_PSSM_SCORE = None

    # end of function __init__


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

        @type  node: *
        @param node: One Node

        @rtype:  *
        @return: organism identifier (presumably a string)
        """
        return node[0]

    # end of function _organism_from_node 


    def organism_by_node(self,node):
        """
        @attention: alias for _organism_from_node()
        """
        return self._organism_from_node(node)

    # end of function organism_by_node

    def __str__(self,formatted_nodestring=""):
        """ Nicely formatted oneliner reprenting this object """
        # get formatted nodestring if not applied as argument
        if not formatted_nodestring:
            formatted_nodestring = _format_nodes_to_txt(self)

        # some inheriting objects have a phase; not all!
        try:    phases = self.phase()
        except: phases = ""

        # format objsettings; if None, float formatting can fail
        try:
            objsettings = "%1.1f d%s" % ( 
                    float(self.MIN_PSSM_SCORE),
                    int(self.ALIGNED_SITE_AA_OFFSET)
                    )
        except:
            objsettings = "%s d%s" % (
                    self.MIN_PSSM_SCORE,
                    self.ALIGNED_SITE_AA_OFFSET
                    )

        # return a nicely formatted string
        return "<%s (%s) N%s-O%s-E%s %s (csat: %1.1f tw:%2.1f) %s >" % (
                self.__class__.__name__,
                objsettings,
                self.node_count(),
                self.organism_set_size(),
                self.edge_count(),
                phases,
                self.connectivitysaturation(),
                self.total_weight(),
                formatted_nodestring
                )

    # end of function __str__


    def cumulative_score(self):
        """
        Return None because a Collection of sites has no cumulative_score!!
        """
        return None

    # end of function cumulative_score


    def find_fully_connected_subgraphs(self,edges=None,max_missing_edges=2):
        """
        See subclass_sitealignment.BasalSiteAlignmentFunctions.find_fully_connected_subgraphs
        """
        return BasalSiteAlignmentFunctions.find_fully_connected_subgraphs(
                    self,edges=edges,
                    max_missing_edges=max_missing_edges
                    )

    # end of find_fully_connected_subgraphs


    def find_perfectly_conserved_sites(self,edges=None):
        """
        See subclass_sitealignment.BasalSiteAlignmentFunctions.find_perfectly_conserved_sites
        """
        return BasalSiteAlignmentFunctions.find_perfectly_conserved_sites(
                    self,edges=edges
                    )

    # end of function find_perfectly_conserved_sites


    def find_conserved_sites(self,edges=None):
        """
        See subclass_sitealignment.BasalSiteAlignmentFunctions.find_conserved_sites
        """
        # if edges is not applied get by definition form the OrganismGraph
        if not edges: edges = self.organism_set_size() - 1

        return BasalSiteAlignmentFunctions.find_conserved_sites(
                    self,edges=edges
                    )

    # end of function find_conserved_sites


    def align_all_sites(self,edges=None,minimal_edges=2):
        """
        See subclass_sitealignment.BasalSiteAlignmentFunctions.align_all_sites
        """
        return BasalSiteAlignmentFunctions.align_all_sites(
                    self,edges=edges,
                    minimal_edges=minimal_edges
                    )

    # end of function align_all_sites


    def collection2alignedsites(self,edges=None,minimal_edges=2):
        """
        @attention: MUST BE OVERRIDDEN IN THE SUBCLASSES, because action depends on what type of site is to be aligned
        """
        pass

    # end of function collection2alignedsites


    def get_optimal_single_site(self,organism,bestalignedsite=None):
        """
        Get the optimal site from the collection

        @attention: MUST BE OVERRIDDEN IN THE SUBCLASSES, if the returned value is not correct

        @type  organism: *
        @param organism: Organism identifier

        @rtype:  *
        @return: optimal PSSM site
        """
        return None

    # end of function get_optimal_single_site 


    def togff(self,organism=None,gff={}):
        """
        Create gff tuple for ALL objects in the collection of a specific organism

        @type  gff: dictionary
        @param gff: overwrite default gff data, keys: ('fstrand','fphase','fref',etc...)

        @type  organism: * (presumably string)
        @param organism: Organism identifier to make the gff for

        @rtype:  list of tuples
        @return: list of gff tuples with 9 elements
        """
        gffdata = []
        if organism in self.organism_set():
            for object in ordering.order_list_by_attribute( self.get_organism_objects(organism), order_by='pos'):
               gffdata.append( object.togff( gff=gff ) )
        return gffdata

    # end of function togff

# end of class PssmObjectCollectionGraph


class TranslationalStartSiteCollectionGraph(PssmObjectCollectionGraph):
    """
    """
    def __init__(self,tcode_5p_windowsize=TCODE_TSS_5P_WINDOW,tcode_3p_windowsize=TCODE_TSS_3P_WINDOW):
        """
		Initialize a TranslationalStartSiteCollectionGraph
        """
        PssmObjectCollectionGraph.__init__(self)
        self._tcode5pscore = {}
        self._tcode3pscore = {}
        self._TCODE_5P_WINDOWSIZE = tcode_5p_windowsize
        self._TCODE_3P_WINDOWSIZE = tcode_3p_windowsize

    # end of function __init__


    def _update_after_changes(self):
        """
        @attention: MUST BE OVERRIDDEN IN THE SUBCLASSES, if the returned value is not correct
        """
        PssmObjectCollectionGraph._update_after_changes(self)
        # and update the tcode data
        self._update_tcode_data()

    # end of function _update_after_changes


    def _update_tcode_data(self):
        """
        Remove edges from self._tcode5pscore and self._tcode3pscore that
        are no longer present in self.weights (after deletion of nodes/edges)

        @attention: USE THIS FUNCTION AFTER DELETION OF NODES FROM THE GRAPH
        """
        keys_nodes = self.get_nodes()
        # check if tcode5p key is present in weights
        keys_tcode5p  = self._tcode5pscore.keys()
        for ktc in keys_tcode5p:
            if ktc not in keys_nodes:
                del( self._tcode5pscore[ktc] )
        # check if tcode3p key is present in weights
        keys_tcode3p  = self._tcode3pscore.keys()
        for ktc in keys_tcode3p:
            if ktc not in keys_nodes:
                del( self._tcode3pscore[ktc] )

    # end of function _update_tcode_data


    def collection2alignedsites(self,edges=None,minimal_edges=2):
        """
        TODO!!!
        """
        # handle edges argument (error check etc.)
        edges = self._handle_edges_argument(edges)

        # do the first basal alignment
        self.alignedsites = sort_by_cumulative_score(
                [ conversion.TranslationalStartSiteCollectionGraph2AlignedTranslationalStartSiteGraph(algsite,max_node_count=self.organism_set_size()) for algsite in self.find_conserved_sites(edges=edges) ]
                )
        # and order all the sites on cumulative score
        self.alignedsites = sort_by_cumulative_score(self.alignedsites)

        # no keep on aligning the remaining fraction of non-aligned sites
        for current_edges in range(edges-1,minimal_edges-1,-1):
            if self.alignedsites and self.alignedsites[-1].__class__.__name__ == 'TranslationalStartSiteCollectionGraph':
                # redo this non-aligned part
                gra = self.alignedsites.pop()
                self.alignedsites.extend(
                        [ conversion.TranslationalStartSiteCollectionGraph2AlignedTranslationalStartSiteGraph(algsite,max_node_count=self.organism_set_size()) for algsite in gra.find_conserved_sites(edges=current_edges) ]
                        )
                self.alignedsites = sort_by_cumulative_score(self.alignedsites)

    # end of function collection2alignedsites


    def get_optimal_single_site(self,organism,bestalignedsite=None):
        """
        Get the optimal TranslationalStartSite from the collection

        @type  organism: *
        @param organism: Organism identifier

        @type  bestalignedsite: AlignedTranslationalStartSiteGraph 
        @param bestalignedsite: AlignedTranslationalStartSiteGraph instance or None

        @rtype:  TranslationalStartSite
        @return: optimal TranslationalStartSite instance 
        """
        # get the most frontal (5p) TSS; phase check is not needed (all phases==0)
        sites = ordering.order_list_by_attribute( self.get_organism_objects(organism), order_by='pos')

        # TODO: now just the most frontal site is chosen; not the site most nearby the OMSR,
        # or the most frontal site > a certain score. This might be implemented after studying
        # ample examples
        if bestalignedsite:
            if sites:
                return sites[0]
            else:
                return None
        else:
            if sites:
                return sites[0]
            else:
                return None

    # end of function get_optimal_single_site

# end of class TranslationalStartSiteCollectionGraph


def _format_nodes_to_shortened_txt(graphobj):
    """
    """
    nodes = graphobj.get_nodes()
    nodes.sort()
    # try some string shortening here
    try:    nodes = [ '%s_%s' % (node[0],node[-1]) for node in nodes ]
    except: pass
    return " ".join([ str(node) for node in nodes ])

# end of function _format_nodes_to_shortened_txt


def _format_nodes_to_txt(graphobj):
    """
    """
    nodestxt = []
    organisms = list(graphobj.organism_set())
    organisms.sort()
    for org in organisms:
        orgtxt = [ "%s:" % org ]
        nodes = graphobj.get_organism_nodes(org)
        nodes.sort()
        for node in nodes:
             orgtxt.append( "%s" % (node[-1]) )
        nodestxt.extend( orgtxt )
    # return text string
    return " ".join([ str(item) for item in nodestxt ])

# end of function _format_nodes_to_txt



def _format_splicesite_nodes_to_txt(graphobj):
    """
    """
    nodestxt = []
    organisms = list(graphobj.organism_set())
    organisms.sort()
    for org in organisms:
        orgtxt = [ "%s:" % org ]
        for node in graphobj.get_organism_nodes(org):
             obj = graphobj.get_node_object(node)
             symbol = "?"
             if obj.__class__.__name__[0:9] == 'Projected':
                 symbol = 'P'
             elif obj.__class__.__name__ in ['CodingBlockStart','CodingBlockEnd']:
                 symbol = 'X'
             elif graphobj._type == 'Donor':
                 if obj.is_canonical(): symbol = 'D'
                 else:                  symbol = 'N'
             elif graphobj._type == 'Acceptor':
                 symbol = 'A'
             else:
                 pass
             # make obj.phase string and then first letter (None->N)
             orgtxt.append( "%s%s%s" % (node[-1],symbol,str(obj.phase)[0] ) )
        nodestxt.extend( orgtxt )
    # return text string
    return " ".join([ str(item) for item in nodestxt ])

# end of function _format_splicesite_nodes_to_txt 


class SpliceSiteCollectionGraph(PssmObjectCollectionGraph):
    """
    """
    def __init__(self):
        """
        Initialize a SpliceSiteCollectionGraph
        """
        PssmObjectCollectionGraph.__init__(self)
        # which type? donors or acceptors?
        self._type          = None    # to be overwritten in the subclasses

    # end of function __init__


    def add_edge(self, u, v, wt=1, entropy=1.0):
        """
        """
        PssmObjectCollectionGraph.add_edge(self,u,v,wt=wt)
        self._edge_binary_entropies[(u,v)] = (entropy,entropy)
        self._edge_binary_entropies[(v,u)] = (entropy,entropy)
 
    # end of function add_edge


    def __str__(self):
        """ Nicely formatted oneliner reprenting this object """

        return PssmObjectCollectionGraph.__str__(
               self,
               formatted_nodestring=\
               _format_splicesite_nodes_to_txt(self)
               ) 

        ## get list of present phases
        #phases = self.phase()
        ## and return a nicely formatted string
        #return "<%s N%s-O%s-E%s %s (csat: %1.1f tw:%2.1f) %s >" % (
        #        self.__class__.__name__,
        #        self.node_count(),
        #        self.organism_set_size(),
        #        self.edge_count(),
        #        phases,
        #        self.connectivitysaturation(),
        #        self.total_weight(),
        #        _format_splicesite_nodes_to_txt(self),
        #        )

    # end of function __str__


    def omsr_entropy_of_position(self,site,omsraacoord,window=15):
        """
        @attention: MUST BE OVERWRITTEN IN Donor- and Acceptor Collection graphs!
        """
        pass

    # end of function omsr_entropy_of_position


    def set_consideredsplicesiterange(self,organism,start,end):
        """
        Set the considered splice site range of this organism (on its orf)

        @type  organism: * (presumably string)
        @param organism: Organism identifier to return the splice site range for

        @type  start: positive integer
        @param start: start coordinate (on current orf) from where splice sites are collected

        @type  end: positive integer
        @param end: end coordinate (on current orf) untill where splice sites are collected
        """
        self._organism_consideredrange[organism] = (start,end)

    # end of function set_consideredsplicesiterange


    def get_consideredsplicesiterange(self,organism):
        """
        Get the considered splice site range of this organism (on its orf)

        @type  organism: * (presumably string)
        @param organism: Organism identifier to return the splice site range for

        @rtype:  tuple
        @return: tuple with (python) start and end coordinates of the nt position on the orf
        """
        if self._organism_consideredrange:
            if self._organism_consideredrange.has_key(organism):
                return self._organism_consideredrange[organism]
            elif organism not in self.organism_set():
                print "OrganismNotPresentInGraph:%s" % organism
                print self._organism_consideredrange
                raise OrganismNotPresentInGraph
            else:
                return ( None, None )
        else:
            return ( None, None )

    # end of function get_consideredsplicesiterange


    def phase(self):
        """
        Phase function of SpliceSiteCollectionGraph; return unique list of phases or single phase
        """
        if len( Set([ o.phase for o in self._node_object.values() ]) ) != 1:
            return list( Set( [ o.phase for o in self._node_object.values() ] ) )
        elif len(self._node_object) == 0:
            # no objects in graph (yet)
            return None
        else:
            return self._node_object.values()[0].phase

    # end of function phase


    def collection2alignedsites(self,edges=None,minimal_edges=2):
        """
        """
        # handle edges argument (error check etc.)
        edges = self._handle_edges_argument(edges)

        # do the first basal alignment
        self.alignedsites = sort_by_cumulative_score(
                [ conversion.SpliceSiteCollectionGraph2AlignedSpliceSiteGraph(algsite,max_node_count=self.organism_set_size()) for algsite in self.find_conserved_sites(edges=edges) ]
                )

        # now keep on aligning the remaining fraction of non-aligned sites
        for current_edges in range(edges-1,minimal_edges-1,-1):
            if self.alignedsites and self.alignedsites[-1].__class__.__name__ in\
            ['DonorSiteCollectionGraph','AcceptorSiteCollectionGraph','SpliceSiteCollectionGraph']:
                # redo this non-aligned part
                gra = self.alignedsites.pop()
                self.alignedsites.extend(
                        [ conversion.SpliceSiteCollectionGraph2AlignedSpliceSiteGraph(algsite,max_node_count=self.organism_set_size()) for algsite in gra.find_conserved_sites(edges=current_edges) ]
                        )
                self.alignedsites = sort_by_cumulative_score(self.alignedsites)

        # if there are AlignedSpliceSiteWithPhaseShiftGraphs, remove all that have not all organisms represented
        lenrange = range(len(self.alignedsites)-1,-1,-1) 
        for pos in lenrange: 
            if self.alignedsites[pos].__class__.__name__ in ['AlignedDonorSiteWithPhaseShiftGraph',
            'AlignedAcceptorSiteWithPhaseShiftGraph','AlignedSpliceSiteWithPhaseShiftGraph']:
                if self.alignedsites[pos].organism_set_size() != self.organism_set_size():
                    _removed = self.alignedsites.pop(pos)

        # finally, merge AlignedSites that are separated due possible erroneous alignments
        # around (aligned) inframe-introns. Due to ALIGNED_DONOR_MAX_TRIPLET_DISTANCE,
        # situations like the following can occur. Suppose organism A-E, E having an inframe intron
        # that is in some Pacbps splitted, and in some aligned. Due to the differences in the location
        # where BLAST places the gaps, offset can arrise, resulting in 2 i.s.o. 1 AlignedSite:
        # A(x)-B(x)-C(x)-D(x) and A(x)-B(x)-C(x)-E(x), where A(x)-B(x)-C(x) are the same sites!
        if self.alignedsites:
            currentpos = 0
            while True:
                for pos in range(currentpos,len(self.alignedsites)):
                    site_merged = False
                    site = self.alignedsites[pos]
                    if site.organism_set_size() == self.organism_set_size():
                        continue
                    if site.__class__.__name__ not in\
                    ['AlignedDonorSiteGraph','AlignedAcceptorSiteGraph','AlignedSpliceSiteGraph']:
                        continue
                    for otherpos in range(pos+1,len(self.alignedsites)):
                        othersite = self.alignedsites[otherpos]
                        if othersite.organism_set_size() == self.organism_set_size():
                            continue
                        if othersite.__class__.__name__ not in\
                        ['AlignedDonorSiteGraph','AlignedAcceptorSiteGraph','AlignedSpliceSiteGraph']:
                            continue
                        if site.phase() == othersite.phase():
                            mutual_nodes = graphPlus.comparison.mutual_nodes(site,othersite)
                            if not mutual_nodes: continue
                            # check if the difference in nodes completes `site`
                            new_nodes = othersite.node_set().difference( site.get_nodes() )
                            new_orgs  = [ othersite._organism_from_node(node) for node in new_nodes]
                            if not site.organism_set().intersection(new_orgs):
                                # yes, there are mutual nodes! Now update `site` with the nodes
                                # that are only occuring in `othersite`
                                site.nodes.update(othersite.nodes)
                                site.weights.update(othersite.weights)
                                site._node_pssm.update(othersite._node_pssm)
                                site._node_object.update(othersite._node_object)
                                site._edge_binary_entropies.update(othersite._edge_binary_entropies)
                                # set site_merged to True to make shure the outern for loop
                                # is broken, remove `othersite` and break out
                                site_merged = True
                                self.alignedsites.pop(otherpos)
                                break
                    if site_merged:
                        # yes, a succesfull merge; break the outern forloop because the
                        # length of the list self.alignedsites has changed!
                        currentpos = pos
                        break
                else:
                    # eof list; break the while loop
                    break
            # if currentpos>0, >=1 sites are merged -> resort the sites
            self.alignedsites = sort_by_cumulative_score(self.alignedsites)

    # end of function collection2alignedsites


    def get_optimal_single_site(self,organism,phase=None,bestalignedsite=None,use_entropy_around_omsr=None):
        """
        Get the highest PSSM scoring splice site from the collection

        @type  organism: *
        @param organism: Organism identifier

        @type  phase: 0,1,2 or None 
        @param phase: specify phase to get the highest scoring site of 

        @type  bestalignedsite: AlignedSpliceSiteGraph (Donor,Acceptor...)
        @param bestalignedsite: AlignedSpliceSiteGraph or inheriting object

        @type  use_entropy_around_omsr: integer
        @param use_entropy_around_omsr: 

        @rtype:  * (or None)
        @return: optimal PSSM site
        """

        # TODO: take splice site phase shift into account!
        # TODO: take splice site phase shift into account!
        # TODO: take splice site phase shift into account!


        if phase in [0,1,2]:
             # get best object by phase
             if organism in self.organism_set():
                elegiable_sites = []
                for site in self.get_organism_objects(organism):
                    if site.phase != phase: continue
                    if use_entropy_around_omsr:
                        entropyscore = self.omsr_entropy_of_position(site,use_entropy_around_omsr)
                        elegiable_sites.append( ( entropyscore*site.pssm_score, site ) )
                    else:
                        elegiable_sites.append( ( site.pssm_score, site ) )

                if elegiable_sites:
                    elegiable_sites.sort()
                    return elegiable_sites[-1][1]
                else:
                    # no sites of this phase for this organism
                    return None
             else:
                 # no sites at all for this organism
                 return None
        elif bestalignedsite:
            # get the best site assuming the aligned site is the optimal one
            # use only when the organism of interest is ABSENT in the best aligned site
            for aligned in self.alignedsites[1:]:
                if bestalignedsite.phase() == aligned.phase() and organism in aligned.organism_set():
                    return aligned.get_organism_objects(organism)[0]
            else:
                # check the unaligned remaining collection for sites of the proper phase
                elegiable_sites = []
                possible_phases = self.alignedsites[-1].phase()
                if type(possible_phases) == type(int()):
                    possible_phases = [ possible_phases ]
                if bestalignedsite.phase() in possible_phases:
                    if organism in self.alignedsites[-1].organism_set():
                        for site in self.alignedsites[-1].get_organism_objects(organism):
                            if site.phase == bestalignedsite.phase():
                                elegiable_sites.append( ( site.pssm_score, site ) )
                if elegiable_sites:
                    elegiable_sites.sort()
                    return elegiable_sites[-1][1]

                # if we reach this point, no site can be found in accordance
                # with the bestalignedsite. Just get best site without this knowledge
                return self.get_optimal_single_site(organism)
        else:
            if organism in self.organism_set():
                elegiable_sites = []
                for site in self.get_organism_objects(organism):
                    if use_entropy_around_omsr:
                        entropyscore = self.omsr_entropy_of_position(site,use_entropy_around_omsr)
                        elegiable_sites.append( ( entropyscore*site.pssm_score, site ) )
                    else:
                        elegiable_sites.append( ( site.pssm_score, site ) )

                if elegiable_sites:
                    elegiable_sites.sort()
                    return elegiable_sites[-1][1]
                else:
                    return None
            else:
                return None

    # end of function get_optimal_single_site 


# end of class SpliceSiteCollectionGraph


class DonorSiteCollectionGraph(SpliceSiteCollectionGraph):
    """
    DonorSiteCollectionGraph (DSCOLG) class, inheriting from SpliceSiteCollectionGraph class.
    """
    def __init__(self):
        """
        Initialize as a SpliceSiteCollectionGraph
        """
        SpliceSiteCollectionGraph.__init__(self)
        # which type donors or acceptors?
        self._type = "Donor"

    # end of function __init__


    def omsr_entropy_of_position(self,site,omsraacoord,window=15):
        """
        """
        # get the binary entropy array
        bea = zeros(window*2)
        try:    coord = site.pos / 3
        except: return 0.0
        omsrinbeacoord = max([ 0, window+(coord-omsraacoord) ])
        # set omsr area to 1 in bea array
        bea[0:omsrinbeacoord] = 1
        # calculate maximum entropy score
        max_entropy = float(sum(range(1,window+1)))

        # now calculate entropy of position around omsr
        score_left  = sum(range(1,sum(bea[0:window])+1))
        score_rigth = sum(range(1,sum(bea[window:window*2])+1))
        entropy     = score_left - score_rigth
        entropy     = float(entropy) / max_entropy

        # return entropy value
        return entropy

    # end of function omsr_entropy_of_position

# end of class DonorSiteCollectionGraph


class AcceptorSiteCollectionGraph(SpliceSiteCollectionGraph):
    """
    AcceptorSiteCollectionGraph (ASCOLG) class, inheriting from SpliceSiteCollectionGraph class.
    """
    def __init__(self):
        """
        Initialize as a SpliceSiteCollectionGraph
        """
        SpliceSiteCollectionGraph.__init__(self)
        # which type donors or acceptors?
        self._type = "Acceptor"

    # end of function __init__


    def omsr_entropy_of_position(self,site,omsraacoord,window=15):
        """
        """
        # get the binary entropy array
        bea = zeros((window*2))
        try:    coord = site.pos / 3
        except: return 0.0
        omsrinbeacoord = max([ 0, window+(omsraacoord-coord) ])
        # set omsr area to 1 in bea array
        bea[omsrinbeacoord:window*2] = 1
        # calculate maximum entropy score
        max_entropy = float(sum(range(1,window+1)))

        # now calculate entropy of position around omsr
        score_left  = sum(range(1,sum(bea[0:window])+1))
        score_rigth = sum(range(1,sum(bea[window:window*2])+1))
        entropy     = score_rigth - score_left
        entropy     = float(entropy) / max_entropy
        # return entropy value
        return entropy

    # end of function omsr_entropy_of_position

# end of class AcceptorSiteCollectionGraph


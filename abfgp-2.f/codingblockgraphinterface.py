"""
The CodingBlockGraphInterface class described the (intronic) interface
between two adjacent CodingBlockGraphs in a GenestructureOfCodingBlockGraphs
instance.
"""

# Module metadata
__authors__ = "Ate van der Burgt"
__license__ = "MIT"

# ABGP Imports
from graphAbgp.subclass_sitealignment import sort_by_cumulative_score
from graphAbgp.ordering import order_list_by_attribute
from graphAbgp import conversion
from graphAbgp.exceptions import OrganismNotPresentInGraph

# Other imports
from lib_stopwatch import StopWatch
import graphPlus
from intronprojection import (
    project_missing_introns,
    projectedintrons2projectedsites,
    )

# Gene Imports
from gene.codingblock import CodingBlockStart, CodingBlockEnd
from gene.donor import SpliceDonor, ProjectedSpliceDonor
from gene.acceptor import SpliceAcceptor, ProjectedSpliceAcceptor

# Python imports
from sets import Set

# Global variables
from settings.gff.cbginterface import *
from settings.inframeintron import (
    INFRAME_INTRON_MIN_AA_LENGTH,
    )
from settings.genestructure import (
    MIN_INTRON_NT_LENGTH
    )
from settings.sitealignment import (
    ALIGNED_DONOR_MAX_TRIPLET_DISTANCE,
    ALIGNED_ACCEPTOR_MAX_TRIPLET_DISTANCE,
    PROPOSE_IMPROVED_DONOR_MAX_AA_DIST,
    PROPOSE_IMPROVED_DONOR_MIN_PSSM_RATIO,
    PROPOSE_IMPROVED_ACCEPTOR_MAX_AA_DIST,
    PROPOSE_IMPROVED_ACCEPTOR_MIN_PSSM_RATIO,
    # variables in the cbgIF.optimize() function
    OPTIMIZE_CBGIF_DONOR_MAX_TRIPLET_DISTANCE,
    OPTIMIZE_CBGIF_DONOR_ENLARGE_3P_BOUNDARY_BY,
    OPTIMIZE_CBGIF_ACCEPTOR_MAX_TRIPLET_DISTANCE,
    OPTIMIZE_CBGIF_ACCEPTOR_ENLARGE_5P_BOUNDARY_BY,
    # variables in the cbgIF.optimizetinyexon() function
    OPTIMIZE_CBGIF_TINYEXON_MAX_AA_SIZE,
    OPTIMIZE_CBGIF_TINYEXON_DONOR_OMSR_3P_NT_OFFSET,
    OPTIMIZE_CBGIF_TINYEXON_ACCEPTOR_OMSR_5P_NT_OFFSET,
    )



class CodingBlockGraphInterface:
    """
    CodingBlockGraphInterface, describes the introns at the interfac of 2 CodingBlockGraphs
    
    @attention: the cbgIF is stored in the CBG as: cbgD._CBGinterface3p = cbgIF
    @attention: the cbgIF is stored in the CBG as: cbgA._CBGinterface5p = cbgIF
    """
    def __init__(self,donorCBG,acceptorCBG):
        """
        Initialize a CodingBlockGraphInterface between two adjacent CBGs
        
        @type  donorCBG: CodingBlockGraph
        @param donorCBG: 5' side CodingBlockGraph of the interface (donor)
        
        @type  acceptorCBG: CodingBlockGraph
        @param acceptorCBG: 3' side CodingBlockGraph of the interface (acceptor)
        """
        self.donorCBG                   = donorCBG 
        self.acceptorCBG                = acceptorCBG
        self.verbose                    = False
        self._HAS_INTRONS_PROJECTED     = False
        self._USE_ENTROPY               = True  
        self._interface_is_intron       = {}
        
        # splitted/lsrCBG ends for the DONOR CBG
        self._forced_3p_ends            = {}
        # splitted/lsrCBG ends for the ACCEPTOR CBG
        self._forced_5p_ends            = {}
        
        # dict attributes for storing the projected splice sites
        self._projected_acceptor_sites  = {}
        self._projected_donor_sites     = {}
        
        # attributes for splice site collection graph settings
        self._donor_allow_phase_shift      = False
        self._donor_allow_non_canonical    = False
        self._acceptor_allow_phase_shift   = False
        self._acceptor_allow_non_canonical = False
        
        # attributes that contains the optimal AlignedSpliceSiteGraphs
        self._optimal_aligned_donor     = None
        self._optimal_aligned_acceptor  = None
        
        # attributes that contains the SpliceSiteGraphs
        self._spliceacceptorgraph       = None
        self._splicedonorgraph          = None
        
        # list that stores optimization actions that have been done
        self._optimizations             = []
        
        # obtain enforced CBG ends in case of an lsrCBG interface
        self._obtain_forced_ends()
        
        # obtain interface_is_intron data into cache
        for organism in self.organism_set():
            self._interface_is_intron[organism] =\
                    self.interface_is_intron(organism)

    # end of function __init__


    def __str__(self):
        """ """
        dorgsetsize, dnodes,dedges,doffset = 0,0,0, "N"
        aorgsetsize, anodes,aedges,aoffset = 0,0,0, "N"
        try:
            # get properties of the SpliceDonorGraph
            doffset     = self._splicedonorgraph.ALIGNED_SITE_AA_OFFSET 
            dnodes      = self._splicedonorgraph.node_count()
            dedges      = self._splicedonorgraph.edge_count()
            dorgsetsize = self._splicedonorgraph.organism_set_size()
        except:
            pass
        try:
            # get properties of the SpliceAcceptorGraph
            aoffset     = self._spliceacceptorgraph.ALIGNED_SITE_AA_OFFSET
            anodes      = self._spliceacceptorgraph.node_count()
            aedges      = self._spliceacceptorgraph.edge_count()
            aorgsetsize = self._spliceacceptorgraph.organism_set_size()
        except:
            pass
        
        # create substrings of the cbgIF
        dsettings = "%s(%s%s%s)" % (
            self.donorCBG._short_name,
            doffset,
            str(self._donor_allow_phase_shift)[0],
            str(self._donor_allow_non_canonical)[0], 
            )
        asettings = "%s(%s%s%s)" % (
            self.acceptorCBG._short_name,
            aoffset,
            str(self._acceptor_allow_phase_shift)[0],
            str(self._acceptor_allow_non_canonical)[0],
            )
        dgrainfo = "D: %s,%s,%s" % ( dorgsetsize, dnodes, dedges )
        agrainfo = "A: %s,%s,%s" % ( aorgsetsize, anodes, aedges )
        mutnodes = []
        try:
            for (org,orfid) in graphPlus.comparison.mutual_nodes(
            self.donorCBG,self.acceptorCBG):
                mutnodes.append( "%s:%s%s" % (
                    org,orfid,str(self._interface_is_intron[org])[0] )
                    )
        except:
            pass

        # return string representation of the cbgIF
        return "<CBGinterface %s %s %s [%s,%s] mutual:(%s) %s %s >" % (
                self.is_compatible(),
                self.is_optimal_donor(),
                self.is_optimal_acceptor(),
                dsettings,
                asettings,
                " ".join(mutnodes),
                dgrainfo,
                agrainfo,
                )

    # end of function __str__


    def interfaceproperties(self,detailed=True):
        """
        Print a detailed report on the properties of this cbgIF

        @type  detailed: Boolean
        @param detailed: extra information is printed when detailed is True
        """
        print self._interfaceproperties(detailed=detailed)

    # end of function interfaceproperties


    def _interfaceproperties(self,detailed=True):
        """
        Return detailed report on the properties of this cbgIF

        @type  detailed: Boolean
        @param detailed: extra information is returned when detailed is True

        @rtype:  string
        @return: multi-line, somewhat layouted overview of cbgIF properties
        """
        if self.donorCBG._short_name != "CBG":
            return "" 
        if self.acceptorCBG._short_name != "CBG":
            return "" 

        # gather distance per organism between codingblocks
        aadistomsr  = self.donorCBG.distance_between_codingblocks(
                self.acceptorCBG)
        aadistmaxsr = self.donorCBG.maxsr_distance_between_codingblocks(
                self.acceptorCBG)
        aadistorfD  = self.donorCBG.orfenddistance('right')
        aadistorfA  = self.acceptorCBG.orfenddistance('left')
        omsrD       = self.donorCBG.overall_minimal_spanning_range()
        omsrA       = self.acceptorCBG.overall_minimal_spanning_range()
        maxsrD      = self.donorCBG.maximal_spanning_range()
        maxsrA      = self.acceptorCBG.maximal_spanning_range()

        # gather mutual nodes & organisms
        mutual_nodes = Set(self.donorCBG.get_nodes()).intersection(
                self.acceptorCBG.get_nodes())
        mutual_orgs  = [ self.donorCBG.organism_by_node(node) for\
                node in mutual_nodes ]

        # info line for each organism identifier is appended to list of lines
        lines = []

        # loop over the organism identifiers and gather info
        for org in self.organism_set():
            # if aadist* dict is missing keys, create None values
            # this occurs when cbgD or cbgA is a K(s-x) CBG
            if not aadistomsr.has_key(org):  aadistomsr[org]  = None
            if not aadistmaxsr.has_key(org): aadistmaxsr[org] = None
            if not aadistorfD.has_key(org):  aadistorfD[org]  = None
            if not aadistorfA.has_key(org):  aadistorfA[org]  = None
            try:
                donorrange =\
                    self._splicedonorgraph.get_consideredsplicesiterange(org)
            except:
                donorrange = (None,None)
            try:
                acceptorrange =\
                    self._spliceacceptorgraph.get_consideredsplicesiterange(org)
            except:
                acceptorrange = (None,None)

            # gather info concerning the acceptor CBG
            if org in self.acceptorCBG.organism_set():
                accepNode = self.acceptorCBG.node_by_organism(org)
                acceptxt  = "orf:%s<-%s-%s:omsr          " % (
                    aadistorfA[org],
                    min(maxsrA[accepNode]),
                    min(omsrA[accepNode])
                    )
                acceptxt  = acceptxt[0:24]
            else:
                accepNode = None
                acceptxt  = " "*24

            # gather info concerning the acceptor CBG
            if org in self.donorCBG.organism_set():
                donorNode = self.donorCBG.node_by_organism(org)
                donortxt  = "omsr:%s-%s->%s:orf          " % (
                    max(omsrD[donorNode]),
                    max(maxsrD[donorNode]),
                    aadistorfD[org]
                    )
                donortxt  = donortxt[0:24]
            else:
                donorNode = None
                donortxt  = " "*24

            # gather info concerning the distance between donor & acceptor CBG
            if donorNode and accepNode:
                disttxt = "dist(omsr:%s,maxsr:%s)     " % (
                    aadistomsr[org],
                    aadistmaxsr[org]
                    )
                disttxt = disttxt[0:24]
            else:
                disttxt = " "*24
            if org in mutual_orgs: mutual = "M"
            else:                  mutual = " "

            # make properties line of the cbgIF for this Organism indentifier
            line =  "%s %s\t%s \t%s\t %s\t%s\t%s\t%s" % (
                mutual, org, donortxt, disttxt,
                acceptxt, org, donorrange, acceptorrange )
            # append to list of lines
            lines.append(line)

        # if detailed => more data!
        if detailed:
            lines.append( str( self._splicedonorgraph ) )
            lines.append( str( self._optimal_aligned_donor ) )
            if self._splicedonorgraph:
                for org in self._splicedonorgraph.organism_set(): 
                    line = "BEST-d: %s %s optphase: %s intron?: %s improvable?: %s %s" % (
                        org, self._splicedonorgraph.get_optimal_single_site(org),
                        self.highest_pssm_phase(organism=org),
                        self._interface_is_intron[org],
                        self.is_donor_improveable(org),
                        self.propose_improved_donor(org)
                        )
                    lines.append( line )
            else:
                lines.append( "# NO splicedonorgraph in cbgIF !!" )
            if self._spliceacceptorgraph:
                for org in self._spliceacceptorgraph.organism_set():
                    line = "BEST-a: %s %s improvable?: %s %s" % (
                        org, self._spliceacceptorgraph.get_optimal_single_site(org),
                        self.is_acceptor_improveable(org),
                        self.propose_improved_acceptor(org)
                        )
                    lines.append( line )
            else:
                lines.append( "# NO spliceacceptorgraph in cbgIF !!" )
            lines.append( str( self._optimal_aligned_acceptor ) )
            lines.append( str( self._spliceacceptorgraph ) )

        # return one long txt string
        return "\n".join(lines)

    # end of function _interfaceproperties


    def togff(self,organism,FREF=None):
        """
        Create gff for this CodingBlockGraphInterface object
        """
        gffdata = []
        try:
            ( (dSta,dEnd), (aSta,aEnd) ) = self.consideredsplicesiterange()[organism] 
        except:
            print self.donorCBG
            print self.acceptorCBG
            print self._splicedonorgraph
            print self._spliceacceptorgraph
            print self
            ( (dSta,dEnd), (aSta,aEnd) ) = self.consideredsplicesiterange()[organism]

        # no splcie site range listed -> a lsrCBG interface!
        if (dSta,dEnd) == (None,None) or (aSta,aEnd) == (None,None):
            return []

        # make gff for the elegiable donor range that is taken into account
        if ( dSta, dEnd ) == ( None, None ):
            pass  # no splice site range set (yet)
        else:
            groupD = "%s %s-%s" % (GFF_DONORRANGE_GCLASS, dSta+1, dEnd)
            gffD = ( FREF, GFF_DONORRANGE_FSOURCE, GFF_DONORRANGE_FMETHOD,
                     dSta+1, dEnd, '.', '+', '.', groupD )
            gffdata.append(gffD)

        # make gff for the donor sites in this collection
        if GFF_DONOR_OUTPUT and self._splicedonorgraph:
            _gff_donor = {'fref': FREF, 'fsource': GFF_DONOR_FSOURCE,
                          'fmethod': GFF_DONOR_FMETHOD,
                          'gclass': GFF_DONOR_GCLASS  }
            gffdata.extend( self._splicedonorgraph.togff( organism=organism, gff=_gff_donor ) )

        # make gff for the AlignedDonorSiteGraphs
        if GFF_ALIGNED_DONOR_OUTPUT and self._splicedonorgraph:
            _gff_ads = {'fref': FREF, 'fsource': GFF_ALIGNED_DONOR_FSOURCE,
                        'fmethod': GFF_ALIGNED_DONOR_FMETHOD,
                        'gclass': GFF_ALIGNED_DONOR_GCLASS  }
            for adg in self._splicedonorgraph.alignedsites[0:GFF_ALIGNED_DONOR_REPORTBEST]:
                if adg.__class__.__name__ == 'AlignedDonorSiteGraph':
                    gffdata.append( adg.togff( organism=organism, gff=_gff_ads ) )

        # make gff for the elegiable acceptor range that is taken into account
        if ( aSta, aEnd ) == ( None, None ):
            pass  # no splice site range set (yet)
        else:
            groupA = "%s %s-%s" % (GFF_ACCEPTORRANGE_GCLASS, aSta+1, aEnd)
            gffA = ( FREF, GFF_ACCEPTORRANGE_FSOURCE, GFF_ACCEPTORRANGE_FMETHOD,
                     aSta+1, aEnd, '.', '+', '.', groupA )
            gffdata.append(gffA)

        # make gff for the acceptor sites in this collection
        if GFF_ACCEPTOR_OUTPUT and self._spliceacceptorgraph:
            _gff_acceptor = {'fref': FREF, 'fsource': GFF_ACCEPTOR_FSOURCE,
                             'fmethod': GFF_ACCEPTOR_FMETHOD,
                             'gclass': GFF_ACCEPTOR_GCLASS  }
            gffdata.extend( self._spliceacceptorgraph.togff( organism=organism, gff=_gff_acceptor ) )

        # make gff for the AlignedAcceptorSiteGraphs
        if GFF_ALIGNED_ACCEPTOR_OUTPUT and self._spliceacceptorgraph:
            _gff_aas = {'fref': FREF, 'fsource': GFF_ALIGNED_ACCEPTOR_FSOURCE,
                        'fmethod': GFF_ALIGNED_ACCEPTOR_FMETHOD,
                        'gclass': GFF_ALIGNED_ACCEPTOR_GCLASS  }
            for aag in self._spliceacceptorgraph.alignedsites[0:GFF_ALIGNED_ACCEPTOR_REPORTBEST]:
                if aag.__class__.__name__ == 'AlignedAcceptorSiteGraph':
                    gffdata.append( aag.togff( organism=organism, gff=_gff_aas ) )

        # make gff track for the cbgIF itself
        if ( dSta, aEnd ) == ( None, None ):
            pass  # no splice site range set (yet)
        else:
            # create c9tv_txt data
            c9tv_txt = [ ( '; %sObject' % self.__class__.__name__, str(self) ) ]
            c9tv_txt.append( ( "is_optimal", self.is_optimal() ) )
            c9tv_txt.append( ( "is_compatible", self.is_compatible() ) )
            c9tv_txt.append( ( "is_acceptable", self.is_acceptable() ) )
            c9tv_txt.append( ( "is_acceptable(organism)", self.is_acceptable(organism=organism) ) )
            for optimizationstep in self._optimizations:
                c9tv_txt.append( ( "optimization", optimizationstep ) )
            for org in self.organism_set():
                c9tv_txt.append( ( "interface_is_intron", "%s %s" % (org, self._interface_is_intron[org] ) ) )
 
            if self._splicedonorgraph:
                c9tv_txt.append( (GFF_DONORRANGE_FMETHOD, str( ( dSta+1, dEnd ) ) ) )
                c9tv_txt.append( (self._splicedonorgraph.__class__.__name__, str(self._splicedonorgraph) ) )
                if self._optimal_aligned_donor and self._optimal_aligned_donor.__class__.__name__ != 'DonorSiteCollectionGraph':
                    c9tv_txt.append( (self._optimal_aligned_donor.__class__.__name__, str(self._optimal_aligned_donor) ) )
                c9tv_txt.append( ("AlignedDonorSiteOffset", "%s AA" % self._splicedonorgraph.ALIGNED_SITE_AA_OFFSET ) )
                c9tv_txt.append( ( "is_optimal_donor", self.is_optimal_donor() ) )
                c9tv_txt.append( ( "is_optimal_donor(organism)", self.is_optimal_donor(organism=organism) ) )
                phases = {}
                for org in self._splicedonorgraph.organism_set():
                    c9tv_txt.append( ( "optimal_donor", "%s %s" % (org, self._splicedonorgraph.get_optimal_single_site(org) ) ) )
                    phases[org] = self.highest_pssm_phase(organism=org)
                c9tv_txt.append( ( "optimal_phases", str(phases) ) )
                c9tv_txt.append( ( "optimal_overall_phase", self.highest_pssm_phase() ) )


            if self._spliceacceptorgraph:
                c9tv_txt.append( (GFF_ACCEPTORRANGE_FMETHOD, str( ( aSta+1, aEnd ) ) ) )
                c9tv_txt.append( (self._spliceacceptorgraph.__class__.__name__, str(self._spliceacceptorgraph) ) )
                if self._optimal_aligned_acceptor and self._optimal_aligned_acceptor.__class__.__name__ != 'AcceptorSiteCollectionGraph':
                    c9tv_txt.append( (self._optimal_aligned_acceptor.__class__.__name__, str(self._optimal_aligned_acceptor) ) )
                c9tv_txt.append( ("AlignedAcceptorSiteOffset", "%s AA" % self._spliceacceptorgraph.ALIGNED_SITE_AA_OFFSET ) )
                c9tv_txt.append( ( "is_optimal_acceptor", self.is_optimal_acceptor() ) )
                c9tv_txt.append( ( "is_optimal_acceptor(organism)", self.is_optimal_acceptor(organism=organism) ) )
                for org in self._spliceacceptorgraph.organism_set():
                    c9tv_txt.append( ( "optimal_acceptor", "%s %s" % (org, self._spliceacceptorgraph.get_optimal_single_site(org) ) ) )

            # join c9tv_txt into a string
            c9tv_txt = "; ".join( [ "%s '%s'" % (a,b) for (a,b) in c9tv_txt ] )
            # replace '<' '>' symbols for '&lt;' and '&gt;' to enable proper html visualisation
            c9tv_txt = c9tv_txt.replace('<','[ ').replace('>',' ]')
            # obtain score;
            # is_optimal:    2
            # is_acceptable: 1
            # other:         0
            if self.is_optimal():                       score = 2
            elif self.is_acceptable(organism=organism): score = 1 
            else:                                       score = 0
            # create gclass/gname, gff tuple and insert into gffdata list
            groupIF = "%s %s-%s" % (GFF_CBGIF_GCLASS, dSta+1, aEnd)
            gffIF = ( FREF, GFF_CBGIF_FSOURCE, GFF_CBGIF_FMETHOD,
                     dSta+1, aEnd, score, '+', '.', "%s%s" % (groupIF,c9tv_txt) )
            gffdata.insert(0,gffIF)

        # return list with gff tuples
        return gffdata

    # end of function togff

    #############################################################
    #### Some functions that are accesed via one of the CBGs ####
    #############################################################

    def organism_by_node(self,node):
        """
        Get organism identifier by CBG node identifier

        @type  node: tuple
        @param node: Node identifier in (org,orfid) style

        @rtype:  string
        @return: Organism identifier
        """
        return self.donorCBG.organism_by_node(node)

    # end of function organism_by_node 


    def organism_set(self):
        """
        Get Set of organism identifiers from both cbgA and cbgD

        @rtype:  sets.Set
        @return: Set with organism identifiers
        """
        return self.donorCBG.organism_set().union(self.acceptorCBG.organism_set())

    # end of function organism_set


    def organism_set_size(self):
        """
        Get number of distinct organism identifiers in both cbgA and cbgD

        @rtype:  integer
        @return: size of self.organism_set()
        """
        return len(self.organism_set())

    # end of function organism_set_size

    #############################################################
    #### Functions to obtain **the** donor and acceptor      ####
    #############################################################

    def get_assigned_donor_site(self,org):
        """
        Get *the* assigned donor object of this Organism identifier

        @type  org: * (string)
        @param org: Organism Identifier

        @rtype  donor: mixed
        @return donor: None, DonorSite, ProjectedDonorSite or CodingBlockEnd
        """
        if org not in self.organism_set():
            raise OrganismNotPresentInGraph
        try:
            donor = self._optimal_aligned_donor.get_organism_objects(org)[0]
        except AttributeError:
            donor = None
        except OrganismNotPresentInGraph:
            donor = None
        # return donor site (or None)
        return donor

    # end of function get_assigned_donor_site


    def get_assigned_acceptor_site(self,org):
        """
        Get *the* assigned acceptor object of this Organism identifier

        @type  org: * (string)
        @param org: Organism Identifier

        @rtype  accep: mixed
        @return accep: None, AcceptorSite, ProjectedAcceptorSite or CodingBlockStart
        """
        if org not in self.organism_set():
            raise OrganismNotPresentInGraph
        try:
            accep = self._optimal_aligned_acceptor.get_organism_objects(org)[0]
        except AttributeError:
            accep = None
        except OrganismNotPresentInGraph:
            accep = None
        # return acceptor site (or None)
        return accep
    
    # end of function get_assigned_acceptor_site

    #############################################################
    #### Compatibility / Optimal check function (Booleans)   ####
    #############################################################

    def is_compatible(self):
        """
        Is this CodingBlockInterface a compatible one, i.e. assigned slice sites
        and sites of identical phase?

        @rtype:  Boolean
        @return: True or False
        """
        if self.donorCBG._short_name == "CBG" and\
        self.acceptorCBG._short_name == "lsrCBG":
            if self._forced_3p_ends and self._forced_5p_ends:
                return True
        elif self.donorCBG._short_name == "lsrCBG" and\
        self.acceptorCBG._short_name == "CBG":
            if self._forced_3p_ends and self._forced_5p_ends:
                return True
        elif self.donorCBG._short_name == "CBG" and\
        self.acceptorCBG._short_name == "CBG":
            if self._forced_3p_ends or self._forced_5p_ends:
                if self._isa_cbgif_junction():
                    return True
                else:
                    # some of the sites are fixed, but potentially not all
                    # therefor, return a False for the complete interface
                    # THIS CASE SHOULD NOT OCCUR AT ALL ANYMORE !!!!!!
                    return False 
            # check if (1) _optimal_aligned_donor
            # and (2) _optimal_aligned_acceptor
            # and (3) AlignedSites (not Collections)
            # and (4) all nodes/organisms and (5) phase equity
            if self._optimal_aligned_donor and self._optimal_aligned_acceptor:
                if self._optimal_aligned_donor.__class__.__name__ ==\
                'DonorSiteCollectionGraph':
                    return False
                elif self._optimal_aligned_acceptor.__class__.__name__ ==\
                'AcceptorSiteCollectionGraph':
                    return False
                elif self._optimal_aligned_donor.node_count() !=\
                self.organism_set_size():
                    return False
                elif self._optimal_aligned_acceptor.node_count() !=\
                self.organism_set_size():
                    return False
                elif self._optimal_aligned_donor.phase() ==\
                self._optimal_aligned_acceptor.phase() and\
                self._optimal_aligned_donor.phase() in [0,1,2]:
                    return True
                elif self._optimal_aligned_donor.phase() ==\
                self._optimal_aligned_acceptor.phase() and\
                self._optimal_aligned_donor.phase():
                    # check for each organism if the D & A objects have
                    # compatible phases
                    for org in self.organism_set():
                        d = self._optimal_aligned_donor.get_organism_objects(org)[0]
                        a = self._optimal_aligned_acceptor.get_organism_objects(org)[0]
                        if d.phase != a.phase:
                            return False
                    else:
                        # all objects okay!
                        return True
                else:
                    return False
            else:
                return False
        else:
            pass

        # if this point is reached -> return a False
        return False

    # end of function is_compatible


    def _isa_cbgif_junction(self):
        """
        Is this a cbgIFjunction; i.e. 2 CBGs fitting perfectly together
        without an intermediate lsrCBG

        @rtype:  Boolean
        @return: True or False
        """
        if self.donorCBG and self.acceptorCBG and\
        self.donorCBG._short_name == "CBG" and\
        self.acceptorCBG._short_name == "CBG" and\
        self._forced_3p_ends and self._forced_5p_ends and\
        self.donorCBG.get_ordered_nodes() ==\
        self.acceptorCBG.get_ordered_nodes():
            keys = self._forced_3p_ends.keys()
            coords3 = [ (key, self._forced_3p_ends[key].pos) for key in keys ]
            coords5 = [ (key, self._forced_5p_ends[key].pos) for key in keys ]
            if coords3==coords5:
                return True
            else:
                return True # SHOULD BE FALSE MAYBE!?!?! 
        else:
            return False

    # end of function _isa_cbgif_junction


    def is_lacking_organism_node(self,organism=None):
        """
        Does the AlignedDonorSiteGraph or the AlignedAcceptorsSiteGraph
        miss a node (i.e. is a K(s-x) graph)

        @type  organism: string
        @param organism: Organism identifier (or None)

        @attention: when Organism identifier is given, only presence of the
                    respective organism node is checked.

        @rtype:  Boolean
        @return: True or False
        """
        if self.is_donor_lacking_organism_node(organism=organism):
            return True
        elif self.is_acceptor_lacking_organism_node(organism=organism):
            return True
        else:
            return False

    # end of function is_lacking_organism_node


    def is_acceptor_lacking_organism_node(self,organism=None):
        """
        Does the AlignedAcceptorsSiteGraph miss a node (i.e. is a K(s-x) graph)

        @type  organism: string
        @param organism: Organism identifier (or None)

        @attention: when Organism identifier is given, only presence of the
                    respective organism node is checked.

        @rtype:  Boolean
        @return: True or False
        """
        if self._optimal_aligned_acceptor:
            if self._optimal_aligned_acceptor.node_count() !=\
            self.organism_set_size():
                if organism:
                    if organism not in\
                    self._optimal_aligned_acceptor.organism_set():
                        return True
                    else:
                        return False
                else:
                    return True
            else:
                return False
        else:
            # TODO: what do we need/want to return in this case?
            return False

    # end of function is_acceptor_lacking_organism_node 


    def is_donor_lacking_organism_node(self,organism=None):
        """
        Does the AlignedDonorSiteGraph miss a node (i.e. is a K(s-x) graph)

        @type  organism: string
        @param organism: Organism identifier (or None)

        @attention: when Organism identifier is given, only presence of the
                    respective organism node is checked.

        @rtype:  Boolean
        @return: True or False
        """
        if self._optimal_aligned_donor:
            if self._optimal_aligned_donor.node_count() !=\
            self.organism_set_size():
                if organism:
                    if organism not in\
                    self._optimal_aligned_donor.organism_set():
                        return True
                    else:
                        return False
                else:
                    return True
            else:
                return False
        else:
            # TODO: what do we need/want to return in this case?
            return False
        
    # end of function is_donor_lacking_organism_node


    def is_forced_for_organism(self,organism):
        """
        Is the CBGInterface for this organism defined by forced sites?

        @type  organism: string
        @param organism: Organism Identifier

        @rtype:  Boolean 
        @return: True or False
        """
        has3pforced = self._forced_3p_ends.has_key(organism)
        has5pforced = self._forced_5p_ends.has_key(organism)
        if has3pforced and has5pforced:
            return True
        elif not has3pforced and not has5pforced:
            return False
        elif has3pforced and self.acceptorCBG._short_name == "lsrCBG" and\
        organism not in self.acceptorCBG.organism_set():
            return True
        elif has5pforced and self.donorCBG._short_name == "lsrCBG" and\
        organism not in self.donorCBG.organism_set():
            return True
        else:
            # unexpected event; either none or both sides should be defined!
            print "NOT BOTH SIDES DEFINED BY FORCED ENDS", organism
            print self.donorCBG
            print [ (k,str(v)) for k,v in self._forced_3p_ends.iteritems() ]
            print self
            print [ (k,str(v)) for k,v in self._forced_5p_ends.iteritems() ]
            print self.acceptorCBG
            raise "NOT BOTH SIDES DEFINED BY FORCED ENDS: 5p:%s 3p:%s" % (
                    has5pforced, has3pforced )

    # end of function is_forced_for_organism


    def is_meaningfull_intron(self,organism=None):
        """
        Is the assigned intron a biological possible one?
        
        @type  organism: * (string)
        @param organism: Organism Identifier or None
        
        @attention: when Organism identifier is given, the check is performed
                    only for the respective organism
        
        @rtype:  Boolean 
        @return: True or False
        """
        ARE_MEANINGFULL_INTRONS = True
        for org in self.organism_set():
            # only perform check on organism identifier when applied
            if organism and org != organism: continue
            # get assigned donor / acceptor object
            donor    = self.get_assigned_donor_site(org)
            acceptor = self.get_assigned_acceptor_site(org)
            # check if both donor and acceptor sites are available
            if not donor or not acceptor:
                ARE_MEANINGFULL_INTRONS = False
                break
            # get splice sites classes
            dclass = donor.__class__.__name__
            aclass = acceptor.__class__.__name__
            if dclass == 'ProjectedSpliceDonor' and\
            aclass == 'ProjectedSpliceAcceptor':
                # This should be fine, howeever in exceptional cases the
                # projected sites could be non-identical positioned.
                # This fixed in _fix_gap_between_projected_splice_sites()
                # lateron. But, a large gap in an earlier stage could hint on
                # a bogus projection and thus bogus intron(s).
                if acceptor.pos - donor.pos >= MIN_INTRON_NT_LENGTH:
                    ARE_MEANINGFULL_INTRONS = False
                    break
            elif dclass == 'SpliceDonor' and\
            aclass == 'SpliceAcceptor':
                # No phase check, pssm check, entropy check or what soever.
                # That is all taken care for in is_optimal(), is_compatible(),
                # is_optimal_donor() and is_optimal_acceptor() functions.
                # Here, we just check the intron.
                # TODO: check on branch point, tcode score, coding propensity??
                if acceptor.pos - donor.pos < MIN_INTRON_NT_LENGTH:
                    ARE_MEANINGFULL_INTRONS = False
                    break
            else:
                # all other combinations -> no meaningfull intron
                # TODO: combinations of a ProjectedSpliceSite and a Splicesit
                # TODO: can occur in cases of an inframe intron, where the
                # TODO: projected site falls exactly on the actual site ;-((
                # TODO: Here, False is returned -> which is okay!
                # TODO: But, these cases must be recognized in the functions:
                # TODO: self.find_conserved_splice_sites()
                # TODO: self.optimizetinyexoninterface()
                ARE_MEANINGFULL_INTRONS = False
                break

        # return function's outcome
        return ARE_MEANINGFULL_INTRONS
    
    # end of function is_meaningfull_intron


    def is_optimal(self):
        """
        Is this CodingBlockInterface the most optimal possible one?

        @rtype:  Boolean 
        @return: True or False
        """
        if self.is_optimal_donor():
            if self.is_optimal_acceptor():
                if self.is_meaningfull_intron():
                    return True
                else:
                    return False
            else:
                return False
        else:
            return False

    # end of function is_optimal


    def is_acceptable(self,organism=None):
        """
        Is this CodingBlockGraphInterface acceptable (but not per se optimal)
        as a CodingBlockGraphInterface?
        
        @type  organism: * (string)
        @param organism: Organism Identifier or None
        
        @attention: when Organism identifier is given, the check is performed
                    only for the respective organism
        
        @rtype:  Boolean
        @return: True or False
        """
        if self.is_optimal():
            return True
        elif organism and self._forced_3p_ends.has_key(organism) and\
        self._forced_5p_ends.has_key(organism):
            # a forced Non-intron that is not recognized as a lsrCBG interface.
            # But, forced ends overrule site alignment!
            return True
        elif not self.is_meaningfull_intron():
            return False
        elif not self._splicedonorgraph:
            return False
        elif not self._spliceacceptorgraph:
            return False
        elif self.organism_set_size() ==\
        self._splicedonorgraph.organism_set_size() and\
        self.organism_set_size() ==\
        self._spliceacceptorgraph.organism_set_size():
            # for each organism there is at least a single site present in both
            # donor and acceptor graph at the interface
            if self.is_compatible() and not self.is_donor_improveable() and\
            not self.is_acceptor_improveable():
                # not optimizable any further
                return True
            elif self.is_highest_pssm_phase(self.phase()) and\
            self.is_highest_pssm_phase(self.phase(),organism=organism):
                # very unlikely that it is optimizable
                return True
            elif self.is_compatible() and\
            (self.is_donor_improveable() or self.is_acceptor_improveable()):
                # TODO: can we improve this site indeed?? At least, it is a
                # severe warning that for the given organism there is doubt
                # in splice site choice
                if not organism:
                    return False
                else:
                    if self.is_donor_improveable(organism=organism) or\
                    self.is_acceptor_improveable(organism=organism):
                        return None
                    else:
                        return True
            else:
                return False
        else:
            return False 

    # end of function is_acceptable 


    def optimalitycheck(self):
        """
        Return list that represents the optimality of this interface

        @attention: used by checking if cbgIF.optimalitycheck().count(True) >= 2

        @rtype:  list 
        @return: [ Boolean, Boolean, Boolean ]
        """
        return [ self.is_compatible(), self.is_optimal_donor(),
                 self.is_optimal_acceptor() ]

    # end of function optimalitycheck


    def is_optimal_donor(self,organism=None,verbose=False):
        """
        Are the individual DonorSites in the optimal AlignedDonorSiteGraph the
        optimal donor sites in the elegiable range of donor sites?
    
        @type  organism: * (string)
        @param organism: Organism Identifier or None

        @type  verbose: Boolean
        @param verbose: in case True, print logging/debugging messages to STDOUT

        @attention: lsrCBGs in the CodingBlockInterface will return True
        @attention: incomplete/empty CodingBlockInterfaces will return False
        @attention: when Organism identifier is given, the check is
                    performed only for the respective organism
    
        @rtype:  Boolean      
        @return: True or False
        """
        # set verbose parameter to True in case cbgIF.verbose is True
        if self.verbose == True: verbose = True

        if self.donorCBG._short_name == "CBG" and\
        self.acceptorCBG._short_name == "lsrCBG":
            return True
        elif self.donorCBG._short_name == "lsrCBG" and\
        self.acceptorCBG._short_name == "CBG":
            return True
        elif self._isa_cbgif_junction():
            return True
        elif not self.is_compatible() and not organism:
            return False
        elif not self._optimal_aligned_donor or\
        self._optimal_aligned_donor.__class__.__name__ ==\
        'DonorSiteCollectionGraph':
            return False
        elif self.donorCBG._short_name == "CBG" and\
        self.acceptorCBG._short_name == "CBG":
            # because is_compatible(), no phase check needed anymore
            # just check if it is the highest scoring available pssm site
            # see self._optimal_donor_analyses() for description of the
            # structure of the dict 'data'
            data = self._optimal_donor_analyses(use_entropy=self._USE_ENTROPY)

            if organism and data.has_key(organism):
                return data[organism][1][0] == data[organism][1][1]
            elif organism and organism in self.organism_set():
                # no site available in the Collection at all for this organism!?
                return False
            elif organism:
                # organism identifier provided but not present as a dict key
                # in self._optimal_donor_analyses() -> Error
                raise OrganismNotPresentInGraph
            else:
                check = [ data[org][1][0] == data[org][1][1] for org in data.keys() ]
                if False in check: return False
                else:              return True
        else:
            pass

        # if this point is reached -> return an (unexpected) False
        return False

    # end of function is_optimal_donor


    def is_highest_pssm_phase(self,phase,organism=None):
        """
        Is this phase the highest scoring summed pssm (for this organism)

        @type  phase: 0,1,2
        @param phase: intron/cbgIF phase 
 
        @type  organism: * (string)
        @param organism: Organism Identifier or None

        @attention: lsrCBGs in the CodingBlockInterface will return False
        @attention: incomplete/empty CodingBlockInterfaces will return False
        @attention: when Organism identifier is given, the highest pssm phase
                    is checked only for the respective organism

        @rtype:  Boolean
        @return: True or False
        """

        if phase == self.highest_pssm_phase(organism=organism):
            return True
        else:
            return False

    # end of function is_highest_pssm_phase


    def highest_pssm_phase(self,organism=None):
        """
        Get the highest scoring summed pssm phase (for this organism)

        @type  organism: * (string)
        @param organism: Organism Identifier or None

        @attention: lsrCBGs in the CodingBlockInterface will return None
        @attention: incomplete/empty CodingBlockInterfaces will return None
        @attention: when Organism identifier is given, the highest pssm phase
                    is checked only for the respective organism

        @rtype  phase: 0,1,2 or None
        @return phase: intron/cbgIF phase
        """
        # if a splitted CBFinterface, not meaningfull -> return None
        if self._forced_3p_ends and self._forced_5p_ends:
            return None
        # obtain highest pssm score
        data = {}
        for org in self.organism_set():
            if organism and org != organism: continue
            data[org] = {}
            for phase in [0,1,2]:
                try:
                    donor = self._splicedonorgraph.get_optimal_single_site(
                        org,phase=phase)
                    accep = self._spliceacceptorgraph.get_optimal_single_site(
                        org,phase=phase)
                except AttributeError:
                    # when no _splicedonorgraph or _spliceacceptorgraph -> None
                    donor,accep = None,None
                if donor and accep:
                    data[org][phase] = donor.pssm_score + accep.pssm_score
                else:
                    # do not create a dict elem for this phase
                    pass
                
        # now get cumulative scores for each phase
        summation = [ [0.0, 0], [0.0, 1], [0.0, 2] ]
        for phase in [0,1,2]:
            for org,phasedict in data.iteritems():
                if phasedict.has_key(phase):
                    summation[phase][0] += phasedict[phase]

        # remove phases that have not been increased; that means there
        # are no splice sites of that phase available
        for phase in [2,1,0]:
            if summation[phase][0] == 0.0:
                summation.pop(phase)

        if summation:
            # order the summation and return the
            # highest phase (== the final list elem)
            summation.sort()
            return summation[-1][1]
        else:
            # no phase with splice site donors & acceptors
            return None

    # end of function highest_pssm_phase


    def is_highest_pssm_donor(self,organism=None):
        """
        Is the optimal DonorSite (of this organism) the highest PSSM scoring
        site in the elegiable range of DonorSites?

        @type  organism: * (string)
        @param organism: Organism Identifier or None

        @attention: lsrCBGs in the CodingBlockInterface will return False
        @attention: incomplete/empty CodingBlockInterfaces will return False
        @attention: when Organism identifier is given, the check is
                    performed only for the respective organism

        @rtype:  Boolean      
        @return: True or False
        """
        if not self._optimal_aligned_donor:
            return False
        if self._optimal_aligned_donor.__class__.__name__ ==\
        'DonorSiteCollectionGraph':
            return False

        # get data for each of the organisms/genes; use_entropy == False
        # because here the interest is only actual PSSM score, not position
        # see self._optimal_donor_analyses() for description of the
        # structure of the dict 'data'
        data = self._optimal_donor_analyses(use_entropy=False) 
        
        if organism and data.has_key(organism):
            return data[organism][1][0] == data[organism][1][1] 
        elif not data:
            return False
        elif organism:
            raise OrganismNotPresentInGraph 
        else:
            check = [ data[org][1][0]==data[org][1][1] for org in data.keys() ]
            if False in check: return False
            else:              return True

    # end of function is_highest_pssm_donor


    def is_highest_pssm_acceptor(self,organism=None):
        """
        Is the optimal AcceptorSite (of this organism) the highest PSSM scoring
        site in the elegiable range of AcceptorSites?

        @type  organism: * (string)
        @param organism: Organism Identifier or None

        @attention: lsrCBGs in the CodingBlockInterface will return False
        @attention: incomplete/empty CodingBlockInterfaces will return False
        @attention: when Organism identifier is given, the check is
                    performed only for the respective organism

        @rtype:  Boolean      
        @return: True or False
        """
        if not self._optimal_aligned_acceptor: 
            return False
        if self._optimal_aligned_acceptor.__class__.__name__ ==\
        'AcceptorSiteCollectionGraph':
            return False

        # get data for each of the organisms/genes; use_entropy == False
        # because here the interest is only actual PSSM score, not position
        data = self._optimal_acceptor_analyses(use_entropy=False)

        if organism and data.has_key(organism):
            return data[organism][1][0] == data[organism][1][1]
        elif not data:
            return False
        elif organism:
            raise OrganismNotPresentInGraph
        else:
            check = [ data[org][1][0]==data[org][1][1] for org in data.keys() ]
            if False in check: return False
            else:              return True

    # end of function is_highest_pssm_acceptor


    def _optimal_donor_analyses(self,use_entropy=True):
        """
        Get for each organism the *best* site in the DonorSiteCollectionGraph

        @type  use_entropy: Boolean 
        @param use_entropy: Use entropy/position of DonorSite to get the
                            *best* site (True), or use only pssm score (False) 

        @rtype:  dict 
        @return: dictionary with Organism identifiers as keys,
                 tuple structure as values: ( (a,b), (c,d), )
                    a:  PSSM score of the selected DonorSite
                    b:  PSSM score of the optimal DonorSite
                    c:  position (in nt) of the selected DonorSite
                    d:  position (in nt) of the optimal DonorSite
        """
        data = {} 
        if use_entropy:
            omsr = self.donorCBG.overall_minimal_spanning_range()
        for node in self._optimal_aligned_donor.get_nodes():
            org    = self.organism_by_node(node)
            cbgnode= tuple(node[0:2])
            if use_entropy: coord = max(omsr[cbgnode])
            else:           coord = None
            ds_sel = self._optimal_aligned_donor.get_node_object(node)
            ds_opt = self._splicedonorgraph.get_optimal_single_site(
                org,use_entropy_around_omsr=coord)
            if ds_opt:
                data[org] = (
                    ( ds_sel.pssm_score, ds_opt.pssm_score ),
                    ( ds_sel.pos, ds_opt.pos )
                    )
            else:
                # zero DonorSites of this organism/gene in the DSCG
                pass
            
        # return data dict
        return data

    # end of function _optimal_donor_analyses 


    def is_optimal_acceptor(self,organism=None,verbose=False):
        """
        Are the individual AcceptorSites in the optimal AlignedAcceptorSiteGraph
        the optimal acceptor sites in the elegiable range of acceptor sites?
    
        @type  organism: * (string)
        @param organism: Organism Identifier or None

        @type  verbose: Boolean
        @param verbose: in case True, print logging/debugging messages to STDOUT

        @attention: lsrCBGs in the CodingBlockInterface will return True
        @attention: incomplete/empty CodingBlockInterfaces will return False
        @attention: when Organism identifier is given, the check is
                    performed only for the respective organism
    
        @rtype:  Boolean      
        @return: True or False
        """
        # set verbose parameter to True in case cbgIF.verbose is True
        if self.verbose == True: verbose = True

        if self.donorCBG._short_name == "CBG" and\
        self.acceptorCBG._short_name == "lsrCBG":
            return True
        elif self.donorCBG._short_name == "lsrCBG" and\
        self.acceptorCBG._short_name == "CBG":
            return True
        elif self._isa_cbgif_junction():
            return True
        elif not self.is_compatible() and not organism:
            return False
        elif not self._optimal_aligned_acceptor or\
        self._optimal_aligned_acceptor.__class__.__name__ ==\
        'AcceptorSiteCollectionGraph':
            return False
        elif self.donorCBG._short_name == "CBG" and\
        self.acceptorCBG._short_name == "CBG":
            # because is_compatible(), no phase check needed anymore
            # just check if it is the highest scoring available pssm site
            # see self._optimal_donor_analyses() for description of the
            # structure of the dict 'data'
            data = self._optimal_acceptor_analyses(use_entropy=self._USE_ENTROPY)

            if organism and data.has_key(organism):
                return data[organism][1][0] == data[organism][1][1]
            elif organism and organism in self.organism_set():
                # no site available in the Collection at all for this organism!?
                return False
            elif organism:
                raise OrganismNotPresentInGraph
            else:
                check = [ data[org][1][0] == data[org][1][1] for org in data.keys() ]
                if False in check: return False
                else:              return True
        else:
            pass

        # if this point is reached -> return an (unexpected) False
        return False

    # end of function is_optimal_acceptor


    def _optimal_acceptor_analyses(self,use_entropy=True):
        """
        Get for each organism the *best* site in the AcceptorSiteCollectionGraph

        @type  use_entropy: Boolean 
        @param use_entropy: Use entropy/position of AcceptorSite to get the
                            *best* site (True), or use only pssm score (False) 

        @rtype:  dict 
        @return: dictionary with Organism identifiers as keys,
                 tuple structure as values: ( (a,b), (c,d), )
                    a:  PSSM score of the selected AcceptorSite
                    b:  PSSM score of the optimal AcceptorSite
                    c:  position (in nt) of the selected AcceptorSite
                    d:  position (in nt) of the optimal AcceptorSite
        """
        data = {} 
        if use_entropy:
            omsr = self.acceptorCBG.overall_minimal_spanning_range()
        for node in self._optimal_aligned_acceptor.get_nodes():
            org    = self.organism_by_node(node)
            cbgnode= tuple(node[0:2])
            if use_entropy: coord = min(omsr[cbgnode])
            else:           coord = None
            as_sel = self._optimal_aligned_acceptor.get_node_object(node)
            as_opt = self._spliceacceptorgraph.get_optimal_single_site(
                org,use_entropy_around_omsr=coord)
            if as_opt:
                data[org] = (
                    ( as_sel.pssm_score, as_opt.pssm_score),
                    ( as_sel.pos, as_opt.pos )
                    )
            else:
                # zero AcceptorSites of this organism/gene in the ASCG
                pass

        # return data dict
        return data

    # end of function _optimal_acceptor_analyses


    def is_majority_of_donors_optimal(self):
        """
        Is the majority of DonorSites in the optimal ADSG, optimal donors?

        @attention: incomplete/empty CodingBlockInterfaces will return False
        @attention: incompatible CodingBlockInterfaces will return False
        @attention: lsrCBGs in the CodingBlockInterface will return True

        @rtype:  Boolean
        @return: True or False
        """
        # only valid for compatible sites!
        if not self.is_compatible(): return False
        improvability = []
        for org in self.organism_set():
            if self.is_optimal_donor(organism=org):
                improvability.append( False )
            else:
                improvability.append( self.is_donor_improveable(organism=org) )
        # check improvability; if False >= True, majority is optimal
        if improvability.count(False) >= improvability.count(True):
            return True
        else:
            return False

    # end of function is_majority_of_donors_optimal


    def is_majority_of_acceptors_optimal(self):
        """
        Is the majority of AcceptorSites in the optimal AASG, optimal acceptors?

        @attention: incomplete/empty CodingBlockInterfaces will return False
        @attention: incompatible CodingBlockInterfaces will return False
        @attention: lsrCBGs in the CodingBlockInterface will return True

        @rtype:  Boolean
        @return: True or False
        """
        # only valid for compatible sites!
        if not self.is_compatible(): return False
        improvability = []
        for org in self.organism_set():
            if self.is_optimal_acceptor(organism=org):
                improvability.append( False )
            else:
                improvability.append(self.is_acceptor_improveable(organism=org))
        # check improvability; if False >= True, majority is optimal
        if improvability.count(False) >= improvability.count(True):
            return True
        else:
            return False

    # end of function is_majority_of_acceptors_optimal

    #############################################################
    #### Obtain consideredsplicesiterange()                  ####
    #############################################################

    def consideredsplicesiterange(self):
        """
        Get the considered splice site range of this cbgIF

        @rtype:  dict 
        @return: dictionary with Organism identifiers as keys,
                 tuple structure as values: ( (a,b), (c,d), )
                    a:  start nt coordinate of the elegiable donor range
                    b:    end nt coordinate of the elegiable donor range
                    a:  start nt coordinate of the elegiable acceptor range
                    b:    end nt coordinate of the elegiable acceptor range

        @attention: lsrCBGs result in None values in stead of integers (nt)
        """
        data = {}
        for org in self.organism_set():
            if self._splicedonorgraph:
                donorrange = self._splicedonorgraph.get_consideredsplicesiterange(org)
            else:
                donorrange = ( None, None )
            if self._spliceacceptorgraph:
                acceptorrange = self._spliceacceptorgraph.get_consideredsplicesiterange(org)
            else:
                acceptorrange = ( None, None )
            data[org] = ( donorrange, acceptorrange )
        # return data dictionary
        return data

    # end of function consideredsplicesiterange

    ########################################################################
    #### Functions for checking phase (equity) of the splice sites      ####
    ########################################################################

    def phase(self):
        """
        Get the phase of this CodingBlockInterface

        @attention: incompatible donor & acceptor phases will return None
        @attention: incomplete/empty CodingBlockInterfaces will return None

        @rtype  phase: 0,1,2 or None
        @return phase: 0,1,2 or None
        """
        if self.is_compatible():
            return self._optimal_aligned_donor.phase()
        else:
            return None

    # end of function phase


    def donor_phase(self):
        """
        Get the phase of the optimal AlignedDonorSiteGraph

        @attention: incomplete/empty CodingBlockInterfaces will return None

        @rtype  phase: 0,1,2 or None
        @return phase: 0,1,2 or None
        """
        if self._optimal_aligned_donor and\
        self._optimal_aligned_donor.__class__.__name__ !=\
        'DonorSiteCollectionGraph':
            return self._optimal_aligned_donor.phase()
        else:
            return None # no proper phase assigned yet!

    # end of function donor_phase  


    def acceptor_phase(self):
        """
        Get the phase of the optimal AlignedAcceptorSiteGraph

        @attention: incomplete/empty CodingBlockInterfaces will return None

        @rtype  phase: 0,1,2 or None
        @return phase: 0,1,2 or None
        """
        if self._optimal_aligned_acceptor and\
        self._optimal_aligned_acceptor.__class__.__name__ !=\
        'AcceptorSiteCollectionGraph': 
            return self._optimal_aligned_acceptor.phase()
        else:
            return None # no proper phase assigned yet! 

    # end of function acceptor_phase

    ############################################################################
    #### Functions for obtaining the DSCG and ASCG and related functions    ####
    ############################################################################

    def project_missing_introns(self):
        """
        Project missing introns in between donor and acceptor CBG

        @attention: USE WITH CARE! Only called from the harvest_donor_sites()
                    and harvest_acceptor_sites() function.
        """
        # reset projected splice site lists to empty
        self._projected_acceptor_sites  = {}
        self._projected_donor_sites     = {}

        if self.donorCBG.IS_3P_SPLITTED == True and\
        self.acceptorCBG.IS_5P_SPLITTED == True:
            # no intron projection on a splitted interface
            pass
        elif not graphPlus.comparison.mutual_nodes(self.donorCBG,
        self.acceptorCBG):
            # no common nodes shared by both graphs => no projection!
            pass
        else:
            # make intron projections!
            projected_introns = project_missing_introns(
                    ( self.donorCBG, self.acceptorCBG ),
                    interface_intron_status = self._interface_is_intron,
                    inframe_intron_min_aa_length=INFRAME_INTRON_MIN_AA_LENGTH,
                    )

            # remove dict key when organism is not in one of both CBGs
            # this happens in K(s-x) graphs
            for org in self.organism_set():
                if org not in self.donorCBG.organism_set() and\
                projected_introns.has_key(org):
                    del( projected_introns[org] )
                if org not in self.acceptorCBG.organism_set() and\
                projected_introns.has_key(org): 
                    del( projected_introns[org] )

            # only allow projected introns in the proper designated area
            for org,intronlist in projected_introns.iteritems():
                mindonor =\
                    self.donorCBG.minimal_eligable_donor_site_position(org)
                maxaccep =\
                    self.acceptorCBG.maximal_eligable_acceptor_site_position(org)
                prevlen = len(intronlist)
                for pos in range(len(intronlist)-1,-1,-1):
                    intron = intronlist[pos]
                    if intron.projected_donor_pos in\
                    range(mindonor[1],maxaccep[1]+1) and\
                    intron.projected_acceptor_pos in\
                    range(mindonor[1],maxaccep[1]+1):
                        pass
                    else:
                        # remove this projectedintron -> not in proper range!
                        intronlist.pop(pos)

            # get the (unique) projected donor/acceptor sites
            # from the projected introns
            proj_acceptors, proj_donors     =\
                    projectedintrons2projectedsites(projected_introns)
            self._projected_donor_sites     = proj_donors
            self._projected_acceptor_sites  = proj_acceptors

        # check checker for projection to True
        self._HAS_INTRONS_PROJECTED = True

        ########################################################################
        if self.verbose:
            print "IntronProjection: %s donors, %s acceptors" % (
                sum([ len(v) for v in self._projected_donor_sites.values() ]),
                sum([ len(v) for v in self._projected_acceptor_sites.values() ])
                )
        ########################################################################

    # end of function project_missing_introns


    def harvest_donor_sites(self,
        allow_phase_shift=False,
        allow_non_canonical=False,
        enlarge_5p_boundary_by=None,
        enlarge_3p_boundary_by=None,
        aligned_donor_max_triplet_distance=ALIGNED_DONOR_MAX_TRIPLET_DISTANCE,
        store_all_projected_sites=False):
        """
        (Re)create the DonorSiteCollectionGraph (DSCG)

        @type  allow_phase_shift: Boolean
        @param allow_phase_shift: create additional edges in the DSCG for
                    DonorSites of distinct phase (when less than
                    MAX_SPLICE_SITE_PHASE_SHIFT_NT_DISTANCE apart) when True,
                    default False

        @type  allow_non_canonical: Boolean
        @param allow_non_canonical: create nodes for non-canonical DonorSites
                    when set on True, default False

        @type  enlarge_5p_boundary_by: None or integer
        @param enlarge_5p_boundary_by: in- or decrease the 5' side defined by
                    self.consideredsplicesiterange()

        @type  enlarge_3p_boundary_by: None or integer
        @param enlarge_3p_boundary_by: in- or decrease the 3' side defined by
                    self.consideredsplicesiterange()

        @type  aligned_donor_max_triplet_distance: integer
        @param aligned_donor_max_triplet_distance: maximum AA distance between
                    DonorSites (of identical phase) in a PacbPORF for edge
                    creation in the DSCG

        @type  store_all_projected_sites: Boolean
        @param store_all_projected_sites: add ProjectedDonorSites to the DSCG
                    even when located outside self.consideredsplicesiterange()
                    when True, default False
        """
        ########################################################################
        if self.verbose:
            # StopWatch time logger in verbose mode
            stw = StopWatch(name="harvestdonor (aps:%s anc:%s)" % (
                allow_phase_shift,allow_non_canonical) )
            stw.start()
        ########################################################################

        # project introns if not done yet!
        if not self._HAS_INTRONS_PROJECTED:
            self.project_missing_introns()
            ####################################################################
            if self.verbose:
                print "projected sites DONE!"
                print self._projected_donor_sites
                print self._projected_acceptor_sites
            ####################################################################

        if not self.donorCBG.IS_3P_SPLITTED:
            # scan donorCBG for donorsites
            self.donorCBG.harvest_elegiable_donor_sites(
                    projected_donors          = self._projected_donor_sites,
                    forced_codingblock_ends   = self._forced_3p_ends,
                    next                      = self.acceptorCBG,
                    allow_phase_shift         = allow_phase_shift,
                    store_all_projected_sites = store_all_projected_sites,
                    allow_non_canonical       = allow_non_canonical,
                    enlarge_5p_boundary_by    = enlarge_5p_boundary_by,
                    enlarge_3p_boundary_by    = enlarge_3p_boundary_by,
                    aligned_donor_max_triplet_distance =\
                            aligned_donor_max_triplet_distance 
                    )
        else:
            # donorCBG IS_3P_SPLITTED
            pass

        # set attributes for donor site settings
        self._donor_allow_phase_shift = allow_phase_shift 
        self._donor_allow_non_canonical = allow_non_canonical 

        # update splicesitegraph attributes
        self._splicedonorgraph  = self.donorCBG._splicedonorgraph

        # remove donor nodes for organisms/genes that have
        # definately NO intron interface
        for organism,intronstatus in self._interface_is_intron.iteritems():
            if intronstatus == False:
                # remove all the SpliceSites, leave Projections intact! 
                status = self._remove_object_type_of_organism_from_graph(
                        organism,'SpliceDonor','_splicedonorgraph')

        ########################################################################
        if self.verbose:
            print stw.lap(), 
            if self._splicedonorgraph:
                print self._splicedonorgraph.node_count(), "nodes,",
                print self._splicedonorgraph.edge_count(), "edges"
            else:
                print ""
        ########################################################################


    # end of function harvest_donor_sites


    def harvest_acceptor_sites(self,
        allow_phase_shift=False,
        allow_non_canonical=False,
        enlarge_5p_boundary_by=None,
        enlarge_3p_boundary_by=None,
        aligned_acceptor_max_triplet_distance=\
            ALIGNED_ACCEPTOR_MAX_TRIPLET_DISTANCE,
        store_all_projected_sites=False):
        """
        (Re)create the AcceptorSiteCollectionGraph (ASCG)

        @type  allow_phase_shift: Boolean
        @param allow_phase_shift: create additional edges in the ASCG for
                    AcceptorSites of distinct phase (when less than
                    MAX_SPLICE_SITE_PHASE_SHIFT_NT_DISTANCE apart) when True,
                    default False

        @type  allow_non_canonical: Boolean
        @param allow_non_canonical: create nodes for non-canonical AcceptorSites
                    when set on True, default False

        @attention: non-canonical AcceptorSites are not supported (yet)

        @type  enlarge_5p_boundary_by: None or integer
        @param enlarge_5p_boundary_by: in- or decrease the 5' side defined by
                    self.consideredsplicesiterange()

        @type  enlarge_3p_boundary_by: None or integer
        @param enlarge_3p_boundary_by: in- or decrease the 3' side defined by
                    self.consideredsplicesiterange()

        @type  aligned_acceptor_max_triplet_distance: integer
        @param aligned_acceptor_max_triplet_distance: maximum AA distance
                    between AcceptorSites (of identical phase) in a PacbPORF for
                    edge creation in the ASCG

        @type  store_all_projected_sites: Boolean
        @param store_all_projected_sites: add ProjectedAcceptorSites to the ASCG
                    even when located outside self.consideredsplicesiterange()
                    when True, default False
        """
        ########################################################################
        if self.verbose:
            # StopWatch time logger in verbose mode
            stw = StopWatch(name="harvestaccep (aps:%s anc:%s)" % (
                allow_phase_shift,allow_non_canonical) )
            stw.start()
        ########################################################################

        # project introns if not done yet!
        if not self._HAS_INTRONS_PROJECTED: self.project_missing_introns()

        if not self.acceptorCBG.IS_5P_SPLITTED:
            # scan acceptorCBG for acceptorsites
            self.acceptorCBG.harvest_elegiable_acceptor_sites(
                    forced_codingblock_ends      = self._forced_5p_ends,
                    prev                         = self.donorCBG,
                    allow_phase_shift            = allow_phase_shift,
                    store_all_projected_sites    = store_all_projected_sites,
                    enlarge_5p_boundary_by       = enlarge_5p_boundary_by,
                    enlarge_3p_boundary_by       = enlarge_3p_boundary_by,
                    # TODO -> implement this variable in function!
                    # ALLOW_NON_CANONICAL_ACCEPTOR = allow_non_canonical,
                    projected_acceptors          =\
                        self._projected_acceptor_sites,
                    aligned_acceptor_max_triplet_distance =\
                        aligned_acceptor_max_triplet_distance 
                    )
        else:
            # acceptorCBG IS_5P_SPLITTED
            pass

        # set attributes for acceptor site settings
        self._acceptor_allow_phase_shift   = allow_phase_shift
        self._acceptor_allow_non_canonical = allow_non_canonical

        # update splicesitegraph attributes
        self._spliceacceptorgraph = self.acceptorCBG._spliceacceptorgraph

        # remove acceptor nodes for organisms/genes that have
        # definately NO intron interface
        for organism,intronstatus in self._interface_is_intron.iteritems():
            if intronstatus == False:
                # remove all the SpliceSites, leave Projections intact!
                status = self._remove_object_type_of_organism_from_graph(
                        organism,'SpliceAcceptor','_spliceacceptorgraph')

        ########################################################################
        if self.verbose:
            print stw.lap(),
            if self._spliceacceptorgraph: 
                print self._spliceacceptorgraph.node_count(), "nodes,",
                print self._spliceacceptorgraph.edge_count(), "edges"
            else:
                print ""
        ########################################################################

    # end of function harvest_acceptor_sites


    def harvest_splice_sites(self,
        allow_phase_shift=False,
        allow_non_canonical=False,
        store_all_projected_sites=False):
        """
        (Re)create DonorSiteCollectionGraph and AcceptorSiteCollectionGraph

        @type  allow_phase_shift: Boolean
        @param allow_phase_shift: create additional edges in the SSCG for
                    SpliceSites of distinct phase (when less than
                    MAX_SPLICE_SITE_PHASE_SHIFT_NT_DISTANCE apart) when True,
                    default False

        @type  allow_non_canonical: Boolean
        @param allow_non_canonical: create nodes for non-canonical SpliceSites
                    when set on True, default False

        @type  store_all_projected_sites: Boolean
        @param store_all_projected_sites: add ProjectedSpliceSites to the SSCG
                    even when located outside self.consideredsplicesiterange()
                    when True, default False

        @attention: shortcut function for self.harvest_donor_sites() and
                    self.harvest_acceptor_sites().
        """
        self.harvest_donor_sites(allow_phase_shift=allow_phase_shift,
                store_all_projected_sites=store_all_projected_sites,
                allow_non_canonical=allow_non_canonical)
        self.harvest_acceptor_sites(allow_phase_shift=allow_phase_shift,
                store_all_projected_sites=store_all_projected_sites,
                allow_non_canonical=allow_non_canonical)

    # end of function harvest_splice_sites


    def interface_is_intron(self,organism):
        """
        Is an intron required in this CodingBlockInterface for this organism? 

        @type  organism: * (string)
        @param organism: Organism Identifier

        @rtype:  NoneBoolean
        @return: True   when an intron is required
                 None   when an intron is doubtfull
                 False  when an intron is forbidden (e.g. lsrCBG)
        """
        if organism not in self.organism_set():
            raise OrganismNotPresentInGraph

        if self._forced_3p_ends and self._forced_5p_ends:
            # splitted Interfaces do NOT represent introns!
            return False 

        if organism not in self.donorCBG.organism_set():
            # interface that lacks an organism in the donor
            return None

        if organism not in self.acceptorCBG.organism_set():
            # interface that lacks an organism in the acceptor 
            return None

        # check for mutual nodes and convert to organisms
        mutual_nodes = graphPlus.comparison.mutual_nodes(
            self.donorCBG,self.acceptorCBG)
        mutual_orgs  = [ self.donorCBG.organism_by_node(node) for\
            node in mutual_nodes ]

        # if organism not in mutual_orgs, 2 distinct orfs -> intron required
        if organism not in mutual_orgs: return True

        # do omsr/maxsr analyses
        aadistmaxsr = self.donorCBG.maxsr_distance_between_codingblocks(
            self.acceptorCBG)

        # if aadistmaxsr > 0, a (very) doubtfull case. Possible inframe intron. 
        if aadistmaxsr[organism] > 0: return None

        # do a check on the omsr/maxsr region
        omsrD   = self.donorCBG.overall_minimal_spanning_range(
                    organism=organism)
        omsrA   = self.acceptorCBG.overall_minimal_spanning_range(
                    organism=organism)
        orgnode = self.donorCBG.node_by_organism(organism)

        interfacerange = Set(range(max(omsrD)+1,min(omsrA)))
        rangecount = {}
        for pos in interfacerange: rangecount[pos] = 0

        # TODO: we can 'measure' the alignment quality here too!
        # TODO: a first effort was made earlier here.
        # TODO: old code can be accessed in earlier version for future update.

        for node in mutual_nodes:
            if node == orgnode: continue
            pacbporfA = self.acceptorCBG.get_pacbps_by_nodes(
                node1=orgnode,node2=node)[0]
            interAlen = len(interfacerange.intersection(
                pacbporfA.alignment_protein_range_query()))
            pacbporfD = self.donorCBG.get_pacbps_by_nodes(
                node1=orgnode,node2=node)[0]
            interDlen = len(interfacerange.intersection(
                pacbporfD.alignment_protein_range_query()))
            if interDlen >= interAlen:
                for aapos in interfacerange.intersection(
                pacbporfD.alignment_protein_range_query()):
                    rangecount[aapos] += 1
            else:
                for aapos in interfacerange.intersection(
                pacbporfA.alignment_protein_range_query()):
                    rangecount[aapos] += 1

        # obtain scores for the interface
        maxscore     = (len(mutual_nodes)-1)*len(rangecount)
        totalscore   = sum(rangecount.values())
        nonzeroscore = len(rangecount) - rangecount.values().count(0)

        # return the outcome: is this an intron (True, None or False)
        if len(rangecount) == 0 and len(mutual_nodes) >= 2:
            # hmmm... the CBGs glue perfectly together or even overlap
            # that might be a (clear) sign of a sequence error in one
            # of the other sequences! The continuity in the CBGs is
            # supported by at least 2 organisms/genes, so safe to 
            # enforce a non-intron here -> return a False
            return False
        elif maxscore == 0 and len(rangecount) > 0:
            return None  # just a single mutual_node -> no mappable evidence!
        elif totalscore == maxscore:
            return False # no intron here!
        elif nonzeroscore == len(rangecount):
            # check the OMSR coords of the other CBGs
            distances = self.donorCBG.distance_between_codingblocks(
                self.acceptorCBG)
            thisdist  = distances[organism] - min(distances.values())
            if thisdist >= 11 and min(distances.values()) <= 4:
                return True  # Force an INFRAME intron here!
            else:
                return False # no intron here!
        else:
            return None  # doubtfull case....

    # end of function interface_is_intron

    ############################################################################
    #### Function for removing object(type)s from the DSCG and ASCG         ####
    #### THis function is written generic for a graph with node_objects     ####
    ############################################################################

    def _remove_object_type_of_organism_from_graph(self,organism,
        objclassname,graphattributename,verbose=False):
        """
        @type  organism: * (string)
        @param organism: Organism Identifier
        
        @type  objclassname: string
        @param objclassname: literal class name of the to be removed objects
        
        @type  graphattributename: string
        @param graphattributename: literal attribute name where the graph object
                                   from which the objects must be removed
        
        @attention: only applicable on graphs with gra._node_object attribute!
        
        @type  verbose: Boolean
        @param verbose: in case True, print logging/debugging messages to STDOUT
        
        @rtype:  True
        @return: (non-meaningfull) True
        """
        gra = getattr(self,graphattributename)
        if gra and organism in gra.organism_set():
            splicesitenodes = gra.get_organism_nodes(organism)
            ####################################################################
            if verbose:
                print "remove '%s' of '%s' from '%s' (current %s)" % (
                    objclassname, organism, graphattributename,
                    len(splicesitenodes)
                    )
            ####################################################################

            for node in splicesitenodes:
                nodeobj = gra.get_node_object(node)
                if nodeobj.__class__.__name__ == objclassname:
                    gra.del_node(node)
            # run _update_after_changes() function
            gra._update_after_changes()

            ####################################################################
            if verbose:
                cnt = 0
                if organism in gra.organism_set():
                    cnt = len(gra.get_organism_nodes(organism))
                print "removal done!", organism, cnt, "sites",
                print len(gra._node_object)
            ####################################################################

        # return a non-meaningfull True
        return True

    # end of function _remove_object_type_of_organism_from_graph

    ############################################################################
    #### Functions allowing / disallowing introns for a specific Organism   ####
    ############################################################################

    def disallow_intron_in_organism(self,organism):
        """
        Remove the SpliceSites (not projected sites) from the DSCG and ASCG for
        this Organism identifier

        @type  organism: * (string)
        @param organism: Organism Identifier

        @rtype:  True 
        @return: True
        """
        if self._forced_3p_ends and self._forced_5p_ends:
            # this function is not applicable for splitted Interfaces
            # return a non-meaningfull True
            return True

        # disallow intron for this organism
        status = self._remove_object_type_of_organism_from_graph(
                organism,'SpliceDonor','_splicedonorgraph')
        status = self._remove_object_type_of_organism_from_graph(
                organism,'SpliceAcceptor','_spliceacceptorgraph')

        # set interface_is_intron status to False!
        self._interface_is_intron[organism] = False

        # return a non-meaningfull True
        return True

    # end of function disallow_intron_in_organism


    def force_intron_in_organism(self,organism):
        """
        Remove the ProjectedSpliceSites (not the True SpliceSites) from the DSCG
        and ASCG for this Organism identifier

        @type  organism: * (string)
        @param organism: Organism Identifier

        @rtype:  True 
        @return: True
        """
        if self._forced_3p_ends and self._forced_5p_ends:
            # this function is not applicable for splitted Interfaces
            # return a non-meaningfull True
            return True

        # force intron boundary in this organism
        status = self._remove_object_type_of_organism_from_graph(
                organism,'ProjectedSpliceDonor','_splicedonorgraph')
        status = self._remove_object_type_of_organism_from_graph(
                organism,'ProjectedSpliceAcceptor','_spliceacceptorgraph')

        # set interface_is_intron status to True!
        self._interface_is_intron[organism] = True

        # return a non-meaningfull True
        return True

    # end of function force_intron_in_organism


    def force_intron_in_organisms(self,organisms):
        """
        Remove the ProjectedSpliceSites (not the True SpliceSites) from the DSCG
        and ASCG for a list of Organism identifier

        @type  organisms: list
        @param organisms: list with Organism Identifiers

        @rtype:  True 
        @return: True
        """
        if self._forced_3p_ends and self._forced_5p_ends:
            # this function is not applicable for splitted Interfaces
            # return a non-meaningfull True
            return True

        # force intron boundaries for these organisms
        for organism in organisms:
            self.force_intron_in_organism(organism)

        # return a non-meaningfull True
        return True

    # end of function force_intron_in_organisms


    def allow_intron_in_organisms(self,organisms):
        """
        Remove the SpliceSites (not projected sites) from the DSCG and ASCG for
        all Organisms in the cbgIF except for those provided as a list

        @type  organisms: list
        @param organisms: list with Organism Identifiers for which the
                SpliceSites must NOT be removed

        @rtype:  True 
        @return: True

        @attention: used in case of an lsrCBG or other merging situation.
        """
        if self._forced_3p_ends and self._forced_5p_ends:
            # this function is not applicable for splitted Interfaces
            # return a non-meaningfull True
            return True

        for organism in self.acceptorCBG.organism_set():
            if organism not in organisms:
                self.disallow_intron_in_organism(organism)

        # return a non-meaningfull True
        return True

    # end of function allow_intron_in_organisms

    ############################################################################
    #### Functions that find the optimal ADSG and AASG in the               ####
    #### DSCG and ASCG, respectively. There are 2 options, of increasing    ####
    #### accuracy but at the cost of (large) increase of computation effort:####
    ####    (1) find_conserved_splice_sites()                               ####
    ####    (2) collection2alignedsites()                                   ####
    #### In easy solvable cases (the vast majority), the first function     ####
    #### will do.                                                           ####
    ############################################################################

    def find_conserved_splice_sites(self,edges=None):
        """ 
        """
        # reset result variables of optimal splice sites
        self._optimal_aligned_donor    = None
        self._optimal_aligned_acceptor = None

        ##################################################
        if self.verbose:
            # StopWatch time logger in verbose mode
            stw = StopWatch(name="conservedsites")
            stw.start()
        ##################################################

        # get number of expected/required organism nodes
        if not edges: edges = self.organism_set_size()-1

        # Only find splicesites between non-lsrCBGs and for CBGs that have
        # already splicedonor- and acceptor graphs constructed. It can occur
        # that this is performed, but simply no sites have been found ->
        # no object created, still None. This is verified by checking if the
        # graph objects exist or not.
        if self.donorCBG._short_name == 'CBG' and\
        self.acceptorCBG._short_name == 'CBG' and\
        self._splicedonorgraph and self._spliceacceptorgraph and\
        not (self._forced_3p_ends or self._forced_5p_ends):

            # find perfectly (dist=0nt, wt=1.0) aligned donor sites
            self._splicedonorgraph.alignedsites = sort_by_cumulative_score(
                [ conversion.SpliceSiteCollectionGraph2AlignedSpliceSiteGraph(
                    algsite,max_node_count=self.donorCBG.organism_set_size()
                    ) for algsite in\
                 self._splicedonorgraph.find_conserved_sites(edges=edges) ] )

            if self._splicedonorgraph.alignedsites:
                # conserved sites found -> set to optimal site
                self._optimal_aligned_donor =\
                        self._splicedonorgraph.alignedsites[0]

                ################################################################
                if self.verbose:
                    print stw.lap(), "donors:",
                    print len(self._splicedonorgraph.alignedsites),
                    print [ ds.phase() for ds in\
                            self._splicedonorgraph.alignedsites ],
                    print [ ds.cumulative_score() for ds in\
                            self._splicedonorgraph.alignedsites ]
                    print self._optimal_aligned_donor
                ################################################################

            else:
                pass
                ################################################################
                if self.verbose: print stw.lap(), "donors: 0"
                ################################################################

            # find perfectly (dist=0nt, wt=1.0) aligned acceptor sites
            self._spliceacceptorgraph.alignedsites = sort_by_cumulative_score(
                [ conversion.SpliceSiteCollectionGraph2AlignedSpliceSiteGraph(
                    algsite,max_node_count=self.acceptorCBG.organism_set_size()
                    ) for algsite in\
                self._spliceacceptorgraph.find_conserved_sites(edges=edges) ] )

            if self._spliceacceptorgraph.alignedsites:
                # conserved sites found -> set to optimal site
                self._optimal_aligned_acceptor =\
                        self._spliceacceptorgraph.alignedsites[0]

                ################################################################
                if self.verbose:
                    print self._optimal_aligned_acceptor
                    print stw.lap(), "acceps:",
                    print len(self._spliceacceptorgraph.alignedsites),
                    print [ asite.phase() for asite in\
                            self._spliceacceptorgraph.alignedsites ],
                    print [ asite.cumulative_score() for asite in\
                            self._spliceacceptorgraph.alignedsites ]
                ################################################################
            else:
                pass
                ################################################################
                if self.verbose: print stw.lap(), "acceps: 0"
                ################################################################
        else:
            # an interface with enforced ends due to a lsrCBG. No splice sites
            pass

    # end of function find_conserved_splice_sites

  
    def collection2alignedsites(self):
        """
        Perform full site alignment on currently present donors & acceptors in the Interface
        """
        self._spliceacceptorgraph.collection2alignedsites()
        self._splicedonorgraph.collection2alignedsites()
        if self._splicedonorgraph.alignedsites:
            self._optimal_aligned_donor =\
                    self._splicedonorgraph.alignedsites[0]
        if self._spliceacceptorgraph.alignedsites:
            self._optimal_aligned_acceptor =\
                    self._spliceacceptorgraph.alignedsites[0]
        optimal_phase = self._assign_optimal_phase()
        if optimal_phase != None:
            self._adjust_to_optimal_phase(optimal_phase,verbose=True)
        # fix gaps between projected splice sites
        self._fix_gap_between_projected_splice_sites()

    # end of function collection2alignedsites


    def _assign_optimal_phase(self):
        """
        Assign the optimal phase when the highest ADSG & AASG are of distinct
        phase, but near-optimal ASSG have phase argeement and are K(s-x) ASSGs.

        @attention: Helper function only for find_conserved_splice_sites()

        @rtype  phase: 0,1,2 or None
        @return phase: 0,1,2 or None
        """
        # check which phase is the highest scoring
        phase_scores = {0:0.0,1:0.0,2:0.0}
        phase_seen   = {0:0.0,1:0.0,2:0.0}
        for phase in phase_scores.keys():
            for dsg in self._splicedonorgraph.alignedsites:
                if dsg.__class__.__name__ == 'DonorSiteCollectionGraph':
                    break
                if dsg.phase() == phase:
                    phase_scores[phase]+=dsg.cumulative_score()
                    phase_seen[phase]+=1
                    break
            for asg in self._spliceacceptorgraph.alignedsites:
                if asg.__class__.__name__ == 'AcceptorSiteCollectionGraph':
                    break
                if asg.phase() == phase:
                    phase_scores[phase]+=asg.cumulative_score()
                    phase_seen[phase]+=1
                    break
        # get the maximal phase score and the corresponding phase
        max_score = max(phase_scores.values())
        optimal_phase = None
        for phase,phasescore in phase_scores.iteritems():
            if phasescore == max_score:
                if phase_seen[phase] == 2:
                    optimal_phase = phase
                    break
        # return the optimal phase; 0,1,2 or None
        return optimal_phase

    # end of function _assign_optimal_phase 


    def _adjust_to_optimal_phase(self,phase,verbose=False):
        """
        Pick the highest scoring ADSG & AASG of the assigned phase as the
        optimal ADSG & AASG

        @attention: Helper function only for find_conserved_splice_sites()
        @attention: Either the current ADSG or AASG must have the assigned phase

        @type  phase: 0,1,2 or None
        @param phase: 0,1,2 or None

        @type  verbose: Boolean
        @param verbose: in case True, print logging/debugging messages to STDOUT
        """
        if self._optimal_aligned_donor:
            phase_donor = self._optimal_aligned_donor.phase() 
        else:
            phase_donor = None
        if self._optimal_aligned_acceptor: 
            phase_acceptor = self._optimal_aligned_acceptor.phase()    
        else:
            phase_acceptor = None

        # now check phase agreement
        if phase_donor != phase and phase_acceptor != phase:
            ####################################################################
            if verbose:
                print "BOTH SITES DO NOT MATCH OPTIMAL PHASE!!!!!",
                print phase_donor, phase_acceptor, "phase:", phase
            ####################################################################
            pass

        elif phase_donor != phase:
            for ii in range(0,len(self._splicedonorgraph.alignedsites)):
                dsg = self._splicedonorgraph.alignedsites[ii]
                if dsg.phase() == phase:
                    # set this aligned site in the front of the list
                    self._splicedonorgraph.alignedsites.pop(ii)
                    self._splicedonorgraph.alignedsites.insert(0,dsg)
                    # and set as the optimal donor too
                    self._optimal_aligned_donor = dsg
                    ############################################################
                    if verbose: print "DONOR adjusted from %s to %s" % (
                            phase_donor, phase)
                    ############################################################
                    break

        elif phase_acceptor != phase:
            for ii in range(0,len(self._spliceacceptorgraph.alignedsites)):
                asg = self._spliceacceptorgraph.alignedsites[ii]
                if asg.phase() == phase:
                    # set this aligned site in the front of the list
                    self._spliceacceptorgraph.alignedsites.pop(ii)
                    self._spliceacceptorgraph.alignedsites.insert(0,asg)
                    # and set as the optimal acceptor too
                    self._optimal_aligned_acceptor = asg
                    ############################################################
                    if verbose: print "ACCEPTOR adjusted from %s to %s" % (
                            phase_acceptor, phase)
                    ############################################################
                    break
        else:
            # both AADG & ADSG are already of the assigned phase
            pass

    # end of function _adjust_to_optimal_phase

    ############################################################################
    #### Functions for creation of CodingBlockStart & CodingBlockEnd        ####
    #### objects at the interface between the CBGs. Used in case of         ####
    #### lsrCBGs and unusual cases of CBGs that are perfectly glued together####
    ############################################################################

    def _obtain_forced_ends(self):
        """
        Asign hard-set boundaries to the CBGs in the Interface in case of
        lsrCBGs and uncommen splitted CBGs

        @attention: ONLY call this function once from the __init__ function!
        """
        if self.donorCBG.IS_3P_SPLITTED and self.acceptorCBG.IS_5P_SPLITTED:
            # loop over the organisms
            omsrD = self.donorCBG.overall_minimal_spanning_range()
            omsrA = self.acceptorCBG.overall_minimal_spanning_range()
            for org in self.organism_set():
                if org in self.donorCBG.organism_set():
                    orf_of_org = self.donorCBG.get_orfs_of_graph(
                                    organism=org)[0]
                    node       = self.donorCBG.node_by_organism(org)
                    posAA      = max( omsrD[node] ) + 1
                    posDNA     = orf_of_org.proteinpos2dnapos(posAA)
                    cbgEnd     = CodingBlockEnd(posDNA)
                    # update _forced_3p_ends attribute
                    self._forced_3p_ends[org] = cbgEnd

                if org in self.acceptorCBG.organism_set():
                    orf_of_org = self.acceptorCBG.get_orfs_of_graph(
                                    organism=org)[0]
                    node       = self.acceptorCBG.node_by_organism(org)
                    posAA      = min( omsrA[node] )
                    posDNA     = orf_of_org.proteinpos2dnapos(posAA)
                    cbgSta     = CodingBlockStart(posDNA)
                    # update _forced_5p_ends attribute
                    self._forced_5p_ends[org] = cbgSta

            # set forced end attributes as well to the CBGs themselves
            self.donorCBG._forced_3p_ends = self._forced_3p_ends
            self.acceptorCBG._forced_5p_ends = self._forced_5p_ends

    # end of function _obtain_forced_ends


    def _create_forced_ends_for_non_introns(self,verbose=False):
        """
        Create forced CBG Sta & End objects when intron is not allowed and
        no Projected SpliceSites are available.
        
        @attention: USE WITH CARE; only use properly in the optimize() function
        
        @type  verbose: Boolean
        @param verbose: in case True, print logging/debugging messages to STDOUT
        
        @rtype:  Boolean
        @return: True or False, weather or not forced sites are created
        """
        SITES_ARE_CREATED = False

        #################################################################
        if verbose or self.verbose:
            print "# CREATE _create_forced_ends_for_non_introns(self):"
            print self
            print self._interface_is_intron 
        #################################################################

        # add CodingBlockStart & CodingBlockEnd for _interface_is_intron==False
        for org, status in self._interface_is_intron.iteritems():
            if status == False:
                inDonor  = org in self._splicedonorgraph.organism_set()
                inAccep  = org in self._spliceacceptorgraph.organism_set()
                ################################################################
                if verbose: print org, inDonor,inAccep,str(self)
                ################################################################
                if not inDonor or not inAccep:
                    # cbgIF for this organism identifier is NO intron,
                    # but there is no site in the ADSG or AASG. This can occur
                    # when intron projection was insuccesfull. Trust the
                    # cbgIF._interface_is_intron dict and create CBGend and
                    # CBGstart objects in stead if (Projected)SpliceSites
                    omsrDpos = max(
                        self.donorCBG.overall_minimal_spanning_range(
                            organism=org) ) + 1
                    omsrApos = min(
                        self.acceptorCBG.overall_minimal_spanning_range(
                            organism=org) )
                    omsraa   = (omsrDpos + omsrApos) / 2
                    theorf   = self.donorCBG.get_orfs_of_graph(organism=org)[0]
                    omsrdna  = theorf.proteinpos2dnapos(omsraa)
                    cbgEnd   = CodingBlockEnd(omsrdna)
                    cbgSta   = CodingBlockStart(omsrdna)
                    # update _forced_3p_ends attribute
                    self._forced_3p_ends[org] = cbgEnd
                    # update _forced_5p_ends attribute
                    self._forced_5p_ends[org] = cbgSta
                    SITES_ARE_CREATED = True
                    ############################################################
                    if verbose or self.verbose:
                        print "# NON-INTRON forced end created:",
                        print org, status, cbgEnd
                    ############################################################

        # return Boolean weather or not sites are created
        return SITES_ARE_CREATED

    # end of function _create_forced_ends_for_non_introns

    ############################################################################
    #### Functions for optimizing non-optimal cbgIFs                        ####
    ############################################################################

    def optimize(self,verbose=False):
        """
        If the CBGinterface is non-optimal, try to assign a better site.

        If after (quick) SpliceSite alignment (from DSCG and ASCG to the
        optimal ADSG and AASG, respectively) not an is_optimal() cbgIF
        is produced, the cbgIF.optimize() function will step-by-step (with
        increasing computational effort & decreasing likelihood of cbgIF
        structure (e.g. non-canonical sites, phase shift etc)) try to
        make the cbgIF at least is_acceptable()

        @type  verbose: Boolean
        @param verbose: in case True, print logging/debugging messages to STDOUT
        
        @rtype:  NoneBoolean
        @return: None (already optimal), True (optimized) or False (failed)
        """
        # cleanup self._optimizations list
        self._optimizations = []
        
        # if already optimal, no optimization possible!
        if self.is_optimal():
            ####################################################################
            if verbose: print "cbgIF.optimization() ended: is_optimal() == True"
            ####################################################################
            return None

        ####################################################################
        if verbose: print "cbgIF.optimization() start:", self
        ####################################################################

        # check if forced sites available -> no further improvement possible
        if self._forced_3p_ends or self._forced_5p_ends:
            ####################################################################
            if verbose: print "cbgIF.optimization() ended: forced_ends "
            ####################################################################
            return None

        # (1) first fix possible wrong projected sites
        self._fix_gap_between_projected_splice_sites()

        # append message for 1th optimization action to self._optimizations
        self._optimizations.append( 'fix_gap_between_projected_splice_sites' )
        ########################################################################
        if verbose: print "cbgIF.optimization() step:", self._optimizations[-1] 
        ########################################################################

        if self.is_acceptable():
            ####################################################################
            if verbose: print "cbgIF.optimization() ended: is_acceptable()"
            ####################################################################
            return True

        ########################################################################
        if verbose:
            print "before optimization:", self, "optimal", self.is_optimal(),
            print "acceptable", self.is_acceptable()
            print "before optimization:", self._optimal_aligned_donor
            print "before optimization:", self._optimal_aligned_acceptor
        ########################################################################

        # (2) do full sitealignment of all sites
        # align SpliceSites in ADSGs and AASGs and take best scoring ones
        self.collection2alignedsites()

        # append message for 2th optimization action to self._optimizations
        self._optimizations.append( 'collection2alignedsites' )
        ########################################################################
        if verbose: print "cbgIF.optimization() step:", self._optimizations[-1] 
        ########################################################################

        if self.is_acceptable():
            ####################################################################
            if verbose: print "cbgIF.optimization() ended: is_acceptable()"
            ####################################################################
            return True

        # (3) If some _interface_is_intron == False values exist, overrule
        # these by setting these to None (doubtfull intron). Next, try a
        # couple of cbgIF optimization steps and see if it solves the problem
        if False in self._interface_is_intron.values():
            corCbgIf = CodingBlockGraphInterface(self.donorCBG,self.acceptorCBG)
            for org,isintron in corCbgIf._interface_is_intron.iteritems():
                if isintron == False:
                    corCbgIf._interface_is_intron[org] = None
            # harvest splice sites and do
            # quick site alignment, perfectly conserved intron
            corCbgIf.harvest_donor_sites()
            corCbgIf.harvest_acceptor_sites()
            corCbgIf.find_conserved_splice_sites()
            if corCbgIf.is_acceptable():
                self._copy_interface_properties_from_otherinterface(corCbgIf)
                self._optimizations.append( 'adjusted_interface_is_intron 1' )
                return True

            # optimizetinyexoninterface might help too
            corCbgIf.optimizetinyexoninterface()
            if corCbgIf.is_acceptable():
                self._copy_interface_properties_from_otherinterface(corCbgIf)
                self._optimizations.append( 'adjusted_interface_is_intron 2' )
                return True

            # full site alignment as a final option here
            corCbgIf.collection2alignedsites()
            if corCbgIf.is_acceptable():
                self._copy_interface_properties_from_otherinterface(corCbgIf)
                self._optimizations.append( 'adjusted_interface_is_intron 3' )
                return True


        # (4) add CodingBlockStart & CodingBlockEnd for
        # _interface_is_intron == False cbgIFs. This is rarely
        # the case, but it must be checked here
        SITES_ARE_CREATED = self._create_forced_ends_for_non_introns(
            verbose=verbose)
        if verbose: print "non-introns2forcedends:", SITES_ARE_CREATED
        # TODO: if SITES_ARE_CREATED, do a sequence error check!

        # append information about this optimization step too 
        self._optimizations.append( 'adjusted_interface_is_intron' )
        ########################################################################
        if verbose: print "OPTIMIZ. STEP PERFORMED:", self._optimizations[-1] 
        ########################################################################

        if self.is_acceptable():
            ####################################################################
            if verbose: print "cbgIF.optimization() ended: is_acceptable()"
            ####################################################################
            return True

        ########################################################################
        # (5) rescan for splice site but now with a bigger triplet distance
        # for aligned sites and a larger area to be scanned
        ########################################################################
        self.harvest_donor_sites(
                aligned_donor_max_triplet_distance=\
                OPTIMIZE_CBGIF_DONOR_MAX_TRIPLET_DISTANCE,
                enlarge_3p_boundary_by=\
                OPTIMIZE_CBGIF_DONOR_ENLARGE_3P_BOUNDARY_BY,
                )
        self.harvest_acceptor_sites(
                aligned_acceptor_max_triplet_distance=\
                OPTIMIZE_CBGIF_ACCEPTOR_MAX_TRIPLET_DISTANCE,
                enlarge_5p_boundary_by=\
                OPTIMIZE_CBGIF_ACCEPTOR_ENLARGE_5P_BOUNDARY_BY,
                )
        # align SpliceSites in ADSGs and AASGs and take best scoring ones
        self.collection2alignedsites()
        
        # append message for 3th optimization action to self._optimizations
        self._optimizations.append('increase triplet distance & elegiable area')
        ########################################################################
        if verbose: print "cbgIF.optimization() step:", self._optimizations[-1] 
        ########################################################################
        
        if self.is_acceptable():
            ####################################################################
            if verbose: print "cbgIF.optimization() ended: is_acceptable()"
            ####################################################################
            return True
        
        
        ########################################################################
        # (6) allow phase shift & non-canonical donors
        ########################################################################
        self.harvest_donor_sites(
                allow_non_canonical=True,
                allow_phase_shift=True
                )
        self.harvest_acceptor_sites(
                allow_phase_shift=True
                )
        # align SpliceSites in ADSGs and AASGs and take best scoring ones
        self.collection2alignedsites()

        # append message for 4th optimization action to self._optimizations
        self._optimizations.append(
            'non-canonical donors & splice site phase shift' )
        ########################################################################
        if verbose: print "cbgIF.optimization() step:", self._optimizations[-1] 
        ########################################################################

        if self.is_acceptable():
            ####################################################################
            if verbose: print "cbgIF.optimization() ended: is_acceptable()"
            ####################################################################
            return True

        ########################################################################
        if verbose:
            print "cbgIF.optimization() FAILED"
            print self
            print self._splicedonorgraph
            print self._spliceacceptorgraph
            print self._optimal_aligned_donor
            print self._optimal_aligned_acceptor 
            #for site in self._spliceacceptorgraph.alignedsites:
            #    print "@@@@", site
            #for site in self._splicedonorgraph.alignedsites:
            #    print "$$$$", site
        ########################################################################

        # if this point is reached: a very crapy interface!
        return False

    # end of function optimize


    def _copy_interface_properties_from_otherinterface(self,othercbgif):
        """
        Copy all the main attributes from another CbgIF

        @type  othercbgif: CodingBlockGraphInterface 
        @param othercbgif: CodingBlockGraphInterface 
        """
        self._HAS_INTRONS_PROJECTED     = othercbgif._HAS_INTRONS_PROJECTED
        self._USE_ENTROPY               = othercbgif._USE_ENTROPY
        self._forced_3p_ends            = othercbgif._forced_3p_ends
        self._forced_5p_ends            = othercbgif._forced_5p_ends
        self._optimal_aligned_donor     = othercbgif._optimal_aligned_donor
        self._optimal_aligned_acceptor  = othercbgif._optimal_aligned_acceptor
        self._splicedonorgraph          = othercbgif._splicedonorgraph
        self._spliceacceptorgraph       = othercbgif._spliceacceptorgraph
        self._interface_is_intron       = othercbgif._interface_is_intron

    # end of function _copy_interface_properties_from_otherinterface


    def _cleanup_donorsites_outside_omsr_region(self,
        omsr_3p_nt_offset=0,verbose=False):
        """
        Remove all (Projected)SpliceSite nodes 'omsr_3p_nt_offset' away from
        the 3' OMSR of the DonorSiteCollectionGraph

        @attention: function is used (only) in optimizetinyexoninterface

        @type  omsr_3p_nt_offset: integer
        @param omsr_3p_nt_offset: offset to increase the OMSR with on the 3'side
        
        @type  verbose: Boolean
        @param verbose: in case True, print logging/debugging messages to STDOUT
        
        @rtype:  integer
        @return: number of (Projected)DonorSites that are deleted
        """
        # get OMSR, and redefine the consideredsplicesiterange() on the
        # 3' side of the OMSR based on variable omsr_3p_nt_offset
        donoromsr = self.donorCBG.overall_minimal_spanning_range()
        nodes_removed_cnt=0
        for org in self.organism_set():
            if org not in self.donorCBG.organism_set(): continue
            node = self.donorCBG.node_by_organism(org)
            try:
                donorrange =\
                    self._splicedonorgraph.get_consideredsplicesiterange(org)
            except:
                donorrange = ( min(donoromsr[node])*3 , None )
            # take the original consideredsplicesiterange() on the 5' side
            newstart = donorrange[0]
            # enlarge 3p donor OMSR range by applied variable
            newend   = max(donoromsr[node])*3 + omsr_3p_nt_offset
            # create new consideredsplicesiterange()
            newrange = range(newstart,newend+1)
            # update considered splice donor range 
            self._splicedonorgraph.set_consideredsplicesiterange(
                    org,newstart,newend
                    )
            ################################################################
            if verbose:
                print org, donorrange, donorrange[1]-donorrange[0],
                print (newstart, newend), (newend-newstart)
            ################################################################
            delnodes = []
            for donornode, donorobj in\
            self._splicedonorgraph._node_object.iteritems():
                if self._splicedonorgraph.organism_by_node(donornode)==org:
                    if donorobj.pos not in newrange:
                        delnodes.append( donornode )
            # remove the nodes listed in list 'delnodes'
            for donornode in delnodes:
                self._splicedonorgraph.del_node( donornode )
            nodes_removed_cnt += len(delnodes)
            # run _update_after_changes function
            self._splicedonorgraph._update_after_changes()
            
        # return number of removed nodes
        return nodes_removed_cnt
    
    # end of function _cleanup_donorsites_outside_omsr_region


    def _cleanup_acceptorsites_outside_omsr_region(self,
        omsr_5p_nt_offset=0,verbose=False):
        """
        Remove all (Projected)SpliceSite nodes 'omsr_5p_nt_offset' away from
        the 5' OMSR of the AcceptorSiteCollectionGraph

        @attention: function is used (only) in optimizetinyexoninterface

        @type  omsr_5p_nt_offset: integer
        @param omsr_5p_nt_offset: offset to increase the OMSR with on the 5'side
        
        @type  verbose: Boolean
        @param verbose: in case True, print logging/debugging messages to STDOUT
        
        @rtype:  integer
        @return: number of (Projected)AcceptorSites that are deleted
        """
        # get OMSR, and redefine the consideredsplicesiterange() on the
        # 5' side of the OMSR based on variable omsr_5p_nt_offset
        acceptoromsr = self.acceptorCBG.overall_minimal_spanning_range()
        nodes_removed_cnt=0
        for org in self.organism_set():
            if org not in self.acceptorCBG.organism_set(): continue
            node = self.acceptorCBG.node_by_organism(org)
            try:
                acceptorrange =\
                    self._spliceacceptorgraph.get_consideredsplicesiterange(org)
            except:
                acceptorrange = ( None, max(acceptoromsr[node])*3 )

            # take the original consideredsplicesiterange() on the 3' side
            newend   = acceptorrange[1]
            # enlarge 5p donor OMSR range by applied variable
            newstart = min(acceptoromsr[node])*3 - omsr_5p_nt_offset
            # create new consideredsplicesiterange()
            newrange = range(newstart,newend+1)
            # update considered splice donor range 
            self._spliceacceptorgraph.set_consideredsplicesiterange(
                org,newstart,newend)
            ################################################################
            if verbose:
                print org, acceptorrange, acceptorrange[1]-acceptorrange[0],
                print (newstart, newend), (newend-newstart)
            ################################################################
            delnodes = []
            for acceptornode, acceptorobj in\
            self._spliceacceptorgraph._node_object.iteritems():
                if self._spliceacceptorgraph.organism_by_node(acceptornode)==org:
                    if acceptorobj.pos not in newrange:
                        delnodes.append( acceptornode )
            # remove the nodes listed in list 'delnodes'
            for acceptornode in delnodes:
                self._spliceacceptorgraph.del_node( acceptornode )
            nodes_removed_cnt += len(delnodes)
            # run _update_after_changes function
            self._spliceacceptorgraph._update_after_changes()
          
        # return number of removed nodes
        return nodes_removed_cnt
    
    # end of function _cleanup_acceptorsites_outside_omsr_region


    def optimizetinyexoninterface(self,
        tinyexon_max_aa_size=OPTIMIZE_CBGIF_TINYEXON_MAX_AA_SIZE,
        verbose=False):
        """         
        Optimization of the cbgIF especially for small/tiny exons/CBGs
        
        @type  tinyexon_max_aa_size: integer
        @param tinyexon_max_aa_size:
        
        @type  verbose: Boolean
        @param verbose: in case True, print logging/debugging messages to STDOUT
        
        @rtype:  NoneBoolean
        @return: None (already optimal), True (optimized) or False (failed)
        """
        # cleanup self._optimizations list
        self._optimizations = []

        # if already optimal, no optimization possible!
        if self.is_optimal():
            ####################################################################
            if verbose: print "cbgIF.tyexonoptimiz() ended: is_optimal()=True"
            ####################################################################
            return None

        # check if forced sites available -> no further improvement possible
        if self._forced_3p_ends or self._forced_5p_ends:
            ####################################################################
            if verbose: print "cbgIF.tyexonoptimiz() ended: forced_ends "
            ####################################################################
            return None

        # (1) first fix possible wrong projected sites
        self._fix_gap_between_projected_splice_sites()

        # append message for 1th optimization action to self._optimizations
        self._optimizations.append( 'fix_gap_between_projected_splice_sites' )
        ########################################################################
        if verbose: print "cbgIF.tyexonoptimiz() step:", self._optimizations[-1] 
        ########################################################################
        
        if self.is_acceptable():
            ####################################################################
            if verbose: print "cbgIF.tyexonoptimiz() ended: is_acceptable()"
            ####################################################################
            return True
        
        # check if there is a tinyexonCBG in this interface
        # if not -> do not perform this function and return False
        if not (self.donorCBG.omsrlength() <= tinyexon_max_aa_size or\
        tinyexon_max_aa_size <= tinyexon_max_aa_size):
            ####################################################################
            if verbose: print "cbgIF.tyexonoptimiz() ended: not applicable"
            ####################################################################
            return False

        # (2) BRUTE FORCE: allow phase shift & non-canonical donors
        # Then, delimit the cbgIF suitable area for splice sites
        # and only then do site alignment into ADSG & AASG
        self.harvest_donor_sites(
                allow_non_canonical=True,
                allow_phase_shift=True)
        self.harvest_acceptor_sites(allow_phase_shift=True)
        self._optimizations.append(
            'non-canonical donors & splice site phase shift' )

        if self.donorCBG.omsrlength() <= tinyexon_max_aa_size:
            # donorCBG isa tinyexon -> shorten its consideredsplicesiterange()
            deleted = self._cleanup_donorsites_outside_omsr_region(
                omsr_3p_nt_offset=\
                OPTIMIZE_CBGIF_TINYEXON_DONOR_OMSR_3P_NT_OFFSET)
            # append message for 2th optimization action to self._optimizations
            self._optimizations.append(
                'tinyexon donor elegiable range limitation' )
            
        if self.acceptorCBG.omsrlength() <= tinyexon_max_aa_size:
            # accepCBG isa tinyexon -> shorten its consideredsplicesiterange()
            deleted = self._cleanup_acceptorsites_outside_omsr_region(
                omsr_5p_nt_offset=\
                OPTIMIZE_CBGIF_TINYEXON_ACCEPTOR_OMSR_5P_NT_OFFSET)
            # append message for 2th optimization action to self._optimizations
            self._optimizations.append(
                'tinyexon acceptor elegiable range limitation' )
            
        ########################################################################
        if verbose: print "cbgIF.tyexonoptimiz() step:", self._optimizations[-1] 
        ########################################################################
        # find conserved splice sites
        self.find_conserved_splice_sites()
        
        if self.is_acceptable():
            ####################################################################
            if verbose: print "cbgIF.tyexonoptimiz() ended: is_acceptable()"
            ####################################################################
            return True
        
        # do full sitealignment of all sites
        self.collection2alignedsites()
        
        # Done. Check if we have an acceptable interface now
        if self.is_acceptable(): return True
        else:                    return False

    # end of function optimizetinyexoninterface


    def is_seemless_intron_projection(self,organism,aa_offset=0):
        """
        Are projected sites in this interface seemlessly fitting together?
        @type  organism: * (string)
        @param organism: Organism Identifier or None
        @type  aa_offset: integer
        @param aa_offset: maximal aa offset between projected sites
        @rtype:  NoneBoolean
        @return: None (not applicable), True (yes) or False (no)
        """
        if organism not in self.organism_set():
            raise OrganismNotPresentInGraph
        donor    = self.get_assigned_donor_site(organism)
        acceptor = self.get_assigned_acceptor_site(organism)
        if not donor or not acceptor:
            # irrelevant; not a NON-seemless projection ;-)
            return None
        # check splice sites classes
        if donor.__class__.__name__ == 'ProjectedSpliceDonor' and\
        acceptor.__class__.__name__ == 'ProjectedSpliceAcceptor':
            if abs(acceptor.pos - donor.pos)/3 <= aa_offset:
                return True
            else:
                print "NOT SEEMLESS:", donor
                print "NOT SEEMLESS:", acceptor
                return False
        else:
            # irrelevant; not a NON-seemless projection ;-)
            return None
        
    # end of function is_seemless_intron_projection


    def _fix_gap_between_projected_splice_sites(self,verbose=False):
        """
        Fix gaps that are the result from ProjectedSpliceSite offsets

        @type  verbose: Boolean
        @param verbose: in case True, print logging/debugging messages to STDOUT

        @rtype:  Boolean
        @return: True or False, weather or not the `missing` projection was created
        """
        if self._optimal_aligned_donor and self._optimal_aligned_acceptor:
            if self._optimal_aligned_donor.__class__.__name__ == 'DonorSiteCollectionGraph':
                return False
            if self._optimal_aligned_acceptor.__class__.__name__ == 'AcceptorSiteCollectionGraph':
                return False
            # okay, at this stage we have aligned sites. Now check for phase identity
            if self._optimal_aligned_donor.phase() != self._optimal_aligned_acceptor.phase():
                return False
            # loop over the organisms and get the sites from the graph
            is_adjusted = False

            for org in self.organism_set():
                try:
                    donor  = self._optimal_aligned_donor.get_organism_objects(org)[0]
                except OrganismNotPresentInGraph:
                    donor = None
                try:
                    acceptor = self._optimal_aligned_acceptor.get_organism_objects(org)[0]
                except OrganismNotPresentInGraph:
                    acceptor = None
                # check if both donor and acceptor sites are available
                if not donor or not acceptor:
                    # check if we can create the `missing` projected site
                    is_created = self._create_missing_projected_site(org,donor,acceptor)
                    if not is_created:
                        continue
                    else:
                        # projected donor or acceptor site is created; get them again
                        donor    = self._optimal_aligned_donor.get_organism_objects(org)[0]
                        acceptor = self._optimal_aligned_acceptor.get_organism_objects(org)[0]
                else:
                    pass

                # check which type the splice sites are; here, we need only ProjectedSites
                dclass = donor.__class__.__name__
                aclass = acceptor.__class__.__name__
                if not (dclass == 'ProjectedSpliceDonor' and aclass == 'ProjectedSpliceAcceptor'):
                    continue
                # if here, both projected sites. Check if their positions are identical
                # TODO: what if the gap between the projections is HUGE?
                # a gap of 15nt (max projection distance used by default is 5aa) can occur
                if donor.pos != acceptor.pos:
                    # hard-adjust the projected acceptor to the position of the donor
                    acceptor.pos   = donor.pos
                    acceptor.end   = acceptor.pos
                    acceptor.start = acceptor.end-2
                    acceptor.phase = donor.phase # adapt the phase as well; phases can be inequal!
                    ############################################################
                    if verbose:
                        print "ProjectedSites positions matched!",
                        print org, donor.pos, acceptor.pos, acceptor
                    ############################################################
                    is_adjusted = True
            else:
                return is_adjusted
        else:
            return False

    # end of function _fix_gap_between_projected_splice_sites 


    def _create_missing_projected_site(self,organism,donor,acceptor,
        verbose=False):
        """
        Create a `missing` ProjectedSite that fixes an otherwise
        partial AlignedSpliceSiteGraph

        @attention: USE WITH CARE! In fact, not not use except as implemented in
                    self._fix_gap_between_projected_splice_sites()

        @type  organism: * (string)
        @param organism: Organism Identifier or None

        @type  donor: ProjectedSpliceDonor or None 
        @param donor: ProjectedSpliceDonor or None 

        @type  acceptor: ProjectedSpliceAcceptor or None
        @param acceptor: ProjectedSpliceAcceptor or None
    
        @type  verbose: Boolean
        @param verbose: print logging/debugging messages to STDOUT when True
        
        @rtype:  Boolean
        @return: True when a `missing` projection was created, False otherwise
        """
        if not donor and not acceptor:
            return False
        # only perform this function when one of both objects is a ProjectedSite
        if not donor and acceptor.__class__.__name__ !=\
        'ProjectedSpliceAcceptor':
            return False
        if not acceptor and donor.__class__.__name__ !=\
        'ProjectedSpliceDonor':
            return False
        # only perform this function when the other site is
        # `complete` ( a K(s) ADSG or AASG )
        if not donor and self._optimal_aligned_acceptor.organism_set_size() !=\
        self.organism_set_size():
            return False
        if not acceptor and self._optimal_aligned_donor.organism_set_size() !=\
        self.organism_set_size():
            return False

        # if here, then we have either a missing projected donor or acceptor
        if not donor:
            # create projected donor with knowledge from the acceptor
            dummydonor = SpliceDonor(acceptor.pos,"XX",
                    phase       = acceptor.phase,
                    strand      = acceptor.strand,
                    pattern_offset=(0,0)
                    )
            projdonor  = ProjectedSpliceDonor(acceptor.pos,"XX",dummydonor,
                    distance    = acceptor.distance,
                    entropy     = acceptor.entropy,
                    pssm_score  = acceptor.pssm_score,
                    cases       = acceptor.cases
                    )
            ####################################################################
            if verbose:
                print "PROJDONOR CREATED::", projdonor, acceptor
                print "before creation:", self._optimal_aligned_donor
            ####################################################################

            # now the tiresome part: add to the _optimal_aligned_donor graph
            orfobj  = self.donorCBG.get_orfs_of_graph(organism=organism)[0]
            dsqNode = ( organism, orfobj.id, projdonor.pos )
            nodes   = self._optimal_aligned_donor.get_nodes()
            # add this new node to the graph
            self._optimal_aligned_donor.add_node_and_object(dsqNode,projdonor)
            # create edges between the new node and all the other ones
            for node in nodes:
                self._optimal_aligned_donor.add_edge(dsqNode,node,wt=1.0)
                self._optimal_aligned_donor._edge_binary_entropies[(dsqNode,node)] = (projdonor.entropy,projdonor.entropy)
                self._optimal_aligned_donor._edge_binary_entropies[(node,dsqNode)] = (projdonor.entropy,projdonor.entropy)
            if verbose: print "after creation:", self._optimal_aligned_donor
            # return status True for added object
            return True

        if not acceptor:
            # create projected acceptor with knowledge from the donor
            dummyacceptor = SpliceAcceptor(donor.pos,"XX",
                    phase       = donor.phase,
                    strand      = donor.strand,
                    pattern_offset=(0,0)
                    )
            projacceptor  = ProjectedSpliceAcceptor(donor.pos,"XX",dummyacceptor,
                    distance    = donor.distance,
                    entropy     = donor.entropy,
                    pssm_score  = donor.pssm_score,
                    cases       = donor.cases
                    )
            ####################################################################
            if verbose:
                print "PROJACCEPTOR CREATED::", projacceptor, donor
                print "before creation:", self._optimal_aligned_acceptor
            ####################################################################

            # now the tiresome part: add to the _optimal_aligned_acceptor graph
            orfobj  = self.acceptorCBG.get_orfs_of_graph(organism=organism)[0]
            asqNode = ( organism, orfobj.id, projacceptor.pos )
            nodes   = self._optimal_aligned_acceptor.get_nodes()
            # add this new node to the graph
            self._optimal_aligned_acceptor.add_node_and_object(asqNode,projacceptor)
            # create edges between the new node and all the other ones
            for node in nodes:
                self._optimal_aligned_acceptor.add_edge(asqNode,node,wt=1.0)
                self._optimal_aligned_acceptor._edge_binary_entropies[(asqNode,node)] = (projacceptor.entropy,projacceptor.entropy)
                self._optimal_aligned_acceptor._edge_binary_entropies[(node,asqNode)] = (projacceptor.entropy,projacceptor.entropy)
            if verbose: print "after creation:", self._optimal_aligned_acceptor
            # return status True for added object
            return True

    # end of function _create_missing_projected_site

    ########################################################################
    #### Functions for accessing/obtaining *the* donor and acceptor     ####
    #### sites for an Organism Identifier, labelled with a succes status####
    ########################################################################

    def get_acceptor_object(self,organism):
        """
        Get the object that confines the 5' side of the cbgIF of this organism

        @type  organism: * (string)
        @param organism: Organism Identifier

        @attention: lsrCBGs in the cbgIF will return ( CodingBlockStart, True )
        @attention: incomplete/empty cbgIFs will return ( None, False )

        @rtype:  tuple
        @return: tuple of format ( object, Boolean )
                    object:     AcceptorSite, ProjectedAcceptorSite,
                                CodingBlockStart or None
                    boolean:    True when object is part of the optimal AASG
                                False when None or any other discrepancy
        """
        if organism not in self.organism_set():
            raise OrganismNotPresentInGraph

        if self._forced_5p_ends and self._forced_5p_ends.has_key(organism):
            # LowSimilarityRegionCodingBlockGraph; fixed sites to mask
            # a false inframe intron or: a forced Non-intron that failed
            # alignment and is overruled by forced ends
            boundary = self._forced_5p_ends[organism]
            return ( boundary, True )
        elif self._optimal_aligned_acceptor and\
        self._optimal_aligned_acceptor.__class__.__name__ ==\
        'AcceptorSiteCollectionGraph':
            if organism in self._optimal_aligned_acceptor.organism_set():
                acceptor = self._spliceacceptorgraph.get_optimal_single_site(
                    organism)
                return ( acceptor, False )
            else:
                # hmm.. .we can not find an eligiable acceptor site
                return ( None, False )
        elif self._optimal_aligned_acceptor and\
        organism in self._optimal_aligned_acceptor.organism_set():
            # this is the most regular case:
            # a nicely assigned and aligned splice site object in the AASG
            acceptor = self._optimal_aligned_acceptor.get_organism_objects(
                organism)[0]
            return ( acceptor, True )
        elif self._spliceacceptorgraph:
            if organism in self._spliceacceptorgraph.organism_set():
                # return the 'optimal' site
                acceptor = self._spliceacceptorgraph.get_optimal_single_site(
                    organism)
                return ( acceptor, False )
            else:
                # hmm.. .we can not find any eligiable acceptor site
                return ( None, False )
        else:
            # hmm.. .we can not find an eligiable acceptor site
            return ( None, False )

    # end of function get_acceptor_object


    def get_donor_object(self,organism):
        """
        Get the object that confines the 3' side of the cbgIF of this organism

        @type  organism: * (string)
        @param organism: Organism Identifier

        @attention: lsrCBGs in the cbgIF will return ( CodingBlockEnd, True )
        @attention: incomplete/empty cbgIFs will return ( None, False )

        @rtype:  tuple
        @return: tuple of format ( object, Boolean )
                    object:     DonorSite, ProjectedDonorSite,
                                CodingBlockEnd or None
                    boolean:    True when object is part of the optimal ADSG
                                False when None or any other discrepancy
        """
        if organism not in self.organism_set():
            raise OrganismNotPresentInGraph

        if self._forced_3p_ends and self._forced_3p_ends.has_key(organism):
            # LowSimilarityRegionCodingBlockGraph; fixed sites to mask
            # a false inframe intron or: a forced Non-intron that failed
            # alignment and is overruled by forced ends
            boundary = self._forced_3p_ends[organism]
            return ( boundary, True )
        elif self._optimal_aligned_donor and\
        self._optimal_aligned_donor.__class__.__name__ ==\
        'DonorSiteCollectionGraph':
            if organism in self._optimal_aligned_donor.organism_set():
                donor = self._splicedonorgraph.get_optimal_single_site(organism)
                return ( donor, False )
            else:
                # hmm.. .we can not find an eligiable donor site
                return ( None, False )
        elif self._optimal_aligned_donor and organism in\
        self._optimal_aligned_donor.organism_set():
            # this is the most regular case: a nicely assigned and
            # aligned splice site object in the ADSG
            donor = self._optimal_aligned_donor.get_organism_objects(organism)[0]
            return ( donor, True )
        elif self._splicedonorgraph:
            if organism in self._splicedonorgraph.organism_set():
                # return the 'optimal' site
                donor = self._splicedonorgraph.get_optimal_single_site(organism)
                return ( donor, False )
            else:
                # hmm.. .we can not find any eligiable donor site
                return ( None, False )
        else:
            # hmm...we can not find an eligiable donor site
            return ( None, False )

    # end of function get_donor_object


    def get_optimal_donor(self,organism,phase=None):
        """
        Get the optimal donor site for an organism in this CBGinterface

        @type  organism: * (string)
        @param organism: Organism Identifier

        @type  phase: 0,1,2 or None
        @param phase: when phase is applied, the optimal site with that phase
                      is returned

        @attention: *optimal* is a combination of PSSM score and positioning
                    towards the OMSR of the 3' end of the donorCBG
        @attention: use self.get_highest_pssm_donor(organism) for the highest
                    PSSM scoring DonorSite
        @attention: lsrCBGs in the CodingBlockInterface will return None
        @attention: incomplete/empty CodingBlockInterfaces will return None

        @rtype:  DonorSite (or None)
        @return: optimal DonorSite in the DSCG of the applied Organism
        """
        try:
            omsr = self.donorCBG.overall_minimal_spanning_range(
                organism=organism)
            return self._splicedonorgraph.get_optimal_single_site(
                organism,use_entropy_around_omsr=max(omsr))
        except AttributeError:
            return None
        except:
            raise "UnExpectedException, no CBGinterface.get_optimal_donor()"

    # end of function get_optimal_donor


    def get_highest_pssm_donor(self,organism,phase=None):
        """
        Get the highest PSSM DonorSite from the CBGinterface for an organism
        
        @type  organism: * (string)
        @param organism: Organism Identifier

        @type  phase: 0,1,2 or None
        @param phase: when phase is applied, the optimal site with that phase
                      is returned

        @attention: lsrCBGs in the CodingBlockInterface will return None
        @attention: incomplete/empty CodingBlockInterfaces will return None

        @rtype:  DonorSite (or None)
        @return: highest PSSM DonorSite in the DSCG of the applied Organism
        """
        try:
            return self._splicedonorgraph.get_optimal_single_site(organism)
        except AttributeError:
            return None
        except:
            raise "UnExpectedException, no CBGinterface.get_optimal_donor()"

    # end of function get_optimal_donor


    def get_optimal_acceptor(self,organism,phase=None):
        """
        Get the optimal acceptor site for an organism in this CBGinterface

        @type  organism: * (string)
        @param organism: Organism Identifier

        @type  phase: 0,1,2 or None
        @param phase: when phase is applied, the optimal site with that phase
                      is returned

        @attention: *optimal* is a combination of PSSM score and positioning
                    towards the OMSR of the 5' end of the acceptorCBG
        @attention: use self.get_highest_pssm_acceptor(organism) for the highest
                    PSSM scoring AcceptorSite
        @attention: lsrCBGs in the CodingBlockInterface will return None
        @attention: incomplete/empty CodingBlockInterfaces will return None

        @rtype:  AcceptorSite (or None)
        @return: optimal AcceptorSite in the ASCG of the applied Organism
        """
        try:
            omsr = self.acceptorCBG.overall_minimal_spanning_range(
                organism=organism)
            return self._spliceacceptorgraph.get_optimal_single_site(
                organism,use_entropy_around_omsr=min(omsr))
        except AttributeError:
            return None
        except:
            raise "UnExpectedException, no CBGinterface.get_optimal_aceptor()"

    # end of function get_optimal_acceptor


    def get_highest_pssm_acceptor(self,organism,phase=None):
        """
        Get the highest PSSM AcceptorSite from the CBGinterface for an organism
        
        @type  organism: * (string)
        @param organism: Organism Identifier

        @type  phase: 0,1,2 or None
        @param phase: when phase is applied, the optimal site with that phase
                      is returned

        @attention: lsrCBGs in the CodingBlockInterface will return None
        @attention: incomplete/empty CodingBlockInterfaces will return None

        @rtype:  AcceptorSite (or None)
        @return: highest PSSM AcceptorSite in the ASCG of the applied Organism
        """
        try:
            return self._spliceacceptorgraph.get_optimal_single_site(organism)
        except AttributeError:
            return None
        except:
            raise "UnExpectedException, no CBGinterface.get_optimal_aceptor()"

    # end of function get_optimal_acceptor

    ########################################################################
    #### Functions for improving the ADSG and AASG by swapping with     ####
    #### higher PSSM scoring SpliceSites in the near vicinity.          ####
    ########################################################################

    def is_donor_improveable(self,organism=None):
        """
        Are the DonorSite(s), given the optimal ADSG, improvable?

        @type  organism: * (string)
        @param organism: Organism Identifier or None

        @attention: the check is limited to the phase of the optimal ADSG !!
        @attention: lsrCBGs in the CodingBlockInterface will return None
        @attention: incomplete/empty CodingBlockInterfaces will return None
        @attention: when Organism identifier is given, the check is
                    performed only for the respective organism

        @rtype:  NoneBoolean
        @return: True  when a higher PSSM scoring site of the same phase exists
                 False when this is not the case
                 None  when this check is not applicable
        """
        if not self._optimal_aligned_donor or\
        self._optimal_aligned_donor.__class__.__name__ ==\
        'DonorSiteCollectionGraph':
            return None

        if organism and organism not in\
        self._optimal_aligned_donor.organism_set():
            # organism requested for but not in optimal donor 
            return None
        
        # check for each organism/gene if there is a higher PSSM scoring
        # site available OF THE GIVEN PHASE OF THE ADSG!!
        data = []
        for org in self._optimal_aligned_donor.organism_set():
            if organism and org != organism: continue
            answer = False
            current_site = self._optimal_aligned_donor.get_organism_objects(
                org)[0]
            for site in self._splicedonorgraph.get_organism_objects(org):
                if site.pos == current_site.pos: continue
                if site.phase != current_site.phase: continue
                if site.pssm_score > current_site.pssm_score:
                    answer = True
                    break
            data.append(answer)
            
        if True in data: return True
        else:            return False

    # end of function is_donor_improveable


    def is_acceptor_improveable(self,organism=None):
        """
        Are the AcceptorSite(s), given the optimal AASG, improvable?

        @type  organism: * (string)
        @param organism: Organism Identifier or None

        @attention: the check is limited to the phase of the optimal AASG !!
        @attention: lsrCBGs in the CodingBlockInterface will return None
        @attention: incomplete/empty CodingBlockInterfaces will return None
        @attention: when Organism identifier is given, the check is
                    performed only for the respective organism

        @rtype:  NoneBoolean
        @return: True  when a higher PSSM scoring site of the same phase exists
                 False when this is not the case
                 None  when this check is not applicable
        """
        if not self._optimal_aligned_acceptor or\
        self._optimal_aligned_acceptor.__class__.__name__ ==\
        'AcceptorSiteCollectionGraph':
            return None

        if organism and organism not in\
        self._optimal_aligned_acceptor.organism_set():
            # organism requested for but not in optimal acceptor
            return None

        # check for each organism/gene if there is a higher PSSM scoring
        # site available OF THE GIVEN PHASE OF THE AASG!!
        data = []
        for org in self._optimal_aligned_acceptor.organism_set():
            if organism and org != organism: continue
            answer = False
            current_site = self._optimal_aligned_acceptor.get_organism_objects(
                org)[0]
            for site in self._spliceacceptorgraph.get_organism_objects(org):
                if site.pos == current_site.pos: continue
                if site.phase != current_site.phase: continue
                if site.pssm_score > current_site.pssm_score:
                    answer = True
                    break
            data.append(answer)
            
        if True in data: return True
        else:            return False

    # end of function is_acceptor_improveable


    def propose_improved_donor(self,organism,
        max_aa_dist=PROPOSE_IMPROVED_DONOR_MAX_AA_DIST,
        min_pssm_ratio=PROPOSE_IMPROVED_DONOR_MIN_PSSM_RATIO):
        """
        Return - if there is - a higher scoring donor of the same phase AND
        in the vicinity of the currently optimal aligned donor in ADSG
        
        @type  organism: * (string)
        @param organism: Organism Identifier
        
        @type  max_aa_dist: integer
        @param max_aa_dist: maximal distance in AA positions from current site
        
        @type  min_pssm_ratio: float
        @param min_pssm_ratio: relative gain in PSSM score of improved site
        
        @attention: min_pssm_ratio must be greater than 1.0 !
        @attention: max_aa_dist must be greater than 0
        
        @rtype:  * 
        @return: None, SpliceDonor or ProjectedSpliceDonor object
        """
        # correct for bogus applied min_pssm_ratio and max_aa_dist
        min_pssm_ratio = max([min_pssm_ratio,1.0])
        max_aa_dist    = max([max_aa_dist,1])
        
        if not self._optimal_aligned_donor or\
        self._optimal_aligned_donor.__class__.__name__ ==\
        'DonorSiteCollectionGraph':
            return None
        if organism not in self._optimal_aligned_donor.organism_set():
            return None

        # get current site from the SpliceSiteCollectionGraph
        current_site = self._optimal_aligned_donor.get_organism_objects(
            organism)[0]
        sites = []
        for site in self._splicedonorgraph.get_organism_objects(organism):
            if site.pos == current_site.pos: continue
            if site.phase != current_site.phase: continue
            if site.pos < (current_site.pos-(max_aa_dist*3)): continue
            if site.pos > (current_site.pos+(max_aa_dist*3)): continue
            if site.pssm_score >= (current_site.pssm_score * min_pssm_ratio):
                # this is an elegiable optimized splice site
                sites.append( site )

        # if no sites, no higher PSSM scoring site available
        if not sites: return None
        # order list with SpliceSites by pssm_score attribute
        sites = order_list_by_attribute(
            sites,order_by='pssm_score',reversed=True)
        # return first element, which is the highest PSSM scoring SpliceSite
        return sites[0]

    # end of function propose_improved_donor


    def propose_improved_acceptor(self,organism,
        max_aa_dist=PROPOSE_IMPROVED_ACCEPTOR_MAX_AA_DIST,
        min_pssm_ratio=PROPOSE_IMPROVED_ACCEPTOR_MIN_PSSM_RATIO):
        """
        Return - if there is - a higher scoring acceptor of the same phase AND
        in the vicinity of the currently optimal aligned donor in AASG

        @type  organism: * (string)
        @param organism: Organism Identifier

        @type  max_aa_dist: integer
        @param max_aa_dist: maximal distance in AA positions from current site

        @type  min_pssm_ratio: float
        @param min_pssm_ratio: relative gain in PSSM score of improved site
        
        @attention: min_pssm_ratio must be greater than 1.0 !
        @attention: max_aa_dist must be greater than 0

        @rtype:  *
        @return: None, SpliceAcceptor or ProjectedSpliceAcceptor object
        """
        # correct for bogus applied min_pssm_ratio and max_aa_dist
        min_pssm_ratio = max([min_pssm_ratio,1.0])
        max_aa_dist    = max([max_aa_dist,1])

        if not self._optimal_aligned_acceptor or\
        self._optimal_aligned_acceptor.__class__.__name__ ==\
        'AcceptorSiteCollectionGraph':
            return None
        if organism not in self._optimal_aligned_acceptor.organism_set():
            return None

        # get current site from the SpliceSiteCollectionGraph
        current_site = self._optimal_aligned_acceptor.get_organism_objects(
            organism)[0]
        sites = []
        for site in self._spliceacceptorgraph.get_organism_objects(organism):
            if site.pos == current_site.pos: continue
            if site.phase != current_site.phase: continue
            if site.pos < (current_site.pos-(max_aa_dist*3)): continue
            if site.pos > (current_site.pos+(max_aa_dist*3)): continue
            if site.pssm_score >= (current_site.pssm_score * min_pssm_ratio):
                # this is an elegiable optimized splice site
                sites.append( site )

        # if no sites, no higher PSSM scoring site available
        if not sites: return None
        # order list with SpliceSites by pssm_score attribute
        sites = order_list_by_attribute(
            sites,order_by='pssm_score',reversed=True)
        # return first element, which is the highest PSSM scoring SpliceSite
        return sites[0]

    # end of function propose_improved_acceptor


    def optimize_aligned_donorgraph(self,
        max_aa_dist=PROPOSE_IMPROVED_DONOR_MAX_AA_DIST,
        min_pssm_ratio=PROPOSE_IMPROVED_DONOR_MIN_PSSM_RATIO,
        verbose=False):
        """
        (Try to) optimize the ADSG by replacing nearby higher PSSM scoring sites

        @type  max_aa_dist: integer
        @param max_aa_dist: maximal distance in AA positions from current site

        @type  min_pssm_ratio: float
        @param min_pssm_ratio: relative gain in PSSM score of improved site

        @type  verbose: Boolean
        @param verbose: print logging/debugging messages to STDOUT when True

        @rtype  donor_is_optimized: Boolean
        @return donor_is_optimized: True when at least 1 DonorSite improved
        """
        # return Boolean variable
        donor_is_optimized = False

        if self.donorCBG._short_name == "CBG" and\
        self.acceptorCBG._short_name == "CBG" and\
        self._spliceacceptorgraph and\
        self._splicedonorgraph and\
        self._optimal_aligned_acceptor and\
        self._optimal_aligned_donor and\
        self.is_compatible() and\
        self.donor_phase() in [0,1,2] and\
        self.acceptor_phase() in [0,1,2]:
            ####################################################################
            if verbose: print "POTENTIALLY OPTIMIZABLE cbgIF\n", self
            ####################################################################

            # optimize the sites in the AlignedDonorSiteGraph
            for organism in self.organism_set():
                try:
                    # search for an improved splice site node
                    nodeobj = self.propose_improved_donor(organism,
                                max_aa_dist=max_aa_dist,
                                min_pssm_ratio=min_pssm_ratio
                                )
                    if not nodeobj: continue

                    # get data of current splice site node
                    currentnode = self._optimal_aligned_donor.get_organism_nodes(
                        organism)[0]
                    currentobj =\
                        self._optimal_aligned_donor._node_object[currentnode]

                    ########################################################
                    if verbose:
                        print "ALTERNATIVE BETTER SPLICE SITE NODE!"
                        print organism
                        print "new:", nodeobj
                        print "cur:", currentobj
                    ########################################################

                    # replace and reset new current node and object
                    self._optimal_aligned_donor.replace_splice_site_object(
                            currentnode, currentobj, nodeobj )
                    currentnode = ( currentnode[0], currentnode[1], nodeobj.pos )
                    currentobj  = nodeobj
                    donor_is_optimized = True

                    ############################################################
                    if verbose: print "eof", organism
                    ############################################################
                except OrganismNotPresentInGraph:
                    pass
            ####################################################################
            if verbose: print "EOF", "donor"
            ####################################################################

        # return the counter if anything was optimized
        return donor_is_optimized

    # end of function optimize_aligned_donorgraph


    def optimize_aligned_acceptorgraph(self,
        max_aa_dist=PROPOSE_IMPROVED_ACCEPTOR_MAX_AA_DIST,
        min_pssm_ratio=PROPOSE_IMPROVED_ACCEPTOR_MIN_PSSM_RATIO,
        verbose=False):
        """
        (Try to) optimize the AASG by replacing nearby higher PSSM scoring sites
        
        @type  max_aa_dist: integer
        @param max_aa_dist: maximal distance in AA positions from current site
        
        @type  min_pssm_ratio: float
        @param min_pssm_ratio: relative gain in PSSM score of improved site
        
        @type  verbose: Boolean
        @param verbose: print logging/debugging messages to STDOUT when True
        
        @rtype  acceptor_is_optimized: Boolean
        @return acceptor_is_optimized: True when at least 1 AcceptorSite improved
        """
        # return Boolean variable
        acceptor_is_optimized = False
        
        if self.donorCBG._short_name == "CBG" and\
        self.acceptorCBG._short_name == "CBG" and\
        self._spliceacceptorgraph and\
        self._splicedonorgraph and\
        self._optimal_aligned_acceptor and\
        self._optimal_aligned_donor and\
        self.is_compatible() and\
        self.donor_phase() in [0,1,2] and\
        self.acceptor_phase() in [0,1,2]:
            ####################################################################
            if verbose: print "POTENTIALLY OPTIMIZABLE cbgIF\n", self
            ####################################################################
            
            # optimize the sites in the AlignedDonorSiteGraph
            for organism in self.organism_set():
                try:
                    # search for an improved splice site node
                    nodeobj = self.propose_improved_acceptor(organism,
                                max_aa_dist=max_aa_dist,
                                min_pssm_ratio=min_pssm_ratio
                                )
                    if not nodeobj: continue
                    
                    # get data of current splice site node
                    currentnode = self._optimal_aligned_acceptor.get_organism_nodes(
                        organism)[0]
                    currentobj =\
                        self._optimal_aligned_acceptor._node_object[currentnode]
                    
                    ########################################################
                    if verbose:
                        print "ALTERNATIVE BETTER SPLICE SITE NODE!"
                        print organism
                        print "new:", nodeobj
                        print "cur:", currentobj
                    ########################################################
                    
                    # replace and reset new current node and object
                    self._optimal_aligned_acceptor.replace_splice_site_object(
                                currentnode, currentobj, nodeobj )
                    currentnode = ( currentnode[0], currentnode[1], nodeobj.pos)
                    currentobj  = nodeobj
                    acceptor_is_optimized = True
                    
                    ############################################################
                    if verbose: print "eof", organism
                    ############################################################
                except OrganismNotPresentInGraph:
                    pass
            ####################################################################
            if verbose: print "EOF", "acceptor"
            ####################################################################
            
        # return the counter if anything was optimized
        return acceptor_is_optimized
    
    # end of function optimize_aligned_acceptorgraph
    # end of function optimize_aligned_acceptorgraph


    def optimize_aligned_splicesitegraph(self):
        """
        Shortcut function for self.optimize_aligned_[donor|acceptor]graph()
        
        @rtype:  tuple
        @return: ( Boolean, Boolean ) succes status for donor & acceptor
        """
        donor_is_optimized = self.optimize_aligned_donorgraph()
        acceptor_is_optimized = self.optimize_aligned_acceptorgraph()
        
        return (donor_is_optimized,acceptor_is_optimized)

    # end of function optimize_aligned_splicesitegraph

# end of class CodingBlockGraphInterface

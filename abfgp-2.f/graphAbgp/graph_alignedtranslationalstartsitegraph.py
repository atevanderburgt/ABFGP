################################################################################
### AlignedTranslationalStartSiteGraph class                                ####
################################################################################

# graphAbgp imports
from basal_alignedpssmobjects import BasalAlignedPssmObjectGraph
from subclass_tcodedata import GraphTcodeDataAccesFunctions
from subclass_exteriorcbgoptimality import ExteriorCbgOptimalityAnalyses
from exceptions import *

# Gene Imports
from gene.start import scan_pssm_tss


# Global variable Imports
from settings.alignedtranslationalstartsitegraph import *
from settings.translationalstartsites import *

# Python imports
from sets import Set

class AlignedTranslationalStartSiteGraph(BasalAlignedPssmObjectGraph,GraphTcodeDataAccesFunctions,ExteriorCbgOptimalityAnalyses):
    """
    AlignedTranslationalStartSiteGraph (ATSSG) class, inheriting from AlignedPssmObjectGraph
    """
    def __init__(self,max_node_count=None,min_pssm_score=None,aligned_site_aa_offset=None,tcode_5p_windowsize=TCODE_TSS_5P_WINDOW,tcode_3p_windowsize=TCODE_TSS_3P_WINDOW):
        """
        Initialize asa BasalAlignedPssmObjectGraph, then extend
        """
        BasalAlignedPssmObjectGraph.__init__(self,max_node_count=max_node_count,
                min_pssm_score=min_pssm_score,
                aligned_site_aa_offset=aligned_site_aa_offset ) 
        # attributes for TCODE data
        self._tcode5pscore = {}
        self._tcode3pscore = {}
        self._TCODE_5P_WINDOWSIZE = tcode_5p_windowsize
        self._TCODE_3P_WINDOWSIZE = tcode_3p_windowsize

        # attribute for storing this CBG object
        # when not HARD-SET into this postion,
        # is_optimal_gtgweakestnode() will crash! 
        self._codingblockgraph = None

        # is_optimal_xxx thresholds
        self._optimal_min_tcode      = ALIGNEDTSSGRAPH_OPTIMALITY_MIN_TCODE
        self._optimal_max_tcode      = ALIGNEDTSSGRAPH_OPTIMALITY_MAX_TCODE
        self._optimal_min_weight     = ALIGNEDTSSGRAPH_OPTIMALITY_MIN_WEIGHT
        self._optimal_max_weight     = ALIGNEDTSSGRAPH_OPTIMALITY_MAX_WEIGHT 
        self._optimal_min_gtgweakest = ALIGNEDTSSGRAPH_OPTIMALITY_MIN_GTGWEAKEST
        self._optimal_max_gtgweakest = ALIGNEDTSSGRAPH_OPTIMALITY_MAX_GTGWEAKEST

    # end of function __init__

    ########################################################################
    ### Build-in functions                                               ###
    ########################################################################

    def __str__(self):
        """ Nicely formatted oneliner reprenting this object """
        nodes = self.get_ordered_nodes()
        nodes = [ '%s(%s):%s' % (node[0],node[1],node[3]) for node in nodes ]
        nodes =  " ".join([ n for n in nodes])
        # and return a nicely formatted string
        return "<%s N%s score: %2.1f (tw:%2.1f betw:%2.1f pssmtw: %2.1f tcode:%1.2f-(%1.2f)-%1.2f) [%s]>" % (
                self.__class__.__name__,
                self.node_count(),
                self.cumulative_score(),
                self.total_weight(),
                self.total_binary_entropy(),
                self.total_pssm(),
                self.average_5p_tcode_score(),
                self.average_tcode_entropy(),
                self.average_3p_tcode_score(),
                nodes,
                )

    # end of function __str__


    ########################################################################
    ### Function to asses the optimality of the AlignedTSSGraph    ###
    ########################################################################

    def is_optimal(self,organism=None,node=None):
        """
        Is the AlignedTSSGraph optimal (for the given organism or node) ?

        @type  organism: *
        @param organism: One Organism Identifier (or None)

        @type  node: *
        @param node: One Node Identifier (or None)

        @rtype:  NoneBoolean
        @return: True, None or False
        """
        # get organism from node
        if node: organism = self.organism_by_node(node)

        tcode      = self.nonebooleanintegermapper( self.is_optimal_tcode_ratio() )
        weight     = self.nonebooleanintegermapper( self.is_optimal_weight() )
        gtgweakest = self.nonebooleanintegermapper( self.is_optimal_gtgweakestnode() )
        pssmwindow = self.nonebooleanintegermapper( self.is_optimal_pssm_in_window() )
        scores     = [tcode,weight,gtgweakest,pssmwindow]
        totalsum   = sum(scores)

        ###print "NOORG",[tcode,weight,gtgweakest,pssmwindow], "AVtcode:" , self.average_tcode_entropy(), "TOTALSUM:", totalsum

        # now translate into NoneBoolean outcome
        if totalsum > 1:                                            return True
        elif totalsum >= 0 and min(scores) == 0:                    return True
        elif not organism:                                          return False
        else:                                                       pass


        # switch to individual organism check. If this organism not in graph -> return False
        if organism not in self.organism_set(): return False

        # if this point in the function is reached, do a check for
        # the organism that is asked for
        tcode      = self.nonebooleanintegermapper( self.is_optimal_tcode_ratio(organism=organism) )
        weight     = self.nonebooleanintegermapper( self.is_optimal_weight(organism=organism) )
        gtgweakest = self.nonebooleanintegermapper( self.is_optimal_gtgweakestnode(organism=organism) )
        pssmwindow = self.nonebooleanintegermapper( self.is_optimal_pssm_in_window(organism=organism) )
        scores     = [tcode,weight,gtgweakest,pssmwindow]
        totalsum   = sum(scores)


      
        theNode  = self.node_by_organism(organism)
        ratioorg = self.average_tcode_entropy(
                    tcode5p=self._tcode5pscore[theNode],
                    tcode3p=self._tcode3pscore[theNode] )
        print organism,[tcode,weight,gtgweakest,pssmwindow], theNode, "AVtcode:" , self.average_tcode_entropy(), "TOTALSUM:", totalsum, ratioorg, self._tcode5pscore[theNode], self._tcode3pscore[theNode]

        # now translate into NoneBoolean outcome
        if totalsum >= 2:                                           return True
        elif totalsum >= 1:                                         return None
        else:                                                       return False

    # end of function is_optimal


    def is_optimal_pssm_in_window(self,organism=None,node=None):
        """
        Is the AlignedTranslationalStartSiteGraph optimal for its pssm score in a given window) ?

        @type  organism: *
        @param organism: One Organism Identifier (or None)

        @type  node: *
        @param node: One Node Identifier (or None)

        @rtype:  NoneBoolean
        @return: True, None or False
        """
        # get organism from node
        if node: organism = self.organism_by_node(node)

        if organism and organism not in self.organism_set():
            return False

        status = {} 
        for org in self.organism_set():
            if organism and org!=organism: continue
            # get orf and omsr
            theorf    = self._codingblockgraph.get_orfs_of_graph(organism=org)[0]
            omsr      = self._codingblockgraph.overall_minimal_spanning_range(organism=org)
            theNode   = self.node_by_organism(org)
            minntomsr = min(omsr)*3
            maxntomsr = max(omsr)*3
            thetss    = self.get_organism_objects(org)[0]
            ratio     = average_tcode_entropy_startcodon(tcode5p = self._tcode5pscore[theNode], tcode3p = self._tcode3pscore[theNode] )
            thetss._tcode_ratio = ratio
            thetss._pssm_x_tcode = thetss.pssm_score * thetss._tcode_ratio
            # define from where on to scan the genomesequence for TSSs
            startscanpos = min([thetss.pos,minntomsr])  # define start coord for TSS scanning
            startscanpos-=TSS_IS_OPTIMAL_5P_WINDOW      # decrease this pos with window size
            startscanpos = max([startscanpos,0])        # correct when <0 coord is created
            endscanpos   = min([thetss.pos+TSS_IS_OPTIMAL_3P_WINDOW,maxntomsr+TSS_IS_OPTIMAL_3P_WINDOW])

            # scan the input sequence first ~1600nt for high scoring TSS
            tsslist = scan_pssm_tss(theorf.inputgenomicsequence[startscanpos:endscanpos],min_pssm_score=TSS_IS_OPTIMAL_MIN_PSSM_SCORE)

            # now do some precalcs on the TSSs
            for tss in tsslist:
                # recalculate to absolute positions
                tss.pos+=startscanpos
                tss.start+=startscanpos
                tss.end+=startscanpos

                # get aaPos of TSS and obtain tcode entropy ratio
                aaPos = tss.pos / 3
                ( tcode5p,tcode3p ) = theorf.tcode_entropy_of_pos(
                        aaPos,
                        window_left=TSS_IS_OPTIMAL_5P_WINDOW,
                        window_right=TSS_IS_OPTIMAL_3P_WINDOW,
                        )
                ratio =  average_tcode_entropy_startcodon(tcode5p = tcode5p, tcode3p = tcode3p )
                tss._tcode_ratio = ratio
                # multiply pssm_score with tcode ratio
                tss._pssm_x_tcode = tss.pssm_score * tss._tcode_ratio
                #print "dsTSS", org, tss, tss._pssm_x_tcode, tcode5p, tcode3p, ratio

            # check if selectedtss score is optimal in the given window
            if not tsslist:
                # hmm... no TSS in tsslist, but thetss exists !?
                # Can be caused by a difference in tresholds used for TSS scanning
                # that means that thetss has a low value!
                status[org] = None # or True ?
            elif thetss.pssm_score >= max([ tss.pssm_score for tss in tsslist ]):
                status[org] = True
            elif thetss._pssm_x_tcode >= max([ tss._pssm_x_tcode for tss in tsslist ]):
                status[org] = True
            else:
                # get the highest pssm scoring TSS
                maxtss = None
                for tss in tsslist:
                    if tss.pssm_score == max([ tss.pssm_score for tss in tsslist ]):
                        maxtss = tss
                        break
                if maxtss.pos > thetss.pos:
                    # the highest PSSM scoring TSS is downstream of the selected tss
                    # that means: thetss is earlier encountered by the ribosome, so
                    # probably translation starts at thetss
                    status[org] = None
                else:    
                    ### print "ABOUT OT GET FALSE pssm_in_window", org,
                    ### print (thetss.pos, thetss.pssm_score),
                    ### print (maxtss.pos, maxtss.pssm_score),
                    ### print thetss._pssm_x_tcode,
                    ### print  max([ tss._pssm_x_tcode for tss in tsslist ])
                    # assign status==False for this organism for the TSS
                    status[org] = False

        # now check the status of the startcodons/TSSs
        if list(Set(status.values())) == [True]:   return True
        elif None in status.values() and not False in status.values()\
        and status.values().count(None) <= 2:      return None
        elif organism and status[organism]:        return True
        elif organism:                             return False
        else:
            if status.values().count(True) >= status.values().count(False):
                return None
            else:
                return False

    # end of function is_optimal_pssm_in_window


    def cumulative_score(self):
        """
        Cumulative score of a group of aligned TranslationalStartSites (TSS)

        @attention: See AlignedPssmObjectGraph.cumulative_score, CAN BE OVERRIDDEN HERE IF DESIRED!!
        """
        return BasalAlignedPssmObjectGraph.cumulative_score(self)

    # end of function cumulative_score


    def average_tcode_entropy(self,tcode5p=None,tcode3p=None):
        """
        Calculate ratio between Coding (rigth) and Non-Coding (left) side of Start-Codon
        """
        if tcode5p and tcode3p:
            return average_tcode_entropy_startcodon(tcode5p=tcode5p,tcode3p=tcode3p)
        elif self.get_nodes():
            return average_tcode_entropy_startcodon(
                    tcode5p = sum(self._tcode5pscore.values()) / float(self.node_count()),
                    tcode3p = sum(self._tcode3pscore.values()) / float(self.node_count()),
                    )
        else:
            return average_tcode_entropy_startcodon()

    # end of function average_tcode_entropy


    def tcode_offset_analyses(self,cbg,aa_step_size=33,range_of_steps=range(-5,3),window_left=201,window_right=201,default_window_size=201):
        """
        Analysis left and rigth of aligned start codons for a more likely TSS area

        @type  aa_step_size: Integer (positive)
        @param aa_step_size: Size of step (in AA) to next position to take into account

        @type  range_of_steps: list
        @param range_of_steps: range of steps to calculate; default range(-5,3); 0 == window around `pos`

        @type  window_left: Integer (positive)
        @param window_left: nt length to take as a left/5p window

        @type  window_rigth: Integer (positive)
        @param window_rigth: nt length to take as a rigth/3p window

        @type  default_window_size: Integer (positive)
        @param default_window_size: nt length of default window size when left/rigth window is omitted

        @rtype:  tuple of (positive) floats
        @return: tuple of ratios of averaged tcode window scores around a requested position
        """
        node2data = {}
        for tssNode in self.get_nodes():
            (org,orfid,posAA,posDNA) = tssNode
            theorf = cbg.get_orfs_of_graph(organism=org)[0]
            cbgNode = (org,orfid)
            scores = theorf.tcode_entropy_around_pos(posAA,
                    aa_step_size=aa_step_size,
                    range_of_steps=range_of_steps,
                    window_left=window_left,
                    window_right=window_right,
                    default_window_size=default_window_size
                    )
            node2data[cbgNode] = scores
        # calculate the ratios for the averaged sites
        ratios = {}
        for i in range(0,len(range_of_steps)):
            tcode5p = sum([ vlist[i][0] for vlist in node2data.values() ]) / float(len(node2data))
            tcode3p = sum([ vlist[i][1] for vlist in node2data.values() ]) / float(len(node2data))
            offset = range_of_steps[i]*aa_step_size
            ratios[offset] = average_tcode_entropy_startcodon(
                    tcode5p = tcode5p,
                    tcode3p = tcode3p,
                    )
        # done; return the ratios
        return ratios
    # end of function tcode_analyses_offset


    def togff(self,gff={},organism=None):
        """
        Create gff tuple for aligned splice site of a specific organism

        @attention: See AlignedPssmObjectGraph.togff, CAN BE OVERRIDDEN HERE IF DESIRED!!
        """
        return BasalAlignedPssmObjectGraph.togff(self,organism=organism,gff=gff)

    # end of function togff

# end of class AlignedTranslationalStartSiteGraph


########################################################################
### Helper Functions for AlignedTranslationalStartSiteGraph          ###
########################################################################

def average_tcode_entropy_startcodon(tcode5p=0.0,tcode3p=0.0):
    """
    Ratio between Coding (rigth) and Non-Coding (left) side of Start-Codon

    @type  tcode5p: positive float (0.450 - 1.320)
    @param tcode5p: average tcode score 5'/left/upstream of specific site

    @type  tcode3p: positive float (0.450 - 1.320)
    @param tcode3p: average tcode score 3'/right/downstream of specific site

    @rtype:  float (-0.450 - 2.160)
    @return: score attached to the likelihood that site is a startcodon; >1.0 is likely

    """
    TCODE_MAX_NONCODING = 0.740
    TCODE_MIN_CODING    = 0.950
    return 1.0 + ( TCODE_MAX_NONCODING - tcode5p ) * 2.0 + ( tcode3p -  TCODE_MAX_NONCODING )

# end of function average_tcode_entropy_startcodon



"""
Class with functions concerning LowSimilarityRegionCodingBlockGraphs (lsrCBGs)
in the GeneStructureOfCodingBlocks (GSG) used in Alignment Based Gene Prediction

    cexpander_inframe_intron_search

    cexpander_lsrcbg_intron_search

    search_for_lowsimilarity_regions
        DEPRECATED by cexpander_inframe_intron_search

    search_for_inframe_introns
        DEPRECATED by search_for_lowsimilarity_regions

    check_lsrcbgs_for_inframe_introns
        DEPRECATED by cexpander_lsrcbg_intron_search


"""

# Module metadata
__authors__ = "Ate van der Burgt"
__license__ = "MIT"

# grapgAbgp Imports
from codingblock_operations import deepcopy_with_removed_nodes
from codingblock_optimization import _update_cbg_with_pacbporf_replacements
from codingblock_splitting import create_intermediate_lowsimilarity_region
from exceptions import *

# ABGP Imports
from codingblockgraphinterface import CodingBlockGraphInterface
from lib_stopwatch import StopWatch
from gene.intron import _filter_intron_list
import lib_cexpander
import pacb

# Python Imports
from sets import Set
from copy import deepcopy

# Global variables
from settings.genestructure import MIN_INTRON_NT_LENGTH
MIN_TOTAL_PSSM_INFRAME_INTRON = 3.0
from settings.inframeintron import ALIGNED_INTRON_MIN_AA_LENGTH



class LowSimilarityRegionCodingBlockGraphFunctions:
    """ lsrCBG related functions in a GSG """

    def create_intermediary_lsrcbgs(self,verbose=False):
        """
        Create lsrCBGs in between CBGs with identical sets of nodes

        @attention: function can be called as often as desired

        @rtype:  Boolean
        @return: True or False weather or not lsrCBGs are created
        """
        # return status boolean weather or not a lsrCBG is added
        RETURN_STATUS_LSRCBG_IS_ADDED = False

        # loop BACKWARDS over the CBGs in case of an insert
        for i in range(len(self)-1,0,-1):
            # get combinations of 2 neighbouring CBGs
            (firstCBG,secondCBG) = self.codingblockgraphs[i-1:i+1]
            # ignore if one of them IS_IGNORED, lsrCBG or SPLITTED
            if firstCBG.IS_IGNORED:  continue
            if secondCBG.IS_IGNORED: continue
            if firstCBG._short_name  == "lsrCBG": continue
            if secondCBG._short_name == "lsrCBG": continue
            if firstCBG.IS_3P_SPLITTED:  continue
            if secondCBG.IS_5P_SPLITTED: continue
            # ignore if not all mutual nodes
            if firstCBG.node_set().symmetric_difference(secondCBG.get_nodes()):
                # reset possible bogus IS_SPLITTED variable settings
                # in CBGS. Can be instantiated by deletion of CBGs
                firstCBG.IS_3P_SPLITTED = False
                secondCBG.IS_5P_SPLITTED = False
                if not firstCBG.IS_5P_SPLITTED:
                    firstCBG.IS_SPLITTED = False
                if not secondCBG.IS_3P_SPLITTED:
                    secondCBG.IS_SPLITTED = False
                continue

            # If this point is reached, firstCBG and secondCBG are CBGs with
            # exactly the same nodes
            # create intermediate lsrCBG
            lsrCBG = create_intermediate_lowsimilarity_region(firstCBG,secondCBG)

            # check if this lsrCBG got any nodes (lsromsr not added!)
            if lsrCBG.get_nodes():
                ################################################################
                if verbose:
                    print lsrCBG
                    print "potential inframe intron:",
                    print lsrCBG.potentially_contains_inframe_intron()
                ################################################################
                # update the status of CBG firstCBG and secondCBG
                firstCBG.IS_SPLITTED     = True
                firstCBG.IS_3P_SPLITTED  = True
                secondCBG.IS_SPLITTED    = True
                secondCBG.IS_5P_SPLITTED = True
                # insert the LowSimilarityRegionCodingBlockGraph
                # at the proper position
                self.codingblockgraphs.insert(i,lsrCBG)
                RETURN_STATUS_LSRCBG_IS_ADDED = True

        # return the status weather or not a lsrCBG is added
        return RETURN_STATUS_LSRCBG_IS_ADDED

    # end of function create_intermediary_lsrcbgs


    def remove_incorrect_lsrcbgs(self,verbose=False):
        """
        Remove incorrectly placed lsrCBGs and fix IS_SPLITTED status of all CBGs

        @attention: function can be called as often as desired
        @attention: function must be called after CBG removal functions

        @rtype:  Boolean
        @return: True or False weather or not lsrCBGs are removed
        """
        # return status boolean weather or not a lsrCBG is removed
        RETURN_STATUS_LSRCBG_IS_REMOVED = False

        for pos in self.cbgpositerator(reversed=True):
            lsrCBG = self.codingblockgraphs[pos]
            if lsrCBG._short_name != "lsrCBG": continue
            # get next and prev CBG for node comparison
            prevCBG, nextCBG = None, None
            if pos > 0:
                prevCBG = self.codingblockgraphs[pos-1]
            if pos < len(self)-1:
                nextCBG = self.codingblockgraphs[pos+1]

            # check if prevCBG - lsrCBG - nextCBG have mutual nodes
            if (nextCBG and not prevCBG) or (prevCBG and not nextCBG) or\
            (prevCBG and lsrCBG.node_set().difference(prevCBG.get_nodes())) or\
            (nextCBG and lsrCBG.node_set().difference(nextCBG.get_nodes())) or\
            (prevCBG and nextCBG and\
            prevCBG.node_set().symmetric_difference(nextCBG.get_nodes())):
                # remove this lsrCBG and clean neighboring CBGinterfaces
                popped = self.codingblockgraphs.pop(pos)
                RETURN_STATUS_LSRCBG_IS_REMOVED = True
                if prevCBG:
                    prevCBG._CBGinterface3p = None
                    prevCBG.IS_3P_SLITTED   = False
                    if not prevCBG.IS_5P_SPLITTED:
                        prevCBG.IS_SPLITTED = False
                if nextCBG:
                    nextCBG._CBGinterface5p = None
                    nextCBG.IS_5P_SLITTED   = False
                    if not nextCBG.IS_3P_SPLITTED:
                        nextCBG.IS_SPLITTED = False

        # now check the IS_SPLITTED status of the CBGs
        # only remove WRONG IS_SPLITTED status, do not
        # create new IS_SPLITTED = True, becauuse this
        # interfears with the create_intermediary_lsrcbgs() function 

        for pos in self.cbgpositerator(reversed=True):
            cbg = self.codingblockgraphs[pos]
            if cbg._short_name != "CBG": continue
            # get next and prev CBG for node comparison
            prevCBG, nextCBG = None, None
            if pos > 0:
                prevCBG = self.codingblockgraphs[pos-1]
                if prevCBG.node_set().difference(cbg.node_set()):
                    # different nodes in prevCBG -> no split
                    cbg.IS_5P_SPLITTED     = False
                    prevCBG.IS_3P_SPLITTED = False
            else:
                # first CBG -> no 5' split
                cbg.IS_5P_SPLITTED = False
            if pos < len(self)-1:
                nextCBG = self.codingblockgraphs[pos+1]
                if nextCBG and nextCBG.node_set().difference(cbg.node_set()):
                    # different nodes in nextCBG -> no split
                    cbg.IS_3P_SPLITTED     = False
                    nextCBG.IS_5P_SPLITTED = False
            else:
                # final CBG -> no 3' split
                cbg.IS_3P_SPLITTED = False
            # update the main IS_SPLITTED argument
            cbg.IS_SPLITTED = True in [ cbg.IS_5P_SPLITTED, cbg.IS_3P_SPLITTED ]

        # return the status weather or not a lsrCBG is removed
        return RETURN_STATUS_LSRCBG_IS_REMOVED

    # end of function remove_incorrect_lsrcbgs


    def gsg_cexpander_enlarge_lsrcbgs(self,verbose=False):
        """
        """
        lsr_coords_changed = 0

        for pos in range(1,len(self)-1):
            if self.codingblockgraphs[pos].__class__.__name__ !=\
            'LowSimilarityRegionCodingBlockGraph':
                continue
            # get previous and next CBG
            prevCBG = self.codingblockgraphs[pos-1]
            nextCBG = self.codingblockgraphs[pos+1]

            # obtain current CBG data for logging when something fails
            strreprPrevCbg = str(prevCBG)
            strreprLsrCbg  = str(self.codingblockgraphs[pos])
            strreprNextCbg = str(nextCBG)

            # deepcoy Pacbps in case cexpander omsr border gaps
            # operations mingles the CBG(s)
            bckp_prevcbg_pacbps = deepcopy(prevCBG.pacbps)
            bckp_nextcbg_pacbps = deepcopy(nextCBG.pacbps)

            try:
                # optimize the CBGs around the lsrCBG with cexpander data
                statusP = lib_cexpander.cexpander_checkCBG4omsrbordergaps(
                        prevCBG, omit5pside = True )
                statusN = lib_cexpander.cexpander_checkCBG4omsrbordergaps(
                        nextCBG, omit3pside = True )
                if statusP or statusN:
                    # if one or both CBGs changed -> new lsrCBG
                    if statusP: prevCBG.create_cache()
                    if statusN: nextCBG.create_cache()
                    newLsrCBG = create_intermediate_lowsimilarity_region(
                            prevCBG,nextCBG)
                    prepare_lsrcbg_and_cbg_for_gsg_insertion(prevCBG,newLsrCBG)
                    prepare_lsrcbg_and_cbg_for_gsg_insertion(newLsrCBG,nextCBG)
                    self.codingblockgraphs[pos] = newLsrCBG
                    ############################################################
                    if verbose:
                        print "gsg_cexpander_enlarge_lsrcbgs WAS:"
                        print strreprPrevCbg
                        print strreprLsrCbg
                        print strreprNextCbg
                        print "gsg_cexpander_enlarge_lsrcbgs IS:"
                        print prevCBG
                        print newLsrCBG
                        print nextCBG
                    ############################################################
                    lsr_coords_changed += 1

            except NoOverallMinimalSpanningRange:
                # NoOverallMinimalSpanningRange Exception;
                # that is - normally - the signal for deleting this CBG.
                # However, here it is a SEVERE problem. The CBG is 'lost' due to
                # the cexpander optimization. This will result in a later crash
                ########################################################################
                if verbose:
                    print "SeriousWarning: CBG lost due to gsg_cexpander_enlarge_lsrcbgs"
                    print "NoOverallMinimalSpanningRange"
                    print strreprPrevCbg
                    print strreprLsrCbg
                    print strreprNextCbg
                ########################################################################
                # Restore CBGs and lsrCBG in state as before this operation
                prevCBG.pacbps = bckp_prevcbg_pacbps
                prevCBG.create_cache()
                nextCBG.pacbps = bckp_nextcbg_pacbps
                nextCBG.create_cache()
                restoredLsrCBG = create_intermediate_lowsimilarity_region(
                            prevCBG,nextCBG)
                prepare_lsrcbg_and_cbg_for_gsg_insertion(prevCBG,restoredLsrCBG)
                prepare_lsrcbg_and_cbg_for_gsg_insertion(restoredLsrCBG,nextCBG)
                self.codingblockgraphs[pos] = restoredLsrCBG

            except lib_cexpander.ZeroUniformlyAlignedPositions:
                # due to optimization, the multiple alignment collapsed
                # that is - normally - the signal for deleting this CBG.
                # However, here it is a SEVERE problem. The CBG is 'lost' due to
                # the cexpander optimization. This will result in a later crash
                ########################################################################
                if verbose:
                    print "SeriousWarning: CBG lost due to gsg_cexpander_enlarge_lsrcbgs"
                    print "lib_cexpander.ZeroUniformlyAlignedPositions"
                    print strreprPrevCbg
                    print strreprLsrCbg
                    print strreprNextCbg
                ########################################################################
                # Restore CBGs and lsrCBG in state as before this operation
                prevCBG.pacbps = bckp_prevcbg_pacbps
                prevCBG.create_cache()
                nextCBG.pacbps = bckp_nextcbg_pacbps
                nextCBG.create_cache()
                restoredLsrCBG = create_intermediate_lowsimilarity_region(
                            prevCBG,nextCBG)
                prepare_lsrcbg_and_cbg_for_gsg_insertion(prevCBG,restoredLsrCBG)
                prepare_lsrcbg_and_cbg_for_gsg_insertion(restoredLsrCBG,nextCBG)
                self.codingblockgraphs[pos] = restoredLsrCBG

            except:
                # unexpected exception -> raise!
                raise "UnExpectedException in checkCBGs4omsrbordergaps"

        # return the counter how much lsrCBGs are changed
        return lsr_coords_changed

    # end of function gsg_cexpander_enlarge_lsrcbgs


    def gsg_cexpander_inframe_intron_search(self,**kwargs):
        """
        Scan the CBGs in the genestructure for hidden inframe introns

        @attention: see cbg_cexpander_inframe_intron_search function for
                    argument documentation

        @rtype:  tuple ( int, int )
        @return: ( number of created lsrCBGs, number of created CBGs )
        """

        # return counters
        created_lsrcbg_cnt, created_cbg_cnt = 0, 0

        # Loop reversed through genestructure: after CBG splitting the
        # remainder of the GSG.codingblockgraphs list positions stay intact.
        for pos in range(len(self)-1,-1,-1):
            cbg = self.codingblockgraphs[pos]
            # skip IGNORED and lsrCBG codingblocks
            if cbg.__class__.__name__ == 'LowSimilarityRegionCodingBlockGraph':
                continue
            if cbg.IS_IGNORED:
                continue

            # do the cbg inframe intron cexpander search
            result = cbg_cexpander_inframe_intron_search(cbg,**kwargs)
            if not result:
                pass
            else:
                # update splitted CBGs to genestructure GSG
                self.codingblockgraphs.__setslice__(pos,pos+1,result)

                # update the result counters
                lsrcbg_cnt = [ _cbg.__class__.__name__  for _cbg in\
                    result ].count('LowSimilarityRegionCodingBlockGraph')
                created_lsrcbg_cnt += lsrcbg_cnt
                created_cbg_cnt += ( len(result) - lsrcbg_cnt - 1 )

        # return both counters as a result
        return ( created_lsrcbg_cnt, created_cbg_cnt )

    # end of function gsg_cexpander_inframe_intron_search


    def search_for_lowsimilarity_regions(self,aligned_intron_min_aa_length=ALIGNED_INTRON_MIN_AA_LENGTH,verbose=False):
        """
        Search CBGs in genestructure for lowsimilarity regions
        """

        ################################################################
        if verbose:
            stw = StopWatch(name='lsrCBGsearch')
            stw.start()
        ################################################################

        # Loop reversed through genestructure to make sure that once
        # a CBG is splitted, the positions of the remainder of the
        # list stay intact.
        for posinGSG in range(len(self)-1,-1,-1):
            sg = self.codingblockgraphs[posinGSG]
            # skip IGNORED, lsrCBG and CBGs that are incomplete (still await HMM completion) 
            if sg.IS_IGNORED: continue
            if sg.__class__.__name__ == 'LowSimilarityRegionCodingBlockGraph': continue
            if sg.node_count() < self.EXACT_SG_NODE_COUNT: continue

            if verbose: print stw.lap(), posinGSG, "start"

            # check for potential aligned intron
            if sg.potentially_contains_aligned_intron(window_aa_size=aligned_intron_min_aa_length):
                ########################################################
                if verbose:
                    print stw.lap(), posinGSG, "found"
                    for k,v in sg.getomsrproteinsequences().iteritems():
                        print ">%s\n%s\n" % (k,v)
                    print "ABOUT TO SPLIT:", sg
                    print sg._cexpander.binarystring,
                    print sg._cexpander.projected_on
                    sg.printmultiplealignment()
                    for k,pacbp in sg.pacbps.iteritems(): print k, pacbp
                ########################################################
                # now actually split by inframe intron
                res = sg.split_codingblock_by_inframe_intron()
                if len(res) == 1:
                    # no inframe intron found here
                    pass
                else:
                    # prepare the CBGs for insertion 
                    for pos in range(0,len(res)):
                        splittedCBG = res[pos]
                        splittedCBG.extend_pacbporfs(self.input)
                        splittedCBG.update_edge_weights_by_minimal_spanning_range()
                        splittedCBG.IS_SPLITTED = True
                        if pos > 0:
                            splittedCBG.IS_5P_SPLITTED = True
                            splittedCBG.IS_FIRST = False
                        if pos < len(res)-1:
                            splittedCBG.IS_3P_SPLITTED = True
                            splittedCBG.IS_LAST = False
                        # (re)create the cache for the splitted CBGs
                        splittedCBG.create_cache()
                        ################################################
                        if verbose:
                            print stw.lap(), posinGSG, "done!"
                            print "SUCCESFULLY SPLITTED:", splittedCBG
                            splittedCBG.printmultiplealignment()
                            print splittedCBG._cexpander.binarystring, 
                            print splittedCBG._cexpander.projected_on
                            print splittedCBG._omsr
                            for trf in splittedCBG._cexpander._transferblocks:
                                print trf.binarystring, trf.projected_on
                            for k,v in splittedCBG._cexpander.inputsequences.iteritems():
                                print v,"\t",k
                            for _org,orflist in splittedCBG.get_orfs_of_graph().iteritems():
                                print orflist[0], _org
                            for pacbp in splittedCBG.pacbps.values():
                                print pacbp
                                pacbp.print_protein(_linesize=100)
                        ################################################

                    # create lsrCBGs and cbgIFs between them by looping in reversed
                    # order over all pairs of CBGs (because lsrCBG insertion in list)
                    for pos in range(len(res)-2,-1,-1):
                        cbgL,cbgR = res[pos:pos+2]
                        lsrCBG = create_intermediate_lowsimilarity_region(cbgL,cbgR)
                        res.insert(pos+1,lsrCBG)
                        # create cbgIF between the CBGs and the lsrCBG
                        # just create -> cbgIF with lsrCBG is immediately is_optimal()
                        cbgIFa = CodingBlockGraphInterface(cbgL,lsrCBG)
                        cbgIFb = CodingBlockGraphInterface(lsrCBG,cbgR)
                        # set cbgIF objects to the CBGs and the lsrCBG
                        cbgL._CBGinterface3p   = cbgIFa
                        lsrCBG._CBGinterface5p = cbgIFa
                        lsrCBG._CBGinterface3p = cbgIFb
                        cbgR._CBGinterface5p   = cbgIFb

                    # update the first and last CBG in this list with the
                    # cbgIFs of the parental CBG (variable sg)
                    res[0]._CBGinterface5p =  sg._CBGinterface5p
                    res[-1]._CBGinterface3p = sg._CBGinterface3p
                    # update the original IS_FIRST/IS_LAST status
                    res[0].IS_FIRST = sg.IS_FIRST
                    res[-1].IS_LAST = sg.IS_LAST

                    # and set splittedCBGs to genestructure
                    # by replacing the existing CBG (variable sg) on the
                    # position posinGSG with the list op splitted CBGs
                    self.codingblockgraphs.__setslice__(posinGSG,posinGSG+1,res)

            else:
                # nope, no potential inframe intron; just append
                ###print sg.total_weight(), False
                pass
    
    # end of function search_for_lowsimilarity_regions


    def search_for_inframe_introns(self,aligned_intron_min_aa_length=ALIGNED_INTRON_MIN_AA_LENGTH):
        """
        Search CBGs in genestructure for remaining inframe introns
        """
        # Loop reversed through genestructure to make sure that once
        # a CBG is splitted, the positions of the remainder of the
        # list stay intact.
        for pos in range(len(self)-1,-1,-1):
            sg = self.codingblockgraphs[pos]
            if sg.IS_IGNORED: continue    # skip IGNORED codingblocks
            print "TRYING TO SPLIT:", sg
            if sg.potentially_contains_aligned_intron(window_aa_size=aligned_intron_min_aa_length):
                print "ABOUT TO SPLIT:", sg
                for k,pacbp in sg.pacbps.iteritems(): print k, pacbp
                sg.printmultiplealignment()
                res = sg.split_codingblock_by_inframe_intron()
                if len(res) == 1:
                    # no inframe intron found here
                    pass
                else:
                    for splittedCBG in res:
                        splittedCBG.extend_pacbporfs(self.input)
                        splittedCBG.update_edge_weights_by_minimal_spanning_range()
                        print "SUCCESFULLY SPLITTED:", splittedCBG
                        splittedCBG.printmultiplealignment()
                        print splittedCBG._omsr
                    # and set splittedCBGs to genestructure
                    self.codingblockgraphs.__setslice__(pos,pos+1,res)

            else:
                # nope, no potential inframe intron; just append
                pass
    
    # end of function search_for_inframe_intron



    def check_lsrcbgs_for_inframe_introns(self,verbose=False):
        """
        Check the lsrCBGs in the GSG and see if these regions can better be explained by an inframe intron
        """
        INFRAME_INTRONS_PREDICTED = 0
        LSR_RECREATED             = 0
        for cbgpos in range(len(self)-1,-1,-1):
            cbg = self.codingblockgraphs[cbgpos]
            if cbg.__class__.__name__ != 'LowSimilarityRegionCodingBlockGraph':
                continue
            # do the inframe intron analyses on a lsrCBG
            inframeintrons = cbg.potentially_contains_inframe_intron(verbose=verbose)
            # aparantly it seems possible to create one or more introns in the lsrCBG
            if inframeintrons:
                # get the bordering CBGs
                prev = self.codingblockgraphs[cbgpos-1]
                next = self.codingblockgraphs[cbgpos+1]
                # make CBGInterface between prev and next;
                # reset the _IS_SPLITTED tags!
                prev._splicedonorgraph = None
                prev._CBGinterface3p   = None
                prev._forced_3p_ends   = {}
                prev.IS_3P_SPLITTED    = False
                prev.IS_SPLITTED       = prev.IS_5P_SPLITTED
                next._spliceacceptorgraph = None
                next._CBGinterface5p   = None
                next._forced_5p_ends   = {}
                next.IS_5P_SPLITTED    = False
                next.IS_SPLITTED       = next.IS_3P_SPLITTED
        
                # create an actual CBGInterface of both CBGs around the lsrCBG
                cbgIF = CodingBlockGraphInterface(prev,next)
                if verbose: print cbgIF
                # re-harvest splice sites; store ALL the intron-projected sites
                cbgIF.harvest_splice_sites(allow_phase_shift=False,store_all_projected_sites=True)
                if verbose: print cbgIF
                # now remove all non-projected splice-sites in organisms that
                # are not reported to have a potential inframe intron
                cbgIF.allow_intron_in_organisms(inframeintrons)
                cbgIF.find_conserved_splice_sites()
                if verbose:
                    print cbgIF
                    print "compatible:", cbgIF.is_compatible(), "optimal:", cbgIF.is_optimal()
                    print cbgIF._optimal_aligned_donor
                    print cbgIF._optimal_aligned_acceptor
                # yes, this is what we expect; a compatible CBGInterface!
                # this very likely represents an inframe intron!
                if cbgIF.is_compatible():
                    # remove the lsrCBG from the GSG
                    lsrCBG = self.codingblockgraphs.pop(cbgpos)
                    # set the CBGInterface object in next and prev CBG
                    prev._CBGinterface3p = cbgIF
                    next._CBGinterface5p = cbgIF
                    # increase the counter of number of inframe introns predicted
                    INFRAME_INTRONS_PREDICTED+=1
                    ############################################################
                    if verbose: print "INFRAME INTRON PREDICTED!!"
                    ############################################################

                else:
                    # nope, this does not seem like a proper inframe intron
                    # reset the CBGs and the lsrCBG objects as they were!

                    # If this point is reached, `first` and `second` are CBGs with exactly the same nodes
                    # create intermediate lsrCBG
                    prev.IS_SPLITTED    = True
                    prev.IS_3P_SPLITTED = True
                    next.IS_SPLITTED    = True
                    next.IS_5P_SPLITTED = True
                    lsrCBG = create_intermediate_lowsimilarity_region(prev,next)
                    self.codingblockgraphs[cbgpos]   = lsrCBG

                    # recreate the CBGInterfaces (I)
                    cbgIFa = CodingBlockGraphInterface(prev,lsrCBG)
                    cbgIFa.harvest_splice_sites()
                    cbgIFa.find_conserved_splice_sites()
                    # set the interface object to the CBGs in GSG
                    prev._CBGinterface3p   = cbgIFa
                    lsrCBG._CBGinterface5p = cbgIFa

                    # recreate the CBGInterfaces (II)
                    cbgIFb = CodingBlockGraphInterface(lsrCBG,next)
                    cbgIFb.harvest_splice_sites()
                    cbgIFb.find_conserved_splice_sites()
                    # set the interface object to the CBGs in GSG
                    lsrCBG._CBGinterface3p = cbgIFb
                    next._CBGinterface5p   = cbgIFb

                    ############################################################
                    if verbose: print "NO COMPATIBLE SITE!"
                    ############################################################

                    ###for org in inframeintrons:
                    ###    print org, "NO COMPATIBLE SITES FOUND!"
                    ###    print prev
                    ###    print cbgIF
                    ###    print next
                    ###    theorf = next.get_orfs_of_graph(organism = org )[0]
                    ###    print theorf
                    ###    theorf.printproteinanddna()
                    ###    for donor in theorf._donor_sites: print donor
                    ###    for acceptor in theorf._acceptor_sites: print acceptor
        
        # return number of found inframe introns
        return INFRAME_INTRONS_PREDICTED

    # end of function check_lsrcbgs_for_inframe_introns

# end of class LowSimilarityRegionCodingBlockGraphFunctions


def _filter_putative_inframe_intron_list(introns,org,inframe_intron_criteria):
    """
    @type  introns: list
    @param introns: list with all inframe introns in a certain Orf

    @attention: only use from within cbg_cexpander_inframe_intron_search function

    @rtype  introns: list
    @return introns: largely reduced list with suitable inframe introns
    """
    if inframe_intron_criteria[org].has_key('max_nt_length'):
        introns = _filter_intron_list(introns, filter_by='length',
            criterion=inframe_intron_criteria[org]['max_nt_length'],
            operator='<=')
    if inframe_intron_criteria[org].has_key('min_nt_length'):
        introns = _filter_intron_list(introns, filter_by='length',
            criterion=inframe_intron_criteria[org]['min_nt_length'],
            operator='>=')
    if inframe_intron_criteria[org].has_key('min_donor_pos'):
        introns = _filter_intron_list(introns,filter_by='donor_pos',
            criterion=inframe_intron_criteria[org]['min_donor_pos'],
            operator='>')
    if inframe_intron_criteria[org].has_key('max_acceptor_pos'):
        introns = _filter_intron_list(introns,filter_by='acceptor_pos',
            criterion=inframe_intron_criteria[org]['max_acceptor_pos'],
            operator='>')
    if inframe_intron_criteria[org].has_key('min_total_pssm'):
        introns = _filter_intron_list(introns,filter_by='total_pssm',
            criterion=inframe_intron_criteria[org]['min_total_pssm'],
            operator='>=')
    return introns

# end of function _filter_putative_inframe_intron_list


def prepare_lsrcbg_and_cbg_for_gsg_insertion(cbgL,cbgR):
    """
    """
    cbgL.IS_SPLITTED    = True
    cbgL.IS_3P_SPLITTED = True
    cbgR.IS_SPLITTED    = True
    cbgR.IS_5P_SPLITTED = True
    # an lsrCBG / splitted interface is always is_optimal!
    cbgIF = CodingBlockGraphInterface(cbgL,cbgR)
    cbgL._CBGinterface3p = cbgIF
    cbgR._CBGinterface5p = cbgIF

# end of prepare_lsrcbg_and_cbg_for_gsg_insertion


def cbg_cexpander_inframe_intron_search(self,
        min_total_pssm_score = MIN_TOTAL_PSSM_INFRAME_INTRON,
        min_intron_nt_length = MIN_INTRON_NT_LENGTH,
        verbose=False):
        """
        @type  self: CodingBlockGraph
        @param self: CodingBlockGraph instance

        @type  min_total_pssm_score: float
        @param min_total_pssm_score: MIN_TOTAL_PSSM_INFRAME_INTRON

        @type  min_intron_nt_length: integer
        @param min_intron_nt_length: MIN_INTRON_NT_LENGTH

        @type  verbose: Boolean
        @param verbose: print status/debugging messages to STDOUT

        @rtype:  list or False
        @return: list with new (sub)CBGs or False when not splitted
        """
        ########################################################################
        if verbose:
            stw = StopWatch(name="cexpCbgIfIntron")
            stw.start()
        ########################################################################

        # return variable; list of splitted CBGs.
        return_cbg_list = [ self ]

        # create cexpander multiplealignment blocks
        cbgMA = lib_cexpander.cexpander2multiplealignment(self._cexpander,
                verbose=verbose)

        # In freak-accident cases (one in thousends of times), cexpander produces
        # unequal amount of 1's in the binarystrings. This is theoretically impossible.
        # Problem is worked on; in the meanwhile, cexpander2multiplealignment returns
        # False in these cases. Catch this here by quiting current 
        # cbg_cexpander_inframe_intron_search() function call and return False
        TODO=True
        if not cbgMA: return False

        ########################################################################
        if verbose:
            print stw.lap()
            blockscnt = len( cbgMA[ cbgMA.keys()[0] ] )
            print self
            print "BLOCKS:", blockscnt, self._cexpander.binarystring,
            print self._cexpander.projected_on
            for org in cbgMA.keys():
                print org, "\t", 
                for blockid in range(0,blockscnt):
                    if cbgMA[org][blockid].count("1") >= 1:
                        print len(cbgMA[org][blockid]), 
                    else:
                        print cbgMA[org][blockid], 
                print ""
        ########################################################################

        # loop over the aligned cexpander blocks and check the 
        # non-uniformly aligned blocks for length variation
        blockscnt  = len( cbgMA[ cbgMA.keys()[0] ] )
        oricbgomsr = self.overall_minimal_spanning_range()

        for blockid in range(0,blockscnt):
            # obtain non-uniformly aligned AA lengths for this block
            lengths = {}
            for org in cbgMA.keys():
                lengths[org] = cbgMA[org][blockid].count("0")
            # skip the uniformly aligned blocks
            if list(Set(lengths.values())) == [0]: continue
            ####################################################################
            if verbose: print stw.lap(), "lengths:", lengths
            ####################################################################

            # obtain coordinates for this area
            lsrcoords = {}
            for org in cbgMA.keys():
                node = self.node_by_organism(org)
                coordSta = min(oricbgomsr[node])
                # make summation of length of preceeding (non)aligned blocks
                for i in range(0,blockid):
                    coordSta += cbgMA[org][i].count("1") +\
                                cbgMA[org][i].count("0")
                # end coord is start coord + length of current block
                coordEnd = coordSta + lengths[org]
                lsrcoords[org] = ( coordSta, coordEnd )

            ####################################################################
            if verbose: print stw.lap(), "lsrcoords:", lsrcoords
            ####################################################################

            # translate AA lengths to NT lengths
            for k in lengths.keys(): lengths[k] = lengths[k]*3

            # check lenght discrepancy and assign putative inframe introns
            putative_inframe_intron_orgs =\
                _length_discrepancy_to_potential_inframe_introns(lengths)

            if not putative_inframe_intron_orgs:
                # no length discrepancy that can represent an inframe intron
                continue

            # organisms/genes for which an inframe intron can be an improvement
            # data dictionary. Keys: 'max_nt_length', 'min_nt_length', 
            # 'min_donor_pos', 'max_acceptor_pos', 'min_total_pssm'
            inframe_intron_criteria = {}

            # find putative inframe introns in assigned genes/organisms
            putative_inframe_introns = {}
            for org in putative_inframe_intron_orgs:
                # assign inframe intron criteria for this organism
                inframe_intron_criteria[org] = {
                    'min_nt_length'     : min_intron_nt_length,
                    'min_total_pssm'    : min_total_pssm_score,
                    'min_donor_pos'     : (min(lsrcoords[org]) - 5) * 3,
                    'max_acceptor_pos'  : (max(lsrcoords[org]) + 5) * 3,
                    }

                # search for potential introns that can be responsible for this event
                theorf = self.get_orfs_of_graph(organism=org)[0]
                introns = pacb.connecting.merge_orfs_with_intron( theorf,theorf,
                            min_intron_nt_length=min_intron_nt_length
                            )

                ################################################################
                if verbose: print "introns:", org, len(introns), "raw"
                ################################################################

                # filter introns for all outside the OMSR, to short, to long,
                # total pssm_score etc
                introns = _filter_putative_inframe_intron_list(
                        introns,org,inframe_intron_criteria)
                putative_inframe_introns[org] = introns
                ################################################################
                if verbose: print "introns:", org, len(introns), "filtered"
                ################################################################

            # check if all putative_inframe_intron_orgs have indeed introns
            # and check if all have at least a single intron phase in common
            if 0 in [ len(ill) for ill in putative_inframe_introns.values() ]:
                # no introns in one or more organisms/genes -> continue
                continue
            if len( putative_inframe_introns )> 1:
                # do phase check in all organisms/genes
                phases = Set([0,1,2])
                for org, intronlist in putative_inframe_introns.iteritems():
                    thisphases = Set([ intron.phase for intron in intronlist ])
                    phases.intersection_update(thisphases)
                if len(phases) == 0:
                    ################################################################
                    if verbose: print "no mutual phase -> no cbgIF.is_optimal()"
                    ################################################################
                    # no mutual phase -> no cbgIF.is_optimal() possible lateron
                    continue
            else:
                pass

            # if an intron in at least a single organism is still there,
            # then split the involved pacbps in the `original` cbgL, the last
            # added CBG element in the return_cbg_list, and make a (virtual)
            # deepcopy of a novel cbgL. Both CBGs have actually the SAME pacbps!
            cbgR = self.deepcopy()
            cbgL = self.deepcopy()

            # loop over the organisms/genes with inframe introns split
            # the Pacbps of these orgs in both to-become L and R CBGs 
            inframe_intron_orgs = putative_inframe_introns.keys()
            for org in inframe_intron_orgs:
                ################################################################
                if verbose:
                    print "splitting PACBPs for org:", org
                    print "L", cbgL
                    print "R", cbgL
                ################################################################
                node = self.node_by_organism(org)
                replacementsL = {}
                replacementsR = {}
                for (key,node1,node2), pacbporf in cbgL.pacbps.iteritems():
                    if node in [node1,node2]:
                        # get the pacbp of this pacbporf and split it!
                        pacbp = pacb.conversion.pacbporf2pacbp(pacbporf)
                        org1 = self.organism_by_node(node1)
                        org2 = self.organism_by_node(node2)

                        if org1 in putative_inframe_introns.keys() and\
                        org2 in putative_inframe_introns.keys() and\
                        inframe_intron_orgs.index(org) > 0:
                            # already splitted; both orgs are inframe introns!
                            continue

                        # make split coordinates relative
                        splitL = lsrcoords[org1][0] - pacbp.query_start
                        splitR = lsrcoords[org1][1] - pacbp.query_start

                        pacbpL = pacb.splitting.split_pacb_on_coordinates(
                            pacbp,(splitL,splitL),returnside='left')
                        pacbpR = pacb.splitting.split_pacb_on_coordinates(
                            pacbp,(splitR,splitR),returnside='rigth')

                        # check if both cbgL and cbgR make sence
                        # if not -> return False!
                        if not pacbpL: return False
                        if not pacbpR: return False

                        ########################################################
                        if verbose:
                            print "#", node1, node2, lsrcoords[org1], 
                            print "L:", splitL, "R:", splitR
                            print pacbp
                            print pacbpL
                            print pacbpR
                        ########################################################

                        # pacbpL -> extented pacbporfL -> store to replacementsL
                        newpacbporfL = pacb.conversion.pacbp2pacbporf(pacbpL,
                                       pacbporf.orfQ,pacbporf.orfS)
                        newpacbporfL.extend_pacbporf_after_stops()
                        replacementsL[(key,node1,node2)] = newpacbporfL

                        # pacbpR -> extented pacbporfR -> store to replacementsR
                        newpacbporfR = pacb.conversion.pacbp2pacbporf(pacbpR,
                                       pacbporf.orfQ,pacbporf.orfS)
                        newpacbporfR.extend_pacbporf_after_stops()
                        replacementsR[(key,node1,node2)] = newpacbporfR


                # do the pacbporf replacements in both CBGs
                statusL = _update_cbg_with_pacbporf_replacements(
                            cbgL,replacementsL)
                statusR = _update_cbg_with_pacbporf_replacements(
                            cbgR,replacementsR)

                # check if both cbgL and cbgR make sence
                if not statusL or not statusR:
                    # return unchanged cbg status -> False
                    return False
                    


            # Verify the interface between cbgL and cbgR.
            # Most likely, the sites are nicely alignable.
            cbgIF = CodingBlockGraphInterface(cbgL,cbgR)
            cbgIF.force_intron_in_organisms( putative_inframe_introns.keys() )
            cbgIF.allow_intron_in_organisms( putative_inframe_introns.keys() )
            cbgIF.harvest_splice_sites()
            cbgIF.find_conserved_splice_sites()

            ####################################################################
            if verbose:
                print cbgL
                print cbgIF
                print cbgR
                cbgIF.interfaceproperties()
            ####################################################################
            # check the properties of the CBGinterface
            if cbgIF.optimalitycheck().count(True) >= 2:
                # yes; is_compatible and donor and/or acceptor is optimal
                cbgL._CBGinterface3p = cbgIF
                cbgR._CBGinterface5p = cbgIF
                cbgL.copy_5pcbginterface_from_othercbg(self)
                cbgR.copy_3pcbginterface_from_othercbg(self)
                return_cbg_list = [ cbgL, cbgR ]
                ################################################################
                if verbose: print "INFRAME INTRON CONFIRMED!!"
                ################################################################
            else:
                # no compatible interface... although intron(s) was/were found!
                # (at least) two options are now open:
                # 1. enforce the intron(s) and create cbgIF with _forced_ends
                # 2. ignore the intron(s) and create an intermediate lsrCBG

                # 1. is `tricky`. First, how sure is this inframe intron,
                # what type of criteria do we assume etc etc.
                # second, how to create a coorect cbgIF? It must be an
                # IS_SPLITTED interface, of which the boundaries might fall
                # outside the OMSR's of the CBGs.

                # 2. ignore the intron(s) and create an intermediate lsrCBG
                lsrCBG = create_intermediate_lowsimilarity_region(cbgL,cbgR)
                prepare_lsrcbg_and_cbg_for_gsg_insertion(cbgL,lsrCBG)
                prepare_lsrcbg_and_cbg_for_gsg_insertion(lsrCBG,cbgR)
                cbgL.copy_5pcbginterface_from_othercbg(self)
                cbgR.copy_3pcbginterface_from_othercbg(self)
                return_cbg_list = [ cbgL, lsrCBG, cbgR ]
                ################################################################
                if verbose:
                    print "no INFRAME INTRON -> lsrCBG"
                    print cbgL
                    print " ", lsrCBG._CBGinterface5p
                    print " ", lsrCBG
                    print " ", lsrCBG._CBGinterface3p
                    print cbgR
                    self.printmultiplealignment()
                    print cbgL
                    cbgL.printmultiplealignment()
                    print cbgR
                    cbgR.printmultiplealignment()
                ################################################################

        # EOF this function.
        # return False if this CBG remained intact, list of splits when splitted
        if len(return_cbg_list) == 1:
            return False
        else:
            return return_cbg_list
                
# end of function cbg_cexpander_inframe_intron_search


def cexpander_lsrcbg_intron_search(self,prev,next,verbose=False):
        """
        Scan the lsrCBGs in the genestructure for hidden inframe introns

        @type  self: LowSimilarityRegionCodingBlockGraph
        @param self: LowSimilarityRegionCodingBlockGraph (lsrCBG) instance

        @type  prev: CodingBlockGraph
        @param prev: CodingBlockGraph before / 5p this lsrCBG

        @type  next: CodingBlockGraph
        @param next: CodingBlockGraph after / 3p this lsrCBG

        @type  verbose: Boolean
        @param verbose: print status/debugging messages to STDOUT

        @rtype  lsrcbgs_removed_cnt: integer
        @return lsrcbgs_removed_cnt: number of lsrCBGS that have inframe introns
        """
        pass

        # obtain non-uniformly aligned AA lengths for the 
        lengths = self.overall_minimal_spanning_range_sizes()

        ########################################################################
        if verbose: print self, "\n", lengths
        ########################################################################

        # from here on, identical to CBG inframe intron search!
        TODO = True



# end of function cexpander_lsrcbg_intron_search




def _length_discrepancy_to_potential_inframe_introns(
    datadict, min_intron_nt_length = 30):
    """
    Return the organisms/genes (dict keys) for which an inframe introns
        improves the overall gene lengths

    @type  min_intron_nt_length: integer
    @param min_intron_nt_length: MIN_INTRON_NT_LENGTH

    @rtype:  list
    @return: list of the keys of datadict that are potential inframe introns
    """

    tmp = [ (v,k) for k,v in datadict.iteritems() ]
    tmp.sort()
    tmp.reverse()
    keys    = [ k for (v,k) in tmp ]
    lengths = [ v for (v,k) in tmp ]

    inframe_intron_cnt_range = range(0,4)
    datacontainer = []
    for inframe_intron_cnt in inframe_intron_cnt_range:
        curlengths = deepcopy(lengths)
        for i in range(0,inframe_intron_cnt):
            cur_intron_length = curlengths[i] - sum(curlengths[i+1:])/(len(curlengths)-i)
            curlengths[i] = curlengths[i] - max([cur_intron_length,min_intron_nt_length])
        maxdif = float(max(curlengths)) - min(curlengths)
        avdif  = float( sum(curlengths) - min(curlengths)*len(curlengths) ) / len(curlengths)
        if not datacontainer:
            # inframe_intron_cnt == 0 -> first entry
            datacontainer.append( ( curlengths, maxdif, avdif, inframe_intron_cnt ) )

        else:
            # check if it is an improvement
            improvement = True
            for entry in datacontainer:
                if not (maxdif < entry[1] and avdif < entry[2]):
                    improvement = False
                    break
            if improvement:
                # potential inframe intron here!!
                datacontainer.append( ( curlengths, maxdif, avdif, inframe_intron_cnt ) )

    # return the inframe_intron_cnt of the last added row in 
    # the datacontainer; return the corresponding headers
    return keys[0:datacontainer[-1][-1]]

# end of function _length_discrepancy_to_potential_inframe_introns

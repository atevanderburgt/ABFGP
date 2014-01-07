"""
Functions for adding CBGs to the GSG, gathered in a class
used in Aligment Based Gene Predictions
"""

# Module metadata
__authors__ = "Ate van der Burgt"
__license__ = "MIT"

# Python Imports
from sets import Set
from copy import deepcopy

# graphAbgp imports
import ordering
from codingblock_ordering import relatively_positioned_towards
import codingblock_splitting

# Global variables
from settings.codingblockgraph import *
from settings.genestructure import *
from settings.inframeintron import *


class AddCodingBlockGraphFunctions:
    """
    """

    def add_ksminx_codingblocks(self,ksminxcbglist,
        max_cbg_gtg_topo_dif=MAX_KSMINX_CBG_GTG_TOPO_DIF,
        max_cbg_gtg_abs_dif=MAX_KSMINX_CBG_GTG_ABS_DIF,
        min_cbg_gtg_id_ratio=MIN_KSMINX_CBG_GTG_ID_RATIO,
        min_tcode_omsr=None,
        omit_conditional_addition=False,
        max_num_upstream_ksminxcbgs   = 3, # 5p/left  additional K(s-x) CBGs
        max_num_downstream_ksminxcbgs = 2, # 3p/right additional K(s-x) CBGs
        verbose=False):
        """
        Place K(s-x) CBGs into the GSG

        @type  ksminxcbglist: list
        @param ksminxcbglist: list of K(s-x) CBGs to be added to the genestructure

        @type  max_cbg_gtg_topo_dif: float (or None)
        @param max_cbg_gtg_topo_dif:

        @type  max_cbg_gtg_abs_dif: float (or None)
        @param max_cbg_gtg_abs_dif:

        @type  max_cbg_gtg_id_ratio: float (or None)
        @param max_cbg_gtg_id_ratio:

        @type  min_tcode_omsr: float (or None)
        @param min_tcode_omsr:

        @type  omit_conditional_addition: Boolean
        @param omit_conditional_addition: omit all CBG quality checks; just
                verify if the CBG is placeable based on the OMSR coordinates

        @type  max_num_upstream_ksminxcbgs: integer
        @param max_num_upstream_ksminxcbgs: number of K(s-x) CBGs to allow in
            front of the first K(s) CBGS

        @type  max_num_downstream_ksminxcbgs: integer
        @param max_num_downstream_ksminxcbgs: number of K(s-x) CBGs to allow
            after of the final K(s) CBGS

        @type  verbose: Boolean
        @param verbose: print status/debugging messages to STDOUT

        @rtype:  integer
        @return: number of K(s-x) CBGS that is added to the GSG
        """
        # current number of CBGs
        current_cbg_count = len(self)

        # get Set -> list of K(s-x) CBG node counts
        ksmin_cbg_sizes = list(Set([cbg.node_count() for cbg in ksminxcbglist]))
        ksmin_cbg_sizes.sort()
        ksmin_cbg_sizes.reverse()

        # loop over the ksminxsize; the more nodes in a K(s-x) CBG graph,
        # the more likely it is to be a bona fide one!
        for ksminxsize in ksmin_cbg_sizes:
            ksminx_sized_cbg_list = []
            for cbg in ksminxcbglist:
                if cbg.node_count() == ksminxsize:
                    ksminx_sized_cbg_list.append(cbg)

            # place the K(s-x) CBGs of this ksminxsize into the GSG
            placed = self.add_codingblocks(ksminx_sized_cbg_list,
                    max_cbg_gtg_topo_dif=max_cbg_gtg_topo_dif,
                    max_cbg_gtg_abs_dif=max_cbg_gtg_abs_dif,
                    min_cbg_gtg_id_ratio=min_cbg_gtg_id_ratio,
                    min_tcode_omsr=min_tcode_omsr,
                    omit_conditional_addition=omit_conditional_addition,
                    log_debug = verbose
                    )

        ########################################################################
        if verbose:
            print "add_ksminx_codingblocks:", current_cbg_count, "K(s)",
            print len(self)-current_cbg_count, "K(s-x)", self, 
            print len(ksminxcbglist), "K(s-x) candidates"
        ########################################################################

        # find first K(s) CBG
        firstKsCBGpos = None
        for pos in self.cbgpositerator():
            if self.codingblockgraphs[pos].node_count() ==\
            self.EXACT_SG_NODE_COUNT:
                firstKsCBGpos = pos
                break

        # find final K(s) CBG
        finalKsCBGpos = None
        for pos in self.cbgpositerator(reversed=True):
            if self.codingblockgraphs[pos].node_count() ==\
            self.EXACT_SG_NODE_COUNT:
                finalKsCBGpos = pos
                break

        ########################################################################
        if verbose: print "K(s-x) CBG pos:", firstKsCBGpos, finalKsCBGpos, len(self)
        ########################################################################

        # remove K(s-x) CBGs that exceed max_num_downstream_ksminxcbgs
        if len(self)-1-finalKsCBGpos > max_num_downstream_ksminxcbgs:
            # remove all K(s-x) CBGs after finalKsCBGpos + correction
            self.codingblockgraphs =\
            self.codingblockgraphs[0:finalKsCBGpos+max_num_downstream_ksminxcbgs+1]

        # remove K(s-x) CBGs that exceed max_num_upstream_ksminxcbgs
        if firstKsCBGpos > max_num_upstream_ksminxcbgs:
            # remove all K(s-x) CBGs before firstKsCBGpos - correction
            self.codingblockgraphs =\
            self.codingblockgraphs[firstKsCBGpos-max_num_upstream_ksminxcbgs:]

        ########################################################################
        if verbose: print "K(s-x) + K(s) GSG:", self
        ########################################################################

        if len(self) > current_cbg_count:
            # set IS_FIRST / IS_LAST correctly
            self.finalize_genestructure()

        # return added number of K(s-x) CBGs
        return len(self) -  current_cbg_count

    # end of function add_ksminx_codingblocks


    def add_codingblocks(self,cbglist,log_debug=False,
        max_cbg_gtg_topo_dif=None,
        max_cbg_gtg_abs_dif=None,
        min_cbg_gtg_id_ratio=None,
        omit_conditional_addition=False,
        min_tcode_omsr=None ):
        """
        (Try to) add multiple CodingBlockGraphs to the genestructure

        @type  cbglist: list
        @param cbglist: list of CodingBlockGraph objects to be added to the genestructure

        @type  max_cbg_gtg_topo_dif: float (or None)
        @param max_cbg_gtg_topo_dif:

        @type  max_cbg_gtg_abs_dif: float (or None)
        @param max_cbg_gtg_abs_dif:

        @type  max_cbg_gtg_id_ratio: float (or None)
        @param max_cbg_gtg_id_ratio:

        @type  min_tcode_omsr: float (or None)
        @param min_tcode_omsr:

        @type  omit_conditional_addition: Boolean
        @param omit_conditional_addition: omit all CBG quality checks; just
                verify if the CBG is placeable based on the OMSR coordinates

        @rtype:  Boolean
        @return: True or False, weather or not at least 1 CodingBlockGraph was added
        """
        # add cbg's after ordering on total weight
        statusses = []
        # update edge weights to minimal spanning range
        for sg in cbglist:
            sg.update_edge_weights_by_minimal_spanning_range()
        
        for cbg in ordering.order_cbgraphlist(cbglist):
            status = self.add_codingblock(cbg,log_debug=log_debug,
                    max_cbg_gtg_topo_dif=max_cbg_gtg_topo_dif,
                    max_cbg_gtg_abs_dif=max_cbg_gtg_abs_dif,
                    min_cbg_gtg_id_ratio=min_cbg_gtg_id_ratio,
                    min_tcode_omsr=min_tcode_omsr,
                    omit_conditional_addition=omit_conditional_addition,
                    )
            statusses.append( status ) 
        if True in statusses:
            return True
        else:
            return False
        
    # end of function add_codingblocks


    def add_codingblock(self,new,log_debug=False,only_try_adding=False,
        max_cbg_gtg_topo_dif=None,
        max_cbg_gtg_abs_dif=None,
        min_cbg_gtg_id_ratio=None,
        min_tcode_omsr=None,
        omit_conditional_addition=False ):
        """
        (Try to) add a CodingBlockGraph to the genestructure

        @type  new: CodingBlockGraph
        @param new: CodingBlockGraph object to be added to the genestructure

        @type  only_try_adding: Boolean
        @param only_try_adding: only try to add the CBG and return succes status

        @type  omit_conditional_addition: Boolean
        @param omit_conditional_addition: omit all CBG quality checks; just
                verify if the CBG is placeable based on the OMSR coordinates

        @type  max_cbg_gtg_topo_dif: float (or None)
        @param max_cbg_gtg_topo_dif:

        @type  max_cbg_gtg_abs_dif: float (or None)
        @param max_cbg_gtg_abs_dif:

        @type  max_cbg_gtg_id_ratio: float (or None)
        @param max_cbg_gtg_id_ratio:

        @type  min_tcode_omsr: float (or None)
        @param min_tcode_omsr:


        @rtype:  Boolean
        @return: True or False, weather or not adding was succesfull
        """

        verbose = log_debug # log_debug must be replaced by verbose....

        # update edge weights by overall minimal spanning range
        if not only_try_adding:
            new.update_edge_weights_by_minimal_spanning_range()

        # check for difference with GeneTreeGraph
        if new.__class__.__name__ == 'CodingBlockGraph' and\
        not omit_conditional_addition and\
        self.genetree() and len(self) >= 1:
            newgtg = new.genetree()

            if new.node_count() == self.EXACT_SG_NODE_COUNT:
                gtg = self.genetree()
            else:
                # new to-be-placed graph misses certain organism node(s)
                completegtg = self.genetree()
                gtg  = deepcopy(completegtg)
                for missingorg in completegtg.organism_set().difference(new.organism_set()):
                    gtg.del_node(missingorg)


            # calculate identity, topological and absolute GTG differences
            cbg_gtg_topo_dif = gtg.graphalignmentdifference( newgtg )
            cbg_gtg_abs_dif  = gtg.absolutegraphalignmentdifference( newgtg )
            cbg_gtg_id_ratio = newgtg.identity() / gtg.identity()

            ####################################################################
            if verbose: print "cbg2gsg", new
            ####################################################################

            # check the identity ratio
            if min_cbg_gtg_id_ratio:
                threshold_min_cbg_gtg_id_ratio = MIN_GTG_ID_RATIO_FUNCTION(min_cbg_gtg_id_ratio,gtg,new)
                ################################################################
                if verbose:
                    print "CUSTOM", threshold_min_cbg_gtg_id_ratio,
                    print min_cbg_gtg_id_ratio
                ################################################################
            else:
                threshold_min_cbg_gtg_id_ratio = MIN_GTG_ID_RATIO_FUNCTION(self.MIN_CBG_GTG_ID_RATIO,gtg,new)
                ################################################################
                if verbose:
                    print "NORMAL", threshold_min_cbg_gtg_id_ratio,
                    print self.MIN_CBG_GTG_ID_RATIO
                ################################################################

            if cbg_gtg_id_ratio < threshold_min_cbg_gtg_id_ratio:
                ################################################################
                if verbose:
                    print "rejected on ID ratio", threshold_min_cbg_gtg_id_ratio,
                    print ">", cbg_gtg_id_ratio
                ################################################################
                return False
            else:
                pass


            # check the relative topological difference
            if max_cbg_gtg_topo_dif:
                threshold_max_cbg_gtg_topo_dif = MAX_GTG_TOPO_DIF_FUNCTION(max_cbg_gtg_topo_dif,gtg,new)
                ################################################################
                if verbose:
                    print "CUSTOM", threshold_max_cbg_gtg_topo_dif,
                    print max_cbg_gtg_topo_dif, gtg.identity(), new.omsrlength()
                ################################################################
            else:
                threshold_max_cbg_gtg_topo_dif = MAX_GTG_TOPO_DIF_FUNCTION(self.MAX_CBG_GTG_TOPO_DIF,gtg,new)
                ################################################################
                if verbose:
                    print "NORMAL", threshold_max_cbg_gtg_topo_dif,
                    print self.MAX_CBG_GTG_TOPO_DIF, gtg.identity(), new.omsrlength()
                ################################################################

            if cbg_gtg_id_ratio >= 1.10:
                # ignore TOPO_DIF check when newgtg.id% >> gtg.id%
                pass 
            elif cbg_gtg_topo_dif > threshold_max_cbg_gtg_topo_dif: 
                ################################################################
                if verbose:
                    print "rejected on TOPO_DIF", cbg_gtg_topo_dif,
                    print ">", threshold_max_cbg_gtg_topo_dif
                ################################################################
                return False
            else:
                pass

            # check the absolute topological difference
            if max_cbg_gtg_abs_dif:
                threshold_max_cbg_gtg_abs_dif = MAX_GTG_ABS_DIF_FUNCTION(max_cbg_gtg_abs_dif,gtg,new)
                ################################################################
                if verbose:
                    print "CUSTOM", threshold_max_cbg_gtg_abs_dif,
                    print max_cbg_gtg_abs_dif, gtg.identity(), new.omsrlength()
                ################################################################
            else:
                threshold_max_cbg_gtg_abs_dif = MAX_GTG_ABS_DIF_FUNCTION(self.MAX_CBG_GTG_ABS_DIF,gtg,new)
                ################################################################
                if verbose:
                    print "NORMAL", threshold_max_cbg_gtg_abs_dif,
                    print self.MAX_CBG_GTG_ABS_DIF
                ################################################################


            if cbg_gtg_id_ratio >= 1.10:
                # ignore TOPO_DIF check when newgtg.id% >> gtg.id%
                pass
            elif cbg_gtg_abs_dif > threshold_max_cbg_gtg_abs_dif:
                ################################################################
                if verbose:
                    print "rejected on ABS_DIF", threshold_max_cbg_gtg_abs_dif,
                    print "<", cbg_gtg_abs_dif
                ################################################################

                return False
            else:
                pass


            # check the Tcode score
            if min_tcode_omsr > new.msr_tcode_score():
                ################################################################
                if verbose:
                    print "rejected on MIN_TCODE_OMSR", min_tcode_omsr, ">",
                    print new.msr_tcode_score() 
                ################################################################
                return False
            else:
                pass 

        else:
            # probably omit_conditional_addition==True
            pass


        # check if exactly this one is already in the genestructure
        if self.is_codingblockgraph_already_in_genestructure(new):
            ####################################################################
            if verbose: print "already in genestructure!!"
            ####################################################################
            return False

        # if LowSimilarityRegionCodingBlockGraph, do a OMSR check
        if new.__class__.__name__ == 'LowSimilarityRegionCodingBlockGraph':
            # get the OMSR of the complete GeneStructure
            omsrGSG = self.overall_minimal_spanning_range()

            for node,omsr in new._omsr.iteritems():
                org = new.organism_by_node(node)
                if omsrGSG[org].intersection(omsr):
                    ############################################################
                    if verbose:
                        print "ALREADY IN GENESTRUCTURE!", "\n", omsr
                        print omsrGSG[org].intersection(omsr)
                    ############################################################
                    return False
            else:
                ################################################################
                if verbose: print "lsrCBG NOT in genestructure's OMSR"
                ################################################################
                # new lsrCBG! just continue with this function
                pass


        for pos in range(0,len(self)):
            cbg = self.codingblockgraphs[pos]
            # do not compare position towards a LowSimilarityRegionCodingBlockGraph
            if cbg.__class__.__name__ == 'LowSimilarityRegionCodingBlockGraph':
                continue

            # get relative positioning data of 2 CBGs
            ( absPosCbg, absPosNew, binPosCbg, binPosNew, orfIdent, posRel, posBin) = relatively_positioned_towards(cbg,new)

            if verbose: print "eval cbg pos", pos, "binarytuples:", binPosCbg, binPosNew

            # Check the positioning in binaryCbgPositioning
            # Required positioning of new codingblock `new` is cbg-->new
            if binPosCbg == (1,0,0) and binPosNew == (0,0,1):
                # The order is new-->cbg; continue
                continue

            elif new.__class__.__name__ == 'LowSimilarityRegionCodingBlockGraph' and new.node_set().difference(cbg.node_set()):
                # binPos comparison assumes correct position, but it is not because node intersection
                # not all nodes in lsrCBG 'new'  are shared in next CBG 'cbg' -> ignore here!
                continue

            elif new.__class__.__name__ == 'LowSimilarityRegionCodingBlockGraph' and binPosCbg == (1,0,0) and binPosNew == (0,0,0):
                # Final check -> is the lsrCBG directly ajacent to the cbg?
                # In case a single orf set is split up in several CBGs (and several lsrCBGs),
                # insertion can give problems when not verifying the distance
                distances = cbg.distance_between_codingblocks(new)
                if list(Set(distances.values())) == [0]:
                    # yes, this is the position where we want to insert
                    pass
                else:
                    # nope, not the correct, diretly adjacent position
                    if verbose:
                        print "# NOT ADDING HERE a lsrCBG????"
                        print cbg
                        print new
                        print distances
                    # continue to the next cbg position
                    continue

                # if only_try_adding, return a succesfull True!
                if only_try_adding: return True

                if verbose:
                    print new, binPosCbg, binPosNew, "inserting on pos %s+1" % pos

                # Add into the ordered GeneStructure!
                tobeadded = new
                tobeadded.create_cache()
                self.codingblockgraphs.insert( pos+1, tobeadded )
                # succesfull insert; return True
                return True



            elif binPosCbg == (0,0,1) and binPosNew == (1,0,0):
                # The order is cbg-->new; this is the required location!


                if new.__class__.__name__ == 'LowSimilarityRegionCodingBlockGraph':
                    # check if directly neighboring this cbg
                    distances = cbg.distance_between_codingblocks(new)
                    print new
                    print distances

                # if only_try_adding, return a succesfull True!
                if only_try_adding: return True

                # Make deepcopy of this CBG in orde to prevent abnormalities
                # when splitting/merging in lateron function
                #tobeadded = deepcopy(new)
                tobeadded = new
                tobeadded.create_cache()
                self.codingblockgraphs.insert( pos, tobeadded )
                if verbose:
                    print "ADDING:", pos, new
                    print  binPosCbg, binPosNew

                # succesfull insert; return True
                return True

            elif binPosCbg == (0,1,0) and binPosNew == (0,1,0):
                # Exactly this codingblock exists already in the genestructure
                return False

            elif binPosCbg in [(1,0,0),(1,1,0)] and binPosNew in [(0,0,1),(0,1,1)]:
                # Check for: binPosCbg ~= (1,?,0) and binPosNew ~= (0,?,1)
                # Some overlap, and the order is new-->cbg
                # So, potential accepted new CBG 1 position AFTER this one!
                # Check if the overlap is compatible with Orf order
                for i in range(0,len(posBin)):
                    (a1,a2,a3), (b1,b2,b3)       = posRel[i]
                    (Ba1,Ba2,Ba3), (Bb1,Bb2,Bb3) = posBin[i]
                    identicalorfs                = orfIdent[i]
                    if _is_compatible_overlap( (a1,a2,a3), (b1,b2,b3), (Ba1,Ba2,Ba3), (Bb1,Bb2,Bb3), identicalorfs ):
                        pass
                    else:
                        if verbose:
                            print "ABOUT TO BREAK THE FORLOOP I:"
                            print (a1,a2,a3), (b1,b2,b3)
                            print (Ba1,Ba2,Ba3), (Bb1,Bb2,Bb3)
                            cbg_overlap_ratio            = float(a2) / float(a1+a2+a3)
                            new_overlap_ratio            = float(b2) / float(b1+b2+b3)
                            overlap_ratio                = max([cbg_overlap_ratio, new_overlap_ratio])
                            print overlap_ratio
                        # No! incorrect positioning -> break
                        break
                else:
                    # EOF forloop nicely reached; compatible new codingblock
                    # The order is cbg-->new; this is the required location!

                    # check how it is positioned towards the other; NO OVERLAP ALOWED HERE!
                    if pos < len(self)-1:
                        next = self.codingblockgraphs[pos+1]
                        ( absPosNext, absPosNewNext, binPosNext, binPosNewNext,
                          orfIdentNext, posRelNext, posBinNext ) = relatively_positioned_towards(new,next)
                        # now check binPosNext & binPosNewNext: should be (1, 0, 0) & (0, 0, 1)
                        # That means -> no overlap with the next CBG
                        if binPosNext == (1, 0, 0) and binPosNewNext == (0, 0, 1):
                            # no overlap -> ready to store!
                            pass
                        else:
                            # there is overlap! Do not allow addition of this CBG
                            ###########################################################
                            if verbose:
                                print "OVERLAP WITH NEXT CBG!:", pos, new
                                print binPosNext, binPosNewNext 
                            ###########################################################
                            return False
                    else:
                        # no further CBGs in GSG, so no check possible (and no overlap possible ;-)
                        pass

                    # if only_try_adding, return a succesfull True!
                    if only_try_adding: return True

                    # Make deepcopy of this CBG in orde to prevent abnormalities
                    # when splitting/merging in lateron function
                    #tobeadded = deepcopy(new)
                    tobeadded = new
                    tobeadded.create_cache()
                    self.codingblockgraphs.insert( pos+1, tobeadded )
                    if verbose:
                        print "ADDING:", pos+1, new
                        print  binPosCbg, binPosNew
                    # succesfull insert; return True
                    return True

                # forloop broken -> incompatible new codingblock
                return False

            elif binPosCbg in [(0,0,1),(0,1,1)] and binPosNew in [(1,0,0),(1,1,0)]:
                # Check for: binPosCbg ~= (0,?,1) and binPosNew ~= (1,?,0)
                # The order is new-->cbg; this is the required location!
                # Check if the overlap is compatible with Orf order
                for i in range(0,len(posBin)):
                    (a1,a2,a3), (b1,b2,b3)       = posRel[i]
                    (Ba1,Ba2,Ba3), (Bb1,Bb2,Bb3) = posBin[i]
                    identicalorfs                = orfIdent[i]
                    if _is_compatible_overlap( (a1,a2,a3), (b1,b2,b3), (Ba1,Ba2,Ba3), (Bb1,Bb2,Bb3), identicalorfs ):
                        pass
                    else:
                        if verbose:
                            print "ABOUT TO BREAK THE FORLOOP II:"
                            print (a1,a2,a3), (b1,b2,b3)
                            print (Ba1,Ba2,Ba3), (Bb1,Bb2,Bb3)
                            cbg_overlap_ratio            = float(a2) / float(a1+a2+a3)
                            new_overlap_ratio            = float(b2) / float(b1+b2+b3)
                            overlap_ratio                = max([cbg_overlap_ratio, new_overlap_ratio])
                            print overlap_ratio
                        # No! incorrect positioning -> break
                        break
                else:
                    # EOF forloop nicely reached; compatible new codingblock
                    # The order is new-->cbg; this is the required location!

                    # if only_try_adding, return a succesfull True!
                    if only_try_adding: return True

                    # Make deepcopy of this CBG in orde to prevent abnormalities
                    # when splitting/merging in lateron function
                    #tobeadded = deepcopy(new)
                    tobeadded = new
                    tobeadded.create_cache()
                    self.codingblockgraphs.insert( pos, tobeadded )
                    # succesfull insert; return True
                    return True

                # forloop broken -> incompatible new codingblock
                return False

            else:
                # A more messy positioning (overlaps/repetitive etc.)
                if verbose:
                    print "WEIRD BINARY SUMMED ORDER!", absPosCbg, absPosNew, binPosCbg, binPosNew
                    print cbg
                    for i in range(0,len(posBin)):
                        print posRel[i], posBin[i], orfIdent[i]
                # Reject this new codingblockgrap
                return False
        else:
            # If eof for loop is reached, append to the end
            # of current genestructure

            # if only_try_adding, return a succesfull True!
            if only_try_adding: return True

            # Make deepcopy of this CBG in orde to prevent abnormalities
            # when splitting/merging in lateron function
            #tobeadded = deepcopy(new)
            tobeadded = new
            tobeadded.create_cache()
            self.codingblockgraphs.append( tobeadded )
            # check if this is the first added cbg to GeneStructure object
            if len(self.codingblockgraphs) == 1:
                # cache the genetreegraph object
                gtg = self.set_genetree()
                # set threshold values as a function of gtg
                self.initialize_first_added_cbg()

            # succesfull insert; return True
            return True

    # end of function add_codingblock


    def is_codingblockgraph_already_in_genestructure(self,cbgraph):
        """
        """
        # make ordered list of novel pacbp nodes
        novel_pacbps = cbgraph.pacbps.keys()
        novel_pacbps.sort()
        novel_nodes = cbgraph.get_ordered_nodes()
        # Stringifying the cbgraph can fail because it has PacbPs iso PacbPORFS
        # In that case, Tcode calculation fails with a raise
        try:    stringified = str(cbgraph)
        except: stringified = None
        for sg in self.codingblockgraphs:
            # if strings of CBGs are identical -> identical CBGs ;-)
            if str(sg) == stringified: return True

            # check if node sets are identical; if not -> not identical ;-)
            if novel_nodes != sg.get_ordered_nodes(): continue

            # do somewhat more elaborate check(s)
            if cbgraph.__class__.__name__ == sg.__class__.__name__:
                if cbgraph.__class__.__name__ == 'CodingBlockGraph':
                    # make ordered list of pacbp nodes
                    known_pacbps = sg.pacbps.keys()
                    known_pacbps.sort()
                    if novel_pacbps == known_pacbps:
                        return True
                else:
                    # not a CodingBlockGraph but a LowSimilarityRegionCodingBlockGraph
                    omsrNovel = [ (k,cbgraph._omsr[k]) for k in cbgraph._omsr.keys() ]
                    omsrKnown = [ (k,sg._omsr[k]) for k in sg._omsr.keys() ]
                    omsrNovel.sort()
                    omsrKnown.sort()
                    if omsrNovel == omsrKnown:
                        return True
            else:
                # unequal graph types can never be identical!
                pass
        else:
            # if eof is reached, return False
            return False

    # end of function is_codingblockgraph_already_in_genestructure


    def is_pacbp_conflicting_with_genestructure(self,pacbp,orgQ=None,orgS=None):
        """
        """
        MIN_PACBP_GENESTRUCTURE_EXTENTION = 6 
        rQ   = Set(range(pacbp.query_start,pacbp.query_end))
        rS   = Set(range(pacbp.query_start,pacbp.query_end))
        pacbplen = len(pacbp)
        if pacbplen - len( self.overall_minimal_spanning_range(orgQ).intersection( rQ ) ) < MIN_PACBP_GENESTRUCTURE_EXTENTION:
            return len( self.overall_minimal_spanning_range(orgQ).intersection( rQ ) )
            return True
        if pacbplen - len( self.overall_minimal_spanning_range(orgS).intersection( rS ) ) < MIN_PACBP_GENESTRUCTURE_EXTENTION:
            return len( self.overall_minimal_spanning_range(orgS).intersection( rS ) )
            return True
        # if we reach this point, compatible position
        return False

    # end of function is_pacbp_conflicting_with_genestructure


# end of class AddCodingBlockGraphFunctions 


################################################################
#### Helper functions                                       ####
################################################################

def _is_compatible_overlap( tupA, tupB, tupBinA, tupBinB, identicalorfs ):
    """
    @attention: global variable CBG_MAX_ALOWED_GSGINSERTION_OVERLAP_AA_LENGTH must be assigned
    @attention: global variable CBG_MAX_ALOWED_GSGINSERTION_OVERLAP_RATIO must be assigned
    """
    overlap_ratio_A = float(tupA[1]) / float(sum(tupA))
    overlap_ratio_B = float(tupB[1]) / float(sum(tupB))
    overlap_ratio   = max( [ overlap_ratio_A, overlap_ratio_B ] )
    overlap_length  = tupA[1] # remind, tupA[1] == tupB[1] !!

    if ( tupBinA, tupBinB ) == ( (1,0,0), (0,0,1) ):
        # no overlap for this organism
        return True
    elif ( tupBinA, tupBinB ) == ( (0,0,1), (1,0,0) ):
        # no overlap for this organism
        return True
    elif (0,1,0) in [ tupBinA, tupBinB ]:
        # a full inclusion -> reject
        return False
    elif identicalorfs and overlap_ratio > CBG_MAX_ALOWED_GSGINSERTION_OVERLAP_RATIO:
        # overlap is way to large!
        return False
    elif identicalorfs and overlap_length <= CBG_MAX_ALOWED_GSGINSERTION_OVERLAP_AA_LENGTH:
        # overlap is fairly small
        return True
    elif not identicalorfs and overlap_length <= CBG_MAX_ALOWED_GSGINSERTION_OVERLAP_AA_LENGTH:
        # Overlap is fairly small BUT not identicalorfs
        # dangerous.....how to join these codingblocks?
        # TODO: think more in detail about this case, maybe break here
        return True
    else:
        # No! incorrect positioning
        return False

# end of function _is_compatible_overlap




def findmostlikelyCBG2GSGinsert(partialGSG,cbglist,
    order_cbgs_by_total_weight=True,
    reorder_cbgs_on_node_occurrence=True,
    optimizetinyexoninterface=True,
    verbose=False):
    """
    """
    # if no CBG in list -> done!
    if not cbglist: return partialGSG

    # make interfaces in the partialGSG
    created = partialGSG.create_cbginterfaces()
    gsgsize = len(partialGSG)

    if order_cbgs_by_total_weight:
        # order graphs by total weight
        cbglist = ordering.order_cbgraphlist(cbglist)

    if reorder_cbgs_on_node_occurrence:
        # get prev and next cbg from partialGSG
        if len(partialGSG) >= 2:
            prev = partialGSG.codingblockgraphs[0]
            next = partialGSG.codingblockgraphs[-1]
        elif len(partialGSG) == 1:
            prev, next = None, partialGSG.codingblockgraphs[0]
        else:
            prev, next = None, None
        # re-order on node occurrence.
        cbglist  = ordering.reorder_cbgs_on_node_occurrence(cbglist,prev=prev,next=next)

    ####################################################################
    if verbose and partialGSG and cbglist:
        print "elegiable CBGs for addition in partialGSG"
        for cbg in cbglist: print cbg
    ####################################################################

    # now try to insert CBGs in the partialGSG
    # we try several insertions, starting with the most stringest one

    ############################################################
    # (1)   only do both options for the top listed CBG     ####
    #       perfect cbgIF(s), topologically sound CBGs      ####
    #       perfect cbgIF(s), CBG2GTG topology not tested   ####
    ############################################################
    is_added = _place_cbg_in_partialgsg([ cbglist[0] ],partialGSG,
            omit_conditional_addition=False,
            optimizetinyexoninterface=optimizetinyexoninterface,
            verbose=verbose )
    is_added = _place_cbg_in_partialgsg([ cbglist[0] ],partialGSG,
            omit_conditional_addition=True,
            optimizetinyexoninterface=optimizetinyexoninterface,
            verbose=verbose )

    ############################################################
    # (2)   perfect cbgIF(s), topologically sound CBGs      ####
    ############################################################
    is_added = _place_cbg_in_partialgsg(cbglist,partialGSG,
            omit_conditional_addition=False,
            optimizetinyexoninterface=optimizetinyexoninterface,
            verbose=verbose )

    ############################################################
    # (3)   perfect cbgIF(s), CBG2GTG topology not tested   ####
    ############################################################
    is_added = _place_cbg_in_partialgsg(cbglist,partialGSG,
            omit_conditional_addition=True,
            optimizetinyexoninterface=optimizetinyexoninterface,
            verbose=verbose )

    # when partialGSG is not enlarged, its cbgIFs are likely
    # distorted in the proces of trying novel CBG inserts.
    # Recreate the cbgIFs here
    if len(partialGSG) == gsgsize:
        if gsgsize >=2:
            partialGSG.clear_central_cbginterfaces()
            partialGSG.create_cbginterfaces()
            if gsgsize == 2:
                # most obvious case of 2 CBGs in the GSG
                partialGSG.codingblockgraphs[0].IS_3P_SPLITTED = False
                partialGSG.codingblockgraphs[1].IS_5P_SPLITTED = False
                if not partialGSG.codingblockgraphs[0].IS_5P_SPLITTED:
                    partialGSG.codingblockgraphs[0].IS_SPLITTED = False
                if not partialGSG.codingblockgraphs[1].IS_3P_SPLITTED:
                    partialGSG.codingblockgraphs[1].IS_SPLITTED = False
        else:
            # a partialGSG of only a single CBG (first or last CBG)
            partialGSG.codingblockgraphs[0]._CBGinterface5p = None
            partialGSG.codingblockgraphs[0]._CBGinterface3p = None

    # Done with this function. Return the (incremented) partialGSG
    return partialGSG

# end of function findmostlikelyCBG2GSGinsert



def _remove_placed_cbgs_from_list(placed,cbglist):
    """
    Remove the CBGs from cbglist by the integer list positions given in a list

    @type  placed: []
    @param placed: list with position IDs referring to the positions in cbglist

    @type  cbglist: []
    @param cbglist: list with CodingBlockGraphs

    @attention: used only from within function _place_cbg_in_partialgsg()
    """ 
    # remove the CBGs that are placed in the partialGSG
    if placed:
        placed.sort()
        placed.reverse()
        for pos in placed:
            cbglist.pop(pos)

# end of function _remove_placed_cbgs_from_list


class CBGnotinGSG(Exception):
    pass

def _cbg_position_in_gsg(cbg,gsg):
    """
    Find the (integer) list position of a cbg in a GeneStructureOfCodingBlockGraphs

    @type  cbg: CodingBlockGraph
    @param cbg: CodingBlockGraph that might be present in the gsg

    @type  gsg: GeneStructureOfCodingBlockGraphs
    @param gsg: GeneStructureOfCodingBlockGraphs that might contain the cbg 

    @rtype:  integer or CBGnotinGSG(Exception) 
    @return: position of the cbg in the gsg.codingblockgraphs list
    """
    cbgstr = str(cbg)
    for i in range(0,len(gsg)):
        if cbgstr == str(gsg.codingblockgraphs[i]):
            return i
    else:
        return CBGnotinGSG

# end of function _cbg_position_in_gsg



def _place_cbg_in_partialgsg(cbglist,partialGSG,
    optimizetinyexoninterface=True,
    omit_conditional_addition=False,
    verbose=False):
    """
    @type  cbglist: [] 
    @param cbglist: list with CodingBlockGraphs that might be palced in the (partial)GSG

    @type  partialGSG: GeneStructureOfCodingBlockGraphs
    @param partialGSG: partial GeneStructureOfCodingBlockGraphs in which the CBGs are tried to be inserted into

    @rtype:  Boolean 
    @return: are any CBGs from cbglist placed into partialGSG?
    """
    # import function here to prevent circular import
    # TODO: make correct import!
    from genestructure_intermediatecbg import intermediateCBG_node_comparison

    placed_in_partialGSG = []
    curgsglen = len(partialGSG)
    while cbglist:
        for i in range(0,len(cbglist)):
            cbg = cbglist[i]
            # placeability check in the GSG with function's settings
            # for topological check or not
            placeability = partialGSG.add_codingblock(cbg,
                    only_try_adding=True,
                    omit_conditional_addition=omit_conditional_addition
                    )
            ############################################################
            if verbose: print i, cbg, placeability
            ############################################################
            if not placeability: continue

            # place in the partialGSG and find the inserted position
            added = partialGSG.add_codingblock(cbg,omit_conditional_addition=omit_conditional_addition)
            cbgposingsg = _cbg_position_in_gsg(cbg,partialGSG)

            # do intermediateCBG_node_comparison() in the insert position
            if cbgposingsg > 0 and cbgposingsg < len(partialGSG)-1:
                prevCBG = partialGSG.codingblockgraphs[cbgposingsg-1]
                nextCBG = partialGSG.codingblockgraphs[cbgposingsg+1]
                if False == intermediateCBG_node_comparison(prevCBG,cbg,nextCBG):
                    # erroneou CBG insert -> continue
                    continue

            # replace proper pacbporfs from the parents
            if cbgposingsg > 0:
                prevCBG = partialGSG.codingblockgraphs[cbgposingsg-1]
                replacements1 = partialGSG.codingblockgraphs[cbgposingsg]._recrute_pacbporfs_from_parental_cbg(prevCBG,verbose=verbose)
            if cbgposingsg < len(partialGSG)-1:
                nextCBG = partialGSG.codingblockgraphs[cbgposingsg+1]
                replacements2 = partialGSG.codingblockgraphs[cbgposingsg]._recrute_pacbporfs_from_parental_cbg(nextCBG,verbose=verbose)

            # create cbgIFs
            created = partialGSG.create_cbginterfaces()

            # check if one of the direct neighbouring CBGs has
            # the same set of nodes -> signal for a lsrCBG
            cbginterface_isa_lsrcbg = False
            cbginterface_isa_lsrcbg_asses_cbgIFa = None 
            cbginterface_isa_lsrcbg_asses_cbgIFb = None
            lsrCBG = None
            if cbgposingsg > 0 and len(cbg.node_set().intersection(
            partialGSG.codingblockgraphs[cbgposingsg-1].node_set())) == cbg.node_count():
                # left/5p of cbg is a CBG in the partialGSG with identical node set
                cbginterface_isa_lsrcbg = True
                cbginterface_isa_lsrcbg_asses_cbgIFa = False
                cbginterface_isa_lsrcbg_asses_cbgIFb = True 
                lsrCBG = codingblock_splitting.create_intermediate_lowsimilarity_region(
                        partialGSG.codingblockgraphs[cbgposingsg-1], cbg )
            if cbgposingsg < len(partialGSG)-1 and len(cbg.node_set().intersection(
            partialGSG.codingblockgraphs[cbgposingsg+1].node_set())) == cbg.node_count():
                # right/3p of cbg is a CBG in the partialGSG with identical node set
                cbginterface_isa_lsrcbg = True
                cbginterface_isa_lsrcbg_asses_cbgIFa = True 
                cbginterface_isa_lsrcbg_asses_cbgIFb = False 
                lsrCBG = codingblock_splitting.create_intermediate_lowsimilarity_region(
                        cbg, partialGSG.codingblockgraphs[cbgposingsg+1] )

            # assess the created CBGinterfaces
            cbgIFa = partialGSG.codingblockgraphs[cbgposingsg]._CBGinterface5p
            cbgIFb = partialGSG.codingblockgraphs[cbgposingsg]._CBGinterface3p
            if cbgIFa:
                if optimizetinyexoninterface: cbgIFa.optimizetinyexoninterface()
                cbgIFaCheck = cbgIFa.optimalitycheck()
            else:
                cbgIFaCheck = [ None, None, None ]
            if cbgIFb:
                if optimizetinyexoninterface: cbgIFb.optimizetinyexoninterface()
                cbgIFbCheck = cbgIFb.optimalitycheck()
            else:
                cbgIFbCheck = [ None, None, None ]


            # check if this freshly placed CBG makes sense to place in the partialGSG
            if cbginterface_isa_lsrcbg and lsrCBG:
                is_lsrcbg_addable = True
                if cbgIFa and cbginterface_isa_lsrcbg_asses_cbgIFa and\
                cbgIFaCheck.count(True) < 2:
                    is_lsrcbg_addable = False
                if cbgIFb and cbginterface_isa_lsrcbg_asses_cbgIFb and\
                cbgIFbCheck.count(True) < 2:
                    is_lsrcbg_addable = False
                if is_lsrcbg_addable:
                    # addable in the partialGSG; add the lsrCBG too
                    added = partialGSG.add_codingblock(lsrCBG,omit_conditional_addition=True)
                    # create cbgInterfaces for the novel added lsrCBG
                    partialGSG.create_cbginterfaces()
                    placed_in_partialGSG.append(i)
                    # continue trying adding the next CBG
                    continue
                else:
                    # nope, not addable; pass here and solve
                    # removal of this (lsr)CBG lateron
                    pass

            elif cbginterface_isa_lsrcbg and not lsrCBG:
                is_lsrcbg_addable = True
                if cbgIFa and cbginterface_isa_lsrcbg_asses_cbgIFa and\
                cbgIFaCheck.count(True) < 2:
                    is_lsrcbg_addable = False
                if cbgIFb and cbginterface_isa_lsrcbg_asses_cbgIFb and\
                cbgIFbCheck.count(True) < 2:
                    is_lsrcbg_addable = False
                if is_lsrcbg_addable:
                    # addable in the partialGSG as a tight fit to existing CBGs
                    # without a lsrCBG. Weird & rare case but can happen!
                    # re-create cbgInterface for the novel added CBG because
                    # it is not recognized asa splitted interface yet
                    partialGSG.codingblockgraphs[cbgposingsg]._CBGinterface5p = None
                    partialGSG.codingblockgraphs[cbgposingsg]._CBGinterface3p = None
                    if cbgposingsg > 0:
                        partialGSG.codingblockgraphs[cbgposingsg-1]._CBGinterface3p = None
                    if cbgposingsg < len(partialGSG)-1: 
                        partialGSG.codingblockgraphs[cbgposingsg+1]._CBGinterface5p = None
                    # now recreate cbgInterfaces
                    partialGSG.create_cbginterfaces()
                    placed_in_partialGSG.append(i)
                    # continue trying adding the next CBG
                    continue
                else:
                    # nope, not addable; pass here and solve
                    # removal of this (lsr)CBG lateron
                    pass

            elif cbgIFa and cbgIFb:
                if cbgIFaCheck.count(True) >= 2 or cbgIFbCheck.count(True) >= 2:
                    ############################################################
                    if verbose: print "PLACEDab\n", cbgIFa, "\n", cbg, "\n", cbgIFb
                    ############################################################
                    # succesfully placed; leave in place
                    placed_in_partialGSG.append(i)
                    # continue trying adding the next CBG
                    continue

            elif cbgIFa:
                if cbgIFaCheck.count(True) >= 2:
                    ############################################################
                    if verbose: print "PLACEDa\n", cbgIFa, "\n", cbg
                    ############################################################
                    # succesfully placed; leave in place
                    placed_in_partialGSG.append(i)
                    # continue trying adding the next CBG
                    continue

            elif cbgIFb:
                if cbgIFbCheck.count(True) >= 2:
                    ############################################################
                    if verbose: print "PLACEDb\n", cbg, "\n", cbgIFb
                    ############################################################
                    # succesfully placed; leave in place
                    placed_in_partialGSG.append(i)
                    # continue trying adding the next CBG
                    continue

            else:
                # what else!?
                raise "No cbgIFs at all in partialGSG %s for cbg %s" % ( partialGSG, cbg )


            ############################################################
            if verbose:
                print i, "NOTPLACABLE!",
                print cbg
                print "cbgIFa:", cbgIFa
                print "cbgIFb:", cbgIFb
            ############################################################

            # Remove the falsely placed CBG and recreate original cbgIFs
            partialGSG.codingblockgraphs.pop(cbgposingsg)
            created = partialGSG.create_cbginterfaces()
            # done with trying to add this CBGS. Do next...

        if placed_in_partialGSG:
            # remove the CBGs that are placed in the partialGSG
            _remove_placed_cbgs_from_list(placed_in_partialGSG,cbglist)
            # reset placed_in_partialGSG to empty list
            placed_in_partialGSG = []
        else:
            # no cbgs placed in the GSG -> break the while loop
            break

    # check if a CBG was added 
    if len(partialGSG) > curgsglen:
        return True
    else:
        return False

# end of function _place_cbg_in_partialgsg


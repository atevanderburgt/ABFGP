"""
Functions for removing CBGs from the GSG
"""

# graphAbgp Imports
import codingblock_ordering
from graph_codingblock import estimate_bitscores_ksminxcbg

# Global variable Imports
from settings.codingblockgraph import (
    CBG_MAX_ALOWED_GSGREMOVAL_OVERLAP_AA_LENGTH,
    CBG_MAX_ALOWED_GSGREMOVAL_OVERLAP_RATIO,
    KSMINX_CBG_SINGLE_MISSING_UNLIKELY_EDGE_BITSCORE, 
    KSMINX_CBG_MULTIPLE_MISSING_UNLIKELY_EDGE_BITSCORE,
    )


class CodingBlockGraphRemovalFunctions:

    def separate_ds_and_us_gsg(self):
        """
        Remove all CBGs that have status IS_IGNORED in separate GSGs

        @attention: LowSimilarityRegionCodingBlockGraph that are not required
                    anymore are deleted!

        @attention: CodingBlockGraphInterface objects around removed CBGs are set to None!

        @rtype:  tuple
        @return: tuple of 3 (empty) GenestructureOfCodingBlockGraphs
        """
        # create empty GSG to place IS_IGNORED CBGs in
        from graph_genestructure import GenestructureOfCodingBlockGraphs
        dsGSG  = GenestructureOfCodingBlockGraphs(self.input)
        usGSG  = GenestructureOfCodingBlockGraphs(self.input)
        etcGSG = GenestructureOfCodingBlockGraphs(self.input)

        # if no CBGs in GSG -> return all empty ones
        if len(self) == 0: return dsGSG, usGSG, etcGSG

        # check if there is any LowSimilarityRegionCodingBlockGraph directly
        # next to a CBG that has status IS_IGNORED. If so, set the lsrCBG
        # status to IS_IGNORED too!
        for pos in range(0,len(self)):
            if not self.codingblockgraphs[pos].IS_IGNORED: continue
            thiscbg = self.codingblockgraphs[pos]
            if pos > 0:
                prevclass = self.codingblockgraphs[pos-1].__class__.__name__
                if prevclass == 'LowSimilarityRegionCodingBlockGraph':
                    # set to IS_IGNORED too!
                    self.codingblockgraphs[pos-1].IS_IGNORED = True
            if pos < len(self)-1:
                nextclass = self.codingblockgraphs[pos+1].__class__.__name__
                if nextclass == 'LowSimilarityRegionCodingBlockGraph':
                    # set to IS_IGNORED too!
                    self.codingblockgraphs[pos+1].IS_IGNORED = True

        # Separate a potential downsteam GSG from this main GSG in dsGSG
        # that means, all IS_IGNORED CBGs until the first that is not IS_IGNORED
        for pos in range(0,len(self)):
            if pos == 0 and not self.codingblockgraphs[pos].IS_IGNORED:
                # first CBG is not IS_IGNORED -> no dsGSG
                break
            if not dsGSG and not self.codingblockgraphs[pos].IS_IGNORED and pos >= 1:
                # place all CBGs ds of this first CBG that is not IS_IGNORED in dsGSG
                for delpos in range(0,pos):
                    dsGSG.codingblockgraphs.append( self.codingblockgraphs.pop(0) )
                for cbg in dsGSG.codingblockgraphs:
                    cbg.IS_IGNORED = False
                break

        # if dsGSG is not empty, fix cbgIFs in the GSG itself and check for
        # now non-sense lsrCBGs in the dsGSG
        if len(dsGSG) > 0:
            # whipe out the cbgIF object on the 5p side of the first CBG
            self.codingblockgraphs[0]._CBGinterface5p = None
            self.codingblockgraphs[0]._forced_5p_ends = {}

            # check dsGSG; if on of its exterior CBGs is an lsrCBG, remove it!
            if dsGSG.codingblockgraphs[0].__class__.__name__ ==\
            'LowSimilarityRegionCodingBlockGraph':
                removed_lsrCBG = dsGSG.codingblockgraphs.pop(0)
            if dsGSG.codingblockgraphs[-1].__class__.__name__ ==\
            'LowSimilarityRegionCodingBlockGraph':
                removed_lsrCBG = dsGSG.codingblockgraphs.pop()


        # separate a potential upstream GSG from this main GSG in usGSG
        # that means, all IS_IGNORED CBGs until the first that is not IS_IGNORED
        for pos in range(len(self)-1,-1,-1):
            if pos == len(self)-1 and not self.codingblockgraphs[pos].IS_IGNORED:
                # first CBG is not IS_IGNORED -> no usGSG
                break
            if not usGSG and not self.codingblockgraphs[pos].IS_IGNORED and pos < len(self)-1:
                # place all CBGs ds of this first CBG that is not IS_IGNORED in usGSG
                for delpos in range(pos+1,len(self)):
                    usGSG.codingblockgraphs.insert(0, self.codingblockgraphs.pop() )
                for cbg in usGSG.codingblockgraphs:
                    cbg.IS_IGNORED = False
                break

        # if usGSG is not empty, fix cbgIFs in the GSG itself and check for
        # now non-sense lsrCBGs in the usGSG
        if len(usGSG) > 0:
            # whipe out the cbgIF object on the 3p side of the last CBG
            self.codingblockgraphs[len(self)-1]._CBGinterface3p = None
            self.codingblockgraphs[len(self)-1]._forced_3p_ends = {}

            # check usGSG; if on of its exterior CBGs is an lsrCBG, remove it!
            if usGSG.codingblockgraphs[0].__class__.__name__ ==\
            'LowSimilarityRegionCodingBlockGraph':
                removed_lsrCBG = usGSG.codingblockgraphs.pop(0)
            if usGSG.codingblockgraphs[-1].__class__.__name__ ==\
            'LowSimilarityRegionCodingBlockGraph':
                removed_lsrCBG = usGSG.codingblockgraphs.pop()


        # check for intermediate IS_IGNORED CBGs in the main GSG and place in etcGSG
        for pos in range(len(self)-1,-1,-1):
            if self.codingblockgraphs[pos].IS_IGNORED:
                # place the IS_IGNORED one in etcGSG
                removed_cbg = self.codingblockgraphs.pop(pos)
                classname   = removed_cbg.__class__.__name__
                if classname != 'LowSimilarityRegionCodingBlockGraph':
                    # only place non-lsrCBGs in the etcGSG
                    # intermediate lsrCBGs that are IS_IGNORED are now nonsense
                    etcGSG.codingblockgraphs.insert(0, removed_cbg )

                # whipe out potential cbgIF objects surrounding this CBG
                if pos > 0:
                    self.codingblockgraphs[pos-1]._CBGinterface3p = None
                    self.codingblockgraphs[pos-1]._forced_3p_ends = {}
                if pos <= len(self):
                    self.codingblockgraphs[pos]._CBGinterface5p = None
                    self.codingblockgraphs[pos]._forced_5p_ends = {}

        # and return dsGSG, usGSG, etcGSG
        return dsGSG, usGSG, etcGSG

    # end of function separate_ds_and_us_gsg


    def remove_cbg_by_pos(self,remove_pos):
        """
        """
        # pop the requested CBG from the GSG
        deletedCBG = self.codingblockgraphs.pop(remove_pos)

        # update cbgIF objects -> set to None
        if remove_pos > 0:
            # there is a CBG 5' of the just deleted cbg
            self.codingblockgraphs[remove_pos-1]._CBGinterface3p = None
        if remove_pos < len(self):
            # there is a CBG 3' of the just deleted cbg
            # check for len(self), not len(self)-1 !!! because
            # the `remove_pos` could have been the final cbg...
            self.codingblockgraphs[remove_pos]._CBGinterface5p = None

        # check if neighboring CBG was overlapping with identical nodes
        if deletedCBG.IS_5P_SPLITTED and remove_pos > 0:
            prevCBG = self.codingblockgraphs[remove_pos-1]
            if prevCBG._short_name == 'CBG':
                prevCBG.IS_3P_SPLITTED = False
                if not prevCBG.IS_5P_SPLITTED:
                    prevCBG.IS_SPLITTED = False
        if deletedCBG.IS_3P_SPLITTED and remove_pos < len(self):
            nextCBG = self.codingblockgraphs[remove_pos]
            if nextCBG._short_name == 'CBG':
                nextCBG.IS_5P_SPLITTED = False
                if not nextCBG.IS_3P_SPLITTED:
                    nextCBG.IS_SPLITTED = False

        # check if the neighboring CBG was a lsrCBG -> cascade deletion!
        self.remove_incorrect_lsrcbgs()

        # return the removed CBG
        return deletedCBG

    # end of function remove_cbg_by_pos


    def remove_noncoding_cbgs(self,
        max_cbg_aa_length=15,
        max_omsr_tcode=0.70,
        cbg_optimality_check = True,
        perform_cbgif_check = True,
        recreate_cbgifs = True,
        verbose=True): 
        """
        Remove CBGs that are highly likely non-coding alignments

        @type  verbose: Boolean
        @param verbose: report information/debugging messages to STDOUT

        @type  max_cbg_aa_length: integer (or None)
        @apram max_cbg_aa_length: maximal OMSR length to allow removal

        @rtype:  integer
        @return: number of CBGs that is removed from the GSG
        """
        # loop REVERSED because elements can be removed during this for loop
        removed_cnt = 0
        for pos in range(len(self)-1,-1,-1):
            cbg = self.codingblockgraphs[pos]
            # no not take lsrCBGs into account
            if cbg.__class__.__name__ == 'LowSimilarityRegionCodingBlockGraph': continue
            # check if CBG has properties that allow a deletion
            if max_cbg_aa_length and cbg.omsrlength() > max_cbg_aa_length: continue
            if cbg.omsr_tcode_score() > max_omsr_tcode: continue

            # only allow incompatible CBGs to be deleted
            if cbg_optimality_check:
                statuslist = [ self._codingblock_prediction_status(cbg,org) for org in cbg.organism_set() ] 
                if None in statuslist:
                    # no data on cbg optimality yet (no cbgIFs) -> omit deletion
                    continue
                elif not False in statuslist:
                    # this CBG is optimal (no grounds to have it deleted)
                    continue
                else:
                    # CBG still elegiable for deletion...
                    pass

            # only allow CBGs with poor cbgIFs to be deleted
            if perform_cbgif_check:
                # perform_cbgif_check requires presence of cbgIFs.
                # When an interface without a cbgIF is encountered,
                # it is omitted for deletion. Initialize None values
                # for the cbgIF properties. When None is maintained,
                # it means that it is the first / last CBG in the GSG
                c5p,d5p,a5p = None, None, None
                c3p,d3p,a5p = None, None, None
                if pos > 0:
                    try:
                        c5p,d5p,a5p = cbg._CBGinterface5p.cbg_optimality_check()
                    except:
                        # no cbgIF given!  -> omit deletion
                        c5p,d5p,a5p = True, True, True
                if pos < len(self)-1:
                    try:
                        c3p,d3p,a3p = cbg._CBGinterface3p.cbg_optimality_check()
                    except:
                        # no cbgIF given!  -> omit deletion
                        c3p,d3p,a3p = True, True, True
                # now check the Booleans...
                if None in [c5p,c3p] and [c5p,d5p,a5p,c3p,d3p,a5p].count(True)==3:
                    # first/last CBG and its interface is optimal -> do not delete
                    continue
                elif [c5p,c3p] == [True,True] and (cbg._CBGinterface5p.is_optimal() or\
                cbg._CBGinterface3p.is_optimal()):
                    # intermediate CBG and compatible cbgIF on both sides -> do not delete
                    continue
                else:
                    # all other cases: elegiable for deletion!
                    pass 

            # if here, then remove this CBG!
            removed = self.remove_cbg_by_pos(pos)
            removed_cnt += 1
            ########################################################
            if verbose:
                print "noncoding CBG removed:"
                print removed
            ########################################################

        if recreate_cbgifs and removed_cnt:
            self.create_cbginterfaces()

        # return the counter of how many CBGs were removed
        return removed_cnt

    # end of function remove_noncoding_cbgs


    def remove_cexpander_allzeroscbgs(self,verbose=False,max_cbg_aa_length=20,
        omit_final_lsrcbg=True,omit_first_lsrcbg=True):
        """
        Remove CBGs from the GSG that have a negative cexpander analyses result

        @type  verbose: Boolean
        @param verbose: report information/debugging messages to STDOUT

        @type  max_cbg_aa_length: integer (or None)
        @param max_cbg_aa_length: maximal OMSR length to allow removal 

        @type  omit_final_lsrcbg: Boolean
        @param omit_final_lsrcbg: do not remove the final CBG with a lsrCBG left of it

        @type  omit_first_lsrcbg: Boolean
        @param omit_first_lsrcbg: do not remove the first CBG with a lsrCBG rigth of it

        @rtype:  integer
        @return: number of CBGs that is removed from the GSG
        """
        # loop REVERSED because elements can be removed during this for loop
        removed_cnt = 0
        for pos in range(len(self)-1,-1,-1):
            # check if position exists; if a lsrCBG is removed too
            # GSG can have shortened by >1 position
            if pos >= len(self): continue
            cbg = self.codingblockgraphs[pos]
            # no not take lsrCBGs and to long CBGs into account
            if cbg.__class__.__name__ == 'LowSimilarityRegionCodingBlockGraph': continue
            if max_cbg_aa_length and cbg.omsrlength() > max_cbg_aa_length: continue
            # get cexpander data of this CBG (if not available yet)
            if not cbg._cexpander: cbg.cexpanderanalyses()
            # if cexpanderstring is all zeros -> remove this (tiny) CBG 
            if cbg._cexpander.binarystring.count("1") == 0:
                if omit_final_lsrcbg:
                    if pos > 0 and self.codingblockgraphs[pos-1].__class__.__name__ ==\
                    'LowSimilarityRegionCodingBlockGraph':
                        ########################################################
                        if verbose:
                            print "CBG remove_cexpander_allzeroscbgs:",
                            print "omit_final_lsrcbg=True"
                            print cbg
                        ########################################################
                        # do NOT delete this finalCBG in the GSG
                        continue
                if omit_first_lsrcbg:
                    if pos == 0 and len(self) >= 2 and\
                    self.codingblockgraphs[pos+1].__class__.__name__ ==\
                    'LowSimilarityRegionCodingBlockGraph':
                        ########################################################
                        if verbose:
                            print "CBG remove_cexpander_allzeroscbgs:",
                            print "omit_first_lsrcbg=True"
                            print cbg
                        ########################################################
                        # do NOT delete this finalCBG in the GSG
                        continue

                # if here, then delete this CBG without cexpander evidence
                removed = self.remove_cbg_by_pos(pos)
                removed_cnt += 1
                ########################################################
                if verbose:
                    print "CBG removed on only zeros"
                    print removed
                ########################################################

        # return the counter of how many CBGs were removed
        return removed_cnt

    # end of function remove_cexpander_allzeroscbgs


    def remove_cexpander_lowmatchratiocbgs(self,verbose=False,
        min_match_ratio=0.21,
        max_cbg_aa_length=20):
        """
        Remove CBGs from the GSG that have a negative cexpander analyses result

        @type  verbose: Boolean
        @param verbose: report information/debugging messages to STDOUT

        @type  min_match_ratio: float 
        @apram min_match_ratio: threshold cexpander match ratio value 
        
        @type  max_cbg_aa_length: integer (or None)
        @apram max_cbg_aa_length: maximal OMSR length to allow removal
        
        @rtype:  integer
        @return: number of CBGs that is removed from the GSG
        """
        # loop REVERSED because elements can be removed during this for loop
        removed_cnt = 0 
        for pos in range(len(self)-1,-1,-1):
            cbg = self.codingblockgraphs[pos]
            # no not take lsrCBGs and to long CBGs into account
            if cbg.__class__.__name__ == 'LowSimilarityRegionCodingBlockGraph': continue
            if max_cbg_aa_length and cbg.omsrlength() > max_cbg_aa_length: continue
            # get cexpander data of this CBG (if not available yet)
            if not cbg._cexpander: cbg.cexpanderanalyses() 
            match_ratio = cbg._cexpander.uniformly_matched_ratio()
            if match_ratio == 0.0:
                # if cexpanderstring is all zeros -> ignore here.
                # allzeros is dealth with in another function
                continue
            elif match_ratio >= min_match_ratio:
                # high enough match ratio
                continue
            else:
                # remove this CBG from the GSG
                removed = self.remove_cbg_by_pos(pos)
                removed_cnt += 1
                ########################################################
                if verbose:
                    print "CBG removed on cexpander low match ratio",
                    print "%1.3f < %1.3f" % (min_match_ratio,match_ratio)
                    print removed
                ########################################################

        # return the counter of how many CBGs were removed
        return removed_cnt

    # end of function remove_cexpander_allzeroscbgs


    def remove_overlapping_cbgs(self,verbose=False,
        ignore_is_optimal_cbgif=True,
        ignore_is_compatible_cbgif=False,
        cbg_max_alowed_overlap_aa_length=CBG_MAX_ALOWED_GSGREMOVAL_OVERLAP_AA_LENGTH,
        cbg_max_alowed_overlap_ratio=CBG_MAX_ALOWED_GSGREMOVAL_OVERLAP_RATIO):
        """
        Remove overlapping CBGs in the GSG

        @type  ignore_is_optimal_cbgif: Boolean 
        @param ignore_is_optimal_cbgif: if True, leave is_optimal() cbgIFs intact 

        @type  ignore_is_compatible_cbgif: Boolean
        @param ignore_is_compatible_cbgif: if True, leave is_compatible() cbgIFs intact

        @type  cbg_max_alowed_aa_length: integer
        @param cbg_max_alowed_aa_length: maximal overlap between CBGs in AA's

        @type  cbg_max_alowed_overlap_ratio: float
        @param cbg_max_alowed_overlap_ratio: ratio between overlap and omsr (AA)

        @type  verbose: Boolean
        @param verbose: print status/debugging messages to STDOUT

        @rtype:  Integer
        @return: Number of CBGs that are removed from the GSG
        """
        removed = True
        removed_cnt = 0
        while removed:
            removed = False
            for pos in range(1,len(self)):
                # get concerned CBGs
                (cbg1,cbg2) = self.codingblockgraphs[pos-1:pos+1]
                if cbg1.__class__.__name__ == 'LowSimilarityRegionCodingBlockGraph':
                    continue
                if cbg2.__class__.__name__ == 'LowSimilarityRegionCodingBlockGraph':
                    continue
                prevCBG = None
                nextCBG = None
                if pos-2 >= 0:
                    prevCBG = self.codingblockgraphs[pos-2]
                if pos+1 < len(self):
                    nextCBG = self.codingblockgraphs[pos+1]

                # check cbgIF between these CBGs
                if cbg1._CBGinterface3p:
                    if ignore_is_optimal_cbgif and\
                    cbg1._CBGinterface3p.is_optimal():
                        continue
                    if ignore_is_compatible_cbgif and\
                    cbg1._CBGinterface3p.is_compatible():
                        print "COMPATBLE!!!"
                        continue

                # get overlap data between these 2 CBGs
                (  absPosCbg1, absPosCbg2,
                   binPosCbg1, binPosCbg2, orfIdent,
                   posRel, posBin ) = codingblock_ordering.relatively_positioned_towards(cbg1,cbg2)

                # evaluate the overlap data
                if binPosCbg1 != (1, 0, 0) or binPosCbg2 != (0, 0, 1):
                    distances = cbg1.distance_between_codingblocks(cbg2)
                    overlaps  = [ min([0,value]) for value in distances.values() ]
                    # remove all non-overlaps
                    while 0 in overlaps: overlaps.remove(0)
                    # now sum them and calculate an average
                    summed    = abs(sum(overlaps))
                    if summed:  average = float(summed) / len(overlaps)
                    else:       average = 0.0

                    ############################################################
                    #if verbose:
                    #    print cbg1
                    #    print cbg2
                    #    print pos-1, pos+1,
                    #    print absPosCbg1, absPosCbg2,
                    #    print binPosCbg1, binPosCbg2
                    #    print distances, summed, average,
                    #    print cbg1.total_weight(), cbg2.total_weight()
                    ############################################################

                    # check if the overlap ratio is not to large
                    ratioCbg1 = ( average / cbg1.omsrlength() ) > cbg_max_alowed_overlap_ratio
                    ratioCbg2 = ( average / cbg2.omsrlength() ) > cbg_max_alowed_overlap_ratio
                    remove_pos = None
                    if average > cbg_max_alowed_overlap_aa_length and ratioCbg1 and ratioCbg2:
                        # hmmm both have a high overlap ratio -> remove the lowest scoring
                        if cbg1.total_weight() < cbg2.total_weight():
                            remove_pos = pos-1
                        else:
                            remove_pos = pos
                    elif average > cbg_max_alowed_overlap_aa_length and ratioCbg1:
                        remove_pos = pos-1
                    elif average > cbg_max_alowed_overlap_aa_length and ratioCbg2:
                        remove_pos = pos
                    elif ratioCbg1 and ratioCbg2:
                        # hmmm both have a high overlap ratio -> remove the lowest scoring
                        if cbg1.total_weight() < cbg2.total_weight():
                            remove_pos = pos-1
                        else:
                            remove_pos = pos
                    elif ratioCbg1:
                        remove_pos = pos-1
                    elif ratioCbg2:
                        remove_pos = pos
                    else:
                        # omit removal!
                        pass

                    if remove_pos != None:
                        # final check: if node_count() not identical ->
                        # then remove the one that has lowest number of nodes!
                        if cbg1.node_count() > cbg2.node_count():
                            remove_pos = pos
                        elif cbg1.node_count() < cbg2.node_count():
                            remove_pos = pos-1
                        else:
                            # identical node_count -> stick to the oppointed one!
                            pass

                        # get positional pointer to CBG with which the to-be-deleted
                        # CBG is overlapping with (for verbose logginf)
                        if remove_pos == pos:
                            overlapping_cbg_pos = pos-1
                            theCBGif = cbg2._CBGinterface5p
                            # final check: does to-be-deleted CBG
                            # fill a gap in the genestructure scaffold?
                            if _is_intermediate_overlapping_cbg_a_gsg_scaffold_enrichment(
                            self,cbg1,cbg2,nextCBG):
                                # this is a perfect example of a scaffold enrichment,
                                # caused by a small exon in >= 1 genes, compared to
                                # continious exons in >= 1 other genes.
                                #######################################
                                if verbose:
                                    print "SCAFFOLD ENRICHMENT"
                                    print cbg1
                                    print cbg2
                                    print nextCBG,"NEXT"
                                #######################################
                                continue
                        else:
                            overlapping_cbg_pos = remove_pos 
                            theCBGif = cbg1._CBGinterface3p
                            # final check: does to-be-deleted CBG
                            # fill a gap in the genestructure scaffold?
                            if _is_intermediate_overlapping_cbg_a_gsg_scaffold_enrichment(
                            self,prevCBG,cbg1,cbg2):
                                # this is a perfect example of a scaffold enrichment,
                                # caused by a small exon in >= 1 genes, compared to
                                # continious exons in >= 1 other genes.
                                #######################################
                                if verbose:
                                    print "SCAFFOLD ENRICHMENT"
                                    print prevCBG,"PREV"
                                    print cbg1
                                    print cbg2
                                #######################################
                                continue

                        # remove this codingblock!
                        deletedCBG = self.remove_cbg_by_pos(remove_pos)
                        ############################################################
                        if verbose:
                            print "REMOVED!", pos, overlapping_cbg_pos, remove_pos
                            print deletedCBG
                            deletedCBG.printmultiplealignment()
                            for (key,n1,n2),pacbp in deletedCBG.pacbps.iteritems():
                                print pacbp,n1,n2
                            print "IN VIOLENCE WITH:"
                            print self.codingblockgraphs[overlapping_cbg_pos] 
                            self.codingblockgraphs[overlapping_cbg_pos].printmultiplealignment()
                            for (key,n1,n2),pacbp in\
                            self.codingblockgraphs[overlapping_cbg_pos].pacbps.iteritems():
                                print pacbp,n1,n2
                            if theCBGif: 
                                print "INTERFACE"
                                print theCBGif
                                print theCBGif._interface_is_intron
                        ############################################################
                        # set removed variable to True for the next iteration!
                        removed = True
                        removed_cnt += 1
                        break
                    else:
                        pass
                        ############################################################
                        #if verbose: print "overlap oke...."
                        ############################################################

        # return the number of removed CBGs from the GSG
        return removed_cnt

    # end of function remove_overlapping_cbgs


    def remove_unlikely_ksminxcbgs(self,
        ksminx_cbg_single_missing_unlikely_edge_bitscore=\
        KSMINX_CBG_SINGLE_MISSING_UNLIKELY_EDGE_BITSCORE,
        ksminx_cbg_multiple_missing_unlikely_edge_bitscore=\
        KSMINX_CBG_MULTIPLE_MISSING_UNLIKELY_EDGE_BITSCORE,
        verbose = False ):
        # shortcut names for parameters
        single_bitscore   = ksminx_cbg_single_missing_unlikely_edge_bitscore
        multiple_bitscore = ksminx_cbg_multiple_missing_unlikely_edge_bitscore
        # list with positions to remove from the GSG
        remove_unlikely_ksminxcbgs_positions = []

        for cbgpos in self.cbgpositerator():
            cbg = self.codingblockgraphs[cbgpos]
            if cbg.__class__.__name__ == 'LowSimilarityRegionCodingBlockGraph':
                continue
            if cbg.node_count() == self.EXACT_SG_NODE_COUNT:
                continue
            # get (estimated) bitscores of missing edges
            bitscores = estimate_bitscores_ksminxcbg(cbg,self.genetree()).values()
            single_unlikely_cnt   = [ bs > single_bitscore for bs in bitscores ]
            multiple_unlikely_cnt = [ bs > multiple_bitscore for bs in bitscores ]
            #######################################################################
            if verbose:
                print cbg, "\n", sum(single_unlikely_cnt),
                print sum(multiple_unlikely_cnt), bitscores
            #######################################################################
            # check if estimated bitscores are small enough
            if sum(single_unlikely_cnt) == 0 and sum(multiple_unlikely_cnt) < 2: continue

            # if here, this ksminxCBG is likely to be deleted (if no N-tracks in DNA sequences)
            diforgnodes = self.genetree().different_organisms(cbg)
            dnaseqoffsets = {}
            ksminx_is_removable = True
            for org in diforgnodes:
                dnaseqoffsets[org] = { 'sta': 0, 'end': len(self.input[org]['genomeseq']) }
                for stapos in range(cbgpos-1,-1,-1):
                    compareCBG = self.codingblockgraphs[stapos]
                    if org in compareCBG.organism_set():
                        _omsr = compareCBG.overall_minimal_spanning_range(organism=org)
                        dnaseqoffsets[org]['sta'] = max(_omsr) * 3
                        break
                for endpos in range(cbgpos+1,len(self)):
                    compareCBG = self.codingblockgraphs[endpos]
                    if org in compareCBG.organism_set():
                        _omsr = compareCBG.overall_minimal_spanning_range(organism=org)
                        dnaseqoffsets[org]['end'] = min(_omsr) * 3
                        break
                # check if this organisms' DNA sequences has N's in this range
                if self.input[org]['genomeseq'][ dnaseqoffsets[org]['sta']: dnaseqoffsets[org]['end'] ].upper().count("N") > 0:
                    ksminx_is_removable = False
                    break
            #######################################################################
            if verbose: print "IS_REMOVABLE:", ksminx_is_removable
            #######################################################################
            if ksminx_is_removable:
                remove_unlikely_ksminxcbgs_positions.append( cbgpos )


        # remove the remove_unlikely_ksminxcbgs_positions
        remove_unlikely_ksminxcbgs_positions.reverse()
        for cbgpos in remove_unlikely_ksminxcbgs_positions:
            self.codingblockgraphs.pop( cbgpos )

        # call remove_incorrect_lsrcbgs(); should not make a difference,
        # but just to be certain
        self.remove_incorrect_lsrcbgs()

        # return number of ksminxCBGs that are removed
        return len(remove_unlikely_ksminxcbgs_positions)

    # end of function remove_unlikely_ksminxcbgs


    def replace_scaffold_breaking_cbgs(self,verbose=False):
        """
        (Try) to replace CBG that break the GSG scaffold by other CBGs

        @type  verbose: Boolean
        @param verbose: print status/debugging messages to STDOUT

        @rtype:  Boolean
        @return: Is any CBG replaced?
        """
        # Boolean return value
        scaffold_breaking_cbg_replaced = False

        for cbgpos in self.cbgpositerator(reversed=True)[1:]:
            cbg = self.codingblockgraphs[cbgpos]
            if cbg.__class__.__name__ == 'LowSimilarityRegionCodingBlockGraph':
                continue
            if cbg.node_count() < self.EXACT_SG_NODE_COUNT:
                continue

            # loop forwards through the GSG and look for mutual nodes
            identical_nodes = []
            for backwardspos in range(cbgpos+1,len(self)):
                comparecbg = self.codingblockgraphs[backwardspos]
                if cbg.mutual_nodes(comparecbg):
                    identical_nodes.append(True)
                    break
                else:
                    identical_nodes.append(False)
            if not identical_nodes:
                continue    # final CBG -> continue
            elif identical_nodes == [True]:
                continue    # neighboring node has mutual nodes -> continue
            elif identical_nodes.count(True) == 0:
                continue    # no mutual nodes at all -> continue
            else:
                # this is what we are looking for, a list like
                # [ False, ... True ] with >1 False
                # get total_weights of the intermediate CBGs
                tws = [ self.codingblockgraphs[_pos].total_weight() for\
                        _pos in range(cbgpos+1,backwardspos) ]
                ################################################################
                if verbose:
                    print cbgpos, backwardspos, identical_nodes
                    print cbg
                    print tws
                    print comparecbg
                ################################################################

                # first, check if the CBGs are already partially overlapping
                # if so, get rid of this intermediate CBG
                omsrdist = cbg.omsr_distance_between_codingblocks(comparecbg)
                if max(omsrdist.values()) <= 1:
                    # yes, all organisms glue these CBGs perfectly together
                    # just remove this one without further checks
                    cbg._CBGinterface3p        = None
                    comparecbg._CBGinterface5p = None
                    self.codingblockgraphs.__setslice__(
                            cbgpos+1,
                            backwardspos+1,
                            [ comparecbg ] )
                    scaffold_breaking_cbg_replaced = True
                    # go to the next cbg in the list
                    continue

                # do a more eleborate check by trying to create a CBG
                # in this large_scaffold_gap
                from graph_genestructure import GenestructureOfCodingBlockGraphs
                partialGSG = GenestructureOfCodingBlockGraphs(self.input)
                partialGSG.codingblockgraphs = [ cbg, comparecbg ]
                partialGSG._GENETREE = self._GENETREE
                partialGSG.create_large_intermediate_cbg_for_scaffold_gap(
                        sprdif_min_node_count = 2,
                        cbg_min_node_count = self.EXACT_SG_NODE_COUNT,
                        verbose = verbose
                        )
                if len(partialGSG) == 2:
                    ############################################################
                    if verbose: print "NO scaffold CBGs found!"
                    ############################################################
                    pass
                else:
                    new_tws = [ _cbg.total_weight() for _cbg in\
                            partialGSG.codingblockgraphs[1:-1] ]
                    if sum(new_tws) > sum(tws):
                        # replace!
                        self.codingblockgraphs.__setslice__(
                            cbgpos+1,
                            backwardspos,
                            partialGSG.codingblockgraphs[1:-1] )
                        scaffold_breaking_cbg_replaced = True
                    else:
                        pass
                    ############################################################
                    if verbose:
                        if sum(new_tws) > sum(tws):
                            print "REPLACING scaffold-breaking CBGs!!!"
                        else:
                            print "MAINTAINING scaffold-breaking CBGs!!!"
                        for cbg in partialGSG: print cbg
                    ############################################################

        # return if CBGs are removed
        return scaffold_breaking_cbg_replaced

    # end of function replace_scaffold_breaking_cbgs


    def remove_scaffold_breaking_cbgs(self,verbose=False):
        """
        Remove any CBG that breaks the GSG scaffold

        @type  verbose: Boolean
        @param verbose: print status/debugging messages to STDOUT

        @attention: USE WITH CARE!!!

        @rtype:  Integer
        @return: Number of CBGs that have been removed from the GSG
        """
        scaffold_breaking_cbg_removed_cnt = 0

        for cbgpos in self.cbgpositerator(reversed=True)[1:]:
            cbg = self.codingblockgraphs[cbgpos]
            if cbg.__class__.__name__ == 'LowSimilarityRegionCodingBlockGraph':
                continue
            if cbg.node_count() < self.EXACT_SG_NODE_COUNT:
                continue

            # loop forwards through the GSG and look for mutual nodes
            identical_nodes = []
            for backwardspos in range(cbgpos+1,len(self)):
                comparecbg = self.codingblockgraphs[backwardspos]
                if cbg.mutual_nodes(comparecbg):
                    identical_nodes.append(True)
                    break
                else:
                    identical_nodes.append(False)
            if not identical_nodes:
                continue    # final CBG -> continue
            elif identical_nodes == [True]:
                continue    # neighboring node has mutual nodes -> continue
            elif identical_nodes.count(True) == 0:
                continue    # no mutual nodes at all -> continue
            else:
                # this is what we are looking for, a list like
                # [ False, ... True ] with >1 False
                # remove these CBGs
                for _tobepopped in range(backwardspos-1,cbgpos,-1):
                    self.codingblockgraphs.pop( _tobepopped )
                    scaffold_breaking_cbg_removed_cnt += 1

        # if cbgs are removed -> call the lsr correction function
        if scaffold_breaking_cbg_removed_cnt:
            self.remove_incorrect_lsrcbgs()

        # return if CBGs are removed
        return scaffold_breaking_cbg_removed_cnt

    # end of function remove_scaffold_breaking_cbgs

# end of class CodingBlockGraphRemovalFunctions


def _is_intermediate_overlapping_cbg_a_gsg_scaffold_enrichment(gsg,
    cbgA,cbgB,cbgC,minimal_scaffold_aa_enrichment=5):
    """
    """
    if not (gsg and cbgA and cbgB and cbgC):
        # (most likely) cbgA or cbgC is not defined
        # function behavious should be ti return False
        return False
    if not cbgA.mutual_nodes(cbgC):
        # series of CBGs does not represent an suitable
        # gene structure scaffold -> return False
        return False

    # perform this check
    from graph_genestructure import GenestructureOfCodingBlockGraphs
    partGSG = GenestructureOfCodingBlockGraphs(gsg.input)
    partGSG.codingblockgraphs = [ cbgA,cbgB,cbgC ]
    partGSG._GENETREE = gsg._GENETREE
    partOMSRa = partGSG.overall_minimal_spanning_range()
    partGSG.codingblockgraphs = [ cbgA,cbgC ]
    partOMSRb = partGSG.overall_minimal_spanning_range()
    scaffold_enrichments = []
    for node in cbgA.mutual_nodes(cbgC):
        org = gsg.organism_by_node(node)
        scaffold_enrichments.append(
            len(partOMSRa[org]) - len(partOMSRb[org]) >=\
            minimal_scaffold_aa_enrichment
            )
    # check if True in scaffold_enrichments
    if True in scaffold_enrichments:
        return True
    else:
        return False

# end of function _is_intermediate_overlapping_cbg_a_gsg_scaffold_enrichment

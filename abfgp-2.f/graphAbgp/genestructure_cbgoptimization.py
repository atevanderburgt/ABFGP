"""
Functions for improving CBGs in the GSG
"""

# graphAbgp Imports
from exceptions import NoOverallMinimalSpanningRange

# Abgp Imports
import lib_cexpander

# Global variable Imports

# Python Imports



class CbgInGeneStructureOptimizationFunctions:
    """  """
    def optimization_cexpandercheckCBGs4omsrbordergaps(self,**kwargs):
        """
        Improve (shorten) the CBGs when non-uniformly aligned areas are present at the OMSR's exterior

        @attention: see lib_cexpander.checkCBGs4omsrbordergaps for documentation
        """
        return lib_cexpander.checkCBGs4omsrbordergaps(self,**kwargs)

    # end of function optimization_cexpandercheckCBGs4omsrbordergaps


    def optimization_improvecbgalignment(self,verbose=False,
        ignore_ksminx_cbgs=False,ignore_ks_cbgs=False):
        """
        Improve the PacbPORF alignments of the CBGs in the GSG

        @type  ignore_ks_cbgs: Boolean
        @param ignore_ks_cbgs: if True, K(s) CBGs are skipped

        @type  ignore_ksminx_cbgs: Boolean
        @param ignore_ksminx_cbgs: if True, K(s-x) CBGs are skipped

        @type  verbose: Boolean
        @param verbose: print status/debugging messages to STDOUT

        @rtype:  list
        @return: list with Booleans for each of the CBGs in the GSG
        """
        statusses = []
        for pos in range(0,len(self)):
            cbg = self.codingblockgraphs[pos]
            if ignore_ks_cbgs and cbg.node_count() == self.EXACT_SG_NODE_COUNT: 
                continue
            if ignore_ksminx_cbgs and cbg.node_count() < self.EXACT_SG_NODE_COUNT:
                continue
            # check which side(s) of the CBG can be optimized
            if pos==0:
                allow_5p_optimization = True
            elif self.codingblockgraphs[pos-1].__class__.__name__ ==\
            'LowSimilarityRegionCodingBlockGraph':
                allow_5p_optimization = False
            else:
                allow_5p_optimization = True
            if pos==len(self)-1:
                allow_3p_optimization = True
            elif self.codingblockgraphs[pos+1].__class__.__name__ ==\
            'LowSimilarityRegionCodingBlockGraph':
                allow_3p_optimization = False
            else:
                allow_3p_optimization = True

            # call the cbg.improvealignment() function
            status = cbg.improvealignment(verbose=verbose,
                allow_5p_optimization=allow_5p_optimization,
                allow_3p_optimization=allow_3p_optimization )

            # append the status to the return list
            statusses.append(status)

            # WARNING: this function is NOT called on CBGs
            # that neighbour on a lsrCBG. When this is changed
            # (maybe in the future), do not forget to call the
            # self.correct_lsrcbgs_after_optimization() function 

        # return list with Boolean statusses
        return statusses

    # end of function  optimization_improvecbgalignment


    def optimization_correct_pacbpgaps_nearby_omsr(self,omsr_offset=3,gap_size=3,
        ignore_ksminx_cbgs=False, ignore_ks_cbgs=False,
        ignore_lsrcbg_boundaries=True, verbose=False):
        """
        Correct PacbPORFs in the GSG for gaps directly around the OMSR

        @type  omsr_offset: integer 
        @param omsr_offset: ...

        @type  gap_size: integer
        @param gap_size: ...

        @type  ignore_lsrcbg_boundaries: Boolean
        @param ignore_lsrcbg_boundaries: what to do with a CBG that is neighbored by an lsrCBG?

        @type  ignore_ks_cbgs: Boolean
        @param ignore_ks_cbgs: if True, K(s) CBGs are skipped

        @type  ignore_ksminx_cbgs: Boolean
        @param ignore_ksminx_cbgs: if True, K(s-x) CBGs are skipped

        @type  verbose: Boolean
        @param verbose: print status/debugging messages to STDOUT

        @attention: see codingblock_optimization.correct_pacbpgaps_nearby_omsr for documentation

        @rtype:  list
        @return: list with NoneBooleans for each of the CBGs in the GSG
        """
        statusses = []
        # loop REVERSED because elements might be added/removed during this for loop
        for pos in range(len(self)-1,-1,-1):
            cbg = self.codingblockgraphs[pos]
            if cbg.IS_IGNORED:
                continue
            if cbg.__class__.__name__ == 'LowSimilarityRegionCodingBlockGraph':
                continue
            if ignore_ks_cbgs and cbg.node_count() == self.EXACT_SG_NODE_COUNT: 
                continue
            if ignore_ksminx_cbgs and cbg.node_count() < self.EXACT_SG_NODE_COUNT:
                continue

            # vars stating if or not a side must be omitted from optimization
            omit5pside, omit3pside = False, False

            # check the 5' side / previous CBG if it is a lsrCBG
            if ignore_lsrcbg_boundaries and pos > 0:
                prev = self.codingblockgraphs[pos-1]
                if prev.__class__.__name__ == 'LowSimilarityRegionCodingBlockGraph':
                    omit5pside = True

            # check the 3' side / next CBG if it is a lsrCBG
            if ignore_lsrcbg_boundaries and pos < len(self)-1:
                next = self.codingblockgraphs[pos+1]
                if next.__class__.__name__ == 'LowSimilarityRegionCodingBlockGraph':
                    omit3pside = True

            try:
                status = cbg.correct_pacbpgaps_nearby_omsr(
                        omit5pside = omit5pside, omit3pside = omit3pside,
                        omsr_offset=omsr_offset,gap_size=gap_size,verbose=verbose)
                statusses.append( status )
                if status:
                    # CBG is changed, so it's GeneTreeGraph is likely changed too
                    gtg = cbg.set_genetree()
            except NoOverallMinimalSpanningRange:
                # mark this cbg as to-be-ignored!
                cbg.IS_IGNORED = True
                statusses.append( None )
            except:
                # hmmmm.... an very unexpected error. This should be raised,
                # resulting in hard-termination of the program!
                status = cbg.correct_pacbpgaps_nearby_omsr(verbose=True)
                raise "Error occurred previously (no OMSR)"

        # Check if there are CBGs that are IS_IGNORED
        # If so, remove these from the GSG by the microsyntheny function
        if None in statusses:
            # there is at least 1 CBG marked for deletion (set to IS_IGNORED)
            # here, we have to check if one of is neighbours happens to be an
            # lsrCBG. If so, mark this one for deletion too!
            for pos in range(0,len(statusses)):
                if statusses[pos] == None:
                    if pos > 0:
                        prevcbgclass = self.codingblockgraphs[pos-1].__class__.__name__
                        if prevcbgclass == 'LowSimilarityRegionCodingBlockGraph':
                            self.codingblockgraphs[pos-1].IS_IGNORED = True
                    if pos < len(statusses)-1:
                        nextcbgclass = self.codingblockgraphs[pos+1].__class__.__name__
                        if nextcbgclass == 'LowSimilarityRegionCodingBlockGraph':
                            self.codingblockgraphs[pos+1].IS_IGNORED = True
            # remove the IS_IGNORED cbgs from the GSG
            dsGSG, usGSG, etcGSG = self.separate_ds_and_us_gsg()
            etcGSG.codingblockgraphs.extend(dsGSG) # when IS_IGNORED at the start of the GSG list
            etcGSG.codingblockgraphs.extend(usGSG) # when IS_IGNORED at the end   of the GSG list
        else:
            etcGSG = []

        # Check if there are CBGs changed by correct_pacbpgaps_nearby_omsr
        # If so, recreate lsrCBGs because OMSR coords might be changed
        # And, re-create the GeneTreeGraph of the GSG, because some of
        # the CBGs are changed, potentially causing a change in *the* GTG
        if True in statusses:
            LSRCBG_CORRECTED = self.correct_lsrcbgs_after_optimization()
            thegtg = self.set_genetree()
        else:
            LSRCBG_CORRECTED = 0

        ################################################################
        if verbose:
            print "SUMMARY: correct_pacbpgaps_nearby_omsr"
            print "CBGS removed beause no OMSR:", len(etcGSG)
            print "CBGs succesfully optimized:", statusses.count(True)
            print "lsrCBGs corrected after CBG optimization", LSRCBG_CORRECTED
        ################################################################

        # return list of statusses
        return statusses

    # end of function optimization_correct_pacbpgaps_nearby_omsr

# end of class CbgInGeneStructureOptimizationFunctions

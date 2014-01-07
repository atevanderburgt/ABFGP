"""
Alignment Based Gene Predictions settings: settings for AligedPSSMSiteGraphs for
selection of optimal splice sites, TSS, StopCodons, StartCodons etc.
"""

################################################################################
# Thresholds for length range to allow a search for aligned splicesites and TSS
# Recommended to set spatiously large!!
################################################################################
ALIGNED_ACCEPTOR_MAX_TRIPLET_DISTANCE               = 4
ALIGNED_DONOR_MAX_TRIPLET_DISTANCE                  = 4
ALIGNED_STOP_MAX_AA_DISTANCE                        = 50    # NOT IMPLEMENTED YET!!!!
ALIGNED_TSS_MAX_AA_DISTANCE                         = 20    # How far can TSS-Methionines be appart
                                                            # in order to be taken into account as
                                                            # alignable?
ALIGNED_ATG_MAX_AA_DISTANCE                         = 20    # How far can ATG-Methionines be appart
                                                            # in order to be taken into account as
                                                            # alignable?


ELIGABLE_ALIGNED_START_SITES_AA_OFFSET = ALIGNED_ATG_MAX_AA_DISTANCE    # TO BE DEPRECATED!!



ELIGABLE_DONOR_SITE_MINIMAL_AA_OFFSET           = 5
ELIGABLE_ACCEPTOR_SITE_MINIMAL_AA_OFFSET        = 5


ELIGABLE_DONOR_SITE_LEFT_OF_OMSR_AA_OFFSET      = 20    # Threshold for how far LEFT of the max(OMSR)
                                                        # donor sites may be harvested. So, this position
                                                        # is located IN the OMSR of the CBG.
                                                        # In the function minimal_eligable_donor_site_position
                                                        # in graphAbgp.codingblock_collectionharvesting
                                                        # this number is further restricted by the identity of
                                                        # the CBG and possibly limited by min(OMSR)
ELIGABLE_ACCEPTOR_SITE_RIGTH_OF_OMSR_AA_OFFSET  = 20    # Threshold for how far RIGTH of the min(OMSR)
                                                        # acceptor sites may be harvested. So, this position
                                                        # is located IN the OMSR of the CBG.
                                                        # In the function maximal_eligable_acceptor_site_position
                                                        # in graphAbgp.codingblock_collectionharvesting
                                                        # this number is further restricted by the identity of 
                                                        # the CBG and possibly limited by max(OMSR)

# ------------------------------------------------------------------------------------------------
# ELIGABLE_DONOR_SITE_LEFT_OF_OMSR_AA_OFFSET   ID%   ELIGABLE_DONOR_SITE_RIGTH_OF_OMSR_AA_OFFSET
# ------------------------------------------------------------------------------------------------
#                                         ..   100   ............
#                                        ...         ..........
#                                       ....    75   .......
#                                      .....         .....
#                                    .......    50   ....
#                                 ..........         ...
#                               ............    25   ..
# ------------------------------------------------------------------------------------------------


ELIGABLE_DONOR_SITE_RIGTH_OF_OMSR_AA_OFFSET         = 100   # threshold for the current
                                                            # PacbpOrf/graph
ELIGABLE_ACCEPTOR_SITE_LEFT_OF_OMSR_AA_OFFSET       = 100    # threshold for the current
                                                            # PacbpOrf/graph


################################################################################
# Maximal offset (in AA) to allow for a projected splice site.
# Used when organism A has an intron and B lacks this intron.
# When the aligned protein sequence is 100% continious, donor and acceptor
# sites can be projected on exactly the same sequence coordinate.
# However, a small offset is likely to arise over evolutionary time.
# This threshold defines the maximal distance to occur.
# The best projection possible has an offset of zero (offset is absolute value)
# The H0-hypothesis (in fungal genes) is conservation of protein sequence
# and therefor as little offset as possible.
################################################################################
MAX_PROJECTED_INTRON_AA_OFFSET  = 5


################################################################################
# Thresholds for alignment of TSS 
################################################################################
ELIGABLE_ALIGNED_TSS_3P_AA_OFFSET = 50
ELIGABLE_ALIGNED_TSS_5P_AA_OFFSET = None # (no limit)

################################################################################
# Thresholds for proposing a higher PSSM scoring SpliceSite in favour of
# the SpliceSite in the optimal AlignedSpliceSiteGraph in a cbgIF
################################################################################
PROPOSE_IMPROVED_DONOR_MAX_AA_DIST      = 4  # in older code 6 was chosen
PROPOSE_IMPROVED_ACCEPTOR_MAX_AA_DIST   = 4  # in older code 6 was chosen
PROPOSE_IMPROVED_DONOR_MIN_PSSM_RATIO   = 1.5
PROPOSE_IMPROVED_ACCEPTOR_MIN_PSSM_RATIO= 1.5


################################################################################
# Thresholds for creating the DSCG & ASCG in a cbgIF in the optimize() function
################################################################################
OPTIMIZE_CBGIF_DONOR_MAX_TRIPLET_DISTANCE       = 8
OPTIMIZE_CBGIF_DONOR_ENLARGE_3P_BOUNDARY_BY     = 0 # TODO!!
OPTIMIZE_CBGIF_ACCEPTOR_MAX_TRIPLET_DISTANCE    = 8    
OPTIMIZE_CBGIF_ACCEPTOR_ENLARGE_5P_BOUNDARY_BY  = 15

################################################################################
# Thresholds used in the optimizetinyexon() function in the cbgIF
################################################################################
OPTIMIZE_CBGIF_TINYEXON_MAX_AA_SIZE                = 10
OPTIMIZE_CBGIF_TINYEXON_DONOR_OMSR_3P_NT_OFFSET    = 15
OPTIMIZE_CBGIF_TINYEXON_ACCEPTOR_OMSR_5P_NT_OFFSET = 15


################################################################################
# Thresholds when allowing splice site phase shift in Aligned Site Graphs 
################################################################################
MAX_SPLICE_SITE_PHASE_SHIFT_NT_DISTANCE = 5
MIN_DONOR_SITE_PHASE_SHIFT_PSSM_SCORE   = 0.0     # example NCU02366 tiny exon donor
                                                  # ('ncu', 61, 2852)
                                                  # <SpliceDonor 2849-2858 (2), score=0.927, catGTacag>
MIN_ACCEP_SITE_PHASE_SHIFT_PSSM_SCORE   = 0.0


################################################################################
# Threshold values aligned donor site graphs
################################################################################
ALIGNED_DONOR_MIN_TOTAL_WEIGTH              = 0.0       # NOT IMPLEMENTED (YET)
ALIGNED_DONOR_MIN_TOTAL_PSSM                = 0.0       # NOT IMPLEMENTED (YET)
ALIGNED_DONOR_MIN_BINARY_ENTROPY            = 0.0       # NOT IMPLEMENTED (YET)
ALIGNED_DONOR_MIN_CUMULATIVE_SCORE          = 0.0       # NOT IMPLEMENTED (YET)


################################################################################
# Threshold values aligned acceptor site graphs
################################################################################
ALIGNED_ACCEPTOR_MIN_TOTAL_WEIGTH          = 0.0       # NOT IMPLEMENTED (YET)
ALIGNED_ACCEPTOR_MIN_TOTAL_PSSM            = 0.0       # NOT IMPLEMENTED (YET)
ALIGNED_ACCEPTOR_MIN_BINARY_ENTROPY        = 0.0       # NOT IMPLEMENTED (YET)
ALIGNED_ACCEPTOR_MIN_CUMULATIVE_SCORE      = 0.0       # NOT IMPLEMENTED (YET)


################################################################################
# Threshold values aligned TranslationalStartSite graphs
################################################################################
ALIGNED_TSS_MIN_TOTAL_WEIGTH                = 0.0       # NOT IMPLEMENTED (YET)
ALIGNED_TSS_MIN_TOTAL_PSSM                  = 0.0       # NOT IMPLEMENTED (YET)
ALIGNED_TSS_MIN_BINARY_ENTROPY              = 0.0       # NOT IMPLEMENTED (YET)
ALIGNED_TSS_MIN_CUMULATIVE_SCORE            = 0.0       # NOT IMPLEMENTED (YET)


################################################################################
# Threshold values aligned ATG-Methionine graphs
################################################################################
ALIGNED_ATG_MIN_TOTAL_WEIGTH                = 0.0       # NOT IMPLEMENTED (YET)
ALIGNED_ATG_MIN_TOTAL_PSSM                  = None      # DOES NOT EXIST!!!
ALIGNED_ATG_MIN_BINARY_ENTROPY              = 0.0       # NOT IMPLEMENTED (YET)
ALIGNED_ATG_MIN_CUMULATIVE_SCORE            = 0.0       # NOT IMPLEMENTED (YET)


################################################################################
# Threshold values aligned STOP-codon graphs
################################################################################
ALIGNED_STOP_MIN_TOTAL_WEIGTH               = 0.0       # NOT IMPLEMENTED (YET)
ALIGNED_STOP_MIN_TOTAL_PSSM                 = None      # DOES NOT EXIST!!!
ALIGNED_STOP_MIN_BINARY_ENTROPY             = 0.0       # NOT IMPLEMENTED (YET)
ALIGNED_STOP_MIN_CUMULATIVE_SCORE           = 0.0       # NOT IMPLEMENTED (YET)

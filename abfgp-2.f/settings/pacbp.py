"""
Alignment Based Gene Predictions settings: settings for PacbP acceptance/rejection
"""

################################################################################
#
################################################################################
MINIMUM_OVERLAPPING_ORF_NT_LENGTH   = 21 
MINIMUM_PACB_ALIGNMENT_AA_LENGTH    = MINIMUM_OVERLAPPING_ORF_NT_LENGTH/3
GTGADJUSTED_PACBPS_SPLIT_ON_GAPSIZE = None # calculated from PACBPS_SPLIT_ON_GAPSIZE and GTG
PACBPS_SPLIT_ON_GAPSIZE             = 8    # This number is corrected by the identity
                                           # of the GTG: PACBPS_SPLIT_ON_GAPSIZE / id%
                                           # Example values of actual used gapsize (8):
                                           # id% 100     8
                                           # id%  90     9
                                           # id%  80    10 
                                           # id%  70    11
                                           # id%  60    13
                                           # id%  50    16


################################################################################
# Linearization acceptance margin thresholds
# Recommended to set spatiously large and small!
################################################################################
LINEARIZATION_ACCEPTANCE_MARGIN     = 1000      # in AA coordinates; so 3kb (was 1kb...)
                                                # when much smaller, this can cause
                                                # bona fide pacbps to be rejected in the
                                                # extreme case of a gene with many and a
                                                # gene with little exons
LINEARIZATION_STARTWITHBEST         = 1         # take the n best pacbs as a start
                                                # best performance is garantied by 1
                                                # >=2 incorporates the risk of bridging
                                                # the intergenic distance in case of
                                                # a region of microsyntheny
LINEARIZATION_WEIGTHED              = True      # weigth the linearization by the
                                                # alignment scores of the pacb's


################################################################################
PACBPORF_HIGH_GAP_RATIO_THRESHOLD = 0.40


################################################################################
# Some thresholds for discarding overlapping PACB's
# TODO: are all still in use??
################################################################################
WEAKER_SIMILARITY_OTHER_FRAME_RATIO = 2.0
ALTERNATIVE_ALIGNMENT_OVERLAP_RATIO = 0.8
REPETITIVE_ALIGNMENT_OVERLAP_RATIO  = 0.8
ALWAYS_SPLIT_HSPS                   = False # WAS ORIGINALLY TRUE ....
SPLIT_INDENTITY_CUTOFF              = 0.90
SPLIT_ON_GAP_SIZE                   = 2



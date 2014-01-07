"""
Alignment Based Gene Predictions settings: settings for (PSSM) Translational
Start Sites (TSS)
"""
from settings.abgp import MAIN_ABGP_PATH
from os.path import join as osPathJoin


################################################################################
# PSSM_IC file for Translational Start Sites
################################################################################
IC_TSS_DATA_FILE       = osPathJoin(MAIN_ABGP_PATH,"datafiles/ic_start_5fungi.txt")
IC_TSS_PATTERN_OFFSET  = (9,7)

################################################################################
# Threshold values for TranslationalStartSite (ATG) PSSM cutoffs
# Recommended to set spatiously low!
# Non-cannonical splice sites are not supported yet
################################################################################
TSS_MIN_PSSM_SCORE                  = float(-1)
TSS_ALLOW_NON_CANONICAL             = False
TSS_NON_CANONICAL_MIN_PSSM_SCORE    = float(0)


################################################################################
# Threshold for defining if a TSS is optimal in a given range
################################################################################
TCODE_TSS_5P_WINDOW                 = 201     # nt coords
TCODE_TSS_3P_WINDOW                 = 201     # nt coords
TCODE_AVERAGE_SCORE                 = 0.845


TSS_IS_OPTIMAL_5P_WINDOW            = 200     # nt coords, was 250 in the juli2009 abgp_tssanalyses test
TSS_IS_OPTIMAL_3P_WINDOW            = 100     # nt coords
TSS_IS_OPTIMAL_MIN_PSSM_SCORE       = 2.0

################################################################################
# Threshold values aligned TSS graphs
################################################################################
ALIGNED_TSS_MIN_TOTAL_WEIGTH        = 0.0       # NOT IMPLEMENTED (YET)
ALIGNED_TSS_MIN_TOTAL_PSSM          = 0.0       # NOT IMPLEMENTED (YET)
ALIGNED_TSS_MIN_BINARY_ENTROPY      = 0.0       # NOT IMPLEMENTED (YET)
ALIGNED_TSS_MIN_CUMULATIVE_SCORE    = 0.0       # NOT IMPLEMENTED (YET)

"""
Alignment Based Gene Predictions settings: settings for iterative blastp
"""

from settings.abgp import MAIN_ABGP_PATH
from settings import ParameterOptions
from settings.emboss import (
    BLOSUM62_PATH,
    BLOSUM80_PATH,
    BLOSUM45_PATH,
    PAM30_PATH,
    PAM70_PATH,  
    )

class CrossBlastOptions(ParameterOptions):
    """ Empty class to assemble blast-related settings """
    def __init__(self): pass
# end of class CrossBlastOptions

################################################################################
# Thresholds defining Pacb characteristics (~= aligned protein sequences) 
# I.e. limitations for blastp hits.
# Take care that when choosing a blastp matrix (BLASTP_MATRIX_NAME),
# some BLASTP_EXTRA_PARAMS (gap opening and extension!) are required to be altered as well.
# Below, valid combinations are listed (matrix (M), opening (G), extension (E):
#   BLOSUM62     6-11        2              BLOSUM45    10-13        3
#   BLOSUM62     9-13        1              BLOSUM45    12-16        2
#   BLOSUM80     6-9,13,25   2              BLOSUM45    16-19        1
#   BLOSUM80     9-11        1
# Recommended combinations for usage in ABGP are:
#   BLOSUM45, 10, 3
#   BLOSUM62,  7, 2
#   BLOSUM80,  7, 2
################################################################################

BLASTP_MATRIX_NAME                = "BLOSUM62"  # default matrix name, mostly overwritten.
BLASTP_DIRECTLY_IGNORE_TINY_HITS  = 4
BLASTP_LONGEST_ORFS               = 20


# blastoptions for the first iteration
blastoptions_iter1 = CrossBlastOptions()
blastoptions_iter1.BLASTP_DIRECTLY_IGNORE_TINY_HITS = BLASTP_DIRECTLY_IGNORE_TINY_HITS
blastoptions_iter1.BLASTP_LONGEST_ORFS              = BLASTP_LONGEST_ORFS 
blastoptions_iter1.BLASTP_MATRIX_NAME               = "BLOSUM62"
blastoptions_iter1.BLASTP_HSP_MINIMAL_LENGTH        = 12
blastoptions_iter1.BLASTP_HSP_MINIMAL_BITS          = 60
blastoptions_iter1.BLASTP_HSP_MAXIMAL_EXPECT        = 0.1   # 0.1 means little filtering
blastoptions_iter1.BLASTP_EXTRA_PARAM_F             = 'F'   # F   blastp low-complexity filtering
blastoptions_iter1.BLASTP_EXTRA_PARAM_G             = 11    # 11  blastp gap opening cost
blastoptions_iter1.BLASTP_EXTRA_PARAM_E             = 1     # 1   blastp gap extension cost
blastoptions_iter1.BLASTP_EXTRA_PARAM_W             = 3     # 3   !!blastp word size!! 
blastoptions_iter1.BLASTP_EXTRA_PARAM_e             = blastoptions_iter1.BLASTP_HSP_MAXIMAL_EXPECT
blastoptions_iter1.extra_blastp_params = {
                        'F': blastoptions_iter1.BLASTP_EXTRA_PARAM_F,
                        'e': blastoptions_iter1.BLASTP_EXTRA_PARAM_e,
                        'G': blastoptions_iter1.BLASTP_EXTRA_PARAM_G,
                        'E': blastoptions_iter1.BLASTP_EXTRA_PARAM_E,
                        'W': blastoptions_iter1.BLASTP_EXTRA_PARAM_W,
                        'M': blastoptions_iter1.BLASTP_MATRIX_NAME,
                      }


# blastoptions for the second iteration
blastoptions_iter2 = CrossBlastOptions()
blastoptions_iter2.BLASTP_MATRIX_NAME               = "BLOSUM45"
blastoptions_iter2.BLASTP_DIRECTLY_IGNORE_TINY_HITS = BLASTP_DIRECTLY_IGNORE_TINY_HITS
blastoptions_iter2.BLASTP_LONGEST_ORFS              = BLASTP_LONGEST_ORFS 
blastoptions_iter2.BLASTP_MATRIX_NAME               = "BLOSUM45"
blastoptions_iter2.BLASTP_HSP_MINIMAL_LENGTH        = 6 
blastoptions_iter2.BLASTP_HSP_MINIMAL_BITS          = 40
blastoptions_iter2.BLASTP_HSP_MAXIMAL_EXPECT        = 50.0  # 50  means little or no filtering at all 
blastoptions_iter2.BLASTP_EXTRA_PARAM_F             = 'F'   # F   blastp low-complexity filtering
blastoptions_iter2.BLASTP_EXTRA_PARAM_G             = 10    # 10  blastp gap opening cost
blastoptions_iter2.BLASTP_EXTRA_PARAM_E             = 3     # 3   blastp gap extension cost
blastoptions_iter2.BLASTP_EXTRA_PARAM_W             = 2     # 2   !!blastp word size!!
blastoptions_iter2.BLASTP_EXTRA_PARAM_e             = blastoptions_iter2.BLASTP_HSP_MAXIMAL_EXPECT
blastoptions_iter2.extra_blastp_params = {
                        'F': blastoptions_iter2.BLASTP_EXTRA_PARAM_F,
                        'e': blastoptions_iter2.BLASTP_EXTRA_PARAM_e,
                        'G': blastoptions_iter2.BLASTP_EXTRA_PARAM_G,
                        'E': blastoptions_iter2.BLASTP_EXTRA_PARAM_E,
                        'W': blastoptions_iter2.BLASTP_EXTRA_PARAM_W,
                        'M': blastoptions_iter2.BLASTP_MATRIX_NAME,
                      }


# list with the iterative blast settings
ITERATIVE_BLAST_SETTINGS = [ blastoptions_iter1, blastoptions_iter2 ] 


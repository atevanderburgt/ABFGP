"""
Alignment Based Gene Predictions settings: settings for (PSSM) splice sites
"""
from settings.abgp import MAIN_ABGP_PATH
from settings import ParameterOptions
# defaults for splice site thresholds
from settings.genestructure import (
    MIN_INTRON_NT_LENGTH,
    MAX_INTRON_NT_LENGTH,
    TINYEXON_MAX_NT_LENGTH,
    TINYEXON_MIN_NT_LENGTH,
    TINYEXON_MAX_INTRON_NT_LENGTH,
    TINYEXON_MIN_INTRON_NT_LENGTH,
    OPTIMAL_BRACNHPOINT_TO_ACCEPTOR_DISTANCE,
    )

from os.path import join as osPathJoin

class DonorPSSMOptions(ParameterOptions):
    """ Empty class to assemble DonorSite PSSM settings """
    def __init__(self): pass
# end of class DonorPSSMOptions 

class AcceptorPSSMOptions(ParameterOptions):
    """ Empty class to assemble AcceptorSite PSSM settings """
    def __init__(self): pass
# end of class AcceptorPSSMOptions


################################################################################
# Which (Non)canonical splice sites are recognized/accepted?
################################################################################
CANONICAL_DONOR_SITES     = ['GT']
NON_CANONICAL_DONOR_SITES = ['GC']
CANONICAL_ACCEPTOR_SITES  = ['AG']


################################################################################
# PSSM_IC files for donors & acceptors
################################################################################
IC_DONOR_DATA_FILE         = osPathJoin(MAIN_ABGP_PATH,
                                        "datafiles/ic_donor_5fungi.txt")
IC_DONOR_PATTERN_OFFSET    = (3,4)
IC_ACCEPTOR_DATA_FILE      = osPathJoin(MAIN_ABGP_PATH,
                                        "datafiles/ic_acceptor_5fungi.txt")
IC_ACCEPTOR_PATTERN_OFFSET = (6,3)
IC_DONOR_NCGC_DATA_FILE    = osPathJoin(MAIN_ABGP_PATH,
                                        "datafiles/ic_ncdonorgc_5fungi.txt")

################################################################################
# Threshold values for donor site PSSM cutoffs
# Recommended to set spatiously low!
# Non-cannonical donor sites are supported yet, but disabled by default
################################################################################
MIN_DONOR_PSSM_SCORE                    = 0.0 
ALLOW_NON_CANONICAL_DONOR               = True
NON_CANONICAL_MIN_DONOR_PSSM_SCORE      = float(3.0)


################################################################################
# Threshold values for acceptor site PSSM cutoffs
# Recommended to set spatiously low!
# Non-cannonical acceptor splice sites are not supported yet
################################################################################
MIN_ACCEPTOR_PSSM_SCORE                 = 0.0       # 0.20 tinyexon acceptor in MYCGR_102837 tiny exon
ALLOW_NON_CANONICAL_ACCEPTOR            = False
NON_CANONICAL_MIN_ACCEPTOR_PSSM_SCORE   = float(0)


################################################################################
# Assemble DonorPSSM Thresholds into a single object
################################################################################
DONORPSSMOPTIONS = DonorPSSMOptions()
DONORPSSMOPTIONS.MIN_DONOR_PSSM_SCORE               = MIN_DONOR_PSSM_SCORE
DONORPSSMOPTIONS.ALLOW_NON_CANONICAL_DONOR          = ALLOW_NON_CANONICAL_DONOR
DONORPSSMOPTIONS.NON_CANONICAL_MIN_DONOR_PSSM_SCORE = NON_CANONICAL_MIN_DONOR_PSSM_SCORE
DONORPSSMOPTIONS.IC_DONOR_DATA_FILE                 = IC_DONOR_DATA_FILE
DONORPSSMOPTIONS.IC_DONOR_PATTERN_OFFSET            = IC_DONOR_PATTERN_OFFSET


################################################################################
# Assemble AcceptorPSSM Thresholds into a single object
################################################################################
ACCEPTORPSSMOPTIONS = AcceptorPSSMOptions()
ACCEPTORPSSMOPTIONS.MIN_ACCEPTOR_PSSM_SCORE               = MIN_ACCEPTOR_PSSM_SCORE
ACCEPTORPSSMOPTIONS.ALLOW_NON_CANONICAL_ACCEPTOR          = ALLOW_NON_CANONICAL_ACCEPTOR
ACCEPTORPSSMOPTIONS.NON_CANONICAL_MIN_ACCEPTOR_PSSM_SCORE = NON_CANONICAL_MIN_ACCEPTOR_PSSM_SCORE
ACCEPTORPSSMOPTIONS.IC_ACCEPTOR_DATA_FILE                 = IC_ACCEPTOR_DATA_FILE
ACCEPTORPSSMOPTIONS.IC_ACCEPTOR_PATTERN_OFFSET            = IC_ACCEPTOR_PATTERN_OFFSET



################################################################################
# global variables for stopless 3n introns 
################################################################################
STOPLESS_3N_INTRON_MIN_DONOR_PSSM_SCORE                  = 0.0 
STOPLESS_3N_INTRON_ALLOW_NON_CANONICAL_DONOR             = ALLOW_NON_CANONICAL_DONOR
STOPLESS_3N_INTRON_NON_CANONICAL_MIN_DONOR_PSSM_SCORE    = 1.5 
STOPLESS_3N_INTRON_MIN_ACCEPTOR_PSSM_SCORE               = 1.5 # increased acceptor affinity 
STOPLESS_3N_INTRON_ALLOW_NON_CANONICAL_ACCEPTOR          = ALLOW_NON_CANONICAL_ACCEPTOR
STOPLESS_3N_INTRON_NON_CANONICAL_MIN_ACCEPTOR_PSSM_SCORE = NON_CANONICAL_MIN_ACCEPTOR_PSSM_SCORE
STOPLESS_3N_INTRON_MIN_INTRON_NT_LENGTH                  = MIN_INTRON_NT_LENGTH
STOPLESS_3N_INTRON_MAX_INTRON_NT_LENGTH                  = MAX_INTRON_NT_LENGTH



################################################################################
# global variables for tinyexons and their surrounding introns
################################################################################
TINYEXON_MIN_NT_LENGTH                          = 8 
TINYEXON_MAX_NT_LENGTH                          = 21
TINYEXON_MIN_DONOR_PSSM_SCORE                   = 3.0
TINYEXON_MIN_ACCEPTOR_PSSM_SCORE                = 2.5
TINYEXON_ALLOW_NON_CANONICAL_DONOR              = True
TINYEXON_ALLOW_NON_CANONICAL_ACCEPTOR           = False
TINYEXON_NON_CANONICAL_MIN_DONOR_PSSM_SCORE     = 3.5
TINYEXON_NON_CANONICAL_MIN_ACCEPTOR_PSSM_SCORE  = 0.0
TINYEXON_MIN_TOTAL_PSSM_SCORE                   = TINYEXON_MIN_DONOR_PSSM_SCORE +\
                                                  TINYEXON_MIN_ACCEPTOR_PSSM_SCORE +\
                                                  1.0
TINYEXON_MIN_INTRON_NT_LENGTH                   = MIN_INTRON_NT_LENGTH
TINYEXON_MAX_INTRON_NT_LENGTH                   = 140


################################################################################
# global variables for detection of falsely PACBPORF alignments
# spanning stopless 3n introns 
################################################################################
PACBPORF_ALIGNED_STOPLESS_3N_INTRON_MIN_DONOR_PSSM_SCORE                  = 3.0 
PACBPORF_ALIGNED_STOPLESS_3N_INTRON_ALLOW_NON_CANONICAL_DONOR             = ALLOW_NON_CANONICAL_DONOR
PACBPORF_ALIGNED_STOPLESS_3N_INTRON_NON_CANONICAL_MIN_DONOR_PSSM_SCORE    = 3.0 
PACBPORF_ALIGNED_STOPLESS_3N_INTRON_MIN_ACCEPTOR_PSSM_SCORE               = 2.5
PACBPORF_ALIGNED_STOPLESS_3N_INTRON_ALLOW_NON_CANONICAL_ACCEPTOR          = False
PACBPORF_ALIGNED_STOPLESS_3N_INTRON_NON_CANONICAL_MIN_ACCEPTOR_PSSM_SCORE = NON_CANONICAL_MIN_ACCEPTOR_PSSM_SCORE
PACBPORF_ALIGNED_STOPLESS_3N_INTRON_MIN_INTRON_NT_LENGTH                  = 35
PACBPORF_ALIGNED_STOPLESS_3N_INTRON_MAX_INTRON_NT_LENGTH                  = 90
PACBPORF_ALIGNED_STOPLESS_3N_INTRON_HAS_BRANCHPOINT                       = True
PACBPORF_ALIGNED_STOPLESS_3N_INTRON_HAS_POLYPYRIMIDINE                    = True 
PACBPORF_ALIGNED_STOPLESS_3N_INTRON_HAS_AT_INCREASE                       = True
PACBPORF_ALIGNED_STOPLESS_3N_INTRON_ALLOW_PHASE_SHIFT                     = False # hard-coded in function
PACBPORF_ALIGNED_STOPLESS_3N_INTRON_MAX_NT_OFFSET                         = 3
PACBPORF_ALIGNED_STOPLESS_3N_INTRON_MAX_AA_OFFSET                         = PACBPORF_ALIGNED_STOPLESS_3N_INTRON_MAX_NT_OFFSET / 3


################################################################################
# global variables for projected introns
################################################################################
PROJECTED_INTRON_MIN_DONOR_PSSM_SCORE                   = MIN_DONOR_PSSM_SCORE
PROJECTED_INTRON_ALLOW_NON_CANONICAL_DONOR              = ALLOW_NON_CANONICAL_DONOR
PROJECTED_INTRON_NON_CANONICAL_MIN_DONOR_PSSM_SCORE     = NON_CANONICAL_MIN_DONOR_PSSM_SCORE
PROJECTED_INTRON_MIN_ACCEPTOR_PSSM_SCORE                = MIN_ACCEPTOR_PSSM_SCORE
PROJECTED_INTRON_ALLOW_NON_CANONICAL_ACCEPTOR           = ALLOW_NON_CANONICAL_ACCEPTOR
PROJECTED_INTRON_NON_CANONICAL_MIN_ACCEPTOR_PSSM_SCORE  = NON_CANONICAL_MIN_ACCEPTOR_PSSM_SCORE
PROJECTED_INTRON_MIN_INTRON_NT_LENGTH                   = MIN_INTRON_NT_LENGTH
PROJECTED_INTRON_MAX_INTRON_NT_LENGTH                   = MAX_INTRON_NT_LENGTH
PROJECTED_INTRON_ALLOW_PHASE_SHIFT                      = None                  # n.a.!
PROJECTED_INTRON_MAX_NT_OFFSET                          = 6                     # keep small
PROJECTED_INTRON_MAX_AA_OFFSET                          = PROJECTED_INTRON_MAX_NT_OFFSET / 3


################################################################################
# global variables for mapped introns
################################################################################
MAPPED_INTRON_MIN_DONOR_PSSM_SCORE                  = MIN_DONOR_PSSM_SCORE
MAPPED_INTRON_ALLOW_NON_CANONICAL_DONOR             = ALLOW_NON_CANONICAL_DONOR
MAPPED_INTRON_NON_CANONICAL_MIN_DONOR_PSSM_SCORE    = NON_CANONICAL_MIN_DONOR_PSSM_SCORE
MAPPED_INTRON_MIN_ACCEPTOR_PSSM_SCORE               = MIN_ACCEPTOR_PSSM_SCORE
MAPPED_INTRON_ALLOW_NON_CANONICAL_ACCEPTOR          = ALLOW_NON_CANONICAL_ACCEPTOR
MAPPED_INTRON_NON_CANONICAL_MIN_ACCEPTOR_PSSM_SCORE = NON_CANONICAL_MIN_ACCEPTOR_PSSM_SCORE
MAPPED_INTRON_MIN_INTRON_NT_LENGTH                  = MIN_INTRON_NT_LENGTH
MAPPED_INTRON_MAX_INTRON_NT_LENGTH                  = MAX_INTRON_NT_LENGTH
MAPPED_INTRON_ALLOW_PHASE_SHIFT                     = False                     # keep False
MAPPED_INTRON_MAX_NT_OFFSET                         = 6                         # keep small
MAPPED_INTRON_MAX_AA_OFFSET                         = MAPPED_INTRON_MAX_NT_OFFSET / 3


################################################################################
# global variables for mapped introns with strict acceptor, loose donor
################################################################################
MAPPED_INTRON_STRACC_MIN_DONOR_PSSM_SCORE                  = 3.5 # 4.16 in CFU_832284
MAPPED_INTRON_STRACC_ALLOW_NON_CANONICAL_DONOR             = ALLOW_NON_CANONICAL_DONOR
MAPPED_INTRON_STRACC_NON_CANONICAL_MIN_DONOR_PSSM_SCORE    = 3.5 
MAPPED_INTRON_STRACC_MIN_ACCEPTOR_PSSM_SCORE               = 5.0
MAPPED_INTRON_STRACC_ALLOW_NON_CANONICAL_ACCEPTOR          = ALLOW_NON_CANONICAL_ACCEPTOR
MAPPED_INTRON_STRACC_NON_CANONICAL_MIN_ACCEPTOR_PSSM_SCORE = NON_CANONICAL_MIN_ACCEPTOR_PSSM_SCORE
MAPPED_INTRON_STRACC_MIN_INTRON_NT_LENGTH                  = MIN_INTRON_NT_LENGTH
MAPPED_INTRON_STRACC_MAX_INTRON_NT_LENGTH                  = MAX_INTRON_NT_LENGTH
MAPPED_INTRON_STRACC_ALLOW_PHASE_SHIFT                     = False              # keep False
MAPPED_INTRON_STRACC_MAX_NT_OFFSET                         = 15                 # fairly keep small
MAPPED_INTRON_STRACC_MAX_AA_OFFSET                         = MAPPED_INTRON_STRACC_MAX_NT_OFFSET / 3


################################################################################
# global variables for mapped introns with strict donor, loose acceptor 
################################################################################
MAPPED_INTRON_STRDON_MIN_DONOR_PSSM_SCORE                  = 6.0
MAPPED_INTRON_STRDON_ALLOW_NON_CANONICAL_DONOR             = ALLOW_NON_CANONICAL_DONOR
MAPPED_INTRON_STRDON_NON_CANONICAL_MIN_DONOR_PSSM_SCORE    = 6.0
MAPPED_INTRON_STRDON_MIN_ACCEPTOR_PSSM_SCORE               = 3.5
MAPPED_INTRON_STRDON_ALLOW_NON_CANONICAL_ACCEPTOR          = ALLOW_NON_CANONICAL_ACCEPTOR
MAPPED_INTRON_STRDON_NON_CANONICAL_MIN_ACCEPTOR_PSSM_SCORE = NON_CANONICAL_MIN_ACCEPTOR_PSSM_SCORE
MAPPED_INTRON_STRDON_MIN_INTRON_NT_LENGTH                  = MIN_INTRON_NT_LENGTH
MAPPED_INTRON_STRDON_MAX_INTRON_NT_LENGTH                  = MAX_INTRON_NT_LENGTH
MAPPED_INTRON_STRDON_ALLOW_PHASE_SHIFT                     = False              # keep False
MAPPED_INTRON_STRDON_MAX_NT_OFFSET                         = 15                 # fairly keep small
MAPPED_INTRON_STRDON_MAX_AA_OFFSET                         = MAPPED_INTRON_STRDON_MAX_NT_OFFSET / 3


################################################################################
# global variables for phase shifted (mapped) introns
################################################################################
PHASE_SHIFT_INTRON_MIN_DONOR_PSSM_SCORE                 = 4.0
PHASE_SHIFT_INTRON_ALLOW_NON_CANONICAL_DONOR            = True
PHASE_SHIFT_INTRON_NON_CANONICAL_MIN_DONOR_PSSM_SCORE   = NON_CANONICAL_MIN_DONOR_PSSM_SCORE
PHASE_SHIFT_INTRON_MIN_ACCEPTOR_PSSM_SCORE              = 2.0
PHASE_SHIFT_INTRON_ALLOW_NON_CANONICAL_ACCEPTOR         = False
PHASE_SHIFT_INTRON_NON_CANONICAL_MIN_ACCEPTOR_PSSM_SCORE= 0.0                   # n.a.
PHASE_SHIFT_INTRON_MIN_INTRON_NT_LENGTH                 = MIN_INTRON_NT_LENGTH
PHASE_SHIFT_INTRON_MAX_INTRON_NT_LENGTH                 = 200
PHASE_SHIFT_INTRON_ALLOW_PHASE_SHIFT                    = True                  # keep True
# DISTANCE vs. OFFSET: Distance is measured distance between phase shifted splice sites
# DISTANCE vs. OFFSET: Offset is max length of the gap due to the phase shift
PHASE_SHIFT_INTRON_MAX_NT_OFFSET                        = 3
PHASE_SHIFT_INTRON_MAX_AA_OFFSET                        = 1
PHASE_SHIFT_INTRON_MAX_NT_DISTANCE                      = 5                     # keep short!
PHASE_SHIFT_INTRON_MAX_AA_DISTANCE                      = PHASE_SHIFT_INTRON_MAX_NT_DISTANCE / 3
if (PHASE_SHIFT_INTRON_MAX_NT_DISTANCE % 3) > 0:
    PHASE_SHIFT_INTRON_MAX_AA_DISTANCE += 1



################################################################################
# global variables for closeby independant intron gain
################################################################################
CLOSEBY_INDEPENDANT_INTRON_GAIN_MIN_DONOR_PSSM_SCORE                    = MIN_DONOR_PSSM_SCORE
CLOSEBY_INDEPENDANT_INTRON_GAIN_ALLOW_NON_CANONICAL_DONOR               = False
CLOSEBY_INDEPENDANT_INTRON_GAIN_NON_CANONICAL_MIN_DONOR_PSSM_SCORE      = 0.0   # n.a.
CLOSEBY_INDEPENDANT_INTRON_GAIN_MIN_ACCEPTOR_PSSM_SCORE                 = MIN_ACCEPTOR_PSSM_SCORE
CLOSEBY_INDEPENDANT_INTRON_GAIN_ALLOW_NON_CANONICAL_ACCEPTOR            = False
CLOSEBY_INDEPENDANT_INTRON_GAIN_NON_CANONICAL_MIN_ACCEPTOR_PSSM_SCORE   = 0.0   # n.a.
CLOSEBY_INDEPENDANT_INTRON_GAIN_MIN_INTRON_NT_LENGTH                    = MIN_INTRON_NT_LENGTH
CLOSEBY_INDEPENDANT_INTRON_GAIN_MAX_INTRON_NT_LENGTH                    = 200
CLOSEBY_INDEPENDANT_INTRON_GAIN_ALLOW_PHASE_SHIFT                       = True  # keep True
CLOSEBY_INDEPENDANT_INTRON_GAIN_MAX_AA_OFFSET                           = 1     # 0 or 1 AA!!
CLOSEBY_INDEPENDANT_INTRON_GAIN_MAX_NT_OFFSET                           = ( CLOSEBY_INDEPENDANT_INTRON_GAIN_MAX_AA_OFFSET * 3 ) + 2
# Length of the resulting intermediate PacbPORF due to C.I.G.
CLOSEBY_INDEPENDANT_INTRON_GAIN_MIN_AA_LENGTH                           = 2     # keep short!
CLOSEBY_INDEPENDANT_INTRON_GAIN_MAX_AA_LENGTH                           = 7     # keep short!



################################################################################
# global variables for projected tinyexons - 2 introns
################################################################################
PROJECTED_TINYEXON_MIN_DONOR_PSSM_SCORE                   = MIN_DONOR_PSSM_SCORE
PROJECTED_TINYEXON_ALLOW_NON_CANONICAL_DONOR              = ALLOW_NON_CANONICAL_DONOR
PROJECTED_TINYEXON_NON_CANONICAL_MIN_DONOR_PSSM_SCORE     = NON_CANONICAL_MIN_DONOR_PSSM_SCORE
PROJECTED_TINYEXON_MIN_ACCEPTOR_PSSM_SCORE                = MIN_ACCEPTOR_PSSM_SCORE
PROJECTED_TINYEXON_ALLOW_NON_CANONICAL_ACCEPTOR           = ALLOW_NON_CANONICAL_ACCEPTOR
PROJECTED_TINYEXON_NON_CANONICAL_MIN_ACCEPTOR_PSSM_SCORE  = NON_CANONICAL_MIN_ACCEPTOR_PSSM_SCORE
PROJECTED_TINYEXON_MIN_INTRON_NT_LENGTH                   = TINYEXON_MIN_INTRON_NT_LENGTH
PROJECTED_TINYEXON_MAX_INTRON_NT_LENGTH                   = 250 # TINYEXON_MAX_INTRON_NT_LENGTH
PROJECTED_TINYEXON_MIN_EXON_NT_LENGTH                     = TINYEXON_MIN_NT_LENGTH
PROJECTED_TINYEXON_MAX_EXON_NT_LENGTH                     = TINYEXON_MAX_NT_LENGTH
PROJECTED_TINYEXON_ALLOW_PHASE_SHIFT                      = None                  # n.a.!
PROJECTED_TINYEXON_MAX_NT_OFFSET                          = 3                     # keep small
PROJECTED_TINYEXON_MAX_AA_OFFSET                          = PROJECTED_INTRON_MAX_NT_OFFSET / 3



KWARGS_SPLICESITES = {
    'min_donor_pssm_score':                 MIN_DONOR_PSSM_SCORE,
    'allow_non_canonical_donor':            ALLOW_NON_CANONICAL_DONOR,
    'non_canonical_min_donor_pssm_score':   NON_CANONICAL_MIN_DONOR_PSSM_SCORE,
    'min_acceptor_pssm_score':              MIN_ACCEPTOR_PSSM_SCORE,
    'allow_non_canonical_acceptor':         ALLOW_NON_CANONICAL_ACCEPTOR,
    'non_canonical_min_acceptor_pssm_score':NON_CANONICAL_MIN_ACCEPTOR_PSSM_SCORE,
}



KWARGS_STOPLESS_3N_INTRONS = {
    'min_donor_pssm_score':                 STOPLESS_3N_INTRON_MIN_DONOR_PSSM_SCORE,
    'allow_non_canonical_donor':            STOPLESS_3N_INTRON_ALLOW_NON_CANONICAL_DONOR,
    'non_canonical_min_donor_pssm_score':   STOPLESS_3N_INTRON_NON_CANONICAL_MIN_DONOR_PSSM_SCORE,
    'min_acceptor_pssm_score':              STOPLESS_3N_INTRON_MIN_ACCEPTOR_PSSM_SCORE,
    'allow_non_canonical_acceptor':         STOPLESS_3N_INTRON_ALLOW_NON_CANONICAL_ACCEPTOR,
    'non_canonical_min_acceptor_pssm_score':STOPLESS_3N_INTRON_NON_CANONICAL_MIN_ACCEPTOR_PSSM_SCORE,
    'min_intron_nt_length':                 STOPLESS_3N_INTRON_MIN_INTRON_NT_LENGTH,
    'max_intron_nt_length':                 STOPLESS_3N_INTRON_MAX_INTRON_NT_LENGTH,
}


KWARGS_PACBPORF_ALIGNED_STOPLESS_3N_INTRONS = {
    'min_donor_pssm_score':                 PACBPORF_ALIGNED_STOPLESS_3N_INTRON_MIN_DONOR_PSSM_SCORE,
    'allow_non_canonical_donor':            PACBPORF_ALIGNED_STOPLESS_3N_INTRON_ALLOW_NON_CANONICAL_DONOR,
    'non_canonical_min_donor_pssm_score':   PACBPORF_ALIGNED_STOPLESS_3N_INTRON_NON_CANONICAL_MIN_DONOR_PSSM_SCORE,
    'min_acceptor_pssm_score':              PACBPORF_ALIGNED_STOPLESS_3N_INTRON_MIN_ACCEPTOR_PSSM_SCORE,
    'allow_non_canonical_acceptor':         PACBPORF_ALIGNED_STOPLESS_3N_INTRON_ALLOW_NON_CANONICAL_ACCEPTOR,
    'non_canonical_min_acceptor_pssm_score':PACBPORF_ALIGNED_STOPLESS_3N_INTRON_NON_CANONICAL_MIN_ACCEPTOR_PSSM_SCORE,
    'min_intron_nt_length':                 PACBPORF_ALIGNED_STOPLESS_3N_INTRON_MIN_INTRON_NT_LENGTH,
    'max_intron_nt_length':                 PACBPORF_ALIGNED_STOPLESS_3N_INTRON_MAX_INTRON_NT_LENGTH,
    'has_branchpoint':                      PACBPORF_ALIGNED_STOPLESS_3N_INTRON_HAS_BRANCHPOINT,
    'has_polypyrimidine':                   PACBPORF_ALIGNED_STOPLESS_3N_INTRON_HAS_POLYPYRIMIDINE,
    'has_at_increase':                      PACBPORF_ALIGNED_STOPLESS_3N_INTRON_HAS_AT_INCREASE,
    'max_nt_offset':                        PACBPORF_ALIGNED_STOPLESS_3N_INTRON_MAX_NT_OFFSET,
    'max_aa_offset':                        PACBPORF_ALIGNED_STOPLESS_3N_INTRON_MAX_AA_OFFSET,
}


KWARGS_TINYEXON_PAIRWISE = {
    'min_donor_pssm_score':                 TINYEXON_MIN_DONOR_PSSM_SCORE,
    'allow_non_canonical_donor':            TINYEXON_ALLOW_NON_CANONICAL_DONOR,
    'non_canonical_min_donor_pssm_score':   TINYEXON_NON_CANONICAL_MIN_DONOR_PSSM_SCORE,
    'min_acceptor_pssm_score':              TINYEXON_MIN_ACCEPTOR_PSSM_SCORE,
    'allow_non_canonical_acceptor':         TINYEXON_ALLOW_NON_CANONICAL_ACCEPTOR,
    'non_canonical_min_acceptor_pssm_score':TINYEXON_NON_CANONICAL_MIN_ACCEPTOR_PSSM_SCORE,
    'min_total_pssm_score':                 TINYEXON_MIN_TOTAL_PSSM_SCORE,
    'max_tinyexon_nt_length':               TINYEXON_MAX_NT_LENGTH,
    'min_tinyexon_nt_length':               TINYEXON_MIN_NT_LENGTH,
    'min_intron_nt_length':                 TINYEXON_MIN_INTRON_NT_LENGTH,
    'max_intron_nt_length':                 TINYEXON_MAX_INTRON_NT_LENGTH,
}



KWARGS_PROJECTED_INTRON = {
    'min_donor_pssm_score':                 PROJECTED_INTRON_MIN_DONOR_PSSM_SCORE,
    'allow_non_canonical_donor':            PROJECTED_INTRON_ALLOW_NON_CANONICAL_DONOR,
    'non_canonical_min_donor_pssm_score':   PROJECTED_INTRON_NON_CANONICAL_MIN_DONOR_PSSM_SCORE,
    'min_acceptor_pssm_score':              PROJECTED_INTRON_MIN_ACCEPTOR_PSSM_SCORE,
    'allow_non_canonical_acceptor':         PROJECTED_INTRON_ALLOW_NON_CANONICAL_ACCEPTOR,
    'non_canonical_min_acceptor_pssm_score':PROJECTED_INTRON_NON_CANONICAL_MIN_ACCEPTOR_PSSM_SCORE,
    'min_intron_nt_length':                 PROJECTED_INTRON_MIN_INTRON_NT_LENGTH,
    'max_intron_nt_length':                 PROJECTED_INTRON_MAX_INTRON_NT_LENGTH,
    # phase shift is not relevant for intron projection
    'max_nt_offset':                        PROJECTED_INTRON_MAX_NT_OFFSET,
    'max_aa_offset':                        PROJECTED_INTRON_MAX_AA_OFFSET,
}


KWARGS_MAPPED_INTRON_STRACC = {
    'min_donor_pssm_score':                 MAPPED_INTRON_STRACC_MIN_DONOR_PSSM_SCORE,
    'allow_non_canonical_donor':            MAPPED_INTRON_STRACC_ALLOW_NON_CANONICAL_DONOR,
    'non_canonical_min_donor_pssm_score':   MAPPED_INTRON_STRACC_NON_CANONICAL_MIN_DONOR_PSSM_SCORE,
    'min_acceptor_pssm_score':              MAPPED_INTRON_STRACC_MIN_ACCEPTOR_PSSM_SCORE,
    'allow_non_canonical_acceptor':         MAPPED_INTRON_STRACC_ALLOW_NON_CANONICAL_ACCEPTOR,
    'non_canonical_min_acceptor_pssm_score':MAPPED_INTRON_STRACC_NON_CANONICAL_MIN_ACCEPTOR_PSSM_SCORE,
    'min_intron_nt_length':                 MAPPED_INTRON_STRACC_MIN_INTRON_NT_LENGTH,
    'max_intron_nt_length':                 MAPPED_INTRON_STRACC_MAX_INTRON_NT_LENGTH,
    'allow_phase_shift':                    MAPPED_INTRON_STRACC_ALLOW_PHASE_SHIFT,
    'max_nt_offset':                        MAPPED_INTRON_STRACC_MAX_NT_OFFSET,
    'max_aa_offset':                        MAPPED_INTRON_STRACC_MAX_AA_OFFSET,
}


KWARGS_MAPPED_INTRON_STRDON = {
    'min_donor_pssm_score':                 MAPPED_INTRON_STRDON_MIN_DONOR_PSSM_SCORE,
    'allow_non_canonical_donor':            MAPPED_INTRON_STRDON_ALLOW_NON_CANONICAL_DONOR,
    'non_canonical_min_donor_pssm_score':   MAPPED_INTRON_STRDON_NON_CANONICAL_MIN_DONOR_PSSM_SCORE,
    'min_acceptor_pssm_score':              MAPPED_INTRON_STRDON_MIN_ACCEPTOR_PSSM_SCORE,
    'allow_non_canonical_acceptor':         MAPPED_INTRON_STRDON_ALLOW_NON_CANONICAL_ACCEPTOR,
    'non_canonical_min_acceptor_pssm_score':MAPPED_INTRON_STRDON_NON_CANONICAL_MIN_ACCEPTOR_PSSM_SCORE,
    'min_intron_nt_length':                 MAPPED_INTRON_STRDON_MIN_INTRON_NT_LENGTH,
    'max_intron_nt_length':                 MAPPED_INTRON_STRDON_MAX_INTRON_NT_LENGTH,
    'allow_phase_shift':                    MAPPED_INTRON_STRDON_ALLOW_PHASE_SHIFT,
    'max_nt_offset':                        MAPPED_INTRON_STRDON_MAX_NT_OFFSET,
    'max_aa_offset':                        MAPPED_INTRON_STRDON_MAX_AA_OFFSET,
}


KWARGS_MAPPED_INTRON = {
    'min_donor_pssm_score':                 MAPPED_INTRON_MIN_DONOR_PSSM_SCORE,
    'allow_non_canonical_donor':            MAPPED_INTRON_ALLOW_NON_CANONICAL_DONOR,
    'non_canonical_min_donor_pssm_score':   MAPPED_INTRON_NON_CANONICAL_MIN_DONOR_PSSM_SCORE,
    'min_acceptor_pssm_score':              MAPPED_INTRON_MIN_ACCEPTOR_PSSM_SCORE,
    'allow_non_canonical_acceptor':         MAPPED_INTRON_ALLOW_NON_CANONICAL_ACCEPTOR,
    'non_canonical_min_acceptor_pssm_score':MAPPED_INTRON_NON_CANONICAL_MIN_ACCEPTOR_PSSM_SCORE,
    'min_intron_nt_length':                 MAPPED_INTRON_MIN_INTRON_NT_LENGTH,
    'max_intron_nt_length':                 MAPPED_INTRON_MAX_INTRON_NT_LENGTH,
    'allow_phase_shift':                    MAPPED_INTRON_ALLOW_PHASE_SHIFT,
    'max_nt_offset':                        MAPPED_INTRON_MAX_NT_OFFSET,
    'max_aa_offset':                        MAPPED_INTRON_MAX_AA_OFFSET,
}

KWARGS_CLOSEBY_INDEPENDANT_INTRON_GAIN = {
    'min_donor_pssm_score':                 CLOSEBY_INDEPENDANT_INTRON_GAIN_MIN_DONOR_PSSM_SCORE,
    'allow_non_canonical_donor':            CLOSEBY_INDEPENDANT_INTRON_GAIN_ALLOW_NON_CANONICAL_DONOR,
    'non_canonical_min_donor_pssm_score':   CLOSEBY_INDEPENDANT_INTRON_GAIN_NON_CANONICAL_MIN_DONOR_PSSM_SCORE,
    'min_acceptor_pssm_score':              CLOSEBY_INDEPENDANT_INTRON_GAIN_MIN_ACCEPTOR_PSSM_SCORE,
    'allow_non_canonical_acceptor':         CLOSEBY_INDEPENDANT_INTRON_GAIN_ALLOW_NON_CANONICAL_ACCEPTOR,
    'non_canonical_min_acceptor_pssm_score':CLOSEBY_INDEPENDANT_INTRON_GAIN_NON_CANONICAL_MIN_ACCEPTOR_PSSM_SCORE,
    'min_intron_nt_length':                 CLOSEBY_INDEPENDANT_INTRON_GAIN_MIN_INTRON_NT_LENGTH,
    'max_intron_nt_length':                 CLOSEBY_INDEPENDANT_INTRON_GAIN_MAX_INTRON_NT_LENGTH,
    'allow_phase_shift':                    CLOSEBY_INDEPENDANT_INTRON_GAIN_ALLOW_PHASE_SHIFT,
    'max_nt_offset':                        CLOSEBY_INDEPENDANT_INTRON_GAIN_MAX_NT_OFFSET,
    'max_aa_offset':                        CLOSEBY_INDEPENDANT_INTRON_GAIN_MAX_AA_OFFSET,
    # closeby independant intron gain specific attributes
    'cig_min_aa_length':                    CLOSEBY_INDEPENDANT_INTRON_GAIN_MIN_AA_LENGTH,
    'cig_max_aa_length':                    CLOSEBY_INDEPENDANT_INTRON_GAIN_MAX_AA_LENGTH,
}



KWARGS_PHASE_SHIFT_INTRON = {
    'min_donor_pssm_score':                 PHASE_SHIFT_INTRON_MIN_DONOR_PSSM_SCORE,
    'allow_non_canonical_donor':            PHASE_SHIFT_INTRON_ALLOW_NON_CANONICAL_DONOR,
    'non_canonical_min_donor_pssm_score':   PHASE_SHIFT_INTRON_NON_CANONICAL_MIN_DONOR_PSSM_SCORE,
    'min_acceptor_pssm_score':              PHASE_SHIFT_INTRON_MIN_ACCEPTOR_PSSM_SCORE,
    'allow_non_canonical_acceptor':         PHASE_SHIFT_INTRON_ALLOW_NON_CANONICAL_ACCEPTOR,
    'non_canonical_min_acceptor_pssm_score':PHASE_SHIFT_INTRON_NON_CANONICAL_MIN_ACCEPTOR_PSSM_SCORE,
    'min_intron_nt_length':                 PHASE_SHIFT_INTRON_MIN_INTRON_NT_LENGTH,
    'max_intron_nt_length':                 PHASE_SHIFT_INTRON_MAX_INTRON_NT_LENGTH,
    'allow_phase_shift':                    PHASE_SHIFT_INTRON_ALLOW_PHASE_SHIFT,
    'max_nt_offset':                        PHASE_SHIFT_INTRON_MAX_NT_OFFSET,
    'max_aa_offset':                        PHASE_SHIFT_INTRON_MAX_AA_OFFSET,
    # phase shift intron specific attributes
    'max_nt_distance':                      PHASE_SHIFT_INTRON_MAX_NT_DISTANCE,
    'max_aa_distance':                      PHASE_SHIFT_INTRON_MAX_AA_DISTANCE,
}

KWARGS_PROJECTED_TINYEXON = {
    'min_donor_pssm_score':                 PROJECTED_TINYEXON_MIN_DONOR_PSSM_SCORE,
    'allow_non_canonical_donor':            PROJECTED_TINYEXON_ALLOW_NON_CANONICAL_DONOR,
    'non_canonical_min_donor_pssm_score':   PROJECTED_TINYEXON_NON_CANONICAL_MIN_DONOR_PSSM_SCORE,
    'min_acceptor_pssm_score':              PROJECTED_TINYEXON_MIN_ACCEPTOR_PSSM_SCORE,
    'allow_non_canonical_acceptor':         PROJECTED_TINYEXON_ALLOW_NON_CANONICAL_ACCEPTOR,
    'non_canonical_min_acceptor_pssm_score':PROJECTED_TINYEXON_NON_CANONICAL_MIN_ACCEPTOR_PSSM_SCORE,
    'min_intron_nt_length':                 PROJECTED_TINYEXON_MIN_INTRON_NT_LENGTH,
    'max_intron_nt_length':                 PROJECTED_TINYEXON_MAX_INTRON_NT_LENGTH,
    'allow_phase_shift':                    PROJECTED_TINYEXON_ALLOW_PHASE_SHIFT,
    'max_nt_offset':                        PROJECTED_TINYEXON_MAX_NT_OFFSET,
    'max_aa_offset':                        PROJECTED_TINYEXON_MAX_AA_OFFSET,
    # projected tinyexon specific attributes
    'min_exon_nt_length':                   PROJECTED_TINYEXON_MIN_EXON_NT_LENGTH,
    'max_exon_nt_length':                   PROJECTED_TINYEXON_MAX_EXON_NT_LENGTH,
    'min_tinyexon_nt_length':               PROJECTED_TINYEXON_MIN_EXON_NT_LENGTH,
    'max_tinyexon_nt_length':               PROJECTED_TINYEXON_MAX_EXON_NT_LENGTH,
    'min_tinyexon_intron_nt_length':        PROJECTED_TINYEXON_MIN_INTRON_NT_LENGTH,
    'max_tinyexon_intron_nt_length':        PROJECTED_TINYEXON_MAX_INTRON_NT_LENGTH,
    'min_total_pssm_score':                 None, # not in use as threshold!
}

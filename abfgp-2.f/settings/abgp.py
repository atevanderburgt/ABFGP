"""
Alignment Based Gene Predictions settings: all paths and settings required
for succesful execution of the ABGP suite of algorithms
"""

from os.path import join as osPathJoin

################################################################################
# Main ABGP directory
################################################################################
MAIN_ABGP_PATH         = "/home/avdb/code/abfgp-2.0/"               # ADJUST TO YOUR PERSONAL SYSTEM
ABGP_OUTDIR_PATH       = "/tmp"                                     # ADJUST TO YOUR PERSONAL SYSTEM; e.g. /tmp/abfgp/results

# decorative variable used e.g. in GFF fsource
ABGP_PROGRAM_NAME      = "abfgp"
ABGP_VERSION_NUMBER    = "2.0"
ABGP_VERSION           = ABGP_PROGRAM_NAME + ABGP_VERSION_NUMBER

MAIN_ABFGP_PATH        = MAIN_ABGP_PATH     # alias -> to be DEPRECATED
ABFGP_VERSION          = ABGP_VERSION       # alias -> to be DEPRECATED


################################################################################
# Settings on main data sanity input 
################################################################################
ABGP_MINIMAL_NUM_LOCI       = 4     # minimal required loci to start/continue
ABGP_OPTIMAL_NUM_LOCI       = 5     # optimal number of loci; okay, it is a guess ;-)
ABGP_MAXIMAL_NUM_LOCI       = 5     # more loci -> explosion in calculation time!

ABGP_MINIMAL_NUM_LOCI_LOWER = 2     # 2 loci is not an option (yet)
ABGP_MINIMAL_NUM_LOCI_UPPER = 5     # at least 1 target and 2 informants needed
ABGP_MAXIMAL_NUM_LOCI_LOWER = 3     # more loci -> explosion in calculation time!
ABGP_MAXIMAL_NUM_LOCI_UPPER = 7     # more loci -> explosion in calculation time!

# The same numbers can be expressed in numbers of informants
# Substract -1 from all LOCI numbers (because the target isa locus itself too!) 
ABGP_MINIMAL_NUM_INFORMANTS       = ABGP_MINIMAL_NUM_LOCI -1
ABGP_OPTIMAL_NUM_INFORMANTS       = ABGP_OPTIMAL_NUM_LOCI -1
ABGP_MAXIMAL_NUM_INFORMANTS       = ABGP_MAXIMAL_NUM_LOCI -1
ABGP_MINIMAL_NUM_INFORMANTS_LOWER = ABGP_MINIMAL_NUM_LOCI_LOWER -1
ABGP_MINIMAL_NUM_INFORMANTS_UPPER = ABGP_MINIMAL_NUM_LOCI_UPPER -1
ABGP_MAXIMAL_NUM_INFORMANTS_LOWER = ABGP_MAXIMAL_NUM_LOCI_LOWER -1
ABGP_MAXIMAL_NUM_INFORMANTS_UPPER = ABGP_MAXIMAL_NUM_LOCI_UPPER -1


# minimal and maximal locus lengths.
# Morale:
# ABFGP expects a gene-encoding locus, including  some flanking sequence.
# Flanking sequence is required to detect incorrectly (truncated) annotated genes
# ABFGP is designed to tackle a --single-- gene locus.
# Few to none - fungal - genes will span 20kb in length.
# The very broad length range of 2-20kb allows for virtually all informants to be
# incorporated. It is not recommended to adjust these thresholds. Only do it in
# case you want to reannotate a - SINGLE - gene that is over 20kb in length
ABGP_MINIMAL_NT_LOCUS_LENGTH = 2000
ABGP_MAXIMAL_NT_LOCUS_LENGTH = 20000


################################################################################
# Settings for some main actions to take into account or not
################################################################################
USE_PSSM_SPLICE_SITES       = True  # If False, any GT...AG are scored equivalently
USE_PSSM_TSS                = True  # If False, any ATG is scored equivalently
USE_TINYEXON                = True  # If False, not searched for
                                    # Can only be used if USE_PSSM_SPLICE_SITES=True
USE_SHORT_LEADINGEXON       = True  # If False, not searched for
                                    # Can only be used if USE_PSSM_SPLICE_SITES=True
                                    # Can only be used if USE_PSSM_TSS=True
USE_SHORT_TAILINGEXON       = True  # If False, not searched for
                                    # Can only be used if USE_PSSM_SPLICE_SITES=True
                                    # Can only be used if USE_PSSM_TSS=True
USE_TCODE                   = True  # If False, not used as scoring
USE_SIGNALP                 = False # If False, not used upon scoring a TSS/ATG
                                    # NOT IMPLEMENTED (YET)
USE_SIGNALP_WITH_INTRON     = False # If False, not used upon scoring a TSS/ATG
                                    # NOT IMPLEMENTED (YET)
                                    # Can only be used if USE_SIGNALP=True
                                    # Can only be used if USE_PSSM_SPLICE_SITES=True


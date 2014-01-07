"""
Alignment Based Gene Predictions settings: genestructure requirements for
genes, exons, introns, intergenic and regions. As well settings for
tiny exons (leading, intermediate, tailing).
"""

# import BLOSUM62 matrix path
from emboss import BLOSUM62_PATH

################################################################################
# Thresholds for length range to allow for introns and exons
# Recommended to set spatiously large and small!
################################################################################
ABSOLUTE_MAX_INTRON_NT_LENGTH       = 1500  # larger introns than this are NOT accepted
MAX_INTRON_NT_LENGTH                =  900  # regular max intron size. When a size larger than this
                                    # 1000  # is encountered (in the intergenecity function(s)), the
                                            # cbgIF is checked. Accepted as intron only when:
                                            # (1) the cbgIF is_optimal()
                                            # (2) the 2th CBG does not have uniform TSS sites 
MIN_INTRON_NT_LENGTH                = 21    # smallest intron possible (7AA)
MAX_EXON_NT_LENGTH                  = 10000 # not implemented...
MIN_EXON_NT_LENGTH                  = 3     # not implemented ??
MIN_INTERGENIC_NT_LENGTH            = 300   # smallest intergenic distance possible
MAX_INTERGENIC_MIN_NT_LENGTH        = 550   # if all intergenic distances > MIN_INTERGENIC_NT_LENGTH,
                                            # as soon as a single one > MAX_INTERGENIC_MIN_NT_LENGTH,
                                            # we still threat this as an intergenic distance
AVERAGE_INTERGENIC_MIN_NT_LENGTH    = 400   # if average distance > AVERAGE_INTERGENIC_MIN_NT_LENGTH
                                            # we threat this as an intergenic distance, not an intron 

# according to Kupfer e.a. 2004
OPTIMAL_BRACNHPOINT_TO_ACCEPTOR_DISTANCE = range(13,21)
# maximal distance towards optimal distance to still allow as a perfect branchpoint;
# currently only used to filter for (strong) stopless3n introns
MAXIMAL_OPTIMAL_BRACNHPOINT_TO_ACCEPTOR_DISTANCE = 5



################################################################################
# Thresholds for lower and upper boundaries for a gene
# Recommended to set spatiously large and small!
################################################################################
GENE_MIN_EXON_COUNT                 = 1     # not implemented (yet)
GENE_MAX_EXON_COUNT                 = 20    # not implemented (yet)
GENE_MAX_NT_LENGTH                  = 6000  # not implemented (yet)




################################################################################
# Labels to place in Objects to recognize gene structure elements
################################################################################
FIRST_ANNOTATED_EXON_LABEL = "IS_FIRST"
FINAL_ANNOTATED_EXON_LABEL = "IS_FINAL"
IS_ANNOTATED_EXON_LABEL    = "IS_KNOWN"
ORF_IS_UNIGENE_LABEL       = "IS_UNIGENE"


################################################################################
# Thresholds from at which nt identity of informant alignments special
# threatment is allowed to circumvent out-of-frame alignments
################################################################################
INFORMANT_HIGH_NT_IDENTITY      = 0.94
MIN_QUARANTINE_HIGH_NT_IDENTITY = 0.94



STOPCODONGRAPH_MINIMAL_TOTAL_WEIGHT_RATIO = 0.4


################################################################################
# Thresholds for allowing some extra AA search space in the
# interface of 2 CBGs that are scaffolded by at least 1 mutual Orf
################################################################################
SCAFFOLD_GAP_OMSR_OFFSET            = 1    # KEEP VERY SMALL!! (0 or 1)
                                           # because we are in search for a tiny
                                           # exon, elongation of the aligned
                                           # sequence results in lower sensitivity
                                           # that cause higher chance for nonsense
                                           # alignments



################################################################################
# Thresholds for a short leading exon
# Non-cannonical splice sites are not supported yet
################################################################################
# length range of short leading exon and first intron
SHORT_LEADINGEXON_MIN_NT_LENGTH                         = 2 
SHORT_LEADINGEXON_MAX_NT_LENGTH                         = 70  #60 #30
SHORT_LEADINGEXON_MAX_INTRON_NT_LENGTH                  = 700 #500; increasement done because of mgg0528 MGG_06685
SHORT_LEADINGEXON_MIN_INTRON_NT_LENGTH                  = MIN_INTRON_NT_LENGTH
# thresholds for TSS of short leading exon
SHORT_LEADINGEXON_TSS_MIN_PSSM_SCORE                    = 1.0
SHORT_LEADINGEXON_TSS_ALLOW_NON_CANONICAL               = False     # not implemented yet
SHORT_LEADINGEXON_TSS_NON_CANONICAL_MIN_PSSM_SCORE      = float(0)  # not implemented yet
# thresholds for splice site signals of short leading exon
SHORT_LEADINGEXON_ELEGIABLE_ACCEPTOR_OMSR_NT_OFFSET     = 21
SHORT_LEADINGEXON_MIN_TOTAL_PSSM_SCORE                  = 5.0
SHORT_LEADINGEXON_MIN_DONOR_PSSM_SCORE                  = 1.0
SHORT_LEADINGEXON_MIN_ACCEPTOR_PSSM_SCORE               = 1.0
SHORT_LEADINGEXON_ALLOW_NON_CANONICAL_DONOR             = False     # not implemented yet
SHORT_LEADINGEXON_ALLOW_NON_CANONICAL_ACCEPTOR          = False     # not implemented yet
SHORT_LEADINGEXON_NON_CANONICAL_MIN_PSSM_SCORE          = None      # DEPRECATED
SHORT_LEADINGEXON_NON_CANONICAL_MIN_DONOR_PSSM_SCORE    = float(0)  # not implemented yet
SHORT_LEADINGEXON_NON_CANONICAL_MIN_ACCEPTOR_PSSM_SCORE = float(0)  # not implemented yet


################################################################################
# Thresholds for a short TAILING exon
# Non-cannonical splice sites are not supported yet
################################################################################
# length range of short TAILING exon and first intron
SHORT_TAILINGEXON_MIN_NT_LENGTH                         = 2 # was 5 nt 
SHORT_TAILINGEXON_MAX_NT_LENGTH                         = 80 # was 50 nt; longer exon chance is (very) small; chance of problems too
SHORT_TAILINGEXON_MAX_INTRON_NT_LENGTH                  = 200
SHORT_TAILINGEXON_MIN_INTRON_NT_LENGTH                  = MIN_INTRON_NT_LENGTH
SHORT_TAILINGEXON_TAKE_MAX_BEST_ACCEPTORS               = 10        # do not take to little;
                                                                    # 4-7 used to give problems in
                                                                    # exceptional cases
SHORT_TAILINGEXON_TAKE_MAX_BEST_ECGS                    = 50        # very high number; top scoring is not
                                                                    # alwaysthe best after clustalw alignments!
SHORT_TAILINGEXON_TAKE_MAX_BEST_CBGS                    = 15        # not so high number; top scoring CBGs 
                                                                    # should contain the bona fide final CBG 


# thresholds for splice site signals of short TAILING exon
SHORT_TAILINGEXON_MIN_TOTAL_PSSM_SCORE                  = 3.0
SHORT_TAILINGEXON_MIN_DONOR_PSSM_SCORE                  = 0.1 
SHORT_TAILINGEXON_MIN_ACCEPTOR_PSSM_SCORE               = 0.1 
SHORT_TAILINGEXON_ALLOW_NON_CANONICAL_DONOR             = False     # not implemented yet
SHORT_TAILINGEXON_ALLOW_NON_CANONICAL_ACCEPTOR          = False     # not implemented yet
SHORT_TAILINGEXON_NON_CANONICAL_MIN_DONOR_PSSM_SCORE    = float(0)  # not implemented yet
SHORT_TAILINGEXON_NON_CANONICAL_MIN_ACCEPTOR_PSSM_SCORE = float(0)  # not implemented yet


################################################################################
# Thresholds for a tiny intermediate exon
# Non-cannonical splice sites are not supported yet
################################################################################
TINYEXON_MATRIX_PATH                            = BLOSUM62_PATH
TINYEXON_MAX_NT_LENGTH                          = 32        # 40; old value was not implemented
TINYEXON_MIN_NT_LENGTH                          = 5         # 15; old value was not implemented
TINYEXON_MAX_INTRON_NT_LENGTH                   = 200       # much stricter as MAX_INTRON_LENGTH
TINYEXON_MIN_INTRON_NT_LENGTH                   = 35        # much stricter as MIN_INTRON_LENGTH

# TODO -> FUTURE DEPRECATED/MOVED to splicesites.py
TINYEXON_MIN_TOTAL_PSSM_SCORE                   = 5.0
TINYEXON_MIN_DONOR_PSSM_SCORE                   = 0.5       # othwise a case in example9 is missed; site at ncu:2852
TINYEXON_MIN_ACCEPTOR_PSSM_SCORE                = 0.5
TINYEXON_ALLOW_NON_CANONICAL_DONOR              = False     # not implemented yet 
TINYEXON_ALLOW_NON_CANONICAL_ACCEPTOR           = False     # not implemented yet
TINYEXON_NON_CANONICAL_MIN_DONOR_PSSM_SCORE     = float(0)  # not implemented yet 
TINYEXON_NON_CANONICAL_MIN_ACCEPTOR_PSSM_SCORE  = float(0)  # not implemented yet 


################################################################################
# Thresholds for removing a intermediary CBG with poor CBGinterfaces
################################################################################
INTERMEDIATE_CBG_REMOVAL_MAX_AA_LENGTH          = 20
INTERMEDIATE_CBG_REMOVAL_MAX_GTG_ID_RATIO       = 0.900 # 0.800 failed for some

"""
Alignment Based Gene Predictions settings: settings for iterative blastp
"""


################################################################################
# Thresholds for detection of `inframe introns`.
# An inframe intron is a case where two exons are located on the same ORF.
# TODO: thresholds are now rather intuitively chosen ;-)
# improve this by a more supported method/numbers
################################################################################

# maybe to be DEPRECATED by INFRAME_INTRON_MIN_AA_LENGTH
IGNORE_INTRON_IN_ORF_MIN_AA_GAP_SIZE          = 12      # ~minimal intron AA size
IGNORE_INTRON_IN_ORF_MIN_AVERAGE_SCORE_CUTOFF = 0.20    # 0.2 ``identity`` score per aligned AA position
IGNORE_INTRON_IN_ORF_MIN_WINDOW_SUMMED_SCORE  = 1.5     # total ``identity`` score summed over a window of IGNORE_INTRON_IN_ORF_MIN_AA_GAP_SIZE
                                                        # this number must be fairly low; already 2.0 predicted a non-existing intron in
                                                        # orf 73 of MGG_10961_dna (first junction, not the second)

INFRAME_INTRON_MIN_AA_LENGTH                  = 10      # ~minimal intron AA size
                                                        # Somewhat larger than
                                                        # INTRON_MIN_AA_LENGTH



# DEPRECATED !?!?
MINIMAL_PROJECTED_INTRON_SUMMEDSITE_ENTROPY = 0.25  # A measure for sequence similarity
                                                    # at the interface of the projected
# DEPRECATED !?!?
ALIGNED_INTRON_MIN_AA_LENGTH                = 12    # Minimal length (in AA) of the gained/lost
                                                    # intron in one of the species. When set to small,
                                                    # e.g. near MAX_PROJECTED_INTRON_AA_OFFSET,
                                                    # the difference between a poorly conserved
                                                    # area of the protein and a (small) gained/lost intron,
                                                    # disappears



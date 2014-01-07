"""
Alignment Based Gene Predictions settings: AlignedTranslationalStartSiteGraph
requirements for accepting the predicted first CodingBlockGraph (CBG) as the
True first CBG and some other variables.
"""

################################################################################
# Thresholds for lower and upper boundaries for an AlignedTranslationalStartSiteGraph 
# to be accepted as the first CodingBlockGraph in the GeneStructureGraph 
################################################################################

ALIGNEDTSSGRAPH_OPTIMALITY_MIN_TCODE      = 1.00
ALIGNEDTSSGRAPH_OPTIMALITY_MAX_TCODE      = 1.15 # not 1.20 as in StopCodons
ALIGNEDTSSGRAPH_OPTIMALITY_MIN_WEIGHT     = 0.30
ALIGNEDTSSGRAPH_OPTIMALITY_MAX_WEIGHT     = 0.80
ALIGNEDTSSGRAPH_OPTIMALITY_MIN_GTGWEAKEST = 0.725 # the default value of
                                                  # OPTIONS.minimal_gtgweakestnode
ALIGNEDTSSGRAPH_OPTIMALITY_MAX_GTGWEAKEST = 0.900



################################################################################
# Threshold for replacing an (optimal) aligned TSS with a more 5' located
# TSS in the same ORF. This followes the logic that he Ribosome can not
# `read-forward` in search for a marginally higher (PSSM) scoring TSS, but just
# accepts the first encountered TSS
################################################################################
OPTIMIZE_ALIGNEDTSSGRAPH_MINIMAL_PSSM_RATIO = 0.70

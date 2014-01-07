"""
Alignment Based Gene Predictions settings: AlignedStopCodonGraph requirements
for accepting the predicted final CodingBlockGraph (CBG) as the True final CBG.
"""

################################################################################
# Thresholds for lower and upper boundaries for a AlignedStopCodonGraph 
# to be accepted as the final CodingBlockGraph in the GeneStructureGraph 
################################################################################

ALIGNEDSTOPCODONGRAPH_OPTIMALITY_MIN_TCODE      = 1.00
ALIGNEDSTOPCODONGRAPH_OPTIMALITY_MAX_TCODE      = 1.20
ALIGNEDSTOPCODONGRAPH_OPTIMALITY_MIN_WEIGHT     = 0.30
ALIGNEDSTOPCODONGRAPH_OPTIMALITY_MAX_WEIGHT     = 0.80
ALIGNEDSTOPCODONGRAPH_OPTIMALITY_MIN_GTGWEAKEST = 0.725 # the default value of
                                                        # OPTIONS.minimal_gtgweakestnode
ALIGNEDSTOPCODONGRAPH_OPTIMALITY_MAX_GTGWEAKEST = 0.900


def MAX_OMSR_DIST_FUNCTION(stopcodongraph):
    """
    Get maximal alowed distance between OMSR and StopCodon based on threshold and GTG identity
    """
    MAX_OMSR_DIST = 8
    cbgident   = stopcodongraph._codingblockgraph.genetree().identity()
    return MAX_OMSR_DIST  / pow(cbgident,2)

# end of function MAX_OMSR_DIST_FUNCTION



"""
Alignment Based Gene Predictions settings: CodingBlockGraph (CBG) settings & requirements
"""
from math import sqrt


################################################################################
# Thresholds for CodingBlockGraphs
################################################################################
CBG_MIN_AA_LENGTH                       = 5
CBG_MINIMAL_OVERAL_SPANNING_RANGE_SIZE  = 3
CBG_SPRDIF_MIN_NODE_COUNT               = 3


################################################################################
# Thresholds for CodingBlockGraph overlap in adding them to the GSG
# and subsequently removing them again from the GSG
################################################################################
CBG_MAX_ALOWED_GSGINSERTION_OVERLAP_AA_LENGTH      = 10 # ???
CBG_MAX_ALOWED_GSGINSERTION_OVERLAP_RATIO          = 0.5
CBG_MAX_ALOWED_GSGREMOVAL_OVERLAP_AA_LENGTH        = 10 # ???
CBG_MAX_ALOWED_GSGREMOVAL_OVERLAP_RATIO            = 0.50


################################################################################
# Thresholds for CodingBlockGraph extension on spanningrange difference
################################################################################
CBG_EXTEND_LEFT_SPRDIF_MIN_AA_LENGTH        = 6
CBG_EXTEND_RIGTH_SPRDIF_MIN_AA_LENGTH       = 6
CBG_EXTEND_PACBPORF_MIN_LENGTH_RATIO        = 1.0  # Extended PacbPORF must be longer
CBG_EXTEND_PACBPORF_MIN_IDENTITYSCORE_RATIO = 0.7  # Extended PacbPORF must not have dropped
                                                   # to much in identityscore

################################################################################
# Thresholds and parameters for CodingBlockGraph (ClustalW) optimization
################################################################################
CBG_OPTIMIZE_MAXIMAL_IDENTITY            = 0.80 # optimize CBG only when identity < 0.80
                                                # optimization is somewhat time-consuming and
                                                # blastp has no problem correctly assigning
                                                # high-identity gene alignments
CBG_OPTIMIZE_CLUSTALW_GAP_SIZE           = 6    # split ClustalW pacbps on gapsize of 6
                                                # this enshures breaking up the often whay to
                                                # long ClustalW alignments
CBG_OPTIMIZE_MINIMAL_BITSCORE_RATIO      = 0.95 # allow longer ClustalW-PacbP when bitscore
                                                # is at least 95% of the shorter (blastp
                                                # derived) PacbP 
CBG_OPTIMIZE_MINIMAL_IDENTITY_RATIO      = 0.90 # allow longer ClustalW-PacbP when identityscore
                                                # is at least 90% of the shorter (blastp
                                                # derived) PacbP


# CBG_PACBPGAP_NEARBY_OMSR_OFFSET and CBG_PACBPGAP_NEARBY_OMSR_GAP_SIZE
# relate to each other! The area around the OMSR start/end which is checked is
# of checkAAsize = ( 2 * CBG_PACBPGAP_NEARBY_OMSR_OFFSET + 1 ).
# So, CBG_PACBPGAP_NEARBY_OMSR_GAP_SIZE must be <(<) as checkAAsize!
CBG_PACBPGAP_NEARBY_OMSR_OFFSET     = 3 # offset around OMSR to check for gaps
CBG_PACBPGAP_NEARBY_OMSR_GAP_SIZE   = 3 # minimal required AA gap size



# thresholds for the cexpander_omsrbordergaps function
CBG_CEXPANDER_OMSRBORDERGAPS_NONUNIFORM_MAX_AA_LENGTH    = 15
CBG_CEXPANDER_OMSRBORDERGAPS_MAX_BITSCORERATIO_THRESHOLD = 0.0
CBG_CEXPANDER_OMSRBORDERGAPS_NONUNIFORM_AA_OFFSET        = 10
CBG_CEXPANDER_OMSRBORDERGAPS_GAP_SIZE                    = 3




################################################################################
# Thresholds for CodingBlockGraph splitting on spanningrange difference
################################################################################
# splitting on LARGE left sprdif
CBG_LARGE_LEFT_SPRDIF_MIN_AA_LENGTH         = 20
CBG_LARGE_LEFT_SPRDIF_MIN_NODE_COUNT        = 2     #EXACT_SG_NODE_COUNT-2

# splitting on LARGE rigth sprdif
CBG_LARGE_RIGTH_SPRDIF_MIN_AA_LENGTH        = 20
CBG_LARGE_RIGTH_SPRDIF_MIN_NODE_COUNT       = 2

CBG_FIRST_SPRDIF_MIN_AA_LENGTH              = 8
CBG_FIRST_SPRDIF_MIN_NODE_COUNT             = 2
CBG_SPRDIF_VS_OMSR_MIN_IDENTITYSCORE_RATIO  = 0.60

# splitting on SMALL left sprdif
CBG_SMALL_LEFT_SPRDIF_MIN_AA_LENGTH         = 5
CBG_SMALL_LEFT_SPRDIF_MIN_NODE_COUNT        = 2
CBG_SMALL_LEFT_SPRDIF_SPRDIF_MIN_IDENTITY   = 0.65
CBG_SMALL_LEFT_SPRDIF_CBG_MIN_IDENTITY      = 0.70

# splitting on SMALL rigth sprdif
CBG_SMALL_RIGTH_SPRDIF_MIN_AA_LENGTH        = 5
CBG_SMALL_RIGTH_SPRDIF_MIN_NODE_COUNT       = 2
CBG_SMALL_RIGTH_SPRDIF_SPRDIF_MIN_IDENTITY  = 0.65
CBG_SMALL_RIGTH_SPRDIF_CBG_MIN_IDENTITY     = 0.70


################################################################################
# Thresholds for inframe introns in CodingBlockGraphs
################################################################################
INFRAME_INTRON_IN_CBG_GAP_SIZE                      = 5
INFRAME_INTRON_IN_CBG_GAP_SIZE_ALONE                = 10
INFRAME_INTRON_IN_CBG_LENGTH_DISCREPANCY            = 10
INFRAME_INTRON_IN_CBG_WINDOW_AA_SIZE                = 12
INFRAME_INTRON_IN_CBG_WINDOW_MIN_SIMILARITY_SCORE   = 2.0
INFRAME_INTRON_IN_CBG_OBSERVED_VS_EXPECTED          = 1
INFRAME_INTRON_IN_CBG_MIN_TOTAL_PSSM                = 3.0



################################################################################
# Thresholds for accepting CodingBlockGraph into the GeneStructureOfCodingBlocks
# These thresholds are valid for a GeneTreeGRaph with identity == 1.0
# Threshold values are adjusted as a function of genetreegraph's identity 
################################################################################
MAX_CBG_GTG_TOPO_DIF                    = 0.075 # adjusted by  gtg.identity()    # was 0.075
MAX_CBG_GTG_ABS_DIF                     = 0.125 # adjusted by  gtg.identity()    # was 0.080 <20/11/2009
MIN_CBG_GTG_ID_RATIO                    = 0.700 # adjusted by  gtg.identity()    # was 0.850 <20/11/2009

# Threshold values for functions below as as function of GTG.identity()
#	TOPO	ABS	RATIO
# id%	0,075	0,750	0,850
# 1,00	0,075	0,750	0,850
# 0,95	0,079	0,831	0,825
# 0,90	0,083	0,926	0,799
# 0,85	0,088	1,038	0,772
# 0,80	0,094	1,172	0,744
# 0,75	0,100	1,333	0,716
# 0,70	0,107	1,531	0,687
# 0,60	0,125	2,083	0,625
# 0,50	0,150	3,000	0,557
# 0,40	0,188	4,688	0,482


#def MAX_GTG_TOPO_DIF_FUNCTION(value,gtg):
#    """ """
#    return value / gtg.identity() 
#
## end of function MAX_GTG_TOPO_DIF_FUNCTION
#
#def MAX_GTG_ABS_DIF_FUNCTION(value,gtg):
#    """ """
#    return value / ( pow(gtg.identity(),2) )
#
## end of function MAX_GTG_ABS_DIF_FUNCTION
#
#def MIN_GTG_ID_RATIO_FUNCTION(value,gtg):
#    """ """
#    return value - (1.0 - sqrt(gtg.identity()))
#
## end of function MIN_GTG_ID_RATIO_FUNCTION


def MAX_GTG_TOPO_DIF_FUNCTION(value,gtg,cbg):
    """ """
    MIN_OMSR_LEN = 20
    MAX_OMSR_LEN = 100 
    threshold = value / gtg.identity()
    omsrlength = cbg.omsrlength() 
    if omsrlength <= MIN_OMSR_LEN:
        return threshold
    elif omsrlength >= MAX_OMSR_LEN:
        return threshold * 2.0
    else:
        ratio = float(omsrlength - MIN_OMSR_LEN) / float(MAX_OMSR_LEN - MIN_OMSR_LEN)
        return threshold * (1.0 + ratio)

# end of function MAX_GTG_TOPO_DIF_FUNCTION


def MAX_GTG_ABS_DIF_FUNCTION(value,gtg,cbg):
    """ """
    MIN_OMSR_LEN = 20
    MAX_OMSR_LEN = 100
    threshold = value / ( pow(gtg.identity(),2) )
    omsrlength = cbg.omsrlength()
    if omsrlength <= MIN_OMSR_LEN:
        return threshold
    elif omsrlength >= MAX_OMSR_LEN:
        return threshold * 2.0
    else:
        ratio = float(omsrlength - MIN_OMSR_LEN) / float(MAX_OMSR_LEN - MIN_OMSR_LEN)
        return threshold * (1.0 + ratio)

# end of function MAX_GTG_ABS_DIF_FUNCTION


def MIN_GTG_ID_RATIO_FUNCTION(value,gtg,cbg):
    """ """
    MIN_OMSR_LEN = 20
    MAX_OMSR_LEN = 100
    threshold = value - (1.0 - sqrt(gtg.identity())) 
    omsrlength = cbg.omsrlength()
    if omsrlength <= MIN_OMSR_LEN:
        return threshold
    elif omsrlength >= MAX_OMSR_LEN:
        return threshold / 2.0
    else:
        ratio = float(omsrlength - MIN_OMSR_LEN) / float(MAX_OMSR_LEN - MIN_OMSR_LEN)
        return threshold / (1.0 + ratio)

# end of function MIN_GTG_ID_RATIO_FUNCTION



################################################################################
# Thresholds for accepting CodingBlockGraph into the GeneStructureOfCodingBlocks
# for the IS_FIRST small sprdif splitting
################################################################################
MAX_CBG_FIRST_SPRDIF_GTG_TOPO_DIF           = 0.110
MAX_CBG_FIRST_SPRDIF_GTG_ABS_DIF            = 0.150
MIN_CBG_FIRST_SPRDIF_GTG_ID_RATIO           = 0.600


################################################################################
# Thresholds for accepting CodingBlockGraph into the GeneStructureOfCodingBlocks
# for splitting on final sprdif difference, final cbg 
################################################################################
MAX_CBG_FINAL_SPRDIF_GTG_TOPO_DIF           = 0.110
MAX_CBG_FINAL_SPRDIF_GTG_ABS_DIF            = 0.150
MIN_CBG_FINAL_SPRDIF_GTG_ID_RATIO           = 0.700

# values for when & how to use this function
CBG_FINAL_SPRDIF_MIN_AA_LENGTH              = 7 # was 10, but migth be somewhat to large!
CBG_FINAL_SPRDIF_MIN_NODE_COUNT             = 2 
CBG_FINAL_SPRDIF_ONLY_IF_STOP_TW_RATIO_LTE  = 0.70
CBG_FINAL_SPRDIF_ONLY_IF_CBG_ID_GTE         = 0.80


################################################################################
# Thresholds for accepting CodingBlockGraph into the GeneStructureOfCodingBlocks
# after splitting on spanningrange difference
################################################################################
MAX_CBG_LARGE_SPRDIF_GTG_TOPO_DIF           = 0.160
MAX_CBG_LARGE_SPRDIF_GTG_ABS_DIF            = 0.200
MIN_CBG_LARGE_SPRDIF_GTG_ID_RATIO           = 0.700


################################################################################
# Thresholds for accepting CodingBlockGraph into the 
# GeneStructureOfCodingBlocks after HMMsearch on K(s-x) CBGs.
# There are two sets of thresholds: one for normal sized CBGs
# and one for small CBGs. `small` is more or less arbitrary,
# but it is consistently used in fuctions that indeed search
# for small or tiny CBGs (<10AA)
################################################################################
MAX_CBG_HMM_COMPLETION_GTG_TOPO_DIF         = 0.100
MAX_CBG_HMM_COMPLETION_GTG_ABS_DIF          = 0.200
MIN_CBG_HMM_COMPLETION_GTG_ID_RATIO         = 0.800

MAX_SMALL_CBG_HMM_COMPLETION_GTG_TOPO_DIF   = 0.200
MAX_SMALL_CBG_HMM_COMPLETION_GTG_ABS_DIF    = 0.500
MAX_SMALL_CBG_HMM_COMPLETION_GTG_ID_RATIO   = 0.800



################################################################################
# Thresholds for accepting CodingBlockGraph into the GeneStructureOfCodingBlocks
# after splitting on lowconnected spanningrange difference 
################################################################################
MAX_CBG_LOWCONNECTED_SPRDIF_GTG_TOPO_DIF    = 0.275
MAX_CBG_LOWCONNECTED_SPRDIF_GTG_ABS_DIF     = 0.300
MIN_CBG_LOWCONNECTED_SPRDIF_GTG_ID_RATIO    = 0.350
MIN_CBG_LOWCONNECTED_TCODE_OMSR             = 0.780

################################################################################
# Thresholds for removing CodingBlockGraph from the GeneStructureOfCodingBlocks
# from the START of the GSG when there is a better one further up in line 
################################################################################
MAX_CBG_REMOVE_NONSENSE_FIRST_AA_LENGTH     = 8
MAX_CBG_REMOVE_NONSENSE_FIRST_MUTUAL_NODES  = 1
MAX_CBG_REMOVE_NONSENSE_FIRST_GTG_TOPO_DIF  = 0.110  # example: 0.123 in abfgp_mgg0157
MAX_CBG_REMOVE_NONSENSE_FIRST_GTG_ABS_DIF   = 0.150
MIN_CBG_REMOVE_NONSENSE_FIRST_GTG_ID_RATIO  = 0.720
MIN_CBG_REMOVE_NONSENSE_FIRST_TCODE_OMSR    = 0.780
MIN_CBG_REMOVE_NONSENSE_FIRST_CRITERIA_COUNT= 2

################################################################################
# Thresholds for removing CodingBlockGraph from the GeneStructureOfCodingBlocks
# from the END of the GSG when there is a better one further down in line
################################################################################
MAX_CBG_REMOVE_NONSENSE_FINAL_AA_LENGTH     = 8
MAX_CBG_REMOVE_NONSENSE_FINAL_GTG_TOPO_DIF  = 0.110
MAX_CBG_REMOVE_NONSENSE_FINAL_GTG_ABS_DIF   = 0.150
MIN_CBG_REMOVE_NONSENSE_FINAL_GTG_ID_RATIO  = 0.720
MIN_CBG_REMOVE_NONSENSE_FINAL_TCODE_OMSR    = 0.780



################################################################################
# Thresholds for accepting K(s-x) CodingBlockGraph into the GSG
# These thresholds are valid for a GeneTreeGraph with identity == 1.0
# Threshold values are adjusted as a function of genetreegraph's identity
################################################################################
MAX_KSMINX_CBG_GTG_TOPO_DIF                    = 0.065 # adjusted by gtg.identity()
MAX_KSMINX_CBG_GTG_ABS_DIF                     = 0.090 # adjusted by gtg.identity()
MIN_KSMINX_CBG_GTG_ID_RATIO                    = 0.700 # adjusted by gtg.identity()
MAXIMAL_KSMINX_PACBP_GSGOMSR_OVERLAP           = 6     # required in 
                                                       # PacbpCollectionGraph.obtain_ksminxgraphs()
                                                       # was originally 4

################################################################################
# Thresholds for removing K(s-x) CBGs from the GSG
################################################################################
KSMINX_CBG_SINGLE_MISSING_UNLIKELY_EDGE_BITSCORE   = 160.0
KSMINX_CBG_MULTIPLE_MISSING_UNLIKELY_EDGE_BITSCORE = 140.0


################################################################################
# Thresholds for PacbPgraph analyses
# Depends largely of the number of homologous sequences are taken into account
# Number of homologous sequences == EXACT_SG_NODE_COUNT
# A fully connected graph of EXACT_SG_NODE_COUNT has a certain EXACT_SG_EDGE_COUNT
################################################################################
EXACT_SG_NODE_COUNT                 = None  # to be assigned on-the-fly
EXACT_SG_EDGE_COUNT                 = None  # to be assigned on-the-fly
MINIMAL_SG_EDGE_COUNT               = None  # to be assigned on-the-fly
MAXIMAL_SG_MISSING_EDGE_COUNT       = 2     # Threshold when using 5 nodes!!
                                            # When a different number of nodes are used,
                                            # this value is likely to perform suboptimal.
MINIMAL_SG_NODE_CONNECTIVITY        = 2     # Minimal required number of outgoing edges
                                            # for each node in the graph.
                                            # A value of 1 means nothing is removed (each node
                                            # has at least 1 edge, otherwise absent in graph...).
                                            # 2 is a meaningfull threshold when using 5 nodes!!
                                            # When a different number of nodes are used,
                                            # this value is likely to perform suboptimal.
                                            # TODO: make this number dependant on EXACT_SG_NODE_COUNT
                                            # e.g.: [ (4,1), (5,2), (6,2), (7,3), (8,4) ]

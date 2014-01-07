from settings.executables import *
from settings.abgp import ABGP_PROGRAM_NAME, ABGP_VERSION 

# import from all the subclasses/files
from currentannotation import *
from unigene import *
from cbginterface import *

################################################################################
# GFF Feature names for abfgp outcome
################################################################################

GFF_ABFGP_OUTCOME_FSOURCE   = 'abfgp'
GFF_ABFGP_OUTCOME_FMETHOD   = ''                    # defined in subclasses
GFF_ABFGP_OUTCOME_GCLASS    = 'AbfgpCodingBlock'
GFF_ABFGP_OUTCOME_GNAME     = ''                    # auto assigned
GFF_ABFGP_OUTCOME_OUTPUT    = True

GFF_CODINGBLOCK_FSOURCE     = ABGP_VERSION          # use ABGP_VERSION
GFF_CODINGBLOCK_FMETHOD     = 'AbfgpCodingBlock'    # do not change!
GFF_CODINGBLOCK_GCLASS      = 'AbfgpCodingBlock'    # do not change!
GFF_CODINGBLOCK_GNAME       = ''                    # auto assigned
GFF_CODINGBLOCK_OUTPUT      = True                  # do not change!

# track data for AbfgpExon, AbfgpFirstExon, AbfgpLastExon, AbfgpSingleExon
GFF_ABFGP_EXON_FSOURCE      = ABGP_VERSION          # use ABGP_VERSION 
GFF_ABFGP_EXON_FMETHOD      = 'AbfgpExon'           # can be overridden
GFF_ABFGP_EXON_GCLASS       = 'AbfgpExon'           # do not change!
GFF_ABFGP_EXON_GNAME        = ''                    # auto assigned
GFF_ABFGP_EXON_OUTPUT       = True                  # do not change!

# track data for AbfgpIntrons
GFF_ABFGP_INTRON_FSOURCE    = ABGP_VERSION          # use ABGP_VERSION
GFF_ABFGP_INTRON_FMETHOD    = 'AbfgpIntron'         # can be overridden
GFF_ABFGP_INTRON_GCLASS     = 'AbfgpIntron'         # do not change!
GFF_ABFGP_INTRON_GNAME      = ''                    # auto assigned
GFF_ABFGP_INTRON_OUTPUT     = True                  # do not change!

################################################################################
# GFF Feature names for abfgp performance benchmark
################################################################################

GFF_PERFORMANCE_FSOURCE         = ABGP_VERSION              # use ABGP_VERSION 
GFF_PERFORMANCE_FMETHOD         = ''                        # defined in subclasses
GFF_PERFORMANCE_GCLASS          = 'AbfgpPerformance'
GFF_PERFORMANCE_GNAME           = ''                        # auto assigned
GFF_PERFORMANCE_OUTPUT          = True                      # do not change!

# track data for AbfgpCorrectPredicted, TP (True  Positives)
GFF_PERFORMANCE_TP_FSOURCE     = GFF_PERFORMANCE_FSOURCE    # do not change!
GFF_PERFORMANCE_TP_FMETHOD     = 'AbfgpPerformanceTP'       # can be overridden
GFF_PERFORMANCE_TP_GCLASS      = GFF_PERFORMANCE_GCLASS     # do not change!
GFF_PERFORMANCE_TP_GNAME       = ''                         # auto assigned
GFF_PERFORMANCE_TP_OUTPUT      = GFF_PERFORMANCE_OUTPUT     # do not change!

# track data for AbfgpOverPredicted,    FP (False Positives)
GFF_PERFORMANCE_FP_FSOURCE     = GFF_PERFORMANCE_FSOURCE    # do not change!
GFF_PERFORMANCE_FP_FMETHOD     = 'AbfgpPerformanceFP'       # can be overridden
GFF_PERFORMANCE_FP_GCLASS      = GFF_PERFORMANCE_GCLASS     # do not change!
GFF_PERFORMANCE_FP_GNAME       = ''                         # auto assigned
GFF_PERFORMANCE_FP_OUTPUT      = GFF_PERFORMANCE_OUTPUT     # do not change!

# track data for AbfgpUnderPredicted,   FN (False Negatives)
GFF_PERFORMANCE_FN_FSOURCE     = GFF_PERFORMANCE_FSOURCE    # do not change!
GFF_PERFORMANCE_FN_FMETHOD     = 'AbfgpPerformanceFN'       # can be overridden
GFF_PERFORMANCE_FN_GCLASS      = GFF_PERFORMANCE_GCLASS     # do not change!
GFF_PERFORMANCE_FN_GNAME       = ''                         # auto assigned
GFF_PERFORMANCE_FN_OUTPUT      = GFF_PERFORMANCE_OUTPUT     # do not change!


################################################################################
# GFF Feature names for abfgp orf performance benchmark
################################################################################

GFF_ORF_PERFORMANCE_FSOURCE     = ABGP_VERSION                  # use ABGP_VERSION 
GFF_ORF_PERFORMANCE_FMETHOD     = ''                            # defined in subclasses
GFF_ORF_PERFORMANCE_GCLASS      = 'AbgpOrfPerformance'
GFF_ORF_PERFORMANCE_GNAME       = ''                            # auto assigned
GFF_ORF_PERFORMANCE_OUTPUT      = True                          # do not change!

# track data for AbfgpCorrectPredicted, TP (True  Positives)
GFF_ORF_PERFORMANCE_TP_FSOURCE  = GFF_ORF_PERFORMANCE_FSOURCE   # do not change!
GFF_ORF_PERFORMANCE_TP_FMETHOD  = 'AbgpOrfPerformanceTP'        # can be overridden
GFF_ORF_PERFORMANCE_TP_GCLASS   = GFF_ORF_PERFORMANCE_GCLASS    # do not change!
GFF_ORF_PERFORMANCE_TP_GNAME    = ''                            # auto assigned
GFF_ORF_PERFORMANCE_TP_OUTPUT   = GFF_ORF_PERFORMANCE_OUTPUT    # do not change!

# track data for AbfgpOverPredicted,    FP (False Positives)
GFF_ORF_PERFORMANCE_FP_FSOURCE  = GFF_ORF_PERFORMANCE_FSOURCE   # do not change!
GFF_ORF_PERFORMANCE_FP_FMETHOD  = 'AbgpOrfPerformanceFP'        # can be overridden
GFF_ORF_PERFORMANCE_FP_GCLASS   = GFF_ORF_PERFORMANCE_GCLASS    # do not change!
GFF_ORF_PERFORMANCE_FP_GNAME    = ''                            # auto assigned
GFF_ORF_PERFORMANCE_FP_OUTPUT   = GFF_ORF_PERFORMANCE_OUTPUT    # do not change!

# track data for AbfgpUnderPredicted,   FN (False Negatives)
GFF_ORF_PERFORMANCE_FN_FSOURCE  = GFF_ORF_PERFORMANCE_FSOURCE   # do not change!
GFF_ORF_PERFORMANCE_FN_FMETHOD  = 'AbgpOrfPerformanceFN'        # can be overridden
GFF_ORF_PERFORMANCE_FN_GCLASS   = GFF_ORF_PERFORMANCE_GCLASS    # do not change!
GFF_ORF_PERFORMANCE_FN_GNAME    = ''                            # auto assigned
GFF_ORF_PERFORMANCE_FN_OUTPUT   = GFF_ORF_PERFORMANCE_OUTPUT    # do not change!

# track data for AbfgpUnderPredicted,   FN (False Negatives)
GFF_ORF_PERFORMANCE_IGNORED_FSOURCE  = GFF_ORF_PERFORMANCE_FSOURCE   # do not change!
GFF_ORF_PERFORMANCE_IGNORED_FMETHOD  = 'AbgpOrfPerformanceIGNORED'   # can be overridden
GFF_ORF_PERFORMANCE_IGNORED_GCLASS   = GFF_ORF_PERFORMANCE_GCLASS    # do not change!
GFF_ORF_PERFORMANCE_IGNORED_GNAME    = ''                            # auto assigned
GFF_ORF_PERFORMANCE_IGNORED_OUTPUT   = GFF_ORF_PERFORMANCE_OUTPUT    # do not change!

################################################################################
# GFF Feature names for StartCodon & StopCodon estimated quality 
################################################################################
GFF_STARTCODON_OPTIMALITY_FSOURCE = ABGP_VERSION 
GFF_STARTCODON_OPTIMALITY_FMETHOD = 'GSGstartPerf'
GFF_STARTCODON_OPTIMALITY_GCLASS  = 'GSGstartPerf'
GFF_STARTCODON_OPTIMALITY_GNAME   = ''
GFF_STARTCODON_OPTIMALITY_OUTPUT  = True

GFF_STOPCODON_OPTIMALITY_FSOURCE = ABGP_VERSION
GFF_STOPCODON_OPTIMALITY_FMETHOD = 'GSGstopPerf'
GFF_STOPCODON_OPTIMALITY_GCLASS  = 'GSGstopPerf'
GFF_STOPCODON_OPTIMALITY_GNAME   = ''
GFF_STOPCODON_OPTIMALITY_OUTPUT  = True


################################################################################
# GFF Feature names for abfgp introduced/predicted intermediates
################################################################################

# track for GETORF predicted ORFs
GFF_ORF_FSOURCE             = 'GETORF'
GFF_ORF_FMETHOD             = 'orf'
GFF_ORF_GCLASS              = 'Orf'
GFF_ORF_GNAME               = ''                    # auto assigned
GFF_ORF_OUTPUT              = True

# track for PacbP / aligned orfs
GFF_ORFSIMILARITY_FSOURCE   = 'getorf-TBLASTX'
GFF_ORFSIMILARITY_FMETHOD   = 'PacbPORF'
GFF_ORFSIMILARITY_GCLASS    = 'Similarity'
GFF_ORFSIMILARITY_GNAME     = ''                    # auto assigned
GFF_ORFSIMILARITY_OUTPUT    = True

# track for PacbP / aligned UNIGENES 
GFF_UNIGENESIMILARITY_FSOURCE   = 'unigene-TBLASTX'
GFF_UNIGENESIMILARITY_FMETHOD   = GFF_ORFSIMILARITY_FMETHOD
GFF_UNIGENESIMILARITY_GCLASS    = GFF_ORFSIMILARITY_GCLASS
GFF_UNIGENESIMILARITY_GNAME     = GFF_ORFSIMILARITY_GNAME
GFF_UNIGENESIMILARITY_OUTPUT    = True

# track for ProjectedLeadingStop of a PacbPORF
GFF_ORFSIM_PLS_FSOURCE      = GFF_ORFSIMILARITY_FSOURCE # use GFF_ORFSIMILARITY_FSOURCE
GFF_ORFSIM_PLS_FMETHOD      = 'ProjectedLeadingStop'    # free naming
GFF_ORFSIM_PLS_GCLASS       = GFF_ORFSIMILARITY_GCLASS  # use GFF_ORFSIMILARITY_GCLASS
GFF_ORFSIM_PLS_GNAME        = ''                        # auto assigned
GFF_ORFSIM_PLS_OUTPUT       = True

# track for ProjectedTailingStop of a PacbPORF
GFF_ORFSIM_PTS_FSOURCE      = GFF_ORFSIMILARITY_FSOURCE # use GFF_ORFSIMILARITY_FSOURCE
GFF_ORFSIM_PTS_FMETHOD      = 'ProjectedTailingStop'    # free naming
GFF_ORFSIM_PTS_GCLASS       = GFF_ORFSIMILARITY_GCLASS  # use GFF_ORFSIMILARITY_GCLASS
GFF_ORFSIM_PTS_GNAME        = ''                        # auto assigned
GFF_ORFSIM_PTS_OUTPUT       = True

# track for PacbpORFs of ReversecomplementCodingBlockGraphs
GFF_REVCBG_SIMILARITY_FSOURCE   = 'revCBG-getorf-TBLASTX'
GFF_REVCBG_SIMILARITY_FMETHOD   = 'PacbPORF'
GFF_REVCBG_SIMILARITY_GCLASS    = 'Similarity'
GFF_REVCBG_SIMILARITY_GNAME     = ''                    # auto assigned
GFF_REVCBG_SIMILARITY_OUTPUT    = True

# track for Transmembranehelices
GFF_TMHMM_FSOURCE             = 'TMHMM2.0'
GFF_TMHMM_FMETHOD             = 'TMhelix'
GFF_TMHMM_GCLASS              = ''                    # auto assigned
GFF_TMHMM_GNAME               = ''                    # auto assigned
GFF_TMHMM_OUTPUT              = True

# track for SignalPeptides
GFF_SIGNALP_FSOURCE             = 'SignalP3.0'
GFF_SIGNALP_FMETHOD             = 'predSignalPeptide'
GFF_SIGNALP_GCLASS              = 'SignalPSignalPeptide'  # auto assigned
GFF_SIGNALP_GNAME               = ''                      # auto assigned
GFF_SIGNALP_OUTPUT              = True

# track for TranslationalStartSite
GFF_TSS_FSOURCE             = 'tssPSSM4FUNGI'
GFF_TSS_FMETHOD             = 'tsspssm'
GFF_TSS_GCLASS              = 'TranslationalStartSite'  
GFF_TSS_GNAME               = ''                    # auto assigned
GFF_TSS_OUTPUT              = True

# track for AlignedTranslationalStartSite
GFF_ALIGNED_TSS_FSOURCE     = 'tssPSSM4FUNGI-PACBP'
GFF_ALIGNED_TSS_FMETHOD     = 'aligned_tss'
GFF_ALIGNED_TSS_GCLASS      = 'AlignedTranslationalStartSite'
GFF_ALIGNED_TSS_GNAME       = ''                    # auto assigned
GFF_ALIGNED_TSS_OUTPUT      = True
GFF_ALIGNED_TSS_REPORTBEST  = 2

# track for AbfgpMultipleAlignment
GFF_AMA_FSOURCE             = ABGP_VERSION          # use ABGP_VERSION 
GFF_AMA_FMETHOD             = 'AbfgpMultipleAlignment'
GFF_AMA_GCLASS              = 'AMA'
GFF_AMA_GNAME               = ''                    # auto assigned
GFF_AMA_OUTPUT              = True

# track for AbfgpMultipleOccurrence
GFF_AMO_FSOURCE             = ABGP_VERSION          # use ABGP_VERSION 
GFF_AMO_FMETHOD             = 'AbfgpMultipleOccurrence'
GFF_AMO_GCLASS              = 'AMO'
GFF_AMO_GNAME               = ''                    # auto assigned
GFF_AMO_OUTPUT              = True

# track for AbfgpMultipleCombination
GFF_AMC_FSOURCE             = ABGP_VERSION          # use ABGP_VERSION 
GFF_AMC_FMETHOD             = 'AbfgpMultipleCombination'
GFF_AMC_GCLASS              = 'AMC'
GFF_AMC_GNAME               = ''                    # auto assigned
GFF_AMC_OUTPUT              = True


# track for tcode tracks
GFF_TCODE_FSOURCE             = "%s-step%s-window%s" % (
                                    EXECUTABLE_TCODE_VERSION,
                                    EXECUTABLE_TCODE_STEP,
                                    EXECUTABLE_TCODE_WINDOW,
                                    )
GFF_TCODE_FMETHOD             = 'tcode'
GFF_TCODE_GCLASS              = 'TcodeWindowCenter'
GFF_TCODE_GNAME               = ''                    # auto assigned
GFF_TCODE_OUTPUT              = True


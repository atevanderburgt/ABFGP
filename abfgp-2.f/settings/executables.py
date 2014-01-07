"""
Executables and command line arguments for 3th party software used in ABGP

Valid at dev1.ab.wurnet.nl 
"""

from os.path import join as osPathJoin
from settings.abgp import MAIN_ABGP_PATH as BASEPATH

# python
PYTHON_PATH             = "/usr/bin/python2.6" # tested are 2.4, 2.6 and 2.7.4
PYTHON_VERSION          = "python2.6.2"
# emboss
EMBOSS_EXECUTABLES_PATH = "/usr/bin"
EMBOSS_DATA_DIRECTORY   = "/usr/share/EMBOSS/data"
EMBOSS_VERSION          = "4.0.0"
# perl
PERL_PATH               = "perl"    # only required in case a MySQL GGB database is
                                    # coupled to store and visualize ABFGP output.
PERL_VERSION            = "?"       # not used at all.

################################################################################
# Full paths to executables
################################################################################
EXECUTABLE_BLASTALL     = osPathJoin(BASEPATH,"software/blast-2.2.8/blastall")
EXECUTABLE_FORMATDB     = osPathJoin(BASEPATH,"software/blast-2.2.8/formatdb")
EXECUTABLE_GETORF       = osPathJoin(EMBOSS_EXECUTABLES_PATH,"getorf")
EXECUTABLE_TCODE        = osPathJoin(EMBOSS_EXECUTABLES_PATH,"tcode")
EXECUTABLE_TRANSEQ      = osPathJoin(EMBOSS_EXECUTABLES_PATH,"transeq")
EXECUTABLE_CLUSTALW     = osPathJoin(BASEPATH,"software/clustalw-1.83/clustalw")
EXECUTABLE_SIGNALP      = osPathJoin(BASEPATH,"software/signalp-3.0/signalp")
EXECUTABLE_TMHMM        = osPathJoin(BASEPATH,"software/tmhmm-2.0/bin/tmhmm")
EXECUTABLE_SFM          = osPathJoin(BASEPATH,"software/ScanForMatches/scan_for_matches")
EXECUTABLE_HMMPATH      = osPathJoin(BASEPATH,"software/hmmer-2.3.2")
EXECUTABLE_HMMSEARCH    = osPathJoin(EXECUTABLE_HMMPATH,"hmmsearch")
EXECUTABLE_HMMBUILD     = osPathJoin(EXECUTABLE_HMMPATH,"hmmbuild")

EXECUTABLE_CEXPANDER_PATH       = osPathJoin(BASEPATH,"software/cexpander-1.0")
EXECUTABLE_CEXPANDER_ALLVSALL   = osPathJoin(EXECUTABLE_CEXPANDER_PATH,
                                             "prep_launch.py")
EXECUTABLE_CEXPANDER_CBALIGNP   = osPathJoin(EXECUTABLE_CEXPANDER_PATH,
                                             "cbalignp_test_2/src/cbalignp")
EXECUTABLE_CEXPANDER_CEXPANDER  = osPathJoin(EXECUTABLE_CEXPANDER_PATH,
                                             "cexpander_dr")

EXECUTABLE_GFF2FASTA             = osPathJoin(BASEPATH,"software/get_gff_sequence_from_fasta.py")
EXECUTABLE_UNIGENEANNOTATION     = osPathJoin(BASEPATH,"unigeneannotation.py")


EXECUTABE_SCAN_BRANCHPOINT       = osPathJoin(BASEPATH,"software/gsregex/scan_branchpoint.sh")
EXECUTABE_SCAN_POLYPYRIMIDINES   = osPathJoin(BASEPATH,"software/gsregex/scan_polypyrimidines.py")
EXECUTABLE_FASTALENGTH           = osPathJoin(BASEPATH,"software/fastalength.sh")
EXECUTABLE_GFFLENGTH             = osPathJoin(BASEPATH,"software/gfflength.sh")
EXECUTABLE_LOAD_GFF              = osPathJoin(BASEPATH,"software/load_gff.pl")


################################################################################
# Versions of executables - used for GffFeature fsource and meta-data
################################################################################
EXECUTABLE_BLAST_VERSION        = "blast-2.2.8"
EXECUTABLE_BLASTALL_VERSION     = EXECUTABLE_BLAST_VERSION
EXECUTABLE_FORMATDB_VERSION     = EXECUTABLE_BLAST_VERSION
EXECUTABLE_GETORF_VERSION       = "getorf-"+EMBOSS_VERSION
EXECUTABLE_TCODE_VERSION        = "tcode-"+EMBOSS_VERSION
EXECUTABLE_TRANSEQ_VERSION      = "transeq-"+EMBOSS_VERSION
EXECUTABLE_CLUSTALW_VERSION     = "clustalw-2.0.9"
EXECUTABLE_SIGNALP_VERSION      = "SignalP3.0"
EXECUTABLE_TMHMM_VERSION        = "TMHMM2.0"
EXECUTABLE_HMM_VERSION          = "hmm-2.3.2"
EXECUTABLE_HMMSEARCH_VERSION    = EXECUTABLE_HMM_VERSION
EXECUTABLE_HMMBUILD_VERSION     = EXECUTABLE_HMM_VERSION
EXECUTABLE_CEXPANDER_VERSION    = "cexpander-1.0"

################################################################################
# tcode command line parameters
################################################################################
EXECUTABLE_TCODE_STEP           = 1     # >=200, see `tcode -help`
EXECUTABLE_TCODE_WINDOW         = 200   # >=200, see `tcode -help`
EXECUTABLE_TCODE_GFFOUTPUT_STEP = 9     # >=EXECUTABLE_TCODE_STEP
EXECUTABLE_TCODE_N_CORRECTION   = 100   # <<200 nt, long N-tracks cause errors
TCODE_MAX_NONCODING             = 0.740 # <= TCODE calls Non-coding           
TCODE_MIN_CODING                = 0.950 # >= TCODE calls Coding

################################################################################
# SignalP thresholds
################################################################################
SIGNALP_MIN_D                   = 0.430 # minimal SignalP D value; column 13
SIGNALP_MIN_SPROB               = 0.500 # minimal SignalP Sprob value; column 20

################################################################################
# getorf command line parameters
################################################################################
EXECUTABLE_GETORF_MINSIZE       = 9     # minimum orf nt length; was 15 nt
                                        # BC1G_15776 has tiny final ORF: ..agNTGA
                                        # of in total 9nt/3AA length *AAA*



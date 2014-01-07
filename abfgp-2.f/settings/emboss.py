"""
Alignment Based Gene Predictions settings: EMBOSS data files 
"""

from os.path import join as osPathJoin
from abgp import MAIN_ABGP_PATH

################################################################################
# Full paths to Amino Acid similarity matrices
################################################################################
BLOSUM62_PATH   = osPathJoin(MAIN_ABGP_PATH,"datafiles/EMBOSS/EBLOSUM62")
BLOSUM80_PATH   = osPathJoin(MAIN_ABGP_PATH,"datafiles/EMBOSS/EBLOSUM80")
BLOSUM45_PATH   = osPathJoin(MAIN_ABGP_PATH,"datafiles/EMBOSS/EBLOSUM45")
PAM30_PATH      = osPathJoin(MAIN_ABGP_PATH,"datafiles/EMBOSS/EPAM30")
PAM70_PATH      = osPathJoin(MAIN_ABGP_PATH,"datafiles/EMBOSS/EPAM70")


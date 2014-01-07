"""
Alignment Based Gene Predictions settings: all paths and settings required
for succesfull execution of the ABGP suite of algorithms
"""

class ParameterOptions:
    """ Empty class to assemble related settings """
    def __init__(self): pass
# end of class ParameterOptions

from abgp import *
from dbwarehouse import *
from genetreegraph import *
from genestructure import *
from codingblockgraph import *
from executables import *
from splicesites import *
from translationalstartsites import *
from sitealignment import *
from alignedstopcodongraph import *
from pacbp import *
from blastp import *
from inframeintron import *
from gff import *
from ggb import *

### ETCETERA; SOME THINGIES THAT MUST BE PLACES SOMEWHERE ELSE OR ARE ABOUT TO BE DEPRECATED!


################################################################################            # DEPRECATED; NOT USED (ANYMORE)
# Score values to assign to identity/similarity on an aligned amino acid                    # DEPRECATED; NOT USED (ANYMORE)
# this is used to construct the data structure that can output the VISTA-like track         # DEPRECATED; NOT USED (ANYMORE)
# TODO: implement these global variables in graphAbgp                                       # DEPRECATED; NOT USED (ANYMORE)
# TODO: use similarity maxtrix in stead of (semi)binary matrix                              # DEPRECATED; NOT USED (ANYMORE)
################################################################################            # DEPRECATED; NOT USED (ANYMORE)
POSITIONAL_ALIGNMENT_SIMILARITY_SCORES = ( 1.0, 0.5, -1.0 ) # identity_score   = 1.0        # DEPRECATED; NOT USED (ANYMORE)
                                                            # similarity_score = 0.5        # DEPRECATED; NOT USED (ANYMORE)
                                                            # gap_score        = -1.0       # DEPRECATED; NOT USED (ANYMORE)

################################################################################            # DEPRECATED; NOT USED (ANYMORE)
# Thresholds for (expected) flanking sequence length.                                       # DEPRECATED; NOT USED (ANYMORE)
# TODO: this dependancy should be discarded in the future!                                  # DEPRECATED; NOT USED (ANYMORE)
# It is used to make sure the central gene of interest is evaluated                         # DEPRECATED; NOT USED (ANYMORE)
# in case of a region of microsyntheny.                                                     # DEPRECATED; NOT USED (ANYMORE)
################################################################################            # DEPRECATED; NOT USED (ANYMORE)
FLANKING_SEQUENCE_NT_LENGTH         = 1500                                                  # DEPRECATED; NOT USED (ANYMORE)
FLANKING_SEQUENCE_AA_LENGTH         = FLANKING_SEQUENCE_NT_LENGTH/3                         # DEPRECATED; NOT USED (ANYMORE)
START_SITE_OVERLAP_SEARCH_OFFSET    = 200                                                   # DEPRECATED; NOT USED (ANYMORE)


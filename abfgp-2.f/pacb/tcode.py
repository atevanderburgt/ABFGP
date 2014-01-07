"""
TCODE accessibility functions of PacbpORFs
"""

# Module metadata
__authors__ = "Ate van der Burgt"
__license__ = "MIT"


def tcode_query(pacbporf):
    """ """
    aa_start= pacbporf._get_original_alignment_pos_start().query_pos
    aa_end  = pacbporf._get_original_alignment_pos_end().query_pos + 1
    return pacbporf.orfQ.tcode_score(aa_start,aa_end)

# end of function tcode_query


def tcode_sbjct(pacbporf):
    """ """
    aa_start= pacbporf._get_original_alignment_pos_start().sbjct_pos
    aa_end  = pacbporf._get_original_alignment_pos_end().sbjct_pos + 1
    return pacbporf.orfS.tcode_score(aa_start,aa_end)

# end of function tcode_sbjct

def tcode(pacbporf):
    """ """
    tcQ,tcS = tcode_query(pacbporf),tcode_sbjct(pacbporf)
    if tcQ==0.0 and tcS==0.0: return 0.0
    elif tcQ==0.0:            return tcS
    elif tcS==0.0:            return tcQ
    else:                     return (tcQ+tcS)/2

# end of function tcode



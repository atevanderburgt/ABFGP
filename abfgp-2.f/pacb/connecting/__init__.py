"""
Functions for connecting PacPORFs into a genemodel 
"""
# Module metadata
__authors__ = "Ate van der Burgt"
__license__ = "MIT"


import connecting

from mapping import (
    merge_pacbporfs_with_introns,
    merge_pacbporfs_with_closeby_independant_introns,
    merge_pacbporfs_with_phase_shift_introns,
    )


from projecting import (
    merge_pacbporfs_by_intron_tinyexon_intron_in_query,
    merge_pacbporfs_by_intron_tinyexon_intron_in_sbjct,
    merge_pacbporfs_by_intron_in_query,
    merge_pacbporfs_by_intron_in_sbjct,
    merge_pacbporfs_by_inframe_intron_in_query,
    merge_pacbporfs_by_inframe_intron_in_sbjct,
    # DEPRECATED FUNCTION NAMES
    project_splicesites_on_pacbporfs_with_lacking_intron_in_sbjct,
    project_splicesites_on_pacbporfs_with_lacking_intron_in_query,
    )



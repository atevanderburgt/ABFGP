"""
PacbPORF connection by sequence error
"""

# Module metadata
__authors__ = "Ate van der Burgt"
__license__ = "MIT"

# Imports
from pacb.validators import IsPacbPORF, IsIdenticalOrfs

# Other Imports
from gene.sequenceerror import (
    SequenceErrorConnectingOrfs,
    merge_orfs_with_sequenceerror,
    )

# Python Imports
from sets import Set

# Global variable Imports


def merge_pacbporf_with_sequenceerror_in_query(pacbporfD,pacbporfA,verbose=False,**kwargs):
    """
    Merge 2 PacbPORFs with a sequence error

    @type  pacbporfD: PacbPORF object
    @param pacbporfD: PacbPORF object that has to deliver PSSM donor objects

    @type  pacbporfA: PacbPORF object
    @param pacbporfA: PacbPORF object that has to deliver PSSM acceptor objects

    @type  verbose: Boolean
    @param verbose: print status/debugging messages to STDOUT

    @rtype:  list
    @return: list with ( intron, intron ), in query and sbjct
    """

    # at least 6, preferably somewhat larger (but not very large!)
    MAX_SEQUENCE_ERROR_PACBP_CONNECTION_DISTANCES = (
        ( 21, 1 ), # 7AA offset, maximal 1nt difference (insert/deletion)
        (  7, 2 ), # 2AA offset, maximal 2nt difference due to 1 AA loss
        )


    # input validation
    IsPacbPORF(pacbporfD)
    IsPacbPORF(pacbporfA)

    # Orfs of SBJCT must be identical
    IsIdenticalOrfs(pacbporfD.orfS,pacbporfA.orfS)

    # Orfs of SBJCT must be non-identical
    if pacbporfD.orfQ.id == pacbporfA.orfQ.id: return False

    # check if PacbPORFs are (nearly) seamlessly connected
    qDend = pacbporfD._get_original_alignment_pos_end().query_dna_end
    qAsta = pacbporfA._get_original_alignment_pos_start().query_dna_start
    sDend = pacbporfD._get_original_alignment_pos_end().sbjct_dna_end
    sAsta = pacbporfA._get_original_alignment_pos_start().sbjct_dna_start

    # calculate distance between aligments
    distQ = max([ 0 , qAsta - qDend ])
    distS = max([ 0 , sAsta - sDend ])

    ############################################################################
    if verbose:
        orfids = (pacbporfD.orfQ.id,pacbporfA.orfQ.id,pacbporfD.orfS.id)
        print "SEtest::",distQ,distS, orfids
    ############################################################################

    # small distances && identical length? -> check for potential sequence error
    is_accepted_sequence_error_offset = False
    for (max_pacbp_nt_distance, max_nt_gap) in\
    MAX_SEQUENCE_ERROR_PACBP_CONNECTION_DISTANCES:
        if abs(distQ-distS) <= max_nt_gap and\
        max([distQ,distS])  <= max_pacbp_nt_distance:
            is_accepted_sequence_error_offset = True
            break

    if is_accepted_sequence_error_offset:
        # get SequenceErrorConnectingOrfs objects & store metadata
        se = merge_orfs_with_sequenceerror(pacbporfD.orfQ,pacbporfA.orfQ)
        if se:
            se._linked_to_introns   = []
            se._linked_to_pacbporfs = []
            # _distance in ABFGP introns es expressed in NTs.
            # Because the (putative) sequence errors cannot be aligned itself,
            # correct with 3nt and correct to a minimum of 0
            se._distance            = ( max([ 0, (max([distQ,distS])-3) ]) / 3 ) + abs(distQ-distS)
            se._apps_donor          = 1.0
            se._apps_accep          = 1.0
            se._gff['fsource']      = "AbfgpSequenceError"
            se.pssm_score           = 1
            # when the SequenceError is a Deletion/Insertion,
            # overwrite se._estimated_nt_coord
            if se.typeof in ['insertion','deletion']:
                se._estimated_nt_coord = qDend+(distQ/2)
                # make shure estimated nt coord is in the proper frame!
                while qAsta % 3 != se._estimated_nt_coord % 3:
                    se._estimated_nt_coord += 1 
            ####################################################################
            if verbose: print "SEQUENCE-ERROR",distQ,distS,se
            if verbose: print pacbporfD, "\n", pacbporfA
            ####################################################################
            return se
        else:
            # not an obvious sequence error
            return False
    else:
        return False

# end of function merge_pacbporf_with_sequenceerror_in_query


def merge_pacbporf_with_sequenceerror_in_sbjct():
    """ todo.... """
    pass

# end of function merge_pacbporf_with_sequenceerror_in_sbjct

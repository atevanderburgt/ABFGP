"""
Overlap functions for python-pacb.
"""

# Module metadata
__authors__ = "Ate van der Burgt"
__license__ = "MIT"


from exceptions import NoPacbP, InproperlyAppliedArgument, ZeroSizedPacb
from comparison import (
    IsIdenticalPacbPORF,
    IsPacbPORFSubstringOf,
    )



def correct_overlap_for_sbjct(prevPACBP,nextPACBP,**kwargs):
    """ """
    return _correct_overlap("sbjct",prevPACBP,nextPACBP,**kwargs)
# end of function correct_overlap_for_sbjct


def correct_overlap_for_query(prevPACBP,nextPACBP,**kwargs):
    """ """
    return _correct_overlap("query",prevPACBP,nextPACBP,**kwargs)
# end of function correct_overlap_for_query


def _correct_overlap(sbjctorquery,prevPACBP,nextPACBP,max_aa_overlap=2,
    verbose=False):
    """
    """
    # minimal value for overlap is 2 AAs -> otherwise no overlap ;-)
    if max_aa_overlap < 2: max_aa_overlap = 2

    # check sbjctorquery argument
    if sbjctorquery not in ['sbjct','query']:
        message = "'sbjctorquery' is '%s', not 'sbjct' or 'query'" % sbjctorquery
        raise InproperlyAppliedArgument, message

    # check if pacbps are PacbP or PacbPDNA, PacbPORF
    prevObjType = prevPACBP.__class__.__name__
    nextObjType = nextPACBP.__class__.__name__
    if prevObjType not in ['PacbP','PacbPDNA','PacbPORF']:
        raise NoPacbP, prevPACBP
    if nextObjType not in ['PacbP','PacbPDNA','PacbPORF']:
        raise NoPacbP, nextPACBP


    # check if PacbP(ORFs) are zero length (due to previous truncation)
    if len(prevPACBP) == 0 or len(nextPACBP) == 0:
        return ( prevPACBP, nextPACBP, None )

    if IsIdenticalPacbPORF(prevPACBP,nextPACBP):
        # identical objects -> quit!
        return ( prevPACBP, nextPACBP, None )

    if IsPacbPORFSubstringOf(prevPACBP,nextPACBP):
        # included PacbPs !?!? please, exit here, this can ge horribly wrong!
        return ( prevPACBP, nextPACBP, None )

    if sbjctorquery == 'sbjct':
        # get all coordinated from the sbjct in the alignment objects
        if prevObjType == 'PacbP':
            # coordinate minus 1 -> sbjct_end isa list slice coordinate
            end      = prevPACBP.sbjct_end - 1
            prev_sta = prevPACBP.sbjct_start
        else:
            end      = prevPACBP._get_original_alignment_pos_end().sbjct_pos
            prev_sta = prevPACBP._get_original_alignment_pos_start().sbjct_pos
        if nextObjType == 'PacbP':
            sta      = nextPACBP.sbjct_start
            next_end = nextPACBP.sbjct_end
        else:
            sta      = nextPACBP._get_original_alignment_pos_start().sbjct_pos
            # coordinate plus 1 -> correct to a list slice coordinate
            next_end = nextPACBP._get_original_alignment_pos_end().sbjct_pos+1
    else:
        # get all coordinated from the query in the alignment objects
        if prevObjType == 'PacbP':
            # coordinate minus 1 -> query_end isa list slice coordinate
            end      = prevPACBP.query_end - 1
            prev_sta = prevPACBP.query_start
        else:
            end      = prevPACBP._get_original_alignment_pos_end().query_pos
            prev_sta = prevPACBP._get_original_alignment_pos_start().query_pos
        if nextObjType == 'PacbP':
            sta      = nextPACBP.query_start
            next_end = nextPACBP.query_end
        else:
            sta      = nextPACBP._get_original_alignment_pos_start().query_pos
            # coordinate plus 1 -> correct to a list slice coordinate
            next_end = nextPACBP._get_original_alignment_pos_end().query_pos+1

    if sta-end >= -max_aa_overlap:
        # no overlap -> return input objects
        return ( prevPACBP, nextPACBP, None )

    # if here, pacbps are overlapping
    overlap_length = abs(sta-end)
    # list arrays to store calculated alignment bitscore differences into
    score_array_prev = []
    score_array_next = []

    # calculate bitscore of NON-overlapping fraction of prevPACBP 
    if sbjctorquery == 'sbjct':
        prev_iter_prev_bit = prevPACBP.bitscore_slice_by_abs_protein_sbjct(
                prev_sta,sta)
    else:
        prev_iter_prev_bit = prevPACBP.bitscore_slice_by_abs_protein_query(
                prev_sta,sta)
    # get bitscore of COMPLETE nextPACBP
    prev_iter_next_bit = nextPACBP.bitscore

    ############################################################################
    if verbose:
        print prevPACBP, "OVERLAP", overlap_length
        print nextPACBP, "OVERLAP", overlap_length
    ############################################################################

    # now loop over the overlapping region.
    # start with NON-overlapping fraction of prevPACBP and complete prevPACBP
    # by iterating over the offset
    # prevPACBP is extended  by 1 AA -> calc bitscore -> calc dif with prev iter
    # nextPACBP is shortened by 1 AA -> calc bitscore -> calc dif with prev iter
    # the differences are stored into score_array_prev and score_array_next
    for offset in range(1,overlap_length):
        if sbjctorquery == 'sbjct':
            prevBit = prevPACBP.bitscore_slice_by_abs_protein_sbjct(
                        prev_sta,sta+offset)
            nextBit = nextPACBP.bitscore_slice_by_abs_protein_sbjct(
                        sta+offset,next_end)
        else:
            prevBit = prevPACBP.bitscore_slice_by_abs_protein_query(
                        prev_sta,sta+offset)
            nextBit = nextPACBP.bitscore_slice_by_abs_protein_query(
                        sta+offset,next_end)

        ########################################################################
        if verbose:
            print offset, 
            if prevObjType != 'PacbP':
                if sbjctorquery == 'sbjct':
                    print prevPACBP.alignmentobject_by_sbjctaa(sta+offset-1).sbjct,
                else:
                    print prevPACBP.alignmentobject_by_queryaa(sta+offset-1).query,
            print prevBit-prev_iter_prev_bit, prev_iter_next_bit-nextBit
        ########################################################################


        # append bitscore differences to score arrays
        score_array_prev.append( prevBit - prev_iter_prev_bit )
        score_array_next.append( prev_iter_next_bit - nextBit )

        # update alignemnt bitscore for next iteration
        prev_iter_prev_bit = prevBit
        prev_iter_next_bit = nextBit

    # now, loop over the same interface, -1 corrected
    # correction is because overlap_length -1 == len(score_array_xxxx)
    summed_score = []
    for pos in range(1,overlap_length-1):
        bit_dif_prev_pacbp = sum(score_array_prev[0:pos])
        bit_dif_next_pacbp = sum(score_array_next[pos:])
        summed_score.append( bit_dif_prev_pacbp + bit_dif_next_pacbp )
        ########################################################################
        if verbose:
            print pos, ":", bit_dif_prev_pacbp,bit_dif_next_pacbp,
            print summed_score[-1]
        ########################################################################

    # now find the position where the summed score is optimal,
    # which corresponds to the position where to correct the overlap
    max_score = max(summed_score)

    # find position to correct the nextPACBP for
    cor_offset_next = summed_score.index(max_score)

    # find position to correct the prevPACBP for
    if summed_score.count(max_score) > 1:
        # take the most 3' occurring position; max_score is observed > 1x
        offset = 0
        while True:
            try:
                cor_offset_prev = summed_score.index(max_score,offset)
                offset = cor_offset_prev+1
            except ValueError:
                break
    else:
        cor_offset_prev = summed_score.index(max_score)

    # correct prevPACBP
    if prevObjType == 'PacbPORF':
        if prevPACBP.is_extended():
            prevPACBP.unextend_pacbporf()
            was_extended = True
        else:
            was_extended = False
        try:
            prevPACBP.pop_alignment_object(aas=overlap_length-cor_offset_prev)
            if prevPACBP._positions and\
            prevPACBP._positions[-1].match in list(" .+:*"):
                prevPACBP.pop_alignment_object()
            prevPACBP.strip_unmatched_ends()
            # get corrected end coordinate
            if sbjctorquery == 'sbjct':
                corrected_end = prevPACBP.sbjct_end - 1
            else:
                corrected_end = prevPACBP.query_end - 1
            if was_extended:
                prevPACBP.extend_pacbporf_after_stops()
        except ZeroSizedPacb:
            # PacbP is completely stripped (zerosized)
            corrected_end = 0
    else:
        try:
            for iter in range(0,(overlap_length-cor_offset_prev)):
                prevPACBP._strip_one_rigth_position()
            if prevPACBP.match and prevPACBP.match[-1] in list(" .+:*"):
                prevPACBP._strip_one_rigth_position()
            prevPACBP.strip_unmatched_ends()
            # get corrected end coordinate
            if sbjctorquery == 'sbjct':
                corrected_end = prevPACBP.sbjct_end - 1
            else:
                corrected_end = prevPACBP.query_end - 1
        except ZeroSizedPacb:
            # PacbP is completely stripped (zerosized)
            corrected_end = 0

    # correct nextPACBP
    if nextObjType == 'PacbPORF':
        try:
            if nextPACBP.is_extended():
                nextPACBP.unextend_pacbporf()
                was_extended = True
            else:
                was_extended = False
            nextPACBP.shift_alignment_object(cor_offset_next)
            if nextPACBP._positions and\
            nextPACBP._positions[0].match in list(" .+:*"):
                nextPACBP.shift_alignment_object()
            nextPACBP.strip_unmatched_ends()
            # get corrected start coordinate
            if sbjctorquery == 'sbjct':
                corrected_start = nextPACBP.sbjct_start
            else:
                corrected_start = nextPACBP.query_start
            if was_extended:
                nextPACBP.extend_pacbporf_after_stops()
        except ZeroSizedPacb:
            # PacbP is completely stripped (zerosized)
            corrected_end = 0
    else:
        try:
            for iter in range(0,cor_offset_next):
                nextPACBP._strip_one_left_position()
            if nextPACBP.match and nextPACBP.match[0] in list(" .+:*"):
                prevPACBP._strip_one_left_position()
            nextPACBP.strip_unmatched_ends()
            # get corrected start coordinate
            if sbjctorquery == 'sbjct':
                corrected_start = nextPACBP.sbjct_start
            else:
                corrected_start = nextPACBP.query_start
        except ZeroSizedPacb:
            # PacbP is completely stripped (zerosized)
            corrected_end = 0

    ############################################################################
    if verbose:
        new_overlap_length = abs(corrected_start - corrected_end)
        print prevPACBP, "DONE", overlap_length, "->", new_overlap_length
        print nextPACBP, "DONE", overlap_length, "->", new_overlap_length
        #q,m,s = prevPACBP.get_unextended_aligned_protein_sequences()
        #print q
        #print m
        #print s
        #print "\tOVERLAP"
        #print nextPACBP
        #q,m,s = nextPACBP.get_unextended_aligned_protein_sequences()
        #print q
        #print m
        #print s
        #print ""
    ############################################################################

    return prevPACBP, nextPACBP, True

# end of function _correct_overlap

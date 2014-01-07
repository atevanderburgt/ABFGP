"""
"""
# Module metadata
__authors__ = "Ate van der Burgt"
__license__ = "MIT"

# Imports
from pacb.ordering import order_pacbporf_list

# Make sure Numeric is installed somewhere!!!
# For convenience, it is installed in ./requiredmodules
#import sys, os
#from os.path import join as osPathJoin
#from settings.abgp import MAIN_ABGP_PATH as BASEPATH
#sys.path.append(osPathJoin(BASEPATH,"requiredmodules"))
#from Numeric import zeros, where, greater, greater_equal

import sys,os
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)),"requiredmodules"))
from Numeric import zeros, where, greater, greater_equal


# Global variable Imports
from settings.genestructure import ORF_IS_UNIGENE_LABEL

def PCG2codingarray(PCG,organism,aalength,omit_unigenes=True):
    """ """
    array_algpresence = zeros(aalength)
    for orgS in PCG.organism_set():
        if organism == orgS: continue
        pacbporfs = order_pacbporf_list(PCG.get_pacbps_by_organisms(organism,orgS))
        if pacbporfs and omit_unigenes and hasattr(pacbporfs[0].orfS,ORF_IS_UNIGENE_LABEL):
            continue
        orgPresArray = pacbporflist2codingarray(pacbporfs,"query",aalength)
        array_algpresence+=orgPresArray

    # return coding/presence array
    return array_algpresence

# end of function PCG2codingarray


def PCG2similarityarray(PCG,organism,aalength,omit_unigenes=True):
    """ """
    array_algsimilarity = zeros(aalength)
    for orgS in PCG.organism_set():
        if organism == orgS: continue
        pacbporfs = order_pacbporf_list(PCG.get_pacbps_by_organisms(organism,orgS))
        if pacbporfs and omit_unigenes and hasattr(pacbporfs[0].orfS,ORF_IS_UNIGENE_LABEL):
            continue
        orgSimArray = pacbporflist2similarityarray(pacbporfs,"query",aalength)
        array_algsimilarity+=orgSimArray

    # return similarity array
    return array_algsimilarity

# end of function PCG2similarityarray


def pacbporflist2codingarray(pacbps,queryorsbjct,length):
    """ """
    bea = zeros(length)
    for pacbporf in pacbps:
        spos = pacbporf._get_original_alignment_pos_start()
        epos = pacbporf._get_original_alignment_pos_end()
        if queryorsbjct == 'query':
            bea[spos.query_pos:epos.query_pos+1] = 1
        else:
            bea[spos.sbjct_pos:epos.sbjct_pos+1] = 1
    return bea
# end of function pacbporflist2codingarray


def pacbporflist2similarityarray(pacbps,queryorsbjct,length):
    """ """
    bea = zeros(length)
    for pacbporf in pacbps:
        spos  = pacbporf._get_original_alignment_pos_start()
        epos  = pacbporf._get_original_alignment_pos_end()
        q,m,s = pacbporf.get_unextended_aligned_protein_sequences()
        if queryorsbjct == 'query':
            start = spos.query_pos
            end   = epos.query_pos + 1
            seqa  = list(q)
            ma    = list(m)
        else:
            start = spos.sbjct_pos
            end   = epos.sbjct_pos + 1
            seqa  = list(s)
            ma    = list(m)

        for pos in range(len(seqa)-1,-1,-1):
            if seqa[pos] == '-':
                seqa.pop(pos)
                ma.pop(pos)
        # prepare replacement of match string into match score list
        matcharray = zeros(end-start)
        for pos in range(0,len(ma)):
            symbol = ma[pos]
            if symbol != ' ':
                matcharray[pos] = 1
        # update (binary) array
        bea[spos.query_pos:epos.query_pos+1] += matcharray

    # correct bea for values > 1
    bea = where(greater_equal(bea, 2), 1, bea)
    return bea

# end of function pacbporflist2similarityarray


def get_ordered_array_window_scores(arrayobj,window_aa_size):
    """ """ 
    array_windows = []
    offset_start = 0
    offset_end   = len(arrayobj)-1
    while arrayobj[offset_start] == 0: offset_start+=1
    while arrayobj[offset_end] == 0:   offset_end-=1
    offset_end+=1

    for offset in range(offset_start,offset_end):
        window = arrayobj[offset:offset+window_aa_size]
        total = sum(window)
        if total != 0: array_windows.append(total)
    array_windows.sort()
    array_windows.reverse()
    return array_windows

# end of function get_ordered_array_window_scores


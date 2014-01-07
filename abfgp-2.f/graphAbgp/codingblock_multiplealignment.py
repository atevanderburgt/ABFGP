"""
Multialignment functionality for CodingBlockGraphs
"""

# Module metadata
__authors__ = "Ate van der Burgt"
__license__ = "MIT"

# graphAbgp Imports
from graphAbgp.exceptions import OrganismNotPresentInGraph

# Abgp imports

# Python Imports

# Global Variable Import

def cbgproteinalignmentdata(cbg,organism,extra_aa_offset=10):
    """
    """
    if organism not in cbg.organism_set(): raise OrganismNotPresentInGraph
    omsr = cbg.overall_minimal_spanning_range(organism=organism)

    # aacoord2scores dict has keys (absolute) AA positions
    # and a list as value. This list has as first element
    # a tuple (aa,score); the encountered AA and its MATRIX score
    # All subsequence elements in the list are the MATRIX scores
    # for the aligned AA's in the PacbPORFs
    aacoord2scores = {}
    pacbporf_count = 0

    # get (a) ProteinSimilarityMatrix from the CBG 
    MATRIX = _get_protsim_matrix_from_cbg(cbg)

    for key, pacbporf in cbg.pacbps.iteritems():
        (a,b,c,d),(orgQ,orfQid),(orgS,orfSid) = key
        # continue if pacbporf is not from this identifier
        if organism not in [orgQ,orgS]: continue
        # Boolean variable stating pacbporf orientation (Q/S)
        pacbporf_count+=1
        is_organism_pacbp_query = orgQ == organism
        sta = pacbporf._original_alignment_pos_start - extra_aa_offset
        end = pacbporf._original_alignment_pos_end + extra_aa_offset
        sta = max([0, sta ])
        end = min([len(pacbporf._positions), end ])
        for pos in range(sta,end):
            algPos = pacbporf._positions[pos]
            if not algPos.isa_gap and is_organism_pacbp_query:
                aapos, aa, otheraa = algPos.query_pos, algPos.query, algPos.sbjct
            elif not algPos.isa_gap and not is_organism_pacbp_query:
                aapos, aa, otheraa = algPos.sbjct_pos, algPos.sbjct, algPos.query
            elif orgQ == organism and\
            (algPos.sbjct_dna_start,algPos.sbjct_dna_end) == (0,0):
                aapos, aa, otheraa = algPos.query_pos, algPos.query, "*"
            elif orgS == organism and\
            (algPos.query_dna_start,algPos.query_dna_end) == (0,0):
                aapos, aa, otheraa = algPos.sbjct_pos, algPos.sbjct, "*"
            else:
                # gap in the organism identifiers AA sequence
                continue

            # check for exceptiobnal cases of <0 aapos;
            # can happen for aligned sequences that fall at the
            # extreme beginning of the input DNA sequences
            if aapos < 0: continue

            # init aapos when it does not exist yet
            if not aacoord2scores.has_key(aapos):
                score = MATRIX.scoreaapair(aa,aa)
                aacoord2scores[aapos] = [ ( aa, score ) ]

            # update alignment data
            score = MATRIX.scoreaapair(aa,otheraa)
            aacoord2scores[aapos].append( score )

    # get lowest possible value from the matrix
    lowestscore = min(MATRIX.matrix.values())

    # correct for alignments that extend outside this orf -> non-coding
    orf = cbg.get_orfs_of_graph(organism=organism)[0]
    for aapos in range(min(aacoord2scores.keys()),orf.protein_startPY): 
       while len(aacoord2scores[aapos])-1 < pacbporf_count:
           aacoord2scores[aapos].append(lowestscore)
    for aapos in range(orf.protein_endPY,max(aacoord2scores.keys())+1):
       while len(aacoord2scores[aapos])-1 < pacbporf_count:
           aacoord2scores[aapos].append(lowestscore)

    # fill all `empty` positions. Non-complete lists in aacoord2scores[aapos]
    # means a stop codon is observed somewhere and therefor the
    # extended PacbPORF is ended
    for aapos in aacoord2scores.keys():
        while len(aacoord2scores[aapos])-1 < pacbporf_count:
           aacoord2scores[aapos].append(lowestscore)

    aapositions = aacoord2scores.keys()
    aapositions.sort()
    for aapos in aapositions:
        aa = aacoord2scores[aapos][0][0]
        positions = len(aacoord2scores[aapos])-1
        maxscore = aacoord2scores[aapos][0][1] * positions
        minscore = lowestscore * positions
        score = sum( aacoord2scores[aapos][1:] )
        try:
            ratio = float(score-minscore) / (maxscore-minscore)
        except:
            # minscore == maxscore: a gap as reference position !?
            # fixed in lib_hmm.py 03/04/2010... this error should
            # now be repaired.
            print cbg
            print organism, len(aapositions), aapos, aacoord2scores[aapos], score, minscore, maxscore
            ratio = float(score-minscore) / (maxscore-minscore)

        # store this ratio to this aapos
        aacoord2scores[aapos] = ratio

    # create leading and tailing 0.0 values
    for i in range(0,2):
       aacoord2scores[min(aacoord2scores)-1] = 0.0
       aacoord2scores[max(aacoord2scores)+1] = 0.0

    # return the aacoord2scores dict
    return aacoord2scores 

# end of function cbgproteinalignmentdata


def _get_protsim_matrix_from_cbg(cbg):
    """ """
    MATRIX = None
    for pacbporf in cbg.pacbps.values():
        if pacbporf.MATRIX:
            MATRIX = pacbporf.MATRIX
            break
    if not MATRIX:
        from settings.blastp import BLASTP_MATRIX_NAME
        from lib_proteinsimilaritymatrix import ProteinSimilarityMatrix
        MATRIX = ProteinSimilarityMatrix(name=BLASTP_MATRIX_NAME)
    # return the ProtSimMatrix
    return MATRIX

# end of function _get_protsim_matrix_from_cbg

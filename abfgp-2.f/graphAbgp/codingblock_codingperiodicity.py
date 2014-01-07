"""
Coding periodicity functions for CodingBlockGraphs
"""

# graphAbgp Imports
from exceptions import NoOverallMinimalSpanningRange

# ABFGP Imports
from lib_sequenceperiodicity import codon_nucleotide_conservation

# Python Imports


def codingperiodicitystatistics(cbg):
    """
    Obtain coding periodicity statistics of this CodingBlockGraph

    @attention: CodingBlockGraph *MUST* have an Overall Minimal Spanning Range

    @rtype  cncSummed:  list
    @return cncSummed: list-of-lists-of-lists; see cncSummed variable
    """
    # check if this CBG hasa OMSR
    if not cbg.has_overall_minimal_spanning_range():
        raise NoOverallMinimalSpanningRange, cbg

    cncSummed = [
        [ [ 0.0, 0.0, 0.0, 0 ] , [ 0.0, 0.0, 0.0, 0 ] , [ 0.0, 0.0, 0.0, 0 ] ],
        [ [ 0.0, 0.0, 0.0, 0 ] , [ 0.0, 0.0, 0.0, 0 ] , [ 0.0, 0.0, 0.0, 0 ] ],
        [ [ 0.0, 0.0, 0.0, 0 ] , [ 0.0, 0.0, 0.0, 0 ] , [ 0.0, 0.0, 0.0, 0 ] ],
        ]
    for k,pacbporf in cbg.pacbps.iteritems():
        # get dnaQuery and dnaSbjct sequences of the OMSR range of this Pacbporf
        aaStart   = min( cbg.overall_minimal_spanning_range( node = k[1]) )
        aaEnd     = max( cbg.overall_minimal_spanning_range( node = k[1]) )
        posStart  = pacbporf.alignmentposition_by_query_pos(aaStart)
        posEnd    = pacbporf.alignmentposition_by_query_pos(aaEnd) + 1
        dnaQ,dnaS = pacbporf.get_aligned_dna_sequences_by_algpositions(
                    posStart,posEnd)

        # obtain the pairwise codon_nucleotide_conservation data
        cnc = codon_nucleotide_conservation(dnaQ,dnaS)

        # add this data to the overall list-of-lists-of-lists `cncSummed`
        for phase in [0,1,2]:
            cncphaserow = cnc[phase]
            for i in range(0,len(cncphaserow)):
                for j in range(0,len(cncphaserow[i])):
                    cncSummed[phase][i][j]+=cncphaserow[i][j]

    # average all ratios and the number of positions taken into account
    for phase in [0,1,2]:
        for item in cncSummed[phase]:
            for i in range(0,len(item)): item[i] = item[i] / cbg.edge_count()

    # return the cncSummed list-of-lists-of-lists
    return cncSummed

# end of function codingperiodicitystatistics


def printcodingperiodicitystatistics(cbg):
    """
    Print nicely layouted coding periodicity statistics of this CodingBlockGraph

    @attention: CodingBlockGraph *MUST* have an Overall Minimal Spanning Range
    """
    periodicitydata = codingperiodicitystatistics(cbg)
    for phase in [0,1,2]:
        for item in periodicitydata[phase]:
            print "[%1.2f %1.2f %1.2f (%s)]" % tuple(item),
        print ""

# end of function printcodingperiodicitystatistics

"""
Ordering functions for CodingBlockGraphs
"""

# Module metadata
__authors__ = "Ate van der Burgt"
__license__ = "MIT"

# imports
from pacb.ordering import relpos2binaryrelpos, rel_pos_sets


def relatively_positioned_towards(cbg1,cbg2):
    """
    Return positional ordering data between 2 CodingBlockGraphs
    
    @type  cbg1: CodingBlockGraph
    @param cbg1: CodingBlockGraph instance

    @type  cbg2: CodingBlockGraph
    @param cbg2: CodingBlockGraph instance

    @rtype:  tuple
    @return: tuple with 7 data structures
    """
    # get union of organisms between two CBGs
    allorgs = cbg1.organism_set().union( cbg2.organism_set() )
    posRel  = []
    posBin  = []
    orfIdent= []
    # get overall minimal spanning range data structures for comparison
    omsrCbg1 = cbg1.overall_minimal_spanning_range()
    omsrCbg2 = cbg2.overall_minimal_spanning_range()
    for org in allorgs:
        if org in cbg1.organism_set() and org in cbg2.organism_set():
            nodeCbg1 = cbg1.node_by_organism(org)
            nodeCbg2 = cbg2.node_by_organism(org)
            if len(omsrCbg1[nodeCbg1]) == 0 or len(omsrCbg2[nodeCbg2]) == 0:
                continue
            # get relative positioning of the OMSRs
            (a1,a2,a3), (b1,b2,b3) = rel_pos_sets(
                    omsrCbg1[nodeCbg1],
                    omsrCbg2[nodeCbg2],
                    )
            # get binary relative positioning of the OMSRs
            (Ba1,Ba2,Ba3), (Bb1,Bb2,Bb3) = relpos2binaryrelpos(
                    (a1,a2,a3), (b1,b2,b3) )

            # append to lists with relative and binary position coords
            posRel.append( ( (a1,a2,a3), (b1,b2,b3) ) )
            posBin.append( ( (Ba1,Ba2,Ba3), (Bb1,Bb2,Bb3) ) )

            # check if the orfs are identical
            orfCbg1 = cbg1.get_orfs_of_graph(organism=org)[0]
            orfCbg2 = cbg2.get_orfs_of_graph(organism=org)[0]
            orfIdent.append( orfCbg1.id == orfCbg2.id )

    # make summed binary tuples
    summedPosBin = ( [0,0,0], [0,0,0] )
    for i in range(0,2):
        for j in range(0,3):
            summedPosBin[i][j] = sum([ item[i][j] for item in posBin ] )
    # and make binary tuples of this one
    ( binPosCbg1, binPosCbg2 ) = relpos2binaryrelpos(
            summedPosBin[0],
            summedPosBin[1],
            overlap_offset = 0,
            percentual_overlap_offset = 0.0
            )

    # return summed tuples, binary tuples and orfIdent for cbg1 and cbg2
    return ( summedPosBin[0], summedPosBin[1],
             binPosCbg1, binPosCbg2, orfIdent, posRel, posBin

# end of function relatively_positioned_towards

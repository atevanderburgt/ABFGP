"""
Ordering algorithms for python-pacb.
"""

# Module metadata
__authors__ = "Ate van der Burgt"
__license__ = "MIT"

# Imports
import pacb

from sets import Set

def relatively_positioned_towards(pacbp1,pacbp2):
        """
        What is the relative positioning of two pacbp's?
    
        Return data structure is a dictionary:
        {   'Q1': ( p,q,r ),
            'Q2': ( x,y,z ),
            'S1': ( p,q,r ),
            'S2': ( x,y,z ),
        }
    
        p   5p extention of `pacbp1` compared to `pacbp2`
        q   overlapping positions between `pacbp1` and `pacbp2`
        r   3p extention of  `pacbp1` compared to `pacbp2`
    
        x   5p extention of `pacbp2` compared to `pacbp1`
        y   overlapping positions between `pacbp2` and `pacbp1`
        z   3p extention of `pacbp2` compared to `pacbp1`
    
        Some examples, between brackets the length of the alignment:
        Coordinates are as from Query, same logic accounts for Sbjct
    
        (1) ================ (16)
        (2)                     ========== (10)
            p,q,r = 16,0,0
            x,y,z = 0,0,10
    
        (1) ================ (16)
        (2)               ========== (10)
            p,q,r = 14,2,0
            x,y,z = 0,2,8
    
        (1) ================ (16)
        (2)       ========== (10)
            p,q,r = 6,10,0
            x,y,z = 0,10,0
    
        (1) ================ (16)
        (2)  ========== (10)
            p,q,r = 1,10,5
            x,y,z = 0,10,0
    
        # pacbp's in reverse order
        (1)                  ================ (16)
        (2) ========== (10)
            p,q,r = 0,0,16
            x,y,z = 10,0,0
        """

        data = {}
        coords = {
            'query': {
                1: {
                    'startPY':  pacbp1.query_start,
                    'endPY':    pacbp1.query_end,
                },
                2: {
                    'startPY':  pacbp2.query_start,
                    'endPY':    pacbp2.query_end,
                }
            },
            'sbjct': {
                1: {
                    'startPY':  pacbp1.sbjct_start,
                    'endPY':    pacbp1.sbjct_end,
                },
                2: {
                    'startPY':  pacbp2.sbjct_start,
                    'endPY':    pacbp2.sbjct_end,
                },
            },
        }
        # PacbPORF object which can have extended alignemnt. 
        # In that case, *_start & *_end coords are not relevant.
        if pacbp1.__class__.__name__ == 'PacbPORF':
            start  = pacbp1._original_alignment_pos_start
            end    = pacbp1._original_alignment_pos_end
            staPos = pacbp1._positions[start]
            endPos = pacbp1._positions[end-1]
            # reset start & end coordinates of pacbp1
            coords['query'][1]['startPY'] = staPos.query_pos
            coords['query'][1]['endPY']   = endPos.query_pos
            coords['sbjct'][1]['startPY'] = staPos.sbjct_pos
            coords['sbjct'][1]['endPY']   = endPos.sbjct_pos
        if pacbp2.__class__.__name__ == 'PacbPORF':
            start  = pacbp2._original_alignment_pos_start
            end    = pacbp2._original_alignment_pos_end
            staPos = pacbp2._positions[start]
            endPos = pacbp2._positions[end-1]
            # reset start & end coordinates of pacbp1
            coords['query'][2]['startPY'] = staPos.query_pos
            coords['query'][2]['endPY']   = endPos.query_pos
            coords['sbjct'][2]['startPY'] = staPos.sbjct_pos
            coords['sbjct'][2]['endPY']   = endPos.sbjct_pos


        for strand in ['query','sbjct']:
            data[strand] = {}
            # seqid is from here on alias for pacbp1 (1) and pacbp2 (2)
            for seqid in [1,2]:
                # a coordinate Set enables easy overlap calculation
                coords[strand][seqid]['set'] =\
                    Set( range( coords[strand][seqid]['startPY'],
                    coords[strand][seqid]['endPY']+1 ) )
                # length of the alignments
                coords[strand][seqid]['length'] =\
                     len( coords[strand][seqid]['set'] )
                # coordinates list for 5p-extention, overlap 3p-extention
                data[strand][seqid] = [ 0,0,0 ]

        # store the overlap coordinate length
        data['query'][1][1] =\
            len(coords['query'][1]['set'].intersection(coords['query'][2]['set']))
        data['sbjct'][1][1] =\
            len(coords['sbjct'][1]['set'].intersection(coords['sbjct'][2]['set']))
        data['query'][2][1] = data['query'][1][1]
        data['sbjct'][2][1] = data['sbjct'][1][1]

        # now do the positioning towards each other
        for strand in ['query','sbjct']:
            if coords[strand][1]['endPY'] <= coords[strand][2]['startPY']:
                # pacbp1 FULLY in front of pacbp2
                data[strand][1][0] = coords[strand][1]['length']
                data[strand][2][2] = coords[strand][2]['length']
            elif coords[strand][2]['endPY'] <= coords[strand][1]['startPY']:
                # pacbp2 FULLY in front of pacbp1
                data[strand][1][2] = coords[strand][1]['length']
                data[strand][2][0] = coords[strand][2]['length']
            elif data[strand][1][1] == coords[strand][1]['length']:
                # pacbp1 included in pacbp2
                data[strand][2][0] = coords[strand][1]['startPY'] -\
                    coords[strand][2]['startPY']
                data[strand][2][2] = coords[strand][2]['endPY'] -\
                    coords[strand][1]['endPY']
            elif data[strand][2][1] == coords[strand][2]['length']:
                # pacbp2 included in pacbp1
                data[strand][1][0] = coords[strand][2]['startPY'] -\
                    coords[strand][1]['startPY']
                data[strand][1][2] = coords[strand][1]['endPY'] -\
                    coords[strand][2]['endPY']
            elif coords[strand][1]['startPY'] < coords[strand][2]['startPY']:
                # pacbp1 in front of pacbp2 but partially overlapping
                data[strand][1][0] = coords[strand][2]['startPY'] -\
                    coords[strand][1]['startPY']
                data[strand][2][2] = coords[strand][2]['endPY'] -\
                    coords[strand][1]['endPY']
            elif coords[strand][2]['startPY'] < coords[strand][1]['startPY']:
                # pacbp2 in front of pacbp1 but partially overlapping
                data[strand][2][0] = coords[strand][1]['startPY'] -\
                    coords[strand][2]['startPY']
                data[strand][1][2] = coords[strand][1]['endPY'] -\
                    coords[strand][2]['endPY']
            else:
                # what else !?!?! we have al options!
                raise "UNEXPECTED HSP ORDER!!"
    
        # and return!
        return {
            'Q1': tuple(data['query'][1]),
            'Q2': tuple(data['query'][2]),
            'S1': tuple(data['sbjct'][1]),
            'S2': tuple(data['sbjct'][2]),
        }
    
# end of function relatively_positioned_towards


def distance_towards(pacbp1,pacbp2):
    """
    Give the (AA) distance between two pacbps

    @attention: pacbp1 is expected to by positioned in front of pacbp
    @attention: incorrectly positioned or overlapping pacbps get distance 0

    @rtype:  positive integer or zero
    @return: AA distance between both pacbps
    """ 
    distQ = min(pacbp2.alignment_protein_range_query()) -\
            max(pacbp1.alignment_protein_range_query())
    distS = min(pacbp2.alignment_protein_range_sbjct()) -\
            max(pacbp1.alignment_protein_range_sbjct())
    if distQ <= 0:
        return 0
    elif distS <= 0:
        return 0
    else:
        return min([distQ,distS])

# end of function distance_towards


def overlap(pacbp1,pacbp2):
    """
    Give the overlap between two pacbps

    @rtype:  positive float
    @return: highest overlap ratio of query and sbjct of these 2 pacbps
    """
    relpositioning = relatively_positioned_towards(pacbp1,pacbp2)
    if relpositioning['Q1'][1] == 0 and relpositioning['S1'][1] == 0:
        return 0.0
    # if here, then there is some overlap!
    overlaps = []
    for item in relpositioning.keys():
        if relpositioning[item][1] > 0:
            overlaps.append( float(relpositioning[item][1]) / float(sum(relpositioning[item])) )
    # return the maximum od the overlaps
    return max(overlaps)
    
# end of function overlap


def query_overlap(pacbp1,pacbp2):
    """
    Give the overlap between the query sequences of two pacbps

    @rtype:  positive float
    @return: overlap ratio of query sequences of these 2 pacbps
    """
    relpositioning = relatively_positioned_towards(pacbp1,pacbp2)
    if relpositioning['Q1'][1] == 0:
        return 0.0
    else:
        q1o = relpositioning['Q1'][1] / float(sum(relpositioning['Q1']))
        q2o = relpositioning['Q2'][1] / float(sum(relpositioning['Q2']))
        return max([q1o,q2o])
    
# end of function query_overlap


def sbjct_overlap(pacbp1,pacbp2):
    """
    Give the overlap between the sbjct sequences of two pacbps

    @rtype:  positive float
    @return: overlap ratio of sbjct sequences of these 2 pacbps
    """
    relpositioning = relatively_positioned_towards(pacbp1,pacbp2)
    if relpositioning['S1'][1] == 0:
        return 0.0
    else:
        s1o = relpositioning['S1'][1] / float(sum(relpositioning['S1']))
        s2o = relpositioning['S2'][1] / float(sum(relpositioning['S2']))
        return max([s1o,s2o])
    
# end of function sbjct_overlap


def rel_pos_sets(rA,rB):
    """
    Relative positioning of two Array-like data structures (Sets, Arrays, Lists)

    @type  rA: Set
    @param rA: Array-like Set (supposed to be left of rB)

    @type  rB: Set
    @param rB: Array-like Set (supposed to be right of rA)

    How are 2 continious ranges of coordinates overlapping with each other?
    First argument is supposed to be LEFT / IN FRONT OF / BEFORE second argument.
    rA  Set     continious Set of coordinates; Set(range(A1,A2))
    rB  Set     continious Set of coordinates; Set(range(B1,B2))
    (a1,a2,a3)  tuple with tree absolute integer values >= 0
    (b1,b2,b3)  tuple with tree absolute integer values >= 0
    x1          outsticking coordinates on the LEFT  side of coordinate range
    x2          overlapping coordinates in both coordinate ranges
    x3          outsticking coordinates on the RIGTH side of coordinate range
    Example 1:
        rA = Set(list(range(10,31))
        rB = Set(list(range(28,31))
    Outcome 1:
        (a1,a2,a3) = (18,3,0)
        (b1,b2,b3) = (0,3,0)
    Example 2:
        rA = Set(list(range(10,31))
        rB = Set(list(range(6,28))
    Outcome 2:
        (a1,a2,a3) = (0,18,3)
        (b1,b2,b3) = (4,18,0)
    """
    (a1,a2,a3) = (0,0,0)
    (b1,b2,b3) = (0,0,0)
    if max(rA) < min(rB):
        # order is ra - rB, no overlap
        a1,b3 = len(rA), len(rB)
    elif min(rA) > max(rB):
        # order is rB - ra, no overlap
        a3,b1 = len(rA), len(rB)
    elif rA.symmetric_difference(rB) == Set():
        # rA and rB overlap perfectly
        b2 = len(rA.intersection(rB))
        a2 = b2
    else:
        # rA and rB (partially) overlap
        b2 = len(rA.intersection(rB))
        a2 = b2
        difA = rA.difference(rB)
        difB = rB.difference(rA)
        if difA:
            # rA has outsticking parts in respect to rB
            difAleft   = difA.difference(range(max(rB),max(rA)+1))
            difArigth  = difA.difference(range(min(rA),min(rB)+1))
            a1, a3     = len(difAleft), len(difArigth)
        if difB:
            # rB has outsticking parts in recpect to rA
            difBleft   = difB.difference(range(max(rA),max(rB)+1))
            difBrigth  = difB.difference(range(min(rB),min(rA)+1))
            b1, b3     = len(difBleft), len(difBrigth)
    # return the relative positioning tuples
    return (a1,a2,a3), (b1,b2,b3)

# end of function rel_pos_sets


def relpos2binaryrelpos( (a1,a2,a3), (b1,b2,b3), overlap_offset = 3, percentual_overlap_offset = 1.5 ):
    """
    """
    aret, bret = [], []
    zerorange = range(1,overlap_offset+1)
    for a in (a1,a2,a3):
        if a == 0:                  aret.append(0)        
        elif a > overlap_offset:    aret.append(1)
        elif a in zerorange:        aret.append(0)
        else:                       aret.append(-1)
    for b in (b1,b2,b3):
        if b == 0:                  bret.append(0)
        elif b > overlap_offset:    bret.append(1)
        elif b in zerorange:        bret.append(0)
        else:                       bret.append(-1)

    # Check for creation of an all-zeros tuple (0,0,0)
    # This can happen for absolute sizes <=overlap_offset
    # Correct this here!
    if aret == [0,0,0]:
        for pos in range(0,3):
            if (a1,a2,a3)[pos] in zerorange: aret[pos] = 1
    if bret == [0,0,0]:
        for pos in range(0,3):
            if (b1,b2,b3)[pos] in zerorange: bret[pos] = 1

    # check for the percentual overlap offset
    if not a3 and a1 and 100.0*(float(a2)/float(a1)) <= percentual_overlap_offset: aret[1] = 0
    if not a1 and a3 and 100.0*(float(a2)/float(a3)) <= percentual_overlap_offset: aret[1] = 0
    if not b3 and b1 and 100.0*(float(b2)/float(b1)) <= percentual_overlap_offset: bret[1] = 0
    if not b1 and b3 and 100.0*(float(b2)/float(b3)) <= percentual_overlap_offset: bret[1] = 0

    # and return the binary positioning tuple
    return ( tuple(aret), tuple(bret) )

# end of function relpos2binaryrelpos


def is_identical_binaryrelpos( (a1,a2,a3), (b1,b2,b3) ):
    """
    """
    if (a1,a2,a3) == (b1,b2,b3):
        return True
    else:
        return False

# end of function is_identical_binaryrelpos


def is_included_binaryrelpos( (a1,a2,a3), (b1,b2,b3) ):
    """
    """
    if (0,1,0) in [ (a1,a2,a3), (b1,b2,b3) ]:
        return True
    else:
        return False

# end of function is_included_binaryrelpos


def reverse_binaryrelpos( (a1,a2,a3) ):
    """
    """
    return ( a3,a2,a1 )

# end of function reverse_binaryrelpos


def issubsetorsuperset(pacbp1,pacbp2):
    """
    Is pacbp2 included in pacbp1 or vise versa?
    """
    if issubset(pacbp1,pacbp2) or issubset(pacbp2,pacbp1):
        return True
    else:
        return False

# end of function issubsetorsuperset


def issubset(pacbp1,pacbp2):
    """
    Is pacbp1 included in pacbp2?
    """ 
    setQ1 = Set(pacbp1.alignment_protein_range_query())
    setQ2 = Set(pacbp2.alignment_protein_range_query())
    setS1 = Set(pacbp1.alignment_protein_range_sbjct())
    setS2 = Set(pacbp2.alignment_protein_range_sbjct())
        
    if setQ1.issubset(setQ2) and setS1.issubset(setS2):
        return True 
    else:
        return False
        
# end of function issubset


def issuperset(pacbp1,pacbp2):
    """
    Is pacbp1 fully covering pacbp2?
    """
    return issubset(pacbp2,pacbp1) 

# end of function issuperset


def order_pacbporf_list(pacbps):
    """ """
    ordered = [ ( p._get_original_alignment_pos_start().query_pos, p ) for p in pacbps ]
    ordered.sort()
    return [ p for (coord,p) in ordered ]
# end of function order_pacbporf_list


def order_pacbp_list(pacbps):
    """ """
    ordered = [ ( p.query_start, p ) for p in pacbps ]
    ordered.sort()
    return [ p for (coord,p) in ordered ]
# end of function order_pacbp_list


def order(pacbps):
    """
    """
    ordered_pacbps = []
    for pacbp in pacbps:
        ordered_pacbps.append( ( pacbp.query_start, pacbp ) )
    ordered_pacbps.sort()
    ordered_pacbps = [ pacbp for (order,pacbp) in ordered_pacbps ]
    print 0, ordered_pacbps[0]
    for i in range(1,len(ordered_pacbps)):
        print relatively_positioned_towards(ordered_pacbps[i-1],ordered_pacbps[i])
        print i, ordered_pacbps[i]

# end of function order


def order_list_by_attribute(inputlist,order_by='',reversed=False):
    """
    Helper function to order objects in a list by any given attribute

    @type  inputlist: list
    @param inputlist: list of items to be ordered on a attribute

    @type  order_by: string
    @param order_by: attribute of items in inputlist to be sorted upon

    @type  reversed: Boolean
    @param reversed: reverse the ordered list (True) or not (False)

    @rtype:  list
    @return: ordered list of items, ordered on order_by attribute
    """
    if not inputlist: return inputlist
    if hasattr(inputlist[0],str(order_by)):
        orderedlist = [ ( getattr(o,str(order_by)), o ) for o in inputlist ]
        orderedlist.sort()
        if reversed: orderedlist.reverse()
        return [ o for (attr,o) in orderedlist ]
    else:
        return inputlist

# end of function order_list_by_attribute



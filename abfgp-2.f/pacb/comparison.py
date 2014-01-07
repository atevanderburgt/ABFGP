"""
Comparison functions for Pacb class
"""

# Module metadata
__authors__ = "Ate van der Burgt"
__license__ = "MIT"


def IsIdenticalPacbPORF(objA,objB,quickanddirty=True):
    """
    Are 2 objects PacbPORFs and indentical?

    @type  objA: PacbPORF
    @param objA: PacbPORF object

    @type  objB: PacbPORF
    @param objB: PacbPORF object

    @type  quickanddirty: Boolean
    @param quickanddirty: if True, check only AAs, not DNAs

    @rtype:  Boolean
    @return: Boolean
    """
    try:
        if objA.__class__.__name__  == 'PacbPORF' and objB.__class__.__name__  == 'PacbPORF':
            if objA.query_start == objB.query_start and\
            objA.sbjct_start == objB.sbjct_start and\
            objA.query_end == objB.query_end and\
            objA.sbjct_end == objB.sbjct_end:
                # coordinates are the same; check if aa's are identical too
                if objA.query == objB.query and objA.sbjct == objB.sbjct:
                    # TODO: param quickanddirty not implemented yet at this stage;
                    # it is not really neccesairy because the above checks seem bulletproof
                    return True
                else:
                    return False
            else:
                return False
        else:
            return False
    except:
        # no __class__.__name__ attributes -> cannot be PacbPORFs
        return False

# end of function IsIdenticalPacbPORF


def IsPacbPORFSubstringOf(objA,objB,quickanddirty=True):
    """
    Is a PacbPORF object a `substring` of the other PacbPORF?

    @type  objA: PacbPORF
    @param objA: PacbPORF object

    @type  objB: PacbPORF
    @param objB: PacbPORF object

    @type  quickanddirty: Boolean
    @param quickanddirty: if True, check only AAs, not DNAs

    @attention: function is more or less identical to pacb.ordering.issubset() but much more elaborate

    @rtype:  Boolean
    @return: Boolean
    """
    try:
        if objA.__class__.__name__  == 'PacbPORF' and objB.__class__.__name__  == 'PacbPORF':
            offsetAinBquery = objB.query.find(objA.query)
            offsetAinBsbjct = objB.sbjct.find(objA.sbjct)
            if offsetAinBquery >= 0 and offsetAinBsbjct >=0 and offsetAinBquery == offsetAinB: 
                # alignment of objA is included in objB; check the coordinates
                if objA._positions[0].query_pos == objA._positions[offsetAinBquery].query_pos and\
                objA._positions[0].sbjct_pos == objA._positions[offsetAinBquery].sbjct_pos:
                    # TODO: param quickanddirty not implemented yet at this stage;
                    # it is not really neccesairy because the above checks seem bulletproof
                    return True
                else:
                    return False
            else:
                return False
        else:
            return False
    except:
        # no __class__.__name__ attributes -> cannot be PacbPORFs
        return False

# end of function IsPacbPORFSubstringOf 



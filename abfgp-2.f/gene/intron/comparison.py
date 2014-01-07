"""
Functions for comparing 2 IntronConnectingOrf objects
"""

# Module metadata
__authors__ = "Ate van der Burgt"
__license__ = "MIT"


# Global variable Imports
from settings.splicesites import OPTIMAL_BRACNHPOINT_TO_ACCEPTOR_DISTANCE
MAXIMAL_OPTIMAL_BRACNHPOINT_TO_ACCEPTOR_DISTANCE = 5


def _branchpoint_comparison(intron1,intron2):
    """
    Is the branchpoint of intron2 an improvement over the branchpoint of intron1?

    @type  intron1: IntronConnectingOrfs object
    @param intron1: IntronConnectingOrfs object

    @type  intron2: IntronConnectingOrfs object
    @param intron2: IntronConnectingOrfs object

    @attention: intron1 is regarded as H0, intron2 as a H1 hypothesis

    @rtype:  NoneBoolean
    @return: NoneBoolean (improved,indifferent,deteriorated)
    """
    # return variable
    return_noneboolean = None

    if not intron1.branchpoint and not intron2.branchpoint:
        return_noneboolean = None

    elif intron1.branchpoint and intron2.branchpoint:

        intron1_bp_dist = intron1.get_branchpoint_nt_distance()
        intron1_bp_optimality = min([ abs(offset-intron1_bp_dist) for offset in OPTIMAL_BRACNHPOINT_TO_ACCEPTOR_DISTANCE ])
        intron2_bp_dist = intron2.get_branchpoint_nt_distance()
        intron2_bp_optimality = min([ abs(offset-intron2_bp_dist) for offset in OPTIMAL_BRACNHPOINT_TO_ACCEPTOR_DISTANCE ])

        if intron1_bp_optimality == intron2_bp_optimality:
            return_noneboolean = None
        elif intron2_bp_optimality < intron1_bp_optimality:
            if intron2_bp_optimality <= MAXIMAL_OPTIMAL_BRACNHPOINT_TO_ACCEPTOR_DISTANCE:
                return_noneboolean = True
            else:
                # both H0 and H1 acceptor site are poorly spaced towards bp
                # do not allow an improvement of this intron position
                # based on branchpoint
                return_noneboolean = False
        elif intron1_bp_optimality < intron2_bp_optimality:
            if intron1_bp_optimality <= MAXIMAL_OPTIMAL_BRACNHPOINT_TO_ACCEPTOR_DISTANCE:
                return_noneboolean = False
            else:
                # both H0 and H1 acceptor site are poorly spaced towards bp
                # do not allow an improvement of this intron position
                # based on branchpoint
                return_noneboolean = False
        else:
            # cannot happen
            pass

    elif intron1.branchpoint:

        # no branchpoint in intron2 !!
        intron1_bp_dist = intron1.get_branchpoint_nt_distance()
        intron1_bp_optimality = min([ abs(offset-intron1_bp_dist) for offset in OPTIMAL_BRACNHPOINT_TO_ACCEPTOR_DISTANCE ])
        if intron1_bp_optimality <= MAXIMAL_OPTIMAL_BRACNHPOINT_TO_ACCEPTOR_DISTANCE:
            # perfect brachpoint in intron1
            return_noneboolean = False
        else:
            # although horrible, 2th intron does not deliver any obvious branchpoint
            return_noneboolean = None

    elif intron2.branchpoint:

        # no branchpoint in intron1 !!
        intron2_bp_dist = intron2.get_branchpoint_nt_distance()
        intron2_bp_optimality = min([ abs(offset-intron2_bp_dist) for offset in OPTIMAL_BRACNHPOINT_TO_ACCEPTOR_DISTANCE ])
        if intron2_bp_optimality <= MAXIMAL_OPTIMAL_BRACNHPOINT_TO_ACCEPTOR_DISTANCE:
            return_noneboolean = True
        else:
            # although horrible, 1th intron does not deliver any obvious branchpoint
            return_noneboolean = None
    else:
        # cannot happen!
        pass

    # return the noneboolean outcome
    return return_noneboolean

# end of function _branchpoint_comparison


def _polypyrimidinetract_comparison(intron1,intron2):
    """
    Are the polypyrimidinetract(s) of intron2 an improvement over those of intron1?

    @type  intron1: IntronConnectingOrfs object
    @param intron1: IntronConnectingOrfs object

    @type  intron2: IntronConnectingOrfs object
    @param intron2: IntronConnectingOrfs object

    @rtype:  NoneBoolean
    @return: NoneBoolean (improved,indifferent,deteriorated)
    """

    # return variable
    return_noneboolean = None

    intron1ppt = 0
    intron2ppt = 0
    if intron1.ppt5p: intron1ppt+=intron1.ppt5p.length
    if intron1.ppt3p: intron1ppt+=intron1.ppt3p.length
    if intron2.ppt5p: intron2ppt+=intron2.ppt5p.length
    if intron2.ppt3p: intron2ppt+=intron2.ppt3p.length


    if not intron1ppt and not intron2ppt:
        return_noneboolean = None
    elif intron1ppt and intron2ppt:
        if intron1ppt > intron2ppt and abs(intron1ppt - intron2ppt) >= 2:
            return_noneboolean = False
        elif intron2ppt > intron1ppt and abs(intron1ppt - intron2ppt) >= 2:
            return_noneboolean = True
        else:
            return_noneboolean = None
    elif intron1ppt:
        return_noneboolean = False
    else:
        return_noneboolean = True

    # return the noneboolean outcome
    return return_noneboolean

# end of function _polypyrimidinetract_comparison


def _algsimilarity_comparison(intron1,intron2,prev_exon,next_exon,algsimilarity):
    """
    Is intron2 an improvement over intron1 in terms of AlgimentSimilarity coverage?

    @type  intron1: IntronConnectingOrfs object
    @param intron1: IntronConnectingOrfs object

    @type  intron2: IntronConnectingOrfs object
    @param intron2: IntronConnectingOrfs object

    @type  prev_exon: ExonOnOrf object or None
    @param prev_exon: ExonOnOrf object to test different donor site on

    @type  next_exon: ExonOnOrf object or None
    @param next_exon: ExonOnOrf object to test different acceptor site on

    @type  algsimilarity: Numerical array
    @param algsimilarity: AlignmentSimilarity (or AlignmentPresence) array

    @rtype:  NoneBoolean
    @return: NoneBoolean (improved,indifferent,deteriorated)
    """

    # return variable
    return_noneboolean = None

    algsim_intron1 = sum(algsimilarity[intron1.donor.pos/3:intron1.acceptor.pos/3])
    algsim_intron2 = sum(algsimilarity[intron2.donor.pos/3:intron2.acceptor.pos/3])
    if next_exon:
        algsim_next_exon1 = sum(algsimilarity[intron1.acceptor.pos/3:next_exon.donor.pos/3])
        algsim_next_exon2 = sum(algsimilarity[intron2.acceptor.pos/3:next_exon.donor.pos/3])
        if algsim_next_exon1 == 0 and (intron1.acceptor.pos/3) == (next_exon.donor.pos/3):
            # avoid ZeroDivisionError
            algsim_next_exon1 = sum(algsimilarity[(intron1.acceptor.pos/3)-1:(next_exon.donor.pos/3)+1])
        if algsim_next_exon2 == 0 and (intron2.acceptor.pos/3) == (next_exon.donor.pos/3):
            # avoid ZeroDivisionError
            algsim_next_exon2 = sum(algsimilarity[(intron2.acceptor.pos/3)-1:(next_exon.donor.pos/3)+1])
    if prev_exon:
        algsim_prev_exon1 = sum(algsimilarity[prev_exon.acceptor.pos/3:intron1.donor.pos/3])
        algsim_prev_exon2 = sum(algsimilarity[prev_exon.acceptor.pos/3:intron2.donor.pos/3])
        if algsim_prev_exon1 == 0 and (prev_exon.acceptor.pos/3) == (intron1.donor.pos/3):
            # avoid ZeroDivisionError
            algsim_prev_exon1 = sum(algsimilarity[(prev_exon.acceptor.pos/3)-1:(intron1.donor.pos/3)+1])
        if algsim_prev_exon2 == 0 and (prev_exon.acceptor.pos/3) == (intron2.donor.pos/3):
            # avoid ZeroDivisionError
            algsim_prev_exon2 = sum(algsimilarity[(prev_exon.acceptor.pos/3)-1:(intron2.donor.pos/3)+1])

    if next_exon and prev_exon:
        # avoid ZeroDivisionError
        if algsim_next_exon1 and not algsim_next_exon2: return False
        if algsim_next_exon2 and not algsim_next_exon1: return True
        if algsim_prev_exon1 and not algsim_prev_exon2: return False
        if algsim_prev_exon2 and not algsim_prev_exon1: return True
        print (algsim_prev_exon1,algsim_prev_exon2), (algsim_next_exon1,algsim_next_exon2)
        print intron1
        print intron2

        if (algsim_prev_exon1,algsim_prev_exon2) == (0,0):
            # omit prev_exon -> both introns yield 0 similarity score for this exon
            return _algsimilarity_comparison(intron1,intron2,None,next_exon,algsimilarity)
        if (algsim_next_exon1,algsim_next_exon2) == (0,0):
            # omit next_exon -> both introns yield 0 similarity score for this exon
            return _algsimilarity_comparison(intron1,intron2,prev_exon,None,algsimilarity)


        # calculate alignmentsimilarity ratios
        algsim_ratio1 = float(algsim_intron1)/algsim_prev_exon1 + float(algsim_intron1)/algsim_next_exon1
        algsim_ratio1 = algsim_ratio1/2
        algsim_ratio2 = float(algsim_intron2)/algsim_prev_exon2 + float(algsim_intron2)/algsim_next_exon2
        algsim_ratio2 = algsim_ratio2/2
    elif next_exon:
        # avoid ZeroDivisionError
        if algsim_next_exon1 and not algsim_next_exon2: return False
        if algsim_next_exon2 and not algsim_next_exon1: return True
        if not algsim_next_exon2 and not algsim_next_exon1: return None
        # calculate alignmentsimilarity ratios
        algsim_ratio1 = float(algsim_intron1)/algsim_next_exon1
        algsim_ratio2 = float(algsim_intron2)/algsim_next_exon2
    elif prev_exon:
        # avoid ZeroDivisionError
        if algsim_prev_exon1 and not algsim_prev_exon2: return False
        if algsim_prev_exon2 and not algsim_prev_exon1: return True
        if not algsim_prev_exon2 and not algsim_prev_exon1: return None
        algsim_ratio1 = float(algsim_intron1)/algsim_prev_exon1
        algsim_ratio2 = float(algsim_intron2)/algsim_prev_exon2
    else:
        # no exons applied -> cannot be calculated
        algsim_ratio1 = 1
        algsim_ratio2 = 1

    if algsim_ratio2 < algsim_ratio1:
        return_noneboolean = True
    elif algsim_ratio1 > algsim_ratio2:
        return_noneboolean = False
    else:
        return_noneboolean = None

    # return the noneboolean outcome
    return return_noneboolean

# end of function _algsimilarity_comparison


def _finetune_splicesite_comparison(site1,site2):
    """
    Is the splice site2 (of intron2) an improvement over the site of intron1?

    @type  site1: SpliceSite object
    @param site1: SpliceSite object

    @type  site2: SpliceSite object
    @param site2: SpliceSite object

    @rtype:  NoneBoolean
    @return: NoneBoolean (improved,indifferent,deteriorated)
    """

    # return variable
    return_noneboolean = None

    # indifferently_ratio; 12.5% affinity range to be indifferent
    indifferently_ratio = 0.125

    if site2.pssm_score > site1.pssm_score:
        difference = site2.pssm_score - site1.pssm_score
        try:
            if difference/site1.pssm_score < indifferently_ratio:
                return_noneboolean = None
            else:
                return_noneboolean = True
        except ZeroDivisionError:
            return_noneboolean = True
    else:
        difference = site1.pssm_score - site2.pssm_score
        try:
            if difference/site2.pssm_score < indifferently_ratio:
                return_noneboolean = None
            else:
                return_noneboolean = False
        except ZeroDivisionError:
            return_noneboolean = False

    # return the noneboolean outcome
    return return_noneboolean

# end of function _finetune_splicesite_comparison

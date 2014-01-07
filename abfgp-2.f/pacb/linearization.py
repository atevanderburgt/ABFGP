"""
Linearization algorithms for python-pacb.

``overlaps`` is a dictionary of:
key     tuple( pacbp.bits, pacbp.length, orfpointerQid, orfpointerSid )
value   PacbP object belonging to this key

"""

# Module metadata
__authors__ = "Ate van der Burgt"
__license__ = "MIT"

# Imports
import pacb
import recombination
import ordering

def linearise_pacbs(overlaps,ACCEPTANCE_MARGIN=333,start_with_bests=1,is_weigthed=True,order_by='bits'):
    """
    Filter all overlaps based on linearization
    
    start_with_bests means how many of the top-best
    overlaps/pacbs must be used for the calibration
    of the linearization.

    WEIGHTED is done based on the the bitscore of the pacbps.
    Pacbps with high bitscores contribute more to the average
    position then those with low bitscores.

    distadvantage of start_with_best == 1
    TODO: take central best scoring pacb as starting point 

    @type  order_by: string
    @param order_by: 'length' (DESC), 'bits' (DESC), default 'bits'

    """
    # no overlaps means return the (0) overlaps and (0) removed onces
    if not overlaps: return ( {}, {} )

    # only 1 key? nothing to linearize
    if len(overlaps) == 1: return ( overlaps, {} )

    # convert dict of overlaps to list of pacbplist -> take the values
    pacbplist = overlaps.values()

    # linearize them in the new fuction
    linear, removed = linearise_pacbp_list(
            pacbplist,ACCEPTANCE_MARGIN=ACCEPTANCE_MARGIN,
            start_with_bests=start_with_bests,is_weigthed=is_weigthed,
            order_by=order_by
            )
    
    # and re-map to the original overlap dictionary
    if not removed:
        # all accepted -> just return AS IS
        return ( overlaps, {} )
    else:
        removed_pacbp_strings = [ str(pacbp) for pacbp in removed ]
        removed_overlaps = {}
        # set removed ones to removed_overlaps
        for key, pacbp in overlaps.iteritems():
            if str(pacbp) in removed_pacbp_strings:
                removed_overlaps[key] = pacbp
        # and delete all keys from overlaps that are now in removes
        for key in removed_overlaps.keys():
            del( overlaps[key] )
        # and return the linearized pacbs and the removed ones
        return ( overlaps, removed_overlaps )

# end of function linearise_pacbs


def linearise_pacbp_list(pacbplist,ACCEPTANCE_MARGIN=333,start_with_bests=1,is_weigthed=True,order_by='bits'):
    """
    """
    # check how much pacps are applied in the list
    if len(pacbplist) == 0:
        return ( [], [] )
    elif len(pacbplist) == 1:
        return ( pacbplist, [] )
    else:
        pass

    # check if order_by is properly applied
    if order_by not in ['bits','length']: order_by ='bits'

    # order the pacbplist
    if order_by == 'bits':
        # order pacbplist, highest bitscore first
        pacbplist = pacb.ordering.order_list_by_attribute(pacbplist,order_by='bits',reversed=True)
        # get lowest absolute and non-zero bitscore for weighting
        LOWEST_WEIGTH = max( [ 1, pacbplist[-1].bits ] )
    else:
        # order pacbplist, longest pacbp first
        pacbplist = pacb.ordering.order_list_by_attribute(pacbplist,order_by='length',reversed=True)
        # get lowest absolute and non-zero bitscore for weighting
        LOWEST_WEIGTH = min( [ pacbp.bits for pacbp in pacbplist ] )
        LOWEST_WEIGTH = max( [ 1, LOWEST_WEIGTH ] )

    # check if start_with_bests is properly given
    if start_with_bests not in range(1,6):
        # take default (1) start_with_bests
        start_with_bests = 1


    # return lists with colinear and removed pacbps
    linear_pacbps  = []
    removed_pacbps = []
    # list to store the Linearization Calibration Overlaps into
    lco = []

    # Take the best pacbp(s) to calibrate the linearization proces
    # However, check if we are not about to create a nonsense
    # linearization upon adding the n-th element
    for pacbp in pacbplist[0:start_with_bests]:
        # get offset/weigth data and set update lco
        if _is_compatible_positioned(pacbp,linear_pacbps):
            ( offset, weigth ) = _getpacbpdata( pacbp, LOWEST_WEIGTH, is_weigthed=is_weigthed )
            lco.append( ( offset, weigth ) )
            linear_pacbps.append( pacbp )
        else:
            removed_pacbps.append( pacbp )

    # Now loop over the remaining pacbps
    for pacbp in pacbplist[start_with_bests:]:
        # get Linear Calibration Offset of current overlap
        ( offset, weigth ) = _getpacbpdata( pacbp, LOWEST_WEIGTH, is_weigthed=is_weigthed )
        # get range of accepted Linear Calibration Offsets
        if offset in _get_lco_range(lco,ACCEPTANCE_MARGIN):
            # check the ordering towards the other pacbps
            if _is_compatible_positioned(pacbp,linear_pacbps):
                # compatible positioning; save and set to lco
                linear_pacbps.append( pacbp )
                lco.append( ( offset, weigth ) )
            else:
                # incompatible positioning
                removed_pacbps.append( pacbp )
        else:
            # distance is to far to be accepted
            removed_pacbps.append( pacbp )

    # actually but the pacbps in linear order
    #linear_pacbps = pacb.ordering.order_list_by_attribute( linear_pacbps, order_by='query_start' )

    # and return the linearized pacbs and the removed ones
    return ( linear_pacbps, removed_pacbps )

# end of function linearize_pacbp_list


def remove_repeats_from_pacbp_list(pacbplist,overlap_ratio=0.85):
    """
    """
    if not pacbplist: return pacbplist
    # order pacbplist by `bits` attribute, highest first
    ordered = pacb.ordering.order_list_by_attribute(pacbplist,order_by='bits',reversed=True)
    # make upper cross of element indexes in ordered list of pacbps
    pairs = recombination.pairwise( range(0,len(ordered)) )
    # loop over all pacbp combis, calculate overlap and store the
    # index to `toberemoved` when higher than `overlap_ratio`
    toberemoved = []
    for (posA,posB) in pairs:
        if posA in toberemoved: continue
        if posB in toberemoved: continue
        pacbpA = ordered[posA]
        pacbpB = ordered[posB]
        if ordering.overlap(pacbpA,pacbpB) >= overlap_ratio:
            toberemoved.append(posB)
    # order `toberemoved` and pop the (pacbp) elements from the ordered input list
    toberemoved.sort()
    toberemoved.reverse()
    for pos in toberemoved: ordered.pop(pos)
    # return the remainder of the odered input list
    return ordered 

# end of function remove_repeats_from_pacbp_list


def _is_compatible_positioned(pacbp,acceptedpacbps):
    """
    @rtype:  Boolean
    @return: True or False
    """
    is_compatible = True
    for acceptedPacbp in acceptedpacbps:
        order = pacb.ordering.relatively_positioned_towards(acceptedPacbp,pacbp)
        q1,q2,s1,s2 = order['Q1'], order['Q2'], order['S1'], order['S2']
        q1b, q2b = pacb.ordering.relpos2binaryrelpos(q1,q2)
        s1b, s2b = pacb.ordering.relpos2binaryrelpos(s1,s2)
        q2br = pacb.ordering.reverse_binaryrelpos(q2b)
        s2br = pacb.ordering.reverse_binaryrelpos(s2b)
        # mind taking the correct variables here: q1b and s2b , s1b and q2b
        if pacb.ordering.is_identical_binaryrelpos(q1b,s2b) and pacb.ordering.is_identical_binaryrelpos(s1b,q2b):
            # query & sbjct position contradict each other compared to previously accepted pacbp
            is_compatible = False
            break
        # mind taking the correct variables here: q2br and s2br
        elif sum([ q2br[0] != s2br[0], q2br[1] != s2br[1], q2br[2] != s2br[2] ]) == 3:
            # Q and S are not congruently positioned towards what accepted pabcp
            is_compatible = False
            break
        # mind taking the correct variables here: q1b and q1b2 , s1b and s2br
        elif pacb.ordering.is_identical_binaryrelpos(q1b,q2br) and pacb.ordering.is_identical_binaryrelpos(s1b,s2br):
            # normal order, Q << S or Q >> S
            pass
        # mind taking the correct variables here: q1b and q1b2 , s1b and s2br
        elif pacb.ordering.is_included_binaryrelpos(q1b,q2br) and pacb.ordering.is_included_binaryrelpos(s1b,s2br):
            # Q includes S or S includes Q
            pass
        else:
            # a more messy relation, not quickly resolved
            # likely, most but NOT ALL of these, are false
            pass

    # return the result status
    return is_compatible

# end of function _is_compatible_positioned


def _getpacbpdata(pacbp,LOWEST_WEIGTH,is_weigthed=True):
    """
    """
    q1,q2  = pacbp.query_start, pacbp.query_end 
    s1,s2  = pacbp.sbjct_start, pacbp.sbjct_end 
    ll     = (q2-q1+s2-s1)/2
    offset = ( q1-s1 + q2-s2 ) / 2
    if is_weigthed:
        # Get relatve weigth.
        # In case of negative bits, set to a very low positive weight (~0.1)
        weigth = max( [ 0.1, float(pacbp.bits) / LOWEST_WEIGTH ] )
    else:
        weigth = 1.0
    # return offset and weigth
    return (offset, weigth)

# end of function _getpacbpdata


def _get_average_lco(lco):
    """
    """
    tot  = sum( [ t[0]*t[1] for t in lco ] )
    wgth = sum( [ t[1] for t in lco ] )
    return int(tot/wgth)

# end of function _get_average_lco


def _get_lco_range(lco,ACCEPTANCE_MARGIN):
    """
    """
    avlco = _get_average_lco(lco)
    return range(avlco-ACCEPTANCE_MARGIN,avlco+ACCEPTANCE_MARGIN+1)

# end of function _get_lco_range


def _get_lco_range_from_average(averagelco,ACCEPTANCE_MARGIN):
    """
    """
    return range(averagelco-ACCEPTANCE_MARGIN,averagelco+ACCEPTANCE_MARGIN+1)

# end of function _get_lco_range_from_average


def _get_lco_of_overlaps(overlaps,is_weigthed=True):
    """
    get current average lco
    """
    keys = overlaps.keys()
    keys.sort()
    keys.reverse()
    # get lowest absolute and non-zero bitscore for weighting
    LOWEST_WEIGTH = max( [ 1, float(keys[-1][0]) ] )
    lco = []

    # take the best overlap(s) to calibrate
    # the linearization proces
    # however, check if we are not about to create
    # a nonsense linearization upon adding the n-th element
    for key in keys:
        # get offset/weigth data and set update lco
        ( offset, weigth ) = _getoverlapdata( key, overlaps[key], LOWEST_WEIGTH, is_weigthed=is_weigthed )
        lco.append( ( offset, weigth ) )
    # and return
    return lco

# end of function _get_lco_of_overlaps


def DEPRECATED_linearise_pacbs(overlaps,ACCEPTANCE_MARGIN=333,start_with_bests=1,is_weigthed=True,order_by='bits'):
    """
    Filter all overlaps based on linearization

    start_with_bests means how many of the top-best
    overlaps/pacbs must be used for the calibration
    of the linearization.

    WEIGHTED is done based on the the bitscore of the pacbps.
    Pacbps with high bitscores contribute more to the average
    position then those with low bitscores.

    distadvantage of start_with_best == 1
    TODO: take central best scoring pacb as starting point

    @type  order_by: string
    @param order_by: 'length' (DESC), 'bits' (DESC), default 'bits'

    """
    # check if order_by is properly applied
    if order_by not in ['bits','length']: order_by ='bits'

    # the return dictionairy of overlaps
    # that are co-linear with each other
    linear_overlaps  = {}
    removed_overlaps = {}

    # no overlaps means return the (0) overlaps and (0) removed onces
    if not overlaps: return ( overlaps, removed_overlaps )

    if order_by == 'bits':
        # get order overlap keys, highest bitscore first
        keys = overlaps.keys()
        keys.sort()
        keys.reverse()
    else:
        # get order overlap keys, longest pacbp first
        keys = [ ( o[1], o ) for o in overlaps.keys() ]
        keys.sort()
        keys.reverse()
        keys = [ k[1] for k in keys ]

    # only 1 key? nothing to linearize
    if len(keys) == 1: return ( overlaps, removed_overlaps )

    # check if start_with_bests is properly given
    if start_with_bests not in range(1,6):
        # take default (1) start_with_bests
        start_with_bests = 1

    # get lowest absolute and non-zero bitscore for weighting
    LOWEST_WEIGTH = max( [ 1, float(keys[-1][0]) ] )


    # list to store the Linearization Calibration Overlaps into
    lco = []

    # take the best overlap(s) to calibrate
    # the linearization proces
    # however, check if we are not about to create
    # a nonsense linearization upon adding the n-th element
    for key in keys[0:start_with_bests]:
        # get offset/weigth data and set update lco
        ( offset, weigth ) = _getoverlapdata( key, overlaps[key], LOWEST_WEIGTH, is_weigthed=is_weigthed )
        lco.append( ( offset, weigth ) )
        linear_overlaps[key] = overlaps[key]

    # now loop over the remaining keys
    for key in keys[start_with_bests:]:
        # get Linear Calibration Offset of current overlap
        (offset, weigth ) = _getoverlapdata( key, overlaps[key], LOWEST_WEIGTH, is_weigthed=is_weigthed )
        # get range of accepted Linear Calibration Offsets
        if offset in _get_lco_range(lco,ACCEPTANCE_MARGIN):
            # check the ordering towards the other pacbps
            is_removed = False
            for accepted in linear_overlaps.values():
                order = pacb.ordering.relatively_positioned_towards(accepted,overlaps[key])
                q1,q2,s1,s2 = order['Q1'], order['Q2'], order['S1'], order['S2']
                q1b, q2b = pacb.ordering.relpos2binaryrelpos(q1,q2)
                s1b, s2b = pacb.ordering.relpos2binaryrelpos(s1,s2)
                q2br = pacb.ordering.reverse_binaryrelpos(q2b)
                s2br = pacb.ordering.reverse_binaryrelpos(s2b)
                # mind taking the correct variables here: q1b and s2b , s1b and q2b
                if pacb.ordering.is_identical_binaryrelpos(q1b,s2b) and pacb.ordering.is_identical_binaryrelpos(s1b,q2b):
                    # pacbps in REVERSE order compared to previously accepted pacbp
                    removed_overlaps[key] = overlaps[key]
                    is_removed = True
                    break
                # mind taking the correct variables here: q1b and q1b2 , s1b and s2br
                elif pacb.ordering.is_identical_binaryrelpos(q1b,q2br) and pacb.ordering.is_identical_binaryrelpos(s1b,s2br):
                    # normal order, Q << S or Q >> S
                    pass
                # mind taking the correct variables here: q1b and q1b2 , s1b and s2br
                elif pacb.ordering.is_included_binaryrelpos(q1b,q2br) and pacb.ordering.is_included_binaryrelpos(s1b,s2br):
                    # Q includes S or S includes Q
                    pass
                else:
                    # a more messy relation, not quickly resolved
                    # likely, most but NOT ALL of these, are false
                    pass

            if not is_removed:
                # yep, this overlap is accepted!
                # save it and set it to lco
                linear_overlaps[key] = overlaps[key]
                lco.append( ( offset, weigth ) )
        else:
            # nope, this overlap/pacb is not co-linear
            # with the already accepted overlaps/pacbs
            removed_overlaps[key] = overlaps[key]


    # and return the linearized pacbs and the removed ones
    return ( linear_overlaps, removed_overlaps )

# end of function linearise_pacbs


def _getoverlapdata(key,pacbp,LOWEST_WEIGTH,is_weigthed=True):
    """
    """
    q1,q2  = pacbp.query_start, pacbp.query_end 
    s1,s2  = pacbp.sbjct_start, pacbp.sbjct_end 
    ll     = (q2-q1+s2-s1)/2
    offset = ( q1-s1 + q2-s2 ) / 2
    if is_weigthed: weigth = float(key[0]) / LOWEST_WEIGTH
    else:           weigth = 1.0
    # return offset and weigth
    return (offset, weigth)

# end of function _getoverlapdata

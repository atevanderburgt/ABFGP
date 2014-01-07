"""
Functions for stopless3n introns used in Alignment Based Gene Prediction
"""

# Module metadata
__authors__ = "Ate van der Burgt"
__license__ = "MIT"


def _filter_overlapping_stopless3n_introns(intronlist):
    """
    Filter a list of stopless3n introns on >1x presence of donor/acceptor sites

    @type  intronlist: list
    @param intronlist: list with stopless3n IntronConnectingOrfs

    @attention: make shure only stopless3n IntronConnectingOrfs are applied
    @attention: make shure intron.assign_bp_and_ppts() has been called

    @rtype  intronlist: list
    @return intronlist: list with filtered stopless3n introns
    """
    donor_site_filter = [ (intronlist[pos].acceptor.pos,intronlist[pos].donor.pos ) for pos in range(0,len(intronlist)) ]
    donor_site_filter.sort()
    donor_site_filter.reverse()
    to_be_removed_intron_coords = []
    # make shure acceptor sites are uniquely encountered
    while donor_site_filter:
        if len(donor_site_filter) == 1:
            donor_site_filter = []
            break
        if donor_site_filter[0][0] == donor_site_filter[1][0]:
            # remove this (longer) intron with alternative donor
            to_be_removed_intron_coords.append( donor_site_filter[1] )
            donor_site_filter.pop(1)
        else:
            donor_site_filter.pop(0)

    # remove labeled introns
    for (apos,dpos) in to_be_removed_intron_coords:
        for pos in range(0,len(intronlist)):
            if intronlist[pos].coords() == (dpos,apos):
                intronlist.pop(pos)
                break

    acceptor_site_filter = [ (intronlist[pos].donor.pos,intronlist[pos].acceptor.pos ) for pos in range(0,len(intronlist)) ]
    acceptor_site_filter.sort()
    acceptor_site_filter.reverse()
    to_be_removed_intron_coords = []
    # make shure donor sites are uniquely encountered
    while acceptor_site_filter:
        if len(acceptor_site_filter) == 1:
            acceptor_site_filter = []
            break
        if acceptor_site_filter[0][0] == acceptor_site_filter[1][0]:
            # remove the intron with the poorest branchpoint
            intron0 = None
            intron1 = None
            for intron in intronlist:
                if intron.coords() == acceptor_site_filter[0]:
                    intron0 = intron
                    break
            for intron in intronlist:
                if intron.coords() == acceptor_site_filter[1]:
                    intron1 = intron
                    break
            if intron0.branchpoint and not intron1.branchpoint:
                to_be_removed_intron_coords.append( acceptor_site_filter[1] )
                acceptor_site_filter.pop(1)
            elif not intron0.branchpoint and intron1.branchpoint:
                to_be_removed_intron_coords.append( acceptor_site_filter[0] )
                acceptor_site_filter.pop(0)
            elif not intron0.branchpoint and not intron1.branchpoint:
                # hmmm no branchpoints applied. We expect them here..
                # just remove the longest one
                to_be_removed_intron_coords.append( acceptor_site_filter[0] )
                acceptor_site_filter.pop(0)
            else:
                score0 = intron0.branchpoint._gff['fscore']  
                score1 = intron1.branchpoint._gff['fscore']
                # measure spacing of branchpoint to acceptor site
                dist0  = abs(20-(intron0.acceptor.pos-intron0.branchpoint.end))
                dist1  = abs(20-(intron1.acceptor.pos-intron1.branchpoint.end))
                if score0 > score1:
                    to_be_removed_intron_coords.append( acceptor_site_filter[1] )
                    acceptor_site_filter.pop(1)
                elif score1 > score0:
                    to_be_removed_intron_coords.append( acceptor_site_filter[0] )
                    acceptor_site_filter.pop(0)
                elif dist0 > dist1:
                    to_be_removed_intron_coords.append( acceptor_site_filter[0] )
                    acceptor_site_filter.pop(0)
                elif dist1 > dist0:
                    to_be_removed_intron_coords.append( acceptor_site_filter[1] )
                    acceptor_site_filter.pop(1)
                else:
                    # same score & dist -> delete based on length
                    # just remove the longest one
                    to_be_removed_intron_coords.append( acceptor_site_filter[0] )
                    acceptor_site_filter.pop(0)
        else:
            acceptor_site_filter.pop(0)

    # remove labeled introns
    for (dpos,apos) in to_be_removed_intron_coords:
        for pos in range(0,len(intronlist)):
            if intronlist[pos].coords() == (dpos,apos):
                intronlist.pop(pos)
                break

    # return list of filtered stopless3n introns
    return intronlist

# end of function _filter_overlapping_stopless3n_introns
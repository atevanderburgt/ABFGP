"""
PacbPORF connection by (enforced) bridgeing with intron(s)
"""

# Module metadata
__authors__ = "Ate van der Burgt"
__license__ = "MIT"

# Imports
from pacb.validators import IsPacbPORF
from pacb.ordering import order_list_by_attribute as olba
from pacb.connecting.orfs import merge_orfs_with_intron
from pacb.connecting.functions import (
    _update_kwargs,
    set_apps_intron_query,
    set_apps_intron_sbjct,
    )

# Python Imports
from sets import Set

# Global variable Imports
from settings.splicesites import (
    KWARGS_MAPPED_INTRON,
)


def merge_pacbporfs_with_query_intron_bridgeing(pacbporfD,pacbporfA,verbose=False,**kwargs):
    """
    Merge query Orfs in PacbPORF by **best** intron

    @attention: see orfs.merge_orfs_with_intron for **kwargs

    @type  pacbporfD: PacbPORF object
    @param pacbporfD: PacbPORF object that has to deliver PSSM donor objects

    @type  pacbporfA: PacbPORF object
    @param pacbporfA: PacbPORF object that has to deliver PSSM acceptor objects

    @type  verbose: Boolean
    @param verbose: print status/debugging messages to STDOUT

    @rtype:  list
    @return: list with ( intron, intron ), in query and sbjct
    """
    # input validation
    IsPacbPORF(pacbporfD)
    IsPacbPORF(pacbporfA)

    # edit **kwargs dictionary for some forced attributes
    _update_kwargs(kwargs,KWARGS_MAPPED_INTRON)
    if not kwargs.has_key('aligned_site_max_triplet_distance'):
        kwargs['aligned_site_max_triplet_distance'] = kwargs['max_aa_offset']

    # calculate maximal/minimal donor/acceptor site position based on alignment
    ELEGIABLE_SPLICE_SITE_AA_RANGE = 75

    qdr = pacbporfD.alignment_dna_range_query()
    qar = pacbporfA.alignment_dna_range_query()
    min_donor_query_pos = max([ min(qdr), max(qdr)-(ELEGIABLE_SPLICE_SITE_AA_RANGE*3) ])
    max_accep_query_pos = min([ max(qar), min(qar)+(ELEGIABLE_SPLICE_SITE_AA_RANGE*3) ])

    # get list of introns
    intronlist = merge_orfs_with_intron(pacbporfD.orfQ,pacbporfA.orfQ,
            min_donor_pos   =min_donor_query_pos,
            max_acceptor_pos=max_accep_query_pos,**kwargs)


    # filter on entropy
    # settings for minimal alignment entropy score
    if min([pacbporfD.identityscore,pacbporfA.identityscore]) > 0.55:
        min_donor_site_entropy = 0.01
        min_acceptor_site_entropy = 0.01
        intronlist = _filter_introns_on_entropy(intronlist,pacbporfD,pacbporfA,
                min_donor_site_entropy=min_donor_site_entropy,
                min_acceptor_site_entropy=min_acceptor_site_entropy)
    else:
        # do not filter, but do not forget to store apps data to intron(s)
        for intron in intronlist:
            succes = set_apps_intron_query(intron,pacbporfD,pacbporfA)


    for intron in intronlist:
        intron._distance = 0 # ??
        # set GFF fsource attribute for recognition of intron sources
        intron._gff['fsource'] = 'ABGPbridgeing'

    # get unique list of donors & acceptors
    donor = olba( list(Set([intron.donor for intron in intronlist ])), order_by='pos')
    accep = olba( list(Set([intron.acceptor for intron in intronlist ])), order_by='pos')

    ############################################################################
    if verbose: print "dQ1",[d.pos for d in donor],"aQ1",[a.pos for a in accep]
    ############################################################################

    intronlist = _filter_introns_on_pssm_entropy_combination(intronlist)

    # get unique list of donors & acceptors
    donor = olba( list(Set([intron.donor for intron in intronlist ])), order_by='pos')
    accep = olba( list(Set([intron.acceptor for intron in intronlist ])), order_by='pos')

    ############################################################################
    if verbose: print "dQ1",[d.pos for d in donor],"aQ1",[a.pos for a in accep]
    ############################################################################

    filtered_intron_list = []
    for intron in intronlist:
        intron.assign_bp_and_ppts()
        if intron.branchpoint and (intron.ppt5p or intron.ppt3p):
            filtered_intron_list.append( intron )
        else:
            pass

    # check if list is emptied due to branchpoint filtering
    # in that case, filter for either branchpoint OR polyppt
    if not filtered_intron_list and intronlist:
        for intron in intronlist:
            if intron.branchpoint or (intron.ppt5p or intron.ppt3p):
                filtered_intron_list.append( intron )

    # return list of filtered introns
    return filtered_intron_list

# end of function merge_pacbporfs_with_query_intron_bridgeing



def _filter_introns_on_entropy(intronlist,pacbporfD,pacbporfA,
    min_donor_site_entropy=0.0,min_acceptor_site_entropy=0.0,verbose=False):
    """ """
    remove_list_index_ids = []
    # accelerator. Do no re-obtain entropy when already omitted on this position
    remove_list_donor_positions = []
    remove_list_accep_positions = []
    for pos in range(0,len(intronlist)):
        intron = intronlist[pos]
        if intron.donor.pos in remove_list_donor_positions:
            remove_list_index_ids.append(pos)
            continue
        if intron.acceptor.pos in remove_list_accep_positions:
            remove_list_index_ids.append(pos)
            continue
        is_accepted = True
        succes = set_apps_intron_query(intron,pacbporfD,pacbporfA)
        if intron._apps_donor < min_donor_site_entropy:
            is_accepted = False
            remove_list_donor_positions.append( intron.donor.pos )
        if intron._apps_accep < min_acceptor_site_entropy:
            is_accepted = False
            remove_list_accep_positions.append( intron.acceptor.pos )
        if not is_accepted:
            remove_list_index_ids.append(pos)

    # remove aligned introns
    remove_list_index_ids.reverse()
    for pos in remove_list_index_ids: intronlist.pop(pos)
    # return intronlist
    return intronlist

# end of function _filter_introns_on_entropy


def _filter_introns_on_pssm_entropy_combination(intronlist,verbose=False):
    """
    """
    # return empty list if no intronlist was given
    if not intronlist: return intronlist

    # calculate summed PSSM score
    for pos in range(0,len(intronlist)):
        intron = intronlist[pos]
        summed_pssm = intron.donor.pssm_score + intron.acceptor.pssm_score
        intronlist[pos] = (summed_pssm,intron)
    # order by summed_pssm
    intronlist.sort()
    intronlist.reverse()
    best_pssm_score,bestintron = intronlist[0]
    best_donor_entropy = bestintron._apps_donor
    best_accep_entropy = bestintron._apps_accep

    # select introns to be removed from the list of candidates
    # based on lower summed PSSM AND/OR entropy <= 0.0
    remove_list_index_ids = []
    for pos in range(1,len(intronlist)):
        summed_pssm,intron = intronlist[pos]
        if summed_pssm < (best_pssm_score * 0.8):
            if intron._apps_donor <= 0.0 or intron._apps_accep <= 0.0:
                # introns have both poor PSSM and poor entropy
                remove_list_index_ids.append(pos)
        elif summed_pssm < best_pssm_score:
            if best_donor_entropy > 0.0 and intron._apps_donor <= 0.0:
                # introns have okay PSSM but (very) poor donor entropy
                remove_list_index_ids.append(pos)
            elif best_accep_entropy > 0.0 and intron._apps_accep <= 0.0:
                # introns have okay PSSM but (very) poor acceptor entropy
                remove_list_index_ids.append(pos)
            else:
                pass
        else:
            pass

    # remove aligned introns
    remove_list_index_ids.reverse()
    for pos in remove_list_index_ids: intronlist.pop(pos)
    # translate list back to list of introns only (remove summed_pssm)
    intronlist = [ intron for (summed_pssm,intron) in intronlist ]

    # return intronlist
    return intronlist

# end of function _filter_introns_on_pssm_entropy_combination

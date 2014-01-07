"""
Helper functions for PacPORF connecting
"""
# Module metadata
__authors__ = "Ate van der Burgt"
__license__ = "MIT"

# Pacb Imports
from pacb.exceptions import CoordinateOutOfRange

# Other Imports

# Python Imports

# Global variable Imports
from settings.splicesites import KWARGS_STOPLESS_3N_INTRONS


def _update_kwargs(kwargs,defaults):
    """
    Update **kwargs dict with defaults for non-existing keys
    """
    for k,v in defaults.iteritems():
        if not kwargs.has_key(k):
            kwargs[k] = v
# end of function_update_kwargs


def _filter_for_entropy(sites,pacbporf,method,min_alignment_entropy=0.0):
    """
    Helper function for merge_pacbporfs_with_introns
    """
    # remove sites with to low alignment entropy
    tobedeleted = []
    for sQ,sS in sites:
        # calculate based on query. THIS IS AN ASSUMPTION FOR THE SBJCT
        dpp, phase = pacbporf.dnaposition_query(sQ.pos,forced_return=True)
        algEntropy = pacbporf.alignment_entropy(dpp,method=method)
        if algEntropy < min_alignment_entropy:
            tobedeleted.append( ( sQ,sS ) )
    # remove & return
    for ( sQ,sS ) in tobedeleted: sites.remove( ( sQ,sS ) )
    return sites

# end of function _filter_for_entropy


def _tinyexon_list_2_dict(resultlist):
    """
    Helper function to re-structure data returned by merge_orfs_with_tinyexon()
    """
    resultdict = {}
    key2exon = {}
    for (intronD,tinyexon,intronA) in resultlist:
        exon_key = (tinyexon.orf.id, tinyexon.donor.pos, tinyexon.acceptor.pos)
        if not resultdict.has_key(exon_key):
            resultdict[exon_key] = [ {}, {} ]
            key2exon[exon_key] = tinyexon
        intronDkey = (intronD.orfDonor.id,intronD.donor.pos)
        intronAkey = (intronA.orfAcceptor.id,intronA.acceptor.pos)
        if not resultdict[exon_key][0].has_key(intronDkey):
            resultdict[exon_key][0][intronDkey] = intronD
        if not resultdict[exon_key][1].has_key(intronAkey):
            resultdict[exon_key][1][intronAkey] = intronA
    return resultdict, key2exon
# end of function _tinyexon_list_2_dict


def _score_introns_obtained_by_mapping(intQ,intS,pacbporfD,pacbporfA,source="ABGPmapping"):
    """
    Helper function to score introns obtained by mapping
    """
    distDnt = pacbporfD.get_distance_aligned_nucleotide_positions(
                    query = intQ.donor.pos, sbjct = intS.donor.pos
                    )
    distAnt = pacbporfA.get_distance_aligned_nucleotide_positions(
                    query = intQ.acceptor.pos, sbjct = intS.acceptor.pos
                    )

    # add distance score to introns
    intQ._distance = abs(distDnt) + abs(distAnt)
    intS._distance = abs(distDnt) + abs(distAnt)

    # add Alignment Positional Periphery Score into objects
    succes = set_apps_intron_query(intQ,pacbporfD,pacbporfA)
    succes = set_apps_intron_sbjct(intS,pacbporfD,pacbporfA)

    # set GFF fsource attribute for recognition of intron sources
    intQ._gff['fsource'] = source
    intS._gff['fsource'] = source

# end of function _score_introns_obtained_by_mapping


def _filter_for_alignable_splice_sites(sitesQ,sitesS,pacbporf,
    allow_non_canonical = False,
    allow_phase_shift = False,
    aligned_site_max_triplet_distance = 1,   # ... 4 CBG.collectionharvesting
    **kwargs):
    """
    Helper function for merge_pacbporfs_with_introns
    """
    # filter for alignable splice sites
    sites = []
    for siteQ in sitesQ:
        if not allow_non_canonical and siteQ.is_canonical() == False:
            continue
        for siteS in sitesS:
            if not allow_non_canonical and siteS.is_canonical() == False:
                continue
            if not allow_phase_shift and siteQ.phase != siteS.phase:
                continue


            try:
                # calculate the distance in aligned nt positions
                dist = pacbporf.get_distance_aligned_nucleotide_positions(
                        query = siteQ.pos, sbjct = siteS.pos
                        )
            except CoordinateOutOfRange:
                continue

            # check if sites are alignable
            if dist > aligned_site_max_triplet_distance*3:
                # sites way to far apart, not alignable
                continue
            elif not allow_phase_shift and siteQ.phase != siteS.phase:
                # allow_phase_shift = False, incompatible phases 
               continue
            elif dist <= aligned_site_max_triplet_distance*3:
                pass
            else:
                # all other cases -> failed
                continue


            # if here, alignable splice sites!
            # Final check if one of the sites is already listed as alignable
            # If so, store only the *best* alignable site!
            if siteQ.pos in [ sQ.pos for sQ,sS,d in sites ] or\
            siteS.pos in [ sS.pos for sQ,sS,d in sites ]:
                donotaddnew = False
                replacements = []
                for sQ,sS,d in sites:
                    if (sQ.pos == siteQ.pos or sS.pos == siteS.pos):
                        if d == 0:
                            # already perfectly aligned combi -> omit replacements
                            replacements = []
                            donotaddnew = True
                            break
                        elif dist == 0:
                            # new combi isa perfect one -> omit current listed
                            replacements.append( ( sQ,sS,d ) )
                        else:
                            # both current as existing one are imperfect.
                            # no obvious decision possible -> leave in place
                            pass
                # remove the replacement(s)
                for (sQ,sS,d) in replacements: sites.remove( (sQ,sS,d) )
                # add the novel one if not disallowed
                if not donotaddnew: sites.append( ( siteQ, siteS, dist ) )
            else:
                # append alignable sites & distance
                sites.append( ( siteQ, siteS, dist ) )

    # return list with sites, omit the distances
    return [ (sQ,sS) for (sQ,sS,d) in sites ]

# end of function _filter_for_alignable_splice_sites


def _filter_sites_on_pssm(sitelist,min_pssm_score=None):
    """ """
    if min_pssm_score or min_pssm_score==0.0:
        return_sitelist = []
        for pos in range(0,len(sitelist)):
            site = sitelist[pos]
            if site.pssm_score < min_pssm_score:
                pass
            else:
                return_sitelist.append( site )
        # return list of sites
        return return_sitelist
    else:
        # no filtering value set
        return sitelist

# end of function _filter_sites_on_pssm


def _filter_aligned_sites_on_pssm(alignedsitelist,min_pssm_score=None):
    """ """
    if min_pssm_score:
        remove_list_index_ids = []
        for pos in range(0,len(alignedsitelist)):
            siteQ,siteS = alignedsitelist[pos]
            if siteQ.pssm_score < min_pssm_score:
                remove_list_index_ids.append(pos)
                continue
            if siteS.pssm_score < min_pssm_score:
                remove_list_index_ids.append(pos)
                continue
        # remove sites
        remove_list_index_ids.reverse()
        for pos in remove_list_index_ids: alignedsitelist.pop(pos)

    # return list of aligned sites
    return alignedsitelist

# end of function _filter_aligned_sites_on_pssm


def _filter_stopless_3n_introns(intronlist,verbose=False,**kwargs):
    """ """
    # stopless3n introns are (easily) overpredicted in
    # Ab Initio Gene Prediction, and in ABFGP it is as well the case
    # Therefor, stopless 3n introns have stricter splice site concensus
    # criterion. Apply a list of introns, and the stopless3n introns
    # will be filtered on some **kwargs

    # update **kwargs dictionary with forced attributes
    _update_kwargs(kwargs,KWARGS_STOPLESS_3N_INTRONS)
    returnlist = []
    for intron in intronlist:
        if intron.__class__.__name__ == "SequenceErrorConnectingOrfs":
            returnlist.append( intron )
            continue
        if not intron.is_stopless_3n_intron():
            returnlist.append( intron )
            continue
        # if here, start filtering stopless3nintron
        if intron.donor.is_canonical():
            if intron.donor.pssm_score < kwargs['min_donor_pssm_score']:
                continue
        else:
            if intron.donor.pssm_score < kwargs['non_canonical_min_donor_pssm_score']:
                continue
        if intron.acceptor.is_canonical():
            if intron.acceptor.pssm_score < kwargs['min_acceptor_pssm_score']:
                continue
        else:
            if intron.acceptor.pssm_score < kwargs['non_canonical_min_acceptor_pssm_score']:
                continue
        if intron.length < kwargs['min_intron_nt_length']:
            continue
        if intron.length > kwargs['max_intron_nt_length']:
            continue

        # if here: an accepted stopless3n intron!
        returnlist.append( intron )
        
    # return list of accepted introns
    return returnlist

# end of function _filter_stopless_3n_introns


def set_apps_intron_query(intron,pacbporfD,pacbporfA):
    """ Shortcut function for setting both APPS in the intron object """
    # add Alignment Positional Periphery Score into objects
    intron._apps_donor = get_apps_donor_query(intron.donor,pacbporfD)
    intron._apps_accep = get_apps_acceptor_query(intron.acceptor,pacbporfA)
    return True # succesfully set attributes in intron object

# end of function set_apps_intron_query


def set_apps_intron_sbjct(intron,pacbporfD,pacbporfA):
    """ Shortcut function for setting both APPS in the intron object """
    intron._apps_donor = get_apps_donor_sbjct(intron.donor,pacbporfD)
    intron._apps_accep = get_apps_acceptor_sbjct(intron.acceptor,pacbporfA)
    return True # succesfully set attributes in intron object

# end of function set_apps_intron_sbjct


def get_apps_donor_query(donor,pacbporf):
    """ Get Alignment Positional Periphery Score for donor in the query """
    return _get_apps_query(donor,pacbporf,"donor")

# end of function get_apps_donor_query


def get_apps_donor_sbjct(donor,pacbporf):
    """ Get Alignment Positional Periphery Score for donor in the sbjct """
    return _get_apps_sbjct(donor,pacbporf,"donor")

# end of function get_apps_donor_sbjct


def get_apps_acceptor_query(acceptor,pacbporf):
    """ Get Alignment Positional Periphery Score for acceptor in the query """
    return _get_apps_query(acceptor,pacbporf,"acceptor")

# end of function get_apps_acceptor_query


def get_apps_acceptor_sbjct(acceptor,pacbporf):
    """ Get Alignment Positional Periphery Score for acceptor in the sbjct """
    return _get_apps_sbjct(acceptor,pacbporf,"acceptor")

# end of function get_apps_acceptor_sbjct


def _get_apps_query(obj,pacbporf,method):
    """ Get Alignment Positional Periphery Score ob object in the query """
    posObj, phase = pacbporf.dnaposition_query(obj.pos,forced_return=True)
    entropy       = pacbporf.alignment_entropy(posObj,method=method)
    return entropy

# end of function _get_apps_query


def _get_apps_sbjct(obj,pacbporf,method):
    """ Get Alignment Positional Periphery Score ob object in the sbjct """
    posObj, phase = pacbporf.dnaposition_sbjct(obj.pos,forced_return=True)
    entropy       = pacbporf.alignment_entropy(posObj,method=method)
    return entropy

# end of function _get_apps_sbjct

"""
Functions performed on crossdata dictionary used in ABFGP

In many of these functions three different predefined dictionary structures
are used as input or output variables:

- crossdata

- pacbpsdict

- input


"""

# Module metadata
__authors__ = "Ate van der Burgt"
__license__ = "MIT"


# import ABFGP modules, classes and functions
from graphAbgp import GeneTreeGraph
from lib_sequenceperiodicity import coding_sequence_periodicity_test
import pacb.recombination
import pacb.conversion

# import linearization settings
from settings.pacbp import (
    LINEARIZATION_STARTWITHBEST,
    LINEARIZATION_ACCEPTANCE_MARGIN,
    LINEARIZATION_WEIGTHED,
    ALTERNATIVE_ALIGNMENT_OVERLAP_RATIO,
    REPETITIVE_ALIGNMENT_OVERLAP_RATIO,
    )

# Python Imports
from copy import deepcopy
from sets import Set

# dict for each of the main orgA-orgB keys in crossdata
CROSSDATA_SUBDICT = {
        'accepted_pacbs'    : {},
        'codingblocks'      : {},   # pacbps in CBGS -> NOT 100% covered all through the code!
        'lowscoring_pacbs'  : {},
        'rejected_pacbs_nl' : {},   # non linear
        'rejected_pacbs_nc' : {},   # non-coding
        'rejected_pacbs_wp' : {},   # wrong periodicty
        'rejected_pacbs_ls' : {},   # weaker similarity out-of-frame
        'rejected_pacbs_aa' : {},   # alternative alignment
        'rejected_pacbs_lc' : {},   # low connectivity in graph
        'best_central_pacbp': None,
        }


def createcrossdata(input):
    """
    Create crossdata dict structure from input data structure

    @rtype  crossdata: dict
    @return crossdata: empty crossdata <dict data structure>
    """
    # create crossdata data structure
    crossdata = {}
    for (geneQ,geneS) in pacb.recombination.pairwise(input.keys()):
        crossdata[(geneQ,geneS)] = deepcopy(CROSSDATA_SUBDICT)
        
    # return the crossdata data dictionary
    return crossdata

# end of function createcrossdata


def count_pacbps_in_crossdata(crossdata,keys=['accepted_pacbs']):
    """
    Count total number of Pacbps in crossdata for given key(s)
    
    @type  crossdata: dict
    @param crossdata: crossdata <dict data structure>

    @type  keys: list
    @param ksys: list with string(s) from CROSSDATA_SUBDICT.keys()

    @attention: empty keys list counts all CROSSDATA_SUBDICT.keys()!

    @rtype:  integer
    @return: toal number of pacbps
    """
    cnt = 0
    if not keys: keys = CROSSDATA_SUBDICT.keys()
    for (orgA,orgB) in crossdata.keys():
        for key in keys:
            if crossdata[(orgA,orgB)].has_key(key):
                cnt+=len(crossdata[(orgA,orgB)][key])
            else:
                # key not recognized. Maybe raise error?
                pass
    # return total number
    return cnt

# end of function count_pacbps_in_crossdata


def create_pacbpcollectiongraph_from_crossdata(
    crossdata,
    increment_to_graph=None):
    """
    Create (incremental) PacbpCollectionGraph from crossdata dict structure

    @type  crossdata: dict
    @param crossdata: crossdata <dict data structure>

    @type  increment_to_graph: None or PacbpCollectionGraph
    @param increment_to_graph: when applied, increment crossdata to applied PCG

    @rtype  g: PacbpCollectionGraph
    @return g: PacbpCollectionGraph
    """
    if increment_to_graph:
        # add new nodes to existing graph
        g = increment_to_graph
    else:
        # create a new blank graph
        from graphAbgp import PacbpCollectionGraph
        g = PacbpCollectionGraph()

    for (orgA,orgB) in crossdata.keys():
        keys = crossdata[(orgA,orgB)]['accepted_pacbs'].keys()
        # sort keys in order to start with highest bitscore
        keys.sort()
        keys.reverse()
        for key in keys:
            (bitscore,lenght,pointerA,pointerB) = key
            nodeA = (orgA,pointerA)
            nodeB = (orgB,pointerB)
            # check if (org,ORF) node exist already
            if nodeA not in g.get_nodes(): g.add_node(nodeA)
            if nodeB not in g.get_nodes(): g.add_node(nodeB)
            if g.has_edge(nodeA,nodeB):
                wt = g.get_edge_weight(nodeA,nodeB)
                if bitscore > wt:
                    g.set_edge_weight(nodeA,nodeB,bitscore)
                else:
                    pass
            else:
                # and create a new edge
                g.add_edge(nodeA,nodeB,wt=bitscore)
    # ready!
    return g

# end of function create_pacbpcollectiongraph_from_crossdata


def crossdata2pacbpcollectiongraph(crossdata,**kwargs):
    """
    Create (incremental) PacbpCollectionGraph from crossdata dict structure

    @attention: alias name for create_pacbpcollectiongraph_from_crossdata
    """
    return create_pacbpcollectiongraph_from_crossdata(crossdata,**kwargs)
    
# end of function crossdata2pacbpcollectiongraph


def create_genetree_from_crossdata(crossdata):
    """
    (Try to) create the GeneTreeGraph from crossdata dictionairy

    @type  crossdata: dict
    @param crossdata: crossdata <dict data structure>

    @rtype:  GeneTreeGraph
    @return: estimated! GeneTreeGraph constructed from PabcPs with highest bitscore
    """
    # fill GeneTreeGraph with nodes
    GTG = GeneTreeGraph()
    for (geneA,geneB) in crossdata.keys():
        if geneA not in GTG.get_nodes():
            GTG.add_node(geneA)
        if geneB not in GTG.get_nodes():
            GTG.add_node(geneB)

    # fill GeneTreeGraph with edges
    for (geneA,geneB) in crossdata.keys():
        keys = crossdata[(geneA,geneB)]['accepted_pacbs'].keys()
        if keys:
            keys.sort()
            keys.reverse()
            bestpacbp = crossdata[(geneA,geneB)]['accepted_pacbs'][keys[0]]
            # store this edge to the GeneTreeGraph
            GTG.add_edge(geneA,geneB,bestpacbp.identityscore)
        else:
            # no keys at all, meaning no Pacbps between 2 species meaining
            # GTG can not be created yet! Set GSG back to empty graph and return
            # set GTG back to an empty graph, break out anf return
            GTG = GeneTreeGraph()
            break

    # return the GeneTreeGraph
    return GTG

# end of function create_genetree_from_crossdata


def crossdata2genetreegraph(crossdata):
    """
    (Try to) create the GeneTreeGraph from crossdata dictionairy

    @attention: alias name for create_genetree_from_crossdata
    """
    return create_genetree_from_crossdata(crossdata)
    
# end of function crossdata2genetreegraph


def remove_nonlinear_pacbs(crossdata,
    startwithbest=LINEARIZATION_STARTWITHBEST,
    acceptance_margin=LINEARIZATION_ACCEPTANCE_MARGIN,
    is_weigthed=LINEARIZATION_WEIGTHED):
    """
    Remove PacbP(ORF)s that break a colinear pairwise alignment

    @type  crossdata: dict
    @param crossdata: crossdata <dict data structure>

    @type  startwithbest: integer
    @param startwithbest: start with this number of highest scoring Pacbs

    @type  acceptance_margin: integer
    @param acceptance_margin: offset (in AA coordinates) for overlap used
                              during the linearization proces

    @type  is_weigthed: Boolean
    @param is_weigthed: weight linearization by bitscore (default True)

    @rtype  crossdata: dict
    @return crossdata: crossdata <dict data structure>
    """
    for (geneQ,geneS) in crossdata.keys():
        pacbs = crossdata[(geneQ,geneS)]['accepted_pacbs']
        # filter out all non-linear pacbs
        ( accepted, nonlinear ) = pacb.linearization.linearise_pacbs(
                pacbs,
                ACCEPTANCE_MARGIN=acceptance_margin,
                start_with_bests=startwithbest,
                is_weigthed=is_weigthed
                )
        # update crossdata dictionaries with
        # accepted pacbps & rejected nonlinear pacbps
        crossdata[(geneQ,geneS)]['accepted_pacbs'] = accepted
        crossdata[(geneQ,geneS)]['rejected_pacbs_nl'].update(nonlinear)

    # and return crossdata
    return crossdata

# end of function remove_nonlinear_pacbs


def _update_pacbp_dict_by_accepted_list(pacbp_dict,pacbp_list):
    """ Update *THE* pacbp_dict by a list of ACCEPTED pacbps """
    if not pacbp_list:
        # ALL pacbps are rejected!
        return ( {}, pacbp_dict )
    else:
        accepted_pacbp_strings = [ str(pacbp) for pacbp in pacbp_list ]
        removed_pacbps = {}
        # set removed ones to removed_pacbps
        for key, pacbp in pacbp_dict.iteritems():
            if str(pacbp) in accepted_pacbp_strings:
                pass
            else:
                removed_pacbps[key] = pacbp
        # delete all keys from pacbp_dict that are now in removed_pacbps
        for key in removed_pacbps.keys():
            del( pacbp_dict[key] )
        # and return the updated accepted pacbps and the removed ones
        return ( pacbp_dict, removed_pacbps )

# end of function _update_pacbp_dict_by_accepted_list


def _update_pacbp_dict_by_rejected_list(pacbp_dict,pacbp_list):
    """ Update *THE* pacbp_dict by a list of REJECTED pacbps """
    if not pacbp_list:
        # ALL pacbps are accepted, none are rejected!
        return ( pacbp_dict, {} )
    else:
        rejected_pacbp_strings = [ str(pacbp) for pacbp in pacbp_list ]
        removed_pacbps = {}
        # set removed ones to removed_pacbps
        for key, pacbp in pacbp_dict.iteritems():
            if str(pacbp) in rejected_pacbp_strings:
                removed_pacbps[key] = pacbp

            else:
                pass
        # delete all keys from pacbp_dict that are now in removed_pacbps
        for key in removed_pacbps.keys():
            del( pacbp_dict[key] )
        # and return the updated accepted pacbps and the removed ones
        return ( pacbp_dict, removed_pacbps )

# end of function _update_pacbp_dict_by_rejected_list


def remove_inclusive_pacbps(crossdata,query=True,sbjct=False,both=False,verbose=False):
    """
    Remove PacbP(ORF)s that are fully included in another PacbP(ORF)
    """
    # IMPORTANT: leave sbjct=False!!
    # this will remove overlapping QUERY pacbps, but leave overlapping
    # SBJCT pacbps intact -> in the pairwise mode we are interested in
    # a global QUERY alignment, sbjct doesn't matter. If sbjct = True,
    # cases like CFU_828161 will easily go wrong. This is a gene with 4
    # point mutations, causing the PacbPs to overlap a lot.
    for (geneQ,geneS) in crossdata.keys():
        rejected_keys = Set()
        for keyA in crossdata[(geneQ,geneS)]['accepted_pacbs'].keys():
            for keyB in crossdata[(geneQ,geneS)]['accepted_pacbs'].keys():
                if keyA == keyB: continue
                pacbpA = crossdata[(geneQ,geneS)]['accepted_pacbs'][keyA]
                pacbpB = crossdata[(geneQ,geneS)]['accepted_pacbs'][keyB]
                pdict  = pacbpA.relatively_positioned_towards(pacbpB)
                incQ1  = (pdict['Q1'][0],pdict['Q1'][2])
                incQ2  = (pdict['Q2'][0],pdict['Q2'][2])
                incS1  = (pdict['S1'][0],pdict['S1'][2])
                incS2  = (pdict['S2'][0],pdict['S2'][2])
                orfQid = keyA[2] == keyB[2]
                orfSid = keyA[3] == keyB[3]
                incs   = [incQ1,incQ2,incS1,incS2]
                if (0,0) not in incs:
                    continue
                elif incQ1 == (0,0) and incS1 == (0,0) and pacbpA.bitscore < pacbpB.bitscore:
                    if both: rejected_keys.add(keyA)
                elif incQ2 == (0,0) and incS2 == (0,0) and pacbpB.bitscore < pacbpA.bitscore:
                    if both: rejected_keys.add(keyB)
                elif incQ1 == (0,0) and query and False in (orfQid,orfSid) and\
                pacbpA.bitscore < pacbpB.bitscore:
                    rejected_keys.add(keyA)
                    ############################################################
                    if verbose: print "## inclusive:",pacbpA,incs,keyA,"->",keyB
                    ############################################################
                elif incQ2 == (0,0) and query and False in (orfQid,orfSid) and\
                pacbpB.bitscore < pacbpA.bitscore:
                    rejected_keys.add(keyB)
                    ############################################################
                    if verbose: print "## inclusive:",pacbpB,incs,keyB,"->",keyA
                    ############################################################
                elif incS1 == (0,0) and sbjct and False in (orfQid,orfSid) and\
                pacbpA.bitscore < pacbpB.bitscore:
                    rejected_keys.add(keyA)
                    ############################################################
                    if verbose: print "## inclusive:",pacbpA,incs,keyA,"->",keyB
                    ############################################################
                elif incS2 == (0,0) and sbjct and False in (orfQid,orfSid) and\
                pacbpB.bitscore < pacbpA.bitscore:
                    rejected_keys.add(keyB)
                    ############################################################
                    if verbose: print "## inclusive:",pacbpB,incs,keyB,"->",keyA
                    ############################################################
                else:
                    # no valid PacbP inclusion detected
                    pass

        # check if there are rejected pacbp keys
        if not rejected_keys: continue

        pacbplist = crossdata[(geneQ,geneS)]['accepted_pacbs']
        accepted, rejected = _update_pacbp_dict_by_rejected_list(
                pacbplist, [ crossdata[(geneQ,geneS)]['accepted_pacbs'][k] for k in rejected_keys ])

        # update crossdata dictionaries with
        # accepted pacbps & rejected nonlinear pacbps
        crossdata[(geneQ,geneS)]['accepted_pacbs'] = accepted
        crossdata[(geneQ,geneS)]['rejected_pacbs_aa'].update(rejected)

    # and return crossdata
    return crossdata

# end of function remove_inclusive_pacbps


def remove_repetitive_pacbps(crossdata,
    overlap_ratio=REPETITIVE_ALIGNMENT_OVERLAP_RATIO):
    """
    Remove PacbP(ORF)s caused by repetitiveness and as such breaking the pairwise alignment

    @type  crossdata: dict
    @param crossdata: crossdata <dict data structure>

    @type  overlap_ratio: float
    @param overlap_ratio: max overlap ratio between PacbPs

    @rtype  crossdata: dict
    @return crossdata: crossdata <dict data structure>
    """
    for (geneQ,geneS) in crossdata.keys():
        pacbps = crossdata[(geneQ,geneS)]['accepted_pacbs']
        # filter out all repetitive (and non-linear) pacbs
        accepted_list = pacb.linearization.remove_repeats_from_pacbp_list(
                pacbps.values(), overlap_ratio = overlap_ratio
                )
        
        accepted, rejected = _update_pacbp_dict_by_accepted_list(pacbps,accepted_list)

        # update crossdata dictionaries with
        # accepted pacbps & rejected nonlinear pacbps
        crossdata[(geneQ,geneS)]['accepted_pacbs'] = accepted
        crossdata[(geneQ,geneS)]['rejected_pacbs_aa'].update(rejected)

    # and return crossdata
    return crossdata

# end of function remove_repetitive_pacbps


def remove_wrong_periodicity(input,crossdata):
    """
    Remove PacbPs with wrong periodicity from crossdata

    @type  input: dict
    @param input: input <dict data structure> with lists of Orfs

    @type  crossdata: dict
    @param crossdata: crossdata <dict data structure>
    """
    for (geneQ,geneS) in crossdata.keys():
        keys = crossdata[(geneQ,geneS)]['accepted_pacbs'].keys()
        for key in keys:
            pacbp = crossdata[(geneQ,geneS)]['accepted_pacbs'][key]
            if pacbp.__class__.__name__ == 'PacbPORF':
                # pacbp is already a PacbPORF
                pacbporf = pacbp
            else:
                (bitscore,length,orfQid,orfSid) = key
                orfQ     = input[geneQ]['orfs'].get_orf_by_id(orfQid)
                orfS     = input[geneS]['orfs'].get_orf_by_id(orfSid)
                pacbporf = pacb.conversion.pacbp2pacbporf(pacbp,orfQ,orfS)

            # get aligned DNA sequences from the pacbporf
            (dnaQ,dnaS) = pacbporf.get_aligned_dna_sequences()
            # do coding_sequence_periodicity_test
            (_test1,_test2) = coding_sequence_periodicity_test(dnaQ,dnaS)
            if (_test1,_test2) != (True, True):
                # store to rejected
                crossdata[(geneQ,geneS)]['rejected_pacbs_wp'][key] = pacbporf
                # remove from accepted_pacbps dict
                del( crossdata[(geneQ,geneS)]['accepted_pacbs'][key] )

    # okay, return the updated crossdata
    return crossdata

# end of function remove_wrong_periodicity


def recover_rejected_pacbps(subgraph,crossdata):
    """
    Recover rejected pacbps that seem relevant (as missing edges) in the graph.

    @attention: make sure crossdata dictionary has correct data structure!

    @type  subgraph: PacbpCollectionGrapg or CodingBlockGraph
    @param subgraph: PacbpCollectionGrapg or CodingBlockGraph

    @type  crossdata: {}
    @param crossdata: crossdata dictionary with predefined structure

    @rtype:  Graph object
    @return: Graph object with additional edges recruited from crossdata

    @rtype:  {}
    @return: crossdata dictionary with predefined structure
    """

    # only perform this function if graph is NOT fully connected
    if subgraph.connectivitysaturation() == 1.0: return subgraph

    rejected_keys = [
        'accepted_pacbs',       # yes! as well check here for abandoned edges!!
        'lowscoring_pacbs',
        'rejected_pacbs_wp',    # wrong periodicity
        'rejected_pacbs_ls',    # weaker similarity out-of-frame
        'rejected_pacbs_lc'     # low connectivity ??? TODO is not filled yet!!
        #'rejected_pacbs_nl',   # non-linear
        #'rejected_pacbs_aa',   # alternative alignment

    ]
    # loop over the organism/gene combinations
    for (geneQ,geneS) in crossdata.keys():
        crossdata[(geneQ,geneS)]['potential_pacbs'] = {}
        known_combis = {}
        # make a lookup dictionary with all orf-2-orf combinations
        # this is needed lateron to check if the newly recruited
        # edge is compatible with the already present edges
        for key in crossdata[(geneQ,geneS)]['accepted_pacbs'].keys():
            # `key` is a 4-elem-tuple ( bitscore, length, orf_id_1, orf_id_2 )
            if (geneQ,key[2]) in subgraph.get_nodes() and\
	    (geneS,key[3]) in subgraph.get_nodes():
                known_combis[( key[2], key[3] )] = key

        # (1) loop over all the rejected_key names,
        # (2) then loop over all the individual rejected pacbps,
        # (3) then check if both orfs of the pacbp are present as nodes
	#     in the graph
        # (4) then check if the egde is compatible with already present edges
        for rej_key in rejected_keys:
            # (2) then loop over all the individual rejected pacbps,
            for (a,b,orfQ,orfS) in crossdata[(geneQ,geneS)][rej_key].keys():
                # (3) then check if both orfs of the pacbp are
		#     present as nodes in the graph
                if (geneQ,orfQ) in subgraph.get_nodes() and\
		(geneS,orfS) in subgraph.get_nodes():
                    # yep, be have a potential (missing) edge for this graph!
                    potential_new_pacbp =\
			crossdata[(geneQ,geneS)][rej_key][(a,b,orfQ,orfS)]
                    if rej_key == 'accepted_pacbs':
                        if not subgraph.has_edge( (geneQ,orfQ), (geneS,orfS) ):
                            # add edge to graph and append to crossdata
                            subgraph.add_edge(
				(geneQ,orfQ), (geneS,orfS),
				wt=potential_new_pacbp.bitscore
				)
                        # and continue to the next case
                        continue

                    # (4) then check if the egde is compatible with
		    #     already present edges
                    compatibility = []
                    # loop over the already existing edges gathered in
                    # the lookup dictionary `known_combis` and check for
                    # compatibility. At this point, do NOT check upon pacbps
                    # of this specific edge. If a false pacbp should be
                    # assigned to this edge, re-recruitment of the true
		    # pacbp is nearly impossible
                    for fullkey in known_combis.values():
                        orig_pacbp =\
			    crossdata[(geneQ,geneS)]['accepted_pacbs'][fullkey]
                        if potential_new_pacbp.bitscore > orig_pacbp.bitscore:
                            # it has a higher bitscore; do not compare
                            # its relative position
                            continue
                        # get relative positioning of potential new pacbp
			# with other pacbps (==edges in the graph)
                        relpos =\
			    potential_new_pacbp.relatively_positioned_towards(
				orig_pacbp
				)
                        # check compatibility of the position (True or False)
                        compatibility.append(
			    _is_compatible_pacbp_positioning(relpos)
			    )

                    # now check the compatibility; a single False is a no-go!
                    if False in compatibility:
                        # nope, not a compatible (new) edge. Continue here
                        continue

		    # When here, a compatible edge. IMPORTANT! check if the edge
		    # is already present in the graph. That seems weird,
		    # but is possible! 2 orfs can be connected by >1 pacbp
		    # (e.g. within-orf intron!). In a previous step,
		    # only the highest scoring pacbp is added to the graph
                    # as an edge and loaded in
		    # crossdata[(geneQ,geneS)]['accepted_pacbs'].
		    # Here, the lower scoring, alternative alignment is
		    # assigned and stored to the crossdata dict as well.
                    if not subgraph.has_edge( (geneQ,orfQ), (geneS,orfS) ):
                        # add edge to graph and append to crossdata
                        subgraph.add_edge( (geneQ,orfQ), (geneS,orfS),
			        wt=potential_new_pacbp.bitscore
			        )
                    # append to `accepted_pacbps` if not already there
                    crossdata[(geneQ,geneS)]['accepted_pacbs'][(a,b,orfQ,orfS)] =\
			    potential_new_pacbp

    # and return!
    return subgraph, crossdata

# end of function recover_rejected_pacbps


def split_accepted_pacbps_on_gapsize(crossdata,gapsize=100):
    """
    Split pacbps on (large) gaps.

    @type  crossdata: dict
    @param crossdata: crossdata <dict data structure>

    @type  gapsize: integer
    @param gapsize: minimal AA gap length in a PacbP in order to be splitted

    @rtype  crossdata: dict
    @return crossdata: crossdata <dict data structure>
    """     
    for (geneQ,geneS) in crossdata.keys():
        # try to perform the splits on the dict `accepted_pacbps`
        (new_pacbps, new_keys, old_keys) = _split_pacbps_on_gapsize(
                crossdata[(geneQ,geneS)]['accepted_pacbs'],
                gapsize=gapsize
                )
        # replace old dict with new dict
        crossdata[(geneQ,geneS)]['accepted_pacbs'] = new_pacbps
    
    # return the updated version of crossdata dict
    return crossdata

# end of function split_accepted_pacbps_on_gapsize


def split_lowscoring_pacbps_on_gapsize(crossdata,gapsize=100):
    """ 
    Split lowscoring pacbps on (large) gaps.

    @type  crossdata: dict
    @param crossdata: crossdata <dict data structure>

    @type  gapsize: integer
    @param gapsize: minimal AA gap length in a PacbP in order to be splitted

    @rtype  crossdata: dict
    @return crossdata: crossdata <dict data structure>
    """     
    for (geneQ,geneS) in crossdata.keys():
        # try to perform the splits on the dict `accepted_pacbps`
       
        (new_pacbps, new_keys, old_keys) = _split_pacbps_on_gapsize(
                crossdata[(geneQ,geneS)]['lowscoring_pacbs'],
                gapsize=gapsize
                )
        # Set all that are splitted to 'accepted_pacbs'
        # This was most likely an inframe intron that is now resolved.
        # Maybe one should want to check the scores of these splitted
        # parts, but at this stage it is believed that it is oke.
        for key in new_keys:
            # store to 'accepted_pacbs'
            crossdata[(geneQ,geneS)]['accepted_pacbs'][key] = new_pacbps[key]
            # delete from the newly splitted dict
            del( new_pacbps[key] )

        # replace old dict with new dict
        crossdata[(geneQ,geneS)]['lowscoring_pacbs'] = new_pacbps

    # return the updated version of crossdata dict
    return crossdata

# end of function split_lowscoring_pacbps_on_gapsize


def _split_pacbps_on_gapsize(pacbpdict,gapsize=100):
    """
    Split individual pacbps in a dict on (large) gaps.

    @type  pacbpdict: dict
    @param pacbpdict: dict with keys (a,b,c,d) and PacbP objects as values

    @type  gapsize: integer (positive)
    @param gapsize: continious gap to be split the PacbPs by (in AA length)
    """
    splitted_keys = []
    new_accepted  = {}
    # loop over dict and try to perform the splits
    for key,thepacbp in pacbpdict.iteritems():
        # split the pacbp if there is a too long gapsize
        splitted_pacbps, is_splitted = pacb.splitting.split_pacb_on_gaps(
                thepacbp,gapsize=gapsize)
        if is_splitted:
            # add these to new_accepted
            for pacbp in splitted_pacbps:
                new_key = ( pacbp.bits, pacbp.length, key[2], key[3] )
                new_accepted[new_key] = pacbp
            # and list `old` key for removal
            splitted_keys.append(key)

    # now check if splits have been performed
    if splitted_keys:
        # remove the old pacbps
        for key in splitted_keys: del( pacbpdict[key] )
        # add the new splitted pacbps
        pacbpdict.update( new_accepted )
    
    # return the updated dict of pacbps and list of new & old keys
    return pacbpdict, new_accepted.keys(), splitted_keys
    
# end of function _split_pacbps_on_gapsize


def recover_lowscoring_pacbps(crossdata,GTG,MATRIX):
    """
    based on identityscore and maximal attainable bitscore,
    recover rejected PacbPs.
    PacbPs have been rejected based on BLASTP_MIN_BITSCORE
    With GTG, we have a gene tree with approxination
    of protein identity. Based on this, take a (considerable)
    number of lower scoring PacbPs into account as well.
    """
    cnt_total=0
    if not (GTG and GTG.get_nodes() ):
        # a GTG is not applied yet, so this function is not callable.
        return crossdata

    for (geneQ,geneS) in crossdata.keys():
        # TODO TODO!! these values must become a ratio of the GTG / genetree
        # In the -!- single -!- example checked (example10, problematic first exon)
        # the lowest ratio was 0.41 (ratio1) and 0.23 (ratio2). However, this was
        # only for the WEAKEST connected edge.
        # The other 3 were 0.62 (ratio1) and 0.42 (ratio2).
        # This weakest edge was still quit strong; in the Organism graph 84%
        # Identities of OrganismGraph
        #, ('fgsg', 'fveg'), 0.936170212766
        #, ('mgg', 'fgsg'), 0.857638888889
        #, ('mgg', 'ncu'), 0.836805555556
        #, ('ncu', 'fveg'), 0.855614973262
        #, ('ncu', 'foxg'), 0.860962566845
        #, ('foxg', 'mgg'), 0.840277777778
        #, ('fgsg', 'ncu'), 0.842105263158
        #, ('foxg', 'fgsg'), 0.913265306122
        #, ('fveg', 'foxg'), 0.99481865285
        #, ('mgg', 'fveg'), 0.840277777778
        # relative weight of OrganismGraph
        #, ('fgsg', 'mgg'), 0.0977039242058
        #, ('ncu', 'fveg'), 0.0974733557211
        #, ('fveg', 'fgsg'), 0.106650368467
        #, ('ncu', 'foxg'), 0.0980825641944
        #, ('foxg', 'mgg'), 0.0957261119749
        #, ('ncu', 'mgg'), 0.0953305495287
        #, ('fgsg', 'ncu'), 0.095934302736
        #, ('fgsg', 'foxg'), 0.104040996048
        #, ('fveg', 'foxg'), 0.113331715149
        #, ('mgg', 'fveg'), 0.0957261119749
        # So, look further into weighting these numbers!!

        # For this example, when the complete blast is performed with BLOSUM45,
        # a vast total of >3000..>6000 additional hits are gathered, causing the
        # algorithm to get fully stuck.
        # A smart idea is to include a 3th ratio:
        # ratio1:   bit / maxbit
        # ratio2:  (bit / maxbit) * indentityscore
        # ratio3: ((bit / maxbit) * indentityscore) * length
        # For the example of example10 above, this would result in:
        # 0.41 - 0.23 and 3.68 for the weakest one   (rank ~700/3500)
        # 0.62 - 0.42 and 5.88 for the three others  (rank  ~50/3500)
        # TODO: Have a closer look on this!
        #BLASTP_RECOVERY_RATIO_1 = 0.375 * GTG.weights[(geneQ,geneS)]
        #BLASTP_RECOVERY_RATIO_2 = 0.200 * GTG.weights[(geneQ,geneS)]

        identityscore = GTG.get_edge_weight(geneQ,geneS)

        # create a recovery ratio based on the GTG of this pair of geneQ/geneS
        BLASTP_RECOVERY_RATIO = 7.0 * GTG.weights[(geneQ,geneS)] * GTG.weights[(geneQ,geneS)]

        recovered_keys = []
        for key,pacbp in crossdata[(geneQ,geneS)]['lowscoring_pacbs'].iteritems():
            maxscore = MATRIX.scorealignment(pacbp.query,pacbp.query)
            bitscore = MATRIX.scorealignment(pacbp.query,pacbp.sbjct)
            #ratio1   = float(bitscore)/float(maxscore)
            #ratio2   = ratio1*pacbp.identityscore

            try:
                ratio   = ( float(bitscore)/float(maxscore) ) * pacbp.identityscore * float(pacbp.length)
                if ratio < BLASTP_RECOVERY_RATIO or bitscore <= 0:
                    continue
            except ZeroDivisionError:
                # maxscore variable == 0.0; do not evaluate this pacbp at all
                continue

            # if here, then recover this low scoring pacbp
            crossdata[(geneQ,geneS)]['accepted_pacbs'][key] = pacbp
            # and mark as to be removed from lowscoring
            recovered_keys.append(key)
            cnt_total+=1

        # and delete from lowscoring
        for key in recovered_keys:
            del( crossdata[(geneQ,geneS)]['lowscoring_pacbs'][key] )

    # and return
    return crossdata

# end of function recover_lowscoring_pacbps


def remove_alternative_alignments(crossdata,
    keep_only_highest_scoring=False,
    overlap_ratio=ALTERNATIVE_ALIGNMENT_OVERLAP_RATIO):
    """
    Remove PacbP objects that overlap (to much) with each other

    @type  crossdata: dict
    @param crossdata: crossdata <dict data structure>

    @type  overlap_ratio: float
    @param overlap_ratio: max overlap ratio between PacbPs

    @type  keep_only_highest_scoring: Boolean
    @param keep_only_highest_scoring: when True and alternative PacbP
            alignment(s) occur, remove all except the highest bitscore PacbP

    @attention: overlap_ratio >> 0.0 and < 1.0

    @rtype  crossdata: dict
    @return crossdata: crossdata <dict data structure>
    """
    for (geneQ,geneS) in crossdata.keys():
        # remove alternative alignments for this Organism identifier pair
        (accepted, rejected) = remove_alternatives_from_pacbps_dict(
                crossdata[(geneQ,geneS)]['accepted_pacbs'],
                keep_only_highest_scoring=keep_only_highest_scoring,
                overlap_ratio=overlap_ratio)

        # replace accepted_pacbps
        crossdata[(geneQ,geneS)]['accepted_pacbs'] = accepted

        # update the rejected ones to crossdata
        crossdata[(geneQ,geneS)]['rejected_pacbs_aa'].update( rejected )

    # return crossdata dict
    return crossdata

# end of function remove_alternative_alignments


def TOBECHECKED_WRONG_PACBPDICT_remove_alternatives_from_pacbps_dict(pacbpsdict,
    keep_only_highest_scoring=False,
    overlap_ratio=ALTERNATIVE_ALIGNMENT_OVERLAP_RATIO):
    """
    Remove PacbP objects that overlap (to much) with each other

    @attention: helper function of remove_alternative_alignments()

    @type  pacbpsdict: {}
    @param pacbpsdict: pacbpsdict <dict data structure>

    @type  keep_only_highest_scoring: Boolean
    @param keep_only_highest_scoring: when True and alternative PacbP
            alignment(s) occur, remove all except the highest bitscore PacbP

    @type  overlap_ratio: float
    @param overlap_ratio: max overlap ratio between PacbPs

    @attention: overlap_ratio >> 0.0 and < 1.0

    @rtype  pacbpsdict: dictionary
    @return pacbpsdict: cleaned-up pacbpsdict dictionary

    @rtype  rejected: ( {}, {} )
    @return rejected: cleanedup and rejected pacbpsdict dictionary
    """
    keys = pacbpsdict.keys()
    keys.sort()
    keys.reverse()
    rejected = {}
    # Capital A/B for arbitrarily first/second pacbp
    # a,b,c,d represents the 4-elem-tuple of a pacbp (bitscore,length,orf1,orf2)
    # g1,o1 represents gene/organism and orf of gene/organisms 1
    # g2,o2 represents gene/organism and orf of gene/organisms 2
    # make a full cross of all vs. all pacbps keys in the dict
    for ( (Aa,Ab,Ac,Ad),(Ag1,Ao1),(Ag2,Ao2) ) in keys:
        for ( (Ba,Bb,Bc,Bd),(Bg1,Bo1),(Bg2,Bo2) ) in keys:
            # ignore identical keys
            if (Aa,Ab,Ac,Ad) == (Ba,Bb,Bc,Bd): continue

            # only do upper triangle of square matrix
            if Aa < Ba: continue

            # ignore alignments between non-identical orf sets
            if ((Ag1,Ao1),(Ag2,Ao2)) != ((Bg1,Bo1),(Bg2,Bo2)): continue

            # yep, identical orf sets; get the PacbP objects
            pacbp1 = pacbpsdict[( (Aa,Ab,Ac,Ad),(Ag1,Ao1),(Ag2,Ao2) )]
            pacbp2 = pacbpsdict[( (Ba,Bb,Bc,Bd),(Bg1,Bo1),(Bg2,Bo2) )]

            # check if argument keep_only_highest_scoring is True
            if keep_only_highest_scoring:
                # just discard the 2th pacbp
                rejected[( (Ba,Bb,Bc,Bd),(Bg1,Bo1),(Bg2,Bo2) )] = pacbp2
                # done here; continue to next comparison
                continue

            # check relative position of PacbP objects
            relpos = pacbp1.relatively_positioned_towards(pacbp2)

            # calculate overlap based on pacbp with lowest bitscore
            overlap_q = float( relpos['Q2'][1] ) / float(sum(relpos['Q2']))
            overlap_s = float( relpos['S2'][1] ) / float(sum(relpos['S2']))

            # check if overlap is sufficiently large to remove
            if max([overlap_q,overlap_s]) > overlap_ratio:
                rejected[( (Ba,Bb,Bc,Bd),(Bg1,Bo1),(Bg2,Bo2) )] = pacbp2

    # now separate the alternative alignments
    # from the original overlaps
    for deletethiskey,value in rejected.iteritems():
        del( pacbpsdict[deletethiskey] )

    # and return the (new) pacbpsdict and rejected
    return pacbpsdict, rejected

# end of function remove_alternatives_from_pacbps_dict


def remove_alternatives_from_pacbps_dict(pacbpsdict,
    keep_only_highest_scoring=False,
    overlap_ratio=ALTERNATIVE_ALIGNMENT_OVERLAP_RATIO):
    """
    Remove PacbP objects that overlap (to much) with each other

    @attention: helper function of remove_alternative_alignments()

    @type  pacbpsdict: {}
    @param pacbpsdict: pacbpsdict <dict data structure>

    @type  keep_only_highest_scoring: Boolean
    @param keep_only_highest_scoring: when True and alternative PacbP
            alignment(s) occur, remove all except the highest bitscore PacbP

    @type  overlap_ratio: float
    @param overlap_ratio: max overlap ratio between PacbPs

    @attention: overlap_ratio >> 0.0 and < 1.0

    @rtype  pacbpsdict: dictionary
    @return pacbpsdict: cleaned-up pacbpsdict dictionary

    @rtype  rejected: ( {}, {} )
    @return rejected: cleanedup and rejected pacbpsdict dictionary
    """
    keys = pacbpsdict.keys()
    keys.sort()
    keys.reverse()
    rejected = {}
    # Capital A/B for arbitrarily first/second pacbp
    # a,b,c,d represents the 4-elem-tuple of a pacbp (bitscore,length,orf1,orf2)
    # make a full cross of all vs. all pacbps keys in the dict
    for (Aa,Ab,Ac,Ad) in keys:
        for (Ba,Bb,Bc,Bd) in keys:
            # ignore identical keys
            if (Aa,Ab,Ac,Ad) == (Ba,Bb,Bc,Bd): continue

            # only do upper triangle of square matrix
            if Aa < Ba: continue

            # ignore alignments between non-identical orf sets
            if (Ac,Ad) != (Bc,Bd): continue

            # yep, identical orf sets; get the PacbP objects
            pacbp1 = pacbpsdict[ (Aa,Ab,Ac,Ad) ]
            pacbp2 = pacbpsdict[ (Ba,Bb,Bc,Bd) ]

            # check if argument keep_only_highest_scoring is True
            if keep_only_highest_scoring:
                # just discard the 2th pacbp
                rejected[ (Ba,Bb,Bc,Bd) ] = pacbp2
                # done here; continue to next comparison
                continue

            # check relative position of PacbP objects
            relpos = pacbp1.relatively_positioned_towards(pacbp2)

            # calculate overlap based on pacbp with lowest bitscore
            overlap_q = float( relpos['Q2'][1] ) / float(sum(relpos['Q2']))
            overlap_s = float( relpos['S2'][1] ) / float(sum(relpos['S2']))

            # check if overlap is sufficiently large to remove
            if max([overlap_q,overlap_s]) > overlap_ratio:
                rejected[ (Ba,Bb,Bc,Bd) ] = pacbp2

    # now separate the alternative alignments
    # from the original overlaps
    for deletethiskey,value in rejected.iteritems():
        del( pacbpsdict[deletethiskey] )  

    # and return the (new) pacbpsdict and rejected
    return pacbpsdict, rejected

# end of function remove_alternatives_from_pacbps_dict


def remove_lower_similarities_in_other_frame(crossdata):
    """
    """
    for (geneA,geneB) in crossdata.keys():
        # get currently accepted pacbs
        pacbps = crossdata[(geneA,geneB)]['accepted_pacbs']
        accepted, rejected = separateoverlaps(pacbps,
            ratio_cutoff=WEAKER_SIMILARITY_OTHER_FRAME_RATIO)
        # and set back to crossdata
        crossdata[(geneA,geneB)]['accepted_pacbs'] = accepted
        crossdata[(geneA,geneB)]['rejected_pacbs_ls'].update( rejected )

    # return crossdata
    return crossdata

# end of function remove_lower_similarities_in_other_frame


def separateoverlaps(overlaps,ratio_cutoff=2.0):
    """
    Separate list of PacbPs into two fractions
        overlaps (accepted)
        overlaps (rejected) with stronger similaties in another frame

    An (graphical) example is this:
        hsp-A1[0]-B1[2]   +++++++++++++++++++++++++++   bitscore 1000
        hsp-A2[1]-B2[1]             ====                bitscore 45
        sliceout hsp                ++++                bitscore 150

     The first, long and high scoring PacbP (hsp-A1[0]-B1[0]) is a hit between
        orf A1[0] (organism A, orf 1, frame 0)
        orf B1[0] (organism B, orf 1, frame 2)
    The second, short and low scoring PacbP (hsp-A2[1]-B2[1]) is a hit between
        orf A2[1] (organism A, orf 2, frame 1)
        orf B2[1] (organism B, orf 2, frame 1)
    This second PacbP is just an accidential hit.
    These PacbPs should be discarded based on periodicity analyses too, but in
    exceptional hases (very low evolutionary distances),
    this might be not the case.
    Higher similarity in another frame is proven by slicing out the fraction
    of the first PacbP that overlaps with the first (based on nt coordinates).
    For this partial PacbP, the bitscore is calculated.
    When this calculated bitscore is ``ratio_cutoff`` times higher than
    the bitscore of the second, short PacbP, this second hsp is placed
    in the list of rejected overlaps with stronger similarity in another frame.
    This means: discarded from further analyses!!
    """

    weakersimilarities = {}
    # get overlap keys and sort
    keys = overlaps.keys()
    keys.sort()
    keys.reverse()

    # now loop and find orfoverlappingHsp's that are most
    # likely aligned in the wrong frame
    for key1 in keys:
        #(bitscore, overlaplength, _orfpointerA,_orfpointerB) = key1
        #thiscorrectedhsp = overlaps[(bitscore, overlaplength, _orfpointerA,_orfpointerB)]
        for key2 in keys:
            #(bitscore2, overlaplength2, _orfpointerA2,_orfpointerB2) = key2
            #thiscorrectedhsp2 = overlaps[(bitscore2, overlaplength2, _orfpointerA2,_orfpointerB2)]
            # continue on identical keys
            if key1 == key2: continue
            # check if PacbP 2 is fully included by PacbP 1
            pacbp1 = overlaps[key1]
            pacbp2 = overlaps[key2]
            relpos = pacbp1.relatively_positioned_towards(pacbp2)
            check = ( relpos['Q2'][0], relpos['Q2'][2], relpos['S2'][0], relpos['S2'][2] )
            # check if PacbP 2 is fully included by PacbP 1
            if check != (0,0,0,0): continue

            # if here, then PacbP 2 is fully included by PacbP 1
            try:
                otherframeQ = pacbp1.returnslice(pacbp2.query_start,pacbp2.query_end,coords_on="query")
                if otherframeQ.bitscore == 0:
                    ratioQ = 0.0
                else:
                    ratioQ = float(pacbp2.bitscore) / float(otherframeQ.bitscore)
            except pacb.ZeroSizedPacb:
                otherframeQ = "<ZeroSizedPacb>"
                ratioQ = 0.0
            except:
                print pacbp1
                print pacbp2
                print relpos
                print check
                otherframeQ = pacbp1.returnslice(pacbp2.query_start,pacbp2.query_end,coords_on="query")
                ratioQ = float(pacbp2.bitscore) / float(otherframeQ.bitscore)
                raise "UNEXPECTED query slicing ERROR OCCURRED!"
            try:
                otherframeS = pacbp1.returnslice(pacbp2.sbjct_start,pacbp2.sbjct_end,coords_on="sbjct")
                if otherframeS.bitscore == 0:
                    ratioS = 0.0
                else:
                    ratioS = float(pacbp2.bitscore) / float(otherframeS.bitscore)
            except pacb.ZeroSizedPacb:
                otherframeS = "<ZeroSizedPacb>"
                ratioS = 0.0
            except:
                print pacbp1
                print pacbp2
                print relpos
                print check
                otherframeS = pacbp1.returnslice(pacbp2.sbjct_start,pacbp2.sbjct_end,coords_on="sbjct")
                ratioS = float(pacbp2.bitscore) / float(otherframeS.bitscore)
                print otherframeS
                print ratioS
                raise "UNEXPECTED sbjct slicing ERROR OCCURRED!"

            if max([ratioQ,ratioS]) >= ratio_cutoff:
                # yep, this most likely is an of-frame accidential alignment
                # the otherframe alignment-slice is not far less scoring as the fully included pacbp
                weakersimilarities[key2] = overlaps[key2]

    # now separate the weakersimilarities
    # from the original overlaps
    for deletethiskey,value in weakersimilarities.iteritems():
        del(overlaps[deletethiskey])

    # and return
    return ( overlaps, weakersimilarities )

# end of function separateoverlaps


def _is_compatible_pacbp_positioning(p):
    """
    
    @rtype:  Boolean
    @return: True or False depending if the relative positioning is compatible
    """
    bp = {} 
    for k,tup in p.iteritems():
        bp[k] = []
        for digit in tup:
            if digit >= 1:  bp[k].append(1)
            else:           bp[k].append(0)
    # reverse the 2th keys
    bp['S2'].reverse()
    bp['Q2'].reverse()
    if len(Set([tuple(l) for l in bp.values() ] ) )==1:
        return True
    else:
        # correct for potential overlap
        for k in bp.keys(): bp[k][1] = 0
        if len(Set([tuple(l) for l in bp.values() ] ) )==1:
            return True
        else:
            return False

# end of function _is_compatible_pacbp_positioning



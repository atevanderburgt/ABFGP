"""
Functions required for Intron Projection in a CodingBlockGraphInterface
"""

# ABGP Imports
from graphPlus.comparison import mutual_nodes
# import intron projection functions from pacb.connecting
from pacb.connecting import\
    project_splicesites_on_pacbporfs_with_lacking_intron_in_sbjct as PrjIntLsbj
from pacb.connecting import\
    project_splicesites_on_pacbporfs_with_lacking_intron_in_query as PrjIntLqry

# Global variables Imports
from settings.inframeintron import INFRAME_INTRON_MIN_AA_LENGTH

def project_missing_introns(cbgs,
    interface_intron_status={},
    inframe_intron_min_aa_length=INFRAME_INTRON_MIN_AA_LENGTH,
    verbose=False):
    """
    Create ProjectedIntrons for intron gain/loss events in ordered CBGs

    @type  cbgs: list
    @param cbgs: ordered (!) list of CodingBlockGraphs

    @type  interface_intron_status: dict
    @param interface_intron_status:

    @type  inframe_intron_min_aa_length: integer
    @param inframe_intron_min_aa_length: minimal AA length of an inframe intron

    @type  verbose: Boolean
    @param verbose: print logging/debugging messages to STDOUT when True

    @rtype:  dict
    @return: dict with organism as keys, list of Projected Introns as values
    """
    # fill dictionary projected_introns with keys of each org identifier
    return_projected_introns= {}
    for node in cbgs[0].get_nodes():
        return_projected_introns[cbgs[0]._organism_from_node(node)] = []

    # loop over all other coding blocks
    # and evaluate the junction between `prev` and `this`
    for position in range(1,len(cbgs)):
        this = cbgs[position]
        prev = cbgs[position-1]

        # find intersection nodes between graphs
        intersectnodes = mutual_nodes(prev,this)

        if intersectnodes:
            # yep, there are shared ORFs between the adjacent CBGs
            # this means we can try to project introns
            intersectorgs =\
                [ prev._organism_from_node(node) for node in intersectnodes ]

            # check for identical pacbps shared by both CBGs ->
            # no intron projection here!
            nodes_seen,max_nodes_seen = cbgs_identical_pacbp_analysis(prev,this)

            ####################################################################
            if verbose:
                print "nodes:", nodes_seen
                print "inters", max_nodes_seen, "intersectorgs:", intersectorgs
            ####################################################################

            for orgQ,orgS in prev.pairwisecrosscombinations_organism():

                # check interface_intron_status dict. Do only a projection
                # when orgQ status in (True,None) and orgS status in
                # (False,None) or vice versa
                if interface_intron_status:
                    if interface_intron_status[orgQ] in (True,None) and\
                    interface_intron_status[orgS] in (False,None):
                        pass
                    elif interface_intron_status[orgS] in (True,None) and\
                    interface_intron_status[orgQ] in (False,None):
                        pass
                    else:
                        # no aparantly projection required in this interface!
                        continue
                else:
                    # Interface_intron_status not applied or empty.
                    # Just pass here, will be solved lateron
                    pass

                # ignore by continue if orgQ or orgS identifier is absent
                # in this or prev CBG; this happens when it is a K(s-x) CBG
                if orgQ not in prev.organism_set(): continue
                if orgQ not in this.organism_set(): continue
                if orgS not in prev.organism_set(): continue
                if orgS not in this.organism_set(): continue

                # get the corresponding nodes
                nodeQprev = prev.node_by_organism(orgQ)
                nodeSprev = prev.node_by_organism(orgS)
                nodeQthis = this.node_by_organism(orgQ)
                nodeSthis = this.node_by_organism(orgS)
                
                # check if edge exists in both CBGs (K(s) CBG, no edges missing)
                if not prev.has_edge(nodeQprev,nodeSprev) or\
                not this.has_edge(nodeQthis,nodeSthis):
                    continue
                # if both organisms are not part of the intersection ->
                # 2x non-identical Orfs so, no projection possible here
                if not orgQ in intersectorgs and not orgS in intersectorgs:
                    continue

                # calculate distance between codingblocks
                distQ = prev.distance_between_codingblocks(this,organism=orgQ)
                distS = prev.distance_between_codingblocks(this,organism=orgS)

                # check for potential inframe intron; if found,
                # ignore the intron projection here
                if orgQ in intersectorgs and orgS in intersectorgs:
                    if nodes_seen[nodeQprev] >= max_nodes_seen[nodeQprev]-1 and\
                    nodes_seen[nodeSprev] >= max_nodes_seen[nodeSprev]-1:
                        # The Pacbp of this edge is identical in `prev`
                        # and `this` CBG. Now check the distance beteen the
                        # CBGs (distance between their OMSRs). If this distance
                        # is large or distinct between the 2 organisms/genes,
                        # then intron projection is still required/usefull.
                        # If small, then possibly an inframe intron or a region
                        # of low AA similarity of variable length range.
                        if abs(distQ-distS) < inframe_intron_min_aa_length:
                            # distance discrepancy is fairly small; this
                            # represents a continious alignment.
                            # Do not project introns here!
                            continue

                # get Pacbp of donor/prev CBG of this organism combination
                pfD = prev.get_pacbps_by_nodes(
                            node1=nodeQprev,node2=nodeSprev)[0]
                pfA = this.get_pacbps_by_nodes(
                            node1=nodeQthis,node2=nodeSthis)[0]

                # Now map Projected Introns!
                if nodeQprev == nodeQthis and nodeSprev != nodeSthis:
                    # continious query alignment
                    tmp_projected_introns = PrjIntLqry(pfD,pfA)
                    org_of_projection = orgQ
                elif nodeSprev == nodeSthis and nodeQprev != nodeQthis:
                    # continious sbjct alignment
                    tmp_projected_introns = PrjIntLsbj(pfD,pfA)
                    org_of_projection = orgS
                elif nodeQprev != nodeQthis and nodeSprev != nodeSthis:
                    # this case is recognized a few lines above:
                    # not orgQ in intersectorgs and not orgS in intersectorgs
                    raise "THIS CANNOT HAPPEN!!!"
                else:
                    # intersection nodes/organisms. Several options possible
                    if distS > distQ and nodes_seen[nodeQprev] >=\
                    max_nodes_seen[nodeQprev]-1:
                        # continious query alignment
                        tmp_projected_introns = PrjIntLqry(pfD,pfA)
                        org_of_projection = orgQ
                    elif distQ > distS and nodes_seen[nodeSprev] >=\
                    max_nodes_seen[nodeSprev]-1:
                        # continious sbjct alignment
                        tmp_projected_introns = PrjIntLsbj(pfD,pfA)
                        org_of_projection = orgS
                    else:
                        ########################################################
                        if verbose:
                            # hmm.... seems an inframe introns in two species!?
                            print "WARNING: inframe intron 2 species!!!",
                            print orgQ,orgS, distQ, distS
                            print nodeQprev, nodeQthis,
                            print "seen:", nodes_seen[nodeQprev],
                            print "max:", max_nodes_seen[nodeQthis],
                            print nodeSprev, nodeSthis,
                            print "seen:", nodes_seen[nodeSprev],
                            print "max:", max_nodes_seen[nodeSthis]
                        ########################################################
                        # no mapping of intron possible because
                        # no continious query or sbjct
                        continue

                ################################################################
                if verbose:
                    print orgQ,orgS, "intron count:", len(tmp_projected_introns)
                ################################################################

                # done! now remap/process the `tmp_projected_introns` and
                # store to return_projected_introns dictionary.
                # See if it is a new projected intron, or an intron that
                # coinsides with a previous found one. If it coinsides, it
                # is even stringer proof that this is a valid Projected Intron
                # and as a consequence its score will be increased
                for proj_intron in tmp_projected_introns:
                    current_projected_intron = None

                    # see if there is in current species already an exact
                    # projected intron of this kind; if so,
                    # add_projection_intron() to the earlier found intron.
                    # This will increase the score of this particular intron
                    for known_projected_intron in\
                    return_projected_introns[org_of_projection]:
                        if proj_intron.barcode() ==\
                        known_projected_intron.barcode():
                            known_projected_intron.add_projected_intron(
                                proj_intron.projected_introns[0] )
                            current_projected_intron = known_projected_intron
                            break
                    else:
                        # not seen before; add to return_projected_introns 
                        return_projected_introns[org_of_projection].append(
                            proj_intron )
                        current_projected_intron = proj_intron

                    # Not projected onto the earlier found introns yet,
                    # so do it here by add_projected_intron() to them too.
                    true_intron = proj_intron.projected_introns[0]
                    for _org in return_projected_introns.keys():
                        if _org == org_of_projection: continue
                        for known_projected_intron in\
                        return_projected_introns[_org]:
                            if true_intron.barcode() in\
                            [ i.barcode() for i in\
                            known_projected_intron.projected_introns ]:
                                # append this intron as extra evidence in
                                # this other organism
                                known_projected_intron.add_projected_intron(
                                    true_intron )
                                # append this same intron another time
                                # to this organism
                                current_projected_intron.add_projected_intron(
                                    true_intron )

        ########################################################################
        if verbose:
            print "SUMMED INTRONS PER ORG:", (position-1,position),
            print [ (org,len(return_projected_introns[org])) for org\
                in return_projected_introns.keys() ]
            print interface_intron_status
        ########################################################################

    # sort the projected introns by total_score
    for org in return_projected_introns.keys():
        tmp = []
        for projintron in return_projected_introns[org]:
            projD = projintron.get_projected_donor_site()
            projA = projintron.get_projected_acceptor_site()
            tmp.append( ( projintron.total_score(), projintron ) )
        tmp.sort()
        tmp.reverse()
        return_projected_introns[org] = [ proj for score,proj in tmp ]

    # and return the dictionary of projected_introns
    return return_projected_introns

# end of function project_missing_introns


def projectedintrons2projectedsites(projintrondict):
    """
    Convert the result of project_missing_introns() to dicts of
    ProjectedAcceptors and ProjectedDonors

    @type  verbose: dict
    @param verbose: result dictionary of project_missing_introns()

    @rtype:  tuple
    @return: ( projAcceptors, projDonors ) tuple, projAcceptors and ProjDonors
                dictionaries like projintrondict but with ProjectedSpliceSites
                in stead of ProjectedIntrons
    """
    projected_donors = {}
    projected_acceptors = {}
    for org, proj_intron_list in projintrondict.iteritems():
        projected_donors[org] = {}
        projected_acceptors[org] = {}
        for proj_intron in proj_intron_list:
            projD = proj_intron.get_projected_donor_site()
            projA = proj_intron.get_projected_acceptor_site()
            if projected_donors[org].has_key(projD.pos):
                if projected_donors[org][projD.pos].total_score() <\
                projD.total_score():
                    projected_donors[org][projD.pos] = projD
            else:
                projected_donors[org][projD.pos] = projD
            if projected_acceptors[org].has_key(projA.pos):
                if projected_acceptors[org][projA.pos].total_score() <\
                projA.total_score():
                    projected_acceptors[org][projA.pos] = projA
            else:
                projected_acceptors[org][projA.pos] = projA

    # now get rid of the ``position`` keys,
    # so return just dict with key==organism,
    # value== a list with projected intron positions
    for org, _dict in projected_donors.iteritems():
        projected_donors[org] = _dict.values()
    for org, _dict in projected_acceptors.iteritems():
        projected_acceptors[org] = _dict.values()

    # return projected acceptor sites and projected donor sites
    return projected_acceptors, projected_donors

# end of function projectedintrons2projectedsites


def cbgs_identical_pacbp_analysis(first,second):
    """
    CBG PacbPORF comparision

    @type  first: CodingBlockGraph 
    @param first: CodingBlockGraph 

    @type  second: CodingBlockGraph 
    @param second: CodingBlockGraph 

    @rtype:  {}, {}
    @return: dicts with counts how often pacbps of mutual nodes are identical

    """ 
    # gather list of nodes present in both CBGs
    mutual_node_list = mutual_nodes(first,second)
    # set counter dicts for:
    # nodes_seen    -> nodes with identical PacbPs in both CBGs
    # max_node_seen -> nodes that have edges between them in both CBGs
    nodes_seen    = {}
    max_node_seen = {}
    mutual_pacbp_edges = []
    for node in mutual_node_list:
        nodes_seen[node] = 0
        max_node_seen[node] = 0
    # loop over all node combinations / edges
    for nodeQ,nodeS in first.pairwisecrosscombinations_node():
        # only deal with combinations of mutual nodes
        if nodeQ not in mutual_node_list or nodeS not in mutual_node_list:
            continue
        # get PacbP(ORF) of first (donor) and second (acceptor) CBG
        donorPacbp, accepPacbp = None, None
        if first.has_edge(nodeQ,nodeS):
            donorPacbp = first.get_pacbps_by_nodes(
                    node1=nodeQ,node2=nodeS)[0]
        if second.has_edge(nodeQ,nodeS):
            accepPacbp = second.get_pacbps_by_nodes(
                    node1=nodeQ,node2=nodeS)[0]
        # check if both edges exist -> increase max_node_seen
        if accepPacbp and donorPacbp:
            max_node_seen[nodeQ]+=1
            max_node_seen[nodeS]+=1
            # check if Pacbps are identical -> increase nodes_seen
            if donorPacbp.barcode() == accepPacbp.barcode():
                mutual_pacbp_edges.append( (nodeQ,nodeS) )
                nodes_seen[nodeQ]+=1
                nodes_seen[nodeS]+=1

    # return both dictionaries
    return nodes_seen, max_node_seen

# end of function cbgs_identical_pacbp_analysis


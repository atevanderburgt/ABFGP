################################################################################
#### Helper functions for gathering elegiable Collections of (PSSM) sites   #### 
#### in grpahAbgp CoodingBlockGraph class                                   ####
################################################################################

# graphAbgp imports
from graph_pssmcollections import DonorSiteCollectionGraph, AcceptorSiteCollectionGraph, TranslationalStartSiteCollectionGraph
from exceptions import *

# Python Imports
from copy import deepcopy

# Global Variables
from settings.sitealignment import (
    ELIGABLE_DONOR_SITE_LEFT_OF_OMSR_AA_OFFSET, ELIGABLE_DONOR_SITE_RIGTH_OF_OMSR_AA_OFFSET, 
    ELIGABLE_ACCEPTOR_SITE_LEFT_OF_OMSR_AA_OFFSET, ELIGABLE_ACCEPTOR_SITE_RIGTH_OF_OMSR_AA_OFFSET, 
    ELIGABLE_DONOR_SITE_MINIMAL_AA_OFFSET, ELIGABLE_ACCEPTOR_SITE_MINIMAL_AA_OFFSET,
    ELIGABLE_ALIGNED_TSS_3P_AA_OFFSET, ELIGABLE_ALIGNED_TSS_5P_AA_OFFSET,
    ELIGABLE_ALIGNED_START_SITES_AA_OFFSET, ALIGNED_TSS_MAX_AA_DISTANCE,
    MAX_SPLICE_SITE_PHASE_SHIFT_NT_DISTANCE, MIN_DONOR_SITE_PHASE_SHIFT_PSSM_SCORE,
    MIN_ACCEP_SITE_PHASE_SHIFT_PSSM_SCORE
    )
from settings.translationalstartsites import TSS_MIN_PSSM_SCORE


def harvest_elegiable_donor_sites(self,projected_donors={},forced_codingblock_ends={},next=None,
    store_all_projected_sites=False,
    allow_phase_shift=False,
    enlarge_5p_boundary_by=None, # in AA coordinates
    enlarge_3p_boundary_by=None, # in AA coordinates
    ALIGNED_DONOR_MAX_TRIPLET_DISTANCE=None,
    MIN_DONOR_PSSM_SCORE=None,ALLOW_NON_CANONICAL_DONOR=False,
    NON_CANONICAL_MIN_DONOR_PSSM_SCORE=None ):
    """
    Harvest elegiable donor sites from this CodingBlockGraph into a DonorSiteCollectionGraph
    """

    if next and next.__class__.__name__ not in ["CodingBlockGraph","LowSimilarityRegionCodingBlockGraph"]:
        message = "next must be a CodingBlock graph object, not a %s" % next.__class__.__name__
        raise InproperlyAppliedArgument, message

    # update minimal pssm score to stg collection object
    stg = DonorSiteCollectionGraph()
    stg.MIN_PSSM_SCORE = MIN_DONOR_PSSM_SCORE
    stg.ALIGNED_SITE_AA_OFFSET = ALIGNED_DONOR_MAX_TRIPLET_DISTANCE

    # First, process each individual organism.
    # (A) obtain elegiable splice site range
    # (B) scan for splice sites
    # (C) add the projected sites to the graph
    # (D) add splice sites to the stg collection graph
    for org in self.organism_set():
        # take the first (and only) orf of this organism
        theorf = self.get_orfs_of_graph(organism=org)[0]

        if forced_codingblock_ends.has_key(org):
            # the node that represents this site
            cbgEnd = forced_codingblock_ends[org]
            cbgEndNode = ( org,theorf.id,cbgEnd.pos )
            # add to the collection graph
            stg.add_node_and_object(cbgEndNode,cbgEnd)
            # ready with this organism, no splice site setting!
            #continue

        if next.__class__.__name__ == "LowSimilarityRegionCodingBlockGraph":
            # continue; all `donor` boundaries are hard-set 
            # no splice_site_range or actual site prediction needed
            continue

        ########################################################################
        ### get the considered splice site range
        ########################################################################

        # calculate considered splice site range based on EOF Orf object
        # take theorf.endPY + 2 (two) !, because EOF Orf is the start of the
        # STOP codon. Example:
        # ... tca TAG tac gtc ...
        # ... tca                   EOF Orf
        #         TAG               STOP codon
        #     ..a taG Tac gt.       perfect DONOR Site; PSSM-score ~7.7

        # calculate considered splice site range based on EOF Orf object
        (min_aa_pos, min_nt_pos) = self.minimal_eligable_donor_site_position(org)
        (max_aa_pos, max_nt_pos) = (theorf.endPY+2)/3, theorf.endPY+2

        if next and org in next.organism_set():
            (next_max_aa_pos, next_max_nt_pos) = self.maximal_eligable_donor_site_position(org,nextcbg=next)
        else:
            (next_max_aa_pos, next_max_nt_pos) = self.maximal_eligable_donor_site_position(org)
        if next_max_nt_pos < max_nt_pos:
             # minimal range falls within the orf's start point
             (max_aa_pos, max_nt_pos) = (next_max_aa_pos, next_max_nt_pos)

        if enlarge_5p_boundary_by:
            min_aa_pos = min_aa_pos - enlarge_5p_boundary_by
            min_nt_pos = min_nt_pos - (enlarge_5p_boundary_by*3)
        if enlarge_3p_boundary_by:
            max_aa_pos = max_aa_pos + enlarge_3p_boundary_by
            max_nt_pos = max_nt_pos + (enlarge_3p_boundary_by*3)


        # set range to stg Collection objects
        stg.set_consideredsplicesiterange(org,min_nt_pos,max_nt_pos)

        if forced_codingblock_ends.has_key(org):
            # ready with this organism, no splice site setting!
            continue

        ########################################################################
        ### obtain splice sites for current collection
        ########################################################################

        # scan for splice sites
        theorf.scan_orf_for_pssm_splice_sites(splicetype="donor",
                min_pssm_score=MIN_DONOR_PSSM_SCORE,allow_non_canonical=ALLOW_NON_CANONICAL_DONOR,
                non_canonical_min_pssm_score=NON_CANONICAL_MIN_DONOR_PSSM_SCORE,
                forced=True)

        # first, add the projected splicesites (they overrule true sites)
        if projected_donors.has_key(org):
            for projsite in projected_donors[org]:
                # check if we can ignore this site
                if not store_all_projected_sites:
                    if projsite.pos < min_nt_pos: continue
                    if max_nt_pos and projsite.pos > max_nt_pos: continue
                # create and add this projected site!
                projNode = ( org,theorf.id,projsite.pos )
                stg.add_node_and_object(projNode,projsite)

        # add the splice sites to the graph
        for dsq in theorf._donor_sites:

            if org == 'mgg' and theorf.id == 98: print dsq, dsq.pos, max_nt_pos

            # check if we can ignore this site
            if dsq.pos < min_nt_pos: continue
            if max_nt_pos and dsq.pos > max_nt_pos: continue

            # the node that represents this site
            dsqNode = ( org,theorf.id,dsq.pos )

            # check if this splice site is not already added as a projected site
            if dsqNode not in stg.get_nodes():
                stg.add_node_and_object(dsqNode,dsq)


    # now loop over all aligned combinations of organisms
    for ( (a,b,c,d),(g1,o1),(g2,o2) ), pacbporf in self.pacbps.iteritems():
        # only proces this combination if both organisms have splice sites!
        if g1 not in stg.organism_set(): continue
        if g2 not in stg.organism_set(): continue

        # now loop over all donor sites in Query and Sbjct
        # and align them in a graph; an edge is added if 2 sites
        # are less then ``ALIGNED_DONOR_MAX_TRIPLET_DISTANCE*3`` apart from each other
        for dsq in stg.get_organism_objects(g1):
            # the node that represents this site
            dsqNode  = ( g1,o1,dsq.pos )
            dsqClass = dsq.__class__.__name__

            for dss in stg.get_organism_objects(g2):
                # the node that represents this site
                dssNode  = ( g2,o2,dss.pos )
                dssClass = dss.__class__.__name__

                if 'CodingBlockEnd' in [ dsqClass,dssClass ]:
                    if dsqClass == dssClass:
                        # both CodingBlockEnd objects
                        dist = 0
                    else:
                        # calculate the distance in aligned nt positions
                        dist = pacbporf.get_distance_aligned_nucleotide_positions(
                                query = dsq.pos, sbjct = dss.pos
                                )

                    # check for the distance constrain
                    if dist > ALIGNED_DONOR_MAX_TRIPLET_DISTANCE*3: continue

                else:
                    # Both Donor sites; check for phase compatibility
                    if not allow_phase_shift and dsq.phase != dss.phase: continue

                    # calculate the distance in aligned nt positions
                    dist = pacbporf.get_distance_aligned_nucleotide_positions(
                            query = dsq.pos, sbjct = dss.pos
                            )

                    if dsq.phase == dss.phase:
                        # ignore uniformly aligned sites here
                        pass
                    elif allow_phase_shift and dist <= MAX_SPLICE_SITE_PHASE_SHIFT_NT_DISTANCE and\
                    dsq.phase != dss.phase and min([dsq.pssm_score, dss.pssm_score ]) >= MIN_DONOR_SITE_PHASE_SHIFT_PSSM_SCORE:
                        #print "PhaseShift:", dist, (g1,dsq.pos), (g2,dss.pos), min([dsq.pssm_score, dss.pssm_score ])
                        pass # a potential splice site phase shift
                    else:
                        continue
    
                    # check for the distance constrain for sites with uniform phase
                    if dist > ALIGNED_DONOR_MAX_TRIPLET_DISTANCE*3: continue

                # calculate binary entropies from Query
                if dsqClass == 'SpliceDonor':
                    dsqPositionPos, phaseQ = pacbporf.dnaposition_query(dsq.pos,forced_return=True)
                    entropyQ = pacbporf.alignment_entropy(dsqPositionPos,method='donor')
                elif dsqClass == 'ProjectedSpliceDonor':
                    entropyQ = dsq.entropy
                elif dsqClass == 'CodingBlockEnd':
                    entropyQ = 1.0
                else:
                    raise "NOT in [ SpliceDonor, ProjectedSpliceDonor, CodingBlockEnd ]"

                # calculate binary entropies from Sbjct
                if dssClass == 'SpliceDonor':
                    dssPositionPos, phaseS = pacbporf.dnaposition_query(dss.pos,forced_return=True)
                    entropyS = pacbporf.alignment_entropy(dssPositionPos,method='donor')
                elif dssClass == 'ProjectedSpliceDonor':
                    entropyS = dss.entropy
                elif dssClass == 'CodingBlockEnd':
                    entropyS = 1.0
                else:
                    raise "NOT in [ SpliceDonor, ProjectedSpliceDonor, CodingBlockEnd ]"

                # if here, then we have an aligned splice site!
                # calculate weight from distance, add edge and binary entropy values
                wt = 1.0 / ( 1.0 + float(dist/3) )
                stg.add_edge(dsqNode,dssNode,wt=wt)
                stg._edge_binary_entropies[(dsqNode,dssNode)] = (entropyQ,entropyS)
                stg._edge_binary_entropies[(dssNode,dsqNode)] = (entropyS,entropyQ)

    # return filled splicesitecollection graph
    return stg

# end of function harvest_elegiable_donor_sites


def harvest_elegiable_acceptor_sites(self,projected_acceptors={},forced_codingblock_ends={},prev=None,
    store_all_projected_sites=False,
    allow_phase_shift=False,
    enlarge_5p_boundary_by=None, # in AA coordinates
    enlarge_3p_boundary_by=None, # in AA coordinates
    ALIGNED_ACCEPTOR_MAX_TRIPLET_DISTANCE=None,
    MIN_ACCEPTOR_PSSM_SCORE=None,ALLOW_NON_CANONICAL_ACCEPTOR=None,
    NON_CANONICAL_MIN_ACCEPTOR_PSSM_SCORE=None ):
    """
    Harvest elegiable acceptor sites from this CodingBlockGraph into a AcceptorSiteCollectionGraph
    """
    if prev and prev.__class__.__name__ not in ["CodingBlockGraph","LowSimilarityRegionCodingBlockGraph"]:
        message = "prev must be a CodingBlock graph object, not a %s" % prev.__class__.__name__
        raise InproperlyAppliedArgument, message

    # update minimal pssm score to stg collection object
    stg = AcceptorSiteCollectionGraph()
    stg.MIN_PSSM_SCORE = MIN_ACCEPTOR_PSSM_SCORE
    stg.ALIGNED_SITE_AA_OFFSET = ALIGNED_ACCEPTOR_MAX_TRIPLET_DISTANCE

    # First, proces each individual organism.
    # (A) obtain elegiable splice site range
    # (B) scan for splice sites
    # (C) add the projected sites to the graph
    # (D) add splice sites to the stg collection graph
    for org in self.organism_set():
        # take the first (and only) orf of this organism
        theorf = self.get_orfs_of_graph(organism=org)[0]

        if forced_codingblock_ends.has_key(org):
            # the node that represents this site
            cbgSta = forced_codingblock_ends[org]
            cbgStaNode = ( org,theorf.id,cbgSta.pos )
            # add to the collection graph
            stg.add_node_and_object(cbgStaNode,cbgSta)
            # ready with this organism, no splice site setting!
            #continue

        if prev.__class__.__name__ == "LowSimilarityRegionCodingBlockGraph":
            # continue; all `acceptor` boundaries are hard-set
            # no splice_site_range or actual site prediction needed
            continue

        ########################################################################
        ### get the considered splice site range
        ########################################################################

        (max_aa_pos, max_nt_pos) = self.maximal_eligable_acceptor_site_position(org)
        (min_aa_pos, min_nt_pos) = (theorf.startPY-2)/3, theorf.startPY-2

        if prev and org in prev.organism_set():
            (next_min_aa_pos, next_min_nt_pos) = self.minimal_eligable_acceptor_site_position(org,prevcbg=prev)
        else:
            (next_min_aa_pos, next_min_nt_pos) = self.minimal_eligable_acceptor_site_position(org)
        if next_min_nt_pos > min_nt_pos:
             # minimal range falls within the orf's start point
             (min_aa_pos, min_nt_pos) = (next_min_aa_pos, next_min_nt_pos)
        if enlarge_5p_boundary_by:
            min_aa_pos = min_aa_pos - enlarge_5p_boundary_by
            min_nt_pos = min_nt_pos - (enlarge_5p_boundary_by*3)
        if enlarge_3p_boundary_by:
            max_aa_pos = max_aa_pos + enlarge_3p_boundary_by
            max_nt_pos = max_nt_pos + (enlarge_3p_boundary_by*3)

        # set range to stg Collection objects
        stg.set_consideredsplicesiterange(org,min_nt_pos,max_nt_pos)

        if forced_codingblock_ends.has_key(org):
            # ready with this organism, no splice site setting!
            continue

        ########################################################################
        ### obtain splice sites for current collection
        ########################################################################

        # scan for splice sites
        theorf.scan_orf_for_pssm_splice_sites(splicetype="acceptor",
                min_pssm_score=MIN_ACCEPTOR_PSSM_SCORE,allow_non_canonical=ALLOW_NON_CANONICAL_ACCEPTOR,
                non_canonical_min_pssm_score=NON_CANONICAL_MIN_ACCEPTOR_PSSM_SCORE)

        # first, add the projected splicesites (they overrule true sites)
        if projected_acceptors.has_key(org):
            for projsite in projected_acceptors[org]:
                # check if we can ignore this site
                if not store_all_projected_sites:
                    if projsite.pos < min_nt_pos: continue
                    if max_nt_pos and projsite.pos > max_nt_pos: continue
                # create and add this projected site!
                projNode = ( org,theorf.id,projsite.pos )
                stg.add_node_and_object(projNode,projsite)

        # add the splice sites to the graph
        for asq in theorf._acceptor_sites:
            # check if we can ignore this site
            if asq.pos < min_nt_pos: continue
            if max_nt_pos and asq.pos > max_nt_pos: continue

            # the node that represents this site
            asqNode = ( org,theorf.id,asq.pos )

            # check if this splice site is not already added as a projected site
            if asqNode not in stg.get_nodes():
                stg.add_node_and_object(asqNode,asq)


    # now loop over all aligned combinations of organisms
    for ( (a,b,c,d),(g1,o1),(g2,o2) ), pacbporf in self.pacbps.iteritems():
        # only proces this combination if both organisms have splice sites!
        if g1 not in stg.organism_set(): continue
        if g2 not in stg.organism_set(): continue

        # now loop over all acceptor sites in Query and Sbjct
        # and align them in a graph; an edge is added if 2 sites
        # are less then ``ALIGNED_ACCEPTOR_MAX_TRIPLET_DISTANCE*3`` apart from each other
        for asq in stg.get_organism_objects(g1):
            # the node that represents this site
            asqNode = ( g1,o1,asq.pos )
            asqClass = asq.__class__.__name__

            for ass in stg.get_organism_objects(g2):
                # the node that represents this site
                assNode = ( g2,o2,ass.pos )
                assClass = ass.__class__.__name__

                if 'CodingBlockStart' in [ asqClass,assClass ]:
                    if asqClass == assClass:
                        # both CodingBlockEnd objects
                        dist = 0
                    else:
                        # calculate the distance in aligned nt positions
                        dist = pacbporf.get_distance_aligned_nucleotide_positions(
                                query = asq.pos, sbjct = ass.pos
                                )

                        # check for the distance constrain
                        if dist > ALIGNED_ACCEPTOR_MAX_TRIPLET_DISTANCE*3: continue

                else:
                    # Both Acceptor sites; check for phase compatibility
                    if not allow_phase_shift and asq.phase != ass.phase: continue

                    # calculate the distance in aligned nt positions
                    dist = pacbporf.get_distance_aligned_nucleotide_positions(
                            query = asq.pos, sbjct = ass.pos
                            )

                    if asq.phase == ass.phase:
                        # ignore uniformly aligned sites here
                        pass
                    elif allow_phase_shift and dist <= MAX_SPLICE_SITE_PHASE_SHIFT_NT_DISTANCE and\
                    asq.phase != ass.phase and min([ asq.pssm_score, ass.pssm_score ]) >= MIN_ACCEP_SITE_PHASE_SHIFT_PSSM_SCORE:
                        #print "PhaseShift:", dist, (g1,asq.pos), (g2,ass.pos), min([ asq.pssm_score, ass.pssm_score ]) 
                        pass # a potential splice site phase shift
                    else:
                        continue

                    # check for the distance constrain for sites of uniform phase
                    if dist > ALIGNED_ACCEPTOR_MAX_TRIPLET_DISTANCE*3: continue

                # calculate binary entropies from Query
                if asqClass == 'SpliceAcceptor':
                    asqPositionPos, phaseQ = pacbporf.dnaposition_query(asq.pos,forced_return=True)
                    entropyQ = pacbporf.alignment_entropy(asqPositionPos,method='acceptor')
                elif asqClass == 'ProjectedSpliceAcceptor':
                    entropyQ = asq.entropy
                elif asqClass == 'CodingBlockStart':
                    entropyQ = 1.0
                else:
                    raise "NOT a SpliceAcceptor or a ProjectedSpliceAcceptor"

                # calculate binary entropies from Sbjct
                if assClass == 'SpliceAcceptor':
                    assPositionPos, phaseS = pacbporf.dnaposition_query(ass.pos,forced_return=True)
                    entropyS = pacbporf.alignment_entropy(assPositionPos,method='acceptor')
                elif assClass == 'ProjectedSpliceAcceptor':
                    entropyS = ass.entropy
                elif assClass == 'CodingBlockStart':
                    entropyS = 1.0
                else:
                    raise "NOT a SpliceAcceptor or a ProjectedSpliceAcceptor"

                # if here, then we have an aligned splice site!
                # calculate weight from distance, add edge and binary entropy values
                wt = 1.0 / ( 1.0 + float(dist/3) )
                stg.add_edge(asqNode,assNode,wt=wt)
                stg._edge_binary_entropies[(asqNode,assNode)] = (entropyQ,entropyS)
                stg._edge_binary_entropies[(assNode,asqNode)] = (entropyS,entropyQ)

    # return filled splicesitecollection graph
    return stg

# end of function harvest_elegiable_acceptor_sites


def harvest_elegiable_tss_sites(self,max_aa_distance=ALIGNED_TSS_MAX_AA_DISTANCE,
    tss_min_pssm_score=TSS_MIN_PSSM_SCORE,
    skip_nonelegiable_sites=True):
    """
    """
    # update minimal pssm score to stg collection object
    stg = TranslationalStartSiteCollectionGraph()
    stg.MIN_PSSM_SCORE = tss_min_pssm_score 
    stg.ALIGNED_SITE_AA_OFFSET = max_aa_distance


    # First, proces each individual organism.
    for org in self.organism_set():
        # take the first (and only) orf of this organism
        theorf = self.get_orfs_of_graph(organism=org)[0]
        # ready if there are no potential tss loci (no ATG sequence)
        if not theorf.has_start(): continue
        # scan for tss loci
        theorf.scan_orf_for_pssm_tss(min_pssm_score=tss_min_pssm_score)


        if skip_nonelegiable_sites:
            # get the considered TSS range
            (min_aa_pos, min_nt_pos) = self.minimal_eligable_tss_position(org)
            (max_aa_pos, max_nt_pos) = self.maximal_eligable_tss_position(org)
        else:
            (min_aa_pos, min_nt_pos) = None, None
            (max_aa_pos, max_nt_pos) = None, None

        for tss in theorf._tss_sites:
            # check if we can ignore this site
            if min_nt_pos and tss.pos < min_nt_pos: continue
            if max_nt_pos and tss.pos > max_nt_pos: continue
            # an accepted site; add to TSS Collection Graph
            startpos = tss.pos / 3
            tssNode = ( org, theorf.id, startpos, tss.pos )
            stg.add_node_and_object(tssNode,tss)


    # Second, evaluate all cross combinations
    for ( (a,b,c,d),(g1,o1),(g2,o2) ), pacbporf in self.pacbps.iteritems():
        # only proces this combination if both organisms have splice sites!
        if g1 not in stg.organism_set(): continue
        if g2 not in stg.organism_set(): continue

        # now loop over all TSS in Query and Sbjct
        # and align them in a graph; an edge is added if 2 sites
        # are less then ``max_aa_distance`` apart from each other
        for tssQ in stg.get_organism_objects(g1):
            # the node that represents this site
            startQpos = tssQ.pos / 3
            startQnode  = ( g1, o1, startQpos, tssQ.pos )
            for tssS in stg.get_organism_objects(g2):
                # the node that represents this site
                startSpos = tssS.pos / 3
                startSnode = ( g2, o2, startSpos, tssS.pos )

                # get distance between (aligned) start-codons
                dist = pacbporf.get_distance_aligned_protein_positions(
                        query=startQpos,sbjct=startSpos)

                # continue if distance between start sites is to big
                if dist > max_aa_distance: continue

                # calculate binary entropies from both positions
                startQpositionPos,phaseQ = pacbporf.dnaposition_query(tssQ.pos,forced_return=True)
                startSpositionPos,phaseS = pacbporf.dnaposition_sbjct(tssS.pos,forced_return=True)
                entropyQ = pacbporf.alignment_entropy(startQpositionPos,method='left')
                entropyS = pacbporf.alignment_entropy(startSpositionPos,method='left')

                # calculate a weight from distance between startQpos and startSpos
                wt = 1.0 / ( 1.0 + float(dist) )

                # check if edge already in graph
                if stg.has_edge( startQnode, startSnode ):
                    _wt = stg.weights[( startQnode, startSnode )]
                    if wt > _wt:
                        stg.set_edge_weight( startQnode, startSnode, wt=wt )
                        # and add binary entropy values
                        stg._edge_binary_entropies[(startQnode, startSnode)] = (entropyQ,entropyS)
                        stg._edge_binary_entropies[(startSnode, startQnode)] = (entropyS,entropyQ)
                else:
                    stg.add_edge( startQnode, startSnode, wt=wt )
                    # and add binary entropy values
                    stg._edge_binary_entropies[(startQnode, startSnode)] = (entropyQ,entropyS)
                    stg._edge_binary_entropies[(startSnode, startQnode)] = (entropyS,entropyQ)


    # Get tcode data for these start codon nodes
    # Assuming that this is indeed the start-codon,
    # the stretch of ATG untill max(OMSR) will be coding.
    # Take the length of this stretch (in nt) as right/3p/upstream window size
    omsr = self.overall_minimal_spanning_range()
    for (org,orfid,aaPos,dnaPos) in stg.get_nodes():
        theorf = self.get_orfs_of_graph(organism=org)[0]
        right_window_size = ( max(omsr[(org,orfid)])+1 - aaPos )*3
        # confirm that window size is not < 0; this is possible
        # once the Methionine/TSS is located downstream of the
        # OMSR max site
        if right_window_size <= 0:
            right_window_size = stg._TCODE_3P_WINDOWSIZE
        # calculate the average TCODE scores for the windows
        ( tcode5p,tcode3p ) = theorf.tcode_entropy_of_pos(
                aaPos,
                window_left=stg._TCODE_5P_WINDOWSIZE,
                window_right=right_window_size,
                )
        stg._tcode5pscore[(org,orfid,aaPos,dnaPos)] = tcode5p
        stg._tcode3pscore[(org,orfid,aaPos,dnaPos)] = tcode3p

    # return filled tsscollection graph
    return stg

# end of function harvest_elegiable_tss_sites


#def align_start_codons_NO_TSS(self,max_aa_distance=ELIGABLE_ALIGNED_START_SITES_AA_OFFSET):
#    """
#    """
#    stg = AlignedStartCodonGraph()
#    for ( (a,b,c,d),(g1,o1),(g2,o2) ), pacbporf in self.pacbps.iteritems():
#        if not pacbporf.orfQ.has_start(): continue
#        if not pacbporf.orfS.has_start(): continue
#
#        for startQpos in pacbporf.orfQ.potential_start_aa_positions():
#            startQdnapos = pacbporf.orfQ.aapos2dnapos(startQpos)
#            startQnode   = (g1,o1,startQpos,startQdnapos)
#            for startSpos in pacbporf.orfS.potential_start_aa_positions():
#                startSdnapos = pacbporf.orfS.aapos2dnapos(startSpos)
#                startSnode   = (g2,o2,startSpos,startSdnapos)
#                # get distance between (aligned) start-codons
#                dist = pacbporf.get_distance_aligned_protein_positions(
#                    query=startQpos,sbjct=startSpos)
#
#                # continue if distance between start sites is to big
#                if dist > max_aa_distance: continue
#
#                # calculate binary entropies from both positions
#                startQpositionPos = pacbporf.alignmentposition_by_query_pos(startQpos,forced_return=True)
#                startSpositionPos = pacbporf.alignmentposition_by_sbjct_pos(startSpos,forced_return=True)
#                entropyQ = pacbporf.alignment_entropy(startQpositionPos,method='left')
#                entropyS = pacbporf.alignment_entropy(startSpositionPos,method='left')
#
#                # add these nodes if not in the graph yet
#                if startQnode not in stg.get_nodes(): stg.add_node(startQnode)
#                if startSnode not in stg.get_nodes(): stg.add_node(startSnode)
#
#                # calculate a weight from algQpos and algSpos
#                wt = 1.0 / ( 1.0 + float(dist) )
#
#                # check if edge already in graph
#                if stg.has_edge( startQnode, startSnode ):
#                    _wt = stg.weights[( startQnode, startSnode )]
#                    if wt > _wt:
#                        stg.set_edge_weight( startQnode, startSnode, wt=wt )
#                        # and add binary entropy values
#                        stg._edge_binary_entropies[(startQnode, startSnode)] = (entropyQ,entropyS)
#                        stg._edge_binary_entropies[(startSnode, startQnode)] = (entropyS,entropyQ)
#                else:
#                    stg.add_edge( startQnode, startSnode, wt=wt )
#                    # and add binary entropy values
#                    stg._edge_binary_entropies[(startQnode, startSnode)] = (entropyQ,entropyS)
#                    stg._edge_binary_entropies[(startSnode, startQnode)] = (entropyS,entropyQ)
#
#    # done; set graph object to object
#    self._startcodongraph = stg
#
## end of function align_start_codons_NO_TSS


def align_stop_codons(codingblockgraph,alignedstopcodongraph):
    """
    Align the stop-codons from the aligned orfs
    in the pacbporfs in a graph object.
    Perfectly alignable stop-codons result in a
    graph with total_weight() == 1.0, all below 1.0
    means offsets in the alignment. Non-perfect aligned
    stop-codons are common, but (very) poor alignable
    stop-codons are very rare.
    """
    for ( (a,b,c,d),(g1,o1),(g2,o2) ), pacbporf in codingblockgraph.pacbps.iteritems():
        stopQpos = pacbporf.orfQ.protein_endPY+1
        stopSpos = pacbporf.orfS.protein_endPY+1
        stopQdnapos = pacbporf.orfQ.proteinpos2dnapos(stopQpos)
        stopSdnapos = pacbporf.orfS.proteinpos2dnapos(stopSpos)
        stopQnode = (g1,o1,stopQpos,stopQdnapos)
        stopSnode = (g2,o2,stopSpos,stopSdnapos)
        # get distance between (aligned) stop-codons
        dist = pacbporf.get_distance_aligned_protein_positions(
                query=stopQpos,sbjct=stopSpos)
        # calculate weight from distance
        wt = 1.0 / ( 1.0 + float(dist) )

        # calculate binary entropies from both positions
        stopQpositionPos,phaseQ = pacbporf.dnaposition_query(stopQdnapos,forced_return=True)
        stopSpositionPos,phaseS = pacbporf.dnaposition_sbjct(stopSdnapos,forced_return=True)
        entropyQ = pacbporf.alignment_entropy(stopQpositionPos,method='right')
        entropyS = pacbporf.alignment_entropy(stopSpositionPos,method='right')

        # check if node in graph
        if stopQnode not in alignedstopcodongraph.get_nodes(): alignedstopcodongraph.add_node(stopQnode)
        if stopSnode not in alignedstopcodongraph.get_nodes(): alignedstopcodongraph.add_node(stopSnode)

        # check if edge already in graph
        if alignedstopcodongraph.has_edge( stopQnode, stopSnode ):
            _wt = alignedstopcodongraph.weights[( stopQnode, stopSnode )]
            if wt > _wt:
                alignedstopcodongraph.set_edge_weight( stopQnode, stopSnode, wt=wt )
                # and add binary entropy values
                alignedstopcodongraph._edge_binary_entropies[(stopQnode, stopSnode)] = (entropyQ,entropyS)
                alignedstopcodongraph._edge_binary_entropies[(stopSnode, stopQnode)] = (entropyS,entropyQ)
        else:
            alignedstopcodongraph.add_edge( stopQnode, stopSnode, wt=wt )
            # and add binary entropy values
            alignedstopcodongraph._edge_binary_entropies[(stopQnode, stopSnode)] = (entropyQ,entropyS)
            alignedstopcodongraph._edge_binary_entropies[(stopSnode, stopQnode)] = (entropyS,entropyQ)


    # Get tcode data for these stop codon nodes
    # Assuming that this is indeed the stop-codon,
    # the stretch of min(OMSR) untill TGA will be coding.
    # Take the length of this stretch (in nt) as left/5p/downstream window size
    omsr = codingblockgraph.overall_minimal_spanning_range()
    for (org,orfid,aaPos,dnaPos) in alignedstopcodongraph.get_nodes():
        theorf = codingblockgraph.get_orfs_of_graph(organism=org)[0]
        left_window_size = ( aaPos - min(omsr[(org,orfid)]) )*3
        ( tcode5p,tcode3p ) = theorf.tcode_entropy_of_pos(
                aaPos,
                window_left=left_window_size,
                window_right=alignedstopcodongraph._TCODE_3P_WINDOWSIZE,
                )
        alignedstopcodongraph._tcode5pscore[(org,orfid,aaPos,dnaPos)] = tcode5p
        alignedstopcodongraph._tcode3pscore[(org,orfid,aaPos,dnaPos)] = tcode3p

    # return filled stopcodon graph
    return alignedstopcodongraph

# end of function align_stop_codons


def minimal_eligable_donor_site_position(cbg,organism):
    """
    Get the position on this CodingBlockGraph from where to take donor positions into account

    @type  organism: *
    @param organism: organism identifier (or None)

    @rtype:  tuple
    @return: ( absolute_aa_position, absolute_nt_position )

    @atttention: requires global variable ELIGABLE_DONOR_SITE_LEFT_OF_OMSR_AA_OFFSET
    """
    # take the first (and only) orf of this organism
    orf_of_org = cbg.get_orfs_of_graph(organism=organism)[0]
    omsr = cbg.overall_minimal_spanning_range(organism=organism)
    # adjust ELIGABLE_DONOR_SITE_LEFT_OF_OMSR_AA_OFFSET based on identity of the cbg
    offset = int( ELIGABLE_DONOR_SITE_LEFT_OF_OMSR_AA_OFFSET * cbg.get_genetree().identity() )
    # calculate absolute aa and nt positions from where to take donors into account
    abs_aa_pos = max([ min(omsr), max(omsr)-offset ])
    abs_nt_pos = orf_of_org.proteinpos2dnapos(abs_aa_pos)
    return ( abs_aa_pos, abs_nt_pos )

# end of function minimal_eligable_donor_site_position


def maximal_eligable_acceptor_site_position(cbg,organism):
    """
    Get the position on this CodingBlockGraph up to where to take acceptor positions into account

    @type  cbg: CodingBlockGraph
    @param cbg: CodingBlockGraph instance

    @type  organism: *
    @param organism: organism identifier (or None)

    @rtype:  tuple
    @return: ( absolute_aa_position, absolute_nt_position )

    @attention: requires global variable ELIGABLE_ACCEPTOR_SITE_RIGTH_OF_OMSR_AA_OFFSET
    """
    # take the first (and only) orf of this organism
    orf_of_org = cbg.get_orfs_of_graph(organism=organism)[0]
    omsr = cbg.overall_minimal_spanning_range(organism=organism)
    # adjust ELIGABLE_ACCEPTOR_SITE_RIGTH_OF_OMSR_AA_OFFSET based on identity of the cbg
    offset = int( ELIGABLE_ACCEPTOR_SITE_RIGTH_OF_OMSR_AA_OFFSET * cbg.get_genetree().identity() )
    # calculate absolute aa and nt positions untill where to take acceptors into account
    abs_aa_pos = min([ max(omsr), min(omsr)+offset ])
    abs_nt_pos = orf_of_org.proteinpos2dnapos(abs_aa_pos)
    return ( abs_aa_pos, abs_nt_pos )

# end of function maximal_eligable_acceptor_site_position


def maximal_eligable_donor_site_position(cbg,organism,nextcbg=None):
    """
    Get the position on this CodingBlockGraph untill where to take donor positions into account

    @type  cbg: CodingBlockGraph
    @param cbg: CodingBlockGraph instance

    @type  organism: *
    @param organism: organism identifier (or None)

    @type  nextcbg: None or CodingBlockGraph
    @param nextcbg: the CodingBlockGraph next in line in an ordered CBG series

    @rtype:  tuple
    @return: ( absolute_aa_position, absolute_nt_position )

    @atttention: requires global variable ELIGABLE_DONOR_SITE_RIGTH_OF_OMSR_AA_OFFSET and ELIGABLE_DONOR_SITE_MINIMAL_AA_OFFSET
    """
    # take the first (and only) orf of this organism
    orf_of_org = cbg.get_orfs_of_graph(organism=organism)[0]
    omsr = cbg.overall_minimal_spanning_range(organism=organism)
    # adjust ELIGABLE_DONOR_SITE_RIGTH_OF_OMSR_AA_OFFSET based on identity of the cbg
    # when, due to very high identity, lt ELIGABLE_DONOR_SITE_MINIMAL_AA_OFFSET, this bottom value is taken as offset
    offset = max([ ELIGABLE_DONOR_SITE_MINIMAL_AA_OFFSET, int( ELIGABLE_DONOR_SITE_RIGTH_OF_OMSR_AA_OFFSET * ( 1.0 - cbg.get_genetree().identity() ) ) ])
    # calculate absolute aa and nt positions untill where to take donors into account
    abs_aa_pos = max(omsr)+offset
    # if the next cbg is provided, check if the splice site range does not interfere with its omsr

    if nextcbg:
        nextminomsr = min(nextcbg.overall_minimal_spanning_range(organism=organism))
        if nextminomsr < max(omsr):
            # (slightly) overlapping CBGs -> take spatious overlap to allow the projected sites!
            pass
        else: 
            # correct nextminomsr with 6 AA -> the distance allowed in intron projection is 5 AA
            abs_aa_pos = min([ abs_aa_pos, nextminomsr+6 ])

    # get nt pos of aa pos and return
    abs_nt_pos = orf_of_org.proteinpos2dnapos(abs_aa_pos)
    return ( abs_aa_pos, abs_nt_pos )

# end of function maximal_eligable_donor_site_position


def minimal_eligable_acceptor_site_position(cbg,organism,prevcbg=None):
    """
    Get the position on this CodingBlockGraph from where to take acceptor positions into account

    @type  cbg: CodingBlockGraph
    @param cbg: CodingBlockGraph instance

    @type  organism: *
    @param organism: organism identifier (or None)

    @type  prevcbg: None or CodingBlockGraph
    @param prevcbg: the CodingBlockGraph previous in line in an ordered CBG series

    @rtype:  tuple
    @return: ( absolute_aa_position, absolute_nt_position )

    @atttention: requires global variable ELIGABLE_ACCEPTOR_SITE_LEFT_OF_OMSR_AA_OFFSET and ELIGABLE_ACCEPTOR_SITE_MINIMAL_AA_OFFSET
    """
    # take the first (and only) orf of this organism
    orf_of_org = cbg.get_orfs_of_graph(organism=organism)[0]
    omsr = cbg.overall_minimal_spanning_range(organism=organism)

    # adjust ELIGABLE_ACCEPTOR_SITE_LEFT_OF_OMSR_AA_OFFSET based on identity of the cbg
    # when, due to very high identity, lt ELIGABLE_DONOR_SITE_MINIMAL_AA_OFFSET, this bottom value is taken as offset
    offset = max([ ELIGABLE_ACCEPTOR_SITE_MINIMAL_AA_OFFSET, int( ELIGABLE_ACCEPTOR_SITE_LEFT_OF_OMSR_AA_OFFSET * ( 1.0 - cbg.get_genetree().identity() ) ) ])
    # calculate absolute aa and nt positions from where to take acceptors into account
    abs_aa_pos = min(omsr)-offset
    # if the prev cbg is provided, check if the splice site range does not interfere with its omsr
    if prevcbg:
        prevmaxomsr = max(prevcbg.overall_minimal_spanning_range(organism=organism))
        if prevmaxomsr > min(omsr):
            # (slightly) overlapping CBGs -> take spatious overlap to allow the projected sites!
            pass
        else:
            # correct nextminomsr with 6 AA -> the distance allowed in intron projection is 5 AA
            abs_aa_pos = max([ abs_aa_pos, prevmaxomsr-6 ])

    # get nt pos of aa pos and return
    abs_nt_pos = orf_of_org.proteinpos2dnapos(abs_aa_pos)
    return ( abs_aa_pos, abs_nt_pos )

# end of function minimal_eligable_acceptor_site_position


def minimal_eligable_tss_position(cbg,organism):
    """
    Get the position on the Orf of this organism in the CBG from where to take TSS positions into account

    @type  cbg: CodingBlockGraph
    @param cbg: CodingBlockGraph instance

    @type  organism: *
    @param organism: organism identifier (or None)

    @rtype:  tuple
    @return: ( absolute_aa_position, absolute_nt_position )
    """
    # take the first (and only) orf of this organism
    orf_of_org = cbg.get_orfs_of_graph(organism=organism)[0]
    omsr = cbg.overall_minimal_spanning_range(organism=organism)
    # calculate absolute aa and nt positions from where to take acceptors into account
    if ELIGABLE_ALIGNED_TSS_5P_AA_OFFSET == None:
        abs_aa_pos = orf_of_org.protein_startPY
    else:
        abs_aa_pos = max([ min(omsr)-ELIGABLE_ALIGNED_TSS_5P_AA_OFFSET, orf_of_org.protein_startPY ])
    abs_nt_pos = orf_of_org.proteinpos2dnapos(abs_aa_pos)
    return ( abs_aa_pos, abs_nt_pos )

# end of function minimal_eligable_tss_position


def maximal_eligable_tss_position(cbg,organism):
    """
    Get the position on the Orf of this organism in the CBG untill where to take TSS positions into account

    @type  cbg: CodingBlockGraph
    @param cbg: CodingBlockGraph instance

    @type  organism: *
    @param organism: organism identifier (or None)

    @rtype:  tuple
    @return: ( absolute_aa_position, absolute_nt_position )
    """
    # take the first (and only) orf of this organism
    orf_of_org = cbg.get_orfs_of_graph(organism=organism)[0]
    omsr = cbg.overall_minimal_spanning_range(organism=organism)
    # calculate absolute aa and nt positions from where to take acceptors into account
    if ELIGABLE_ALIGNED_TSS_3P_AA_OFFSET == None:
        abs_aa_pos = orf_of_org.protein_endPY
    else:
        abs_aa_pos = min([ min(omsr)+ELIGABLE_ALIGNED_TSS_3P_AA_OFFSET, orf_of_org.protein_endPY ])
    abs_nt_pos = orf_of_org.proteinpos2dnapos(abs_aa_pos)
    return ( abs_aa_pos, abs_nt_pos )

# end of function maximal_eligable_tss_position


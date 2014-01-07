"""
Functions for optimization of CodingBlockGraphs (CBGs)
"""

# Abgp Imports
from lib_clustalw import clustalw
import pacb

# graphAbgp Imports
from exceptions import NoOverallMinimalSpanningRange

# CodingBlockGraph Imports
from operations import _update_cbg_with_pacbporf_replacements

# Python Imports
from sets import Set
from copy import deepcopy

# Global variables
from settings.codingblockgraph import (
    CBG_OPTIMIZE_MAXIMAL_IDENTITY,
    CBG_OPTIMIZE_CLUSTALW_GAP_SIZE,
    CBG_OPTIMIZE_MINIMAL_BITSCORE_RATIO,
    CBG_OPTIMIZE_MINIMAL_IDENTITY_RATIO,
    CBG_PACBPGAP_NEARBY_OMSR_OFFSET,
    CBG_PACBPGAP_NEARBY_OMSR_GAP_SIZE,
    )


def correct_pacbpgaps_nearby_omsr(cbg,
    omsr_offset=CBG_PACBPGAP_NEARBY_OMSR_OFFSET,
    gap_size=CBG_PACBPGAP_NEARBY_OMSR_GAP_SIZE,
    omit5pside=False, omit3pside=False, verbose=False):
    """
    Correct a CBG when (large) gaps are occuring on the edges of the OMSR

    @type  cbg: CodingBlockGraph
    @param cbg: CodingBlockGraph instance (with cbgIF objects!)

    @type  omsr_offset: positive (small) integer
    @param omsr_offset: offset around the OMSR to check for gap occurrence 

    @type  omit5pside: Boolean
    @param omit5pside: Do not process the 5' side (left) of the CBG

    @type  omit3pside: Boolean
    @param omit3pside: Do not process the 3' side (rigth) of the CBG

    @type  verbose: Boolean
    @param verbose: print status/debugging messages to STDOUT

    @attention: When the CBG's OMSR is lost in the proces ->
                NoOverallMinimalSpanningRange exception is raised

    @rtype:  Boolean
    @return: Status for if the CBG is optimized/changed or not 
    """
    # gather data of the current cbg to compare before/after optimization
    current_cbg_total_weight = cbg.total_weight()
    current_cbg_string       = str(cbg)
    current_cbg_omsr         = cbg.overall_minimal_spanning_range()
    current_cbg_maxsr        = cbg.maximal_spanning_range()

    # first, check 3' side of OMSR
    PACBPS_CORRECTED = 0

    if not omit3pside:
        replacements = {}
        for (currentkey,node1,node2),pacbporf in cbg.pacbps.iteritems():
            # get slice of the pacbporf around the max(OMSR) query value
            aaQpos = max(current_cbg_omsr[node1])
            (q,m,s,coords) = pacbporf.alignmentpart_by_query(
                        aaQpos-omsr_offset, aaQpos+1+omsr_offset )
            # If more gaps in this alignment slice then expected ->
            # a pacbp split will follow. The slice is of size (2*omsr_offset)+1
            if q.count('-') >= gap_size:
                # convert (back) to pacbp and obtain position where to split
                pacbp = pacb.conversion.pacbporf2pacbp(pacbporf)

                # obtain splitpos based on gap in query
                splitpos = _obtain_3p_splitpos( q, coords[0], coords[1],
                                    pacbp.query, pacbp.query_start )

            elif s.count('-') >= gap_size:
                # convert (back) to pacbp and obtain position where to split
                pacbp = pacb.conversion.pacbporf2pacbp(pacbporf)

                # obtain splitpos based on gap in sbjct
                splitpos = _obtain_3p_splitpos( s, coords[2], coords[3],
                                    pacbp.sbjct, pacbp.sbjct_start )

            else:
                # no gap to split on; go to next pacbp
                continue

            ####################################################################
            if verbose:
                print "3p,", node1,node2, aaQpos, (q,m,s,coords)
                print pacbp, "relative splitpos:", splitpos
                pacbp.print_protein(_linesize=125)
            ####################################################################

            # when here, split pacbp on splitpos and recreate PacbPORF
            pacbpL = pacb.splitting.split_pacb_on_coordinates(
                                pacbp,(splitpos,splitpos),returnside='left' )
            if pacbpL:
                # recreate pacbpORF from pacbp
                newpacbporf = pacb.conversion.pacbp2pacbporf(
                                pacbpL, pacbporf.orfQ, pacbporf.orfS )
                newpacbporf.extend_pacbporf_after_stops()
                # store to replacements dict
                replacements[(currentkey,node1,node2)] = newpacbporf
                ############################################################
                if verbose:
                    print pacbpL
                    pacbpL.print_protein(_linesize=125)
                    print newpacbporf
                ############################################################
                # increase counter for how much pacbps are corrected
                PACBPS_CORRECTED+=1
            else:
                # Split failed. Wrong coords and/or no alignment left
                pass

        ########################################################################
        if verbose and replacements:
            print "PACBP REPLACEMENTS RIGTH/3p TO BE DONE:", len(replacements)
            print cbg
            cbg.printmultiplealignment()
            omsr = cbg.overall_minimal_spanning_range()
            for node in omsr.keys(): print node,min(omsr[node]),max(omsr[node])
        ########################################################################
    
        # do the replacements of 3' PacbP corrections
        status = _update_cbg_with_pacbporf_replacements(cbg,replacements)
        if status == True:
            pass    # cbg succesfully updated; still an OMSR
        elif status == False:
            # raise a NoOverallMinimalSpanningRange Exception
            raise NoOverallMinimalSpanningRange, str(cbg)
        else:
            pass 

        ########################################################################
        if verbose and replacements and status:
            print "PACBP REPLACEMENTS RIGTH/3p DONE!!:", len(replacements)
            print cbg
            cbg.printmultiplealignment()
            omsr = cbg.overall_minimal_spanning_range()
            for node in omsr.keys(): print node,min(omsr[node]),max(omsr[node])
        ########################################################################

    else:
        # omit the 3' side of this CBG
        pass


    if not omit5pside:
        # now, repeat this function for the 5' side of the cbg
        replacements = {}
        for (currentkey,node1,node2),pacbporf in cbg.pacbps.iteritems():
            # get slice of the pacbporf around the min(OMSR) query value
            aaQpos = min(current_cbg_omsr[node1])
            (q,m,s,coords) = pacbporf.alignmentpart_by_query(
                    aaQpos-omsr_offset, aaQpos+1+omsr_offset )
            # If more gaps in this alignment slice then expected ->
            # a pacbp split will follow. The slice is of size (2*omsr_offset)+1

            if q.count('-') >= gap_size:
                # convert (back) to pacbp and obtain position where to split
                pacbp = pacb.conversion.pacbporf2pacbp(pacbporf)

                # obtain splitpos based on gap in query
                splitpos = _obtain_5p_splitpos( q, coords[0], coords[1],
                                    pacbp.query, pacbp.query_start )
    
            elif s.count('-') >= gap_size:
                # convert (back) to pacbp and obtain position where to split
                pacbp = pacb.conversion.pacbporf2pacbp(pacbporf)

                # obtain splitpos based on gap in query
                splitpos = _obtain_5p_splitpos( s, coords[2], coords[3],
                                    pacbp.sbjct, pacbp.sbjct_start )
    
            else:
                # no gap to split on; go to next pacbp
                continue
            
            ####################################################################
            if verbose:
                print "5p,", node1,node2, aaQpos, (q,m,s,coords)
                print pacbp, "relative splitpos:", splitpos
                pacbp.print_protein(_linesize=125)
            ####################################################################

            # when here, split pacbp on splitpos and recreate PacbPORF
            pacbpR = pacb.splitting.split_pacb_on_coordinates(
                        pacbp,(splitpos,splitpos),returnside='rigth' )
            if pacbpR:
                # recreate pacbpORF from pacbp
                newpacbporf = pacb.conversion.pacbp2pacbporf(
                            pacbpR, pacbporf.orfQ, pacbporf.orfS )
                newpacbporf.extend_pacbporf_after_stops()
                # store to replacements dict
                replacements[(currentkey,node1,node2)] = newpacbporf
                ############################################################
                if verbose:
                    print pacbpR
                    pacbpR.print_protein(_linesize=125)
                    print newpacbporf
                ############################################################
                # increase counter for how much pacbps are corrected
                PACBPS_CORRECTED+=1

        ########################################################################
        if verbose and replacements:
            print "PACBP REPLACEMENTS LEFT/5p TO BE DONE:", len(replacements)
            print cbg
            cbg.printmultiplealignment()
            omsr = cbg.overall_minimal_spanning_range()
            for node in omsr.keys(): print node,min(omsr[node]),max(omsr[node])
        ########################################################################

        # do the replacements of 5' PacbP corrections
        status = _update_cbg_with_pacbporf_replacements(cbg,replacements)
        if status == True:
            pass    # cbg succesfully updated; still an OMSR
        elif status == False:
            # raise a NoOverallMinimalSpanningRange Exception
            raise NoOverallMinimalSpanningRange, str(cbg)
        else:
            pass 

        ########################################################################
        if verbose and replacements and status:
            print "PACBP REPLACEMENTS LEFT/5p DONE!!:", len(replacements)
            print cbg
            cbg.printmultiplealignment()
            omsr = cbg.overall_minimal_spanning_range()
            for node in omsr.keys(): print node,min(omsr[node]),max(omsr[node])
        ########################################################################

    else:
        # omit the 5' side of this CBG
        pass

    # return if there is something improved
    if PACBPS_CORRECTED: return True
    else:                return False

# end of function correct_pacbpgaps_nearby_omsr


def _obtain_3p_splitpos((pacbpAAsubstr,substrStaCoord,substrEndCoord),
    pacbpAAstring,pacbpStaCoord):
    """
    Get PacbP coordinate on which to split it around its 3' OMSR end
    
    @attention: Helper function in correct_pacbpgaps_nearby_omsr()

    @attention: USE WITH CARE! make sure all variables are either referring
                to the query or the sbjct, not mixed!

    @type  pacbpAAsubstr: string
    @param pacbpAAsubstr: query or sbjct substring around the 3' OMSR end
    
    @type  substrStaCoord: integer
    @param substrStaCoord: absolute AA start coordinate of the pacbpAAsubstr
    
    @type  substrEndCoord: integer
    @param substrEndCoord: absolute AA end coordinate of the pacbpAAsubstr
    
    @type  pacbpAAstring: string
    @param pacbpAAstring: complete (aligned) AA query or sbjct of this PacbP
    
    @type  pacbpStaCoord: integer
    @param pacbpStaCoord: absolute query or sbjct AA start coordinate

    @rtype  splitpos: integer
    @return splitpos: relative PacbP position where to split it
    """
    # obtain absolute split position ( in pacbpAAsubstr )
    if pacbpAAsubstr[0] == '-':
        # (most likely) a string like '----XXX' -> split on start coord
        splitpos = substrStaCoord + 1
    elif pacbpAAsubstr[-1] == '-':
        # (most likely) a string like 'XXX---' -> split on end coord -1
        splitpos = substrEndCoord - 1
    else:
        # (less likely) mixed string -> split on first gap occurrence
        splitpos =  substrStaCoord + pacbpAAsubstr.find("-")

    # absolute2relative split position
    splitpos = splitpos - pacbpStaCoord

    # correct splitpos for gaps in the aligned sequence!
    corrected = deepcopy(splitpos)
    # correct position; run while loop until < / LT !!
    while corrected-pacbpAAstring[0:corrected].count("-") < splitpos:
        corrected+=1

    # return updated splitpos
    return corrected

# end of function _obtain_3p_splitpos


def _obtain_5p_splitpos((pacbpAAsubstr,substrStaCoord,substrEndCoord),
    pacbpAAstring,pacbpStaCoord):
    """
    Get PacbP coordinate on which to split it around its 5' OMSR start
    
    @attention: Helper function in correct_pacbpgaps_nearby_omsr()

    @attention: USE WITH CARE! make sure all variables are either referring
                to the query or the sbjct, not mixed!

    @type  pacbpAAsubstr: string
    @param pacbpAAsubstr: query or sbjct substring around the 5' OMSR start
    
    @type  substrStaCoord: integer
    @param substrStaCoord: absolute AA start coordinate of the pacbpAAsubstr
    
    @type  substrEndCoord: integer
    @param substrEndCoord: absolute AA end coordinate of the pacbpAAsubstr
    
    @type  pacbpAAstring: string
    @param pacbpAAstring: complete (aligned) AA query or sbjct of this PacbP
    
    @type  pacbpStaCoord: integer
    @param pacbpStaCoord: absolute query or sbjct AA start coordinate

    @rtype  splitpos: integer
    @return splitpos: relative PacbP position where to split it
    """
    # obtain absolute split position ( in pacbpAAsubstr )
    if pacbpAAsubstr[-1] == '-':
        # (most likely) a substring like 'XXX---' -> split on first gap
        splitpos = substrStaCoord + pacbpAAsubstr.find("-") -1
    else:
        # other substring --> split on last gap
        splitpos =  substrEndCoord -\
                    ( len(pacbpAAsubstr) - pacbpAAsubstr.rfind("-") ) + 1

    # absolute2relative split position
    splitpos = splitpos - pacbpStaCoord

    # correct splitpos for gaps in the aligned sequence!
    corrected = deepcopy(splitpos)
    # correct position; run while loop until <= / LTE !!
    while corrected-pacbpAAstring[0:corrected].count("-") <= splitpos:
        corrected+=1

    # return updated splitpos
    return corrected

# end of function _obtain_5p_splitpos


def improvealignment(cbg,verbose=False,
    allow_3p_optimization=True,
    allow_5p_optimization=True,
    maximal_cbg_identity=CBG_OPTIMIZE_MAXIMAL_IDENTITY,
    clustalw_gap_size=CBG_OPTIMIZE_CLUSTALW_GAP_SIZE,
    optimization_bitscore_ratio=CBG_OPTIMIZE_MINIMAL_BITSCORE_RATIO,
    optimization_identity_ratio=CBG_OPTIMIZE_MINIMAL_IDENTITY_RATIO):
    """
    (Try to) Improve the multiple alignment of this CBG with clustalw

    @type  cbg: CodingBlockGraph
    @param cbg: CodingBlockGraph instance to optimize 

    @type  verbose: Boolean 
    @param verbose: print debugging/intermediate information to STDOUT 

    @type  allow_3p_optimization: Boolean  
    @param allow_3p_optimization: allow optimization(extension!) on the 3p side

    @type  allow_5p_optimization: Boolean
    @param allow_5p_optimization: allow optimization(extension!) on the 5p side

    @type  clustalw_gap_size: integer
    @param clustalw_gap_size: split ClustalW-multiplealignment obtained
                              PacbPs on this AA gap size

    @type  maximal_cbg_identity: float
    @param maximal_cbg_identity: do not optimize CBG when its
                                 GTG.identity() > this number

    @type  optimization_bitscore_ratio: float
    @param optimization_bitscore_ratio: only allow longer ClustalW-PacbPs when
                        at least this ratio towards the original PacbP

    @type  optimization_identity_ratio: float
    @param optimization_identity_ratio: only allow longer ClustalW-PacbPs when
                        at least this ratio towards the original PacbP

    @attention: when a CBG is flanked by a lsrCBG in the GSG, it advised to set
                        allow_*p_optimization to False

    @rtype:  Boolean 
    @return: is the CBG optimized or not
    """
    IS_IMPROVED = False

    # if both allow_*p_optimization are False -> no optimization!
    if not allow_3p_optimization and not allow_5p_optimization:
        return False

    # check if there is a likely chance that we can optimize this cbg
    # This chance is defined by parameter
    if cbg.get_genetree().identity() > maximal_cbg_identity:
        return False

    # gather current CBG data to compare before/after ClustalW optimization
    current_cbg_total_weight = cbg.total_weight()
    current_cbg_string       = str(cbg)
    current_cbg_omsr         = cbg.overall_minimal_spanning_range()
    current_cbg_maxsr        = cbg.maximal_spanning_range()

    # get the orf's sequences in a dict and do clustalw
    seqs = cbg.getorfproteinsequences()
    (_algseqs,_algm) = clustalw(seqs=seqs)

    # check if there is at least a single aligned position
    if len(_algm) == _algm.count(' '):
        return False
    # get position of the first and last aligned AA in the clustalw alignment
    firstalignedpos = 0
    finalalignedpos  = len(_algm)-1
    while _algm[firstalignedpos] == ' ': firstalignedpos+=1
    while _algm[finalalignedpos] == ' ': finalalignedpos-=1
    # increase finalalignedpos+=1 for compatibility asa list slice
    finalalignedpos+=1

    # translate ClustalW multiple alignment start & end to OMSR coordinates
    # While doing this, check if the current OMSR is fully covered by the
    # ClustalW OMSR. In case of long orf sequences and small CBGs,
    # ClustalW is likely to produce out-of-range alignments!
    newomsr = {}
    OMSR_IS_COMPLETELY_COVERED = True
    for org in seqs.keys():
        orf  = cbg.get_orfs_of_graph(organism=org)[0]
        node = cbg.node_by_organism(org)
        omsrstart = orf.protein_startPY +\
            ( firstalignedpos - _algseqs[org][0:firstalignedpos].count('-') )
        omsrend   = omsrstart + ( finalalignedpos - firstalignedpos -\
            _algseqs[org][firstalignedpos:finalalignedpos].count('-') )
        newomsr[org] = (omsrstart,omsrend)
        # get the union between cirrent CBGs OMSR and this novel OMSR
        omsrunion = current_cbg_omsr[node].intersection(
            Set(range( omsrstart, omsrend+1 ) ) )
        if len(omsrunion) < len(current_cbg_omsr[node]):
            # no, OMSR shrunk in stead of increased
            OMSR_IS_COMPLETELY_COVERED = False 
            ####################################################################
            if verbose:
                print org, len(omsrunion), " < ", len(current_cbg_omsr[node])
            ####################################################################
            # continue here; for this Organism identifier no improvement
            continue

        ########################################################################
        if verbose:
            print org, min(current_cbg_omsr[node]), max(current_cbg_omsr[node]),
            print "new:", (omsrstart,omsrend),
            print "maxsr:", min(current_cbg_maxsr[node]),
            print max(current_cbg_maxsr[node]), node, orf ,
            print orf.protein_startPY, orf.protein_endPY, len(_algseqs[org]),
            print len(_algseqs[org])-_algseqs[org].count('-'), orf.length/3
        ########################################################################

    # Check if current CBG OMSR is overlapping with clustalw OMSR
    if not OMSR_IS_COMPLETELY_COVERED:
        ########################################################################
        if verbose: print "NO improvement, ClustalW out-of-range-alignment"
        ########################################################################
        return False 


    ############################################################################
    if verbose:
        linesize=100
        print "<ClustalW obtained multiple alignment>"
        for offset in range(0,len(_algm),linesize):
            start = firstalignedpos + offset
            end   = start + linesize
            if end > finalalignedpos: end = finalalignedpos
            if offset==0 and finalalignedpos-firstalignedpos < linesize:
                end = finalalignedpos 
            for org in seqs.keys():
                print _algseqs[org][start:end], org
            print _algm[start:end]
            print ""
            if end == finalalignedpos: break
        print current_cbg_string
        cbg.printmultiplealignment()
    ############################################################################


    # loop over the pairwise organism combinations and make new pacbps
    # but only if the new OMSR extends the known OMSR.
    # In this process, split the ClustalW PacbpOrfs for gaps
    # of size clustalw_gap_size
    for orgA,orgB in cbg.pairwisecrosscombinations_organism():
        # get the current/original pacbporf
        pacbporf = cbg.get_pacbp_by_organisms(orgA,orgB)
        # check if there is novel extention and on which side
        extention = _does_clustalw_omsr_extend_pacbporf_omsr(
                        pacbporf,orgA,orgB,newomsr)

        # Check if extention is alowed in this side. This check is recommended
        # to be included for CBGs # that are neigbored/delimited by lsrCBG(s)
        if not allow_3p_optimization and extention in ['both','3p']:
            continue    # not alowed!
        if not allow_5p_optimization and extention in ['both','5p']:
            continue    # not alowed!
        if extention == None:
            continue    # not alowed!

        # get orf objects and aligned sequence parts
        orfA  = cbg.get_orfs_of_graph(organism=orgA)[0]
        orfB  = cbg.get_orfs_of_graph(organism=orgB)[0]
        seqA  = _algseqs[orgA][firstalignedpos:finalalignedpos]
        seqB  = _algseqs[orgB][firstalignedpos:finalalignedpos]
        nodeQ = cbg.node_by_organism(orgA)
        nodeS = cbg.node_by_organism(orgB)
        # make pacbp from this clustalw alignment and extend it
        alignment = ( seqA, _algm[firstalignedpos:finalalignedpos], seqB )
        coords = ( newomsr[orgA][0], newomsr[orgA][1],
                   newomsr[orgB][0], newomsr[orgB][1] )
        newpacbp = pacb.conversion.pacbp_from_clustalw(
                        alignment=alignment,coords=coords)

        # check for gaps in the clustalw alignment; if so, split them and
        # select the PacbP that overlaps with the OMSR
        if newpacbp.alignment_has_gaps(gap_size=clustalw_gap_size):
            splitted,status = pacb.splitting.split_pacb_on_gaps(
                    newpacbp,gapsize=clustalw_gap_size)
            if not status:
                # pacbp cannot be splitted for some reason.
                # Ignore it and continue with the next orgA/orgB comparison
                continue
            split_is_compatible = False
            for splittedpacbp in splitted:
                if splittedpacbp.query_start <= min(current_cbg_omsr[nodeQ]) and\
                splittedpacbp.query_end >= max(current_cbg_omsr[nodeQ]) and\
                splittedpacbp.sbjct_start <= min(current_cbg_omsr[nodeS]) and\
                splittedpacbp.sbjct_end >= max(current_cbg_omsr[nodeS]):
                    # this is the splitted PacbP that overlaps with current OMSR
                    newpacbp = splittedpacbp
                    split_is_compatible = True

                    # update the clustalwOMSR coords (newomsr)
                    newomsr[orgA] = (newpacbp.query_start,newpacbp.query_end)
                    newomsr[orgB] = (newpacbp.sbjct_start,newpacbp.sbjct_end)

                    # check - again - if this clustalw OMSR is an extention
                    # check with the ORIGINAL pacbporf!
                    extention = _does_clustalw_omsr_extend_pacbporf_omsr(
                                    pacbporf,orgA,orgB,newomsr)
                    
                    # if no extention, set split as incompatible
                    if extention == None: split_is_compatible = False

                    # break out of looping over the splits
                    break

            # check if the split was compatible with the OMSR of the current CBG
            if not split_is_compatible:
                # pacbp splits rigth through the relevant OMSR region.
                # Ignore it and continue with the next orgA/orgB comparison
                continue

        # If here, convert the (splitted) pacbp into pacbporf
        newpacbporf = pacb.conversion.pacbp2pacbporf(newpacbp,orfA,orfB)

        # now merge the clustalw pacbporf with the existing blast pacbporf
        status3p, status5p = False,False
        if extention in ['3p','both']:
            merged, status3p = pacb.merging.merge_pacbporfs(
                            pacbporf,newpacbporf,'rigth',verbose=verbose)
        if extention in ['5p','both']:
            if extention == 'both':
                # take `merged` as input pacbporf, not `pacbporf` ->
                # it is changed 4 lines higher up!
                merged, status5p = pacb.merging.merge_pacbporfs(
                            merged,newpacbporf,'left',verbose=verbose)
            else:
                merged, status5p = pacb.merging.merge_pacbporfs(
                            pacbporf,newpacbporf,'left',verbose=verbose)
           
        # Only reset the old (pacbporf) by the new (merged) if:
        #  True in (status3p, status5p) AND 
        #  orf.bitscore ratio >= optimization_bitscore_ratio AND
        #  orf.identityscore  >= optimization_identity_ratio

        try:
            # Be aware of a potential ZeroDivisionError in the bitscore ratio
            bitscore_ratio_check = ( float(merged.bitscore) /\
                    float(pacbporf.bitscore) ) >= optimization_bitscore_ratio 
        except ZeroDivisionError:
            # do not take ratio, just check if bigger.
            # by default, optimization_bitscore_ratio < 1.0, so
            # checking for gte is even a more stringent check
            bitscore_ratio_check = merged.bitscore >= pacbporf.bitscore

        try:
            # ZeroDivisionError in the identityscore can not/hardly be possible.
            # Identityscore == 0 means nothing that is alignable at all!
            # But, safety first ...
            identity_ratio_check = ( merged.identityscore /\
                    pacbporf.identityscore ) >= optimization_identity_ratio
        except ZeroDivisionError:
            # do not take ratio, just check if bigger.
            # by default, optimization_identity_ratio < 1.0, so
            # checking for gte is even a more stringent check
            identity_ratio_check = merged.identityscore >=pacbporf.identityscore 


        if True in (status3p, status5p) and\
        bitscore_ratio_check and identity_ratio_check:
            # reset 'old' pacbporf by 'merged'
            nodeQ = cbg.node_by_organism(orgA)
            nodeS = cbg.node_by_organism(orgB)
            cbg.remove_pacbp(pacbporf,nodeQ,nodeS)
            # and reset the pacbporf into the cbg
            merged.extend_pacbporf_after_stops()
            merged.source="clustalw-OPTIMIZED"
            newkey = merged.construct_unique_key(nodeQ,nodeS)
            cbg.pacbps[(newkey,nodeQ,nodeS)] = merged
            IS_IMPROVED = True
            ####################################################################
            if verbose: print "IMPROVEMENT", orgA, orgB
            ####################################################################
        else:
            ####################################################################
            if verbose: print "DISCARDED", orgA, orgB
            ####################################################################
            continue

    if IS_IMPROVED:
        # CBG is succesfully changed. Recreate cache etc.
        cbg.clear_cache()
        cbg.update_edge_weights_by_minimal_spanning_range()
        cbg.create_cache()
        ############################################################
        if verbose:
            print "### OPTIMIZED CBG", cbg
            cbg.printmultiplealignment()
        ############################################################
        # return status True -> this CBG is optimized!
        return True
    else:
        ############################################################
        if verbose: print "### no CBG optimization"
        ############################################################
        # return status False -> no CBG optimized!
        return False

# end of function improvealignment


def _does_clustalw_omsr_extend_pacbporf_omsr(pacbporf,orgA,orgB,clustalwOMSR)
    """
    Is the ClustalW assigned OMSR bigger as the OMSR of this PacbPORF?
    
    @attention: Helper functiom for improvealignment()

    @type  pacbporf: PacbPORF object
    @param pacbporf: PacbPORF object

    @type  orgA: * (string)
    @param orgA: Organism Query Indentifier of the pacbporf

    @type  orgB: * (string)
    @param orgB: Organism Sbjct Indentifier of the pacbporf

    @type  clustalwOMSR: {}
    @param clustalwOMSR: Organism Indentifiers as keys, (sta,end) value tuples

    @rtype: *
    @return: one of ('both','3p','5p',None)
    """
    # are the new multiplealignment OMSR coords bigger as the current ones?
    spos = pacbporf._get_original_alignment_pos_start()
    epos = pacbporf._get_original_alignment_pos_end()
    isextended5p = ( clustalwOMSR[orgA][0] < spos.query_pos,
                     clustalwOMSR[orgB][0] < spos.sbjct_pos )
    isextended3p = ( clustalwOMSR[orgA][1]-1 > epos.query_pos,
                     clustalwOMSR[orgB][1]-1 > epos.sbjct_pos )

    # check if there is novel extention and on which side
    extention = None
    if isextended5p == (True,True) and isextended3p == (True,True):
        extention = 'both'  # extention on both sides
    elif isextended5p == (True,True):
        extention = '5p'    # extention on 5p side alone
    elif isextended3p == (True,True):
        extention = '3p'    # extention on 3p side alone
    else:
        extention = None    # no compatible extention at all

    # return the extention variable
    return extention

# end of function _does_clustalw_omsr_extend_pabcp_omsr

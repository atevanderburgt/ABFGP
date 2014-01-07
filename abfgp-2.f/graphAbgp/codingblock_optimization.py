"""
Functions for optimization of CodingBlockGraphs (CBGs)
"""

# Abgp Imports
from lib_clustalw import clustalw
import pacb

# graphAbgp Imports
from exceptions import *

# Python Imports
from sets import Set
from copy import deepcopy

# Global variables
from settings.codingblockgraph import (
    CBG_OPTIMIZE_MAXIMAL_IDENTITY, CBG_OPTIMIZE_CLUSTALW_GAP_SIZE,
    CBG_OPTIMIZE_MINIMAL_BITSCORE_RATIO, CBG_OPTIMIZE_MINIMAL_IDENTITY_RATIO,
    )


def correct_pacbpgaps_nearby_omsr(cbg,omsr_offset=3,gap_size=3,
    omit5pside=False, omit3pside=False, verbose=False):
    """
    Correct a CBG when (large) gaps are occuring on the edges of the OMSR

    @type  cbg: CodingBlockGraph
    @param cbg: CodingBlockGraph instance (with cbgIF objects!)

    @type  omsr_offset: positive (small) integer
    @param omsr_offset: offset on both sides of the OMSR to check for gap occurrence 

    @type  omit5pside: Boolean
    @param omit5pside: Do not process the 5' side (left) of the CBG (for any reason)

    @type  omit3pside: Boolean
    @param omit3pside: Do not process the 3' side (rigth) of the CBG (for any reason)

    @type  verbose: Boolean
    @param verbose: print status/debugging messages to STDOUT

    @attention: When the CBG's OMSR is lost in the proces -> NoOverallMinimalSpanningRange exception is raised

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
            (q,m,s,coords) = pacbporf.alignmentpart_by_query( aaQpos-omsr_offset, aaQpos+1+omsr_offset )
            # if more gaps in this alignment slice then expected -> a pacbp split will follow
            # the slice is of size (2*omsr_offset)+1
            if q.count('-') >= gap_size:
                # convert (back) to pacbp and obtain position where to split
                pacbp = pacb.conversion.pacbporf2pacbp(pacbporf)
                if q[0] == '-':
                    # (most likely) a string like '----XXX' -> split on minimal query coordinate
                    splitpos = coords[0] + 1
                elif q[-1] == '-':
                    # (most likely) a string like 'XXX---' -> split on maximal query coordinate
                    # split on end coord, correct -1 because it is a range
                    splitpos = coords[1] - 1
                else:
                    # (less likely) mixed string -> split on aaQpos of OMSR
                    splitpos =  coords[0] + q.find("-")
    
                # correct splitpos by pacbp.query_start
                splitpos = splitpos - pacbp.query_start
                if verbose: print splitpos, pacbp.query_start
    
                # correct splitpos for the gaps in the aligned sequence!
                corrected = deepcopy(splitpos)
                # correct position; run while loop until < / LT !!
                while corrected - pacbp.query[0:corrected].count("-") < splitpos:
                    corrected+=1
                splitpos = corrected
                ####################################################################
                if verbose:
                    print "3p, Q", node1,node2, aaQpos, "QUERY", (q,m,s,coords)
                    print pacbp, "relative splitpos:", splitpos
                    pacbp.print_protein(_linesize=125)
                ####################################################################
                # now split the pacbp on this position and recreate the pacbporf
                pacbpL = pacb.splitting.split_pacb_on_coordinates(pacbp,(splitpos,splitpos),returnside='left')
                if pacbpL:
                    newpacbporf = pacb.conversion.pacbp2pacbporf(pacbpL,pacbporf.orfQ,pacbporf.orfS)
                    newpacbporf.extend_pacbporf_after_stops()
                    # store to replacements dict
                    replacements[(currentkey,node1,node2)] = newpacbporf
                    ####################################################################
                    if verbose:
                        print pacbpL
                        pacbpL.print_protein(_linesize=125)
                        print newpacbporf
                    ####################################################################
                    # increase counter for how much pacbps are corrected
                    PACBPS_CORRECTED+=1
    
            elif s.count('-') >= gap_size:
                # convert (back) to pacbp and obtain position where to split
                pacbp = pacb.conversion.pacbporf2pacbp(pacbporf)
                if s[0] == '-':
                    # (most likely) a string like '----XXX' -> split on minimal subject coordinate
                    splitpos = coords[2] + 1
                elif s[-1] == '-':
                    # (most likely) a string like 'XXX---' -> split on maximal subject coordinate
                    # split on end coord, correct -1 because it is a range 
                    splitpos = coords[3] - 1
                else:
                    # (less likely) mixed string -> split on aaSpos of OMSR
                    splitpos = coords[2] + s.find("-")
    
                # correct splitpos by pacbp.sbjct_start
                splitpos = splitpos - pacbp.sbjct_start
                if verbose: print splitpos, pacbp.sbjct_start
    
                # correct splitpos for the gaps in the aligned sequence!
                corrected = deepcopy(splitpos)
                # correct position; run while loop until < / LT !!
                while corrected - pacbp.sbjct[0:corrected].count("-") < splitpos:
                    corrected+=1
                splitpos = corrected
                ####################################################################
                if verbose:
                    print "3p, S", node1,node2, aaQpos, "SBJCT", (q,m,s,coords)
                    print pacbp, "relative splitpos:", splitpos
                    pacbp.print_protein(_linesize=125)
                ####################################################################
    
                # now split the pacbp on this position and recreate the pacbporf 
                pacbpL = pacb.splitting.split_pacb_on_coordinates(pacbp,(splitpos,splitpos),returnside='left')
                if pacbpL:
                    newpacbporf = pacb.conversion.pacbp2pacbporf(pacbpL,pacbporf.orfQ,pacbporf.orfS)
                    newpacbporf.extend_pacbporf_after_stops()
                    # store to replacements dict
                    replacements[(currentkey,node1,node2)] = newpacbporf
                    ####################################################################
                    if verbose:
                        print pacbpL
                        pacbpL.print_protein(_linesize=125)
                        print newpacbporf
                    ####################################################################
                    # increase counter for how much pacbps are corrected
                    PACBPS_CORRECTED+=1
                else:
                    # hmm.. the split failed. Wrong coords and/or no alignment left
                    pass
            else:
                pass

        ####################################################################
        if verbose and replacements:
            print "PACBP REPLACEMENTS RIGTH/3p TO BE DONE:", len(replacements)
            print cbg
            cbg.printmultiplealignment()
            omsr = cbg.overall_minimal_spanning_range()
            print [ (node,min(omsr[node]),max(omsr[node])) for node in omsr.keys() ]
        ####################################################################
    
        # do the replacements of 3' PacbP corrections
        status = _update_cbg_with_pacbporf_replacements(cbg,replacements)
        if status == True:
            pass    # cbg succesfully updated; still an OMSR
        elif status == False:
            # raise a NoOverallMinimalSpanningRange Exception
            raise NoOverallMinimalSpanningRange, str(cbg)
        else:
            pass 

        ####################################################################
        if verbose and replacements and status:
            print "PACBP REPLACEMENTS RIGTH/3p DONE!!:", len(replacements)
            print cbg
            cbg.printmultiplealignment()
            omsr = cbg.overall_minimal_spanning_range()
            print [ (node,min(omsr[node]),max(omsr[node])) for node in omsr.keys() ]
        ####################################################################



    else:
        # omit the 3' side of this CBG
        pass


    if not omit5pside:
        # now, repeat this function for the 5' side of the cbg
        replacements = {}
        for (currentkey,node1,node2),pacbporf in cbg.pacbps.iteritems():
            # get slice of the pacbporf around the min(OMSR) query value
            aaQpos = min(current_cbg_omsr[node1])
            (q,m,s,coords) = pacbporf.alignmentpart_by_query( aaQpos-omsr_offset, aaQpos+1+omsr_offset )
            # if more gaps in this alignment slice then expected -> a pacbp split will follow
            if q.count('-') >= gap_size:
                # convert (back) to pacbp and obtain position where to split
                pacbp = pacb.conversion.pacbporf2pacbp(pacbporf)
    
                if q[-1] == '-':
                    splitpos =  coords[0] + q.find("-") -1
                else:
                    splitpos =  coords[1] - ( len(q) - q.rfind("-") ) + 1
    
                # correct splitpos by pacbp.query_start
                splitpos = splitpos - pacbp.query_start
                if verbose: print splitpos, pacbp.query_start
    
                # correct splitpos for the gaps in the aligned sequence!
                corrected = deepcopy(splitpos)
                # correct position; run while loop until <= / LTE !!
                while corrected - pacbp.query[0:corrected].count("-") <= splitpos:
                    corrected+=1
                splitpos = corrected
                ####################################################################
                if verbose:
                    print "5p, Q", node1,node2, aaQpos, "QUERY", (q,m,s,coords)
                    print pacbp, "relative splitpos:", splitpos
                    pacbp.print_protein(_linesize=125)
                ####################################################################
    
                # now split the pacbp on this position and recreate the pacbporf
                pacbpR = pacb.splitting.split_pacb_on_coordinates(pacbp,(splitpos,splitpos),returnside='rigth')
                if pacbpR:
                    newpacbporf = pacb.conversion.pacbp2pacbporf(pacbpR,pacbporf.orfQ,pacbporf.orfS)
                    newpacbporf.extend_pacbporf_after_stops()
                    # store to replacements dict
                    replacements[(currentkey,node1,node2)] = newpacbporf
                    ####################################################################
                    if verbose:
                        print pacbpR
                        pacbpR.print_protein(_linesize=125)
                        print newpacbporf
                    ####################################################################
                    # increase counter for how much pacbps are corrected
                    PACBPS_CORRECTED+=1
    
            elif s.count('-') >= gap_size:
                # convert (back) to pacbp and obtain position where to split
                pacbp = pacb.conversion.pacbporf2pacbp(pacbporf)
    
                if s[-1] == '-':
                    splitpos =  coords[2] + s.find("-") -1
                else:
                    splitpos =  coords[3] - ( len(s) - s.rfind("-") ) + 1
    
                # correct splitpos by pacbp.sbjct_start
                splitpos = splitpos - pacbp.sbjct_start
                if verbose: print splitpos, pacbp.sbjct_start
    
                # correct splitpos for the gaps in the aligned sequence!
                corrected = deepcopy(splitpos)
                # correct position; run while loop until <= / LTE !!
                while corrected - pacbp.sbjct[0:corrected].count("-") <= splitpos:
                    corrected+=1
                splitpos = corrected
                ####################################################################
                if verbose:
                    print "5p, S", node1,node2, aaQpos, "SBJCT", (q,m,s,coords)
                    print pacbp, "relative splitpos:", splitpos,
                    print "ori:", coords[3] - ( len(s) - s.rfind("-") ) + 1 - pacbp.sbjct_start
                    pacbp.print_protein(_linesize=125)
                ####################################################################
    
                # now split the pacbp on this position and recreate the pacbporf
                pacbpR = pacb.splitting.split_pacb_on_coordinates(pacbp,(splitpos,splitpos),returnside='rigth')
                if pacbpR:
                    newpacbporf = pacb.conversion.pacbp2pacbporf(pacbpR,pacbporf.orfQ,pacbporf.orfS)
                    newpacbporf.extend_pacbporf_after_stops()
                    # store to replacements dict
                    replacements[(currentkey,node1,node2)] = newpacbporf
                    ####################################################################
                    if verbose:
                        print pacbpR
                        pacbpR.print_protein(_linesize=125)
                        print newpacbporf
                    ####################################################################
                    # increase counter for how much pacbps are corrected
                    PACBPS_CORRECTED+=1
                else:
                    # hmm.. the split failed. Wrong coords and/or no alignment left
                    pass
            else:
                pass
    
        ####################################################################
        if verbose and replacements:
            print "PACBP REPLACEMENTS LEFT/5p TO BE DONE:", len(replacements)
            print cbg
            cbg.printmultiplealignment()
            omsr = cbg.overall_minimal_spanning_range()
            print [ (node,min(omsr[node]),max(omsr[node])) for node in omsr.keys() ]
        ####################################################################

        # do the replacements of 5' PacbP corrections
        status = _update_cbg_with_pacbporf_replacements(cbg,replacements)
        if status == True:
            pass    # cbg succesfully updated; still an OMSR
        elif status == False:
            # raise a NoOverallMinimalSpanningRange Exception
            raise NoOverallMinimalSpanningRange, str(cbg)
        else:
            pass 


        ####################################################################
        if verbose and replacements and status:
            print "PACBP REPLACEMENTS LEFT/5p DONE!!:", len(replacements)
            print cbg
            cbg.printmultiplealignment()
            omsr = cbg.overall_minimal_spanning_range()
            print [ (node,min(omsr[node]),max(omsr[node])) for node in omsr.keys() ]
        ####################################################################

    else:
        # omit the 5' side of this CBG
        pass

    # return if there is something improved
    if PACBPS_CORRECTED: return True
    else:                return False

# end of function correct_pacbpgaps_nearby_omsr



def _update_cbg_with_pacbporf_replacements(cbg,replacements,verbose=False):
    """
    @rtype:  NoneBoolean
    @return: None (no repleacements to do), True (updated) or False (no OMSR!)
    """
    # do the replacements of (corrected) pacbporfs
    for (currentkey,nodeQ,nodeS),newpacbporf in replacements.iteritems():
        newkey = newpacbporf.construct_unique_key(nodeQ,nodeS)
        cbg.set_edge_weight(nodeQ,nodeS,wt=newpacbporf.bitscore)
        ####################################################################
        if verbose:
            print "OLD:", cbg.pacbps[(currentkey,nodeQ,nodeS)], currentkey
            print "NEW:", newpacbporf, newkey, nodeQ, nodeS
        ####################################################################
        del( cbg.pacbps[(currentkey,nodeQ,nodeS)] )
        cbg.pacbps[(newkey,nodeQ,nodeS)] = newpacbporf


    # Check if, after the replacement, still an OMSR is available.
    # In malafide CBGs with very poor identity% and many gaps,
    # this function likely breaks the OMSR criterion, thereby
    # providing us with a clear signal that this CBG must be removed!

    if replacements:
        cbg.clear_cache()
        hasomsr = cbg.has_overall_minimal_spanning_range()
        ####################################################################
        if verbose: print "after replacements OMSR?:", hasomsr
        ####################################################################
        if not hasomsr:
            # no OMSR left; return status False. In the parental
            # function, this is a signal for a NoOverallMinimalSpanningRange
            # Exception. In the parental function even one level higher up,
            # this Exception is the signal for removal of this cbg.
            return False

        # If here, then there is still an OMSR in the CBG
        # The CBG is succesfully optimized!
        cbg.create_cache()
        cbg.update_edge_weights_by_minimal_spanning_range()
        ####################################################################
        if verbose:
            omsr = cbg.overall_minimal_spanning_range()
            print [ (node,min(omsr[node]),max(omsr[node])) for node in omsr.keys() ]
        ####################################################################
        return True
    else:
        return None

# end of function _update_cbg_with_pacbporf_replacements




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

    @type  maximal_cbg_identity: float
    @param maximal_cbg_identity: do not optimize CBG when its GTG.identity() > this number

    @type  clustalw_gap_size: integer
    @param clustalw_gap_size: split ClustalW-multiplealignment obtained PacbPs on gap size

    @type  optimization_bitscore_ratio: float
    @param optimization_bitscore_ratio: only allow longer ClustalW-PacbPs when at least this ratio towards the original PacbP

    @type  optimization_identity_ratio: float
    @param optimization_identity_ratio: only allow longer ClustalW-PacbPs when at least this ratio towards the original PacbP

    @attention: when a CBG is flanked by a lsrCBG in the GSG, it advised to set allow_*p_optimization to False

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

    # gather data of the current cbg to compare before/after clustalw optimization
    current_cbg_total_weight = cbg.total_weight()
    current_cbg_string       = str(cbg)
    current_cbg_omsr         = cbg.overall_minimal_spanning_range()
    current_cbg_maxsr        = cbg.maximal_spanning_range()

    # get the orf's sequences in a dict and do clustalw
    seqs = cbg.getorfproteinsequences()
    (_algseqs,_algm) = clustalw(seqs=seqs)

    # check if there is at least a single aligned position
    if len(_algm) == _algm.count(' '): return False
    # get the position of the first and last aligned AA in the clustalw alignment
    firstalignedpos = 0
    finalalignedpos  = len(_algm)-1
    while _algm[firstalignedpos] == ' ': firstalignedpos+=1
    while _algm[finalalignedpos] == ' ': finalalignedpos-=1
    # increase finalalignedpos+=1 for compatibility asa list slice
    finalalignedpos+=1

    # translate clustalw multiple alignment start & end to OMSR coordinates
    # While doing this, check if the current OMSR is fully covered by the
    # ClustalW OMSR. In case of long orf sequences and small CBGs,
    # ClustalW is likely to produce out-of-range alignments!
    newomsr = {}
    OMSR_IS_COMPLETELY_COVERED = True
    for org in seqs.keys():
        orf  = cbg.get_orfs_of_graph(organism=org)[0]
        node = cbg.node_by_organism(org)
        omsrstart = orf.protein_startPY + ( firstalignedpos - _algseqs[org][0:firstalignedpos].count('-') )
        omsrend   = omsrstart + ( finalalignedpos - firstalignedpos - _algseqs[org][firstalignedpos:finalalignedpos].count('-') )
        newomsr[org] = (omsrstart,omsrend)
        omsrunion = current_cbg_omsr[node].intersection(Set(range( omsrstart, omsrend+1 ) ) )
        if len(omsrunion) < len(current_cbg_omsr[node]):
            OMSR_IS_COMPLETELY_COVERED = False 
            if verbose: print org, len(omsrunion), " < ", len(current_cbg_omsr[node])
            continue

        if verbose:
            print org, min(current_cbg_omsr[node]), max(current_cbg_omsr[node]), "new:", (omsrstart,omsrend),
            print "maxsr:", min(current_cbg_maxsr[node]), max(current_cbg_maxsr[node]),
            print node, orf , orf.protein_startPY, orf.protein_endPY,
            print len(_algseqs[org]), len(_algseqs[org])-_algseqs[org].count('-'), orf.length/3

    # Check if current CBG OMSR is overlapping with clustalw OMSR
    if not OMSR_IS_COMPLETELY_COVERED:
        if verbose: print "NO improvement, ClustalW out-of-range-alignment"
        return False 


    #######################################################################
    if verbose:
        linesize=100
        print "<ClustalW obtained multiple alignment>"
        for offset in range(0,len(_algm),linesize):
            start = firstalignedpos + offset
            end   = start + linesize
            if end > finalalignedpos: end = finalalignedpos
            if offset==0 and finalalignedpos-firstalignedpos < linesize: end = finalalignedpos 
            for org in seqs.keys():
                print _algseqs[org][start:end], org
            print _algm[start:end]
            print ""
            if end == finalalignedpos: break
        print current_cbg_string
        cbg.printmultiplealignment()
    #######################################################################

    # loop over the pairwise organism combinations and make new pacbps
    # but only if the new OMSR extends the known OMSR.
    # In this process, split the ClustalW PacbpOrfs for gaps
    # of size clustalw_gap_size
    for orgA,orgB in cbg.pairwisecrosscombinations_organism():
        # get the current/original pacbporf
        pacbporf = cbg.get_pacbp_by_organisms(orgA,orgB)
        # are the new multiplealignment OMSR coords bigger as the current ones?
        spos = pacbporf._get_original_alignment_pos_start()
        epos = pacbporf._get_original_alignment_pos_end()
        isextended5p = ( newomsr[orgA][0] < spos.query_pos, newomsr[orgB][0] < spos.sbjct_pos )
        isextended3p = ( newomsr[orgA][1]-1 > epos.query_pos, newomsr[orgB][1]-1 > epos.sbjct_pos )

        # check if there is novel extention and on which side
        extention = None
        if isextended5p == (True,True) and isextended3p == (True,True):
            extention = 'both'  # extention on both sides
        elif isextended5p == (True,True):
            extention = '5p'    # extention on 5p side alone
        elif isextended3p == (True,True):
            extention = '3p'    # extention on 3p side alone
        else:
            # no extention at all -> continue
            continue

        # Check if extention is alowed in this side
        # This check is recommended to be included for CBGs
        # that are neigbored/delimited/separated by lsrCBG(s)
        if not allow_3p_optimization and extention in ['both','3p']:
            continue    # not alowed!
        if not allow_5p_optimization and extention in ['both','5p']:
            continue    # not alowed!

        # get orf objects and aligned sequence parts
        orfA  = cbg.get_orfs_of_graph(organism=orgA)[0]
        orfB  = cbg.get_orfs_of_graph(organism=orgB)[0]
        seqA  = _algseqs[orgA][firstalignedpos:finalalignedpos]
        seqB  = _algseqs[orgB][firstalignedpos:finalalignedpos]
        nodeQ = cbg.node_by_organism(orgA)
        nodeS = cbg.node_by_organism(orgB)
        # make pacbp from this clustalw alignment and extend it
        alignment   = ( seqA, _algm[firstalignedpos:finalalignedpos], seqB )
        alignment   = _remove_gaps_from_clustalw_alignment(alignment)
        coords      = ( newomsr[orgA][0], newomsr[orgA][1], newomsr[orgB][0], newomsr[orgB][1] )
        newpacbp    = pacb.conversion.pacbp_from_clustalw(alignment=alignment,coords=coords)

        # check for gaps in the clustalw alignment; if so, split them and select the
        # pacbp that overlaps with the omsr
        if newpacbp.alignment_has_gaps(gap_size=clustalw_gap_size):
            splitted,status = pacb.splitting.split_pacb_on_gaps(newpacbp,gapsize=clustalw_gap_size)
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
                    newpacbp = splittedpacbp
                    split_is_compatible = True
                    # check - again - if this clustalw-obtained pacbp is an extention
                    newomsr[orgA] = (newpacbp.query_start,newpacbp.query_end)
                    newomsr[orgB] = (newpacbp.sbjct_start,newpacbp.sbjct_end)
                    isextended5p = ( newomsr[orgA][0] < spos.query_pos, newomsr[orgB][0] < spos.sbjct_pos )
                    isextended3p = ( newomsr[orgA][1]-1 > epos.query_pos, newomsr[orgB][1]-1 > epos.sbjct_pos )

                    # check if there is novel extention and on which side
                    extention = None
                    if isextended5p == (True,True) and isextended3p == (True,True):
                        extention = 'both'  # extention on both sides
                    elif isextended5p == (True,True):
                        extention = '5p'    # extention on 5p side alone
                    elif isextended3p == (True,True):
                        extention = '3p'    # extention on 3p side alone
                    else:
                        split_is_compatible = False
                    # break out of looping over the splits
                    break

            # check if the split was compatible with the OMSR of the current CBG
            if not split_is_compatible:
                # pacbp splits rigth through the OMSR region we are interested in.
                # Ignore it and continue with the next orgA/orgB comparison
                continue

        # convert (splitted) pacbp into pacbporf
        newpacbporf = pacb.conversion.pacbp2pacbporf(newpacbp,orfA,orfB)

        # now merge the clustalw pacbporf with the existing blast pacbporf
        status3p, status5p = False,False
        if extention in ['3p','both']:
            merged, status3p = pacb.merging.merge_pacbporfs(pacbporf,newpacbporf,'rigth',verbose=verbose)
        if extention in ['5p','both']:
            if extention == 'both':
                # do not merge `pacbporf` but `merged` -> it is changed 4 lines higher up!
                merged, status5p = pacb.merging.merge_pacbporfs(merged,newpacbporf,'left',verbose=verbose)
            else:
                merged, status5p = pacb.merging.merge_pacbporfs(pacbporf,newpacbporf,'left',verbose=verbose)

        if float(pacbporf.bitscore) == 0.0: print "ZeroDivisionError in creation!"

        # Only reset the old (pacbporf) by the new (merged) if:
        #  True in (status3p, status5p) AND 
        #  orf.bitscore ratio >= optimization_bitscore_ratio AND
        #  orf.identityscore  >= optimization_identity_ratio

        # Be aware of a potential ZeroDivisionError in the bitscore ratio
        try:
            bitscore_ratio_check = ( float(merged.bitscore) / float(pacbporf.bitscore) ) >= optimization_bitscore_ratio 
        except ZeroDivisionError:
            # do not take ratio, just check if bigger.
            # by default, optimization_bitscore_ratio < 1.0, so
            # checking for gte is even a more stringent check
            bitscore_ratio_check = merged.bitscore >= pacbporf.bitscore

        # ZeroDivisionError in the identityscore can not/hardly be possible.
        # identityscore == 0 means nothing that is alignable at all!
        # But, just be certain becasue bitscore ratio ZeroDivisionError occurred as well
        try:
            identity_ratio_check = ( merged.identityscore / pacbporf.identityscore ) >= optimization_identity_ratio
        except ZeroDivisionError:
            # do not take ratio, just check if bigger.
            # by default, optimization_identity_ratio < 1.0, so
            # checking for gte is even a more stringent check
            identity_ratio_check = merged.identityscore >= pacbporf.identityscore 


        if True in (status3p, status5p) and bitscore_ratio_check and identity_ratio_check:
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
            if verbose: print "IMPROVEMENT", orgA, orgB
            ###merged.print_protein(_linesize=150)
        else:
            if verbose: print "DISCARDED", orgA, orgB
            continue

    if IS_IMPROVED:
        cbg.clear_cache()
        cbg.update_edge_weights_by_minimal_spanning_range()
        cbg.create_cache()
        if verbose:
            print "### OPTIMIZED CBG", cbg
            cbg.printmultiplealignment()
        # return status True -> this CBG is optimized!
        return True
    else:
        if verbose: print "### no CBG optimization"
        # return status False -> no CBG optimized!
        return False

# end of function improvealignment


def _remove_gaps_from_clustalw_alignment(alignment):
    """ Remove opposing gaps in pairwise ClustalW alignment, originating from >2 sequences """
    (qseq,qsmatch,sseq)    = alignment
    length = len(qseq)
    qseqL,qsmatchL,sseqL = [],[],[]
    for pos in range(len(qseq)-1,-1,-1):
        if qseq[pos] == "-" and sseq[pos] == "-":
            if not (qseqL):
                qseqL    = list(qseq)
                qsmatchL = list(qsmatch)
                sseqL    = list(sseq)
            qseqL.pop(pos)
            qsmatchL.pop(pos)
            sseqL.pop(pos)
    if qseqL:
        return ( "".join(qseqL), "".join(qsmatchL), "".join(sseqL) )
    else:
        return alignment

# end of function _remove_gaps_from_clustalw_alignment

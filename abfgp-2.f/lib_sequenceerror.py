
# Python Imports
from os.path import join as osPathJoin
from copy import deepcopy

# Global Variable Imports
from settings.genestructure import ORF_IS_UNIGENE_LABEL

# ABGP Imports
from pacb.ordering import order_pacbporf_list
from pacb.pacbp import PacbP
from pacb.conversion import pacbp2pacbporf
from pythonlibs.uniqueness import get_random_string_tag
from abgp_logging import Logger
from lib_fasta import writeMultiFasta
from dna2prot import dna2protein
from lib_clustalw import (
    clustalw,
    clustalw_alignment_score_string,
    clustalw_gap_string,
    find_substring_in_clustalw_alignment,
    strip_alignment_for_exterior_gaps,
    )
from gene.sequenceerror import merge_orfs_with_sequenceerror
from lib_introns_pairwise import (
    _get_main_interface,
    _calc_aligned_score,
    )
from lib_cexpander import runcexpander


def check_and_cleanup_sequenceerror_interfaces(observed_introns,inwpcbgs,PCG,OPTIONS,input,verbose=False):
    """
    @attention: variable 'observed_introns' can be changed in this function
    @attention: variable 'PCG' can be changed in this function
    """
    # deal with verbosity parameter and create a logObject
    verbose= True in [OPTIONS.verbose,verbose]
    logObj = Logger(verbose=verbose)
   
    # lists to keep track of (succesfully) processed Orfid interfaces
    interfaces_done         = []
    se_inferfaces_correct   = []
    se_inferfaces_incorrect = []
    for key,vlist in observed_introns[OPTIONS.target].iteritems():
        seObj = vlist[0]
        if seObj.__class__.__name__ != "SequenceErrorConnectingOrfs": continue
        interface = _get_main_interface([ se._label for se in vlist ])
        print seObj,interface, len(interfaces_done)
        if interface in interfaces_done: continue
        # append to interfaces_done; do not perform this interface next time
        interfaces_done.append( interface )

        # get list of informants with recognized sequence-error
        se_informants = [ se._reference for se in vlist ]
        logObj.logf("CHECK SeqError:", seObj)
        logObj.logf("CHECK SeqError:", interface, len(vlist), se_informants)

        # get inwpCBGs surrounding this interface
        prevInwpCBG, nextInwpCBG = get_inwpcbgs_of_interface(interface,inwpcbgs)
        if not prevInwpCBG or not nextInwpCBG: continue

        # analyse this gap for potential sequence errors
        # OMIT all unigene informants -> unigene informants always have identical Orfs
        # and thereby might favour potential sequence error calls

        gapDataObj = analyse_gapsize_around_seqerror(interface,PCG,OPTIONS,omit_unigene_informants=True)
        # no data could have been gathered -> no SequenceError!
        if not gapDataObj: continue
        # logf gapDataObj analyses
        for row in gapDataObj.data: logObj.logf(row, "known:", row[0] in se_informants)
        for row in gapDataObj.data_ortransition: logObj.logf(row)
        logObj.logf(str(gapDataObj))

        if (gapDataObj.putative_intron_count - gapDataObj.likely_mutation_count) >= 2:
            # REJECTED as a sequence error
            se_inferfaces_incorrect.append( interface )
            continue

        # get fasta data surrounding this potential sequence error
        seqsmerged, coordsmerged, spacer_aaseq =\
            get_fasta_around_potential_sequenceerror(prevInwpCBG,nextInwpCBG,seqerrorobj=seObj)

        # to few sequences remaining to base multiple alignment upon
        # AT LEAST 2 informants must be present to allow trustable seq-error call
        if len(seqsmerged) < 3:
            # do not allow this as a seqerror interface!
            se_inferfaces_incorrect.append( interface )
            continue            

        # analyse this sequence area with ClustalW and Cexpander
        clustalwdata, cexpanderdata, scoresCLW1,scoresCLW2,scoresCEXP1,scoresCEXP2 = analyse_fasta_around_seqerror(
                seqsmerged,coordsmerged,spacer_aaseq,OPTIONS,verbose=verbose)
                
        # make a final judgement: is this a True or False sequence error?
        if scoresCLW1 == (None,None,None) and scoresCLW2 == (None,None,None):
            # failed for unknown reasons -> not a sequence error!
            se_inferfaces_incorrect.append( interface )
        elif scoresCLW2[1] > min(scoresCLW2) or scoresCLW1[1] > min(scoresCLW1):
            # accepted as a sequence error
            se_inferfaces_correct.append( interface )
        elif scoresCEXP2[1] == 1.0:
            # accepted as a sequence error
            se_inferfaces_correct.append( interface )
        elif gapDataObj.likely_mutation_ratio - gapDataObj.putative_intron_ratio >= 0.66:
            # accepted as a sequence error
            se_inferfaces_correct.append( interface )
        elif gapDataObj.putative_intron_ratio - gapDataObj.likely_mutation_ratio >= 0.66:
            # REJECTED as a sequence error
            se_inferfaces_incorrect.append( interface )
        elif gapDataObj.likely_mutation_ratio > gapDataObj.putative_intron_ratio and\
        gapDataObj.average_positive_gap_size <= 15 and\
        scoresCLW2[1]/min([scoresCLW2[0],scoresCLW2[2]]) >= 0.60:
            # accepted as a sequence error
            se_inferfaces_correct.append( interface )
        elif gapDataObj.likely_mutation_ratio >= gapDataObj.putative_intron_ratio and\
        gapDataObj.average_positive_gap_size <= 6:
            # accepted as a sequence error
            se_inferfaces_correct.append( interface )
        else:
            # REJECTED as a sequence error
            se_inferfaces_incorrect.append( interface )
            
        # if interface not added to se_inferfaces_correct -> continue with next interface
        if not interface in se_inferfaces_correct:
            # log message; it is rejected asa sequence error
            logObj.logf("SequenceError assesment", interface, interface in se_inferfaces_correct,
                len(vlist), "->", 0)
            # continue with next elem in the for loop
            continue

        # if accepted as a Sequence Error: make missed SequenceError objects 
        se_subdict = get_informant_supported_seqerrors(seObj,gapDataObj,clustalwdata,inwpcbgs,PCG,OPTIONS,input)
        if se_subdict:
            for seobj in se_subdict[se_subdict.keys()[0]]:
                if seobj._reference not in se_informants:
                    # add this SEobj with somewhat less strong support
                    observed_introns[OPTIONS.target][key].append( seobj )
            new_seqerror_cnt = len(observed_introns[OPTIONS.target][key])
        else:
            new_seqerror_cnt = len(vlist)

        # log message; it is accepted asa sequence error
        logObj.logf("SequenceError assesment", interface, interface in se_inferfaces_correct,
            len(vlist), "->", new_seqerror_cnt)

            
    # okay, done with assesment. logf message
    logObj.logf("SequenceError assesment done: ",
         "True", len(se_inferfaces_correct),
         "False", len(se_inferfaces_incorrect) )
    
    ############################################################################
    ### Cleanup the introns in the interfaces of the correct sequence errors
    ############################################################################
    GENE_MODEL_HAS_SEQUENCE_ERRORS = False
    if se_inferfaces_correct:
        # set label GENE_MODEL_HAS_SEQUENCE_ERRORS to True
        GENE_MODEL_HAS_SEQUENCE_ERRORS = True
        status = _cleanup_true_sequence_error_interfaces(
            OPTIONS.target,observed_introns,se_inferfaces_correct)
                
    # logf message of cleanup of True SequenceErrors
    logObj.logf("inferfaces of True SequenceErrors cleaned (%s)" % se_inferfaces_correct )
    
    ############################################################################
    ### cleanup the SEs in the interfaces of the incorrect sequence errors
    ############################################################################
    if se_inferfaces_incorrect:
        status = _cleanup_false_sequence_error_interfaces(
            OPTIONS.target,observed_introns,se_inferfaces_incorrect)
    
    # logf message of cleanup of False SequenceErrors
    logObj.logf("inferfaces of False SequenceErrors cleaned (%s)\n" % se_inferfaces_incorrect )

    # return boolean weather or not sequence errors are remaining
    return GENE_MODEL_HAS_SEQUENCE_ERRORS

# end of function check_and_cleanup_sequenceerror_interfaces



def check_intron_interfaces_for_hidden_sequenceerrors(observed_introns,inwpcbgs,PCG,OPTIONS,input,verbose=True):
    """
    """
    # deal with verbosity parameter and create a logObject
    verbose= True in [OPTIONS.verbose,verbose]
    logObj = Logger(verbose= True in [OPTIONS.verbose,verbose] )

    interfaces_done   = []
    known_se_interfaces = []
    new_se_inferfaces = []
    new_se_cnt = 0
    
    # find all sequence error interfaces in observed_introns
    # IMPORTANT: only interfaces that are *NOT* linked by sequence errors will be assesed
    for key,vlist in observed_introns[OPTIONS.target].iteritems():
        seObj = vlist[0]
        if seObj.__class__.__name__ != "SequenceErrorConnectingOrfs": continue
        interface = _get_main_interface([ se._label for se in vlist ])
        known_se_interfaces.append( interface )
        interfaces_done.append( interface )
    
    for pos in range(1,len(inwpcbgs)):
        # get inwpCBGs surrounding this interface
        prevInwpCBG, nextInwpCBG = inwpcbgs[pos-1:pos+1]
        interface = ( prevInwpCBG._get_target_node()[1], nextInwpCBG._get_target_node()[1] )
        if interface in interfaces_done: continue
        # append to interfaces_done; do not perform this interface next time
        interfaces_done.append( interface )
        # if interface on a single Orf -> it cannot represent a sequence error!
        if interface[0] == interface[1]: continue

        # find highest observed intron score on this interface
        max_intron_score = 0.0
        max_intron_key = (None,None)
        max_intron_len = 0
        for key,intronlist in observed_introns[OPTIONS.target].iteritems(): 
            if interface == _get_main_interface([ obj._label for obj in intronlist ]):
                iscore = _calc_aligned_score(intronlist)
                if iscore > max_intron_score:
                    max_intron_score = iscore
                    max_intron_key = key
                    max_intron_len = key[1]-key[0]
                elif iscore == max_intron_score and key[1]-key[0] < max_intron_len:
                    # take SHORTEST intron
                    max_intron_key = key
                    max_intron_len = key[1]-key[0]


        # obtain number of shared organisms/gene identifiers ~= max intron score
        shared_org_set_size = len(prevInwpCBG.organism_set().intersection(nextInwpCBG.organism_set()))

        # get OrfObjects of target and check if they can be linked
        # with a potential sequence error
        donorOrf = prevInwpCBG.get_orfs_of_graph(OPTIONS.target)[0]
        accepOrf = nextInwpCBG.get_orfs_of_graph(OPTIONS.target)[0]
        # check if this orf transition could represent an sequence error
        se = merge_orfs_with_sequenceerror(donorOrf,accepOrf)
        if not se and donorOrf.frame == accepOrf.frame:
            # not linked by a single but potential multiple seqerrors
            se = True

        prevminsr = prevInwpCBG.minimal_spanning_range(organism=OPTIONS.target)
        nextminsr = nextInwpCBG.minimal_spanning_range(organism=OPTIONS.target)
        logObj.logf("CHECK Interface:", interface, se, "introndata::",
            shared_org_set_size, max_intron_len,
            max_intron_score, "%s-%s -- %s-%s" % (
            min(prevminsr),max(prevminsr),min(nextminsr),max(nextminsr)))
        # if not linkable by a sequence error -> continue
        if se == False: continue

        # analyse this gap for potential sequence errors
        gapDataObj = analyse_gapsize_around_seqerror(interface,PCG,OPTIONS)
        # no data could have been gathered -> no SequenceError!
        if not gapDataObj: continue
        # logf gapDataObj analyses
        for row in gapDataObj.data: logObj.logf(row)
        for row in gapDataObj.data_ortransition: logObj.logf(row)
        logObj.logf(str(gapDataObj))

        # first rough filtering on gapData
        if gapDataObj.putative_intron_ratio > gapDataObj.likely_mutation_ratio:
            logObj.logf("SequenceError assesment", interface, False )
            continue 

        if gapDataObj.likely_mutation_ratio < gapDataObj.ratio_informant_presence and\
        gapDataObj.ratio_informant_absence > gapDataObj.ratio_informant_presence:
            logObj.logf("SequenceError assesment", interface, False )
            continue

        # get fasta data surrounding this potential sequence error
        seqsmerged, coordsmerged, spacer_aaseq =\
            get_fasta_around_potential_sequenceerror(prevInwpCBG,nextInwpCBG)

        # to few sequences remaining to base multiple alignment upon
        # AT LEAST 2 informants must be present to allow trustable seq-error call
        if len(seqsmerged) < 3: continue
            
        # analyse this sequence area with ClustalW and Cexpander
        clustalwdata, cexpanderdata, scoresCLW1,scoresCLW2,scoresCEXP1,scoresCEXP2 = analyse_fasta_around_seqerror(
                seqsmerged,coordsmerged,spacer_aaseq,OPTIONS,verbose=verbose)                
                
        # now make a final judgement: is this a missed sequence error?
        if clustalwdata == ( {}, "", {} ):
            # ClustalW/CEXPANDER anlyses failed -> leave as intron, REJECTED as potential sequence error
            pass
        elif gapDataObj.putative_intron_ratio > 0.0 and max_intron_score >=\
        0.75*float(shared_org_set_size-1) and max_intron_len <= 120 and\
        scoresCLW2[1]/min([scoresCLW2[0],scoresCLW2[2]]) <= 0.70:
            # leave as intron, REJECTED as potential sequence error
            # THIS IS A STOPLESS-3N intron in some species!!
            pass
        elif scoresCLW2[1] > min(scoresCLW2) or scoresCLW1[1] > min(scoresCLW1):
            # rejected as an intron, accepted as a sequence error
            new_se_inferfaces.append( interface )
        elif scoresCEXP2[1] == 1.0:
            # rejected as an intron, accepted as a sequence error
            new_se_inferfaces.append( interface )
        elif gapDataObj.likely_mutation_ratio - gapDataObj.putative_intron_ratio >= 0.66:
            # rejected as an intron, accepted as a sequence error
            new_se_inferfaces.append( interface )
        elif gapDataObj.likely_mutation_ratio >= gapDataObj.putative_intron_ratio and\
        gapDataObj.average_positive_gap_size <= 6:
            # rejected as an intron, accepted as a sequence error
            new_se_inferfaces.append( interface )
        elif gapDataObj.likely_mutation_ratio >= gapDataObj.putative_intron_ratio and\
        gapDataObj.average_positive_gap_size <= 15 and\
        scoresCLW2[1]/min([scoresCLW2[0],scoresCLW2[2]]) >= 0.60:
            # rejected as an intron, accepted as a sequence error
            new_se_inferfaces.append( interface )
        else:
            # leave as intron, REJECTED as potential sequence error
            pass

        # log message for if it is a sequence error or not
        logObj.logf("SequenceError assesment", interface, interface in new_se_inferfaces)

        # if interface not added to new_se_inferfaces -> continue with next interface
        if not interface in new_se_inferfaces: continue

        # if accepted as a missed Sequence Error: make new SequenceError objects 
        se_subdict = get_informant_supported_seqerrors(se,gapDataObj,clustalwdata,inwpcbgs,PCG,OPTIONS,input)
        # update main observed_introns data dict
        observed_introns[OPTIONS.target].update( se_subdict )
        new_se_cnt+= len(se_subdict)

    # logf message of cleanup of True SequenceErrors
    logObj.logf("SUMMARY:: New hidden SequenceErrors found: %s" % new_se_cnt )

    ############################################################################
    ### Cleanup the introns in the interfaces of the novel sequence errors
    ############################################################################
    HIDDEN_SEQUENCE_ERRORS_DISCOVERED = False
    if new_se_inferfaces:
        # set label GENE_MODEL_HAS_SEQUENCE_ERRORS to True
        HIDDEN_SEQUENCE_ERRORS_DISCOVERED = True
        status = _cleanup_true_sequence_error_interfaces(
            OPTIONS.target,observed_introns,new_se_inferfaces)

        # logf message of cleanup of True SequenceErrors
        logObj.logf("SUMMARY:: inferfaces of novel SequenceErrors cleaned (%s)" % new_se_inferfaces )
            
    # return boolean status HIDDEN_SEQUENCE_ERRORS_DISCOVERED
    return HIDDEN_SEQUENCE_ERRORS_DISCOVERED

# end of function check_intron_interfaces_for_hidden_sequenceerrors


################################################################################
### Functions to acces InwpCBGs and sequences surrounding potential SequenceErrors
################################################################################

def get_inwpcbgs_of_interface(target_orfid_interface,inwpcbgs):
    """
    Get InwardsPointingCodingBlockGraphs surrounding an OrfId interface

    @type  target_orfid_interface: tuple
    @param target_orfid_interface: (orfid,orfid) of target species interface

    @type  inwpcbgs: list
    @param inwpcbgs: ordered list of InwardsPointingCodingBlockGraph

    @rtype:  tuple
    @return: InwardsPointingCodingBlockGraph,InwardsPointingCodingBlockGraph

    @attention: returns (None,None) when interface not found
    """
    for pos in range(0,len(inwpcbgs)):
        prev = inwpcbgs[pos-1]
        next = inwpcbgs[pos]
        orfids = (prev._get_target_node()[1],next._get_target_node()[1])
        if orfids == target_orfid_interface:
            targetorg = prev._get_target_organism()
            # check if the next-next inwpCBG has a better node coverage
            # this can happen in case of poor-similarity stretches
            for nextpos in range(pos+1,len(inwpcbgs)):
                if inwpcbgs[nextpos]._get_target_node() == next._get_target_node() and\
                len(inwpcbgs[nextpos].node_set()) > len(next.node_set()) and\
                len(prev.node_set().intersection(inwpcbgs[nextpos].node_set())) >\
                len(prev.node_set().intersection(next.node_set())):
                    next = inwpcbgs[nextpos]
                else:
                    break
            # check if the prev-prev inwpCBG has a better node coverage
            # this can happen in case of poor-similarity stretches
            for prevpos in range(pos-1,-1,-1):
                if inwpcbgs[prevpos]._get_target_node() == prev._get_target_node() and\
                len(inwpcbgs[prevpos].node_set()) > len(prev.node_set()) and\
                len(next.node_set().intersection(inwpcbgs[prevpos].node_set())) >\
                len(next.node_set().intersection(prev.node_set())):
                    prev = inwpcbgs[prevpos]
                else:
                    break
            # break out of the forloop
            break
    else:
        # if EOF forloop is reached -> interface not found!
        prev,next = None, None

    # return neighboring InwardsPointingCodingBlockGraph
    return prev,next

# end of function get_inwpcbgs_of_interface


def get_fasta_around_potential_sequenceerror(prevInwpCBG,nextInwpCBG,aa_offset=25,seqerrorobj=None):
    """
    Get sequences and their coords in between 2 inwpCBGS around a potential SequenceError

    @type  prevInwpCBG: InwardsPointingCodingBlockGraph
    @param prevInwpCBG: InwardsPointingCodingBlockGraph 5' of potential seqerror

    @type  nextInwpCBG: InwardsPointingCodingBlockGraph
    @param nextInwpCBG: InwardsPointingCodingBlockGraph 3' of potential seqerror

    @rtype:  ( dict,dict )
    @return: sequences (node,sequence) and coords (node,(sta,end)) of area
    """
    seqsprev,coordsprev = prevInwpCBG.get_omsr_proteinsequences_and_coords()
    seqsnext,coordsnext = nextInwpCBG.get_omsr_proteinsequences_and_coords()
    if aa_offset:
        # strip sequences & coords untill at most aa_offset (25AA) of OSMR remains
        while min([ len(seq) for seq in seqsprev.values() ]) > aa_offset:
            for k,seq in seqsprev.iteritems():
                if len(seq) > aa_offset:
                    seqsprev[k] = seq[1:]
                    coordsprev[k].remove(min(coordsprev[k]))
        while min([ len(seq) for seq in seqsnext.values() ]) > aa_offset:
            for k,seq in seqsnext.iteritems():
                if len(seq) > aa_offset:
                    seqsnext[k] = seq[0:-1]
                    coordsnext[k].remove(max(coordsnext[k]))
                    
    prevNode = prevInwpCBG._get_target_node()
    nextNode = nextInwpCBG._get_target_node()
    seqsmerged = {}
    coordsmerged = {}
    organism = prevInwpCBG._get_target_organism()
    for nodeQ in seqsprev.keys():
        if nodeQ == prevInwpCBG._get_target_node():
            continue
        if not seqsnext.has_key(nodeQ):
            # only deal with identical informants with identical Orfs
            continue
        orfObj = prevInwpCBG.get_orfs_of_graph(organism=nodeQ[0])[0]
        seqQ = orfObj.getaas(min(coordsprev[nodeQ]),max(coordsnext[nodeQ])+1)
        strNodeQ = nodeQ[0]+"_"+str(nodeQ[1])
        # update sequences & coords
        seqsmerged[strNodeQ] = seqQ
        coordsmerged[strNodeQ] = [ min(coordsprev[nodeQ]),max(coordsnext[nodeQ]) ]

    # make the spacer for the target sequence
    prevOrfObj = prevInwpCBG.get_orfs_of_graph(organism=organism)[0]
    nextOrfObj = nextInwpCBG.get_orfs_of_graph(organism=organism)[0]
    if prevOrfObj.frame == nextOrfObj.frame:
        # separation is only by substitution(s), not indels
        orfObj = prevOrfObj # just take any from prev/next, we only need its frame
        spacer_nt_sta = (max(coordsprev[prevNode]) + 1)*3 + orfObj.frame
        spacer_nt_end = min(coordsnext[nextNode])*3 + orfObj.frame
        spacer_dnaseq = orfObj.inputgenomicsequence[spacer_nt_sta:spacer_nt_end]
        spacer_aaseq = dna2protein(spacer_dnaseq).replace("*","X")
    elif seqerrorobj.__class__.__name__ == "SequenceErrorConnectingOrfs" and\
    seqerrorobj.typeof in ['insertion','deletion'] and\
    seqerrorobj.start/3 >= max(coordsprev[prevNode]) and\
    seqerrorobj.end/3 <= min(coordsnext[nextNode]):
        # sequence error object (indel) provided which marks the ~exact
        # position of the indel
        x_aa_sta = seqerrorobj.start/3
        x_aa_end = seqerrorobj.end/3
        if x_aa_end == x_aa_sta: x_aa_sta-=1
        # get known sequence area in spacer from ORF objects
        spacerSeqPrev = prevOrfObj.getaas(max(coordsprev[prevNode])+1,x_aa_sta)
        spacerSeqNext = nextOrfObj.getaas(x_aa_end,min(coordsnext[nextNode]))
        spacer_aaseq  = spacerSeqPrev + "X" *( x_aa_end - x_aa_sta ) + spacerSeqNext

    else:
        # indel. dunno where it exists exactly, so just mask out with X's
        # do this based on the MAXSR sequence, not the MINSR sequence
        prevMAXSR = prevInwpCBG.maximal_spanning_range(organism=organism)
        nextMAXSR = nextInwpCBG.maximal_spanning_range(organism=organism)
        # check for weird inclusions (Orfs included in other Orfs)
        if not prevMAXSR.difference(nextMAXSR) or not nextMAXSR.difference(prevMAXSR):
            # take MINSR as putative MAXSR
            prevMAXSR = prevInwpCBG.minimal_spanning_range(organism=organism)
            nextMAXSR = nextInwpCBG.minimal_spanning_range(organism=organism)

        # remove overlaps in MAXSR -> make continious sequence
        while prevMAXSR.intersection(nextMAXSR) and len(prevMAXSR) > 1 and len(nextMAXSR) > 1:
            prevMAXSR.remove(max(prevMAXSR))
            nextMAXSR.remove(min(nextMAXSR))
        if max(prevMAXSR)+1 == min(nextMAXSR):
            # seamless connection -> leave room for an indel
            if len(prevMAXSR) > 1:
                prevMAXSR.remove(max(prevMAXSR))
            elif len(nextMAXSR) > 1:
                nextMAXSR.remove(min(nextMAXSR))
            else:
                # do not delete any coordinate -> empty set is the result
                pass
 

        # get known sequence area in spacer from ORF objects
        spacerSeqPrev = prevOrfObj.getaas(max(coordsprev[prevNode]),max(prevMAXSR)+1)
        spacerSeqNext = nextOrfObj.getaas(min(nextMAXSR),min(coordsnext[nextNode]))
        spacer_aa_end = min(nextMAXSR)
        spacer_aa_sta = max(prevMAXSR)+1
        spacer_aaseq  = spacerSeqPrev + "X" *( spacer_aa_end - spacer_aa_sta ) + spacerSeqNext


    # add target organism sequence to the seqsmerged dict
    seqsmerged[organism] = seqsprev[prevNode]+spacer_aaseq+seqsnext[nextNode]
    # IMPORTANT!!! When separation is by indels, this is NOT 100% correct !!!
    coordsmerged[organism] = [ min(coordsprev[prevNode]),max(coordsnext[nextNode]) ]

    # return sequences & coords
    return seqsmerged, coordsmerged, spacer_aaseq

# end of function get_fasta_around_potential_sequenceerror


################################################################################
### Functions to analyse the alignment surrounding potential SequenceErrors
################################################################################


class SeqErrorAnalysesData:
    def __init__(self):
        pass

    def __str__(self):
        lines = []
        lines.append("<%s>" % self.__class__.__name__)
        for attr in dir(self):
            if attr[0] == "_": continue
            if attr == "data": continue
            vizname = "                                   "+attr
            vizname = vizname[-35:]
            lines.append("%s: %s" % ( vizname, getattr(self,attr) ) )
        lines.append("</%s>" % self.__class__.__name__)
        return "\n".join(lines)
    # end of function __str__

# end of class SeqErrorAnalysesData


def analyse_gapsize_around_seqerror(interface,PCG,OPTIONS,omit_unigene_informants=False):
    """
    """
    informant_data = [] # list with data of informants
    absentinf_data = [] # list with informants which lack one or both Query Orfs
                        # in these cases, it is clear why it has failed;-)
    inforftra_data = [] #  list with data of informants which have orf transition
    informants_orftra = [] # list with informants which have orf transition here -> potential intron!
    
    # loop over all informants
    for informant in PCG.organism_set():
        if informant == OPTIONS.target: continue
        # get ordered PacbPORFs
        thepacbporfs = order_pacbporf_list(PCG.get_pacbps_by_organisms(OPTIONS.target,informant))
        is_unigene_informant = hasattr(thepacbporfs[0].orfS,ORF_IS_UNIGENE_LABEL)
        if thepacbporfs and omit_unigene_informants and is_unigene_informant:
            # omit unigenes if requested for
            continue
        for pos in range(1,len(thepacbporfs)):
            prevPF = thepacbporfs[pos-1]
            nextPF = thepacbporfs[pos]
            if (prevPF.orfQ.id,nextPF.orfQ.id) == interface and\
            prevPF.orfS.id == nextPF.orfS.id:
                endP = prevPF._get_original_alignment_pos_end()
                staP = nextPF._get_original_alignment_pos_start()
                gapQnt  = max([0,staP.query_dna_start - endP.query_dna_end])
                gapSnt  = max([0,staP.sbjct_dna_start - endP.sbjct_dna_end])
                absDist = abs(gapQnt-gapSnt)
                gapDist = gapQnt-gapSnt
                row = (informant, prevPF.orfS.id==nextPF.orfS.id,
                      gapQnt, gapSnt, absDist, gapDist,
                      (prevPF.orfS.id, nextPF.orfS.id) )
                informant_data.append(row)
                # okay, done. Break out of forloop
                break
            elif (prevPF.orfQ.id,nextPF.orfQ.id) == interface and\
            prevPF.orfS.id != nextPF.orfS.id:
                endP = prevPF._get_original_alignment_pos_end()
                staP = nextPF._get_original_alignment_pos_start()
                gapQnt  = max([0,staP.query_dna_start - endP.query_dna_end])
                gapSnt  = max([0,staP.sbjct_dna_start - endP.sbjct_dna_end])
                absDist = abs(gapQnt-gapSnt)
                gapDist = gapQnt-gapSnt
                row = (informant, prevPF.orfS.id==nextPF.orfS.id,
                      gapQnt, gapSnt, absDist, gapDist,
                      (prevPF.orfS.id,nextPF.orfS.id) )
                # append to both informants_orftra and inforftra_data
                informants_orftra.append(informant)
                inforftra_data.append(row)
                # okay, done. Break out of forloop
                break
        else:
            # in case it is an unigene, do *NOT* report it as an absent informant
            if not is_unigene_informant:
                absentinf_data.append(informant)

    if not informant_data: return None
                
    # store data to SeqErrorAnalysesData object
    total_cnt = float( len(informant_data) + len(absentinf_data) + len(inforftra_data) )
    seDataObj = SeqErrorAnalysesData()
    seDataObj.data = informant_data
    seDataObj.data_ortransition = inforftra_data
    seDataObj.interface = interface
    seDataObj.absent_informants = absentinf_data
    seDataObj.orftransition_informants = informants_orftra
    if not total_cnt:
        seDataObj.ratio_informant_presence      = 0.0
        seDataObj.ratio_informant_absence       = 0.0
        seDataObj.ratio_informant_orftransition = 0.0
    else:
        seDataObj.ratio_informant_presence      = len(informant_data) / total_cnt
        seDataObj.ratio_informant_absence       = len(absentinf_data) / total_cnt
        seDataObj.ratio_informant_orftransition = len(inforftra_data) / total_cnt
        
    if informant_data:
        seDataObj.min_gap_size = min([ row[5] for row in informant_data ])
        seDataObj.max_gap_size = max([ row[5] for row in informant_data ])
        seDataObj.average_gap_size = int( sum([ float(row[5]) for row in informant_data ]) / len(informant_data) )
        seDataObj.average_positive_gap_size = int( sum([ max([float(row[5]),0.0]) for row in informant_data ]) / len(informant_data) )
        seDataObj.zero_gap_size_count = [ row[5] for row in informant_data ].count(0)
        seDataObj.zero_gap_size_ratio = float(seDataObj.zero_gap_size_count) / total_cnt
        # count incidences of 1 obvious   # single mutation == -1,0,1 nt distance
                                          # row[2]: seamless connection of PacbPORFs in identical Orf informants
                                          # row[2]: seamless connection of PacbPORFs in Orf transition informants
        seDataObj.single_mutation_count = [ row[5] for row in informant_data ].count(0) +\
                                          [ abs(row[5]) for row in informant_data ].count(1) +\
                                          [ row[2] for row in informant_data ].count(0) +\
                                          [ row[2] for row in informant_data ].count(1) +\
                                          [ row[2] for row in inforftra_data ].count(0) +\
                                          [ row[2] for row in inforftra_data ].count(1)
        seDataObj.single_mutation_ratio = float(seDataObj.single_mutation_count) / total_cnt
        # all with maximal abs(gap) of 6nt is considered asa likely mutation
                                          # row[2]: seamless connection of PacbPORFs in identical Orf informants
                                          # row[2]: seamless connection of PacbPORFs in Orf transition informants
        seDataObj.likely_mutation_count = [ max([abs(float(row[5])),6.5]) for row in informant_data ].count(6.5) +\
                                          [ row[2] for row in informant_data ].count(0) +\
                                          [ row[2] for row in informant_data ].count(1) +\
                                          [ row[2] for row in inforftra_data ].count(0) +\
                                          [ row[2] for row in inforftra_data ].count(1)
        seDataObj.likely_mutation_ratio = float(seDataObj.likely_mutation_count) / total_cnt
        # puative interfaces explained by introns (in both target & informant)
        # intron length threshold taken here == POSITIVE 27nt
        # first line:  take all identical Orf informants which are explained by a TARGET intron only
        # second line: all orftransition informants which are explained by an INFORMANT intron + target NON-intron
        seDataObj.putative_intron_count = len(informant_data) - [ max([float(row[5]),26.5]) for row in informant_data ].count(26.5) +\
                                          len(inforftra_data) - [ row[2] <= 15 for row in inforftra_data ].count(True)

        seDataObj.putative_intron_ratio = float(seDataObj.putative_intron_count) / total_cnt
        seDataObj.putative_intron_informants = deepcopy(informants_orftra)
        for row in informant_data:
            if row[5] >= 27:
                seDataObj.putative_intron_informants.append(row[0])
        seDataObj.putative_intron_informants.sort()        

    # return sequence error gap analyses data
    return seDataObj
    
# end of function analyse_gapsize_around_seqerror


def analyse_fasta_around_seqerror_with_cexpander(seqsmerged,coordsmerged,seqerror_spacer_aa_seq,OPTIONS,verbose=False):
    """
    TODO...
    
    @type  seqsmerged: dictionary
    @param seqsmerged: sequences in between 2 inwpCBGS around a potential SequenceError
    
    @type  coordsmerged: dictionary
    @param coordsmerged: coord tuples (aa_sta,aa_end) as values, nodes as keys

    @type  seqerror_spacer_aa_seq: string
    @param seqerror_spacer_aa_seq: aa string around the potential SequenceError

    @type  OPTIONS: OptParse Options object
    @param OPTIONS: containing at least OPTIONS.outdir, OPTIONS.target, OPTIONS.verbose

    @type  verbose: Boolean
    @param verbose: print information to STDOUT or be silent
    """
    # deal with verbosity parameter and create a logObject
    logObj = Logger(verbose= True in [OPTIONS.verbose,verbose] )

    # check if enough sequences were applied.
    if len(seqsmerged) <= 2:
        logObj.logf(" to few sequences applied (%s)" % len(seqsmerged))
        # return None results
        return None,(None,None,None),(None,None,None)
    
    # write unaligned multifasta file
    fname_fasta = osPathJoin(OPTIONS.outdir,"%s_%s_%s.mfa" % (
            OPTIONS.target, get_random_string_tag(10), "seqerrorcexpander") )
    writeMultiFasta(seqsmerged,fname_fasta)
    cexpOutput = runcexpander(fname_fasta,output='float')
    cexpOutput.set_transferblock(OPTIONS.target)
    gap_cexpander_sta   = cexpOutput.sequence.find(seqerror_spacer_aa_seq)
    gap_cexpander_end   = gap_cexpander_sta + len(seqerror_spacer_aa_seq)
    gap_cexpander_data  = cexpOutput.binarystring[gap_cexpander_sta:gap_cexpander_end]
    prev_cexpander_data = cexpOutput.binarystring[0:gap_cexpander_sta]
    next_cexpander_data = cexpOutput.binarystring[gap_cexpander_end:]

    if type(cexpOutput.binarystring) == type(str()) or not gap_cexpander_data or\
    not prev_cexpander_data or not next_cexpander_data:
        # ERROR WILL HAPPEN!!!
        # print this to STDOUT for debugging purposes
        # analyse_fasta_around_seqerror_with_clustalw(seqsmerged,coordsmerged,seqerror_spacer_aa_seq,OPTIONS,verbose=True)
        print "seqerror_spacer_aa_seq:", seqerror_spacer_aa_seq
        print "cexpOutput.header:", cexpOutput.header 
        print "cexpOutput.sequence:", cexpOutput.sequence
        print "cexpOutput.binarystring:", len(cexpOutput.binarystring)
        # return None results
        return None,(None,None,None),(None,None,None)

    gap_uniformly_matched_ratio = float(gap_cexpander_data.count(1.0)) /\
                                    len(gap_cexpander_data)
    gap_cexpander_score_ratio   = sum(gap_cexpander_data) /\
                                    len(gap_cexpander_data)
    prev_uniformly_matched_ratio= float(prev_cexpander_data.count(1.0)) /\
                                    len(prev_cexpander_data)
    prev_cexpander_score_ratio  = sum(prev_cexpander_data) /\
                                    len(prev_cexpander_data)
    next_uniformly_matched_ratio= float(next_cexpander_data.count(1.0)) /\
                                    len(next_cexpander_data)
    next_cexpander_score_ratio  = sum(next_cexpander_data) /\
                                    len(next_cexpander_data)

    logObj.logf("CEXPANDER::", cexpOutput.header, cexpOutput.uniformly_matched_ratio(),
        len(cexpOutput._transferblocks))
    logObj.logf("%s%s" % ( " "*cexpOutput.sequence.find(seqerror_spacer_aa_seq),
        seqerror_spacer_aa_seq.replace("X","x")) )
    logObj.logf(cexpOutput.sequence.replace("X","x"))
    logObj.logf(cexpOutput.get_formatted_binarystring(),
        gap_uniformly_matched_ratio, gap_cexpander_score_ratio)
    logObj.logf("CEXPANDER prev this next :: %1.2f %1.2f %1.2f || %1.2f %1.2f %1.2f" % (
        prev_uniformly_matched_ratio, gap_uniformly_matched_ratio,
        next_uniformly_matched_ratio, prev_cexpander_score_ratio,
        gap_cexpander_score_ratio, next_cexpander_score_ratio ) )
 
    # return cexpOutput object and scores
    scores1 = ( prev_uniformly_matched_ratio,
                gap_uniformly_matched_ratio,
                next_uniformly_matched_ratio )
    scores2 = ( prev_cexpander_score_ratio,
                gap_cexpander_score_ratio,
                next_cexpander_score_ratio )
    return cexpOutput, scores1,scores2

# end of function analyse_fasta_around_seqerror_with_cexpander


def analyse_fasta_around_seqerror_with_clustalw(seqsmerged,coordsmerged,seqerror_spacer_aa_seq,OPTIONS,verbose=False):
    """
    TODO...
    
    @type  seqsmerged: dictionary
    @param seqsmerged: sequences in between 2 inwpCBGS around a potential SequenceError
    
    @type  coordsmerged: dictionary
    @param coordsmerged: coord tuples (aa_sta,aa_end) as values, nodes as keys

    @type  seqerror_spacer_aa_seq: string
    @param seqerror_spacer_aa_seq: aa string around the potential SequenceError

    @type  OPTIONS: OptParse Options object
    @param OPTIONS: containing at least OPTIONS.outdir, OPTIONS.target, OPTIONS.verbose

    @type  verbose: Boolean
    @param verbose: print information to STDOUT or be silent
    """
    # deal with verbosity parameter and create a logObject
    logObj = Logger(verbose= True in [OPTIONS.verbose,verbose] )

    # make clustalw alignment and strip exterior gaps
    (aligned,alignment) = clustalw(seqs=seqsmerged)
    aligned,alignment,coordsmerged = strip_alignment_for_exterior_gaps(
        aligned,alignment,coordsmerged)

    # find the position of seqerror_spacer_aa_seq in the UNPROCESSED alignment
    (orig_match_sta,orig_match_end) = find_substring_in_clustalw_alignment(
            seqerror_spacer_aa_seq,aligned)

    # This should bot be possible (I think), but error check here.
    # if no (orig_match_sta,orig_match_end) -> bogus data. Return None data
    if (orig_match_sta,orig_match_end) == (None,None):
        # return None results
        return ( {}, "", {} ),(None,None,None),(None,None,None)
        
    # make cexpander-like float score of ClustalW output
    clwscorelist,clwscorestring = clustalw_alignment_score_string(aligned)

    # get ClustalW gap analyses score
    clwgaplist, clwgapstring = clustalw_gap_string(aligned)

    # find the position of seqerror_spacer_aa_seq in the alignment
    (match_sta,match_end) = find_substring_in_clustalw_alignment(
            seqerror_spacer_aa_seq,aligned)

    # check if strip_alignment_for_exterior_gaps consumed the required InwpCBG
    # anchors on both sides. If the case, recreate 1AA anchor(s) and store their
    # scores as reliable
    corrected_match_sta = False
    corrected_match_end = False
    if match_end == len(alignment):
        # alignment was chopped quit extensively; no NEXT inwpCBG minsr region remains
        corrected_match_end = True
        match_end -= 1
        seqerror_spacer_aa_seq = seqerror_spacer_aa_seq[0:-1]
    if match_sta == 0:
        # alignment was chopped quit extensively; no PREV inwpCBG minsr region remains
        corrected_match_sta = True
        match_sta = 1
        seqerror_spacer_aa_seq = seqerror_spacer_aa_seq[1:]
            
    ############################################################################
    if logObj.VERBOSE or True:
        window_size = 150
        for i in range(0,len(alignment),window_size):
            for k,seq in aligned.iteritems():
                # mask sequence around the interesting spacer to lowercase
                maskedseq = seq[0:match_sta].lower() + seq[match_sta:match_end] + seq[match_end:].lower()
                logObj.logf(maskedseq[i:i+window_size],k )
            logObj.logf(alignment[i:i+window_size])
            logObj.logf(clwscorestring[i:i+window_size], "CLUSTALW similarity")
            logObj.logf(clwgapstring[i:i+window_size], "CLUSTALW gaps")
    ############################################################################

    # compare scores around the potential sequence error
    gap_clustalw_data  = clwscorelist[match_sta:match_end]
    prev_clustalw_data = clwscorelist[:match_sta]
    next_clustalw_data = clwscorelist[match_end:]

    gap_clustalw_alg_data  = alignment[match_sta:match_end]
    prev_clustalw_alg_data = alignment[:match_sta]
    next_clustalw_alg_data = alignment[match_end:]

    # create 1 AA reliable anchor(s) if alignment was chopped to much
    if corrected_match_sta:
        prev_clustalw_data = [1.0]
        prev_clustalw_alg_data = "*"
    if corrected_match_end:
        next_clustalw_data = [1.0]
        next_clustalw_alg_data = "*"
        
    #gap_uniformly_matched_ratio = float(gap_clustalw_data.count(1.0)) /\
    #                                len(gap_clustalw_data)
    #prev_uniformly_matched_ratio= float(prev_clustalw_data.count(1.0)) /\
    #                                len(prev_clustalw_data)
    #next_uniformly_matched_ratio= float(next_clustalw_data.count(1.0)) /\
    #                                len(next_clustalw_data)

    gap_uniformly_matched_ratio = float(len(gap_clustalw_alg_data)-gap_clustalw_alg_data.count(" ")) /\
                                    len(gap_clustalw_alg_data)
    prev_uniformly_matched_ratio= float(len(prev_clustalw_alg_data)-prev_clustalw_alg_data.count(" ")) /\
                                    len(prev_clustalw_alg_data)
    next_uniformly_matched_ratio= float(len(next_clustalw_alg_data)-next_clustalw_alg_data.count(" ")) /\
                                    len(next_clustalw_alg_data)

    gap_clustalw_score_ratio    = sum(gap_clustalw_data) /\
                                    len(gap_clustalw_data)
    prev_clustalw_score_ratio   = sum(prev_clustalw_data) /\
                                    len(prev_clustalw_data)
    next_clustalw_score_ratio   = sum(next_clustalw_data) /\
                                    len(next_clustalw_data)

    logObj.logf("CLUSTALW    prev this next :: %1.2f %1.2f %1.2f || %1.2f %1.2f %1.2f" % (
        prev_uniformly_matched_ratio, gap_uniformly_matched_ratio,
        next_uniformly_matched_ratio, prev_clustalw_score_ratio,
        gap_clustalw_score_ratio, next_clustalw_score_ratio ) )

                                    
    # return clustalw aligned data and scores
    scores1 = ( prev_uniformly_matched_ratio,
                gap_uniformly_matched_ratio,
                next_uniformly_matched_ratio )
    scores2 = ( prev_clustalw_score_ratio,
                gap_clustalw_score_ratio,
                next_clustalw_score_ratio )
    return (aligned,alignment,coordsmerged), scores1, scores2

# end of function analyse_fasta_around_seqerror_with_clustalw


def analyse_fasta_around_seqerror(seqsmerged,coordsmerged,seqerror_spacer_aa_seq,OPTIONS,verbose=False):
    """
    @type  seqsmerged: dictionary
    @param seqsmerged: sequences in between 2 inwpCBGS around a potential SequenceError
    
    @type  coordsmerged: dictionary
    @param coordsmerged: coord tuples (aa_sta,aa_end) as values, nodes as keys

    @type  seqerror_spacer_aa_seq: string
    @param seqerror_spacer_aa_seq: aa string around the potential SequenceError

    @type  OPTIONS: OptParse Options object
    @param OPTIONS: containing at least OPTIONS.outdir, OPTIONS.target, OPTIONS.verbose

    @type  verbose: Boolean
    @param verbose: print information to STDOUT or be silent
    """
    # deal with verbosity parameter and create a logObject
    logObj = Logger(verbose= True in [OPTIONS.verbose,verbose] )

    # perform ClustalW analyses
    (aligned,alignment,coordsmerged), scoresCLW1, scoresCLW2 =\
        analyse_fasta_around_seqerror_with_clustalw(
            seqsmerged,coordsmerged,seqerror_spacer_aa_seq,OPTIONS,verbose=False)

    # perform Cexpander analyses - sequence space might be stripped by ClustalW!!
    # for this, recreate seqsmerged
    for k,seq in aligned.iteritems(): seqsmerged[k] = seq.replace("-","")
    cexpOutput, scoresCEXP1, scoresCEXP2 = analyse_fasta_around_seqerror_with_cexpander(
            seqsmerged,coordsmerged,seqerror_spacer_aa_seq,OPTIONS,verbose=False)
 
    ############################################################################
    if logObj.VERBOSE:
        # make cexpander-like float score of ClustalW output
        clwscorelist,clwscorestring = clustalw_alignment_score_string(aligned)

        # get ClustalW gap analyses score
        clwgaplist, clwgapstring = clustalw_gap_string(aligned)

        # find the position of seqerror_spacer_aa_seq in the alignment
        (match_sta,match_end) = find_substring_in_clustalw_alignment(
                seqerror_spacer_aa_seq,aligned)

        if cexpOutput and aligned:
            # format the cexpanderbinary string by projecting the ClustalW gaps onto it
            cexpstring = cexpOutput.get_formatted_binarystring()
            for pos in range(0,len(aligned[OPTIONS.target])):
                if aligned[OPTIONS.target][pos] == "-":
                    cexpstring = cexpstring[0:pos]+"."+cexpstring[pos:]
        else:
            # cexpander FAILED for unknwon reasons
            # or ClustalW fialed because aa_spacer was consumed
            cexpstring = "x"*len(alignment)

        if aligned:
            # print multi-fasta alignment
            window_size = 100
            for i in range(0,len(alignment),window_size):
                for k,seq in aligned.iteritems():
                    # mask sequence around the interesting spacer to lowercase
                    maskedseq = seq[0:match_sta].lower() + seq[match_sta:match_end] + seq[match_end:].lower()
                    logObj.logf(maskedseq[i:i+window_size],k )
                logObj.logf(alignment[i:i+window_size])
                logObj.logf(clwscorestring[i:i+window_size], "CLUSTALW similarity")
                logObj.logf(clwgapstring[i:i+window_size], "CLUSTALW gaps")
                logObj.logf(cexpstring[i:i+window_size],"CEXPANDER")

            (   prev_uniformly_matched_ratio,
                gap_uniformly_matched_ratio,
                next_uniformly_matched_ratio ) = scoresCLW1
            (   prev_clustalw_score_ratio,
                gap_clustalw_score_ratio,
                next_clustalw_score_ratio ) = scoresCLW2

            logObj.logf("CLUSTALW  prev this next :: %1.2f %1.2f %1.2f || %1.2f %1.2f %1.2f" % (
                    prev_uniformly_matched_ratio, gap_uniformly_matched_ratio,
                    next_uniformly_matched_ratio, prev_clustalw_score_ratio,
                    gap_clustalw_score_ratio, next_clustalw_score_ratio ) )            
        else:
            logObj.logf("CLUSTALW FAILED!!!!")

        if cexpOutput:    
            (   prev_uniformly_matched_ratio,
                gap_uniformly_matched_ratio,
                next_uniformly_matched_ratio ) = scoresCEXP1
            (   prev_cexpander_score_ratio,
                gap_cexpander_score_ratio,
                next_cexpander_score_ratio ) = scoresCEXP2
            logObj.logf("CEXPANDER prev this next :: %1.2f %1.2f %1.2f || %1.2f %1.2f %1.2f" % (
                prev_uniformly_matched_ratio, gap_uniformly_matched_ratio,
                next_uniformly_matched_ratio, prev_cexpander_score_ratio,
                gap_cexpander_score_ratio, next_cexpander_score_ratio ) )
        else:
            logObj.logf("CEXPANDER FAILED!!!!")
    ############################################################################

    # return clustalw & cexpander data and all 4 score tuples
    clustalwdata = (aligned,alignment,coordsmerged)
    return (clustalwdata,cexpOutput,scoresCLW1,scoresCLW2,scoresCEXP1,scoresCEXP2)

# end of function analyse_fasta_around_seqerror


################################################################################
### Helper functions 
################################################################################


def _cleanup_true_sequence_error_interfaces(organism,observed_introns,interface_list):
    """
    @attention: variable 'observed_introns' will be changed in this function
    """
    if not interface_list: return False
    ############################################################################
    ### Cleanup the introns in the interfaces of the correct sequence errors
    ############################################################################
    keys = observed_introns[organism].keys()
    keys.sort()
    is_any_deleted = None
    orfid_sta_if_list = [ interface[0] for interface in interface_list ]
    orfid_end_if_list = [ interface[1] for interface in interface_list ]
    for key in keys:
        vlist     = observed_introns[organism][key]
        if len(vlist) == 0:
            # empty list... this is not what expected
            # in any case -> delete it here!
            del( observed_introns[organism][key] )
            is_any_deleted = True
            continue
        # take first element as example element
        example   =  vlist[0]
        interface = _get_main_interface([ se._label for se in vlist ])
        if example.__class__.__name__ != "SequenceErrorConnectingOrfs":
            if interface in interface_list:
                del( observed_introns[organism][key] )
                is_any_deleted = True
            elif interface[0] in orfid_sta_if_list:
                del( observed_introns[organism][key] )
                is_any_deleted = True
            elif interface[1] in orfid_end_if_list:
                del( observed_introns[organism][key] )
                is_any_deleted = True

    # return NoneBoolean is_any_deleted
    return is_any_deleted

# end of function _cleanup_true_sequence_error_interfaces


def _cleanup_false_sequence_error_interfaces(organism,observed_introns,interface_list):
    """
    @attention: variable 'observed_introns' will be changed in this function
    """
    if not interface_list: return False
    ############################################################################
    ### cleanup the SEs in the interfaces of the incorrect sequence errors
    ############################################################################
    keys = observed_introns[organism].keys()
    keys.sort()
    is_any_deleted = None
    for key in keys:
        vlist     = observed_introns[organism][key]
        example   =  vlist[0]
        interface = _get_main_interface([ se._label for se in vlist ])
        if interface in interface_list and\
        example.__class__.__name__ == "SequenceErrorConnectingOrfs":
            del( observed_introns[organism][key] )
            is_any_deleted = True
    # return NoneBoolean is_any_deleted
    return is_any_deleted

# end of function _cleanup_false_sequence_error_interfaces



            
def apply_indel_position_correction(seobjlist,position=None):
    """
    Apply putatitive position correction on a list of indels

    @type  seobjlist: list
    @param seobjlist: list with SequenceErrorConnectingOrfs

    @rtype  seobjlist: list
    @return seobjlist: list with SequenceErrorConnectingOrfs
    """
    if seobjlist[0].__class__.__name__ == "SequenceErrorConnectingOrfs" and\
    seobjlist[0].typeof in ['insertion','deletion']:
        coordestimation = sum([ se._estimated_nt_coord for se in seobjlist ])
        if position:
            # hard-applied position correction. Trust it!
            sta = position
        else:
            sta = coordestimation / len(seobjlist)
            # make shure estimated nt coord is in the proper frame!
            while sta % 3 != seobjlist[0]._estimated_nt_coord % 3:
                sta += 1
        end = sta + seobjlist[0].length
        new_key = (sta,end)
        # update all the individually mapped SequenceErrors
        for se in seobjlist:
            se.donor.pos = sta
            se.acceptor.pos = end
            se.start = se.donor.pos
            se.end = se.acceptor.pos
            se._estimated_nt_coord = sta
    # return list of sequence error objects
    return seobjlist

# end of function apply_indel_position_correction



def get_informant_supported_seqerrors(se,gapDataObj,clustalwdata,inwpcbgs,PCG,OPTIONS,input,verbose=True):
    """
    """
    # deal with verbosity parameter and create a logObject
    verbose= True in [OPTIONS.verbose,verbose]
    logObj = Logger(verbose=verbose)

    dict_with_seqerrors = {}
    
    # OrfIds of this target gene structure containing the SequenceError
    prev_orf_id, next_orf_id = gapDataObj.interface
    prevOrfObj = input[OPTIONS.target]['orfs'].get_orf_by_id(prev_orf_id)
    nextOrfObj = input[OPTIONS.target]['orfs'].get_orf_by_id(next_orf_id)

    # if accepted as a missed Sequence Error: make new SequenceError objects 
    if se == True:
        # shit. Multiple (!basecall!) sequence errors connecting these Orfs

        (aligned,alignment,coordsmerged) = clustalwdata
        sequence = aligned[OPTIONS.target].replace("-","").upper()
        seqerror_abs_positions = []
        seqerror_clw_positions = []
        aa_offset = min(coordsmerged[OPTIONS.target])
        offset = 0
        # find relative and absolute positions of the SequenceErrors
        while sequence.find("X",offset) > -1:
            coord = sequence.find("X",offset)
            seqerror_abs_positions.append( coord + aa_offset )
            offset = coord+1
        offset = 0
        while aligned[OPTIONS.target].upper().find("X",offset) > -1:
            coord = aligned[OPTIONS.target].upper().find("X",offset)
            seqerror_clw_positions.append( coord )
            offset = coord+1

        allow_serial_stopcodon = False
        # check for double sequence-errors: ..AAAXXAAA.. (`Orf` of 0 nt!)
        for pos in range(len(seqerror_abs_positions)-1,0,-1):
            if seqerror_abs_positions[pos] - seqerror_abs_positions[pos-1] == 1:
                seqerror_abs_positions.pop(pos)
                seqerror_clw_positions.pop(pos)
                allow_serial_stopcodon = True
    
        # logf these positions
        logObj.logf(seqerror_abs_positions,seqerror_clw_positions,"serial:",allow_serial_stopcodon)

        # OrfIds of this target gene structure containing the SequenceError
        prev_orf_id, next_orf_id = gapDataObj.interface
        
        # ordered list of concerned Orf objects
        concerned_orfs = [ prevOrfObj ]
        # get orf objects in between
        for pos in range(1,len(seqerror_abs_positions)):
            orf_sta = seqerror_abs_positions[pos-1]+1
            orf_end = seqerror_abs_positions[pos] # do NOT apply -1 correction
            clw_sta = seqerror_clw_positions[pos-1]+1
            clw_end = seqerror_clw_positions[pos] # no NOT apply -1 correction
            orf_aa_seq = aligned[OPTIONS.target][clw_sta:clw_end].replace('-','')
            orf_nt_sta = orf_sta*3 + concerned_orfs[0].frame
            orf_nt_end = orf_end*3 + concerned_orfs[0].frame
            orfs = input[OPTIONS.target]['orfs'].get_elegiable_orfs(
                # get_elegiable_orfs() filtering is in nt coords/length!
                min_orf_length = len(orf_aa_seq)*3,
                max_orf_length = len(orf_aa_seq)*3,
                max_orf_start = orf_nt_sta, min_orf_start = orf_nt_sta,
                max_orf_end = orf_nt_end, min_orf_end= orf_nt_end )
            if orfs and orfs[0].protein_sequence == orf_aa_seq:
                orfObj = orfs[0]
                logObj.logf("EXISTS:", orf_sta, orf_end, orfObj, orfObj.protein_sequence, orf_aa_seq)
            else:
                orfObj = input[OPTIONS.target]['orfs'].add_novel_orf(
                    orf_nt_sta + 1,orf_nt_end,orf_aa_seq) # use +1 offset here -> getorf reports +1 offset too
                logObj.logf("NEW:", orf_sta, orf_end, orfObj, orfObj.protein_sequence, orf_aa_seq)
            
            # append existing/novel orfObj to concerned_orfs
            concerned_orfs.append( orfObj )
    
        # append the nextOrf object to concerned_orfs
        concerned_orfs.append( nextOrfObj )
        for corf in concerned_orfs: logObj.logf("CORF::", corf)
        
        # loop over concerned Orfs -> omit first ( already known)
        for pos in range(1,len(concerned_orfs)):
            donorOrf,accepOrf = concerned_orfs[pos-1:pos+1]
            # init seObj by joining Orfs
            seObj = merge_orfs_with_sequenceerror(donorOrf,accepOrf,allow_double_stopcodon=allow_serial_stopcodon)
            seObj._organism      = OPTIONS.target 
            seObj._label         = (donorOrf.id,accepOrf.id) 
            seObj._distance      = 0
            seObj._apps_donor    = 1.0
            seObj._apps_accep    = 1.0
            seObj._gff['fsource']= "AbfgpSequenceError"
            seObj.pssm_score     = 0.5 # normal SEs get score 1 -> here lower score!

            # prepare intron key and add to dict_with_seqerrors
            key = (seObj.start,seObj.end)
            dict_with_seqerrors[key] = []

            logObj.logf(seObj, donorOrf.id,accepOrf.id, donorOrf, accepOrf,
                    donorOrf.frame,accepOrf.frame )
            
            for row in gapDataObj.data:
                # row is like ('anid', True, 82, 78, 4, 4,(orfid,orfid))
                informant = row[0]
                informant_orf_id = row[6][0]
                sbjctOrf = input[informant]['orfs'].get_orf_by_id(informant_orf_id)
                informant_node = (informant,informant_orf_id)
                informant_str_node = "%s_%s" % (informant,informant_orf_id)
                # check if informant is not more likely intron than SeqError
                if informant in gapDataObj.putative_intron_informants: continue
                # check if informant in aligned sequences; freaky PacbPs -> not present
                if informant_str_node not in aligned.keys(): continue
                
                if pos < len(concerned_orfs)-1:
                    # create PacbPORF object in between these SequenceErrors
                    clw_sta = seqerror_clw_positions[pos-1]+1
                    clw_end = seqerror_clw_positions[pos] # no NOT apply -1 correction
                    query_orf_aa_seq = aligned[OPTIONS.target][clw_sta:clw_end]
                    sbjct_orf_aa_seq = aligned[informant_str_node][clw_sta:clw_end]
                    # calculate query AA start: abs AA start + position in ClustalW alignment - gaps in ClustalW alignment
                    #query_aa_start = seqerror_abs_positions[pos-1] + seqerror_clw_positions[pos-1] - aligned[OPTIONS.target][0:clw_sta].count("-")
                    query_aa_start = seqerror_abs_positions[pos-1] + 1
                    sbjct_aa_start = min(coordsmerged[informant_str_node]) + seqerror_clw_positions[pos-1] - aligned[informant_str_node][0:clw_sta].count("-")
                    # create PacbP object
                    pacbp_input = (query_orf_aa_seq,sbjct_orf_aa_seq,query_aa_start,sbjct_aa_start)
                    pacbpObj = PacbP(input=pacbp_input)
                    pacbpObj.strip_consistent_internal_gaps()
                    pacbpObj.strip_unmatched_ends()
                    if len(pacbpObj) == 0: continue
                    # make pacbporf object
                    pacbpOrfObj = pacbp2pacbporf(pacbpObj,accepOrf,sbjctOrf)
                    pacbpOrfObj.extend_pacbporf_after_stops()
                else:
                    pacbpObj    = None
                    pacbpOrfObj = None
                    
                # deepcopy SequenceEror and append to dict_with_seqerrors
                infSeObj = deepcopy(seObj)
                infSeObj._reference = informant
                infSeObj._informant_label = (informant_orf_id,informant_orf_id)
                infSeObj._linked_to_pacbporfs = [ ]
                if pacbpOrfObj: infSeObj._linked_to_pacbporfs = [ pacbpOrfObj ]
                infSeObj._linked_to_introns = [ ]
                # append as evidence to the dict_with_seqerrors
                dict_with_seqerrors[key].append(infSeObj)
                logObj.logf( informant, infSeObj, row )
                logObj.logf( informant, pacbpObj )
                ###logObj.logf(pacbpObj.query)
                ###logObj.logf(pacbpObj.match)
                ###logObj.logf(pacbpObj.sbjct)

                
    else:
        if se.typeof != "basecall" and se.start == prevOrfObj.endPY:
            # apply sequenceerror correction offset based on ClustalW alignment
            (aligned,alignment,coordsmerged) = clustalwdata
            sequence = aligned[OPTIONS.target].replace("-","").upper()
            seqerror_abs_positions = []
            aa_offset = min(coordsmerged[OPTIONS.target])
            offset = 0
            # find relative and absolute positions of the SequenceErrors
            while sequence.find("X",offset) > -1:
                coord = sequence.find("X",offset)
                seqerror_abs_positions.append( coord + aa_offset )
                offset = coord+1
            se_aa_pos = int(round(float(sum(seqerror_abs_positions))/len(seqerror_abs_positions)))
            se_nt_pos = se_aa_pos*3 + prevOrfObj.frame

            logObj.logf("SE indel correction before:", se)
            objlist = apply_indel_position_correction( [ se ],position = se_nt_pos)
            se = objlist[0]
            logObj.logf("SE indel correction after: ", se)
    
        # single SequenceError linking these target Orfs
        se._organism      = OPTIONS.target 
        se._label         = gapDataObj.interface 
        se._distance      = 0
        se._apps_donor     = 1.0
        se._apps_accep     = 1.0
        se._gff['fsource']= "AbfgpSequenceError"
        se.pssm_score     = 0.5 # normal SEs get score 1 -> here lower score!

        # now store to dict_with_seqerrors
        key = (se.start,se.end)
        dict_with_seqerrors[key] = []
        
        for row in gapDataObj.data:
            # row is like ('anid', True, 82, 78, 4, 4, (orfid,orfid))
            informant = row[0]
            informant_interface = row[6]
            if informant not in gapDataObj.putative_intron_informants:
                infSeObj = deepcopy(se)
                infSeObj._reference = informant
                infSeObj._informant_label = informant_interface
                infSeObj._linked_to_pacbporfs = [ ]
                infSeObj._linked_to_introns = [ ]
                # append as evidence to the dict_with_seqerrors[key]
                dict_with_seqerrors[key].append( infSeObj )
                logObj.logf( informant, infSeObj, row )

    # return new sequence error data
    return dict_with_seqerrors 

# end of function get_informant_supported_seqerrors



#if False and GENE_MODEL_HAS_SEQUENCE_ERRORS:
#    dna2protlength      = len(input[OPTIONS.target]['genomeseq'])/3
#    array_algpresence   = PCG2codingarray(PCG,OPTIONS.target,dna2protlength)
#    array_algsimilarity = PCG2similarityarray(PCG,OPTIONS.target,dna2protlength)
#    
#    for key,vlist in observed_introns[OPTIONS.target].iteritems():
#        if vlist[0].__class__.__name__ == "IntronConnectingOrfs":
#            interface = _get_main_interface([ se._label for se in vlist ])
#            intronObj = vlist[0] 
#            intron_informants = [ intron._reference for intron in vlist ]
#            logf("#", "Intron:", intronObj.coords(), len(vlist), interface, intron_informants)
#            logf("#", intronObj )
#            logf("#", intronObj.branchpoint, intronObj.ppt5p!=None, intronObj.ppt3p!=None )
#
#            # make all windows of this intron size
#            window_aa_size = intronObj.length / 3
#
#            array_algsimilarity_windows = get_ordered_array_window_scores(array_algsimilarity,window_aa_size)
#            # calculate this intron window
#            sta, end = intronObj.start/3, intronObj.end/3
#            while end - sta < window_aa_size: end+=1
#            while end - sta > window_aa_size: end-=1
#            intron_total = sum(array_algsimilarity[sta:end])
#            intron_prevT = sum(array_algsimilarity[sta-window_aa_size:sta])
#            intron_nextT = sum(array_algsimilarity[end:end+window_aa_size])
#
#            offset = 0
#            while intron_total < array_algsimilarity_windows[offset]: offset+=1
#            # print info about this intron position
#            print (intron_prevT,intron_total,intron_nextT), offset,
#            print len(array_algsimilarity_windows),
#            print float(offset)/len(array_algsimilarity_windows),
#            print array_algsimilarity_windows[0], array_algsimilarity_windows[-1]
#
#            array_algpresence_windows = get_ordered_array_window_scores(array_algpresence,window_aa_size)
#            # calculate this intron window
#            sta, end = intronObj.start/3, intronObj.end/3
#            while end - sta < window_aa_size: end+=1
#            while end - sta > window_aa_size: end-=1
#            intron_total = sum(array_algpresence[sta:end])
#            intron_prevT = sum(array_algpresence[sta-window_aa_size:sta])
#            intron_nextT = sum(array_algpresence[end:end+window_aa_size])
#            offset = 0
#            while intron_total < array_algpresence_windows[offset]: offset+=1
#            # print info about this intron position
#            print (intron_prevT,intron_total,intron_nextT),
#            print offset, len(array_algpresence_windows),
#            print float(offset)/len(array_algpresence_windows),
#            print array_algpresence_windows[0], array_algpresence_windows[-1]
#
#            for pos in range(1,len(inwpcbgs)):
#                prevInwpCBG = inwpcbgs[pos-1]
#                nextInwpCBG = inwpcbgs[pos]
#                prevNode = prevInwpCBG._get_target_node()
#                nextNode = nextInwpCBG._get_target_node()
#                # continue if not on correct interface
#                if interface != ( prevNode[1], nextNode[1] ): continue
#                # continue if target nodes are identical (no seqerror!)
#                if prevNode == nextNode: continue
#                # compare informant nodes
#                prevInfNodes = [ prevInwpCBG.get_organism_nodes(org)[0] for org in prevInwpCBG.organism_set().intersection(GENE_INFORMANT_SET) ]
#                nextInfNodes = [ nextInwpCBG.get_organism_nodes(org)[0] for org in nextInwpCBG.organism_set().intersection(GENE_INFORMANT_SET) ]
#                # continue if not all informant Orf nodes shared
#                if Set(prevInfNodes).symmetric_difference(nextInfNodes): continue
#                # loop over these informants and check PacbPORFs
#                for informant_node in prevInfNodes:
#                    prevPF = prevInwpCBG.get_pacbps_by_nodes(prevNode,informant_node)[0]
#                    nextPF = nextInwpCBG.get_pacbps_by_nodes(nextNode,informant_node)[0]
#                    endP = prevPF._get_original_alignment_pos_end()
#                    staP = nextPF._get_original_alignment_pos_start()
#                    gapQnt  = max([0,staP.query_dna_start - endP.query_dna_end])
#                    gapSnt  = max([0,staP.sbjct_dna_start - endP.sbjct_dna_end])
#                    gapDist = abs(gapQnt-gapSnt)
#                    row = ( informant_node[0], gapQnt, gapSnt, gapDist )
#                    print row
#


"""

def check_sequenceerrors(observed_introns,OPTIONS,GENE_INFORMANT_SET,PCG,array_algpresence,verbose=False):
    ############################################################################
    ### Check cases where sequence-error call failed and try again with
    ### the knowledge that in some species they are predicted.
    ############################################################################
    VERBOSE = deepcopy(OPTIONS.verbose)
    if verbose: VERBOSE = verbose
    MIN_SEQUENCE_ERROR_OBSERVATION_RATIO = 0.40
    MIN_SEQUENCE_ERROR_POTENTIAL_RATIO   = 0.67
    se_inferfaces_correct   = []
    se_inferfaces_incorrect = []
    for org,data in observed_introns.iteritems():
        if org != OPTIONS.target: continue
        if not data: continue
        for key,vlist in data.iteritems():
            if vlist[0].__class__.__name__ == "SequenceErrorConnectingOrfs":
                interface = _get_main_interface([ se._label for se in vlist ])
                seObj = vlist[0] 
                se_informants = [ se._reference for se in vlist ]
                logf("#", "SeqError:", seObj, len(vlist), interface, se_informants)
                failed_se_informant_data = [] # list with data of informants without sequence error calls
                absent_se_informant_data = [] # list with informants which lack one or both Query Orfs
                                              # in these cases, it is clear why it has failed;-)
                # loop over non-retrieved informants
                for informant in GENE_INFORMANT_SET.difference(se_informants):
                    # get ordered PacbPORFs
                    thepacbporfs = order_pacbporf_list(
                        PCG.get_pacbps_by_organisms(OPTIONS.target,informant))
                    for pos in range(1,len(thepacbporfs)):
                        prevPF = thepacbporfs[pos-1]
                        nextPF = thepacbporfs[pos]
                        if (prevPF.orfQ.id,nextPF.orfQ.id) == interface and\
                        prevPF.orfS.id==nextPF.orfS.id:
                            endP = prevPF._get_original_alignment_pos_end()
                            staP = nextPF._get_original_alignment_pos_start()
                            gapQnt  = max([0,staP.query_dna_start - endP.query_dna_end])
                            gapSnt  = max([0,staP.sbjct_dna_start - endP.sbjct_dna_end])
                            gapDist = abs(gapQnt-gapSnt)
                            se_located_in_gap = seObj.start >= endP.query_dna_end and seObj.end <= staP.query_dna_start
                            row = (informant, prevPF.orfS.id==nextPF.orfS.id,
                                  gapQnt, gapSnt, gapDist, se_located_in_gap )
                            failed_se_informant_data.append(row)
                            logf("#",row)
                            # okay, done. Break out of forloop
                            break
    
                    # check ordered pacbporfs for non-existance of interface
                    if interface not in [ (thepacbporfs[pos-1].orfQ.id,thepacbporfs[pos].orfQ.id) for pos in range(1,len(thepacbporfs)) ]:
                        absent_se_informant_data.append(informant)
                            
                # check if there are ~identical gap lengths in failed_se_informant_data
                if len(failed_se_informant_data) >= 2:
                    near_identical_lengths = 0
                    for row in failed_se_informant_data:
                        (infrmt,is_ident_orf,gapQ,gapS,dist,se_in_gap) = row
                        distances = [ abs(_row[4]-dist) for _row in failed_se_informant_data ]
                        distances.sort()
                        if distances[1] <= 6:
                            near_identical_lengths+=1
                else:
                    near_identical_lengths=0
    
                # check the failed intron data
                se_obs_ratio    = float(len(vlist)) / float(len(GENE_INFORMANT_SET))
                se_id_orf_ratio = float(len(vlist)+len(failed_se_informant_data)) /\
                        float(len(GENE_INFORMANT_SET))
                se_near_ident_len = float(len(vlist)+near_identical_lengths) /\
                        float(len(GENE_INFORMANT_SET))
                # count absence of Orf as 1/2 convincing evidence
                se_absent_ratio = float(len(vlist)+len(absent_se_informant_data)/2) /\
                        float(len(GENE_INFORMANT_SET))
                se_located_in_gap_ratio = float(len(vlist) + [ row[5] for row in failed_se_informant_data ].count(True) ) /\
                        float(len(GENE_INFORMANT_SET))
    
                # check the array_algpresence around the sequence error gap
                has_presence_gaps = list(array_algpresence[seObj.start/3-10:seObj.end/3+10]).count(0) != 0
                accepted_se_id_orf_ratio   = se_id_orf_ratio >= MIN_SEQUENCE_ERROR_POTENTIAL_RATIO
                accepted_se_near_ident_len = se_near_ident_len >= MIN_SEQUENCE_ERROR_POTENTIAL_RATIO
                accepted_se_located_in_gap_ratio = se_located_in_gap_ratio >= 0.50
                
                # now asses the scores; SequenceErrors should have been called uniformly
                if se_obs_ratio >= MIN_SEQUENCE_ERROR_OBSERVATION_RATIO and\
                (se_id_orf_ratio >= MIN_SEQUENCE_ERROR_POTENTIAL_RATIO or\
                se_absent_ratio >= MIN_SEQUENCE_ERROR_POTENTIAL_RATIO or\
                se_near_ident_len >= MIN_SEQUENCE_ERROR_POTENTIAL_RATIO):
                    logf("#", "SEQUENCE ERROR!! %1.2f %1.2f %1.2f %1.2f %1.2f" % (
                            se_obs_ratio, se_id_orf_ratio, se_near_ident_len,
                            se_located_in_gap_ratio, se_absent_ratio))
                    se_inferfaces_correct.append( interface )
                elif not has_presence_gaps and accepted_se_id_orf_ratio and\
                accepted_se_near_ident_len and accepted_se_located_in_gap_ratio:
                    logf("#", "LESS OBVIOUS SEQUENCE ERROR!! %1.2f %1.2f %1.2f %1.2f %1.2f" % (
                            se_obs_ratio, se_id_orf_ratio, se_near_ident_len,
                            se_located_in_gap_ratio, se_absent_ratio))
                    se_inferfaces_correct.append( interface )
                else:
                    logf("#", "FALSE SEQUENCE ERROR!! %1.2f %1.2f %1.2f %1.2f %1.2f" % (
                            se_obs_ratio, se_id_orf_ratio, se_near_ident_len,
                            se_located_in_gap_ratio, se_absent_ratio))
                    se_inferfaces_incorrect.append( interface )
    
    
    logf("#", "SequenceError assesment done: True",
        len(se_inferfaces_correct), "False", len(se_inferfaces_incorrect) )
    
    ############################################################################
    ### cleanup the introns in the interfaces of the correct sequence errors
    ############################################################################
    GENE_MODEL_HAS_SEQUENCE_ERRORS = False
    if se_inferfaces_correct:
        # set label GENE_MODEL_HAS_SEQUENCE_ERRORS to True
        GENE_MODEL_HAS_SEQUENCE_ERRORS = True
        keys = observed_introns[OPTIONS.target].keys()
        keys.sort()
        for key in keys:
            vlist = observed_introns[OPTIONS.target][key]
            example   =  vlist[0]
            interface = _get_main_interface([ se._label for se in vlist ])
            if interface in se_inferfaces_correct and\
            example.__class__.__name__ != "SequenceErrorConnectingOrfs":
                logf("DELETING::", key, len(vlist), interface, example)
                del( observed_introns[OPTIONS.target][key] )
                
    
    logf("#", "inferfaces of True SequenceErrors cleaned (%s)" % se_inferfaces_correct )
    
    ############################################################################
    ### cleanup the SEs in the interfaces of the incorrect sequence errors
    ############################################################################
    if se_inferfaces_incorrect:
        keys = observed_introns[OPTIONS.target].keys()
        keys.sort()
        for key in keys:
            vlist = observed_introns[OPTIONS.target][key]
            example   =  vlist[0]
            interface = _get_main_interface([ se._label for se in vlist ])
            if interface in se_inferfaces_incorrect and\
            example.__class__.__name__ == "SequenceErrorConnectingOrfs":
                logf("DELETING::", key, len(vlist), interface, example)
                del( observed_introns[OPTIONS.target][key] )
    
    logf("#", "inferfaces of False SequenceErrors cleaned (%s)" % se_inferfaces_incorrect )

    return GENE_MODEL_HAS_SEQUENCE_ERRORS

# end of function ...


"""

"""
################################################################################
### In case GENE_MODEL_HAS_SEQUENCE_ERRORS, re-check all introns
### TODO TODO: leave for now, nice for version 2.
### Perfect example: CFU_833820, NRPS-like enzyme, putative
### Perfect example: DOTSE_062006
################################################################################
if GENE_MODEL_HAS_SEQUENCE_ERRORS:
    dna2protlength      = len(input[OPTIONS.target]['genomeseq'])/3
    array_algpresence   = PCG2codingarray(PCG,OPTIONS.target,dna2protlength)
    array_algsimilarity = PCG2similarityarray(PCG,OPTIONS.target,dna2protlength)

    for key,vlist in observed_introns[OPTIONS.target].iteritems():
        if True or vlist[0].__class__.__name__ == "IntronConnectingOrfs":
            interface = _get_main_interface([ se._label for se in vlist ])
            intronObj = vlist[0] 
            intron_informants = [ intron._reference for intron in vlist ]
            intron_score = (intron.length/3)*len(vlist)
            algprs_score = sum(array_algpresence[intron.start/3:intron.end/3])
            algsim_score = sum(array_algsimilarity[intron.start/3:intron.end/3])

            logf("#", "INTRONorSE", intronObj.coords(), len(vlist), interface, intron_informants)
            logf("#", intronObj )
            if intronObj.__class__.__name__ == "IntronConnectingOrfs":
                logf("#", intronObj.branchpoint, intronObj.ppt5p!=None, intronObj.ppt3p!=None, intron_score, algprs_score, algsim_score )

"""

"""

            failed_se_informant_data = [] # list with data of informants without sequence error calls
            absent_se_informant_data = [] # list with informants which lack one or both Query Orfs
                                          # in these cases, it is clear why it has failed;-)

            # loop over all informants
            for informant in GENE_INFORMANT_SET:
                # get ordered PacbPORFs
                thepacbporfs = pacb.ordering.order_pacbporf_list(
                    PCG.get_pacbps_by_organisms(OPTIONS.target,informant))
                for pos in range(1,len(thepacbporfs)):
                    prevPF = thepacbporfs[pos-1]
                    nextPF = thepacbporfs[pos]
                    if (prevPF.orfQ.id,nextPF.orfQ.id) == interface and\
                    prevPF.orfS.id==nextPF.orfS.id:
                        endP = prevPF._get_original_alignment_pos_end()
                        staP = nextPF._get_original_alignment_pos_start()
                        gapQnt  = max([0,staP.query_dna_start - endP.query_dna_end])
                        gapSnt  = max([0,staP.sbjct_dna_start - endP.sbjct_dna_end])
                        gapDist = abs(gapQnt-gapSnt)
                        ##se_located_in_gap = seObj.start >= endP.query_dna_end and seObj.end <= staP.query_dna_start
                        se_located_in_gap = None
                        row = (informant, prevPF.orfS.id==nextPF.orfS.id,
                              gapQnt, gapSnt, gapDist, se_located_in_gap )
                        failed_se_informant_data.append(row)
                        logf("#",row)
                        # okay, done. Break out of forloop
                        break
                else:
                    absent_se_informant_data.append(informant)
            # print info on absent_se_informant_data
            if absent_se_informant_data: logf("#","absent informants:",absent_se_informant_data)

"""
""" 
            for pos in range(0,len(inwpcbgs)):
                prev = inwpcbgs[pos-1]
                next = inwpcbgs[pos]
                if (prev._get_target_node()[1],next._get_target_node()[1]) == interface:
                    # check if the next-next inwpCBG has a better node coverage
                    # this can happen in case of poor-similarity stretches
                    for nextpos in range(pos+1,len(inwpcbgs)):
                        if inwpcbgs[nextpos]._get_target_node() == next._get_target_node() and\
                        len(inwpcbgs[nextpos].node_set()) > len(next.node_set()) and\
                        len(prev.node_set().intersection(inwpcbgs[nextpos].node_set())) >\
                        len(prev.node_set().intersection(next.node_set())):
                            next = inwpcbgs[nextpos]
                            print "TAKING NEXT!!"
                        else:
                            break
                    # check if the prev-prev inwpCBG has a better node coverage
                    # this can happen in case of poor-similarity stretches
                    for prevpos in range(pos-1,-1,-1):
                        if inwpcbgs[prevpos]._get_target_node() == prev._get_target_node() and\
                        len(inwpcbgs[prevpos].node_set()) > len(prev.node_set()) and\
                        len(next.node_set().intersection(inwpcbgs[prevpos].node_set())) >\
                        len(next.node_set().intersection(prev.node_set())):
                            prev = inwpcbgs[prevpos]
                            print "TAKING PREV!!"
                        else:
                            break
"""
   
"""             
                    seqsprev,coordsprev = prev.get_omsr_proteinsequences_and_coords()
                    seqsnext,coordsnext = next.get_omsr_proteinsequences_and_coords()
                    seqsmerged = {}
                    for nodeQ in seqsprev.keys():
                        if nodeQ[0] == OPTIONS.target:  continue
                        if not seqsnext.has_key(nodeQ): continue
                        orfObj = prev.get_orfs_of_graph(organism=nodeQ[0])[0]
                        seqQ = orfObj.getaas(min(coordsprev[nodeQ]),max(coordsnext[nodeQ])+1)
                        seqsmerged[nodeQ[0]+"_"+str(nodeQ[1])] = seqQ
                    # make the spacer for the target sequence
                    spacer_aaseq = "X"*( min(coordsnext[next._get_target_node()])-max(coordsprev[prev._get_target_node()]) )
                    if prev.get_orfs_of_graph(organism=OPTIONS.target)[0].frame == next.get_orfs_of_graph(organism=OPTIONS.target)[0].frame:
                        # separation is only by substitution(s), not indels
                        orfObj = prev.get_orfs_of_graph(organism=OPTIONS.target)[0]
                        spacer_nt_sta = (max(coordsprev[prev._get_target_node()]) + 1)*3 + orfObj.frame
                        spacer_nt_end = min(coordsnext[next._get_target_node()])*3 + orfObj.frame
                        spacer_dnaseq = orfObj.inputgenomicsequence[spacer_nt_sta:spacer_nt_end]
                        spacer_aaseq = dna2protein(spacer_dnaseq).replace("*","X")
                        print "spacer_aaseq::", spacer_aaseq
                    seqsmerged[OPTIONS.target] = seqsprev[prev._get_target_node()] + spacer_aaseq +  seqsnext[next._get_target_node()]
                    (aligned,alignment) = clustalw(seqs=seqsmerged)
                    window_size = 150
                    for i in range(0,len(alignment),window_size):
                        for k,seq in aligned.iteritems():
                            print seq[i:i+window_size].replace("X","x"),"\t",k 
                        print alignment[i:pos+window_size]
                        print ""

                    fname_fasta = osPathJoin(OPTIONS.outdir,"intronorseqerror.%s.%s.%s.mfa" % (OPTIONS.target,pos-1,pos))
                    writeMultiFasta(seqsmerged,fname_fasta)
                    from parsers.cexpander import runcexpander
                    cexpOutput = runcexpander(fname_fasta,output='float')
                    cexpOutput.set_transferblock(OPTIONS.target)
                    gap_cexpander_sta   = cexpOutput.sequence.find(spacer_aaseq)
                    gap_cexpander_end   = gap_cexpander_sta + len(spacer_aaseq)
                    gap_cexpander_data  = cexpOutput.binarystring[gap_cexpander_sta:gap_cexpander_end]
                    prev_cexpander_data = cexpOutput.binarystring[0:gap_cexpander_sta]
                    next_cexpander_data = cexpOutput.binarystring[gap_cexpander_end:]
                    gap_uniformly_matched_ratio = float(gap_cexpander_data.count(1.0)) / len(gap_cexpander_data)
                    gap_cexpander_score_ratio   = sum(gap_cexpander_data) / len(gap_cexpander_data)
                    prev_uniformly_matched_ratio= float(prev_cexpander_data.count(1.0)) / len(prev_cexpander_data)
                    prev_cexpander_score_ratio  = sum(prev_cexpander_data) / len(prev_cexpander_data)
                    next_uniformly_matched_ratio= float(next_cexpander_data.count(1.0)) / len(next_cexpander_data)
                    next_cexpander_score_ratio  = sum(next_cexpander_data) / len(next_cexpander_data)

                    logf("#","CEXPANDER::", cexpOutput.header, cexpOutput.get_uniform_positions(), cexpOutput.uniformly_matched_ratio(), len(cexpOutput._transferblocks))
                    logf("#","%s%s" % ( " "*cexpOutput.sequence.find(spacer_aaseq), spacer_aaseq.replace("X","x")) )
                    logf("#",cexpOutput.sequence.replace("X","x"))
                    logf("#",cexpOutput.get_formatted_binarystring(), gap_uniformly_matched_ratio, gap_cexpander_score_ratio)
                    logf("#","prev this next :: %1.2f %1.2f %1.2f || %1.2f %1.2f %1.2f" % (
                            prev_uniformly_matched_ratio, gap_uniformly_matched_ratio, next_uniformly_matched_ratio,
                            prev_cexpander_score_ratio, gap_cexpander_score_ratio, next_cexpander_score_ratio ) )
                    break
"""






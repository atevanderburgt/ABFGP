"""
Functions for checking the given (annotated) unigene structure on a DNA sequence of Orfs

    def correct_unigene_for_utrs(unigene_gff_list,start_codon_gff=(),stop_codon_gff=(),dnaseqfname=None):
    def get_annotated_genes_start_codon(gfflist):
    def get_annotated_genes_stop_codon(gfflist):
    def unigeneconfirmation(input,verbose=False):

Imported from abgp_geneconfirmation.py:

    def order_gff_list(gfflist):
    def dna2protein(_sequence,_frame):
    def genestructuregff2orfs(gfftracks,orfs,verbose=False):
    def genestructurehasproperstartcodon(gfflist,dnasequence,fmethod,verbose=False):
    def genestructurehasproperstopcodon(gfflist,dnasequence,fmethod,verbose=False):
    def rungetorf(input):
    def parseinputgff(input):
    def confirmcanonicalsplicesites(sequence,gfflist,exon_fmethod=None,verbose=False):
    def geneconfirmation(input,verbose=False):

"""

# import Python things
import os
from re import finditer
from copy import deepcopy

# import Abgp functions
from abgp_geneconfirmation import *
from abgp_warnings import UniGeneStructureIsNotMappableOnOrfsWarning, GeneStructureIsNotMappableOnOrfsWarning
from gff import gffs2coordset, filtergffs4fmethod, gffs2txt
from lib_fasta import parseSingleFastaHeaderFromFile, parseSingleFasta
from gene.start import score_tss

# import ABGP global variables
from settings.gff import GFF_UGEXON_FMETHOD, GFF_UG3UTREXON_FMETHOD, GFF_UG5UTREXON_FMETHOD
from settings.gff import GFF_GENESTART_FMETHOD, GFF_GENESTOP_FMETHOD, GFF_CDS_FMETHOD
from settings.executables import EXECUTABLE_UNIGENEANNOTATION, PYTHON_PATH
from settings.translationalstartsites import IC_TSS_PATTERN_OFFSET


def correct_unigene_for_utrs(unigene_gff_list,start_codon_gff=(),stop_codon_gff=(),
    minimal_likely_tss_pssm_score=3.0,
    shift_tss_pssm_score_ratio=4.0,
    dnaseqfname=None, verbose=False ):
    """
    Check if unigene contains evidence for non-coding UTRs and if so, correct

    @type  unigene_gff_list: list
    @param unigene_gff_list: list with uncorrected unigene gff tuples

    @type  start_codon_gff: tuple
    @param start_codon_gff: tuple representing the (annotated) protein's start codon
    
    @type  stop_codon_gff: tuple
    @param stop_codon_gff: tuple representing the (annotated) protein's stop codon

    @type  dnaseqfname: string (or None)
    @param dnaseqfname: filename of DNA sequence corresponding to the unigene's GFF

    @type  minimal_likely_tss_pssm_score: float
    @param minimal_likely_tss_pssm_score: minimal (likely) PSSM score of the TSS

    @type  shift_tss_pssm_score_ratio: float
    @param shift_tss_pssm_score_ratio: shift TSS downstream when ratio is exceeded

    @type  verbose: Boolean
    @param verbose: print discrepancies to STDOUT (True) or not (False, default)

    @rtype:  list of (gff) tuples + typeofunigene string
    @return: list with corrected gff tuples + string

    @attention: Global variable GFF_UGEXON_FMETHOD  is required for this function
    @attention: Global variable GFF_UG3UTREXON_FMETHOD  is required for this function
    @attention: Global variable GFF_UG5UTREXON_FMETHOD  is required for this function
    """
    # return list with corrected unigene tracks
    return_unigene_gff_list = []
    start_codon_pos = None
    stop_codon_pos  = None
    typeofunigene   = None
    # make sets of unigene coordinates
    unigene_coordinate_set = gffs2coordset(unigene_gff_list,fmethod=[GFF_UGEXON_FMETHOD])

    if dnaseqfname:
        # print unigene structure annotation
        unigeneexons = filtergffs4fmethod(unigene_gff_list,GFF_UGEXON_FMETHOD)
        unigeneexons.sort()
        # replace fasta header for correct recognition
        #header,descr = parseSingleFastaHeaderFromFile(dnaseqfname)
        header,dnaseq,descr = parseSingleFasta(open(dnaseqfname).readlines())
        for i in range(0,len(unigeneexons)):
            gff = list(unigeneexons[i])
            gff[0] = header
            # correct for negative coordinate. This can happen in case
            # the unigene sticks out of the genelocus
            if gff[3] <= 0: gff[3] = 1
            unigeneexons[i] = tuple(gff)

        # run unigeneannotation command
        command = "%s %s %s" % (PYTHON_PATH,EXECUTABLE_UNIGENEANNOTATION,dnaseqfname)
        ci,co = os.popen2(command)
        ci.write( gffs2txt(unigeneexons) )
        ci.close()
        ugannotation = co.read().strip().split("\t")
        co.close()
        typeofunigene = ugannotation[0]
        # abstract coordinates of start and stop codon
        # from unigene annotation
        try:    start_codon_pos = int(ugannotation[5])
        except: start_codon_pos = None
        try:    stop_codon_pos = int(ugannotation[6])
        except: stop_codon_pos = None

        ################################################################
        if verbose:
            for track in unigene_gff_list: print track
            print ugannotation, start_codon_pos, stop_codon_pos
            print "given ATG:", start_codon_gff
            print "given TGA:", stop_codon_gff
        ################################################################

        if start_codon_pos:
            # check if the PythonRegex obtained Methionine is the most
            # likely TSS. When a far better one is available -> shift
            # the TSS downstream (5p->3p) to this better TSS.
            startcodons = []
            for gffpos in range(start_codon_pos,int(unigeneexons[0][4]),3):
                if dnaseq[gffpos-1:gffpos-1+3].upper() == 'ATG':
                    tssSta = gffpos-1-IC_TSS_PATTERN_OFFSET[0]
                    tssEnd = gffpos-1+3+IC_TSS_PATTERN_OFFSET[1]
                    tssSeq = dnaseq[tssSta:tssEnd]
                    tssSco = score_tss(tssSeq)
                    #print 'ATG', gffpos, "%1.2f" % tssSco
                    startcodons.append( ( tssSco, gffpos ) )
            # check if there are >1 start codon posibilities
            if len(startcodons) > 1 and startcodons[0][0] <\
            minimal_likely_tss_pssm_score:
                for score, gffpos in startcodons[1:]:
                    if score >= minimal_likely_tss_pssm_score and\
                    abs( score / startcodons[0][0] ) > shift_tss_pssm_score_ratio:
                        start_codon_pos = gffpos
                        #print "TSS pos SHIFTED", startcodons[0][1], "->", gffpos
                        # break out after first shift; this is now *THE* TSS
                        break
        elif start_codon_gff:
            # unigene is a fragment or other transcript without
            # likely ATG. Fortunately, ATG is applied from the given
            # gene structure. Take this one.
            start_codon_pos = int(start_codon_gff[3])
        else:
            # NO start_codon_pos available -> unigene fragment!
            pass

    elif start_codon_gff or stop_codon_gff:
        typeofunigene = None # unknown -> no unigeneannotation
        # No dna sequence is applied to verify the ATG/TGA
        # positions of the unigene by unigene annotation.
        # Abstract coordinates of start and/or stop codons
        # from the given coordinates (from the gene's annotation)
        if start_codon_gff:
            start_codon_pos = int(start_codon_gff[3])
        if stop_codon_gff:
            stop_codon_pos = int(stop_codon_gff[4])
    else:
        typeofunigene = None # unknown -> no unigeneannotation
        ########################################################
        if verbose:
            print "NONE GIVEN seq/sta/end:", dnaseqfname 
            print "gff ATG:", start_codon_gff
            print "gff TGA:", stop_codon_gff
        ########################################################
        # no anchors applied in terms of start/stop sites
        # TODO future update: find or predict the putative orf
        # of this unigene. That specific functionallity should
        # NOT be placed in this function!
        # for the time being, just return the input gff list.
        return unigene_gff_list, typeofunigene 

    # create an unigene stop codon track when in unigene_coordinate_set
    if stop_codon_pos and stop_codon_pos in unigene_coordinate_set:
        # make a deepcopy of the first unigene exon track and make a list of it
        newgff = list( deepcopy( unigene_gff_list[0] ) )
        # update the coordinates
        newgff[2] = 'UGstop'
        newgff[3] = stop_codon_pos - 2
        newgff[4] = stop_codon_pos
        return_unigene_gff_list.append( tuple(newgff) )

    # CORRECT the unigene_coordinate_set for 5p nucleotides
    ignore_5p_coords = []
    if start_codon_pos != None and start_codon_pos in\
    unigene_coordinate_set and min(unigene_coordinate_set) < start_codon_pos:
        if verbose: print "CREATE 5pUTR for ug:", typeofunigene
        # yes, there is a 5p unigene alignment part
        for coord in unigene_coordinate_set:
            if coord < start_codon_pos:
                # append to the ignore_5p_coords list
                ignore_5p_coords.append( coord )
        # remove from the unigene coord set
        for coord in ignore_5p_coords:
            unigene_coordinate_set.remove(coord)

    # CORRECT the unigene_coordinate_set for 3p nucleotides
    ignore_3p_coords = []
    if stop_codon_pos != None and stop_codon_pos in\
    unigene_coordinate_set and max(unigene_coordinate_set) > stop_codon_pos:
        if verbose: print "CREATE 3pUTR for ug:", typeofunigene
        # yes, there is a 3p unigene alignment part
        for coord in unigene_coordinate_set:
            if coord > stop_codon_pos:
                # append to the ignore_5p_coords list
                ignore_3p_coords.append( coord )
        # remove from the unigene coord set
        for coord in ignore_3p_coords:
            unigene_coordinate_set.remove(coord)
        #### remove the stop codon position too
        ###unigene_coordinate_set.remove(stop_codon_pos-2)
        ###unigene_coordinate_set.remove(stop_codon_pos-1)
        ###unigene_coordinate_set.remove(stop_codon_pos)

    # make (new) UGExon tracks, corrected for UTRS, if needed
    if not (ignore_5p_coords or ignore_3p_coords) and unigene_coordinate_set:
        # no utrs available; just set the input to the output list
        return_unigene_gff_list.extend(unigene_gff_list)

    elif (ignore_5p_coords or ignore_3p_coords) and unigene_coordinate_set:
        # create new gff tracks for unigene exons
        unigene_exon_coords = list(unigene_coordinate_set)
        unigene_exon_coords.sort()
        track_coords  = [ [ unigene_exon_coords[0] ] ]
        for coord in unigene_exon_coords[1:]:
            if coord == max(track_coords[-1])+1:
                track_coords[-1].append(coord)
            else:
                track_coords.append( [ coord ] )
        for track in track_coords:
            # make a deepcopy of the first unigene exon track and make a list of it
            newgff = list( deepcopy( unigene_gff_list[0] ) )
            # update the coordinates
            newgff[3] = min(track)
            newgff[4] = max(track)
            # and append to the new return unigene gff list
            return_unigene_gff_list.append( tuple(newgff) )

        # make UTR5UGExon track if it exists
        if ignore_5p_coords:
            ignore_5p_coords.sort()
            tracks  = [ [ ignore_5p_coords[0] ] ]
            for coord in ignore_5p_coords[1:]:
                if coord == max(tracks[-1])+1:
                    tracks[-1].append(coord)
                else:
                    tracks.append( [ coord ] )
            # reverse tracks; if there are >1, inserting in the
            # return list will guarantee the correct order
            tracks.reverse()
            for track in tracks:
                # make a deepcopy of the first unigene exon track and make a list of it
                newgff = list( deepcopy( unigene_gff_list[0] ) )
                # update the coordinates
                newgff[2] = GFF_UG5UTREXON_FMETHOD
                newgff[3] = min(track)
                newgff[4] = max(track)
                # and insert as the first the new return unigene gff list
                return_unigene_gff_list.insert( 0, tuple(newgff) )
    
        # make UTR3UGExon track
        if ignore_3p_coords:
            ignore_3p_coords.sort()
            tracks  = [ [ ignore_3p_coords[0] ] ]
            for coord in ignore_3p_coords[1:]:
                if coord == max(tracks[-1])+1:
                    tracks[-1].append(coord)
                else:
                    tracks.append( [ coord ] )
            for track in tracks:
                # make a deepcopy of the first unigene exon track and make a list of it
                newgff = list( deepcopy( unigene_gff_list[0] ) )
                # update the coordinates
                newgff[2] = GFF_UG3UTREXON_FMETHOD
                newgff[3] = min(track)
                newgff[4] = max(track)
                # and append to the new return unigene gff list
                return_unigene_gff_list.append( tuple(newgff) )

    else:
        # hmm... not really expected. There are UniGene tracks,
        # but no UniGene exons are recognized. Probably a wrong setting
        # applied for GFF_UGEXON_FMETHOD (not identical to the naming in
        # the input gff.
        pass

    # order the unigene gff list (stop codon potentially on the front
    return_unigene_gff_list = order_gff_list(return_unigene_gff_list)
    ################################################################
    if verbose and (ignore_5p_coords or ignore_3p_coords):
        for track in return_unigene_gff_list: print track
    ################################################################

    # done! return the new list
    return return_unigene_gff_list, typeofunigene

# end of function correct_unigene_for_utrs


def get_annotated_genes_start_codon(gfflist):
    """
    Find the start codon in a list of gff tuples of a gene's features

    @type  gfflist: list
    @param gfflist: list with gff tuples of a gene structure (exons, introns, etc.)

    @rtype:  tuple
    @return: tuple representing the (annotated) gene's start codon

    @attention: Global variables GFF_GENESTART_FMETHOD and GFF_CDS_FMETHOD are required for this function
    """
    if GFF_GENESTART_FMETHOD in [ gff[2] for gff in gfflist ]:
        for gff in gfflist:
            if gff[2] == GFF_GENESTART_FMETHOD:
                return gff
    else:
        # no start codon track; take the first codon of the CDS
        for gff in gfflist:
            if gff[2] == GFF_CDS_FMETHOD:
                track = list( deepcopy(gff) )
                track[2] = GFF_GENESTART_FMETHOD
                track[4] = track[3]+2
                track[7] = 0
                return tuple( track )
        else:
            # Not expected; no cds track present at all
            # Just return empty tuple; putative erroneous later!
            return ()

# end of function get_annotated_genes_start_codon


def get_annotated_genes_stop_codon(gfflist):
    """
    Find the start codon in a list of gff tuples of a gene's features

    @type  gfflist: list
    @param gfflist: list with gff tuples of a gene structure (exons, introns, etc.)

    @rtype:  tuple
    @return: tuple representing the (annotated) gene's stop codon

    @attention: Global variables GFF_GENESTOP_FMETHOD and GFF_CDS_FMETHOD are required for this function
    """
    if GFF_GENESTOP_FMETHOD in [ gff[2] for gff in gfflist ]:
        for gff in gfflist:
            if gff[2] == GFF_GENESTOP_FMETHOD:
                return gff
    else:
        # no start codon track; take the first codon of the CDS
        for gff in gfflist:
            if gff[2] == GFF_CDS_FMETHOD:
                track = list( deepcopy(gff) )
                track[2] = GFF_GENESTOP_FMETHOD
                track[3] = track[4]-2
                track[7] = 0
                return tuple( track )
        else:
            # Not expected; no cds track present at all
            # Just return empty tuple; putative erroneous later!
            return ()

# end of function get_annotated_genes_stop_codon


def unigeneconfirmation(input,verbose=False):
    """
    Confirm given unigene GFF structure with Orfs of the DNA sequence

    @type  input: dict
    @param input: input dictionary data structure

    @type  verbose: Boolean
    @param verbose: print discrepancies to STDOUT (True) or not (False, default)

    @rtype:  tuple 
    @return: tuple of ( input data structure, TRACKS_ARE_PROPERLY_MATCHED boolean )
    """
    input = rungetorf(input)
    input = parseinputgff(input)

    UNIGENE_TRACKS_ARE_PROPERLY_MATCHED = True

    for org,data in input.iteritems():
        # check the splice sites of the complete unigene
        # BEFORE we chunked the unigene on utr5p and utr3p
        if input[org]['gff-unigene']:
            status,warnings = confirmcanonicalsplicesites(input[org]['genomeseq'],
                input[org]['gff-unigene'],exon_fmethod=GFF_UGEXON_FMETHOD,verbose=False)
            # set the warnings to the warnings list
            input[org]['warnings'].extend(warnings)
            if not status and verbose:
                UNIGENE_TRACKS_ARE_PROPERLY_MATCHED = False
                print "# WARNING: Non-canonical splice site(s) in UniGene for organism/fref:", org, input[org]['FREF']
                for warn in warnings: print warn

    for org,data in input.iteritems():
        # (A) check if both gene and unigene is applied
        # in that case, correct the unigene alignment for
        # putative UTRs that can mess up the performance measurement
        if input[org]['gff-unigene'] and input[org]['gff-gene']:
            # get start and stop codon of current gene annotation
            known_atg_gff = get_annotated_genes_start_codon(input[org]['gff-gene'])
            known_tga_gff = get_annotated_genes_stop_codon(input[org]['gff-gene'])
            # correct the unigene alignment for potential utrs
            input[org]['gff-unigene'],typeofunigene = correct_unigene_for_utrs(
                        input[org]['gff-unigene'],
                        start_codon_gff=known_atg_gff,
                        stop_codon_gff=known_tga_gff,
                        dnaseqfname=input[org]['genomeseqfile'],
                        verbose=verbose
                        )
            # store unigene annotation info into input dict
            input[org]['unigene-annotation'] = typeofunigene

        # (B) (re)order the unigene gff lists
        input[org]['gff-unigene'] = order_gff_list(input[org]['gff-unigene'])

        # (C) find the genestructure_orfmodel of the unigene structure
        # in the correct_unigene_for_utrs() function, unigene UTRs are
        # created when applicable. These UTRs are not taken into
        # account when mapping the unigene structure on the ORFs
        input[org]['orfid-unigenestructure'] = []
        if input[org]['gff-unigene']:
            unigenetracks = []
            for gfftrack in input[org]['gff-unigene']:
                if gfftrack[2] == GFF_UNIGENE_FMETHOD:
                    unigenetracks.append(gfftrack)
            unigenetracks.sort()
            # get the orfs of the genestructure
            orfids,_UNIGENE_TRACKS_ARE_PROPERLY_MATCHED = genestructuregff2orfs(
                    unigenetracks,input[org]['orfs'],
                    typeofunigene=typeofunigene,
                    verbose=verbose
                    )
            input[org]['orfid-unigenestructure'] = orfids
            # check if the unigene was properly matched.
            # If not -> remove the unigene orfid structure
            if not _UNIGENE_TRACKS_ARE_PROPERLY_MATCHED:
                # set unigene orfid structure to empty because it is bogus
                input[org]['orfid-unigenestructure'] = []
                UNIGENE_TRACKS_ARE_PROPERLY_MATCHED = False
                # append warning to warnings list
                warn = UniGeneStructureIsNotMappableOnOrfsWarning("not mappable on ORFs")
                input[org]['warnings'].append( warn )

        # (D) check start & stop codon of the unigene structure
        if input[org]['gff-unigene']:
            # Check - prior to confirmation of ATG & TGA of unigene
            # if the annotated tracks do not exceed EOF sequence.
            # yes, this can happen (locus not far enough extended,
            # improper unigene assembly, etc...)
            if input[org]['unigene-annotation'] in ['utr3p','fragm']:
                # partial unigene, not available at start codon
                UNIGENE_HAS_PROPER_START_CODON = None
                pass
            elif unigenetracks[0][3] >= 1:
                UNIGENE_HAS_PROPER_START_CODON = genestructurehasproperstartcodon(
                    unigenetracks,
                    input[org]['genomeseq'],
                    GFF_UGEXON_FMETHOD,
                    verbose=verbose)
            else:
                UNIGENE_HAS_PROPER_START_CODON = False

            if input[org]['unigene-annotation'] in ['utr5p','fragm']:
                # partial unigene, not available at stop codon
                UNIGENE_HAS_PROPER_STOP_CODON = None
                pass
            elif unigenetracks[-1][4] <= len(input[org]['genomeseq']):
                UNIGENE_HAS_PROPER_STOP_CODON = genestructurehasproperstopcodon(
                    unigenetracks,
                    input[org]['genomeseq'],
                    GFF_UGEXON_FMETHOD,
                    verbose=verbose)
            else:
                UNIGENE_HAS_PROPER_STOP_CODON = False

            # create warning messages if needed
            if UNIGENE_HAS_PROPER_START_CODON == False:
                # append warning to warnings list
                warn = UniGeneStructureIsNotMappableOnOrfsWarning("UniGene lacks proper start codon")
                input[org]['warnings'].append( warn ) 
                UNIGENE_TRACKS_ARE_PROPERLY_MATCHED = False
            if UNIGENE_HAS_PROPER_STOP_CODON == False:
                # append warning to warnings list
                warn = UniGeneStructureIsNotMappableOnOrfsWarning("UniGene lacks proper stop codon")
                input[org]['warnings'].append( warn )
                UNIGENE_TRACKS_ARE_PROPERLY_MATCHED = False
            ####################################################################
            if verbose and (UNIGENE_HAS_PROPER_START_CODON == False or\
            UNIGENE_HAS_PROPER_STOP_CODON == False):
                for track in input[org]['gff-gene']:
                    print "GENE:", track[0:7], "DNA:",
                    print len(input[org]['genomeseq'])
            ####################################################################

    # print the negative outcome in verbose mode
    if not UNIGENE_TRACKS_ARE_PROPERLY_MATCHED and verbose:
        print UniGeneStructureIsNotMappableOnOrfsWarning()

    # return the incremented input dict and error status (True or False)
    return input, UNIGENE_TRACKS_ARE_PROPERLY_MATCHED

# end of function unigeneconfirmation


def geneandunigeneconfirmation(input,verbose=False):
    """
    Carrier function for both gene- and unigene-confirmation

    @type  input: dict
    @param input: input dictionary data structure

    @type  verbose: Boolean
    @param verbose: print discrepancies to STDOUT (True) or not (False, default)

    @rtype:  tuple 
    @return: tuple of ( input data structure, 
                        GENE_TRACKS_ARE_PROPERLY_MATCHED boolean,
                        UNIGENE_TRACKS_ARE_PROPERLY_MATCHED boolean
                        )
    """
    # confirm gene & unigene structure
    input,gene_status = geneconfirmation(input,verbose=verbose)
    input,unigene_status = unigeneconfirmation(input,verbose=verbose)

    # return results
    return input, gene_status, unigene_status

# end of function geneandunigeneconfirmation


if __name__ == "__main__":
    # AbgpGeneLocus and AbgpGeneLocusDirectory imports
    from lib_optparse import (
            abgpoptparser,
            genelocustatsoptions,
            validate_genelocustatsoptions
            )
    from abgpgenelocusdirectory import AbgpGeneLocusDirectory

    # construct the command line option parser
    parser = abgpoptparser()
    genelocustatsoptions(parser)
    # parse the command line & validate
    (OPTIONS, args) = parser.parse_args()
    validate_genelocustatsoptions(parser,OPTIONS)

    # change directory to the one specified by optparse option
    os.chdir(OPTIONS.DIRECTORY)

    try:
        locus = AbgpGeneLocusDirectory(OPTIONS.DIRECTORY)
        input = locus.toinputdict()
    except:
        # TEMPORARILY backwards-compatibility with old input_data_struct.txt file
        input = eval(open('input_data_struct.txt').read().strip())

    # do splice site confirmation for unigene
    for org in input.keys():
        if not input[org].has_key('gff-unigene'): continue
        print org
        status,warnings = confirmcanonicalsplicesites(input[org]['genomeseq'],
                input[org]['gff-unigene'],exon_fmethod=GFF_UGEXON_FMETHOD,verbose=False)
        if not status:
            print "# WARNING: Non-canonical splice site(s) in UniGene for organism/fref:", org, input[org]['FREF']
            for warn in warnings: print warn

    # do unigeneconfirmation
    input,unigene_status = unigeneconfirmation(input,verbose=True)


    # go back to OPTIONS.CURRENTWORKINGDIR
    os.chdir(OPTIONS.CURRENTWORKINGDIR)

    if OPTIONS.verbose:
        print "# EOF unigeneconfirmation; if nothing printed, no abnormalities observed"

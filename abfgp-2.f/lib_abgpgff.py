"""
Functions for creating GFF outcome of ABGP results
"""

# Module metadata
__authors__ = "Ate van der Burgt"
__license__ = "MIT"

# Python Imports
from sets import Set
from copy import deepcopy
import os
from os import system

# Abgp Imports
from lib_stopwatch import StopWatch
from lib_signalp import create_projected_signalp_for_informants
from abgp_etc import abgpsysexit
from lib_pcg2blocks import PCG2inwpCBGS
from lib_codingarray import (
    PCG2codingarray,
    PCG2similarityarray,
    )
from lib_rip import RIPindex
from lib_signalp import ProjectedSignalPSignalPeptide
import pacb
import graphPlus
import dna2prot
from pythonlibs.uniqueness import get_random_string_tag
from gff import (
    gff2db, dbcleanup, gffs2txt, getseqbyfeaturetuple, 
    filtergffs4fmethods, filtergffs4fmethod, mirrorgfflist,
    gffs2coordset, exons2introns, parsegfffile
    )
from gene.start import scan_pssm_tss, TranslationalStartSite, IC_TSS_PATTERN_OFFSET
from gene.splicesite import scan_pssm_splice_site
from gene.intron import IntronConnectingOrfs
from gene.intron.comparison import (
    _branchpoint_comparison,
    _polypyrimidinetract_comparison,
    _algsimilarity_comparison,
    _finetune_splicesite_comparison,
    )
from gene.sequenceerror import (
    SequenceErrorConnectingOrfs,
    SequenceErrorCoordinate,
    ExonStart,
    ExonEnd,
    )

from gene.exon import (
    ExonOnOrf,
    FirstExonOnOrf,
    FinalExonOnOrf,
    SingleExonOnOrf
    )
from gene.gene_exceptions import IncompatibleSpliceSitePhases
    
    

# Import Settings
from settings.abgp import *
from settings.gff import *
from settings.ggb import *
from settings.dbwarehouse import (
    ABGP_ORGANISM2FULLSPECIESNAME_MAPPING, 
    _get_organism_full_name,
    )


def _get_main_interface(labels):
    """ Simple majority voting of list elements """
    ordered = [ ( labels.count(label), label ) for label in Set(labels) ]
    ordered.sort()
    return ordered[-1][1]
# end of function _get_main_interface


def performance_comparison_of_sets( predicted, known ):
    """
    Do a performance comparison on two sets of coordinates

    @type  predicted: sets.Set
    @param predicted: Set with all predicted coordinates

    @type  known: sets.Set
    @param known: Set with all known/correct coordinates

    @rtype:  list of 3 lists of lists of coordinates
    @return: list of 3 list representing the correct (TP), over (FP) and under (FN) predicted coordinates, separated by gaps
    """
    ### Example:
    ### predicted = Set(   2,3,4,5,10,11,12,13,14,15,20,21,22,23,24,25,26,27,28 )
    ### known     = Set( 1,2,3,4,5,10,11,12,13,14,15,   21,22,23,24,25,26       )
    ### Return structure:
    ### [
    ###   [ [2,3,4,5,10,11,12,13,14,15], [21,22,23,24,25,26] ], # TP correct predicted
    ###   [ [20], [27,28] ],                                    # FP over predicted
    ###   [ [1] ],                                              # FN under predicted
    ### ]

    # make CORRECT, UNDER and OVER predicted list
    correct        = predicted.intersection( known )
    overpredicted  = predicted.difference( known )
    underpredicted = known.difference( predicted )

    returnlists    = []
    for item in ( correct, overpredicted, underpredicted ):
        if item:
            item = list(item)
            item.sort()
            tracks = [ [ item[0] ] ]
            for coord in item[1:]:
                if coord == max(tracks[-1])+1:
                    tracks[-1].append(coord)
                else:
                    tracks.append( [ coord ] )
            returnlists.append( tracks )
        else:
            # no overlap of this kind!
            returnlists.append( [] )

    # return the data structure
    return returnlists

# end of function performance_comparison_of_sets


def abgp_performance_measurement_of_gff(gff_predicted, gff_known, fmethod=[], gclass=[], gffdata={}):
    """
    Output gff tracks that measure the performance of the gff of abgp predicted gene structure

    @type  fmethod: list
    @param fmethod: list of strings of fmethods (column 3) to take into account

    @type  gclass: list
    @param gclass: list of strings of glcass's (column 9 first item) to take into account

    @type  gffdata: dictionary
    @param gffdata: fictionairy with at least the exact keys 'fref' and 'gname'

    @attention: gffdata dictionary MUST contain keys 'fref' and 'gname'
    @attention: Global variables GFF_PERFORMANCE_* (20 distinct variables) must be set
    """

    predicted = gffs2coordset(gff_predicted,fmethod=fmethod,gclass=glcass)
    known     = gffs2coordset(gff_known,fmethod=fmethod,gclass=glcass)
    return abgp_performance_measurement_of_sets(predicted,known,gffdata=gffdata)

# end of function abgp_performance_measurement_of_gff


def abgp_performance_measurement_of_sets(predicted,known,gffdata={}):
    """
    Output gff tracks that measure the performance of the gff of abgp predicted gene structure

    @type  predicted: sets.Set
    @param predicted: Set with all predicted coordinates

    @type  known: sets.Set
    @param known: Set with all known/correct coordinates

    @type  gffdata: dictionary
    @param gffdata: dictionary with at least the exact keys 'fref' and 'gname'

    @attention: gffdata dictionary MUST contain keys 'fref' and 'gname'
    @attention: Global variables GFF_PERFORMANCE_* (20 distinct variables) must be set
    """
    # return list of gff tuples
    retgff = []
    fref   = gffdata['fref']
    gname  = gffdata['gname']

    # do the performance comparison of the sets
    (correct,over,under) = performance_comparison_of_sets(predicted,known)

    # make True Positives (TP) CORRECT predicted track
    for track in correct:
        gff = ( fref, GFF_PERFORMANCE_TP_FSOURCE, GFF_PERFORMANCE_TP_FMETHOD,
                min(track), max(track), '.', '+', '.', "%s %s" % (GFF_PERFORMANCE_TP_GCLASS, gname)
                )
        retgff.append(gff)

    # make False Positives (FP) OVER predicted track
    for track in over:
        gff = ( fref, GFF_PERFORMANCE_FP_FSOURCE, GFF_PERFORMANCE_FP_FMETHOD,
                min(track), max(track), '.', '+', '.', "%s %s" % (GFF_PERFORMANCE_FP_GCLASS, gname)
                )
        retgff.append(gff)

    # make False Negatives (FN) UNDER predicted track
    for track in under:
        gff = ( fref, GFF_PERFORMANCE_FN_FSOURCE, GFF_PERFORMANCE_FN_FMETHOD,
                min(track), max(track), '.', '+', '.', "%s %s" % (GFF_PERFORMANCE_FN_GCLASS, gname)
                )
        retgff.append(gff)

    # and return the gff
    return retgff

# end of function abgp_performance_measurement_of_sets


def _tcode2gff(tcodedata,report_step_size=1,FREF=None):
    """
    Convert list to tcode output tuples to gff quantitative data track
    
    @type  tcodedata: sets.Set
    @param tcodedata: Set with all predicted coordinates

    @type  report_step_size: positive integer
    @param report_step_size: defines the granularity of the track (1==each nt)

    @rtype:  list of gff tuples
    @return: list with gff tuples with tcode data

    @attention: Global variables GFF_TCODE_FMETHOD,GFF_TCODE_FSOURCE,GFF_TCODE_GCLASS and FREF must be available
    """
    gffdata = []
    for pos in range(0,len(tcodedata),report_step_size):
        item = tcodedata[pos]
        if item[3] == 'No opinion': continue
        center = ( item[0]-1 + item[1] ) / 2
        if item[3] == 'Coding':
            fmethod = GFF_TCODE_FMETHOD + '_coding'
        else:
            fmethod = GFF_TCODE_FMETHOD + '_noncoding'
        gff = ( FREF, GFF_TCODE_FSOURCE, fmethod, center, center, item[2], '+', '.', "%s %s; Note '%s'" % (
            GFF_TCODE_GCLASS, FREF, item[3] ) )
        gffdata.append( gff )
    # return gffdata with tcode gff tracks
    return gffdata

# end of function _tcode2gff


def _tcodeprofile2gff(tcodedata,report_step_size=9,FREF=None,tcode_5p_windowsize=201,tcode_3p_windowsize=201):
    """
    """
    gffdata = []
    fmethod = GFF_TCODE_FMETHOD + '_profile'
    for pos in range(tcode_5p_windowsize,len(tcodedata)-tcode_3p_windowsize,report_step_size):
        (start,stop,score,status) = tcodedata[pos]
        list5p = [ tup[2] for tup in tcodedata[pos-tcode_5p_windowsize:pos] ]
        list3p = [ tup[2] for tup in tcodedata[pos:pos+tcode_3p_windowsize] ]
        while None in list5p: list5p.remove(None)
        while None in list3p: list3p.remove(None)
        if len(list5p) == 0 or len(list3p) == 0:
            score = 1.0
        else:
            score = ( sum(list5p) / float(len(list5p)) ) / ( sum(list3p) / float(len(list3p)) )
        center = ( start-1 + stop ) / 2
        gff = ( FREF, GFF_TCODE_FSOURCE, fmethod, center, center, score, '+', '.', "%s %s" % (
            GFF_TCODE_GCLASS, FREF ) )
        gffdata.append( gff )
    # return gffdata with tcode gff tracks
    return gffdata

# end of function _tcodeprofile2gff


def _get_abgp_file_basename(OPTIONS):
    """
    Get basename for for the output files (based on target/informants used) 

    @type  OPTIONS: optparse options object
    @param OPTIONS: optparse options object

    @rtype:  string
    @return: basename for the output files
    """
    if OPTIONS.target:
        try:
            num_loci = OPTIONS.selected_num_loci
        except:
            num_loci = len(OPTIONS.loci) + len(OPTIONS.dnafiles)
        return "%s.%s%sSL." % (ABFGP_VERSION,OPTIONS.target,num_loci)
    else:
        return "%s." % ABGP_VERSION

# end of function _get_abgp_file_basename


def _create_unique_filename_with_integer_suffix(fullpath):
    """
    Create an unique and not yet existing filename for this output file

    @type  fullpath: string
    @param fullpath: (absolute) path to arbitrary output file

    @rtype:  string
    @return: unique and full file path for this output file
    """
    # create an unique filename
    suffix = None
    suffix_cnt=1
    while os.path.exists(fullpath):
        if suffix: fullpath = fullpath[0:-len(suffix)]
        suffix = ".%s" % suffix_cnt
        suffix_cnt+=1
        fullpath = fullpath + suffix
    return fullpath

# end of function _create_unique_filename_with_integer_suffix


def _get_abs_fref(organism,input,OPTIONS):
    """
    Get fref string for the absolute GFF output

    @type  organism: string
    @param organism: Organism / Gene Identifier

    @type  input: dict
    @param input: input <dict data structure> with data from input files

    @type  OPTIONS: optparse options object
    @param OPTIONS: optparse options object

    @rtype:  string
    @return: fref gff identifier

    @attention: only truely valid for AbgpGeneLocusDirectory(s) as input
    """
    if input[organism]['locusfile']:
        # fref == parsed from locus file!
        fref = input[organism]['genomefref']
    else:
        # No locusfile, but a DNA file.
        # This means absolute GFF cannot be created!
        # But, that must be tackled somewhere else, not here!
        fref = input[organism]['genomefref']

    # return the absolute fref from the input file/AbgpGeneLocusDirectory
    return fref

# end of function _get_abs_fref


def _get_rel_fref(organism,input,OPTIONS):
    """
    Get fref string for the relative GFF output

    @type  organism: string
    @param organism: Organism / Gene Identifier

    @type  input: dict
    @param input: input <dict data structure> with data from input files

    @type  OPTIONS: optparse options object
    @param OPTIONS: optparse options object

    @rtype:  string
    @return: fref gff identifier

    @attention: informant of --target ABGP_<informant>_for_<target>_<geneloci>SL
    @attention: target identifier     ABGP_<target>_<geneloci>SL
    @attention: no --target           ABGP_<proteinfref>_<geneloci>AL
    """
    protfref = input[organism]['proteinfref']
    if not protfref: protfref = organism

    fref = None
    if OPTIONS.target:
        if organism == OPTIONS.target:
            fref = "%s_%s_%sSL" % ( ABGP_PROGRAM_NAME, protfref, len(input) )
        else:
            prottargetfref = input[OPTIONS.target]['proteinfref']
            if not prottargetfref: prottargetfref = OPTIONS.target
            fref = "%s_%s_for_%s_%sSL" % ( ABGP_PROGRAM_NAME, protfref,
                   prottargetfref, len(input) )
    else:
        fref = "%s_%s_%sAL" % ( ABFGP_VERSION, input[organism]['proteinfref'] , len(input) )

    # return the constructed relative fref
    return fref

# end of function _get_rel_fref


def create_detailed_gff_file(GSofCBG,OPTIONS,organism=None,verbose=False):
    """
    @type  GSofCBG: GenestructureOfCodingBlockGraphs
    @param GSofCBG: GenestructureOfCodingBlockGraphs instance to create gff for

    @type  OPTIONS: optparse options object
    @param OPTIONS: optparse options object

    @type  organism: *
    @param organism: Organism identifier (presumably string)

    @type  verbose: Boolean 
    @param verbose: print (timing) messages to STDOUT (True) or not (False) 

    @rtype:  tuple
    @return: tuple of 2 created file names with relative paths ( gff_file_name, fasta_file_name )

    @attention: Dozens/hundreds of Global variables needed (imported from settings file)
    """
    if not organism: raise "organism not specified"

    stw = StopWatch(name="CreateGffStopWatch")
    stw.start()

    # get old thingies as pointers from object ...
    subgraphs = GSofCBG.codingblockgraphs
    input     = GSofCBG.input

    # get all gff lines that are already gather for this organism
    gffs   = input[organism]['gffs']
    METHOD = input[organism]['METHOD']
    FREF   = _get_rel_fref(organism,GSofCBG.input,OPTIONS)

    # get gff lines for each ORF
    for orf in input[organism]['orfs'].orfs: gffs.append( orf.togff(fref=FREF) )

    ####################################################################
    if verbose: print stw.lap(), "preparation"
    ####################################################################

    # get gff lines for tcode tracks
    if USE_TCODE and GFF_TCODE_OUTPUT:
        # report not all steps, but only a subset of them
        report_step_size = EXECUTABLE_TCODE_GFFOUTPUT_STEP / EXECUTABLE_TCODE_STEP
        gffs.extend( _tcode2gff( input[organism]['tcode'], report_step_size=report_step_size, FREF=FREF ) )
        profile_report_step_size = report_step_size*2
        gffs.extend( _tcodeprofile2gff( input[organism]['tcode'],
                tcode_5p_windowsize=301,
                tcode_3p_windowsize=301,
                report_step_size=profile_report_step_size,
                FREF=FREF ) )

    ####################################################################
    if verbose: print stw.lap(), "tcode"
    ####################################################################

    # scan the input sequence first ~1700nt for high scoring TSS
    endscanpos = 1600
    # obtain endscanpos from aligned startcodon graph
    firstcbg   = subgraphs[0]
    if firstcbg:
        endscanpos = min( firstcbg.overall_minimal_spanning_range(organism=organism) )*3 + 150
    else:
        pass
    if endscanpos < 1600: endscanpos = 1600
    tsslist = scan_pssm_tss(input[organism]['genomeseq'][0:endscanpos],min_pssm_score=0.0)
    for tss in tsslist:
        _gff = {'fref':FREF, 'fsource': 'tssPSSM4FUNGI', 'fmethod': 'tsspssm' }
        # format tss.pssm_score as a string-formatted float
        tss.pssm_score = "%2.1f" % tss.pssm_score
        gffs.append( tss.togff( gff=_gff ) )

    ####################################################################
    if verbose: print stw.lap(), "tssscan"
    ####################################################################

    # dictionaries for projected stop codons of a pacbp
    projected_leading_stops = {}
    projected_tailing_stops = {}

    # get gff lines for each similarity/pacbp
    if GFF_ORFSIMILARITY_OUTPUT:
        for sg in subgraphs:
            if organism not in sg.organism_set(): continue
            theorf = sg.get_orfs_of_graph(organism=organism)[0]
            for (k,(org1,orf1),(org2,orf2)), pacbporf in sg.pacbps.iteritems():
                if organism == org1:
                    thepacbporf = pacbporf
                    group       = "%s-%s-orfid-%s" % (FREF,org2,orf2)
                    # get the TargetFref of the other organism.
                    # This c9tv is used for interlinking in the GGB!
                    targetfref = _get_rel_fref(org2,GSofCBG.input,OPTIONS)
                elif organism == org2:
                    # swap query and sbjct
                    thepacbporf = pacb.swap_query_and_sbjct( pacbporf )
                    group       = "%s-%s-orfid-%s" % (FREF,org1,orf1)
                    # get the TargetFref of the other organism
                    # This c9tv is used for interlinking in the GGB!
                    targetfref = _get_rel_fref(org1,GSofCBG.input,OPTIONS)
                else:
                    continue
    
                # get length of pacbporf
                length = len(thepacbporf ._positions)
                if hasattr(thepacbporf ,'_original_alignment_pos_start'):
                    length = thepacbporf ._original_alignment_pos_end -\
                             thepacbporf ._original_alignment_pos_start
    
                # gff dict for this feature          
                gff = { 'fref'        : FREF, 
                        'fsource'     : GFF_ORFSIMILARITY_FSOURCE,
                        'fmethod'     : GFF_ORFSIMILARITY_FMETHOD,
                        'gclass'      : GFF_ORFSIMILARITY_GCLASS,
                        'column9data' : {
                            'Note'      : "Similarity on ORF %s with ORF %s" % (thepacbporf.orfQ.id,thepacbporf.orfS.id),
                            'Length'    : length,
                            'Bitscore'  : thepacbporf.bitscore,
                            'Similarity': thepacbporf.similarity,
                            'Identity'  : thepacbporf.identity,
                            'TargetFref': targetfref,
                            },
                        'gname' : group,
                    }
                # and make the GFF line
                gffs.append( thepacbporf.togff(gff=gff,with_alignment=False) )

                if GFF_ORFSIM_PLS_OUTPUT or GFF_ORFSIM_PTS_OUTPUT:
                    # fill gff dictionary from setting
                    _gff_pls = { 'fref'    : FREF,
                                 'fsource' : GFF_ORFSIM_PLS_FSOURCE,
                                 'fmethod' : GFF_ORFSIM_PLS_FMETHOD,
                                 'column9data': {'TargetFref': targetfref },
                                 'gclass'  : GFF_ORFSIM_PLS_GCLASS, }
                    _gff_pts = { 'fref'    : FREF,
                                 'fsource' : GFF_ORFSIM_PTS_FSOURCE,
                                 'fmethod' : GFF_ORFSIM_PTS_FMETHOD,
                                 'column9data': {'TargetFref': targetfref },
                                 'gclass'  : GFF_ORFSIM_PTS_GCLASS, }
                    # check if dict key `group` exists already
                    if not projected_leading_stops.has_key(group):
                        projected_leading_stops[group] = []
                        projected_tailing_stops[group] = []
        
                    # make gffs of projected ORF stop sites
                    projected_leading_stops[group].append( thepacbporf.gff_projected_leading_stop(gff=_gff_pls) )
                    projected_tailing_stops[group].append( thepacbporf.gff_projected_tailing_stop(gff=_gff_pts) )

            # done with this CBG
            if verbose: print stw.lap(), "orfsim", sg

    ####################################################################
    if verbose: print stw.lap(), "similarity"
    ####################################################################

    if GFF_ORFSIM_PLS_OUTPUT:
        # take the projection of the most LEFT cbg
        for group, gfflist in projected_leading_stops.iteritems():
            gffs.append(gfflist[0]) 
    if GFF_ORFSIM_PTS_OUTPUT:
        # take the projection of the most RIGTH cbg
        for group, gfflist in projected_tailing_stops.iteritems():
            gffs.append(gfflist[-1])

    ####################################################################
    if verbose: print stw.lap(), "projected * codons"
    ####################################################################

    ## make gff lines for (aligned) stop sites
    #for sg in subgraphs:
    #    if organism not in sg.organism_set(): continue
    #    if not sg._stopcodongraph: continue
    #    gffs.append( sg._stopcodongraph.togff( gff={'fref': FREF }, organism=organism ) )


    for sg in subgraphs[1:]:
        if organism not in sg.organism_set(): continue
        gffs.extend( sg._CBGinterface5p.togff(organism,FREF=FREF) )

    ####################################################################
    if verbose: print stw.lap(), "cbgIFs"
    ####################################################################

    # make gff for the abgp CodingBlocks of the genestructure
    if GFF_CODINGBLOCK_OUTPUT:
        _cb_gff = { 'fref'   : FREF,
                    'fmethod': GFF_CODINGBLOCK_FMETHOD,
                    'fsource': GFF_CODINGBLOCK_FSOURCE,
                    'gclass' : GFF_CODINGBLOCK_GCLASS }
        cbggffdata = GSofCBG.get_codingblock_gff(organism=organism,gff=_cb_gff)
        gffs.extend( cbggffdata )

    ####################################################################
    if verbose: print stw.lap(), "CodingBlockGraphs"
    ####################################################################

    # get multiple alignment data
    gffs.extend( get_gff_of_gsg2multiplealignment(GSofCBG,organism,gff={'fref':FREF}) )

    # get cexpander data
    gffs.extend( get_gff_of_gsg2cexpander(GSofCBG,organism,gff={'fref':FREF}) )

    # coordinate Set of the predicted gene structure
    abgp_genestructure_coordinate_set = Set()

    if GFF_ABFGP_EXON_OUTPUT:
        _outcome_gff = {'fref'   : FREF,
                        'fmethod': GFF_ABFGP_EXON_FMETHOD,
                        'fsource': GFF_ABFGP_EXON_FSOURCE,
                        'gclass' : GFF_ABFGP_EXON_GCLASS }
        outcome = GSofCBG.get_genestructure_gff(organism=organism,gff=_outcome_gff)
        # make intron list of exons
        introns = exons2introns(outcome,fmethod=GFF_ABFGP_INTRON_FMETHOD,
                                        gclass=GFF_ABFGP_INTRON_GCLASS)
        # make the genestructure coordinate set for benchmarking
        abgp_genestructure_coordinate_set = gffs2coordset(
                    outcome,
                    fmethod=[GFF_ABFGP_EXON_FMETHOD]
                    )
        # and extend to gffs
        gffs.extend( outcome )

    ####################################################################
    if verbose: print stw.lap(), "abgp genestructure"
    ####################################################################


    # analyses of the (aligned) stop codons
    finalCBG = GSofCBG.get_final_cbg()
    allperf  = finalCBG._stopcodongraph.is_optimal()
    orgperf  = finalCBG._stopcodongraph.is_optimal(organism=organism)

    # map orgperf/allperf to a binary score
    stopcodongffscoremapper = {
        (True,True):   2,
        (True,False):  1,
        (None,False):  -1, # future deprecated !?
        (False,False): -1,
        (False,True):  -1,
    }
    stopcodongffntsize = 15
    start = max(abgp_genestructure_coordinate_set)+1
    stop  = start + stopcodongffntsize -1

    gfftrack = ( FREF, GFF_STOPCODON_OPTIMALITY_FSOURCE,
        GFF_STOPCODON_OPTIMALITY_FMETHOD, start,stop,
        stopcodongffscoremapper[(orgperf,allperf)],'+','.',
        "%s %s" % (GFF_STOPCODON_OPTIMALITY_GCLASS,FREF)
        )
    gffs.append( gfftrack )

    # analyses of the aligned TSS / start codons
    firstCBG = GSofCBG.get_first_cbg()
    if firstCBG._startcodongraph and firstCBG._startcodongraph.alignedsites:
        besttssgraph = firstCBG._startcodongraph.alignedsites[0]
        if besttssgraph.__class__.__name__ == 'TranslationalStartSiteCollectionGraph':
            orfperf = False
            allperf = False
        else:
            besttssgraph._codingblockgraph = firstCBG
            orfperf = besttssgraph.is_optimal(organism=organism)
            allperf = besttssgraph.is_optimal()
    else:
        orfperf = False
        allperf = False

    # map orgperf/allperf to a binary score
    startcodongffscoremapper = {
        (True,True):   2,
        (True,False):  1,
        (None,False):  0, # future deprecated !?
        (None,True):   0, # future deprecated !?
        (None,None):   0, # future deprecated !?
        (False,False): -1,
        (False,True):  -1,
    }

    startcodongffntsize = 15
    stop  = min(abgp_genestructure_coordinate_set)
    start = stop - startcodongffntsize + 1

    gfftrack = ( FREF, GFF_STARTCODON_OPTIMALITY_FSOURCE, 
        GFF_STARTCODON_OPTIMALITY_FMETHOD, start,stop,
        startcodongffscoremapper[(orgperf,allperf)],'-','.',
        "%s %s" % (GFF_STARTCODON_OPTIMALITY_GCLASS,FREF) 
        )
    gffs.append( gfftrack )


    ####################################################################
    if verbose: print stw.lap(), "performance start & stop codon"
    ####################################################################

    # check if gene structure tracks must be outputted
    if GFF_GENE_OUTPUT:     gffs.extend( input[organism]['gff-gene'] )

    # check if unigene structure tracks must be outputted
    if GFF_UNIGENE_OUTPUT:  gffs.extend( input[organism]['gff-unigene'] )

    # check if performance output tracks must be outputted
    # and take into account if an unigene was applied
    if GFF_PERFORMANCE_OUTPUT and input[organism]['gff-unigene']:
        # make unigene_coordinate_set
        unigene_coordinate_set = gffs2coordset(
                    input[organism]['gff-unigene'],
                    fmethod=[GFF_UGEXON_FMETHOD]
                    )

        # do performance comparison
        performance_gff_unigene = abgp_performance_measurement_of_sets(
                abgp_genestructure_coordinate_set,
                unigene_coordinate_set,
                gffdata = {
                    'gname':   '%s-unigene' % organism,
                    'fref':     FREF,  }  )

        # extend to the gffs list
        gffs.extend( performance_gff_unigene )
        # some printing
        for gff in performance_gff_unigene: print gff

    ####################################################################
    if verbose: print stw.lap(), "performance unigene"
    ####################################################################

    # check if performance output tracks must be outputted
    # and take into account if an annotated gene structure was applied
    if GFF_PERFORMANCE_OUTPUT and input[organism]['gff-gene']:
        # make gene_coordinate_set
        gene_coordinate_set = gffs2coordset(
                    input[organism]['gff-gene'],
                    fmethod=[GFF_CDS_FMETHOD]
                    )

        # correct for STOP codon that is NOT included in the gene structure
        # TODO: make this formal; ik can be that there are gene structures
        # given that DO include a stop codon....
        stop_codon_pos = max(gene_coordinate_set)
        gene_coordinate_set.update(range(stop_codon_pos+1,stop_codon_pos+4))

        # do performance comparison
        performance_gff_gene = abgp_performance_measurement_of_sets(
                abgp_genestructure_coordinate_set,
                gene_coordinate_set,
                gffdata = {
                    'gname':   '%s-gene' % organism,
                    'fref':     FREF,  }  )

        # extend to the gffs list
        gffs.extend( performance_gff_gene )

    ####################################################################
    if verbose: print stw.lap(), "performance unigene"
    ####################################################################

    # check if orf performance output tracks must be outputted
    # and take into account if an annotated gene structure was applied
    if GFF_ORF_PERFORMANCE_OUTPUT and input[organism]['gff-gene']:
        orfperformance = []
        gene_orfids = Set(input[organism]['orfid-genestructure'])
        abgp_orfids = Set(GSofCBG.create_genestructure_orfmodel()[organism])
        # get rid of dashes in the genestructure (`missing` orfs in CBGs of this organism)
        if '-' in abgp_orfids: abgp_orfids.remove('-')
        for orfid_TP in gene_orfids.intersection(abgp_orfids):
            theorf = input[organism]['orfs'].get_orf_by_id(orfid_TP)
            orfperformance.append( theorf.togff(
                    fref=FREF,
                    fmethod=GFF_ORF_PERFORMANCE_TP_FMETHOD,
                    fsource=GFF_ORF_PERFORMANCE_TP_FSOURCE,
                    gclass= GFF_ORF_PERFORMANCE_TP_GCLASS,
                    one_group_per_frame=False
                    ) )
        for orfid_FN in gene_orfids.difference(abgp_orfids):
            theorf = input[organism]['orfs'].get_orf_by_id(orfid_FN)
            orfperformance.append( theorf.togff(
                    fref=FREF,
                    fmethod=GFF_ORF_PERFORMANCE_FN_FMETHOD,
                    fsource=GFF_ORF_PERFORMANCE_FN_FSOURCE,
                    gclass= GFF_ORF_PERFORMANCE_FN_GCLASS,
                    one_group_per_frame=False
                    ) )
        for orfid_FP in abgp_orfids.difference(gene_orfids):
            theorf = input[organism]['orfs'].get_orf_by_id(orfid_FP)
            orfperformance.append( theorf.togff(
                    fref=FREF,
                    fmethod=GFF_ORF_PERFORMANCE_FP_FMETHOD,
                    fsource=GFF_ORF_PERFORMANCE_FP_FSOURCE,
                    gclass= GFF_ORF_PERFORMANCE_FP_GCLASS,
                    one_group_per_frame=False
                    ) )

        # extend the orfperformance to main gffs
        gffs.extend( orfperformance )

    ####################################################################
    if verbose: print stw.lap(), "performance orfs"
    ####################################################################

    if hasattr(GSofCBG,'_reversecomplementcbgs'):
        for revCBG in GSofCBG._reversecomplementcbgs:
            gffs.extend(revCBG.togff(organism,gff={'fref':FREF}))

    ####################################################################
    if verbose: print stw.lap(), "metadata tracks"
    ####################################################################

    # TODO: replace FREF in known gff lines...
    # TODO: only in annotated gene structure & unigene alignment, a wrong `fref` is used.
    while [] in gffs: gffs.remove([])
    while () in gffs: gffs.remove(())
    for i in range(0,len(gffs)):
        try:
            line = list(gffs[i])
            line[0] = FREF
            gffs[i] = tuple(line)
        except:
            print "# WARNING: unexpected gff line: ", line

    # add the SEQUENCE track to the gffs
    seqlen = len(input[organism]['genomeseq'])
    gff = ( FREF, METHOD, 'Sequence', 1, seqlen, '.', '+', '.', 'Sequence %s ; Organism %s' % (FREF, organism) )
    gffs.insert( 0, gff)

    ####################################################################
    if verbose: print stw.lap(), "finalisation"
    ####################################################################

    # Ready with gathering gff data! Create file
    fname  = "%sfullresult.%s.gff" % ( _get_abgp_file_basename(OPTIONS), organism )
    # make fullpath of fname
    fullpath = os.path.join(OPTIONS.outdir,fname)
    # check if file exists and overwriting is aloued
    if not OPTIONS.force:
        fullpath = _create_unique_filename_with_integer_suffix(fullpath)

    # write to GFF file
    fh = open(fullpath, 'w')
    fh.write(gffs2txt(gffs))
    fh.close()

    ####################################################################
    if verbose: print stw.lap(), "fullresult gff file written->done!"
    ####################################################################

    # done! return the freshly created gff file name
    return fullpath 

# end of function create_detailed_gff_file


def create_abgp_prediction_gff_files(GSofCBGs,OPTIONS,organism=None,verbose=False):
    """ """
    performance_data_dict = {}
    created_files = {}
    for thisorg in GSofCBGs.organism_set():
        if organism and thisorg != organism: continue

        gffdata = []
        fref_abs  = _get_abs_fref(thisorg,GSofCBGs.input,OPTIONS)
        fref_rel  = _get_rel_fref(thisorg,GSofCBGs.input,OPTIONS)
        basefname = _get_abgp_file_basename(OPTIONS)

        # extend the currently annotated gene structure    
        gffdata.extend( GSofCBGs.input[thisorg]['gff-gene'] )
    
        # extend the unigene structure tracks
        gffdata.extend( GSofCBGs.input[thisorg]['gff-unigene'] )

        # create some defaults for gff tracks
        METHOD = GSofCBGs.input[thisorg]['METHOD']
        _outcome_gff = {'fref'      : fref_rel,
                        'fmethod'   : GFF_ABFGP_EXON_FMETHOD,
                        'fsource'   : GFF_ABFGP_EXON_FSOURCE,
                        'gclass'    : GFF_ABFGP_EXON_GCLASS }
    
        # obtain the outcome from the GSG
        outcome = GSofCBGs.get_genestructure_gff(
                    organism=thisorg,gff=_outcome_gff)
        # make intron list of exons
        introns = exons2introns(outcome,fmethod=GFF_ABFGP_INTRON_FMETHOD,
                                        gclass=GFF_ABFGP_INTRON_GCLASS)

        # make the genestructure coordinate set for benchmarking
        abgp_genestructure_coordinate_set = gffs2coordset(
                outcome,
                fmethod=[GFF_ABFGP_EXON_FMETHOD]
                ) 
        # and extend to gffs
        gffdata.extend( outcome )
        gffdata.extend( introns )

        ################################################################
        #### Performance analyses of prediction towards unigene     ####
        ################################################################
        if GSofCBGs.input[thisorg]['gff-unigene']:
            # make unigene_coordinate_set
            unigene_coordinate_set = gffs2coordset(
                    GSofCBGs.input[thisorg]['gff-unigene'],
                    fmethod=[GFF_UGEXON_FMETHOD]
                    )
            # do performance comparison
            performance_gff_unigene = abgp_performance_measurement_of_sets(
                    abgp_genestructure_coordinate_set,
                    unigene_coordinate_set,
                    gffdata = {
                        'gname':   '%s-unigene' % thisorg,
                        'fref':     fref_rel,  }  )
            # extend to the gffs list
            gffdata.extend( performance_gff_unigene )
            # update performance_data_dict
            outcome = list(Set([ line[2] for line in performance_gff_unigene ]))
            if not performance_data_dict.has_key(thisorg):
                performance_data_dict[thisorg] = {}
            if outcome == [ GFF_PERFORMANCE_TP_FMETHOD ]:
                performance_data_dict[thisorg]['unigene'] = True
            else:
                performance_data_dict[thisorg]['unigene'] = False


        ################################################################
        #### Performance analyses of prediction towards genestructure###
        ################################################################
        if GSofCBGs.input[thisorg]['gff-gene']:
            # make gene_coordinate_set
            gene_coordinate_set = gffs2coordset(
                    GSofCBGs.input[thisorg]['gff-gene'],
                    fmethod=[GFF_CDS_FMETHOD]
                    )
            # correct for STOP codon that is NOT included in the gene structure
            # TODO: make this formal; ik can be that there are gene structures
            # given that DO include a stop codon....
            stop_codon_pos = max(gene_coordinate_set)
            gene_coordinate_set.update(range(stop_codon_pos+1,stop_codon_pos+4))
            # do performance comparison
            performance_gff_gene = abgp_performance_measurement_of_sets(
                    abgp_genestructure_coordinate_set,
                    gene_coordinate_set,
                    gffdata = {
                        'gname':   '%s-gene' % thisorg,
                        'fref':     fref_rel,  }  )
            # extend to the gffs list
            gffdata.extend( performance_gff_gene )
            # update performance_data_dict
            outcome = list(Set([ line[2] for line in performance_gff_gene ]))
            if not performance_data_dict.has_key(thisorg):
                performance_data_dict[thisorg] = {}
            if outcome == [ GFF_PERFORMANCE_TP_FMETHOD ]:
                performance_data_dict[thisorg]['gene'] = True
            else:
                performance_data_dict[thisorg]['gene'] = False


        # add the SEQUENCE track as first item to the gffs
        seqlen = len(GSofCBGs.input[thisorg]['genomeseq'])
        col9   = 'Sequence %s ; Organism %s' % ( fref_rel, thisorg)
        gfftrack = ( fref_rel, METHOD, 'Sequence', 1, seqlen, '.', '+', '.', col9 )
        gffdata.insert( 0, gfftrack )

        ################################################################
        #### Similarity (PacbPORF) data                             ####
        ################################################################
        projected_leading_stops = {}
        projected_tailing_stops = {}
        for cbg in GSofCBGs.codingblockgraphs:
            if thisorg not in cbg.organism_set(): continue
            theorf = cbg.get_orfs_of_graph(organism=thisorg)[0]
            for (k,(org1,orf1),(org2,orf2)), pacbporf in cbg.pacbps.iteritems():
                if thisorg == org1:
                    thepacbporf = pacbporf
                    group       = "%s-%s-orfid-%s" % (fref_rel,org2,orf2)
                    # get the TargetFref of the other organism.
                    # This c9tv is used for interlinking in the GGB!
                    targetfref = _get_rel_fref(org2,GSofCBGs.input,OPTIONS)
                elif thisorg == org2:
                    # swap query and sbjct
                    thepacbporf = pacb.swap_query_and_sbjct( pacbporf )
                    group       = "%s-%s-orfid-%s" % (fref_rel,org1,orf1)
                    # get the TargetFref of the other organism
                    # This c9tv is used for interlinking in the GGB!
                    targetfref = _get_rel_fref(org1,GSofCBGs.input,OPTIONS)
                else:
                    continue

                # get length of pacbporf
                length = len(thepacbporf ._positions)
                if hasattr(thepacbporf ,'_original_alignment_pos_start'):
                    length = thepacbporf ._original_alignment_pos_end -\
                             thepacbporf ._original_alignment_pos_start

                # gff dict for this feature
                gffdict = { 'fref'        : fref_rel,
                        'fsource'     : GFF_ORFSIMILARITY_FSOURCE,
                        'fmethod'     : GFF_ORFSIMILARITY_FMETHOD,
                        'gclass'      : GFF_ORFSIMILARITY_GCLASS,
                        'column9data' : {
                            'Note'      : "Similarity on ORF %s with ORF %s" % (thepacbporf.orfQ.id,thepacbporf.orfS.id),
                            'Length'    : length,
                            'Bitscore'  : thepacbporf.bitscore,
                            'Similarity': thepacbporf.similarity,
                            'Identity'  : thepacbporf.identity,
                            'TargetFref': targetfref,
                            },
                        'gname' : group,
                    }
                # and make the GFF line
                gffdata.append( thepacbporf.togff(gff=gffdict,with_alignment=False) )

                if GFF_ORFSIM_PLS_OUTPUT or GFF_ORFSIM_PTS_OUTPUT:
                    # fill gff dictionary from setting
                    _gff_pls = { 'fref'    : fref_rel,
                                 'fsource' : GFF_ORFSIM_PLS_FSOURCE,
                                 'fmethod' : GFF_ORFSIM_PLS_FMETHOD,
                                 'column9data': {'TargetFref': targetfref },
                                 'gclass'  : GFF_ORFSIM_PLS_GCLASS, }
                    _gff_pts = { 'fref'    : fref_rel,
                                 'fsource' : GFF_ORFSIM_PTS_FSOURCE,
                                 'fmethod' : GFF_ORFSIM_PTS_FMETHOD,
                                 'column9data': {'TargetFref': targetfref },
                                 'gclass'  : GFF_ORFSIM_PTS_GCLASS, }
                    # check if dict key `group` exists already
                    if not projected_leading_stops.has_key(group):
                        projected_leading_stops[group] = []
                        projected_tailing_stops[group] = []

                    # make gffs of projected ORF stop sites
                    projected_leading_stops[group].append( thepacbporf.gff_projected_leading_stop(gff=_gff_pls) )
                    projected_tailing_stops[group].append( thepacbporf.gff_projected_tailing_stop(gff=_gff_pts) )

        # take the projection of the most LEFT cbg
        for group, gfflist in projected_leading_stops.iteritems():
            gffdata.append(gfflist[0])

        # take the projection of the most RIGTH cbg
        for group, gfflist in projected_tailing_stops.iteritems():
            gffdata.append(gfflist[-1])

        # (A) prepare the RELATIVE GFF file
        # for sake of check-double-check, set fref_rel to all gff tuples
        for i in range(0,len(gffdata)):
            gffl       = list(gffdata[i])
            gffl[0]    = fref_rel
            gffdata[i] = tuple(gffl)

        # now write the result file
        fname  = "%sresult.rel.%s.gff" % ( basefname, thisorg ) 
        # make fullpath of fname 
        fullpath = os.path.join(OPTIONS.outdir,fname)
        # check if file exists and overwriting is aloued
        if not OPTIONS.force:
            fullpath = _create_unique_filename_with_integer_suffix(fullpath)

        # write output file
        fh = open(fullpath, 'w')
        fh.write(gffs2txt(gffdata))
        fh.close()

        # update created_files dictionary
        created_files[thisorg] = fullpath

        # (B) prepare the ABSOLUTE GFF file when applicable
        if GSofCBGs.input[thisorg]['locusfile']:
            # get locusgff line
            locusgff = parsegfffile(GSofCBGs.input[thisorg]['locusfile'])[0]

            if locusgff[6] == '-':
                # mirror the gff data; negative strand!
                locuslength= ( locusgff[4]-locusgff[3]+1 )
                gffdata = mirrorgfflist(gffdata,length=locuslength)

            # set absolute fref and correct for locus coordinate offset
            for i in range(0,len(gffdata)):
                gffl       = list(gffdata[i])
                gffl[0]    = fref_abs
                gffl[3]    = gffl[3]+locusgff[3]
                gffl[4]    = gffl[4]+locusgff[3]
                # convert back from list to tuple
                gffdata[i] = tuple(gffl)

            # write the absolute result file
            fname  = "%sresult.abs.%s.gff" % ( basefname, thisorg ) 
            # make fullpath of fname 
            fullpath = os.path.join(OPTIONS.outdir,fname)
            # check if file exists and overwriting is aloued
            if not OPTIONS.force:
                fullpath = _create_unique_filename_with_integer_suffix(fullpath)

            # write output file
            fh = open(fullpath, 'w')
            fh.write(gffs2txt(gffdata))
            fh.close()

    # return performance outcome & created_files 
    return performance_data_dict, created_files

# end of function create_abgp_prediction_gff_files


def create_abgp_prediction_protein_fasta_files(GSofCBGs,OPTIONS,organism=None,verbose=False,_linesize=60):
    """ """
    allproteins = Set([GSofCBGs.input[thisorg]['proteinfref'] for thisorg in GSofCBGs.organism_set()]) 
    for thisorg in GSofCBGs.organism_set():
        if organism and thisorg != organism: continue
        # get the gff outcome from the GSG
        outcome = GSofCBGs.get_genestructure_gff(organism=thisorg)
        # get dna sequence belonging to these tracks
        dna = []
        for gfftuple in outcome:
            dna.append( getseqbyfeaturetuple(GSofCBGs.input[thisorg]['genomeseq'],gfftuple) )
        # join into dnaseq and translate into protein
        dnaseq  = "".join(dna)
        protseq = dna2prot.dna2protein(dnaseq)

        # create absolute filename
        fref  = GSofCBGs.input[thisorg]['proteinfref']
        fname = "%sprotein.%s.fa" % ( _get_abgp_file_basename(OPTIONS), fref )
        fname = os.path.join(OPTIONS.outdir,fname)
        # check if file exists and overwriting is aloued
        if not OPTIONS.force:
            fname = _create_unique_filename_with_integer_suffix(fname)

        # create used informants string
        informants = allproteins.difference([ GSofCBGs.input[thisorg]['proteinfref'] ]) 
        # create fasta protein file
        fh = open(fname,'w')
        fh.write(">%s [%s informants (%s): %s]\n" % (
                GSofCBGs.input[thisorg]['proteinfref'],
                ABGP_VERSION,
                len(informants),
                " ".join([str(informant) for informant in informants]),
                ))
        # and append the protein sequence to the file
        for i in range(0,len(protseq),_linesize):
            fh.write(protseq[i:i+_linesize]+"\n")
        fh.close()

# end of function create_abgp_prediction_protein_fasta_files


def create_abgp_outcome_txt_file(GSG,OPTIONS,performance_data_dict):
    """ """
    fname = "%soutcome.txt" % ( _get_abgp_file_basename(OPTIONS) )
    fname = os.path.join(OPTIONS.outdir,fname) 
    # check if file exists and overwriting is aloued
    if not OPTIONS.force:
        fname = _create_unique_filename_with_integer_suffix(fname)
    fh = open(fname,'w')
    fh.write("# %s\n" % GSG.genetree())
    # print information about targets and informants used
    if OPTIONS.target:
        # print information about targets and informants used
        informants = GSG.input.keys()
        informants.remove( OPTIONS.target )
        informants.sort()
        fh.write("# --target %s informants used: %s\n" % (OPTIONS.target,informants) )

    # print the ABGP_predictor for its own prediction and confirmation on genes/unigenes
    maxlenorg  = max([ len(org) for org in GSG.input.keys() ])
    maxlenprot = max([ len(GSG.input[org]['proteinfref']) for org in GSG.input.keys() ])
    for organism in GSG.input.keys():
        # get list of ABGP prediction outcome on each of the interfaces, S and E
        ABGP_PREDICTION_LIST   = GSG.abgp_prediction_status(organism,verbose=False)
        abgp_prediction_string = [ str(item)[0] for item in ABGP_PREDICTION_LIST] 
        for i in range(1,len(abgp_prediction_string),2):
            abgp_prediction_string[i] = abgp_prediction_string[i].lower()
        abgp_prediction_string = "".join(abgp_prediction_string)
        if False in ABGP_PREDICTION_LIST:           ABGP_PREDICTION_STATUS = False
        elif ABGP_PREDICTION_LIST.count(None) <= 1: ABGP_PREDICTION_STATUS = True
        elif ABGP_PREDICTION_LIST.count(None) <= 2: ABGP_PREDICTION_STATUS = None
        else:                                       ABGP_PREDICTION_STATUS = False

        gene_status    = "     "
        unigene_status = "     "
        if performance_data_dict.has_key(organism) and performance_data_dict[organism].has_key('gene'):
            gene_status = performance_data_dict[organism]['gene']
        if performance_data_dict.has_key(organism) and performance_data_dict[organism].has_key('unigene'):
            unigene_status = performance_data_dict[organism]['unigene']
        # write results to fh
        orgstring  = organism+" "*(maxlenorg*2)
        orgstring  = orgstring[0:maxlenorg+4]
        protstring = GSG.input[organism]['proteinfref']+" "*(maxlenprot*2)
        protstring = protstring[0:maxlenprot+4]
        txtstring = "%s%sABGP_PREDICTION: %s\tgeneTP: %s\tunigeneTP: %s \t%s" % (
                orgstring, protstring, ABGP_PREDICTION_STATUS, gene_status, unigene_status,abgp_prediction_string)
        fh.write(txtstring+"\n")
    # close filehandle and return filename
    fh.close()
    return fname

# end of function create_abgp_outcome_txt_file


def get_gff_of_gsg2cexpander(gsg,organism,gff={}):
    """ """
    GFF_CBG_CEXPANDER_FSOURCE = 'cexpander'
    GFF_CBG_CEXPANDER_FMETHOD = 'cexpander'
    GFF_CBG_CEXPANDER_GCLASS  = 'Cexpander'

    if not gff.has_key('fref'):    gff['fref']    = None
    if not gff.has_key('fsource'): gff['fsource'] = GFF_CBG_CEXPANDER_FSOURCE
    if not gff.has_key('fmethod'): gff['fmethod'] = GFF_CBG_CEXPANDER_FMETHOD
    if not gff.has_key('gclass'):  gff['gclass']  = GFF_CBG_CEXPANDER_GCLASS

    frame2ggbcoordcor = {0: 3+1, 1:1+1, 2:2+1 }

    gffdata = []
    for cbg in gsg:
        if organism not in cbg.organism_set(): continue

        # get Orf of this cbg do get it phase/frame
        orf = cbg.get_orfs_of_graph(organism=organism)[0]
        omsr = list(cbg.overall_minimal_spanning_range(organism=organism))
        if cbg.__class__.__name__ == 'LowSimilarityRegionCodingBlockGraph':
            omsr.sort()
            for aacoord in omsr:
                fstart = aacoord*3+frame2ggbcoordcor[orf.frame]
                fstop  = fstart + 2
                gfftrack = ( gff['fref'], gff['fsource'], gff['fmethod'],
                             fstart,fstop,0,'+','.',
                             "%s %s" % (gff['gclass'],gff['fref']) )
                gffdata.append(gfftrack)
        else:
            trfblock = cbg._cexpander.get_transferblock(organism)
            if trfblock:
                binarystring = trfblock.binarystring
            else:
                # cexpander FAILED on this batch of sequences!
                binarystring = "0"*len(cbg._cexpander.sequences[organism])

            # calculate fstart & associated fstop and make tracks
            fstart = min(omsr)*3 + frame2ggbcoordcor[orf.frame]
            fstop  = fstart + 2
            for binarysymbol in binarystring:
                gfftrack = ( gff['fref'], gff['fsource'], gff['fmethod'],
                             fstart,fstop,int(binarysymbol),'+','.',
                             "%s %s" % (gff['gclass'],gff['fref']) )
                gffdata.append(gfftrack)
                fstart+=3
                fstop+=3

    # return the gffdata
    return gffdata

# end of function get_gff_of_gsg2cexpander


def get_gff_of_gsg2multiplealignment(gsg,organism,gff={}):
    """ """
    GFF_CBG_MULTIPLEALIGNMENT_FSOURCE = 'cbg-abfgp'
    GFF_CBG_MULTIPLEALIGNMENT_FMETHOD = 'cbgMAPR'
    GFF_CBG_MULTIPLEALIGNMENT_GCLASS  = 'cbgMAPR'

    from graphAbgp.codingblock_multiplealignment import cbgproteinalignmentdata
    if not gff.has_key('fref'):    gff['fref']    = None
    if not gff.has_key('fsource'): gff['fsource'] = GFF_CBG_MULTIPLEALIGNMENT_FSOURCE 
    if not gff.has_key('fmethod'): gff['fmethod'] = GFF_CBG_MULTIPLEALIGNMENT_FMETHOD
    if not gff.has_key('gclass'):  gff['gclass']  = GFF_CBG_MULTIPLEALIGNMENT_GCLASS

    ntdict = {}
    for cbg in gsg:
        if cbg.__class__.__name__ == 'LowSimilarityRegionCodingBlockGraph':
            continue

        # get Orf of this cbg to get it phase/frame
        orf = cbg.get_orfs_of_graph(organism=organism)[0]
        aadict = cbgproteinalignmentdata(cbg,organism=organism)
        # convert aadict to ntdict and store to main ntdict
        for aapos,ratio in aadict.iteritems():
            ntpos = aapos*3 + orf.frame             
            if ntdict.has_key(ntpos):
                if ntdict[ntpos] < ratio:
                    ntdict[ntpos] = ratio
            else:
                # check if 2-nt offset has a higher ratio
                # This can happend when Orfs overlap....
                for offset in range(ntpos+1,ntpos+3):
                    if ntdict.has_key(offset):
                        if ntdict[offset] < ratio:
                            del(ntdict[offset])
                            ntdict[ntpos] = ratio
                            break
                        else:
                            break
                else:
                    # if here, no position in this region known
                    ntdict[ntpos] = ratio
                    

    # convert all positions into GFF
    ntpositions = ntdict.keys()
    ntpositions.sort()
    gffdata = []
    for ntpos in ntpositions:
        fstart,fstop = ntpos, ntpos+2
        gfftrack = ( gff['fref'], gff['fsource'], gff['fmethod'],
                     fstart,fstop,
                     "%1.3f" % ntdict[ntpos],'+','.',
                     "%s %s" % (gff['gclass'],gff['fref']) )

        gffdata.append(gfftrack)
    # return the gffdata
    return gffdata

# end of function get_gff_of_gsg2multiplealignment



def gff2database(fnamegff,fnamefasta):
    """
    Load a gff file (and its fasta) into a GGB database

    @type  fnamegff: string
    @param fnamegff: full path to gff file to store

    @type  fnamefasta: *
    @param fnamefasta: full path to gff file to store or None

    @rtype:  Boolean
    @return: True if succesfull, False otherwise
    """
    return gff2db(
            dbname=DB_NAME,
            dbhost=DB_HOST,
            dbconnection=DB_CONNECTION,
            dbuser=DB_USER,
            dbpass=DB_PASS,
            fnamefasta=fnamefasta,
            fnamegff=fnamegff,
            pathperl=PERL_PATH,
            pathloadgff=EXECUTABLE_LOAD_GFF)

# end of function gff2database

def gff_orf_check(graphlist,input={},org="mgg",basename='orfcheck'):
    """
    """

    FREF   = ".".join(input[org]['proteinseqfile'].split(".")[0:-1])+"_dna"
    FREF   = input[org]['FREF']
    METHOD = input[org]['METHOD']

    gffs = input[org]['gffs']
    orfs = {}
    for sg in graphlist:
        for ( (a,b,c,d),(org1,orf1),(org2,orf2) ), pacbporf in sg.pacbps.iteritems():
            if org1 == org:
                if orfs.has_key(orf1): pass
                else: orfs[orf1] = pacbporf.orfQ
            if org2 == org:
                if orfs.has_key(orf2): pass
                else: orfs[orf2] = pacbporf.orfS


    knownorfs = []
    for gff in gffs:
        if gff[2] == 'known_orf':
            print "\t".join([str(item) for item in gff])
            ###fh.write("%s\n" % "\t".join([str(item) for item in gff]) )
            knownorfs.append( gff[-1].split(" ",1)[1] )

    thisorfs = []
    for orfid,orfObject in orfs.iteritems():
        gff = orfObject.togff(fref=FREF,
                        fmethod='abfgp_orf',
                        fsource='abfgp',
                        gclass='OrfWithSimilarity',
                        one_group_per_frame=False
                        )
        gffs.append(gff)
        print "\t".join([str(item) for item in gff])
        ###fh.write("%s\n" % "\t".join([str(item) for item in gff]) )
        thisorfs.append( gff[-1].split(" ",1)[1] )


    # make sets of orfs
    knownorfs = Set(knownorfs)
    thisorfs  = Set(thisorfs)

    # and append to a specific file
    fname = "/home/avdb/data/abfgp-saas/.%s_%s.gff" % (basename,org)
    txt = open(fname).read()
    fh = open(fname, 'w')
    fh.write(txt)
    fh.write("%s\t%s\t%s\t%s\t%s\t%s\tmissed: %s\textra: %s\n" % (
            os.getcwd(),
            gffs[0][0],
            len(knownorfs),
            knownorfs.difference(thisorfs) == Set([]),
            len( knownorfs.difference(thisorfs) ),
            len( thisorfs.difference(knownorfs) ),
            ", ".join([ orf.split("-")[-1] for orf in knownorfs.difference(thisorfs) ]),
            ", ".join([ orf.split("-")[-1] for orf in thisorfs.difference(knownorfs) ])
            )
        )
    fh.close()
    

# end of function gff_orf_check



def create_intron_mapping_gff_file(input,PCG,introndata,assessed_interfaces,OPTIONS,FREF,organism=None,verbose=False):
    """
    @type  OPTIONS: optparse options object
    @param OPTIONS: optparse options object

    @type  organism: *
    @param organism: Organism identifier (presumably string)

    @type  verbose: Boolean 
    @param verbose: print (timing) messages to STDOUT (True) or not (False) 

    @rtype:  tuple
    @return: tuple of 2 created file names with relative paths ( gff_file_name, fasta_file_name )

    @attention: Dozens/hundreds of Global variables needed (imported from settings file)
    """
    if not organism: raise "organism not specified"

    stw = StopWatch(name="CreateGffStopWatch")
    stw.start()

    # get all gff lines that are already gathered for this organism
    gffs   = input[organism]['gffs']
    METHOD = input[organism]['METHOD']

    PROJECT_FIRSTEXONONORF_MIN_NT_LENGTH = 15
    PROJECT_FINALEXONONORF_MIN_NT_LENGTH = 60

    OPTIONS.force = True

    # data about the annotated gene structure
    geneobjects = input[organism]['gldobj'].as_exons()
    known_introns = [ geneobjects[pos] for pos in range(1,len(geneobjects),2) ]
    known_exons   = [ geneobjects[pos] for pos in range(0,len(geneobjects),2) ]
    known_orf_ids = [ exon.orf.id for exon in known_exons ] 

    ####################################################################
    if verbose: print stw.lap(), "preparation"
    ####################################################################

    # get gff lines for each ORF
    for orf in input[organism]['orfs'].orfs:
        if orf.id in known_orf_ids:
            gfflist = list( orf.togff(fref=FREF) )
            gfflist[-1] = gfflist[-1] + "; knownOrf True"
            gffs.append( tuple(gfflist) )
        else:
            gffs.append( orf.togff(fref=FREF) )

    ####################################################################
    if verbose: print stw.lap(), "Orfs"
    ####################################################################

    # dictionaries for projected stop codons of a pacbp
    projected_leading_stops = {}
    projected_tailing_stops = {}
    projected_eof_protein = {}
    projected_sof_protein = {}
    projected_missed_final_exons = {}
    projected_missed_first_exons = {}

    # get gff lines for each similarity/pacbp

    for orgS in PCG.organism_set():
        if organism == orgS: continue

        # check if informant isa unigene; slightly different treatment here
        # compared to a genomic DNA informant
        IS_UNIGENE = False
        if input[orgS].has_key('is_unigene') and input[orgS]['is_unigene']:
            IS_UNIGENE = True

        if IS_UNIGENE:
            orgSfullname = "%s[%s-%s]" % (
                    _get_organism_full_name(input[orgS]['is_unigene_of'],truncate=True), 
                    orgS.split("[")[0],
                    input[orgS]['unigene-annotation'] )
        else:
            orgSfullname = "%s [%s]" % (
                _get_organism_full_name(orgS,truncate=True),
                input[orgS]['proteinfref'] )

        # get ordered PacbPORF list of this informant
        orderedpacbporfs = pacb.ordering.order_pacbporf_list(PCG.get_pacbps_by_organisms(organism,orgS))

        for _pos in range(0,len(orderedpacbporfs)):
            thepacbporf = orderedpacbporfs[_pos]
            group       = "%s-%s-orfid-%s" % (FREF,orgS,thepacbporf.orfS.id)
            # get the TargetFref of the other organism.
            # This c9tv is used for interlinking in the GGB!
            targetfref = _get_rel_fref(orgS,input,OPTIONS)

            if IS_UNIGENE and _pos >= 1:
                # Check if there is a gap in the chained UniGene.
                # In theory, this should not occur. It is a strong signal that
                # part of the alignment is missing. A biological
                # explanation for a unigene gap is intron retention.
                prevpacbporf = orderedpacbporfs[_pos-1]
                prev_end    = prevpacbporf._get_original_alignment_pos_end().sbjct_pos
                this_start  = thepacbporf._get_original_alignment_pos_start().sbjct_pos
                gap_aa_dist = this_start - prev_end
                if gap_aa_dist >= 8:
                    q_prev_end    = prevpacbporf._get_original_alignment_pos_end().query_pos
                    q_this_start  = thepacbporf._get_original_alignment_pos_start().query_pos
                    gffline = list(gffs[-1])
                    
                    # create a dummy similarity track in the middle!
                    query_nt_gap = thepacbporf._get_original_alignment_pos_start().query_dna_start - gffline[4]
                    dummy_start  = gffline[4] + (( query_nt_gap - (gap_aa_dist*3) ) / 2)
                    dummy_stop   = dummy_start + (gap_aa_dist*3) - 1
                    dummy_score  = 0.0
                    gffline[2] = "ProjectedUnigeneGap"
                    gffline[3] = dummy_start
                    gffline[4] = dummy_stop
                    gffline[5] = dummy_score
                    gffline[8] = gffline[8].split("; ")[0]
                    gffs.append(tuple(gffline))

            # get length of pacbporf
            length = len(thepacbporf._positions)
            if hasattr(thepacbporf ,'_original_alignment_pos_start'):
                length = thepacbporf._original_alignment_pos_end -\
                         thepacbporf._original_alignment_pos_start
    

            column9data =  {
                        #'Note'      : "Similarity on ORF %s with ORF %s" % (thepacbporf.orfQ.id,thepacbporf.orfS.id),
                        'Length'    : length,
                        'Bitscore'  : thepacbporf.bitscore,
                        'Similarity': thepacbporf.similarity,
                        'Identity'  : thepacbporf.identity,
                        'TargetFref': targetfref,
                        'SbjctPos'  : "%s-%s" % (
                                thepacbporf._get_original_alignment_pos_start().sbjct_pos,
                                thepacbporf._get_original_alignment_pos_end().sbjct_pos ),
                        'Informant' : "%s orf(%s)" % (
                                        orgSfullname,
                                        thepacbporf.orfS.id,
                                        )
                        }

            # gff dict for this feature          
            gff = { 'fref'        : FREF, 
                    ###'fsource'     : GFF_ORFSIMILARITY_FSOURCE,
                    'fmethod'     : GFF_ORFSIMILARITY_FMETHOD,
                    'gclass'      : GFF_ORFSIMILARITY_GCLASS,
                    'gname'       : group,
                    'column9data' : column9data,
                    }
            sim_source = GFF_ORFSIM_PLS_FSOURCE

            # check if it is an unigene informant -> other fmethod
            if IS_UNIGENE:
                gff = { 'fref'        : FREF, 
                        'fsource'     : GFF_UNIGENESIMILARITY_FSOURCE,
                        'fmethod'     : GFF_UNIGENESIMILARITY_FMETHOD,
                        'gclass'      : GFF_UNIGENESIMILARITY_GCLASS,
                                            'gname' : group,
                        'column9data' : column9data,
                        }
                # switch projected stop codons as well
                sim_source = GFF_UNIGENESIMILARITY_FSOURCE

            # and make the GFF line
            gffs.append( thepacbporf.togff(gff=gff,with_alignment=False) )

            # fill gff dictionary from setting
            _gff_pls = { 'fref'    : FREF,
                         'fsource' : sim_source,
                         'fmethod' : GFF_ORFSIM_PLS_FMETHOD,
                         'column9data' : column9data,                        
                         'gclass'  : GFF_ORFSIM_PLS_GCLASS, }
            _gff_pts = { 'fref'    : FREF,
                         'fsource' : sim_source,
                         'fmethod' : GFF_ORFSIM_PTS_FMETHOD,
                         'column9data' : column9data,                        
                         'gclass'  : GFF_ORFSIM_PTS_GCLASS, }
            # check if dict key `group` exists already
            if not projected_leading_stops.has_key(group):
                projected_leading_stops[group] = []
                projected_tailing_stops[group] = []

            # make gffs of projected ORF stop sites
            projected_leading_stops[group].append( thepacbporf.gff_projected_leading_stop(gff=_gff_pls) )
            projected_tailing_stops[group].append( thepacbporf.gff_projected_tailing_stop(gff=_gff_pts) )


            # Check if this represents an alignment of the FirstExonOnOrf
            # of the annotated gene structure of the informant
            # Make shure only the FIRST encountered occurrence of this Orf
            # can give a projection (projected_sof_protein.has_key(group))
            if not IS_UNIGENE and input[orgS]['orfid-genestructure'] and\
            thepacbporf.orfS.id == input[orgS]['orfid-genestructure'][0] and\
            not projected_sof_protein.has_key(group) and\
            not thepacbporf.orfS.has_dna_n_symbols():
                firstexon = input[orgS]['gldobj'].as_exons()[0]
                staPosObj = thepacbporf._get_original_alignment_pos_start()
                tss_stapos_dist   = staPosObj.sbjct_dna_start - firstexon.acceptor.pos
                projected_tss_pos = staPosObj.query_dna_start - tss_stapos_dist
                if tss_stapos_dist > PROJECT_FIRSTEXONONORF_MIN_NT_LENGTH:
                    # >5 AAs missing!
                    _gff_sof = { 'fref'    : FREF,
                                 'fsource' : sim_source,
                                 'gclass'  : 'ProjectedFirstExonOnOrf' }
                    dummy = list(thepacbporf.gff_projected_leading_stop(gff=_gff_sof))
                    dummy[2] = 'ProjectedFirstExonOnOrf'
                    dummy[3] = projected_tss_pos
                    dummy[4] = staPosObj.query_dna_start
                    # store to projected_sof_protein dict
                    projected_sof_protein[group] = tuple(dummy)
                else:
                    # start site of informants' FirstExon more or less
                    # covered -> do NOT create projection here!
                    projected_sof_protein[group] = None

            # Check if this represents an alignment of the FinalExonOnOrf
            # of the annotated gene structure of the informant
            if not IS_UNIGENE and input[orgS]['orfid-genestructure'] and\
            thepacbporf.orfS.id == input[orgS]['orfid-genestructure'][-1] and\
            not thepacbporf.orfS.has_dna_n_symbols():
                endPosObj = thepacbporf._get_original_alignment_pos_end()
                endpos_stop_dist = thepacbporf.orfS.endPY - endPosObj.sbjct_dna_end
                projected_stop_start = endPosObj.query_dna_end + endpos_stop_dist
                if endpos_stop_dist > PROJECT_FINALEXONONORF_MIN_NT_LENGTH:
                    # >20 AAs missing!
                    _gff_eof = { 'fref'    : FREF,
                                 'fsource' : sim_source,
                                 'gclass'  : 'ProjectedFinalExonOnOrf' }
                    dummy = list(thepacbporf.gff_projected_tailing_stop(gff=_gff_eof))
                    dummy[2] = 'ProjectedFinalExonOnOrf'
                    dummy[3] = endPosObj.query_dna_end
                    dummy[4] = projected_stop_start
                    if projected_eof_protein.has_key(group):
                        projected_eof_protein[group].append( tuple(dummy) )
                    else:
                        projected_eof_protein[group] = [ tuple(dummy) ]
                else:
                    # Nope; no AAs missing -> whipe out existing
                    # projections from more 5' located PacbPORF alignments.
                    if projected_eof_protein.has_key(group):
                        del( projected_eof_protein[group] )


        # check if there are FirstExons/FinalExons completely missed
        if input[orgS]['orfid-genestructure']:
            known_orfids = input[orgS]['orfid-genestructure']
            informant_orfids = [ pf.orfS.id for pf in orderedpacbporfs ]

            if not informant_orfids:
                # no annotated Orfs known/given for this informant
                pass
            elif informant_orfids[-1] in known_orfids:
                final_index = known_orfids.index(informant_orfids[-1])
                if final_index < len(known_orfids)-1:
                    # get exons of this informant genestructure
                    inf_geneobjects = input[orgS]['gldobj'].as_exons()
                    inf_known_exons   = [ inf_geneobjects[pos] for pos in range(0,len(inf_geneobjects),2) ]

                    # get final pacbporf of this informant for
                    # making relative projection
                    finalpf = orderedpacbporfs[-1]
                    projection_offset = finalpf.orfQ.endPY + 10
                    # make dict element for gff projection tracks
                    projected_missed_final_exons[orgS] = []

                    # loop over the exons and make tracks for those
                    # who are not included in the ABFGP structure
                    for _pos in range(len(inf_known_exons)-1,-1,-1):
                        infe = inf_known_exons[_pos]
                        if infe.orf.id == informant_orfids[-1]: break
                        print "COMPLETELY_MISSED FINAL:", orgS, infe

                        # define a group for this failed exon
                        group = "%s-%s-orfid-%s" % (FREF,orgS,infe.orf.id)
                        # fill gff dictionary from setting
                        _gff_pls = { 'fref'    : FREF,
                                     'fmethod' : GFF_ORFSIM_PLS_FMETHOD,
                                     'gname'   : group,
                                     'gclass'  : "Similarity", }
                        _gff_pts = { 'fref'    : FREF,
                                     'fmethod' : GFF_ORFSIM_PTS_FMETHOD,
                                     'gname'   : group,
                                     'gclass'  : "Similarity", }
            
                        # make gffs of projected ORF stop sites
                        dummy2 = list(finalpf.gff_projected_leading_stop(gff=_gff_pls))
                        dummy3 = list(finalpf.gff_projected_tailing_stop(gff=_gff_pts))
                        dummy3[3] = projection_offset + (dummy3[3]-dummy2[3])
                        dummy3[4] = dummy3[3]+2
                        dummy3[8] = dummy3[8].split("; ")[0]
                        dummy2[3] = projection_offset + 1
                        dummy2[4] = projection_offset + 3
                        dummy2[8] = dummy3[8].split("; ")[0]
                        _gff_eof = { 'fref'    : FREF,
                                     'gclass'  : "Similarity",
                                     'gname'   : group }
                        dummy1 = list( infe.togff(gff=_gff_eof) )
                        dummy1[1] = 'getorf-TBLASTX'
                        dummy1[2] = 'ProjectedFinalExonOnOrf'
                        dummy1[3] = projection_offset + (infe.acceptor.pos-infe.orf.startPY)+1
                        dummy1[4] = dummy1[3] + infe.length
                        projected_missed_final_exons[orgS].append(tuple(dummy1))
                        projected_missed_final_exons[orgS].append(tuple(dummy2))
                        projected_missed_final_exons[orgS].append(tuple(dummy3))
                        #print dummy2
                        #print dummy1
                        #print dummy3
            else:
                # no missing Orf for this informant
                pass


            if not informant_orfids:
                # no annotated Orfs known/given for this informant
                pass
            elif informant_orfids[0] in known_orfids:
                first_index = known_orfids.index(informant_orfids[0])
                if first_index == 1:
                    # get exons of this informant genestructure
                    inf_geneobjects = input[orgS]['gldobj'].as_exons()
                    inf_known_exons   = [ inf_geneobjects[pos] for pos in range(0,len(inf_geneobjects),2) ]

                    # get first pacbporf of this informant for
                    # making relative projection
                    firstpf = orderedpacbporfs[0]
                    projection_offset = firstpf.orfQ.startPY - 10
                    # make dict element for gff projection tracks
                    projected_missed_first_exons[orgS] = []

                    # loop over the exons and make tracks for those
                    # who are not included in the ABFGP structure
                    for _pos in range(0,len(inf_known_exons)):
                        infe = inf_known_exons[_pos]
                        if infe.orf.id == informant_orfids[0]: break
                        print "COMPLETELY_MISSED FIRST:", orgS, infe
                        # define a group for this failed exon
                        group = "%s-%s-orfid-%s" % (FREF,orgS,infe.orf.id)
                        # fill gff dictionary from setting
                        _gff_pls = { 'fref'    : FREF,
                                     'fmethod' : GFF_ORFSIM_PLS_FMETHOD,
                                     'gname'   : group,
                                     'gclass'  : "Similarity", }
                        _gff_pts = { 'fref'    : FREF,
                                     'fmethod' : GFF_ORFSIM_PTS_FMETHOD,
                                     'gname'   : group,
                                     'gclass'  : "Similarity", }
            
                        # make gffs of projected ORF stop sites
                        dummy2 = list(firstpf.gff_projected_leading_stop(gff=_gff_pls))
                        dummy3 = list(firstpf.gff_projected_tailing_stop(gff=_gff_pts))
                        dummy3[3] = projection_offset -3 + 1
                        dummy3[4] = projection_offset
                        dummy3[8] = dummy3[8].split("; ")[0]
                        dummy2[3] = projection_offset - infe.orf.length - 3 -3 + 1
                        dummy2[4] = projection_offset - infe.orf.length - 3
                        dummy2[8] = dummy3[8].split("; ")[0]
                        _gff_sof = { 'fref'    : FREF,
                                     'gclass'  : "Similarity",
                                     'gname'   : group }
                        dummy1 = list( infe.togff(gff=_gff_sof) )
                        dummy1[1] = 'getorf-TBLASTX'
                        dummy1[2] = 'ProjectedFirstExonOnOrf'
                        dummy1[4] = projection_offset - (infe.orf.endPY-infe.donor.pos)
                        dummy1[3] = dummy1[4] - infe.length + 1
                        projected_missed_first_exons[orgS].append(tuple(dummy1))
                        projected_missed_first_exons[orgS].append(tuple(dummy2))
                        projected_missed_first_exons[orgS].append(tuple(dummy3))
                        #print dummy2
                        #print dummy1
                        #print dummy3
                        # check if there are SignalPeptides that can be
                        # projected too
                        if infe.orf._signalp_sites:
                            for spS in infe.orf._signalp_sites:
                                prjQSPstart = (spS.start-infe.start) + dummy1[3]+1
                                prjQSPend   = (spS.start-infe.start) + dummy1[3]+1 + spS.length
                                prjQTSSstart= (spS.start-infe.start) + dummy1[3]+1 - IC_TSS_PATTERN_OFFSET[0]

                                # create a ProjectedSignalPSignalPeptide
                                prjTSS = TranslationalStartSite(
                                    prjQTSSstart,"n"*19,
                                    pssm_score=spS.tss.pssm_score )
                                prjTSS._gff['fmethod'] = 'projectedTSSpssm'
                                prjTSS._gff['column9data'] = {'Informant': orgSfullname}
                                prjSignalP_gff = {
                                    'fmethod'    : 'projectedSignalPeptide',
                                    'column9data': {'Informant': orgSfullname},
                                    'gname'      : "%s_%s_%s" % (orgS,prjQSPstart,prjQSPend) }
                                prjSignalP = ProjectedSignalPSignalPeptide(
                                    prjQSPstart,prjQSPend,
                                    spS.pssm_score,tss=prjTSS,
                                    gff=prjSignalP_gff )
                                prjSignalP._gff['fmethod'] = 'projectedSignalPeptide'
                                prjSignalP._gff['column9data'] = {'Informant': orgSfullname}
                                # append to projected_missed_first_exons
                                projected_missed_first_exons[orgS].extend(
                                    prjSignalP.togff(gff={'fref':FREF}) )

            else:
                # no missing Orf for this informant
                pass


    ####################################################################
    if verbose: print stw.lap(), "similarity"
    ####################################################################

    # take the projection of the most LEFT cbg
    for group, gfflist in projected_leading_stops.iteritems():
        gffs.append(gfflist[0]) 
    # take the projection of the most RIGTH cbg
    for group, gfflist in projected_tailing_stops.iteritems():
        gffs.append(gfflist[-1])

    ####################################################################
    if verbose: print stw.lap(), "projected * codons"
    ####################################################################

    # take projected unaligned FirstExonOnOrf fraction if large enough
    for group, gfftupleornone in projected_sof_protein.iteritems():
        # if gfftupleornone == None -> no FirstExonOnOrf projection!
        if gfftupleornone:
            gffs.append(gfftupleornone)

    ####################################################################
    if verbose: print stw.lap(), "projected unaligned FirstExonOnOrf"
    ####################################################################

    # take projected unaligned FinalExonOnOrf fraction if large enough
    for group, gfflist in projected_eof_protein.iteritems():
        gffs.append(gfflist[-1])

    ####################################################################
    if verbose: print stw.lap(), "projected unaligned FinalExonOnOrf"
    ####################################################################

    for orgS, exongfflist in projected_missed_final_exons.iteritems():
        for inf_gff_exon in exongfflist:
            gffs.append( inf_gff_exon )

    for orgS, exongfflist in projected_missed_first_exons.iteritems():
        for inf_gff_exon in exongfflist:
            gffs.append( inf_gff_exon )


    ####################################################################
    if verbose: print stw.lap(), "projected missed First/FinalExonOnOrf"
    ####################################################################


    seen_donors = []
    seen_acceps = []
    #splicesite_gff = { 'fref': FREF, 'fsource':'PSSM4FUNGI' }
    first_donor_pos = None
    final_donor_pos = None
    for orgS in PCG.organism_set():
        if organism == orgS: continue
        for thepacbporf in pacb.ordering.order_pacbporf_list(PCG.get_pacbps_by_organisms(organism,orgS)):
            for donor in thepacbporf.orfQ._donor_sites:
                if donor.pos not in seen_donors:
                    seen_donors.append( donor.pos )
                    #donor._gff['gname'] = donor.pos
                    splicesite_gff = { 'fref': FREF, 'fsource':'PSSM4FUNGI' }
                    gffs.append( donor.togff( gff=splicesite_gff ) )
                    if not first_donor_pos:
                        first_donor_pos = donor.pos
                    else:
                        if donor.pos < first_donor_pos:
                            first_donor_pos = donor.pos
                    if not final_donor_pos:
                        final_donor_pos = donor.pos
                    else:
                        if donor.pos > final_donor_pos:
                            final_donor_pos = donor.pos
            for acceptor in thepacbporf.orfQ._acceptor_sites:
                if acceptor.pos not in seen_acceps:
                    seen_acceps.append( acceptor.pos )
                    #acceptor._gff['gname'] = acceptor.pos
                    splicesite_gff = { 'fref': FREF, 'fsource':'PSSM4FUNGI' }
                    gffs.append( acceptor.togff( gff=splicesite_gff ) )

    # check if first and final positions have been assigned
    if first_donor_pos == None: first_donor_pos = len(input[organism]['genomeseq'])
    if final_donor_pos == None: final_donor_pos = 0

    # do residual donor scanning 5' of first donor site
    splicesite_gff = { 'fref': FREF, 'fsource':'PSSM4FUNGI' }
    for site in scan_pssm_splice_site(input[organism]['genomeseq'][0:first_donor_pos]):
        gffs.append( site.togff( gff=dict(splicesite_gff) ) )
    # do residual acceptor scanning 3' of final donor site
    for site in scan_pssm_splice_site(input[organism]['genomeseq'][final_donor_pos:]):
        site.pos   = site.pos+final_donor_pos
        site.start = site.start+final_donor_pos
        site.end   = site.end+final_donor_pos
        gffs.append( site.togff( gff=dict(splicesite_gff) ) )

    # do residual donor scanning 5' of first donor site
    for site in scan_pssm_splice_site(input[organism]['genomeseq'][0:first_donor_pos],splicetype="acceptor"):
        gffs.append( site.togff( gff=dict(splicesite_gff) ) )
    # do residual acceptor scanning 3' of final donor site
    for site in scan_pssm_splice_site(input[organism]['genomeseq'][final_donor_pos:],splicetype="acceptor"):
        site.pos   = site.pos+final_donor_pos
        site.start = site.start+final_donor_pos
        site.end   = site.end+final_donor_pos
        gffs.append( site.togff( gff=dict(splicesite_gff) ) )

    ####################################################################
    if verbose: print stw.lap(), "donor & acceptor sites"
    ####################################################################

    for intron in known_introns:
        intron.assign_bp_and_ppts()
        _gff = {'fref':FREF,'fsource':'knownIntron','gclass':'KnownIntron'}
        knownintron_gfflines = intron.togff(gff=_gff,annotated_intron=True)
        # translate the actual intron track to 'IntronSpan'
        # Otherwise, it will pop up in the predicted IntronConnectingOrfs stack
        intron_gffline = list( knownintron_gfflines.pop(0) )
        intron_gffline[2] = 'IntronSpan'
        gffs.append( tuple(intron_gffline) )
        gffs.extend( knownintron_gfflines )

    ####################################################################
    if verbose: print stw.lap(), "known intron data"
    ####################################################################

    # check if gene structure tracks must be outputted
    if GFF_GENE_OUTPUT: gffs.extend( input[organism]['gff-gene'] )

    # check if unigene structure tracks must be outputted
    #if GFF_UNIGENE_OUTPUT:  gffs.extend( input[organism]['gff-unigene'] )

    if GFF_UNIGENE_OUTPUT:
        #gffs.extend( input[organism]['gff-unigene'] )
        gld = input[organism]['gldobj']
        main_ugId = gld.unigene_fref()
        main_ugId_fref = input[organism]['gff-unigene']
        for ugId in gld.get_unigene_frefs():
            gld.set_unigene_by_fref(ugId)
            # blank out the obtained structure -> filled list is a trigger
            # for returning, not novelly obtaining unigene gff data
            input[organism]['gff-unigene'] = []
            gffs.extend( gld._obtain_unigene_gff() )
        # reset `main` unigene ID
        gld.set_unigene_by_fref(main_ugId)
        input[organism]['gff-unigene'] = main_ugId_fref

    ####################################################################
    if verbose: print stw.lap(), "known gene & unigene(s)"
    ####################################################################

    for orf in input[organism]['orfs'].orfs:
        for sp in orf._signalp_sites:
            gffs.extend( sp.togff({'fref':FREF}) )

    ####################################################################
    if verbose: print stw.lap(), "locus SignalP"
    ####################################################################

    _gff = {'fref':FREF, 'fsource': 'tssPSSM4FUNGI', 'fmethod': 'tsspssm' }
    for tss in input[organism]['gldobj'].get_locus_tss_objects():
        # format tss.pssm_score as a string-formatted float
        tss.pssm_score = "%2.1f" % tss.pssm_score
        gffs.append( tss.togff( gff=_gff ) )

    ####################################################################
    if verbose: print stw.lap(), "locus TSS"
    ####################################################################

    gffs.extend( input[organism]['gldobj'].get_locus_tmhmm_gff() )

    ####################################################################
    if verbose: print stw.lap(), "predicted TMHMM domains of target"
    ####################################################################

    gffs.extend( input[organism]['gldobj'].get_locus_proteomics_gff() )

    ####################################################################
    if verbose: print stw.lap(), "proteomics data"
    ####################################################################

    for prjSignalP in create_projected_signalp_for_informants(
    OPTIONS.target,PCG,input,minimal_aa_overlap=0):
        gffs.extend( prjSignalP.togff(gff={'fref':FREF}) )

    ####################################################################
    if verbose: print stw.lap(), "SignalP projections"
    ####################################################################
    dna2protlength      = len(input[organism]['genomeseq'])/3
    array_algpresence   = PCG2codingarray(PCG,organism,dna2protlength)
    array_algsimilarity = PCG2similarityarray(PCG,organism,dna2protlength)

    max_algpresence = float(max(array_algpresence))
    max_algsimilarity = float(max(array_algsimilarity))

    zerotrack = []
    for i in range(0,dna2protlength):
        if array_algpresence[i] == 0:
            zerotrack.append(i)
            if i == dna2protlength-1 or array_algpresence[i+1] != 0:
                start = (zerotrack[0]*3)+1
                end   = (zerotrack[-1]*3)+3
                zerotrack = []
            else:
                continue
        else:
            start = (i*3)+1
            end   = (i*3)+3
        score = "%1.2f" % (float(array_algpresence[i])/max_algpresence)
        gcln  = "AlgPres %s" % FREF
        gffs.append( (FREF, "ABFGP","algPres",start,end,score,"+",".",gcln) )

    zerotrack = []
    for i in range(0,dna2protlength):
        if array_algsimilarity[i] == 0:
            zerotrack.append(i)
            if i == dna2protlength-1 or array_algsimilarity[i+1] != 0:
                start = (zerotrack[0]*3)+1
                end   = (zerotrack[-1]*3)+3
                zerotrack = []
            else:
                continue
        else:
            start = (i*3)+1
            end   = (i*3)+3
        score = "%1.2f" % (float(array_algsimilarity[i])/max_algsimilarity)
        gcln  = "AlgSim %s" % FREF
        gffs.append( (FREF, "ABFGP","algSim",start,end,score,"+",".",gcln) )

    ####################################################################
    if verbose: print stw.lap(), "Exon quantification"
    ####################################################################


    RIP_NT_WINDOW = 500
    RIP_NT_STEP   = 25
    RIP_FSOURCE   = 'RIP_%s_%s' % (RIP_NT_WINDOW,RIP_NT_STEP)
    RIP_FMETHOD   = 'RIPindex'
    for offset in range(0,len(input[organism]['genomeseq'])-RIP_NT_WINDOW,RIP_NT_STEP):
        center = offset+(RIP_NT_WINDOW/2)
        sta    = center - (RIP_NT_STEP/2)
        end    = center + (RIP_NT_STEP/2)
        seqslice = input[organism]['genomeseq'][offset:offset+RIP_NT_WINDOW]
        ripindex = RIPindex(seqslice)
        gff = ( FREF, RIP_FSOURCE, RIP_FMETHOD, sta, end, "%2.2f" % ripindex, '+', '.', 'RIPindex %s' % (FREF) )
        gffs.append( gff)


    ####################################################################
    if verbose: print stw.lap(), "RIPindex"
    ####################################################################


    # finally. The main results of the ABFGP programm. Is the annotated
    # gene structure correct?
    orfids_with_similarity = [ node[1] for node in PCG.get_organism_nodes(organism) ]
    print "PCG orfs:", orfids_with_similarity


    # get ABFGP predicted gene structure
    abfgp_exons, abfgp_introns = get_abfgp_predicted_genestructure(organism,PCG,
            introndata,known_exons,
            array_algpresence,array_algsimilarity,
            verbose=True)

    ####################################################################
    # check if there are NO exons at all returned => bailout!
    ####################################################################
    if not abfgp_exons:
        message="no convincing exons reported"
        if verbose: print stw.lap(), message, "=> bailout" 
        # create a bailout message file
        abgpsysexit(input,OPTIONS,message=message)

    print "ABFGP exons"
    is_first_exon = True
    for exon in abfgp_exons:
        exon._gff['gname'] = FREF

        # get algpresence/algsimilarity scores for this exon
        exon_aa_start = exon.start/3
        exon_aa_stop  = (exon.end/3) + 1
        if exon_aa_stop == exon_aa_start:
            # avoid ZeroDivisionError
            exon_aa_stop+=1
        algprs_score  = sum(array_algpresence[exon_aa_start:exon_aa_stop])
        algsim_score  = sum(array_algsimilarity[exon_aa_start:exon_aa_stop])
        algabs_score  = list(array_algpresence[exon_aa_start:exon_aa_stop]).count(0)
        algprs_ratio  = float(algprs_score) / (max_algpresence*(exon_aa_stop-exon_aa_start))
        algsim_ratio  = float(algsim_score) / (max_algpresence*(exon_aa_stop-exon_aa_start))
        if exon_aa_stop-exon_aa_start == 0:
            algabs_ratio = 0.0
        else:
            algabs_ratio = float(algabs_score) / (exon_aa_stop-exon_aa_start)

        exonNode      = (organism,exon.orf.id)
        if exonNode in PCG.get_organism_nodes(organism):
            informant_set = Set()
            for (nodeQ,nodeS) in PCG.weights.keys():
                if nodeQ == exonNode:
                    informant_set.add( PCG.organism_by_node(nodeS) )
            cnt_informants = len(informant_set)
        else:
            cnt_informants = 0

        ########################################################################
        # check if exon is `confirmed` by algprs_ratio / algsim_ratio /algabs_ratio
        ########################################################################
        if exon.donor.__class__.__name__ == 'SequenceErrorCoordinate' and\
        exon.acceptor.__class__.__name__ == 'SequenceErrorCoordinate':
            # squeezed in between sequence errors -> Trust this one!
            exon.pssm_score = 1
            exon._gff['fscore'] = 1
        elif algprs_ratio == 0.0:
            # already assigned as a non-confirmed exon
            print "KNOWN low exon support:", algprs_ratio, algsim_ratio
            pass
        elif algprs_ratio < 0.30 and algsim_ratio/algprs_ratio < 0.66:
            print "low exon support:", algprs_ratio, algsim_ratio
            exon.pssm_score = 0
            exon._gff['fscore'] = 0
        elif algprs_ratio < 0.15:
            print "low exon support:", algprs_ratio, algsim_ratio
            exon.pssm_score = 0
            exon._gff['fscore'] = 0
        elif algabs_ratio >= 0.20:
            print "low exon support:", algprs_ratio, algsim_ratio
            exon.pssm_score = 0
            exon._gff['fscore'] = 0
        else:
            pass 

        ########################################################################
        # label FirstExon as uncertain when not starting with a Methionine
        ########################################################################
        if is_first_exon and exon.proteinsequence() and\
        exon.proteinsequence()[0].upper() != "M":
            print "not starting with a Methionine!", exon, exon.proteinsequence()
            exon.pssm_score = 0
            exon._gff['fscore'] = 0
            
        # set label to to non-first exon
        is_first_exon = False
            
        # get gffdata list of gff tuple
        gffdata = list(exon.togff())

        col9tv = [
            'Informants %s' % cnt_informants,
            'AlignmentPresence %s' % algprs_score,
            'AlignmentPresenceRatio %1.2f' % algprs_ratio,
            'AlignmentAbsence %s' % algabs_score,
            'AlignmentAbsenceRatio %1.2f' % algabs_ratio,
            'AlignmentSimilarity %s' % algsim_score,
            'AlignmentSimilarityRatio %1.2f' % algsim_ratio,
            'ProteinSequene %s'  % exon.proteinsequence(), 
            ]
        gffdata[8] = gffdata[8] + "; " + "; ".join(col9tv)

        ########################################################################
        if verbose: print exon,exon.pssm_score,exon._gff['fscore'],"\n", gffdata
        ########################################################################
        # store to GFFdata
        gffs.append(tuple(gffdata))

        

    print "ABFGP introns"
    for pos in range(0,len(abfgp_introns)):
        intron = abfgp_introns[pos]
        if not intron:
            # all introns should be introns/seqerrors/connections; not None
            # if None -> create/label as incertain intron
            prevexongff = abfgp_exons[pos].togff()
            nextexongff = abfgp_exons[pos+1].togff()
            score       = 0
            gffdata = (FREF,'abfgp','AbfgpIntron',prevexongff[4]+1,prevexongff[3]-1,score,'+','.',"AbfgpGeneStructure %s" % FREF)
            gffs.append(tuple(gffdata))
            continue

        if intron.__class__.__name__ == "SequenceErrorConnectingOrfs":
            # HARD-LABEL the sequenceerrors as reliable
            intron.pssm_score = 1
            intron._gff['fscore'] = 1

        # Check how reliable this intron is.
        # In introndata some very weak introns are present too!
        if intron.coords() in introndata.keys():
            label = (intron.orfDonor.id,intron.orfAcceptor.id)
            if assessed_interfaces[OPTIONS.target].has_key(label):
                observed_cnt  = float(len(introndata[intron.coords()]))
                interface_cnt = float(len(assessed_interfaces[OPTIONS.target][label]))
                if observed_cnt / interface_cnt < 0.33:
                    # mark intron as NOT confirmed
                    intron.pssm_score = 0
            else:
                # THIS CAN HAPPEN! In case of RELABELLING of introns
                # in case of tinyexons/ cig introns. Not very elegantly
                # solved, but well, that's life ;-(
                pass

        # update to properties
        intron._gff['gname'] = FREF
        gffdata = list(intron.togff({'gname':FREF}))

        intron_aa_start = intron.start/3
        intron_aa_stop  = (intron.end/3) + 1
        algprs_score  = sum(array_algpresence[intron_aa_start:intron_aa_stop])
        algsim_score  = sum(array_algsimilarity[intron_aa_start:intron_aa_stop])
        algprs_ratio  = float(algprs_score) / (max_algpresence*(intron_aa_stop-intron_aa_start))
        algsim_ratio  = float(algsim_score) / (max_algpresence*(intron_aa_stop-intron_aa_start))


        # add individual intron projection data elements
        if intron.coords() in introndata.keys():
            pass
            #assessed_interfaces


        col9tv = [
            'AlignmentPresence %s' % algprs_score,
            'AlignmentPresenceRatio %1.2f' % algprs_ratio,
            'AlignmentSimilarity %s' % algsim_score,
            'AlignmentSimilarityRatio %1.2f' % algsim_ratio,
            ]
        gffdata[8] = gffdata[8] + "; " + "; ".join(col9tv)


        ########################################################################
        if verbose:
            print intron, "SHARED_NT/AA:",
            print (intron.shared_nts,intron.shared_aa), "\n", gffdata
        ########################################################################
        # store to GFFdata
        gffs.append(tuple(gffdata))

    ####################################################################
    if verbose: print stw.lap(), "ABFGP Gene Model"
    ####################################################################


    for key,intronlist in introndata.iteritems():
        for intron in intronlist:
            # organism/gene identifier is stored in intron._referce attribute)
            orgfullname = _get_organism_full_name(intron._reference,truncate=True)
            protfref    = input[intron._reference]['proteinfref']

            # Check if this is an unigene identifier.
            # If so, perform a different layout of orgfullname
            if input.has_key(intron._reference) and\
            input[intron._reference].has_key('is_unigene') and\
            input[intron._reference]['is_unigene'] == True:
                # get rid of the uigeneid[org] organism tag
                protfref = input[intron._reference]['proteinfref'].split("[")[0]     
                org_key  = input[intron._reference]['is_unigene_of']
                orgfullname = _get_organism_full_name(org_key,truncate=True)

            # update supporting intron GFF data
            intron._gff['FREF'] = FREF
            intron._gff['gclass'] = 'intron'
            intron._gff['gname'] = "%s_%s_%s" % (intron._reference,intron.start,intron.end)
            intron._gff['column9data'].update( {
                    'distance': intron._distance,
                    'entropy':  "(%1.2f,%1.2f)" % (intron._apps_donor,intron._apps_accep),
                    'Note':     "%s [%s]" % ( orgfullname, protfref ),
                    } )

            # append gff output of supporting intron to ggfs list
            gffs.append( intron.togff() )

    ####################################################################
    if verbose: print stw.lap(), "predicted intron data"
    ####################################################################

    ERRONEOUSLY_ANNOTATED_EXON_GFF = {
            'fsource':  'abfgp',
            'fmethod':  'ErroneousExon',
            'gname'  :  FREF,
            'gclass' :  'AnnotatedGeneModelError', }
    ERRONEOUSLY_ANNOTATED_INTRON_GFF= {
            'fsource':  'abfgp',
            'fmethod':  'ErroneousIntron',
            'gname'  :  FREF,
            'gclass' :  'AnnotatedGeneModelError', }

    errors = []
        
    for exon in known_exons:
        if exon.coords() not in [ e.coords() for e in abfgp_exons ]:
            erroneous = exon.deepcopy()
            erroneous._gff.update(ERRONEOUSLY_ANNOTATED_EXON_GFF)
            errors.append( erroneous )

    for intron in known_introns:
        abfgp_intron_coords = []
        for i in abfgp_introns:
            # TODO: here all introns should be introns, not None ;-)
            if i:
                abfgp_intron_coords.append( i.coords() )
        if intron.coords() not in abfgp_intron_coords:
            erroneous = intron.deepcopy()
            erroneous._gff.update(ERRONEOUSLY_ANNOTATED_INTRON_GFF)
            errors.append( erroneous )

    if len(errors) > 1:
        # AnnotatedGeneModelError layout class joins all tracks
        # in a single line. 2 unrelated errors are joined/merged by an
        # intron-like representation. To avoid this, ungroup unrelated
        # annotation errors by adding an index to the gname attribute
        errors = [ (error.start,error) for error in errors ]
        errors.sort()
        errors = [ error for (start,error) in errors ]
        index = 0
        errors[0]._gff['gname'] = errors[0]._gff['gname']+"-%s" % index
        for pos in range(1,len(errors)):
            if errors[pos].start != errors[pos-1].end:
                index += 1
            errors[pos]._gff['gname'] = errors[pos]._gff['gname']+"-%s" % index

    # append errors to gffs
    for error in errors: gffs.append( error.togff() )

    ####################################################################
    if verbose: print stw.lap(), "Annotated GeneModel errors"
    ####################################################################

    if OPTIONS.unigenebenchmark:
        # do performance measurement on (full-length) unigene
        abgp_genestructure_coordinate_set = Set()
        for exon in abfgp_exons:
            # update coordinate Set; ignore first coord (python->GGB coords)
            abgp_genestructure_coordinate_set.update( exon.dna_range()[1:] )
        # add final stop codon for unigene benchmarking
        stop = range( abfgp_exons[-1].end+1, abfgp_exons[-1].end+3+1 )
        abgp_genestructure_coordinate_set.update(stop)
        # add SequenceError `introns` to the coordinate set as well
        for intron in abfgp_introns:
            if intron.__class__.__name__ == "SequenceErrorConnectingOrfs":
                abgp_genestructure_coordinate_set.update( intron.dna_range()[1:] )

        # get coordinates of the unigene to benchmark against
        unigene_coordinate_set = gffs2coordset(
                    input[OPTIONS.target]['gff-unigene'],
                    fmethod=[GFF_UGEXON_FMETHOD]
                    )

        performance_gff_unigene = abgp_performance_measurement_of_sets(
                abgp_genestructure_coordinate_set,
                unigene_coordinate_set,
                gffdata = {
                    'gname':   '%s-unigene' % FREF,
                    'fref':     FREF,  }  )
        # extend to the gffs list
        gffs.extend( performance_gff_unigene )


        print abfgp_exons[-1], abfgp_exons[-1].start, abfgp_exons[-1].end
        for tr in performance_gff_unigene:
            print tr

    ####################################################################
    if verbose: print stw.lap(), "unigene performance"
    ####################################################################

    # TODO: replace FREF in known gff lines...
    # TODO: only in annotated gene structure & unigene alignment, a wrong `fref` is used.
    while [] in gffs: gffs.remove([])
    while () in gffs: gffs.remove(())
    for i in range(0,len(gffs)):
        try:
            line = list(gffs[i])
            line[0] = FREF
            gffs[i] = tuple(line)
        except:
            print "# WARNING: unexpected gff line: ", line

    # add the SEQUENCE track to the gffs
    seqlen = len(input[organism]['genomeseq'])
    gff = ( FREF, METHOD, 'Sequence', 1, seqlen, '.', '+', '.', 'Sequence %s ; Organism %s' % (FREF, organism) )
    gffs.insert( 0, gff)

    ####################################################################
    if verbose: print stw.lap(), "metadata tracks"
    ####################################################################


    protfname    = "%s.protein.%sSIL.fa" % (FREF,len(input)-1)
    protfullpath = os.path.join(OPTIONS.outdir,protfname)
    protseq = ""
    for pos in range(0,len(abfgp_introns)):
        protseq+=abfgp_exons[pos].proteinsequence()
        if abfgp_introns[pos]:
            protseq+=abfgp_introns[pos].shared_aa
        else:
            # hmmm... big problem. Stil a NONE intron here !?
            # TODO: check why!!!
            print "STILL A None-Intron in the GeneStructure"
            # append XXXX on this position
            protseq+="XXXX"
    protseq+=abfgp_exons[-1].proteinsequence()

    fh = open(protfullpath, 'w')
    fh.write(">%s.abfgp\n" % FREF)
    fh.write(protseq)
    fh.close()

    ####################################################################
    if verbose: print stw.lap(), "protein single fasta written"
    ####################################################################


    # Ready with gathering gff data! Create file
    fname = "%s.intronmapping.%sSIL.gff" % (OPTIONS.target,len(input)-1)

    # make fullpath of fname
    fullpath = os.path.join(OPTIONS.outdir,fname)
    # check if file exists and overwriting is aloued
    if not OPTIONS.force:
        fullpath = _create_unique_filename_with_integer_suffix(fullpath)

    # write to GFF file
    fh = open(fullpath, 'w')
    fh.write(gffs2txt(gffs))
    fh.close()

    ####################################################################
    if verbose: print stw.lap(), "result gff file: %s " % (fullpath)
    ####################################################################

    # done! return the freshly created gff file name
    return fullpath 

# end of function create_intron_mapping_gff_file


def get_exons_of_inwpcbg(organism,inwpCBG,exons,verbose=True):
    """ """
    orfObj = inwpCBG.get_orfs_of_graph(organism=organism)[0]
    minsr = inwpCBG.minimal_spanning_range(organism=organism)
    maxsr = inwpCBG.maximal_spanning_range(organism=organism)
    ############################################################################
    #if verbose:
    #    print "## inwpCBG MINSR/MAXSR:", ( min(minsr), max(minsr) ),
    #    print ( min(maxsr), max(maxsr) ), inwpCBG
    #    for exon in exons: print "current annotated exon:", exon
    #    print "## inwpCBG"
    ############################################################################
    matched_exons = []
    for exon in exons:
        if exon.orf.id == orfObj.id:
            if len(minsr)>=2 and not minsr.intersection(exon.protein_range()):
                # possibly undetected stopless3n intron; CFU_829297
                ################################################################
                if verbose:
                    print exon, "FAILED A",
                    print min(exon.protein_range()),"..",
                    print max(exon.protein_range())
                ################################################################
                pass

            elif maxsr.intersection(exon.protein_range()):
                matched_exons.append(exon)
                ################################################################
                if verbose:
                    overlap = minsr.intersection(exon.protein_range())
                    print exon,len(minsr), len(overlap), "protrange:",
                    print (min(exon.protein_range()),max(exon.protein_range()))
                ################################################################
            else:
                ################################################################
                if verbose: print exon, "FAILED B"
                ################################################################
                pass
        else:
            pass
    return matched_exons

# end of function get_exons_of_inwpcbg


def _list2dictcounts(ll):
    """ """
    return dict([(item,ll.count(item)) for item in Set(ll)])
# end of function _list2dictcounts


def get_raw_abfgp_genestruture(organism,PCG,introndata,known_exons,verbose=True):
    """ """
    # get InwardsPointingCodingBlockGraphs from PCG
    #inwpcbgs = PCG2inwpCBGS(PCG,omit_singletons=True)
    inwpcbgs = PCG2inwpCBGS(PCG,omit_singletons=False)

    max_inwpCBG_node_count = len(PCG.organism_set())

    # threshold for allowing novel exons into the ABFGP GeneStructure
    NOVEL_ORF_EXON_SUPPORT_RATIO_THRESHOLD = 0.4

    # return list of ABFGP predicted exons
    abfgp_exons = []

    
    # find all sequence error interfaces in observed_introns
    seqerror_interfaces = {} 
    for key,vlist in introndata.iteritems():
        seObj = vlist[0]
        if seObj.__class__.__name__ == "SequenceErrorConnectingOrfs":
            interface = _get_main_interface([ se._label for se in vlist ])
            seqerror_interfaces[key] = interface

    # start looping over the inwpcbgs
    for pos in range(0,len(inwpcbgs)):
        inwpCBG = inwpcbgs[pos]
        exons = get_exons_of_inwpcbg(organism,inwpCBG,known_exons,verbose=verbose)
        #########################################################################
        #if verbose:
        #    print inwpCBG
        #    for exon in exons:
        #        minsr = inwpCBG.overall_minimal_spanning_range(organism=organism)
        #        overlap = minsr.intersection(exon.protein_range())
        #        print exon, len(minsr), len(overlap)
        #########################################################################
        if not exons:
            # build this exon. Get orf, donor & acceptor object
            orfObj = inwpCBG.get_orfs_of_graph(organism=organism)[0]
            ####################################################################
            if verbose: print "building exon on:", orfObj
            ####################################################################
            # get donor object
            donorObj = None
            donorObjSupportCnt = 0
            donorObjSupportTypes = {}
            for nextpos in range(pos+1,len(inwpcbgs)):
                next = inwpcbgs[nextpos]
                nextorfObj = next.get_orfs_of_graph(organism=organism)[0]
                for intron_key in introndata.keys():
                    for thisintron in introndata[intron_key]:
                        if thisintron._label == (orfObj.id,nextorfObj.id):
                            donorObj = thisintron.donor
                            donorObjSupportCnt   = len(introndata[intron_key])
                            ####donorObjSupportTypes = _list2dictcounts([i._gff['fsource'] for i in introndata[intron_key]])
                            break
                    if donorObj: break
                if donorObj: break
            # get acceptor object
            accepObj = None
            accepObjSupportCnt = 0
            accepObjSupportTypes = {}
            for prevpos in range(pos-1,-1,-1):
                prev = inwpcbgs[prevpos]
                prevorfObj = prev.get_orfs_of_graph(organism=organism)[0]
                for intron_key in introndata.keys():
                    for thisintron in introndata[intron_key]:
                        if thisintron._label == (prevorfObj.id,orfObj.id):
                            accepObj = thisintron.acceptor
                            accepObjSupportCnt = len(introndata[intron_key])
                            ####accepObjSupportTypes = _list2dictcounts([i._gff['fsource'] for i in introndata[intron_key]])
                            break
                    if accepObj: break
                if accepObj: break

            # check if NOT donorObj AND NOT accepObj -> in case of tinyexons,
            # CIG exons etc the interface is not found!
            # this is only possible when this happens in the middle of the
            # gene structure, not at the ends. Therefor, positional check.
            if (not donorObj and pos<len(inwpcbgs)-1)and (not accepObj and pos>0):
                # what was the interface if this was tinyexon/cig case
                nextorfObj = inwpcbgs[pos+1].get_orfs_of_graph(organism=organism)[0]
                prevorfObj = inwpcbgs[pos-1].get_orfs_of_graph(organism=organism)[0]
                for intron_key in introndata.keys():
                    for thisintron in introndata[intron_key]:
                        if thisintron._label == (prevorfObj.id,nextorfObj.id):
                            if thisintron.orfAcceptor.id == orfObj.id:
                                accepObj = thisintron.acceptor
                                break
                for intron_key in introndata.keys():
                    for thisintron in introndata[intron_key]:
                        if thisintron._label == (prevorfObj.id,nextorfObj.id):
                            if thisintron.orfDonor.id == orfObj.id:
                                donorObj = thisintron.donor
                                break
                # check if it was succesfull
                if not donorObj or not accepObj:
                    # nope -> both should be mathed, otherwise None
                    donorObj = None
                    accepObj = None
                else:
                    # okay, found. But: RELABEL the introns otherwise
                    # next steps will crash as well!
                    # TODO: RELABELLING should be performed on
                    # TODO: assesed_interfaces too! Now, this label can not
                    # TODO: be found in this data struct....
                    for intron_key in introndata.keys():
                        for thisintron in introndata[intron_key]:
                            if thisintron._label == (prevorfObj.id,nextorfObj.id):
                                if thisintron.orfAcceptor.id == orfObj.id:
                                    thisintron._label = (prevorfObj.id,orfObj.id)
                                if thisintron.orfDonor.id == orfObj.id:
                                    thisintron._label = (orfObj.id,nextorfObj.id)

            ####################################################################
            if verbose: print "donor::", donorObj, "acceptor::", accepObj
            ####################################################################

            if donorObj.__class__.__name__ == 'SequenceErrorCoordinate':
                # TODO: make this more generic
                # temporarily, set pssm_score attribute as float (not None)
                donorObj.pssm_score    = 1.0
            if accepObj.__class__.__name__ == 'SequenceErrorCoordinate':
                # TODO: make this more generic
                # temporarily, set pssm_score attribute as float (not None)
                accepObj.pssm_score    = 1.0

            donorObjSupportRatio = float(donorObjSupportCnt) / (inwpCBG.node_count()-1)
            accepObjSupportRatio = float(accepObjSupportCnt) / (inwpCBG.node_count()-1)


            if (donorObj and donorObjSupportRatio >= NOVEL_ORF_EXON_SUPPORT_RATIO_THRESHOLD\
            and accepObj and accepObjSupportRatio >= NOVEL_ORF_EXON_SUPPORT_RATIO_THRESHOLD) or\
            (donorObj.__class__.__name__ == 'SequenceErrorCoordinate' and\
            accepObj.__class__.__name__ == 'SequenceErrorCoordinate'):
                if donorObj.pos > accepObj.pos:
                    exons = [ ExonOnOrf(accepObj,donorObj,orfObj) ]
                    ############################################################
                    if verbose: print "CREATED::", exons[0]
                    ############################################################
                elif donorObj.__class__.__name__ == 'SequenceErrorCoordinate' and\
                accepObj.__class__.__name__ == 'SequenceErrorCoordinate':
                    # insertion/deletion stopcodon
                    exons = [ ExonOnOrf(accepObj,donorObj,orfObj) ]
                    ############################################################
                    if verbose: print "CREATED::", exons[0]
                    ############################################################
                else:
                    # do not allow this Exon creation
                    print "WARNING:: donorObj.pos <= accepObj.pos", donorObj.pos, accepObj.pos
            elif (donorObj and donorObjSupportRatio >= NOVEL_ORF_EXON_SUPPORT_RATIO_THRESHOLD) or\
            donorObj.__class__.__name__ == 'SequenceErrorCoordinate':
                # donor only -> create FirstExonOnOrf (and hope for the best)
                # do this only when there is not a single exon in the list
                inwpCBGmaxsr = inwpCBG.maximal_spanning_range(organism=organism)
                inwpCBGminsr = inwpCBG.minimal_spanning_range(organism=organism)
                if not abfgp_exons:
                    atg = _propose_atg_on_orf(orfObj,inwpCBG)
                    if atg == False:
                        checks = [  False, None, None ]
                    elif type(atg) == type(int()):
                        # atg is not a TSS object but an integer.
                        # hmm....for now, leave as it is
                        checks = [  atg != False,
                                    donorObj.pos > atg and atg/3,
                                    # ATG in front or in the first 33% percent of the MAXSR
                                    atg/3 <= ( min(inwpCBGmaxsr) + int(0.33*(max(inwpCBGmaxsr)-min(inwpCBGmaxsr))) ) ]

                    else:
                        checks = [  atg != False,
                                    donorObj.pos > atg.pos and atg.pos/3,
                                    # ATG in front or in the first 33% percent of the MAXSR
                                    atg.pos/3 <= ( min(inwpCBGmaxsr) + int(0.33*(max(inwpCBGmaxsr)-min(inwpCBGmaxsr))) ) ]
                    if not False in checks:
                        exons = [ FirstExonOnOrf(atg,donorObj,orfObj) ]
                        exons[0].DO_NOT_CHANGE_ATG = True
                        ########################################################
                        if verbose: print "NOT A -> FirstExon", inwpCBG,exons[0]
                        ########################################################
                    else:
                        if checks[0] == False:
                            print "WARNING:: NO FIRST EXON CAN BE CREATED, NO ATG"
                        if checks[1] == False:
                            print "WARNING:: donorObj.pos <= atg.pos", donorObj.pos, atg
                        if checks[2] == False:
                            print "WARNING:: FirstExon atg.pos > min(maxsr)", atg, (min(inwpCBGmaxsr)*3,max(inwpCBGmaxsr)*3)
                            
                        # have a closer look at this inwpCBG. When it is REALLY
                        # convincing -> store as partial Exon. Most likely an undetected
                        # point mutation / sequence error is present
                        if inwpCBG.count_orfs_labeled_as_annotated_exon() and\
                        float(inwpCBG.count_orfs_labeled_as_first_exon()) / float(inwpCBG.count_orfs_labeled_as_annotated_exon()) >= 0.50 and\
                        float(inwpCBG.count_orfs_labeled_as_first_exon()) / float(inwpCBG.count_genomic_informants()) >= 0.50 and\
                        ( inwpCBG.get_bitscore() > min([inwp.get_bitscore() for inwp in inwpcbgs]) or\
                        inwpCBG.get_identityscore() > min([inwp.get_identityscore() for inwp in inwpcbgs]) ):
                            # yes, store it. Define an artificial start position in between MINSR and MAXSR start
                            dummyatg = min(inwpCBGmaxsr)*3 + orfObj.frame
                            print "CREATING FirstExonoOnOrf with DUMMYATG:", dummyatg
                            exons = [ FirstExonOnOrf(dummyatg,donorObj,orfObj) ]
                            exons[0].DO_NOT_CHANGE_ATG = True
                        elif inwpCBG.count_orfs_labeled_as_annotated_exon() and\
                        inwpCBG.count_orfs_labeled_as_annotated_exon() / float(inwpCBG.count_genomic_informants()) >= 0.50 and\
                        len(inwpCBGminsr) >= 15 and\
                        ( inwpCBG.get_bitscore() > min([inwp.get_bitscore() for inwp in inwpcbgs]) or\
                        inwpCBG.get_identityscore() > min([inwp.get_identityscore() for inwp in inwpcbgs]) ):
                            dummyatg = min(inwpCBGmaxsr)*3 + orfObj.frame
                            print "CREATING FirstExonoOnOrf with DUMMYATG:", dummyatg
                            exons = [ FirstExonOnOrf(dummyatg,donorObj,orfObj) ]
                            # label ATG as unchangable ...
                            exons[0].DO_NOT_CHANGE_ATG = True

                elif donorObj.__class__.__name__ == 'SequenceErrorCoordinate' and\
                abfgp_exons and abfgp_exons[-1].orf.id == orfObj.id:
                    # assume this to be an extention of an already listed exon
                    print "YAHOOOOOOOOO!!!!!!!!!!!!"
                    cur_exon = abfgp_exons[-1]
                    if cur_exon.__class__.__name__ in ['SingleExonOnOrf','FirstExonOnOrf']:
                        new_exon =  FirstExonOnOrf(cur_exon.acceptor,donorObj,orfObj)
                    else:
                        new_exon =  ExonOnOrf(cur_exon.acceptor,donorObj,orfObj)
                    abfgp_exons[-1] = new_exon
                    print "CREATION EXISTING EXON EXTENTION:", new_exon
                    exons = []

                elif donorObj.__class__.__name__ == 'SequenceErrorCoordinate':
                    # SequenceError label -> enforce in the gene structure!
                    exonStartObj = ExonStart( min(inwpCBGmaxsr)*3 + orfObj.frame )
                    exons = [ ExonOnOrf(exonStartObj,donorObj,orfObj) ]
                else:
                    # Something that should become a first exon, but there
                    # is already a first exon. Omit it here!
                    # TODO: this is tricky! When by chance a False FirstExon
                    # is already created, a True first exon could be omitted.
                    # But, vast majority is true FP exon omittings here.
                    pass

            elif (accepObj and accepObjSupportRatio >= NOVEL_ORF_EXON_SUPPORT_RATIO_THRESHOLD) or\
            accepObj.__class__.__name__ == 'SequenceErrorCoordinate':
                # acceptor only -> create FinalExonOnOrf (and hope for the best)
                exons = [ FinalExonOnOrf(accepObj,orfObj.end,orfObj) ]
                ################################################################
                if verbose:
                    print "NOT D -> FinalExon", inwpCBG, exons[0],
                    print orfObj._has_donor_sites_predicted,
                    print len(orfObj._donor_sites)
                ################################################################

            elif accepObj:
                # check if this is for other reasons a serious FinalExon candidate
                is_added_as_exon = False
                if pos > 0:
                    prevInwpCBG = inwpcbgs[pos-1]
                    checkA = inwpCBG.get_projected_tailing_stop_aa_difference() <=\
                             prevInwpCBG.get_projected_tailing_stop_aa_difference()
                    checkB = inwpCBG.get_identityscore() >= prevInwpCBG.get_identityscore()
                    print [ checkA, checkB ], prevInwpCBG
                    if [ checkA, checkB ] == [ True, True ]:
                        is_added_as_exon = True
                        exons = [ FinalExonOnOrf(accepObj,orfObj.end,orfObj) ]
                        ########################################################
                        if verbose:
                            print "NOT D -> FinalExon", inwpCBG, exons[0],
                            print orfObj._has_donor_sites_predicted,
                            print len(orfObj._donor_sites)
                        ########################################################
                   
                if not is_added_as_exon:
                    print "A but intron threshold low:",
                    print accepObjSupportRatio,"<",
                    print NOVEL_ORF_EXON_SUPPORT_RATIO_THRESHOLD
                    pass

            elif donorObj:
                print "D but intron threshold low:",
                print donorObjSupportRatio,"<",
                print NOVEL_ORF_EXON_SUPPORT_RATIO_THRESHOLD
                pass

            elif not donorObj and not accepObj and not exons and\
            len(Set([ _inwpCBG.get_orfs_of_graph(organism=organism)[0].id for _inwpCBG in inwpcbgs ])) == 1 and\
            sum( [ len(_inwpCBG.overall_minimal_spanning_range(organism=organism)) for _inwpCBG in inwpcbgs ] ) >= 50:
                # no donor or acceptor object BUT a fairly proper SingleExonOnOrf thingy
                atg = _propose_atg_on_orf(orfObj,inwpCBG)
                if atg:
                    exons = [ SingleExonOnOrf(atg,orfObj.end,orfObj) ]
                    print "CREATING SingleExonoOnOrf with ATG:", atg
                else:
                    firstInwpCBGmaxsr = inwpcbgs[0].maximal_spanning_range(organism=organism)
                    dummyatg = min(firstInwpCBGmaxsr)*3 + orfObj.frame
                    print "CREATING SingleExonoOnOrf with DUMMYATG:", dummyatg
                    exons = [ SingleExonOnOrf(dummyatg,orfObj.end,orfObj) ]
                    # label ATG as unchangable ...
                    exons[0].DO_NOT_CHANGE_ATG = True
            else:
                print "NOT A & D::", inwpCBG, accepObj, donorObj, orfObj
                pass

                
        if len(exons) == 1:
            # try if we have to change exon end(s) with known SequenceError coordinates
            thisexon = exons[0]
            orfObj   = thisexon.orf
            if orfObj.id in [ itf[0] for itf in seqerror_interfaces.values() ] and\
            thisexon.donor.__class__.__name__ != 'SequenceErrorCoordinate':
                # labelled as a sequence eror interface!
                # get next inwpCBG object and its Orf
                nextorfObj = None
                for nextpos in range(pos+1,len(inwpcbgs)):
                    nextInwpCBG = inwpcbgs[nextpos]
                    nextorfObj  = nextInwpCBG.get_orfs_of_graph(organism=organism)[0]
                    break
                if nextorfObj and (orfObj.id,nextorfObj.id ) in seqerror_interfaces.values():
                    # get SequenceErrorCoord object from SequenceErrorConnectingOrfs object
                    seqErrorCoordObj = None
                    for k,interface in seqerror_interfaces.iteritems():
                        if interface == (orfObj.id,nextorfObj.id ):
                            seqErrorCoordObj = introndata[k][0].donor
                            seqErrorCoordObj.pssm_score = 1.0
                            break
                    # rebuild this exon
                    print "##### SequenceError REBUILDING:", thisexon, 'after',(orfObj.id,nextorfObj.id )
                    if thisexon.__class__.__name__ in ['FirstExonOnOrf','SingleExonOnOrf']:
                        exons = [ FirstExonOnOrf(thisexon.acceptor,seqErrorCoordObj,orfObj) ]
                    else:
                        exons = [ ExonOnOrf(thisexon.acceptor,seqErrorCoordObj,orfObj) ]

            thisexon = exons[0]
            orfObj   = thisexon.orf
            if orfObj.id in [ itf[1] for itf in seqerror_interfaces.values() ] and\
            thisexon.acceptor.__class__.__name__ != 'SequenceErrorCoordinate':
                # labelled as a sequence eror interface!
                # get prev inwpCBG object and its Orf
                prevorfObj = None
                for prevpos in range(pos-1,-1,-1):
                    prevInwpCBG = inwpcbgs[prevpos]
                    prevorfObj  = prevInwpCBG.get_orfs_of_graph(organism=organism)[0]
                    break
                if prevorfObj and ( prevorfObj.id, orfObj.id ) in seqerror_interfaces.values():
                    # get SequenceErrorCoord object from SequenceErrorConnectingOrfs object
                    seqErrorCoordObj = None
                    for k,interface in seqerror_interfaces.iteritems():
                        if interface == ( prevorfObj.id, orfObj.id ):
                            seqErrorCoordObj = introndata[k][0].acceptor
                            seqErrorCoordObj.pssm_score = 1.0
                            break
                    # rebuild this exon
                    print "##### SequenceError REBUILDING:", thisexon, 'before', (orfObj.id,prevorfObj.id )
                    if thisexon.__class__.__name__ in ['FinalExonOnOrf','SingleExonOnOrf']:
                        exons = [ FinalExonOnOrf(seqErrorCoordObj,thisexon.donor,orfObj) ]
                    else:
                        exons = [ ExonOnOrf(seqErrorCoordObj,thisexon.donor,orfObj) ]


        # check if exons(s) already in (ordered) list of abfgp exons
        for exon in exons:
            if not exon:
                continue
            elif exon.coords() in [ e.coords() for e in abfgp_exons]:
                continue
            elif exon.orf.id not in [ e.orf.id for e in abfgp_exons]:
                # in case this is the FIRST ADDED exon:
                # check if it starts with a methionine and if not,
                # try to create it
                if not abfgp_exons:
                    if not exon.proteinsequence():
                        pass # exon of 2 nucleotides!
                    elif hasattr(exon,"DO_NOT_CHANGE_ATG") and getattr(exon,"DO_NOT_CHANGE_ATG") == True:
                        pass # do not change ATG -> a (nonsense) ATG is already chosen!
                    else:
                        atg = _propose_atg_on_orf(exon.orf,inwpCBG)
                        if atg:
                            if hasattr(atg,"pos") and hasattr(exon.acceptor,"pos") and atg.pos == exon.acceptor.pos:
                                # leave as is; TSS/Methionine and unchanged
                                pass
                            elif type(atg) == type(int()) or (hasattr(atg,"pos") and hasattr(exon.acceptor,"pos") and atg.pos != exon.acceptor.pos):
                                print "FIRST EXON START ELEMENT CHANGED:", exon.acceptor,"->",atg
                                if exon.donor.__class__.__name__ == 'StopCodon':
                                    exon = SingleExonOnOrf(atg,exon.donor,exon.orf)
                                else:
                                    exon = FirstExonOnOrf(atg,exon.donor,exon.orf)
                            else:
                                # leave as is; TSS/Methionine and unchanged
                                pass
                        else:
                            # just append.... this will potentially result in an error lateron
                            # but, the chance that this occurs is very very small!
                            pass
                    
                # append to list
                abfgp_exons.append( exon.deepcopy() )

            elif exon.orf.id in [ e.orf.id for e in abfgp_exons] and 'SequenceErrorCoordinate' in\
            [ exon.donor.__class__.__name__, exon.acceptor.__class__.__name__ ]:
                # trust SequenceErrors; append to list
                abfgp_exons.append( exon.deepcopy() )

            elif exon.orf.id in [ e.orf.id for e in abfgp_exons] and len(abfgp_exons)==1 and\
            abfgp_exons[0].__class__.__name__ == 'SingleExonOnOrf' and not introndata:
                # decided that this will become a SingleExonOnOrf ->
                # do not overlay something here
                continue
                
            else:
                # final check: do we not create overlap for exons with:
                # - identical Orfs with stopless3n introns?
                # - FirstExon start site ambiguity?
                print "OVERLAP stopless3n or start site ambiguity!?::", exon, [ e.coords() for e in abfgp_exons]

                exonrange = Set(exon.dna_range())
                for pos in range(0,len(abfgp_exons)):
                    e = abfgp_exons[pos]
                    if e.coords()[1] == exon.coords()[1] and\
                    e.orf.id == exon.orf.id:
                        # start site ambiguity
                        if e.acceptor.pssm_score > exon.acceptor.pssm_score:
                            # poorer score -> do not replace
                            break
                        elif exon.coords()[0] < e.coords()[0]:
                            # replace current exon e with exon
                            print "REPLACING", e, "for", exon
                            abfgp_exons[pos] = exon
                            break
                        elif exon.length > e.length:
                            # replace current exon e with exon
                            print "REPLACING", e, "for", exon
                            abfgp_exons[pos] = exon
                            break
                        else:
                            # smaller first exon -> ignore
                            break
                    if exonrange.intersection(e.dna_range()):
                        break
                else:
                    # if no overlap at all -> append to list
                    abfgp_exons.append( exon.deepcopy() )

    # return the list of ABFGP predicted exons
    return abfgp_exons

# end of function get_raw_abfgp_genestruture


def test_validity_of_final_abfgp_exons(organism,abfgp_genestructure,PCG,
    array_algpresence,array_algsimilarity,verbose=True):
    """ """
    ############################################################################
    if verbose:
        for elem in abfgp_genestructure: print "__GS::", elem
    ############################################################################

    inwpcbgs = PCG2inwpCBGS(PCG,omit_singletons=False)

    exons =   [ abfgp_genestructure[pos] for pos in range(0,len(abfgp_genestructure),2) ]
    introns = [ abfgp_genestructure[pos] for pos in range(1,len(abfgp_genestructure),2) ]

    data = []
    for exclusion_pos in range(0,len(exons)):
        exon            = exons[exclusion_pos]
        orfid           = exon.orf.id
        gene_aa_start   = exon.start / 3 
        gene_aa_stop    = exon.end / 3
        gene_aa_orfend  = exon.orf.protein_endPY
        annotated_orfs  = 0
        num_informants  = 0
        array_range = Set()
        for inwpCBG in inwpcbgs:
            if inwpCBG._get_target_node() == (organism,orfid):
                array_range.update( inwpCBG.maximal_spanning_range(organism=organism) )
                if inwpCBG.count_orfs_labeled_as_annotated_exon() > annotated_orfs:
                    annotated_orfs = inwpCBG.count_orfs_labeled_as_annotated_exon()
                    num_informants = inwpCBG.count_genomic_informants()

        if not array_range:
            # no overlap for this exon AT ALL with ANY PacbPORF.
            # this can happen in case of annotated exons which are added to the abfgp gene model
            score_inwpcbg = 0
        else:
            score_inwpcbg = sum(array_algsimilarity[min(array_range):max(array_range)+1])
        score_exon      = sum(array_algsimilarity[gene_aa_start:gene_aa_stop])
        score_finalexon = sum(array_algsimilarity[gene_aa_start:gene_aa_orfend])
        len_inwpcbg     = len(array_range)
        len_exon        = exon.length / 3
        len_finalexon   = gene_aa_orfend - gene_aa_start
        ########################################################################
        if verbose:
            print "__exon::", annotated_orfs, num_informants, score_inwpcbg,
            print score_exon, score_finalexon,
            print (len_inwpcbg, len_exon, len_finalexon), exon
        ########################################################################
        stat_row = (
                annotated_orfs,
                num_informants,
                score_inwpcbg,
                score_exon,
                score_finalexon,
                (len_inwpcbg, len_exon, len_finalexon) )
        data.append(stat_row)

    # now see if excluding an exon improves the expn explaination score
    best_score  = 0
    best_offset = len(data)-1
    for offset in range(len(data),0,-1):
        current_exon_annot_cnt = data[offset-1][0]
        current_exon_infrm_cnt = data[offset-1][1]
        if offset > 1:
            total_score  = sum( [ data[pos][3] for pos in range(0,offset-1) ] )
            total_score += data[offset-1][4]
            unexpl_score = sum( [ data[pos][2]-data[pos][3] for pos in range(0,offset-1) ] )
            unexpl_score+= ( data[offset-1][2]-data[offset-1][4] )
        else:
            total_score  = data[offset-1][4]
            unexpl_score = data[offset-1][2] - data[offset-1][4]
        # grand total score: explained minus unexplained
        grand_total_score = total_score - unexpl_score
        ########################################################################
        if verbose:
            print offset-1, grand_total_score, total_score,
            print unexpl_score, "annotated::", current_exon_annot_cnt
        ########################################################################
        if grand_total_score > best_score:
            best_score  = grand_total_score
            best_offset = offset-1

        # check if the next iteration is still allowed
        if float(current_exon_annot_cnt) / float(current_exon_infrm_cnt+1) >= 0.33:
            # no, this is an obvious final exon
            # do not allow deletion!
            break

    if best_offset < len(data)-1:
        ########################################################################
        if verbose: print "BEST OFFSET & SCORE::", best_offset, best_score
        ########################################################################
        # truncate this gene model by final exon(s)
        exons =   exons[0:best_offset+1]
        introns = introns[0:best_offset]
        # make stopcodon in stead of donor site in final exon
        wrongexon = exons.pop()
        if exons:
            newexon = FinalExonOnOrf(wrongexon.acceptor,wrongexon.orf.end,wrongexon.orf)
            exons.append(newexon)
        else:
            try:
                newexon = SingleExonOnOrf(wrongexon.acceptor,wrongexon.orf.end,wrongexon.orf)
                exons.append(newexon)
            except:
                # although a Single exon remains, exon does NOT start with an ATG
                # Ignore this here, it must be fixed in another function
                newexon = FinalExonOnOrf(wrongexon.acceptor,wrongexon.orf.end,wrongexon.orf)
                exons.append(newexon)

        # recreate abfgp_genestructure list
        abfgp_genestructure = [ exons.pop(0) ]
        while introns and exons:
            abfgp_genestructure.append(introns.pop(0))
            abfgp_genestructure.append(exons.pop(0))
        # set is_exon_removed pointer to  True
        is_exon_removed = True
    else:
        is_exon_removed = False

    # return the ABFGP genestructure
    return abfgp_genestructure, is_exon_removed

# end of function test_validity_of_final_abfgp_exons


def test_validity_of_first_abfgp_exons(organism,abfgp_genestructure,PCG,
    array_algpresence,array_algsimilarity,verbose=True):
    """ """
    ############################################################################
    if verbose:
        for elem in abfgp_genestructure: print "_first_GS::", elem
    ############################################################################

    inwpcbgs = PCG2inwpCBGS(PCG,omit_singletons=False)

    exons =   [ abfgp_genestructure[pos] for pos in range(0,len(abfgp_genestructure),2) ]
    introns = [ abfgp_genestructure[pos] for pos in range(1,len(abfgp_genestructure),2) ]
    is_exon_removed = False
    
    # do not perform analyses on single-exon genes
    if len(exons) == 1: return abfgp_genestructure, is_exon_removed
    
    data = []
    for exclusion_pos in range(0,len(exons)):
        exon            = exons[exclusion_pos]
        orfid           = exon.orf.id
        gene_aa_start   = exon.start / 3 
        gene_aa_stop    = exon.end / 3
        annotated_orfs  = 0
        num_informants  = 0
        array_range = Set()

        for inwpCBG in inwpcbgs:
            if inwpCBG._get_target_node() == (organism,orfid):
                array_range.update( inwpCBG.maximal_spanning_range(organism=organism) )
                if inwpCBG.count_orfs_labeled_as_annotated_exon() >= annotated_orfs:
                    annotated_orfs = inwpCBG.count_orfs_labeled_as_annotated_exon()
                    num_informants = inwpCBG.count_genomic_informants()
                    
        if not array_range:
            # no overlap for this exon AT ALL with ANY PacbPORF.
            # this can happen in case of annotated exons which are added to the abfgp gene model
            score_inwpcbg = 0
        else:
            score_inwpcbg = sum(array_algsimilarity[min(array_range):max(array_range)+1])

        score_exon      = sum(array_algsimilarity[gene_aa_start:gene_aa_stop])
        score_unexpl    = score_inwpcbg - score_exon
        len_inwpcbg     = len(array_range)

        if exclusion_pos == 0:
            has_alternative = False
            alternativeTSS  = None  
            alternative_score_exon = 0
            alternative_score_unexplained = 0
        else:
            has_alternative = False
            alternativeTSS  = None  
            alternative_score_exon = 0
            alternative_score_unexplained = 0
            if exon.orf.has_methionine():
                exon.orf.scan_orf_for_pssm_tss()
                if exon.orf._tss_sites:
                    score2alternativeTSS = []
                    for tss in exon.orf._tss_sites:
                        _score_exon   = sum(array_algsimilarity[(tss.pos+tss.phase)/3:max(array_range)+1])
                        _score_unexpl = score_inwpcbg - _score_exon
                        _score = _score_exon - _score_unexpl
                        if _score < 0: break
                        score2alternativeTSS.append((_score,_score_exon,_score_unexpl,tss))
                    if score2alternativeTSS:
                        score2alternativeTSS.sort()
                        score2alternativeTSS.reverse()
                        has_alternative = True
                        alternativeTSS  = score2alternativeTSS[0][3]  
                        alternative_score_exon =  score2alternativeTSS[0][1]
                        alternative_score_unexplained =  score2alternativeTSS[0][2]
                        
        stat_row = (
                annotated_orfs,
                num_informants,
                score_inwpcbg,
                score_exon,
                score_unexpl,
                exclusion_pos,
                has_alternative,
                alternativeTSS,
                alternative_score_exon,
                alternative_score_unexplained,
                )
        data.append(stat_row)

    ############################################################################
    if verbose:
        for line in data:
            print "\t".join([ str(elem) for elem in line ])
    ############################################################################
            
    # return the ABFGP genestructure
    return abfgp_genestructure, is_exon_removed

# end of function test_validity_of_first_abfgp_exons



def get_abfgp_predicted_genestructure(organism,PCG,introndata,known_exons,
    array_algpresence,array_algsimilarity,verbose=True):
    """ """
    # get raw list of ABFGP exons
    abfgp_exons = get_raw_abfgp_genestruture(organism,PCG,introndata,known_exons)

    ############################################################################
    if verbose:
        for exon in abfgp_exons: print "raw IE", exon
    ############################################################################

    # Make a list of available exon interfaces
    # The check is required in case of double stopless3n intron occurrence
    abfgp_exon_interfaces = []
    for (prev,next) in [ (_pos-1,_pos) for _pos in range(1,len(abfgp_exons)) ]:
        prevexon = abfgp_exons[prev]
        nextexon = abfgp_exons[next]
        abfgp_exon_interfaces.append( (prevexon.orf.id,nextexon.orf.id) )


    # Make a list of known/current intron coords
    known_intron_coords = []
    for (prev,next) in [ (_pos-1,_pos) for _pos in range(1,len(known_exons)) ]:
        known_intron_coords.append(
            ( known_exons[prev].end, known_exons[next].start ) )

    # find splice site boundary refinements
    for (prev,next) in [ (_pos-1,_pos) for _pos in range(1,len(abfgp_exons)) ]:
        prevexon = abfgp_exons[prev]
        nextexon = abfgp_exons[next]
        thisintronkey = (prevexon.donor.pos,nextexon.acceptor.pos)
        thisintronlabel = (prevexon.orf.id,nextexon.orf.id)
        # simplest case: intron is re-predicted
        if introndata.has_key(thisintronkey): continue
        # if here, look for intron(s) assigned to this label
        exon_is_refined = False
        # check if this exon interface is a stopless3n intron bridgeing
        is_stopless3n_intron_interface = prevexon.orf.id == nextexon.orf.id
        ########################################################################
        if verbose:
            print "start exon refinement:", thisintronkey, thisintronlabel,
            print "stopless3n?:", is_stopless3n_intron_interface
        ########################################################################
        for intron_key in introndata.keys():
            if intron_key in known_intron_coords:
                # A `known` intron from current annotation. This means
                # for this intron no refinement.
                continue
            for thisintron in introndata[intron_key]:
                # Loop over all the individual introns (because is VERY
                # exceptional cases, their interface label is not the same).
                # As soon as a succesfull exon refinement is encountered,
                # this forloop is broken.
                if thisintron._label == thisintronlabel:
                    if is_stopless3n_intron_interface and\
                    abfgp_exon_interfaces.count(thisintron._label) >= 2:
                        # DOUBLE/TRIPLE/... stopless3n intron interfaces
                        ########################################################
                        if verbose: print "DOUBLE stopless3n", thisintron._label
                        ########################################################
                        distDPrev = abs(thisintron.donor.pos-prevexon.donor.pos)
                        distDNext = abs(thisintron.donor.pos-nextexon.donor.pos)
                        distAPrev = abs(thisintron.acceptor.pos-prevexon.acceptor.pos)
                        distANext = abs(thisintron.acceptor.pos-nextexon.acceptor.pos)
                        if distDNext < distDPrev:
                            ####################################################
                            if verbose: print "ignored here",prevexon,nextexon
                            ####################################################
                            # this intron does not belong to this interface
                            break
                        if distAPrev < distANext:
                            ####################################################
                            if verbose: print "ignored here",prevexon,nextexon
                            ####################################################
                            # this intron does not belong to this interface
                            break

                    if thisintron.__class__.__name__ == "SequenceErrorConnectingOrfs":
                        # TODO: make this more generic
                        # temporarily, set pssm_score attribute as float (not None)
                        thisintron.donor.pssm_score    = 1.0
                        thisintron.acceptor.pssm_score = 1.0

                    else:
                        # check if this intron is a strong improvement
                        # based on detailed alignment present / zeros analyses
                        if thisintron.donor.pos != prevexon.donor.pos:
                            # check if this will not result in
                            # an exon with NEGATIVE length
                            if thisintron.donor.pos <= prevexon.acceptor.pos:
                                continue

                            # check in the array_algpresence / array_algsimilarity
                            # if this change makes sence
                            cur_algprs = sum(array_algpresence[prevexon.protein_start():prevexon.protein_end()])
                            cur_algsim = sum(array_algsimilarity[prevexon.protein_start():prevexon.protein_end()])
                            cnt_zeros  = [ p == 0 for p in array_algpresence[prevexon.protein_start():prevexon.protein_end()]].count(True)
    
                            new_algprs = sum(array_algpresence[prevexon.protein_start():thisintron.donor.pos/3])
                            new_algsim = sum(array_algsimilarity[prevexon.protein_start():thisintron.donor.pos/3])
                            new_zeros  = [ p == 0 for p in array_algpresence[prevexon.protein_start():thisintron.donor.pos/3]].count(True)
    
                            # NO IMPROVEMENT, in fact a deterioration!
                            if new_zeros > cnt_zeros: continue
    
                        if thisintron.acceptor.pos != nextexon.acceptor.pos:
                            # check if this will not result in
                            # an exon with NEGATIVE length
                            if thisintron.acceptor.pos >= nextexon.donor.pos:
                                continue

                            # check in the array_algpresence / array_algsimilarity
                            # if this change makes sence
                            cur_algprs = sum(array_algpresence[nextexon.protein_start():nextexon.protein_end()])
                            cur_algsim = sum(array_algsimilarity[nextexon.protein_start():nextexon.protein_end()])
                            cnt_zeros  = [ p == 0 for p in array_algpresence[nextexon.protein_start():nextexon.protein_end()]].count(True)
    
                            new_algprs = sum(array_algpresence[thisintron.acceptor.pos/3:nextexon.protein_end()])
                            new_algsim = sum(array_algsimilarity[thisintron.acceptor.pos/3:nextexon.protein_end()])
                            new_zeros  = [ p == 0 for p in array_algpresence[thisintron.acceptor.pos/3:nextexon.protein_end()]].count(True)
    
                            # NO IMPROVEMENT, in fact a deterioration!
                            if new_zeros > cnt_zeros: continue

                    if thisintron.donor.pos != prevexon.donor.pos:
                        ########################################################
                        if verbose:
                            print "exon donor refinement:", prevexon,
                            print thisintron._label
                            print thisintron
                            for intr in introndata[intron_key]:
                                print "\t",intr.coords(), intr._label
                            print "DONOR COORD:", prevexon.donor.pos, "=>",
                            print thisintron.donor.pos
                            if thisintron.__class__.__name__ !=\
                            "SequenceErrorConnectingOrfs":
                                print "OLD:", (cur_algprs,cur_algsim,cnt_zeros),
                                print "NEW:", (new_algprs,new_algsim,new_zeros)
                        ########################################################
                        # fix prevexon donor
                        if prevexon.IS_FINAL and prevexon.IS_FIRST:
                            prevexon = FirstExonOnOrf(prevexon.acceptor,thisintron.donor,prevexon.orf)
                        elif prevexon.IS_FINAL:
                            prevexon = ExonOnOrf(prevexon.acceptor,thisintron.donor,prevexon.orf)
                        else:
                            prevexon = prevexon.deepcopy() # deepcopy->changes!
                            prevexon.donor      = thisintron.donor
                            prevexon.length     = prevexon.donor.pos -\
                                                  prevexon.acceptor.pos
                            prevexon.end        = prevexon.donor.pos
                            prevexon.pssm_score = prevexon.donor.pssm_score +\
                                                  prevexon.acceptor.pssm_score

                        # place back in list
                        abfgp_exons[prev] = prevexon

                    if thisintron.acceptor.pos != nextexon.acceptor.pos:
                        ########################################################
                        if verbose:
                            print "exon acceptor refinement:", nextexon,
                            print thisintron._label, thisintron,
                            print nextexon.acceptor.pos, "=>",
                            print thisintron.acceptor.pos
                            if thisintron.__class__.__name__ !=\
                            "SequenceErrorConnectingOrfs":
                                print "OLD:", (cur_algprs,cur_algsim,cnt_zeros),
                                print "NEW:", (new_algprs,new_algsim,new_zeros)
                        ########################################################
                        # fix nextexon acceptor
                        if nextexon.IS_FINAL and nextexon.IS_FIRST:
                            nextexon = FinalExonOnOrf(thisintron.acceptor,nextexon.donor,nextexon.orf)
                        elif nextexon.IS_FIRST:
                            nextexon = ExonOnOrf(thisintron.acceptor,nextexon.donor,nextexon.orf)
                        else:
                            nextexon = nextexon.deepcopy() # deepcopy->changes!
                            nextexon.acceptor   = thisintron.acceptor
                            nextexon.length     = nextexon.donor.pos -\
                                                  nextexon.acceptor.pos
                            nextexon.start      = nextexon.acceptor.pos
                            nextexon.pssm_score = nextexon.donor.pssm_score +\
                                                  nextexon.acceptor.pssm_score

                        # place back in list
                        abfgp_exons[next] = nextexon

                    # break out of this loop
                    exon_is_refined = True
                    ############################################################
                    if verbose: print "exon_is_refined = True:", thisintronkey
                    ############################################################
                    break

            # check if refinement was done
            if exon_is_refined: break

    # find OVERPREDICTED stopless 3n introns
    is_replacing = True
    while is_replacing:
        for (prev,next) in [ (pos-1,pos) for pos in range(1,len(abfgp_exons)) ]:
            prevexon = abfgp_exons[prev]
            nextexon = abfgp_exons[next]
            # only asses stopless3n introns -> identical Orfs & phases
            # phase non-indentity is caused by unlikely
            # (but existing) Orf-transitions like a-b-a
            if prevexon.orf.id != nextexon.orf.id: continue
            if prevexon.donor.phase !=  nextexon.acceptor.phase: continue
            # if here, see if this stopless 3n intron is ABFGP confirmed
            stopless3nkey = (prevexon.donor.pos,nextexon.acceptor.pos)
            if not introndata.has_key(stopless3nkey):
                # remove this stopless 3n intron by merging exons
                classPrev = prevexon.__class__.__name__
                classNext = nextexon.__class__.__name__
                if classPrev == 'FirstExonOnOrf' and classNext == 'FinalExonOnOrf':
                    mergedexon = SingleExonOnOrf(prevexon.acceptor,nextexon.donor,prevexon.orf)
                elif classPrev == 'FirstExonOnOrf':
                    mergedexon = FirstExonOnOrf(prevexon.acceptor,nextexon.donor,prevexon.orf)
                elif classNext == 'FinalExonOnOrf':
                    mergedexon = FinalExonOnOrf(prevexon.acceptor,nextexon.donor,prevexon.orf)
                else:
                    mergedexon = ExonOnOrf(prevexon.acceptor,nextexon.donor,prevexon.orf)
                # replace in the list
                abfgp_exons.pop(next)
                abfgp_exons.pop(prev)
                abfgp_exons.insert(prev,mergedexon)
                ########################################################
                if verbose:
                    print "stopless3n intron OVERPREDICTED:",
                    print prevexon, nextexon,
                    print (prevexon.donor.phase, nextexon.acceptor.phase)
                ########################################################
                break    
        else:
            # EOF for loop reached wihtout replacement -> break while loop
            is_replacing = False
                
    for exon in abfgp_exons:
        print "AFTERstopless3n::", exon

    # find UNDERPREDICTED stopless 3n introns
    for intron_key in introndata.keys():
        thisintron = introndata[intron_key][0]
        if thisintron.is_stopless_3n_intron():
            for pos in range(0,len(abfgp_exons)):
                exon = abfgp_exons[pos]
                if thisintron.donor.pos > exon.acceptor.pos and\
                thisintron.acceptor.pos < exon.donor.pos:
                    # inframe intron! split it up...
                    if exon.__class__.__name__ == 'FirstExonOnOrf':
                        exon5p = FirstExonOnOrf(exon.acceptor,thisintron.donor,exon.orf)
                        exon3p = ExonOnOrf(thisintron.acceptor,exon.donor,exon.orf)
                    elif exon.__class__.__name__ == 'FinalExonOnOrf':
                        exon5p = ExonOnOrf(exon.acceptor,thisintron.donor,exon.orf)
                        exon3p = FinalExonOnOrf(thisintron.acceptor,exon.donor,exon.orf)
                    elif exon.__class__.__name__ == 'SingleExonOnOrf':
                        exon5p = FirstExonOnOrf(exon.acceptor,thisintron.donor,exon.orf)
                        exon3p = FinalExonOnOrf(thisintron.acceptor,exon.donor,exon.orf)
                    else:
                        exon5p = ExonOnOrf(exon.acceptor,thisintron.donor,exon.orf)
                        exon3p = ExonOnOrf(thisintron.acceptor,exon.donor,exon.orf)
                    # replace in the list
                    abfgp_exons.pop(pos)
                    abfgp_exons.insert(pos,exon3p)
                    abfgp_exons.insert(pos,exon5p)
                    ########################################################
                    if verbose:
                        print "stopless3n intron UNDERPREDICTED:",
                        print exon5p, exon3p
                        print (exon5p.donor.phase, exon3p.acceptor.phase)
                    ########################################################
                    break    

    # remove FinalExonOnOrfs that should be based on their position ExonOnOrfs
    # BUT have `no` donor sites -> cannot be intermediate exons
    # TODO `no` donor sites must be coded as no donor sites 3' of the PACBP
    # alignment. In the case I am testing now CFU_827898 / MYCGR_083794,
    # by chance this Orf has no donor sites at all ....
    for (prev,next) in [ (pos-1,pos) for pos in range(len(abfgp_exons)-1,0,-1) ]:
        prevexon = abfgp_exons[prev]
        nextexon = abfgp_exons[next]
        if prevexon.IS_FINAL and nextexon.IS_FINAL:
            if prevexon.orf._has_donor_sites_predicted and\
            not prevexon.orf._donor_sites:
                # pop this one out of the GeneStructure
                popped = abfgp_exons.pop(prev)
                print "POPPED non-Final-Final:", popped


    ############################################################################
    # create abfgp_genestructure list by joining the exons with introns
    ############################################################################
    abfgp_genemodel = abfgp_exons2genemodel(abfgp_exons,introndata)


    ############################################################################
    # finetune acceptor boundaries
    ############################################################################
    refined_boundaries = finetune_acceptor_boundaries(abfgp_genemodel,introndata,
            array_algpresence,array_algsimilarity)
    # replace with refined boundaries
    for (old_site,new_site) in refined_boundaries:
        for exon_or_intron in abfgp_genemodel:
            if exon_or_intron == None:
                continue
            if exon_or_intron.__class__.__name__ in ['SequenceErrorConnectingOrfs','SequenceError']: 
                continue
            if exon_or_intron.acceptor.pos == old_site.pos:
                exon_or_intron.acceptor = new_site
                exon_or_intron._init_by_splicesites()
                print "ACCEPTOR finetuned from %s -> %s" % (old_site.pos,new_site.pos) 
        for (donor_pos,acceptor_pos) in introndata.keys():
            if acceptor_pos == old_site.pos:
                new_key = (donor_pos,new_site.pos)
                for thisintron in introndata[(donor_pos,acceptor_pos)]:
                    thisintron.acceptor = new_site
                    thisintron._init_by_splicesites()
                # replace old with new key
                introndata[new_key] = introndata[(donor_pos,acceptor_pos)]
                del( introndata[(donor_pos,acceptor_pos)] )
                break

    ############################################################################
    # finetune donor boundaries
    ############################################################################
    refined_boundaries = finetune_donor_boundaries(abfgp_genemodel,introndata,
            array_algpresence,array_algsimilarity)
    # replace with refined boundaries
    for (old_site,new_site) in refined_boundaries:
        for exon_or_intron in abfgp_genemodel:
            if exon_or_intron == None:
                continue
            if exon_or_intron.__class__.__name__ in ['SequenceErrorConnectingOrfs','SequenceError']:
                continue
            if exon_or_intron.donor.pos == old_site.pos:
                exon_or_intron.donor = new_site
                exon_or_intron._init_by_splicesites()
                print "DONOR finetuned from %s -> %s" % (old_site.pos,new_site.pos) 
        for (donor_pos,acceptor_pos) in introndata.keys():
            if donor_pos == old_site.pos:
                new_key = (new_site.pos,acceptor_pos)
                for thisintron in introndata[(donor_pos,acceptor_pos)]:
                    thisintron.donor = new_site
                    thisintron._init_by_splicesites()
                # replace old with new key
                introndata[new_key] = introndata[(donor_pos,acceptor_pos)]
                del( introndata[(donor_pos,acceptor_pos)] )
                break


    ############################################################################
    # finetune intron boundaries
    ############################################################################
    refined_boundaries = finetune_intron_boundaries(abfgp_genemodel,introndata,
            array_algpresence,array_algsimilarity)
    # replace with refined boundaries
    for (old_site,new_site) in refined_boundaries:
        for exon_or_intron in abfgp_genemodel:
            if exon_or_intron == None:
                continue
            if exon_or_intron.__class__.__name__ in ['SequenceErrorConnectingOrfs','SequenceError']: 
                continue
            if exon_or_intron.donor.pos == old_site.pos:
                exon_or_intron.donor = new_site
                print "SITE (donor on %s) finetuned from %s -> %s (phase:%s->%s)" % (
                        exon_or_intron.__class__.__name__,
                        old_site.pos,new_site.pos,
                        old_site.phase,new_site.phase) 
            if exon_or_intron.acceptor.pos == old_site.pos:
                exon_or_intron.acceptor = new_site
                print "SITE (accep on %s) finetuned from %s -> %s (phase:%s->%s)" % (
                        exon_or_intron.__class__.__name__,
                        old_site.pos,new_site.pos,
                        old_site.phase,new_site.phase)

    # reset attribute values in corrected intron/exon
    # because splice site phase changes, this can be done only
    # all refined sites have been adapted
    new_site_positions = [ new_site.pos for (old_site,new_site) in refined_boundaries ]
    for exon_or_intron in abfgp_genemodel:
        if exon_or_intron.__class__.__name__ in ['SequenceErrorConnectingOrfs','SequenceError']:
            continue
        if exon_or_intron == None:
            continue
        if exon_or_intron.donor.pos in new_site_positions or\
        exon_or_intron.acceptor.pos in new_site_positions:
            # reset attribute values in corrected intron/exon
            print "REINIT OBJECT::",exon_or_intron
            exon_or_intron._init_by_splicesites()

    # replace projected/mapped intron coords with refined boundaries
    for (old_site,new_site) in refined_boundaries:
        is_replaced = True
        # place in while loop because dictionary (introndata) is changed on run-time
        # after each change, while loop is restarted untill all changes have been done
        while is_replaced:
            for (donor_pos,acceptor_pos) in introndata.keys():
                if donor_pos == old_site.pos and acceptor_pos in [ old.pos for old,new in refined_boundaries ]:
                    # adapt donor positions of introns with an OLD to-be-refined acceptor position
                    new_key = (new_site.pos,acceptor_pos)
                    for thisintron in introndata[(donor_pos,acceptor_pos)]:
                        thisintron.donor = new_site
                        # do NOT reinit here -> possible IncompatibleSpliceSitePhases
                    # replace old with new key
                    introndata[new_key] = introndata[(donor_pos,acceptor_pos)]
                    del( introndata[(donor_pos,acceptor_pos)] )
                    break

                if acceptor_pos == old_site.pos and donor_pos in [ new.pos for old,new in refined_boundaries ]:
                    # adapt acceptor positions of introns with an just NEWLY refined donor position
                    new_key = (donor_pos,new_site.pos)
                    for thisintron in introndata[(donor_pos,acceptor_pos)]:
                        thisintron.acceptor = new_site
                        # do NOT reinit here -> possible IncompatibleSpliceSitePhases
                    # replace old with new key
                    introndata[new_key] = introndata[(donor_pos,acceptor_pos)]
                    del( introndata[(donor_pos,acceptor_pos)] )
                    break

            else:
                # EOF for loop reached -> no replacements done
                is_replaced = False

    # re-init introns here; earlier can cause IncompatibleSpliceSitePhases exception
    for (donor_pos,acceptor_pos) in introndata.keys():
        if donor_pos in new_site_positions or\
        acceptor_pos in new_site_positions:
            # reset attribute values in corrected intron/exon
            print "REINIT introns::", (donor_pos,acceptor_pos), introndata[(donor_pos,acceptor_pos)][0]
            for thisintron in introndata[(donor_pos,acceptor_pos)]:
                thisintron._init_by_splicesites()

                

    ############################################################################
    # check if very poor 3' exons are predicted, which cause exons with strong
    # coding signal (array_similarity) to be poorly explained.
    ############################################################################
    abfgp_genemodel, is_final_exon_removed = test_validity_of_final_abfgp_exons(organism,abfgp_genemodel,PCG,
        array_algpresence,array_algsimilarity)

    if is_final_exon_removed:
        ############################################################################
        # replace list with ABFGP exons -> call by reference didn't update these!
        ############################################################################
        abfgp_exons   = [ abfgp_genemodel[pos] for pos in range(0,len(abfgp_genemodel),2) ]
        abfgp_introns = [ abfgp_genemodel[pos] for pos in range(1,len(abfgp_genemodel),2) ]


        
    abfgp_genemodel, is_first_exon_removed = test_validity_of_first_abfgp_exons(organism,abfgp_genemodel,PCG,
        array_algpresence,array_algsimilarity)
        
        
        
    ############################################################################
    # Check if there is any exon predicted; in case very poor similarity is reported,
    # not a single convincing exon is reported
    # Report empty exon/intron lists
    ############################################################################
    if not abfgp_genemodel: return [], []


    ############################################################################
    # GFF defaults. TODO -> move to settings.abfgp / settings.abfgp
    ############################################################################
    ABFGP_EXON_GFF_DEFAULTS = {
            'fsource':  'abfgp',
            'fmethod':  'AbfgpExon',
            'gclass':   'AbfgpGeneStructure' }
    ABFGP_INTRON_GFF_DEFAULTS = {
            'fsource':  'abfgp',
            'fmethod':  'AbfgpIntron',
            'gclass':   'AbfgpGeneStructure' }
    ABFGP_SEQERROR_GFF_DEFAULTS = {
            'fsource':  'abfgp',
            'fmethod':  'AbfgpSequenceError',
            'gclass':   'AbfgpGeneStructure' }

    for exon in abfgp_exons:
        exon._gff.update(ABFGP_EXON_GFF_DEFAULTS)
        exon.pssm_score = 1

    ############################################################################
    # create a Set of DNA positions that are already covered by introns/exons
    ############################################################################
    assigned_coords = Set()
    for elem in abfgp_genemodel:
        if elem: assigned_coords.update(elem.dna_range())
        print "GS ELEM:", elem,
        try:    print elem.pssm_score
        except: print None

        
    ############################################################################
    # add (small) Exons which could not be ABFGP confirmed
    # TODO TODO TODO: this code must be extended much more!!!
    # are potential FinalExonOnOrfs and FirstExonOnOrfs:
    # - identical in length?
    # - have ATGs?
    # - have non-sequence-similarity properties?
    # - have strong splice site consensus?
    ############################################################################
    add_known_exons = _select_known_exons_to_include_in_abfgp_genemodel(
                        known_exons,abfgp_genemodel)
    ############################################################################
    if verbose:
        for ake in add_known_exons:
            print "ABOUT TO ADD GS EXON:", ake
    ############################################################################
    is_any_exon_added = False
    for exon in add_known_exons:
        # if here -> add exon to abfgp_exons (as a not-confirmable one)
        is_any_exon_added = True
        # update features of this NON-ABFGP-CONFIRMED exon
        exon._gff.update(ABFGP_EXON_GFF_DEFAULTS)
        exon.pssm_score = 0

        # place in the correct position towards existing exons
        if not abfgp_exons:
            abfgp_exons = [ exon ]
        elif exon.end < abfgp_exons[0].start:
            abfgp_exons.insert(0,exon)
            ####################################################################
            if verbose: print "GeneModelExon unshifted:", exon
            ####################################################################
        elif exon.start > abfgp_exons[-1].end:
            abfgp_exons.append(exon)
            ####################################################################
            if verbose: print "GeneModelExon appended:", exon
            ####################################################################
        else:
            for pos in range(0,len(abfgp_exons)):
                if exon.end < abfgp_exons[pos].start:
                    ############################################################
                    if verbose: print "GeneModelExon inserted:", exon
                    ############################################################
                    abfgp_exons.insert(pos,exon)
                    # check phase compatibily around hard-insert
                    if not abfgp_exons[pos-1].donor.phase == exon.acceptor.phase:
                        # correct to previous GeneModel exon's donor
                        print abfgp_exons[pos-1].donor.phase == exon.acceptor.phase,
                        print "corrected:", 
                        try:
                            abfgp_exons[pos-1].donor = known_exons[mainpos-1].donor
                        except:
                            pass
                        print abfgp_exons[pos-1].donor.phase == exon.acceptor.phase
                    if not abfgp_exons[pos+1].acceptor.phase == exon.donor.phase:
                        # correct to next GeneModel exon's acceptor
                        print abfgp_exons[pos+1].acceptor.phase == exon.donor.phase,
                        print "corrected:", 
                        try:
                            abfgp_exons[pos+1].acceptor = known_exons[mainpos+1].acceptor
                        except:
                            pass
                        print abfgp_exons[pos+1].acceptor.phase == exon.donor.phase

                    # and break out
                    break

    if is_any_exon_added:
        # recreate abfgp_genestructure list by joining the exons with introns
        abfgp_genemodel = abfgp_exons2genemodel(abfgp_exons,introndata)

        ############################################################################
        # SECOND ROUND for validity_of_final_abfgp_exons
        # this can remove just added annotated exons again from the abfgp genemodel!
        # check if very poor 3' exons are predicted, which cause exons with strong
        # coding signal (array_similarity) to be poorly explained.
        ############################################################################
        abfgp_genemodel, is_final_exon_removed = test_validity_of_final_abfgp_exons(organism,abfgp_genemodel,PCG,
            array_algpresence,array_algsimilarity)

        if is_final_exon_removed:
            print "####################################################"
            print "####################################################"
            print "removed in 2th instance!!!"
            print "####################################################"
            print "####################################################"

            ############################################################################
            # replace list with ABFGP exons -> call by reference didn't update these!
            ############################################################################
            abfgp_exons   = [ abfgp_genemodel[pos] for pos in range(0,len(abfgp_genemodel),2) ]
            abfgp_introns = [ abfgp_genemodel[pos] for pos in range(1,len(abfgp_genemodel),2) ]
        
    # now have a closer look at the abfgp_genemodel. Are First & Final exons there?
    if not abfgp_genemodel[-1].IS_FINAL:
        # make it final!
        # TODO: implement a check if it make sence to skip this one..
        # TODO: write this code at the position where annotated exons are added
        curexon = abfgp_genemodel[-1]
        if len(abfgp_genemodel) == 1:
            if curexon.IS_FIRST:
                newexon = SingleExonOnOrf(curexon.acceptor,curexon.orf.end,curexon.orf)
            else:
                # hmmm... make FirstExon of this one
                atg = _propose_atg_on_orf(curexon.orf,None)
                if atg:
                    newexon = SingleExonOnOrf(atg,curexon.orf.end,curexon.orf)
                else:
                    # cannot be assigned as the first one -> no ATG!
                    # TODO: solve this. Allow fractions!? Or mine in
                    # annotated exon set again for a proper first exon?
                    print "NO TSS/ATG/M on Orf -> no SingleExon"
                    newexon = curexon
        else:
            newexon = FinalExonOnOrf(curexon.acceptor,curexon.orf.end,curexon.orf)
        # update _gff data / pssm_score and place back in the abfgp_genemodel
        newexon._gff.update(curexon._gff)
        newexon.pssm_score = curexon.pssm_score
        abfgp_genemodel[-1] = newexon

    # CHECK for FirstExon
    # TODO: make this code much more extended
    if abfgp_genemodel[0].proteinsequence():
        if abfgp_genemodel[0].proteinsequence()[0] != "M":
            # label first exon - not starting with an Methionine, as uncertain
            print "FIRST EXON not a Methionine -> pssm adjusted::", abfgp_genemodel[0]
            abfgp_genemodel[0].pssm_score = 0
    else:
        # first coding exon of 2 nucleotides!
        pass
           
            
    ############################################################################
    # recreate lists of ABGFP exons & introns
    ############################################################################
    abfgp_exons   = [ abfgp_genemodel[pos] for pos in range(0,len(abfgp_genemodel),2) ]
    abfgp_introns = [ abfgp_genemodel[pos] for pos in range(1,len(abfgp_genemodel),2) ]

    # make 'True' introns on the postions of None elements
    for pos in range(0,len(abfgp_introns)):
        intron = abfgp_introns[pos]
        if intron != None:
            if intron.__class__.__name__ == 'SequenceErrorConnectingOrfs':
                intron._gff.update(ABFGP_SEQERROR_GFF_DEFAULTS)
            else:
                intron._gff.update(ABFGP_INTRON_GFF_DEFAULTS)
            intron._gff['gname'] = 'intron_%s_%s' % intron.coords()
            intron.pssm_score = 1
            # Check how reliable this intron is.
            # In introndata some very weak introns are present too!
            if intron.coords() in introndata.keys():
                # TODO: we need the _calc_aligned_score function here
                # TODO: we need the assessed_interfaces variable here
                pass
            else:
                # not present !?!? hmmmm
                print "WARNING:: intron", intron.coords(), "not in introndata!?" 
                pass

        else:
            # if here, make intron
            prev_exon = abfgp_exons[pos]
            next_exon = abfgp_exons[pos+1]
            try:
                intron = IntronConnectingOrfs(
                    prev_exon.donor,next_exon.acceptor,None,
                    prev_exon.orf,next_exon.orf
                    )

                # set data attributes to intron
                intron._gff.update(ABFGP_INTRON_GFF_DEFAULTS)
                intron._gff['gname'] = 'intron_%s_%s' % intron.coords()
                intron.pssm_score = 0
                # set back to intron list
                abfgp_introns[pos] = intron

            except IncompatibleSpliceSitePhases:
                # dunno why
                print "gene.exceptions.IncompatibleSpliceSitePhases"
                print prev_exon
                print next_exon
                print prev_exon.donor
                print next_exon.acceptor
                for exon in abfgp_exons:
                    print "E", exon
            except:
                # dunno why
                print "Other Exception"
                print prev_exon
                print next_exon
                print prev_exon.donor
                print next_exon.acceptor
                for exon in abfgp_exons:
                    print "E", exon


    # update exon with GFF defaults for proper recognition
    for exon in abfgp_exons: exon._gff.update(ABFGP_EXON_GFF_DEFAULTS)
                    
    # return ABFGP predicted gene structure
    return abfgp_exons, abfgp_introns

# end of function get_abfgp_predicted_genestructure




def finetune_acceptor_boundaries(abfgp_genemodel,introndata,
    array_algpresence,array_algsimilarity,verbose=True):
    """
    Finetune *ACCEPTOR* boundaries of predicted introns
    """

    # Global Variable Imports
    from settings.genestructure import MIN_INTRON_NT_LENGTH
    FINETUNE_ACCEPTOR_NT_OFFSET = 18

    # list with adjusted boundaries
    refined_boundaries = []

    # recreate lists of ABGFP exons & introns
    abfgp_exons   = [ abfgp_genemodel[pos] for pos in range(0,len(abfgp_genemodel),2) ]
    abfgp_introns = [ abfgp_genemodel[pos] for pos in range(1,len(abfgp_genemodel),2) ]

    for intron_pos in range(0,len(abfgp_introns)):
        intron = abfgp_introns[intron_pos]
        if not intron: continue
        if intron.__class__.__name__ == 'SequenceErrorConnectingOrfs': continue

        # assign branchpoint in current intron
        intron.assign_bp_and_ppts()

        has_been_printed = False

        # list of alternatives & associated scores
        alternatives = []
        finetune_range = range(intron.acceptor.pos-FINETUNE_ACCEPTOR_NT_OFFSET,
                intron.acceptor.pos+FINETUNE_ACCEPTOR_NT_OFFSET+1,3)

        for acceptor in intron.orfAcceptor._acceptor_sites:
            if acceptor.pos != intron.acceptor.pos and\
            acceptor.pos in finetune_range:
                # get the next exon (3'of this intron)
                next_exon = abfgp_exons[intron_pos+1]
                if not has_been_printed:
                    has_been_printed = True
                    ############################################################
                    if verbose: print "FINETUNING ACCEPTOR::", intron
                    ############################################################

                # get data on this alternative acceptor position
                test_intron = IntronConnectingOrfs(intron.donor,acceptor,None,intron.orfDonor,intron.orfAcceptor)
                test_intron.assign_bp_and_ppts()

                # test if refinement will result in a long enough intron
                if test_intron.length < MIN_INTRON_NT_LENGTH: continue

                scorelist = []
                # score 1: is acceptor.pssm_score `higher`?
                scorelist.append( _finetune_splicesite_comparison(intron.acceptor,test_intron.acceptor) )
                # score 2: branchpoint comparison?
                scorelist.append( _branchpoint_comparison(intron,test_intron) )
                # score 3: ppt comparison?
                scorelist.append( _polypyrimidinetract_comparison(intron,test_intron) )
                # score 4: is algsimilarity ratio increased (==better)?
                scorelist.append( _algsimilarity_comparison(intron,test_intron,None,next_exon,array_algsimilarity) )

                # evaluate scorelist; improved intron boundary or not?
                # use acceptor, branchpoint & ppt, do *NOT* use algsim score
                if scorelist[0:3].count(False) == 0 and scorelist[0:3].count(True) >= 1:
                    alternatives.append( ( acceptor, scorelist ) )
                    is_accepted = True
                else:
                    is_accepted = False

                ################################################################
                if verbose:
                    print "alternative:", acceptor,
                    print intron.acceptor.pos - acceptor.pos, scorelist,
                    print is_accepted, "BP:",
                    print intron.get_branchpoint_nt_distance(),
                    print "alt:",
                    print test_intron.get_branchpoint_nt_distance()
                ################################################################

        # now evaluate the alternatived and take the best one
        if not alternatives:
            continue
        elif len(alternatives) == 1:
            refined_boundaries.append( ( intron.acceptor, alternatives[0][0] ) )
        else:
            # multiple! again, take the *best* one
            pass

    # return list of refined_boundaries
    return refined_boundaries

# end of function finetune_acceptor_boundaries


def finetune_donor_boundaries(abfgp_genemodel,introndata,
    array_algpresence,array_algsimilarity,verbose=True):
    """
    Finetune *DONOR* boundaries of predicted introns
    """

    # Global Variable Imports
    from settings.genestructure import MIN_INTRON_NT_LENGTH
    FINETUNE_DONOR_NT_OFFSET = 18

    # list with adjusted boundaries
    refined_boundaries = []

    # recreate lists of ABGFP exons & introns
    abfgp_exons   = [ abfgp_genemodel[pos] for pos in range(0,len(abfgp_genemodel),2) ]
    abfgp_introns = [ abfgp_genemodel[pos] for pos in range(1,len(abfgp_genemodel),2) ]

    for intron_pos in range(0,len(abfgp_introns)):
        intron = abfgp_introns[intron_pos]
        if not intron: continue
        if intron.__class__.__name__ == 'SequenceErrorConnectingOrfs': continue
        has_been_printed = False

        # list of alternatives & associated scores
        alternatives = []
        finetune_range = range(intron.donor.pos-FINETUNE_DONOR_NT_OFFSET,
                intron.donor.pos+FINETUNE_DONOR_NT_OFFSET+1,3)

        for donor in intron.orfDonor._donor_sites:
            if donor.pos != intron.donor.pos and\
            donor.pos in finetune_range:
                # get the prev exon (5'of this intron)
                prev_exon = abfgp_exons[intron_pos]
                if not has_been_printed:
                    has_been_printed = True
                    ############################################################
                    if verbose: print "FINETUNING DONOR::", intron
                    ############################################################

                # get data on this alternative donor position
                test_intron = IntronConnectingOrfs(donor,intron.acceptor,None,intron.orfDonor,intron.orfAcceptor)

                # test if refinement will result in a long enough intron
                if test_intron.length < MIN_INTRON_NT_LENGTH: continue

                scorelist = []
                # score 1: is acceptor.pssm_score `higher`?
                scorelist.append( _finetune_splicesite_comparison(intron.donor,test_intron.donor) )
                # score 2: is algsimilarity ratio increased (==better)?
                scorelist.append( _algsimilarity_comparison(intron,test_intron,prev_exon,None,array_algsimilarity) )

                # evaluate scorelist; improved intron boundary or not?
                # use donor and algsim score
                if scorelist.count(True) == 2:
                    alternatives.append( ( donor, scorelist ) )
                    is_accepted = True
                else:
                    is_accepted = False

                ################################################################
                if verbose:
                    print "alternative:", donor,
                    print intron.donor.pos - donor.pos, scorelist, is_accepted
                ################################################################

        # now evaluate the alternatived and take the best one
        if not alternatives:
            continue
        elif len(alternatives) == 1:
            refined_boundaries.append( ( intron.donor, alternatives[0][0] ) )
        else:
            # multiple! again, take the *best* one
            pass

    # return list of refined_boundaries
    return refined_boundaries

# end of function finetune_donor_boundaries



def finetune_intron_boundaries(abfgp_genemodel,introndata,
    array_algpresence,array_algsimilarity,verbose=True):
    """
    Finetune Donor and Acceptor (by changing phase) boundaries of predicted introns
    """

    # Global Variable Imports
    FINETUNE_ACCEPTOR_NT_OFFSET = 12
    FINETUNE_DONOR_NT_OFFSET = 12
    FINETUNE_ACCEPTOR_NT_OFFSET = 18
    FINETUNE_DONOR_NT_OFFSET    = 18
    from settings.genestructure import MIN_INTRON_NT_LENGTH

    # list with adjusted boundaries
    refined_boundaries = []

    # recreate lists of ABGFP exons & introns
    abfgp_exons   = [ abfgp_genemodel[pos] for pos in range(0,len(abfgp_genemodel),2) ]
    abfgp_introns = [ abfgp_genemodel[pos] for pos in range(1,len(abfgp_genemodel),2) ]

    for intron_pos in range(0,len(abfgp_introns)):
        intron = abfgp_introns[intron_pos]
        if not intron: continue
        if intron.__class__.__name__ == 'SequenceErrorConnectingOrfs': continue
        has_been_printed = False
        finetune_acceptor_range = range(intron.acceptor.pos-FINETUNE_ACCEPTOR_NT_OFFSET,
                intron.acceptor.pos+FINETUNE_ACCEPTOR_NT_OFFSET+1)
        finetune_donor_range = range(intron.donor.pos-FINETUNE_DONOR_NT_OFFSET,
                intron.donor.pos+FINETUNE_DONOR_NT_OFFSET+1)

        # assign branchpoint in current intron
        intron.assign_bp_and_ppts()

        # start searching acceptor based
        alternatives = []
        for acceptor in intron.orfAcceptor._acceptor_sites:
            if acceptor.pos != intron.acceptor.pos and\
            acceptor.phase != intron.acceptor.phase and\
            acceptor.pos in finetune_acceptor_range:
                # now see if we can find a donor for this phase too
                for donor in intron.orfDonor._donor_sites:
                    if donor.pos != intron.donor.pos and\
                    donor.phase != intron.donor.phase and\
                    donor.phase == acceptor.phase and\
                    donor.pos in finetune_donor_range:
                        # get the next exon (3'of this intron)
                        next_exon = abfgp_exons[intron_pos+1]
                        prev_exon = abfgp_exons[intron_pos]

                        if not has_been_printed:
                            has_been_printed = True
                            ####################################################
                            if verbose: print "FINETUNING INTRON::", intron
                            ####################################################

                        # get data on this alternative acceptor/donor combination
                        test_intron = IntronConnectingOrfs(donor,acceptor,None,intron.orfDonor,intron.orfAcceptor)
                        test_intron.assign_bp_and_ppts()

                        # test if refinement will result in a long enough intron
                        if test_intron.length < MIN_INTRON_NT_LENGTH: continue

                        scorelist = []
                        # score 1: is donor.pssm_score `higher`?
                        scorelist.append( _finetune_splicesite_comparison(intron.donor,donor) )
                        # score 2: is acceptor.pssm_score `higher`?
                        scorelist.append( _finetune_splicesite_comparison(intron.acceptor,acceptor) )
                        # score 3: branchpoint comparison?
                        scorelist.append( _branchpoint_comparison(intron,test_intron) )
                        # score 4: ppt comparison?
                        scorelist.append( _polypyrimidinetract_comparison(intron,test_intron) )
                        # score 5: is algsimilarity ratio increased (==better)?
                        scorelist.append( _algsimilarity_comparison(intron,test_intron,prev_exon,next_exon,array_algsimilarity) )

                        # evaluate scorelist; improved intron boundary or not?
                        # use donor, acceptor, branchpoint & ppt, do *NOT* use algsim score
                        if scorelist[0:4].count(False) == 0 and scorelist[0:4].count(True) >= 1:
                            alternatives.append( ( donor, acceptor, scorelist ) )
                            is_accepted = True
                        else:
                            is_accepted = False

                        ########################################################
                        if verbose:
                            print "alternatives:", donor, acceptor,
                            print intron.donor.pos - donor.pos,
                            print intron.acceptor.pos - acceptor.pos,
                            print scorelist, is_accepted,
                            print "BPcur:",intron.get_branchpoint_nt_distance(),
                            print "alt:",
                            print test_intron.get_branchpoint_nt_distance()
                        ########################################################

        # now evaluate the alternatived and take the best one
        if not alternatives:
            continue
        elif len(alternatives) == 1:
            refined_boundaries.append( ( intron.donor,    alternatives[0][0] ) )
            refined_boundaries.append( ( intron.acceptor, alternatives[0][1] ) )
        else:
            # multiple! again, take the *best* one
            pass

    # return list of refined_boundaries
    return refined_boundaries

# end of function finetune_intron_boundaries



def correct_intron_boundaries_by_exonsignal(abfgp_genemodel,organism,PCG,introndata,
    array_algpresence,array_algsimilarity,verbose=True):
    """
    Detect & improve discrepancies in exon vs. intron signal
    """

    # Global Variable Imports
    CORRECT_BY_EXONSIGNAL_ACCEPTOR_NT_OFFSET = 200
    CORRECT_BY_EXONSIGNAL_DONOR_NT_OFFSET    = 200
    FINETUNE_ACCEPTOR_NT_OFFSET = 12
    FINETUNE_DONOR_NT_OFFSET = 12
    from settings.genestructure import MIN_INTRON_NT_LENGTH


    # list with adjusted boundaries
    refined_boundaries = []

    # recreate lists of ABGFP exons & introns
    abfgp_exons   = [ abfgp_genemodel[pos] for pos in range(0,len(abfgp_genemodel),2) ]
    abfgp_introns = [ abfgp_genemodel[pos] for pos in range(1,len(abfgp_genemodel),2) ]

    for intron_pos in range(0,len(abfgp_introns)):
        intron = abfgp_introns[intron_pos]
        if not intron: continue
        if intron.__class__.__name__ == 'SequenceErrorConnectingOrfs': continue

        next_exon = abfgp_exons[intron_pos+1]
        prev_exon = abfgp_exons[intron_pos]

        # make elegiable correction range
        correct_acceptor_range = range(
                intron.acceptor.pos-CORRECT_BY_EXONSIGNAL_ACCEPTOR_NT_OFFSET,
                min([ next_exon.donor.pos-6, intron.acceptor.pos+CORRECT_BY_EXONSIGNAL_ACCEPTOR_NT_OFFSET+1 ]) )
        correct_donor_range = range(
                max([ prev_exon.acceptor.pos+6, intron.donor.pos-CORRECT_BY_EXONSIGNAL_DONOR_NT_OFFSET]),
                intron.donor.pos+CORRECT_BY_EXONSIGNAL_DONOR_NT_OFFSET+1)
        finetune_acceptor_range = range(intron.acceptor.pos-FINETUNE_ACCEPTOR_NT_OFFSET,
                intron.acceptor.pos+FINETUNE_ACCEPTOR_NT_OFFSET+1)
        finetune_donor_range = range(intron.donor.pos-FINETUNE_DONOR_NT_OFFSET,
                intron.donor.pos+FINETUNE_DONOR_NT_OFFSET+1)


        staI,endI = intron.donor.pos/3, intron.acceptor.pos/3
        staE,endE = intron.acceptor.pos/3, next_exon.donor.pos/3
        algsim_intron = sum(array_algsimilarity[staI:endI])
        algsim_next_exon = sum(array_algsimilarity[staE:endE])
        algsim_score_abs = float(algsim_intron)/algsim_next_exon
        algsim_score_rel = float(algsim_intron)/algsim_next_exon * (float(endE-staE)/float(endI-staI))

        print "CORRECT", intron
        print intron.acceptor, intron.acceptor.pos, algsim_score_abs, algsim_score_rel
        print [ a.pos for a in intron.orfAcceptor._acceptor_sites ]


        possible_acceptor_boundaries = []
        possible_donor_boundaries = []

        for acceptor in intron.orfAcceptor._acceptor_sites:
            if acceptor.pos - intron.donor.pos < MIN_INTRON_NT_LENGTH: continue
            if acceptor.pos not in correct_acceptor_range: continue
            if acceptor.pos == intron.acceptor.pos: continue
            if _finetune_splicesite_comparison(intron.acceptor,acceptor) == False: continue

            # if here, generate stats
            staI,endI = intron.donor.pos/3, acceptor.pos/3
            algsim_intron_pos    = sum(array_algsimilarity[staI:endI])
            staE,endE = acceptor.pos/3, next_exon.donor.pos/3
            algsim_next_exon_pos = sum(array_algsimilarity[staE:endE]) 
            # check if no similarity left of exon; avoid ZeroDivisionError
            if algsim_next_exon_pos == 0: continue

            # calculate score ratios
            alternative_algsim_score_abs = float(algsim_intron_pos)/algsim_next_exon_pos
            alternative_algsim_score_rel = alternative_algsim_score_abs * (float(endE-staE)/float(endI-staI))

            # test if exon signal is an improvement
            #if alternative_algsim_score_abs > algsim_score_abs: continue
            #if alternative_algsim_score_rel > algsim_score_rel: continue


            if acceptor.phase == intron.donor.phase:
                # get data on this alternative acceptor/donor combination
                test_intron = IntronConnectingOrfs(intron.donor,acceptor,None,intron.orfDonor,intron.orfAcceptor)
                test_intron.assign_bp_and_ppts()
            else:
                test_intron = None

            print acceptor, acceptor.pos, 
            print alternative_algsim_score_abs,
            print alternative_algsim_score_rel,
            if test_intron:
                print _branchpoint_comparison(intron,test_intron),
                print _polypyrimidinetract_comparison(intron,test_intron)
            else:
                print "OTHER PHASE"


        staI,endI = intron.donor.pos/3, intron.acceptor.pos/3
        staE,endE = prev_exon.acceptor.pos/3, intron.donor.pos/3, 
        algsim_intron = sum(array_algsimilarity[staI:endI])
        algsim_prev_exon = sum(array_algsimilarity[staE:endE])
        algsim_score_abs = float(algsim_intron)/algsim_prev_exon
        algsim_score_rel = float(algsim_intron)/algsim_prev_exon * (float(endE-staE)/float(endI-staI))

        print intron.donor, intron.donor.pos, algsim_score_abs, algsim_score_rel
        print [ d.pos for d in intron.orfDonor._donor_sites ]

        for donor in intron.orfDonor._donor_sites:
            if intron.acceptor.pos - donor.pos < MIN_INTRON_NT_LENGTH: continue
            if donor.pos == intron.donor.pos: continue
            if donor.pos not in correct_donor_range: continue
            if _finetune_splicesite_comparison(intron.donor,donor) == False: continue

            # if here, generate stats
            staI,endI = donor.pos/3, intron.acceptor.pos/3
            algsim_intron_pos    = sum(array_algsimilarity[staI:endI])
            staE,endE = prev_exon.acceptor.pos/3, donor.pos/3
            algsim_prev_exon_pos = sum(array_algsimilarity[staE:endE]) 
            # check if no similarity left of exon; avoid ZeroDivisionError
            if algsim_prev_exon_pos == 0: continue

            # calculate score ratios
            alternative_algsim_score_abs = float(algsim_intron_pos)/algsim_prev_exon_pos
            alternative_algsim_score_rel = alternative_algsim_score_abs * (float(endE-staE)/float(endI-staI))

            # test if exon signal is an improvement
            #if alternative_algsim_score_abs > algsim_score_abs: continue
            #if alternative_algsim_score_rel > algsim_score_rel: continue

            if donor.phase == intron.acceptor.phase:
                # get data on this alternative acceptor/donor combination
                test_intron = IntronConnectingOrfs(donor,intron.acceptor,None,intron.orfDonor,intron.orfAcceptor)
                test_intron.assign_bp_and_ppts()
            else:
                test_intron = None

            print donor, donor.pos, 
            print alternative_algsim_score_abs,
            print alternative_algsim_score_rel,
            if test_intron:
                print _branchpoint_comparison(intron,test_intron),
                print _polypyrimidinetract_comparison(intron,test_intron)
            else:
                print "OTHER PHASE"


        #for position in range(intron.acceptor.pos,intron.orfAcceptor.start,-10):
        #    staI,endI = intron.donor.pos/3, position/3
        #    if endI <= staI: break
        #    algsim_intron_pos    = sum(array_algsimilarity[staI:endI])
        #    staE,endE = position/3, next_exon.donor.pos/3
        #    algsim_next_exon_pos = sum(array_algsimilarity[staE:endE]) 
        #    print position, scoreA, scoreB,
        #    print float(algsim_intron_pos)/algsim_next_exon_pos,
        #    print float(algsim_intron_pos)/algsim_next_exon_pos * (float(endE-staE)/float(endI-staI))

# end of function correct_intron_boundaries_by_exonsignal



def abfgp_exons2genemodel(abfgp_exons,introndata):
    """ """
    # create abfgp_genestructure list by joining the exons with introns
    abfgp_genemodel = [  ]

    # if there is just a single exon -> singe exon genemodel!
    if len(abfgp_exons) == 1:
        abfgp_genemodel = abfgp_exons

    for (prev,next) in [ (pos-1,pos) for pos in range(1,len(abfgp_exons)) ]:
        prevexon = abfgp_exons[prev]
        nextexon = abfgp_exons[next]
        # append the first exon in the first iteration through abfgp_exons
        if (prev,next) == (0,1): abfgp_genemodel.append(prevexon)
        intronkey = (prevexon.donor.pos,nextexon.acceptor.pos)
        if introndata.has_key(intronkey):
            example = introndata[intronkey][0] # get any/first intron as example
            #if example.__class__.__name__ == "SequenceErrorConnectingOrfs":
            #    intron = example.deepcopy()
            #else:
            #    intron = IntronConnectingOrfs(
            #        example.donor,example.acceptor,example.shared_nts,
            #        example.orfDonor,example.orfAcceptor
            #        )
            intron = example.deepcopy()
            # remove 'Note' attribute
            if intron._gff['column9data'].has_key("Note"):
                del( intron._gff['column9data']["Note"] )
            # TODO add info on which informants contributed to this intron
        else:
            print "NOINTRON::", intronkey, prevexon, nextexon
            intron = None
        # append intron and exon the genemodel
        abfgp_genemodel.append(intron)
        abfgp_genemodel.append(nextexon)

    # return gene model list
    return abfgp_genemodel

# end of function abfgp_exons2genemodel



def _propose_atg_on_orf(orfObj,inwpCBG):
    """ Propose the most likely TSS/Methionine on this Orf """ 
    # donor only -> create FirstExonOnOrf (and hope for the best)
    if orfObj.has_methionine():
        # RESCAN TSS pssm sites in any case. List might have been edited
        # in exceptional cases somewhere earlier in ABFGP  flow
        orfObj.scan_orf_for_pssm_tss(forced=True)
        if orfObj._tss_sites and inwpCBG:
            maxsr = inwpCBG.maximal_spanning_range(organism=inwpCBG._get_target_organism())
            minsr = inwpCBG.minimal_spanning_range(organism=inwpCBG._get_target_organism())
            firstexon_array_algpresence   = PCG2codingarray(inwpCBG,inwpCBG._get_target_organism(),max(maxsr)+1)
            firstexon_array_algsimilarity = PCG2similarityarray(inwpCBG,inwpCBG._get_target_organism(),max(maxsr)+1)
            score2tssobj = []
            for tss in orfObj._tss_sites:
                tsspos = (tss.pos - orfObj.frame)/3
                exon_estim_end = max(minsr)+1
                len5p = tsspos - min(maxsr)
                len3p = exon_estim_end - tsspos
                score5P = sum(firstexon_array_algpresence[min(maxsr):tsspos])
                score3P = sum(firstexon_array_algpresence[tsspos:exon_estim_end])
                score5S = sum(firstexon_array_algsimilarity[min(maxsr):tsspos])
                score3S = sum(firstexon_array_algsimilarity[tsspos:exon_estim_end])
                if len3p > 0:
                    if len5p > 0:
                        scoreP = float(score3P)/float(len3p) - float(score5P)/float(len5p)
                        scoreS = float(score3S)/float(len3p) - float(score5S)/float(len5p)
                    elif len5p < 0:
                        scoreP = float(score3P)/float(len3p) - 1.0 - float(-len5p)/50.0
                        scoreS = float(score3S)/float(len3p) - 1.0 - float(-len5p)/50.0
                    else:
                        scoreP = float(score3P)/float(len3p)
                        scoreS = float(score3S)/float(len3p)
                else:
                    scoreP = -float(score5P)/float(len5p)
                    scoreS = -float(score5S)/float(len5p)

                print "\t_propose_atg::", scoreP + scoreS + tss.pssm_score, tss 

                # append overall score and tssobj to list ONLY WHEN scoreP + scoreS > 0
                if scoreP + scoreS > 0:
                    score2tssobj.append( ( scoreP + scoreS + tss.pssm_score, tss ) )
                #print "\t", tss, tsspos, scoreP, scoreS, (len5p,len3p), (score5P,score3P),(score5S,score3S),
                #print sum(firstexon_array_algpresence), sum(firstexon_array_algsimilarity), len(firstexon_array_algsimilarity)

            if score2tssobj:
                # get highest scoring TSS
                score2tssobj.sort()
                score2tssobj.reverse()
                atg = score2tssobj[0][1]
            else:
                # no sensible Methionine/StartSite
                atg = False
        elif orfObj._tss_sites:   
            # no inwpCBG applied -> return first tss in line
            atg = orfObj._tss_sites[0]                 
        else:
            # get first Methionine coordinate
            atg = orfObj.potential_start_nt_positions()[0]
        # return this TSS object or ATG coordinate
        return atg
    else:
        # no methionines -> Orf can never become first one...
        return False

# end of function _propose_atg_on_orf



def _select_known_exons_to_include_in_abfgp_genemodel(known_exons,abfgp_genemodel):
    """ """
    # create a Set of DNA positions that are already covered by introns/exons
    assigned_coords = Set()
    for elem in abfgp_genemodel:
        if elem: assigned_coords.update(elem.dna_range())

    # add (small) Exons which could not be ABFGP confirmed
    # TODO TODO TODO: this code must be extended much more!!!

    add_exon_status = []
    for mainpos in range(0,len(known_exons)):
        exon = known_exons[mainpos]
        IS_SMALL_EXON = exon.length <= 40
        IS_FIRST_SMALL_EXON = exon.length <= 100
        IS_FINAL_SMALL_EXON = exon.length <= 100
        ADD_EXON = False
        if exon.end < min(assigned_coords):
            if abfgp_genemodel[0].acceptor.__class__.__name__ in ['SpliceAcceptor','SpliceAcceptorAG']:
                if exon.IS_FIRST and IS_FIRST_SMALL_EXON:
                    ADD_EXON = True
                elif IS_SMALL_EXON:
                    # a second,third,... tiny exon 5' of the evidence
                    ADD_EXON = True
                else:
                    # large exon 5' of the evidence -> do NOT include
                    pass
            else:
                # current ABFGP first exon has a StartCodon/Methionine
                # so, do NOT add this annotated exon in front of it
                pass
        elif exon.start > max(assigned_coords):
            if exon.IS_FINAL and IS_FINAL_SMALL_EXON:
                ADD_EXON = True
            elif IS_SMALL_EXON:
                # a second,third,... tiny exon 3' of the evidence
                ADD_EXON = True
            else:
                # large exon 3' of the evidence -> do NOT include
                pass
        elif not assigned_coords.intersection(exon.dna_range()):
            if IS_SMALL_EXON:
                # annotated exon in the middle of the gene structure
                ADD_EXON = True
            else:
                ADD_EXON = False
        else:
            # exon overlaps with an area which is `closed` by exons & introns
            ADD_EXON = None

        # append statud to add_exon_status 
        add_exon_status.append( ADD_EXON )

    # now have a closer look at add_exon_status
    # [ None, None, True, None, None ] -> central added exon, oke
    # [ True, None, ...., None, None ] -> 5' added exon FirstExon, oke
    # [ None, None, ...., None, True ] -> 3' added exon FinalExon, oke
    # [ True, False, None, None, ... ] -> 5' large exon mist -> do not include!

    # fix additions on the 5' side
    for pos in range(0,len(add_exon_status)-1):
        this = add_exon_status[pos]
        if this != True: continue
        allow_insertion = True
        for nextpos in range(pos+1,len(add_exon_status)):
            next = add_exon_status[nextpos]
            if next == None:
                break
            if next == False:
                allow_insertion = False
                break
        if not allow_insertion:
            # small 5' exon followed by large exon without similarity.
            # overrule addition to False
            add_exon_status[pos] = False

    # fix additions on the 3' side
    for pos in range(len(add_exon_status)-1,0,-1,):
        this = add_exon_status[pos]
        if this != True: continue
        allow_insertion = True
        for prevpos in range(pos-1,-1,-1):
            prev = add_exon_status[prevpos]
            if prev == None:
                break
            if prev == False:
                allow_insertion = False
                break
        if not allow_insertion:
            # small 3' exon preceded by large exon without similarity.
            # overrule addition to False
            add_exon_status[pos] = False

    # translate add_exon_status & known_exons to exon list
    return_addition_exon_list = []
    for pos in range(0,len(add_exon_status)):
        if add_exon_status[pos] == True:
            return_addition_exon_list.append( known_exons[pos] )

    # return return_addition_exon_list
    return return_addition_exon_list

# end of function _select_known_exons_to_include_in_abfgp_genemodel
    

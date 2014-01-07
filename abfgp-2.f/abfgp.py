#!/usr/bin/python
""" Main executable of the Alignment Based Fugal Gene Prediction (ABFGP) Method """
__authors__ = "Ate van der Burgt"
__license__ = "MIT"

# Imports (Python)
from copy import deepcopy
from os.path import (
    isdir  as osPathIsdir,
    join   as osPathJoin,
    exists as osPathExists,
    )
from os import (
    remove as osRemove,
    mkdir  as osMkdir,
    chdir  as osChdir,
    system as osSystem,
    )

    
# Imports (ABGP)
from abgpgenelocusdirectory import (
    AbgpGeneLocusDirectory,
    make_abgpgenelocusdirectory_from_fasta,
    )
from abgp_geneconfirmation import (
    readsequences,
    rungetorf,
    parseinputgff,
    annotatedgeneexonsizeevaluation,
    ntracksindnasequencecheck,
    confirmcanonicalsplicesites,
    )
from abgp_unigeneconfirmation import geneandunigeneconfirmation

from abgp_etc import (
    _blastdb_cleanup,
    _file_cleanup,
    _abgp_safe_identifier,
    recreate_global_informantdata,
    abgpsysexit,
    )
    
# Imports (PACB)
import pacb
from pacb.connecting.orfs import get_potention_first_exons_on_orf
from pacb.connecting.mapping import _filter_aligned_introns_on_pssm_entropy_combination
from pacb.conversion import pacbp_from_clustalw, pacbp2pacbporf

# Imports (graphPlus / graphAbgp)
import graphPlus
import graphAbgp
from graphAbgp import (
    PacbpCollectionGraph,
    GeneTreeGraph,
    GenestructureOfCodingBlockGraphs,
    InwardsPointingCodingBlockGraph,
    OrganismNotPresentInGraph,
    )
from graphAbgp.graph_pacbpcollection import _delete_pacbp

# Pythonlibs imports
from pythonlibs.stdoutmanaging import stdoutManagerClass

from gene.orf import TcodeCodingOrf

# Import form all other  libs
from lib_orfset import OrfSet
from lib_tinyexon_pairwise import (
    discover_tinyexons_simple,
    discover_tinyexons_complex,
    tinyexonpacbporfs2InwardsPointingCodingBlockGraphs,
    )

from lib_synteny_pairwise import (
    detect_and_remove_synteny,
    get_first_and_final_inwpcbg_pos,
    )
from lib_sequenceerror import (
    apply_indel_position_correction,
    check_and_cleanup_sequenceerror_interfaces,
    check_intron_interfaces_for_hidden_sequenceerrors,
    )
from lib_utrornongene_pairwise import (
    detect_and_remove_gtgdiscrepancy,
    detect_and_remove_utrornonegene_inwpcbgs,
    detect_and_remove_single_nonfinal_inwpcbg,
    detect_and_remove_single_nonfirst_inwpcbg,
    )

from lib_introns_pairwise import (
    _update_observed_introns,
    _calc_aligned_score,
    _get_main_interface,
    _filter_introns,
    _filter_introns_for_exon_msr,
    _create_failed_intron_gff,
    )

from lib_proteinsimilaritymatrix import ProteinSimilarityMatrix
from lib_crossdatafunctions import *
from lib_crossblast import (
    iterativecrossblastp,
    createblastdbs,
    )
from lib_tcode import obtaintcodedata
from lib_signalp import (
    obtainlocussignalpdata,
    predict_signalp_exons,
    update_PCG_with_signalpexons,
    )

    
from lib_sequencerepetitiveness import annotatedproteinsequencerepetitivenesscheck
from lib_stopwatch import StopWatch
from lib_fasta import parseSingleFasta, writeMultiFasta
from lib_pcg2blocks import *
from lib_abgpgff import (
    create_intron_mapping_gff_file,
    gff2database, dbcleanup,
    )
from lib_codingarray import (
    PCG2codingarray,
    PCG2similarityarray,
    pacbporflist2codingarray,
    )
from lib_clustalw import clustalw



# HMM Imports
from parsers.hmm import hmmsearch_protein
from lib_hmm_pairwise import (
    _create_hmm_db,
    _create_hmm_profile,
    is_hmmpacbporf_conflicting_with_pacbporflist,
    hmmresults2splittedpacbps,
    hmmpacbporf2PCG
    )

# import settings
from settings import *
from settings.dbwarehouse import _get_organism_full_name

# Import optparse functions and validators
from lib_optparse import (
    abgpoptparser,
    generaloptions,     validate_generaloptions,
    abgpinputoptions,   validate_abgpinputoptions,
    abgppairwisespecialoptions,
    abgppairwiseinformantoptions,
    abgpoutputoptions,  validate_abgpoutputoptions,
    abgpunigeneoptions, validate_abgpunigeneoptions,
    )

################################################################################
# initialize and start the StopWatches
################################################################################
# stwMS will time the MileStones in the algorithm flow
stwMS = StopWatch(name='MileStoneSTW')
# stw   will time all steps in the algorithm flow
stw = StopWatch()
stwMS.start()
stw.start()
stwMS.stepname = "DATALOADING"
stw.stepname = "DATALOADING"

################################################################################
# construct the command line option parser
################################################################################
parser = abgpoptparser()
generaloptions(parser)
abgpinputoptions(parser)
abgppairwisespecialoptions(parser)
abgppairwiseinformantoptions(parser)
abgpoutputoptions(parser)
abgpunigeneoptions(parser)

################################################################################
# parse the command line & validate
################################################################################
(OPTIONS, args) = parser.parse_args()
validate_generaloptions(parser,OPTIONS)
validate_abgpinputoptions(parser,OPTIONS)
validate_abgpoutputoptions(parser,OPTIONS)

# main logging variables
VERBOSE = OPTIONS.verbose
QUIET   = OPTIONS.quiet
SILENT  = OPTIONS.silent

################################################################################
# Helper functions for ABGP
################################################################################

def logf(*args):
    """
    Log results of intermediate steps to STDOUT
    Change this function when more/less information is required/needed
    Major milestone messages start with a '##' symbol
    Milestones messages start with a '#' symbol
    """
    if SILENT: pass
    else:
        x = " ".join([str(item) for item in args])
        if VERBOSE:
            print x
        else:
            if x[0:2] == "##":  print x
            elif x[0] == "#" and x[1] != "#":
                if not QUIET: print x
            else:           pass

# end of function logf


def removeinformantfrominput(informant,input,message=""):
    """ """
    del( input[informant] )
    logf("# INFORMANT DELETED:", informant, message )
    return True

# end of function removeinformantfrominput

################################################################################
# Import & get the blastp Matrices object
################################################################################
from lib_proteinsimilaritymatrix import ProteinSimilarityMatrix
DEFAULT_MATRIX             = ProteinSimilarityMatrix(name=BLASTP_MATRIX_NAME)
BLASTP_SIMILARITY_MATRIX   = DEFAULT_MATRIX
TINYEXON_MATRIX            = ProteinSimilarityMatrix(fname=TINYEXON_MATRIX_PATH)



# init stdoutManager; it will capture all
# STDOUT messages from deeper modules in non-VERBOSE mode
stdoutManager = stdoutManagerClass()

################################################################################
# create outdir for this project/run
################################################################################
if not osPathIsdir(OPTIONS.outdir):
    osMkdir(OPTIONS.outdir)
# chdir to this directory
osChdir(OPTIONS.outdir)

################################################################################
# Parse input datas from/as AbgpGeneLoci 
################################################################################
target_is_seen = False
input = {} # input data structure, empty

if OPTIONS.multifasta:
    for hdr, sequence in OPTIONS.multifasta.iteritems():
        OPTIONS.loci.append( make_abgpgenelocusdirectory_from_fasta(hdr,sequence,OPTIONS.outdir) )
    
if OPTIONS.dnafiles:
    for fname in OPTIONS.dnafiles:
        hdr,sequence,description = parseSingleFasta(open(fname))
        OPTIONS.loci.append( make_abgpgenelocusdirectory_from_fasta(hdr,sequence,OPTIONS.outdir) )

if OPTIONS.loci:
    locus_added = 0
    for locusdir in OPTIONS.loci:
        if target_is_seen and OPTIONS.shortcut and locus_added > OPTIONS.minimal_num_loci:
            # in shortcut mode, only a maximum of OPTIONS.minimal_num_loci
            # are used (means OPTIONS.minimal_num_loci -1 informants)
            continue
            
        if locusdir.__class__.__name__ == "AbgpGeneLocusDirectory":
            # coming from fasta entry
            locus = locusdir
        else:
            locus = AbgpGeneLocusDirectory(locusdir)
            
        for line in locus.filestats().split('\n'): logf(line)
        locuskey = locus._create_auto_key(identifier2organism=ABGP_LOCUS_IDENTIFIER2ORGANISM_MAPPING)

        # check omit_informants, required_informants and maximal_num_loci criterion
        if OPTIONS.target == None or OPTIONS.target not in [ locuskey, locus.protein_fref() ]:
            # check for required informants
            if locuskey in OPTIONS.required_informants:
                pass
            if locus.protein_fref() in OPTIONS.required_informants:
                pass
            # check for omit informants
            elif locuskey in OPTIONS.omit_informants:
                logf("# informant locus removed; in omit_informants", (locuskey,locus.protein_fref()))
                continue
            elif locus.protein_fref() in OPTIONS.omit_informants:
                logf("# informant locus removed; in omit_informants", (locuskey,locus.protein_fref()))
                continue
            # check maximal_num_loci criterion; we assume that
            # loci are listed in potential interesting order.
            # This is definately a FALSE hypothesis, but usually.
            # maximal_num_loci is applied for speedup, so we accept
            # potential loss in accuracy here.
            # In the application as published in the manuscript,
            # --filewithloci were used as input, which were ordered
            # based on (estimated) similarity to the target locus.
            # A such, --maximal_num_loci=15 discarded only the
            # more distantly related informants
            elif locus_added >= OPTIONS.maximal_num_loci:
                logf("# informant locus removed; >= maximal_num_loci", (locuskey,locus.protein_fref()))
                continue
            else:
                # fine; just add as an informant
                pass
        
        else:
            # Yes, we encountered the target locus
            target_is_seen = True
            # log the encountered target
            logf("# TARGET gene locus:", OPTIONS.target )

        if input.has_key(locuskey):
            for suffix in list('abcdefghijklmnopqrstuvwxyz'):
                if not input.has_key(locuskey+suffix):
                    locuskey = locuskey+suffix
                    # check omit_informants criterion
                    if locuskey != OPTIONS.target:
                        if locuskey in OPTIONS.omit_informants: continue
                        # if not -> update input data structure
                        input.update( locus.toinputdict(key=locuskey) )
                        break
        else:
            input.update( locus.toinputdict(key=locuskey) )

        # store the input (sub)dict to locus.input as a pointer
        locus.input = input[locuskey]
        locus_added += 1
        logf("gene locus directory read:", locus.fref)

        
################################################################################
# Start logging! MAIN ABFGP MILESTONE
################################################################################
logf("##", stwMS.lap(), "OUTDIR:", OPTIONS.outdir)

################################################################################
# read input sequence files
################################################################################
input = readsequences(input)

################################################################################
# check if omit_informants is applied
################################################################################
if OPTIONS.omit_informants:
    for informant_org in Set(OPTIONS.omit_informants):
        if informant_org != OPTIONS.target and informant_org in input.keys():
            # delete this informant!
            del( input[informant_org] )
            logf("# --omit_informants:", informant_org)
if len(input) < OPTIONS.minimal_num_loci:
    message = "## %s < %s (--minimal_num_loci)" % (
        len(input), OPTIONS.minimal_num_loci )
    # done -> abgpsysexit()
    logf(message)
    abgpsysexit(input,OPTIONS,message=message)

################################################################################
# check if target name is recognized & specific enough when an abreviation 
################################################################################
if OPTIONS.target:
    if OPTIONS.target in ABGP_LOCUS_IDENTIFIER2ORGANISM_MAPPING.values() and\
    OPTIONS.target+"a" in input.keys():
        informant_names = set()
        informant_names.update(input.keys())
        informant_names.update([ vdict['proteinfref'] for vdict in input.values() ])
        informant_names.update([ vdict['locusfref'] for vdict in input.values() ])
        informant_names = list(informant_names)
        informant_names.sort()
        if None in informant_names: informant_names.remove(None)
        message = [ "## --target %s is not specific enough (paralogs present!)" % OPTIONS.target ]
        for i in range(0,len(informant_names),5):
            message.extend( "##   use one of %s" % (", ".join(informant_names[i:i+5])) )
        message = "\n".join(message)
        # done -> abgpsysexit()
        logf(message)
        abgpsysexit(input,OPTIONS,message=message)
else:
    # --target is not applied
    thisoption = parser.get_option('--target')
    print "Exception.TargetNotApplied", "--target"
    print "--target:", thisoption.help
    informant_names = set()
    informant_names.update(input.keys())
    informant_names.update([ vdict['proteinfref'] for vdict in input.values() ])
    informant_names.update([ vdict['locusfref'] for vdict in input.values() ])
    informant_names = list(informant_names)
    informant_names.sort()
    if None in informant_names: informant_names.remove(None)
    for i in range(0,len(informant_names),5):
        print "## use one of %s" % (", ".join(informant_names[i:i+5]))
        
    # done -> abgpsysexit()
    message = "## --target not applied"
    logf(message)
    abgpsysexit(input,OPTIONS,message=message)

 
if OPTIONS.target and OPTIONS.target not in input.keys():
    # check for auto-dbwarehouse renaming of gene -> org identifier
    target_key = None
    target_key_new = None
    for key in input.keys():
        if OPTIONS.target in [ input[key]['proteinfref'], input[key]['FREF'] ]:
            target_key = key
            if target_key[0:-1] in ABGP_LOCUS_IDENTIFIER2ORGANISM_MAPPING.values() and\
            target_key[0:-1].upper() not in input.keys():
                # paralogous proteins ; ORG[a-z] as target
                target_key_new = target_key[0:-1].upper()
            elif target_key in ABGP_LOCUS_IDENTIFIER2ORGANISM_MAPPING.values():
                # paralogous proteins ; full protein identifier --target
                target_key_new = target_key.upper()
            else:
                # other case -> no hard-assignment of target_key_new
                pass
            # --target recognized -> break out, done here
            break

    if target_key:

        if target_key_new:
            input[target_key_new] = input[target_key]
            del( input[target_key] )
            OPTIONS.target = target_key_new
        elif target_key.upper() not in input.keys():
            input[target_key.upper()] = input[target_key]
            del( input[target_key] )
            OPTIONS.target = target_key.upper()
        elif target_key.lower() not in input.keys():
            input[target_key.lower()] = input[target_key]
            del( input[target_key] )
            OPTIONS.target = target_key.lower() 
        elif 'TARGET' not in input.keys():
            input['TARGET'] = input[target_key]
            del( input[target_key] )
            OPTIONS.target = 'TARGET'
        else:
            raise "UNIQUE --TARGET ASIGNMENT FAILED: %s %s" % ( target_key, input.keys() )
    else:
        # error is raised directly below here..
        pass

if OPTIONS.target and OPTIONS.target not in input.keys():
    message = "## target identifier (%s) not recognized in: %s" % (
        OPTIONS.target, input.keys() )
    # done -> abgpsysexit()
    logf(message)
    abgpsysexit(input,OPTIONS,message=message)

    
logf("number of genes: %s: %s" % (len(input.keys()),input.keys()))
logf("mapping of shortcut/loci FREFs to protein and dna identifiers")
for key in input.keys():
    logf(key, "\t", input[key]['proteinfref'], "\t",
            input[key]['fstrand'], "\t", input[key]['locusfref'],
            "\t", input[key]['genomefref'])

logf("#", stw.lap(),"(commandline) input processed")


################################################################################
# check for request of pre-visualisation of gene locus
################################################################################
if OPTIONS.output_storelocitoggbdb and OPTIONS.target:
    for orgkey in input.keys():
        if orgkey == OPTIONS.target or\
        input[orgkey]['gldobj'].protein_fref() == OPTIONS.target:
            locus = input[orgkey]['gldobj']
            FREF  = "locus_"+locus.protein_fref()

            # store the locus to GGB database     
            locus._gffdata.extend( input[orgkey]['gffs'] )
            locus.storetoggbdb(fref="locus_"+locus.protein_fref())
            logf("# LOCUS STORED TO GGB: ", FREF )


################################################################################
# proces getorf data
################################################################################
input = rungetorf(input)
for key in input.keys():
    logf("number of orfs %s for gene %s" % ( len(input[key]['orfs'].orfs),key))
logf("#", stw.lap(),"getorf processed")


################################################################################
# proces input gff data in the AbgpGeneLocusDirectories (when applied)
################################################################################
input = parseinputgff(input)
logf("#", stw.lap(),"input gff processed")

if VERBOSE:
    for org in input.keys():
        logf(org, 'gff-gene', len(input[org]['gff-gene']))
        #if input[org]['gff-unigene']:
        #    for track in input[org]['gff-unigene']:
        #        logf(track)

################################################################################
# confirm gene & unigene structure on Orfs
################################################################################
# Here ALWAYS verbose=True to obtain (crucial) Warning messages:
# - UniGene not mappable on ORFs
# - Gene not mappable on ORFs 
input, gene_confirmation, unigene_confirmation = geneandunigeneconfirmation(
        input,verbose=VERBOSE)
logf("#", stw.lap(), gene_confirmation, unigene_confirmation,
        "gene and unigene structure(s) confirmed on Orfs")


################################################################################
# perform validate_abgpunigenebenchmarkoptions()
################################################################################
status = validate_abgpunigeneoptions(parser,OPTIONS,input)


################################################################################
# do several sequence / annotation property checks
################################################################################
IS_REPETITIVE           = annotatedproteinsequencerepetitivenesscheck(input)
HAS_SMALL_OR_TINY_EXONS = annotatedgeneexonsizeevaluation(input)
HAS_N_SYMBOLS           = ntracksindnasequencecheck(input)
for org in input.keys():
    status,warnings = confirmcanonicalsplicesites(input[org]['genomeseq'],
        input[org]['gff-gene'],exon_fmethod=GFF_CDS_FMETHOD,verbose=VERBOSE)
    input[org]['warnings'].extend(warnings)

################################################################################
# print/evaluate all the warnings
################################################################################
ORGANISMS = input.keys()
for org in ORGANISMS:
    if input[org]['warnings']:
        logf("#",org,"WARNINGS (%s):" % len(input[org]['warnings']))
        for warn in input[org]['warnings']:
            logf("#\t", warn)
        fatal_warnings = [
            'GeneStructureIsNotMappableOnOrfsWarning',
            #'IncompleteGeneStructureWarning',
            ]
        for fatal_warn in fatal_warnings:
            if fatal_warn in\
            [ warn.__class__.__name__ for warn in input[org]['warnings'] ]:
                # remove this organism/gene, because it is likely
                # to result in bogus predictions
                del( input[org] )
                logf("##",org,"REMOVED:", fatal_warn)
                break


################################################################################
# Check if OPTIONS.minimal_num_loci is still valid
################################################################################
if len(input) < OPTIONS.minimal_num_loci or\
(OPTIONS.target and OPTIONS.target not in input.keys()):
    message = "## %s < %s (--minimal_num_loci) or target (%s) identifier removed" % (
        len(input), OPTIONS.minimal_num_loci, OPTIONS.target
        )
    # done -> abgpsysexit()
    logf(message)
    abgpsysexit(input,OPTIONS,message=message)


################################################################################
# process signalp data
################################################################################
cnt = obtainlocussignalpdata(input)
logf("#", stw.lap(),"SignalP data loaded (%s)" % cnt)


################################################################################
# add unigene(s) to the list of informants
################################################################################
ORGANISMS = input.keys()
unigene_informants_added = []
for org in ORGANISMS:
    if input[org]['orfid-unigenestructure']:
        if OPTIONS.unigenebenchmark and org == OPTIONS.target:
            # do not add unigene of target as an informant (benchmark mode)
            continue
        if OPTIONS.omit_unigene_informants:
            # omit_unigene_informants -> do not add unigenes as informants
            continue
        gldObj = input[org]['gldobj']
        ugannotation = gldObj.unigeneannotation()
        if ugannotation[0] != 'nonc':
            # make single-orf informant from unigene
            cdsseq          = gldObj.get_unigene_cds_sequence()
            protseq, ugframe= gldObj.get_unigene_protein_sequence()
            # correct cdsseq for unigene's frame (in case of utr3p, fragm)

            if ( ugannotation[0], protseq, ugframe ) == ('fragm', None, None ):
                # fragment of a unigene without clear frame
                # ignore this one here as a informant!
                continue
            elif ( protseq, ugframe ) == ( None, None ):
                # most likely a (very) short unigene of other type
                # then 'fragment'. In any case, unigene protein sequence
                # could not be determined. Skip here
                continue
            else:
                pass

            cdsseq = cdsseq[ugframe:]
            # correct protein & cds sequence for stop codon
            if protseq[-1] == '*':
                protseq = protseq[0:-1]
                cdsseq  = cdsseq[0:-3]

            # correct (possible) tailing nucleotides for utr5p UniGenes
            if ugannotation[0] == 'utr5p' and len(cdsseq) % 3 != 0:
                # corect the end of the unigene sequence
                cdsseq = cdsseq[0:-(len(cdsseq) % 3)]
                if len(cdsseq)/3 < len(protseq):
                    protseq = protseq[0:(len(cdsseq)/3)]

            # create CodingOrf object
            orfobj = TcodeCodingOrf("dummyid",protseq,inputgenomicsequence=cdsseq,force_id=1,start=1,end=len(cdsseq))
            # label this Orf as a unigene
            setattr(orfobj,ORF_IS_UNIGENE_LABEL,True)
            orfSetObj = OrfSet()
            orfSetObj.orfs = [ orfobj ]

            # make recognizable and unique (and ABGP safe) unigene identifier
            ugfref = "%s[%s]" % (_abgp_safe_identifier(gldObj.unigene_fref()),org)

            uginputdict = {
                # non-relevant but required GeneLocus attributes
                'dirname':              None,
                'locusfile':            None,
                'proteinseqfile':       None,
                'genomeseqfile':        None,
                'unigenefile':          None,
                'genefile':             None,
                'genomefref':           None, 
                'locusfref':            None,

                # non-relevant but required attributes
                'gffs'                  : [],
                'warnings'              : [],
                'METHOD'                : "",
                'fstrand'               : '+',
                'orfid-unigenestructure': [],

                # required attibutes
                'genomeseq'             : cdsseq,
                'proteinseq'            : protseq,
                'orfs'                  : orfSetObj,
                'tcode'                 : [],
                'blastdb'               : {},
                'FREF'                  : ugfref,
                'proteinfref'           : ugfref,

                # mark this as an unigene
                'is_unigene'            : True,
                'is_unigene_of'         : org,
                'orfid-genestructure'   : [1],
                'unigene-annotation'    : ugannotation[0],

                }
            # add to input dict as an informant!
            input[ugfref] = uginputdict
            unigene_informants_added.append(ugfref)

# log message for unigene_informants_added
logf("#", stw.lap(),"Unigenes added as informants:", unigene_informants_added)

################################################################################
# proces tcode data
################################################################################
input = obtaintcodedata(input)
logf("#", stw.lap(),"tcode processed")

################################################################################
# check if OPTIONS.onlyusegenestructureorfs is applied;
# if so, HARD-REMOVE all orfs except the coding orfs for
# the applied genomes/genes. This means a HUGE speedup!
################################################################################
if OPTIONS.onlyusegenestructureorfs:
    # create global informant data structures
    (   ORGANISMS,
        GENECOMBIS,
        GENE_INFORMANT_SET,
        GENE_IDENTIFIER_SET,
        UNIGENE_INFORMANT_SET )=recreate_global_informantdata(input,OPTIONS.target)
    if OPTIONS.onlyusegenestructureorfs == ['all']:
        # remove non-annotated Orfs for target and ALL informant
        OPTIONS.onlyusegenestructureorfs = ORGANISMS
    if OPTIONS.onlyusegenestructureorfs == ['informants']:
        # remove non-annotated Orfs for ALL informant (but not target)
        OPTIONS.onlyusegenestructureorfs = GENE_INFORMANT_SET
    if OPTIONS.onlyusegenestructureorfs == ['target']:
        # remove non-annotated Orfs of TAREGT only
        OPTIONS.onlyusegenestructureorfs = [ OPTIONS.target ]
    # count current total number of Orfs
    current_orf_count = sum([ len(input[org]['orfs'].orfs) for org in input.keys() ])

    for org in OPTIONS.onlyusegenestructureorfs:
        if input.has_key(org):
            if input[org]['orfid-genestructure']:
                # get the orfs listed in orfid-genestructure
                orfsublist = []
                for orfid in input[org]['orfid-genestructure']:
                    orfsublist.append(
                        input[org]['orfs'].get_orf_by_id(orfid)
                        )
                # and replace this orflist by this sublist
                input[org]['orfs'].orfs = orfsublist
                logf(org, "\tapplied in --onlyusegenestructureorfs:",
                    input[org]['orfid-genestructure'], "remaining")
            else:
                logf(org, "\tapplied in --onlyusegenestructureorfs but",
                    "given annotation not correctly placeable on DNA sequence")
        else:
            logf(org, "\tapplied in --onlyusegenestructureorfs but not",
                    "recognized as a gene/organism/identifier")

    # count new total number of Orfs
    new_orf_count = sum([ len(input[org]['orfs'].orfs) for org in input.keys() ])
    # logf message
    logf("#", stw.lap(), "--onlyusegenestructureorfs: decreased from %s to %s" % (
            current_orf_count, new_orf_count ) )


################################################################################
# create auxiliary lists of (unique) organisms and cross combinations
################################################################################
ORGANISMS = input.keys()
ORGANISMS.sort()
crossdata = createcrossdata(input)
logf("# STEP0.F done;", stw.lap(),"crossdata and other variables created")


################################################################################
# some printing in VERBOSE mode
################################################################################
if VERBOSE:
    for org in ORGANISMS:
        coding,noncoding,other = 0,0,0
        for orf in input[org]['orfs'].orfs:
            if orf.tcode_is_coding(): coding+=1
            elif orf.tcode_is_noncoding(): noncoding+=1
            else: other+=1
        logf("%s\t%s orfs, %s-%s-%s (c-o-n)" % (org,len(input[org]['orfs'].orfs),coding,other,noncoding))


################################################################################
# some printing in none-VERBOSE mode
################################################################################
logf("# annotated genestructure model")
for org in ORGANISMS:
    logf("#\t%s\tgene-orfmodel:" % org, input[org]['orfid-genestructure'], input[org]['proteinfref'] )
logf("# Applied UniGene's orfid-structure")
for org in ORGANISMS:
    if input[org]['orfid-unigenestructure']:
        logf("#\t%s\tUniGene-orfmodel:" % org, input[org]['orfid-unigenestructure'])


# Main Graph objects (initialized as None or empty)
GSG     = GenestructureOfCodingBlockGraphs(input)   # empty Genestructure Object
GTG     = GeneTreeGraph()                           # empty GeneTreeGraph

# initialize the (BLOSUM) similarity matrix
from settings.blastp import blastoptions_iter2 as BLASTOPTIONS
# adjust parameter G; detected in CFU_834418
# TODO: how to solve this more generic !?
# TODO: part of a bona fide - high scoring - alignment is simply NOT
# picked up by blastp due to a large following gap.
# Depending on G=10/G=11 or just default settings, large differences are
# reported by blastp
# Q:71-301 S:154-398    &&   Q:10-173 S:63-231      G=11 [ E=3 W=2 M=BLOSUM45 ]
# Q:71-301 S:154-398                                G=10 [ E=3 W=2 M=BLOSUM45 ]
# Q:10-301 S:63-398                                 G=11 [ E=1 W=? M=BLOSUM45 ]
# The frontal part contains a bona fide, high scoring aligment part
BLASTOPTIONS.BLASTP_EXTRA_PARAM_G  = 11
BLASTOPTIONS.extra_blastp_params['G'] = 11
BLASTP_SIMILARITY_MATRIX = ProteinSimilarityMatrix(name=BLASTOPTIONS.BLASTP_MATRIX_NAME)
BLASTOPTIONS.MATRIX      = BLASTP_SIMILARITY_MATRIX

# create crossdata structure
crossdata = createcrossdata(input) 

# recreate auxiliary lists of (unique) organisms
ORGANISMS = input.keys()
ORGANISMS.sort()

# make pairwise gene combinations for the target organism only
# !!WITH THE TARGET ORGANISM AS QUERY!!
GENECOMBIS = [ (OPTIONS.target,org) for org in ORGANISMS ]
GENECOMBIS.remove((OPTIONS.target,OPTIONS.target))

# create global informant data structures
(   ORGANISMS,
    GENECOMBIS,
    GENE_INFORMANT_SET,
    GENE_IDENTIFIER_SET,
    UNIGENE_INFORMANT_SET )=recreate_global_informantdata(input,OPTIONS.target)



# remove combis in crossdata that are not GENECOMBIS
keys = crossdata.keys()
for (a,b) in keys:
    if (b,a) in GENECOMBIS:
        crossdata[(b,a)] = crossdata[(a,b)]
    if (a,b) not in GENECOMBIS:
        del( crossdata[(a,b)] )


logf("# STEP0.G done;", stw.lap(),"crossdata and other variables created")

################################################################################
### MAIN ABFGP MILESTONE
################################################################################
logf("##", stwMS.lap(), "PREPROCESSING DONE", "gene & unigene informants:",
    len(GENE_INFORMANT_SET), len(UNIGENE_INFORMANT_SET) )


################################################################################
### starting next algorithm step: BLASTP
################################################################################
stwMS.stepname = "BLAST"
stw.stepname = "BLAST"


################################################################################
### make blast database for each of the organisms that serve as SBJCT
################################################################################

for org in Set([orgS for orgQ,orgS in GENECOMBIS]):
    # create blast database for this sbjct organism/gene
    cnt = createblastdbs(input,GSG,OPTIONS,organism=org,acceptorfids=input[org]['orfid-genestructure'])
    logf("blast database created: %s (%s orfs)" % ( org, len(input[org]['orfid-genestructure'])) )

logf("#", stw.lap(), "blastdbs INFORMANTS (known) created")

################################################################################
### do crossblasts of KNOWN orfs of OPTIONS.target against KNOWN orfs of INFORMANTS
################################################################################
crossdata, resultcounts = iterativecrossblastp(
            input,crossdata,OPTIONS,
            dbfraction='annotation', # set on 'annotation' -> only known Exons!
            blastoptions=BLASTOPTIONS,
            )

hit_cnt = sum([ len(crossdata[k]['accepted_pacbs']) for k in crossdata.keys() ])
logf("#", stw.lap(), "blastp TARGET (known) vs. INFORMANTS (known):", hit_cnt)

################################################################################
### now turn around the logics; blast all KNOWN orfs of INFORMANTS against
### all orfs (except the known ones) of OPTIONS.target
################################################################################
# To do so, we first need to create reversed data structures.
# `all` Orfs of OPTIONS.target is limited to an area of ~500nt on both sides
# of currently detected Orfs. (Orfs must start or end ~500nt from current hits)
REVERSED_GENECOMBIS = [ (b,a) for (a,b) in GENECOMBIS ]
REVERSED_GENECOMBIS.sort()
reversed_crossdata = createcrossdata(input) 
# remove combis in reversed_crossdata that are not REVERSED_GENECOMBIS
keys = reversed_crossdata.keys()
for (a,b) in keys:
    if (b,a) in REVERSED_GENECOMBIS:
        reversed_crossdata[(b,a)] = reversed_crossdata[(a,b)]
    if (a,b) not in REVERSED_GENECOMBIS:
        del( reversed_crossdata[(a,b)] )


################################################################################
# make blast database for the query; OMIT all known Orfs of the genestructure
################################################################################
# to do so, set OPTIONS.abinitio TEMPORARILY on True, otherwise only the
# *known* Orfs (exactly the ones we wish to remove here ;-) are taken.
applied_abinitio_option = deepcopy(OPTIONS.abinitio)
OPTIONS.abinitio = True

# make a list of OPTIONS.target orfs which fall in the currently predicted
# gene structure (but in other frames)
if input[OPTIONS.target]['orfid-genestructure']:
    start_orf = input[OPTIONS.target]['orfs'].get_orf_by_id(input[OPTIONS.target]['orfid-genestructure'][0])
    end_orf   = input[OPTIONS.target]['orfs'].get_orf_by_id(input[OPTIONS.target]['orfid-genestructure'][-1])
    ### min_orf_end   = start_orf.startPY - 500 # for GeneModelErrors this was TO RESTRICTED
    ### max_orf_start = end_orf.endPY + 500     # for GeneModelErrors this was TO RESTRICTED
    min_orf_end   = None
    max_orf_start = None
else:
    # No known or incorrect input[OPTIONS.target]['orfid-genestructure'];
    # do no filter on Orf positions
    min_orf_end   = None
    max_orf_start = None
# generate a list of accepted OrfIds to make a blastdb from
accepted_orfs_ids = [ orf.id for orf in input[OPTIONS.target]['orfs'].get_elegiable_orfs(
    min_orf_end   = min_orf_end, max_orf_start = max_orf_start,
    rejectorfids = input[OPTIONS.target]['orfid-genestructure'],
    ) ]

# create the blast database
cnt = createblastdbs(input,GSG,OPTIONS,organism=OPTIONS.target,acceptorfids=accepted_orfs_ids)
logf("#", stw.lap(), "blastdb TARGET (all-known) created: (%s) Orfs" % cnt.values()[0])

# reset OPTIONS.abinitio to applied_abinitio_option
OPTIONS.abinitio = applied_abinitio_option

################################################################################
# do crossblasts of KNOWN orfs of INFORMANTS against ( all - known ) orfs of TARGET
################################################################################
reversed_crossdata, resultcounts = iterativecrossblastp(
            input,reversed_crossdata,OPTIONS,
            dbfraction='annotation', # set on 'annotation' -> only known Exons!
            blastoptions=BLASTOPTIONS,
            )

hit_cnt = sum([ len(reversed_crossdata[k]['accepted_pacbs']) for k in reversed_crossdata.keys() ])
logf("#", stw.lap(), "blastp INFORMANTS (known) vs. TARGET (all-known):", hit_cnt)

################################################################################
# do crossblasts of HCP orfs of INFORMANTS against ( all - known ) orfs of TARGET
################################################################################
reversed_crossdata, resultcounts = iterativecrossblastp(
            input,reversed_crossdata,OPTIONS,
            dbfraction='HCP', # set on 'HCP' -> High Coding Potential
            blastoptions=BLASTOPTIONS,
            )

# file cleanup; remove all the blast databases
_blastdb_cleanup(input)

hit_cnt = sum([ len(reversed_crossdata[k]['accepted_pacbs']) for k in reversed_crossdata.keys() ])
logf("#", stw.lap(), "blastp INFORMANTS (HCP) vs. TARGET (all-known):", hit_cnt)

################################################################################
# store reversed blast results to main crossdata by swapping the Pacbps
################################################################################
for (orgS,orgQ) in REVERSED_GENECOMBIS:
    ###print (orgS,orgQ), len(reversed_crossdata[(orgS,orgQ)]['accepted_pacbs'])
    ###print [ k[0] for k in crossdata[(orgQ,orgS)]['accepted_pacbs'].keys() ]
    ###print [ k[0] for k in reversed_crossdata[(orgS,orgQ)]['accepted_pacbs'].keys() ]
    for (a,b,c,d),pacbp in reversed_crossdata[(orgS,orgQ)]['accepted_pacbs'].iteritems():
        rev_key = (a,b,d,c)
        pacbp._swap_query_and_sbjct()
        crossdata[(orgQ,orgS)]['accepted_pacbs'][rev_key] = pacbp

hit_cnt = sum([ len(crossdata[k]['accepted_pacbs']) for k in crossdata.keys() ])
logf("#", stw.lap(), "blastp INFORMANTS vs. TARGET alignments swapped:", hit_cnt)




################################################################################
### make blast database for each of the informants of all Orfs that are not
### linked to a target Orf yet
################################################################################
# to do so, set OPTIONS.abinitio TEMPORARILY on True, otherwise only the
# *known* Orfs (exactly the ones we wish to remove here ;-) are taken.
applied_abinitio_option = deepcopy(OPTIONS.abinitio)
OPTIONS.abinitio = True

for org in Set([orgS for orgQ,orgS in GENECOMBIS]):
    # create blast database for this sbjct organism/gene
    rejectorfids = list(Set([ orfSid for (a,b,orfQid,orfSid) in crossdata[(OPTIONS.target,org)]['accepted_pacbs'].keys() ]))
    detectedOrfObjs = [ input[org]['orfs'].get_orf_by_id(id) for id in rejectorfids ]
    
    # find Orfs which are FULLY included in these Orfs;
    # append those to rejectorfids too
    for orfObj in input[org]['orfs'].orfs:
        if orfObj.id in rejectorfids: continue
        if orfObj.length <= 30:
            rejectorfids.append(orfObj.id)
            continue
        for detectedOrfObj in detectedOrfObjs:
            if orfObj.start >= detectedOrfObj.start and orfObj.end <= detectedOrfObj.end:
                rejectorfids.append(orfObj.id)
                break

    cnt = createblastdbs(input,GSG,OPTIONS,organism=org,
            #rejectorfids = input[OPTIONS.target]['orfid-genestructure'],
            rejectorfids = rejectorfids,
            )
    logf("# blast database created: %s (%s/%s orfs)" % ( org, cnt.values()[0], len(input[org]['orfs'].orfs) ) )

logf("#", stw.lap(), "blastdbs INFORMANTS (all-detected) created")

# reset OPTIONS.abinitio to applied_abinitio_option
OPTIONS.abinitio = applied_abinitio_option

################################################################################
# do crossblasts of detected orfs of TARGET against ( all - detected ) orfs of INFORMANT
################################################################################
crossdata, resultcounts = iterativecrossblastp(
            input,crossdata,OPTIONS,
            dbfraction='iterative', # set on 'iterative' -> All target Orfids
            blastoptions=BLASTOPTIONS,
            )

hit_cnt = sum([ len(crossdata[k]['accepted_pacbs']) for k in crossdata.keys() ])
logf("#", stw.lap(), "blastp TARGET (detected) vs. TARGET (all-detected):", hit_cnt)

################################################################################
### Blastp steps done; cleanup and MAIN ABFGP MILESTONE
################################################################################
# file cleanup; remove all the blast databases
_blastdb_cleanup(input)
logf("##", stwMS.lap(), "BLASTP steps performed")

################################################################################
### starting next algorithm step: crossdata2PCG
################################################################################
if not VERBOSE: stdoutManager.buffer() # capture all STDOUT messages from deeper modules 
stwMS.stepname = "PCG"
stw.stepname = "PCG"



################################################################################
# Create a PacbpCollectionGraph from crossdata dict structure
# and perform several operations on it for Pacbp filtering.
# Steps are required to translate the local alignments between Orfs
# to a global (gene/exon) alignment
################################################################################

PCG = PacbpCollectionGraph(crossdata=crossdata,
        blastmatrix=BLASTP_SIMILARITY_MATRIX)
logf("#",stw.lap(),PCG,"PCG.__init__()")

PCG.split_pacbps_on_gapsize()
logf(stw.lap(),PCG,"PCG.split_pacbps_on_gapsize()")

PCG.split_lowscoring_pacbps_on_gapsize()
logf(stw.lap(),PCG,"PCG.split_lowscoring_pacbps_on_gapsize()")

PCG.recover_lowscoring_pacbps()
logf(stw.lap(),PCG,"PCG.recover_lowscoring_pacbps()")

PCG.harvest_pacbps_from_crossdata()
logf(stw.lap(),PCG,"PCG.harvest_pacbps_from_crossdata()")

PCG.remove_nonlinear_pacbs()
logf(stw.lap(),PCG,"PCG.remove_nonlinear_pacbs()")

# TEMPORARILY make PacbPORFS of PacbPs
PCG.pacbps2pacbporfs(input)
logf(stw.lap(),PCG,"PCG.pacbps2pacbporfs(input)")

PCG.remove_inclusive_pacbps()
logf(stw.lap(),PCG,"PCG.remove_inclusive_pacbps()")

status = PCG.remove_noncoding_pacbpdnas(verbose=True)
logf(stw.lap(),PCG,status,"PCG.remove_noncoding_pacbpdnas()")

status = PCG.remove_lowerasexpected_ntidentity_pacbps(GENECOMBIS,input)
logf(stw.lap(),PCG,status,"PCG.remove_lowerasexpected_ntidentity_pacbps()")

# TEMPORARILY reset PacbPORFS back to PacbPs
# functions remove_alternative_pacbps(), correct_overlaps()
# remove_repetitive_pacbps() require PacbPs, not PacbPORFs ...
PCG.pacbporfs2pacbps()

PCG.remove_alternative_pacbps(overlap_ratio=0.3) # NOT 0.8...
logf(stw.lap(),PCG,"PCG.remove_alternative_pacbps(overlap_ratio=0.3)")

PCG.correct_overlaps(GENECOMBIS)
logf(stw.lap(),PCG,"PCG.correct_overlaps()")

PCG.remove_repetitive_pacbps(overlap_ratio=0.85) # not 0.8
logf(stw.lap(),PCG,"PCG.remove_repetitive_pacbps(overlap_ratio=0.85)")

PCG.remove_short_pacbps()
logf(stw.lap(),PCG,"PCG.remove_short_pacbps()")

# make EXTENDED pacbporfs from pacbps
PCG.extend_pacbporfs(input)
logf(stw.lap(),PCG,"PCG.extend_pacbporfs(input)")

# extend_matched_ends of the PacbPORFS (caused by pacbp overlap correction)
for pacbporf in PCG.pacbps.values(): pacbporf.extend_matched_ends()
logf(stw.lap(),PCG,"PCG->pacbporf.extend_matched_ends()")

# merge pacbporfs which are probably interrupted by poor identity stretches
PCG_HAS_HIGH_GAP_RATIO_PACBPORFS = PCG.merge_high_gap_ratio_pacbporfs(OPTIONS.target,verbose=True)
logf(stw.lap(),PCG,"PCG.merge_high_gap_ratio_pacbporfs()",PCG_HAS_HIGH_GAP_RATIO_PACBPORFS)


################################################################################
### cleanup nodes without edges
################################################################################
removed = PCG.remove_edgeless_nodes()
if removed:
    # recreate data structures
    ( ORGANISMS, GENECOMBIS,
      GENE_INFORMANT_SET,
      GENE_IDENTIFIER_SET,
      UNIGENE_INFORMANT_SET) = recreate_global_informantdata(PCG,OPTIONS.target)

    # check if an informant is not in ORGANISMS but still in input
    for informant in Set(input.keys()).difference(ORGANISMS):
        status = removeinformantfrominput(informant,input)

logf(stw.lap(), PCG, "nodes without edges removed:", removed)

################################################################################
# PCG operations done. MileStone logging
################################################################################
if not VERBOSE: stdoutManager.reset()
logf("##", stwMS.lap(), PCG, "PCG steps performed")

# print structure in all modes
inwpcbgs = PCG2inwpCBGS(PCG)
if not True in (QUIET,SILENT): print_inwpcbgstructure(inwpcbgs,ORGANISMS)

################################################################################
### starting next algorithm step: QuarantineHighNtIdentityInformants
################################################################################
if not VERBOSE: stdoutManager.buffer() # capture all STDOUT messages from deeper modules 
stwMS.stepname = "QuarantineHighNtIdentityInformants"
stw.stepname = "QuarantineHighNtIdentityInformants"


################################################################################
### Create GeneTreeGraph from PCG
################################################################################
GTG = pcg2gtg_by_identity(PCG,OPTIONS.target,identifier_list=GENE_INFORMANT_SET)
logf("#", stw.lap(), GTG, "created")


################################################################################
### Detect and high-nt-identity informants in GeneTreeGraph and quarantine them
################################################################################
nt_identity_quarantined_input = {}
quarantined_informant_nt_identities = []

for informant in GENE_INFORMANT_SET:
    if not GTG.weights.has_key((OPTIONS.target,informant)):
        # not a single PacbPORF found for this informant
        continue
    ntidentity = GTG.weights[(OPTIONS.target,informant)]
    logf("GTG nt identity:", informant, ntidentity)
    if ntidentity >= MIN_QUARANTINE_HIGH_NT_IDENTITY:
        # if here, an informant with a very high nt identity
        # place this informant in quarantain
        nt_identity_quarantined_input[informant] = input[informant]
        quarantined_informant_nt_identities.append( ( ntidentity, informant ) )
        logf("# Informant Quarantined (%s); nt-identity %1.3f%s" % (
            informant,ntidentity,"%"))

################################################################################
# Check how much informants are placed in quarantine.
# If to much -> recover untill min_num_loci
################################################################################
quarantined_informant_nt_identities.sort()
quarantined_informant_nt_identities.reverse()
while nt_identity_quarantined_input and\
len(GENE_INFORMANT_SET.intersection(input.keys()).difference(nt_identity_quarantined_input.keys()))+1 < OPTIONS.minimal_num_loci:
    # place back an earlier quarantined informant informant
    ntidentity, informant = quarantined_informant_nt_identities.pop()
    del( nt_identity_quarantined_input[informant] )
    logf("# Informant UnQuarantined (%s); nt-identity %1.3f%s" % (
            informant,ntidentity,"%"))


########################################################
if VERBOSE and nt_identity_quarantined_input:
    inwpcbgs = PCG2inwpCBGS(PCG)
    logf(PCG,"prior to hign-nt-identity quarantining")
    print_inwpcbgstructure(inwpcbgs,ORGANISMS)
########################################################

# if there are to-be-quarantined informants, check the inwpCBGs.
# InwpCBGs that lack the to-be-quarantined informants can be (safely?) deleted:
# if not alignable in a >95%nt identity informant,
# it must be a non-significant hit
if nt_identity_quarantined_input:
    delete_pacbp_keys = []
    inwpcbgs = PCG2inwpCBGS(PCG)
    for inwpCBG in inwpcbgs:
        if inwpCBG.count_orfs_labeled_as_annotated_exon() >= 2: continue
        for informant in nt_identity_quarantined_input.keys():
            if informant not in inwpCBG.organism_set():
                for pacbpkey,pacbporf in inwpCBG.pacbps.iteritems():
                    if pacbporf.get_unextended_length() >= 15:
                        delete_pacbp_keys.append(pacbpkey)
                        logf("high ntident conflict to-be-deleted::", inwpCBG, pacbpkey)
                # break out if >1 high identity informant
                break

    # remove listed PacbPORF keys (and their nodes!) from the PCG
    if delete_pacbp_keys:
        for pacbpkey in delete_pacbp_keys:
            _delete_pacbp(PCG,pacbpkey)

        # update global vars after potential informant removal
        # recreate data structures
        ( ORGANISMS, GENECOMBIS,
          GENE_INFORMANT_SET,
          GENE_IDENTIFIER_SET,
          UNIGENE_INFORMANT_SET) = recreate_global_informantdata(PCG,OPTIONS.target)

        # check if an informant is not in ORGANISMS but still in input
        for informant in Set(input.keys()).difference(ORGANISMS):
            status = removeinformantfrominput(informant,input)


# logf message for high-nt-identity quarantined informants
logf("#", stw.lap(), "high-nt-identity discrepancies removed:",
    len(nt_identity_quarantined_input), nt_identity_quarantined_input.keys() )

################################################################################
### place nt_identity_quarantined_input.keys() in quarantine
################################################################################
quarantinePCG = PacbpCollectionGraph(crossdata={},blastmatrix=BLASTP_SIMILARITY_MATRIX)

for informant in nt_identity_quarantined_input.keys():
    delete_pacbp_keys = []
    for (pacbpkey,nodeQ,nodeS),pacbporf in PCG.pacbps.iteritems():
        if PCG.organism_by_node(nodeS) == informant:
            delete_pacbp_keys.append((pacbpkey,nodeQ,nodeS))
    for key in delete_pacbp_keys:
        (pacbpkey,nodeQ,nodeS) = key
        pacbporf = PCG.pacbps[key]
        # store into quarantinePCG
        quarantinePCG.pacbps[key] = pacbporf
        quarantinePCG.add_node(nodeQ)
        quarantinePCG.add_node(nodeS)
        quarantinePCG.add_edge(nodeQ,nodeS,wt=pacbporf.bitscore)
        ############################################
        if VERBOSE:
            print pacbporf
            pacbporf.print_protein_and_dna()
        ############################################
        # remove pacbporf from main PCG
        _delete_pacbp(PCG,key)
    # remove informant (and its) nodes from main PCG
    try:
        for node in PCG.get_organism_nodes(informant): PCG.del_node(node)
    except OrganismNotPresentInGraph:
        pass
    # remove from ORGANISMS and input; recreate other data structures
    del( input[informant] )
    deleted = ORGANISMS.remove(informant)
    deleted = GENECOMBIS.remove((OPTIONS.target,informant))
    deleted = GENE_INFORMANT_SET.remove(informant)
    deleted = GENE_IDENTIFIER_SET.remove(informant)

if nt_identity_quarantined_input:
    # recreate GTG after informant quarantining
    GTG = pcg2gtg_by_identity(PCG,OPTIONS.target,identifier_list=GENE_INFORMANT_SET)

# log message high nt-identity quarantine done
logf("#", stw.lap(), "high nt-identity informants quarantined:", quarantinePCG)

########################################################
if VERBOSE and nt_identity_quarantined_input:
    inwpcbgs = PCG2inwpCBGS(PCG)
    logf(PCG,"after hign-nt-identity quarantining")
    print_inwpcbgstructure(inwpcbgs,ORGANISMS)
########################################################

################################################################################
# High nt identity quarantining done. MileStone logging
################################################################################
if not VERBOSE: stdoutManager.reset()
logf("##", stwMS.lap(), PCG, "high nt-identity quanrantine performed")

################################################################################
### starting next algorithm step: InwardsPointingCodingBlockGraph filtering
################################################################################
if not VERBOSE: stdoutManager.buffer() # capture all STDOUT messages from deeper modules 
stwMS.stepname = "InwpCBGs"
stw.stepname = "InwpCBGs"


################################################################################
### first step: remove poorly covered which have NO known OrfIds
################################################################################
if len(GENE_INFORMANT_SET) > 3:
    inwpcbgs = PCG2inwpCBGS(PCG)
    del_inwpcbg_cnt, del_pacbp_cnt = remove_poorly_covered_unannotated_inwpcbgs(
        inwpcbgs,PCG,GENE_INFORMANT_SET)


    if del_pacbp_cnt > 0:
        # check if an informant is completely deleted from the PCG
        for informant in Set(ORGANISMS).difference(PCG.organism_set()):
            status = removeinformantfrominput(informant,input)

        # recreate data structures
        ( ORGANISMS, GENECOMBIS,
          GENE_INFORMANT_SET,
          GENE_IDENTIFIER_SET,
          UNIGENE_INFORMANT_SET) = recreate_global_informantdata(PCG,OPTIONS.target)

        # check if an informant is not in ORGANISMS but still in input
        for informant in Set(input.keys()).difference(ORGANISMS):
            status = removeinformantfrominput(informant,input)

        # recreate inwpCBGs
        inwpcbgs = PCG2inwpCBGS(PCG)

else:
    del_inwpcbg_cnt, del_pacbp_cnt = 0, 0


logf("#", stw.lap(), PCG, "poorly covered and unannoted removed:",
        del_inwpcbg_cnt, del_pacbp_cnt )


# print structure in all modes
inwpcbgs = PCG2inwpCBGS(PCG)
logf("#", stw.lap(), PCG, "prior to inwpCBGs filtering operations")
if VERBOSE: print_inwpcbgstructure(inwpcbgs,ORGANISMS)


################################################################################
### cleanup nodes without edges
################################################################################
removed = PCG.remove_edgeless_nodes()

if VERBOSE:
    logf(stw.lap(), PCG, "nodes without edges removed:", removed)
    for org in ORGANISMS:
        logf(org, [ orfid for (orgid,orfid) in PCG.get_organism_nodes(org) ],
            [ len(PCG.nodes[node]) for node in PCG.get_organism_nodes(org) ] )

################################################################################
### second step: add (potential, strong) SignalP exon alignemnts
################################################################################
signalpexonseqs = predict_signalp_exons(input,GENE_IDENTIFIER_SET,OPTIONS)

status = update_PCG_with_signalpexons(signalpexonseqs,PCG,OPTIONS)

if status:
    # signalp exon PacbpORFS added to PCG; recreate inwpCBGs
    inwpcbgs = PCG2inwpCBGS(PCG)

# logf message
logf("#", stw.lap(), PCG, "SignalP first exons predicted:", status)



################################################################################
### make blocks & orfidstruct
################################################################################
# translate to blocks
blocks = PCG2blocks(PCG,verbose=False)

# translate to orfidstruct
orfidstruct, blockMSRs = blocks2orfidstruct(blocks,PCG,ORGANISMS,OPTIONS)


################################################################################
### remove blocks that are certainly non-coding alignments
################################################################################

status = remove_noncoding_blocks(blocks,PCG)

# REDO translate to blocks (if a block was deleted)
if status:
    logf("remove_noncoding_blocks")
    ### DO NOT call this function: PCG is not a StarGraph yet!
    ###PCG.remove_low_connectivity_nodes(min_connectivity=1)

    # check if an informant is completely deleted from the PCG
    for informant in Set(ORGANISMS).difference(PCG.organism_set()):
        status = removeinformantfrominput(informant,input)

    # recreate data structures
    ( ORGANISMS, GENECOMBIS,
      GENE_INFORMANT_SET,
      GENE_IDENTIFIER_SET,
      UNIGENE_INFORMANT_SET) = recreate_global_informantdata(PCG,OPTIONS.target)

    # check if an informant is not in ORGANISMS but still in input
    for informant in Set(input.keys()).difference(ORGANISMS):
        status = removeinformantfrominput(informant,input)

    # REDO translation to blocks -> TODO: replace blocks for inwpCBGs
    blocks = PCG2blocks(PCG, verbose = False )

    ####################################################
    if VERBOSE:
        # recreate inwpcbgs
        inwpcbgs = PCG2inwpCBGS(PCG)
        logf(PCG,"remove_noncoding_blocks()")
        print_inwpcbgstructure(inwpcbgs,ORGANISMS)
    ####################################################

# logf message
logf("#",stw.lap(),PCG,status,"status = remove_noncoding_blocks()")

################################################################################
### cleanup blocks of which the parent PacBPORF are probably
### non-coding AND are when the parent is omitted,
### fully contained in another block
################################################################################

status = remove_noncoding_blockorigins(blocks,PCG)

# REDO translate to blocks (if a block was deleted)
if status:
    print "### remove_noncoding_blockorigins"
    ### DO NOT call this function: PCG is not a StarGraph yet!
    ###PCG.remove_low_connectivity_nodes(min_connectivity=1)

    # check if an informant is completely deleted from the PCG
    for informant in Set(ORGANISMS).difference(PCG.organism_set()):
        status = removeinformantfrominput(informant,input)

    # recreate data structures
    ( ORGANISMS, GENECOMBIS,
      GENE_INFORMANT_SET,
      GENE_IDENTIFIER_SET,
      UNIGENE_INFORMANT_SET) = recreate_global_informantdata(PCG,OPTIONS.target)

    # check if an informant is not in ORGANISMS but still in input
    for informant in Set(input.keys()).difference(ORGANISMS):
        status = removeinformantfrominput(informant,input)

    # REDO translation to blocks -> TODO: replace blocks for inwpCBGs
    blocks = PCG2blocks(PCG, verbose = False )

    ####################################################
    if VERBOSE:
        # recreate inwpcbgs
        inwpcbgs = PCG2inwpCBGS(PCG)
        logf(PCG,"remove_noncoding_blockorigins()")
        print_inwpcbgstructure(inwpcbgs,ORGANISMS)
    ####################################################

# logf message
logf("#",stw.lap(),PCG,status,"status = remove_noncoding_blockorigins()")

################################################################################
### remove_unlikely_start_blocks
################################################################################

status = remove_unlikely_start_blocks(blocks,PCG,OPTIONS.target)

# REDO translate to blocks (if a block was deleted)
if status:
    print "### remove_unlikely_start_blocks"
    ### DO NOT call this function: PCG is not a StarGraph yet!
    ###PCG.remove_low_connectivity_nodes(min_connectivity=1)

    # check if an informant is completely deleted from the PCG
    for informant in Set(ORGANISMS).difference(PCG.organism_set()):
        status = removeinformantfrominput(informant,input)

    # recreate data structures
    ( ORGANISMS, GENECOMBIS,
      GENE_INFORMANT_SET,
      GENE_IDENTIFIER_SET,
      UNIGENE_INFORMANT_SET) = recreate_global_informantdata(PCG,OPTIONS.target)

    # check if an informant is not in ORGANISMS but still in input
    for informant in Set(input.keys()).difference(ORGANISMS):
        status = removeinformantfrominput(informant,input)

    # REDO translation to blocks -> TODO: replace blocks for inwpCBGs
    blocks = PCG2blocks(PCG, verbose = False )

    ####################################################
    if VERBOSE:
        # recreate inwpcbgs
        inwpcbgs = PCG2inwpCBGS(PCG)
        logf(PCG,"remove_unlikely_start_blocks()")
        print_inwpcbgstructure(inwpcbgs,ORGANISMS)
    ####################################################

# logf message
logf("#",stw.lap(),PCG,status,"status = remove_unlikely_start_blocks()")

################################################################################
### remove_overlapping_blocks_with_conflicting_orfs
################################################################################

status = remove_overlapping_blocks_with_conflicting_orfs(blocks,PCG)

# REDO translate to blocks (if a block was deleted)
if status:
    print "### remove_overlapping_blocks_with_conflicting_orfs"
    ### DO NOT call this function: PCG is not a StarGraph yet!
    ###PCG.remove_low_connectivity_nodes(min_connectivity=1)

    # check if an informant is completely deleted from the PCG
    for informant in Set(ORGANISMS).difference(PCG.organism_set()):
        status = removeinformantfrominput(informant,input)

    # recreate data structures
    ( ORGANISMS, GENECOMBIS,
      GENE_INFORMANT_SET,
      GENE_IDENTIFIER_SET,
      UNIGENE_INFORMANT_SET) = recreate_global_informantdata(PCG,OPTIONS.target)

    # check if an informant is not in ORGANISMS but still in input
    for informant in Set(input.keys()).difference(ORGANISMS):
        status = removeinformantfrominput(informant,input)

    # REDO translation to blocks -> TODO: replace blocks for inwpCBGs
    blocks = PCG2blocks(PCG, verbose = False )

    ####################################################
    if VERBOSE:
        # recreate inwpcbgs
        inwpcbgs = PCG2inwpCBGS(PCG)
        logf(PCG,"remove_overlapping_blocks_with_conflicting_orfs()")
        print_inwpcbgstructure(inwpcbgs,ORGANISMS)
    ####################################################

# logf message
logf("#",stw.lap(),PCG,status,"status = remove_overlapping_blocks_with_conflicting_orfs()")


################################################################################
### remove_unlikely_orf_transition_blocks
################################################################################

status = remove_unlikely_orf_transition_blocks(blocks,PCG)

# REDO translate to blocks (if a block was deleted)
if status:
    print "### remove_unlikely_orf_transition_blocks"
    ### DO NOT call this function: PCG is not a StarGraph yet!
    ###PCG.remove_low_connectivity_nodes(min_connectivity=1)

    # check if an informant is completely deleted from the PCG
    for informant in Set(ORGANISMS).difference(PCG.organism_set()):
        status = removeinformantfrominput(informant,input)

    # recreate data structures
    ( ORGANISMS, GENECOMBIS,
      GENE_INFORMANT_SET,
      GENE_IDENTIFIER_SET,
      UNIGENE_INFORMANT_SET) = recreate_global_informantdata(PCG,OPTIONS.target)

    # check if an informant is not in ORGANISMS but still in input
    for informant in Set(input.keys()).difference(ORGANISMS):
        status = removeinformantfrominput(informant,input)

    # REDO translation to blocks -> TODO: replace blocks for inwpCBGs
    blocks = PCG2blocks(PCG, verbose = False )

    ####################################################
    if VERBOSE:
        # recreate inwpcbgs
        inwpcbgs = PCG2inwpCBGS(PCG)
        logf(PCG,"remove_unlikely_orf_transition_blocks()")
        print_inwpcbgstructure(inwpcbgs,ORGANISMS)
    ####################################################

# logf message
logf("#",stw.lap(),PCG,status,"status = remove_unlikely_orf_transition_blocks()")


################################################################################
### remove_singleton_blocks_without_intron_signal
################################################################################

status = remove_singleton_blocks_without_intron_signal(blocks,PCG,ORGANISMS,OPTIONS)

# REDO translate to blocks (if a block was deleted)
if status:
    print "### remove_singleton_blocks_without_intron_signal"
    ### DO NOT call this function: PCG is not a StarGraph yet!
    ###PCG.remove_low_connectivity_nodes(min_connectivity=1)

    # check if an informant is completely deleted from the PCG
    for informant in Set(ORGANISMS).difference(PCG.organism_set()):
        status = removeinformantfrominput(informant,input)
        print "JUST DELETED A;:", informant, input.keys()

    # recreate data structures
    ( ORGANISMS, GENECOMBIS,
      GENE_INFORMANT_SET,
      GENE_IDENTIFIER_SET,
      UNIGENE_INFORMANT_SET) = recreate_global_informantdata(PCG,OPTIONS.target)

    # check if an informant is not in ORGANISMS but still in input
    for informant in Set(input.keys()).difference(ORGANISMS):
        status = removeinformantfrominput(informant,input)
        print "JUST DELETED B;:", informant, input.keys()


    # REDO translation to blocks -> TODO: replace blocks for inwpCBGs
    blocks = PCG2blocks(PCG, verbose = False )

    ####################################################
    if VERBOSE:
        # recreate inwpcbgs
        inwpcbgs = PCG2inwpCBGS(PCG)
        logf(PCG,"remove_singleton_blocks_without_intron_signal()")
        print_inwpcbgstructure(inwpcbgs,ORGANISMS)
    ####################################################

# logf message
logf("#",stw.lap(),PCG,status,"status = remove_singleton_blocks_without_intron_signal()")


################################################################################
# inwpCBG filtering operations done. MileStone logging
################################################################################
if not VERBOSE: stdoutManager.reset()
logf("##", stwMS.lap(), "inwpCBG filtering steps performed", PCG)


# recreate/print inwpCBGs in ALL modes
inwpcbgs = PCG2inwpCBGS(PCG)
if VERBOSE: print_inwpcbgstructure(inwpcbgs,ORGANISMS)

################################################################################
### starting next algorithm step: tinyExon InwardsPointingCodingBlockGraph I
################################################################################
if not VERBOSE: stdoutManager.buffer() # capture all STDOUT messages from deeper modules 
stwMS.stepname = "TinyExon-I"
stw.stepname = "TinyExon-I"


################################################################################
#### Find TinyExons - simple
################################################################################

if not OPTIONS.omittinyexon and not PCG_HAS_HIGH_GAP_RATIO_PACBPORFS and\
not nt_identity_quarantined_input:
    # do tinyexonsearches; OMIT when PCG_HAS_HIGH_GAP_RATIO_PACBPORFS
    # or some informants have very high identity
    tinyexonpacbporfs = discover_tinyexons_simple(PCG,input,
            OPTIONS.target,
            GENE_INFORMANT_SET,
            UNIGENE_INFORMANT_SET)
    
    # translate tinyexonpacbporfs to tinyexonInwpCBGs
    tinyexonInwpCBGs = tinyexonpacbporfs2InwardsPointingCodingBlockGraphs(
            OPTIONS.target,
            tinyexonpacbporfs)
    
    is_any_added = False
    for tinyExonInwpCBG in tinyexonInwpCBGs:
        print "TINYEXONCBG::", tinyExonInwpCBG, tinyExonInwpCBG.get_bitscore(), tinyExonInwpCBG.get_nodes()
        # check inconsistencies between existing inwpcbgs
        is_novel_tinyexoncbg_accepted = True
        inconsistent_pacbporf_data = {}
        MINSR_tinyexon = tinyExonInwpCBG.minimal_spanning_range(organism=OPTIONS.target)
        for inwp in inwpcbgs:
            analysesA = inwp.is_positioned_compatibly(tinyExonInwpCBG)
            analysesO = inwp.is_including(tinyExonInwpCBG)
            MINSR_current = inwp.minimal_spanning_range(organism=OPTIONS.target)
            print inwp, inwp.get_bitscore(), analysesA, analysesO
            if analysesO.count(False) == 0:
                # tinyExonInwpCBG is COMPLETELY included in an
                # already existing inwpCBG
                logf("OVERLAP::", inwp.get_organism_nodes(OPTIONS.target),
                    tinyExonInwpCBG.get_organism_nodes(OPTIONS.target))
                # set is_novel_tinyexoncbg_accepted to False and break out
                is_novel_tinyexoncbg_accepted = False
                break

            if not MINSR_tinyexon.difference(MINSR_current) and\
            inwp.get_organism_nodes(OPTIONS.target) ==\
            tinyExonInwpCBG.get_organism_nodes(OPTIONS.target):
                # novel tinyexon CBG is already included as a
                # longer exon in the PCG. Omit it here
                print "OVERLAP::", inwp.get_organism_nodes(OPTIONS.target),
                print tinyExonInwpCBG.get_organism_nodes(OPTIONS.target)
                is_novel_tinyexoncbg_accepted = False
                break
    
            if not analysesA: continue
            if analysesA.count(False)==0: continue
    
            # if current existing inwp is (very) convincing: leave intact
            if inwp.node_count() > tinyExonInwpCBG.node_count():
                is_novel_tinyexoncbg_accepted = False
                break
            # if here, then there are some conflicting pacbporfs
            if analysesA.count(True)==0:
                for key,pacbporf in inwp.pacbps.iteritems():
                    inconsistent_pacbporf_data[key] = pacbporf
            else:
                for (pacbpkey,nodeQ,nodeS),pacbporf in inwp.pacbps.iteritems():
                    dummy = InwardsPointingCodingBlockGraph()
                    dummy.add_node(nodeQ)
                    dummy.add_node(nodeS)
                    dummy.add_edge(nodeQ,nodeS,wt=pacbporf.bitscore)
                    dummy.pacbps[(pacbpkey,nodeQ,nodeS)] = pacbporf
                    if dummy.is_positioned_compatibly(tinyExonInwpCBG) == [False]:
                        inconsistent_pacbporf_data[(pacbpkey,nodeQ,nodeS)] = pacbporf
    
    
        # check if we can add this tinyexon inwpCBG to the main PCG
        if is_novel_tinyexoncbg_accepted:
            print "ADD TO PCG:", tinyExonInwpCBG, tinyExonInwpCBG.get_nodes()
            print "inconsistent:", [ c for a,b,c in inconsistent_pacbporf_data ]
            print PCG
            is_any_added = True
            for key in inconsistent_pacbporf_data.keys(): _delete_pacbp(PCG,key)
            del( inconsistent_pacbporf_data )
            for (pacbpkey,nodeQ,nodeS),pacbporf in tinyExonInwpCBG.pacbps.iteritems():
                if not nodeQ in PCG.get_nodes(): PCG.add_node(nodeQ)
                if not nodeS in PCG.get_nodes(): PCG.add_node(nodeS)
                PCG.add_edge(nodeQ,nodeS,wt=pacbporf.bitscore)
                # label this tinyExonPacbpORF as non-deletable!
                pacbporf._IS_DELETE_PROTECTED = True
                # and store to the PCG
                PCG.pacbps[(pacbpkey,nodeQ,nodeS)] = pacbporf
            print PCG
    
    
    if VERBOSE and is_any_added:
        inwpcbgs = PCG2inwpCBGS(PCG)
        print_inwpcbgstructure(inwpcbgs,ORGANISMS)


################################################################################
# Find TinyExons - simple. MileStone logging
################################################################################
if not VERBOSE: stdoutManager.reset()
logf("##", stwMS.lap(), "Find TinyExons - simple: performed", PCG)

################################################################################
### starting next algorithm step: ORFId transition analyses
################################################################################
if not VERBOSE: stdoutManager.buffer() # capture all STDOUT messages from deeper modules 
stwMS.stepname = "ORFid analyses"
stw.stepname = "ORFid analyses"


################################################################################
### Find gaps in the Orfid structure between informant and target,
### which are convincingly filled by other informants.
### Remove these potentially colinearity-breaking PacbPORFs.
### This must be done prior to HMM analyses!
################################################################################
pacbporf_removal_cnt = 0
while True:
    pacbporf_removal_list = []
    orfid_list_target = [ inwpCBG._get_target_node()[1] for inwpCBG in inwpcbgs ]
    logf("TARGET orfids:", orfid_list_target)
    for informant in GENE_INFORMANT_SET:
        orfid_list_informant = []
        for inwpCBG in inwpcbgs:
            if informant in inwpCBG.organism_set():
                orfid_list_informant.append( inwpCBG.get_organism_nodes(informant)[0][1] )
            else:
                orfid_list_informant.append(None)
        logf("INFORMANT orfids:",orfid_list_informant, informant)

        # find gaps in the structure for this informant which are
        # filled for other informants
        offset = 0
        while True:
            try:
                offset = orfid_list_informant.index(None,offset)
            except ValueError:
                break
            prevorfpos, nextorfpos = None, None
            for pos in range(offset-1,-1,-1):
                if orfid_list_informant[pos] != None:
                    prevorfpos = int(pos)
                    break
            if prevorfpos == None:
                offset+=1
                continue
            for pos in range(offset+1,len(orfid_list_informant),1):
                if orfid_list_informant[pos] != None:
                    nextorfpos = int(pos)
                    break
            if nextorfpos == None:
                offset+=1
                continue

            # get pacbporfs surrounding the gap
            prevpacbporf = inwpcbgs[prevorfpos].get_pacbps_by_organisms(OPTIONS.target,informant)[0]
            nextpacbporf = inwpcbgs[nextorfpos].get_pacbps_by_organisms(OPTIONS.target,informant)[0]

            # continue if prev and next pacbporf are the same (can happen; in fact inwpcbg assembly error!)
            if prevpacbporf.barcode() == nextpacbporf.barcode():
                # continue with the next iteration
                offset = nextorfpos
                # continue with the next iteration
                continue

            # obtain some data on size of the gap and reliability of its filling
            intermediate_inwpcbg_gap = min(inwpcbgs[nextorfpos].minimal_spanning_range(organism=OPTIONS.target)) -\
                                       max(inwpcbgs[prevorfpos].minimal_spanning_range(organism=OPTIONS.target))
            intermediate_inwpcbg_score = sum([ inwpCBG.get_bitscore() for inwpCBG in inwpcbgs[prevorfpos+1:nextorfpos] ])
            intermediate_query_gap = min(nextpacbporf.alignment_protein_range_query()) -\
                                     max(prevpacbporf.alignment_protein_range_query())
            intermediate_sbjct_gap = min(nextpacbporf.alignment_protein_range_sbjct()) -\
                                     max(prevpacbporf.alignment_protein_range_sbjct())
            prev_inwpcbg_score = inwpcbgs[prevorfpos].get_bitscore()
            next_inwpcbg_score = inwpcbgs[nextorfpos].get_bitscore()

            # obtain data on the size of MINSR/MAXSR support in this gap
            intermediate_inwpcbg_maxsr = Set()
            intermediate_inwpcbg_minsr = Set()
            for inwpCBG in inwpcbgs[prevorfpos+1:nextorfpos]:
                intermediate_inwpcbg_maxsr.update(inwpCBG.maximal_spanning_range(organism=OPTIONS.target))
                intermediate_inwpcbg_minsr.update(inwpCBG.minimal_spanning_range(organism=OPTIONS.target))
            intermediate_inwpcbg_maxsr.difference_update(inwpcbgs[nextorfpos].minimal_spanning_range(organism=OPTIONS.target))
            intermediate_inwpcbg_maxsr.difference_update(inwpcbgs[prevorfpos].minimal_spanning_range(organism=OPTIONS.target))
            intermediate_inwpcbg_minsr_len = len(intermediate_inwpcbg_minsr)
            intermediate_inwpcbg_maxsr_len = len(intermediate_inwpcbg_maxsr)
            
            # print info
            logf(OPTIONS.target,"\t",orfid_list_target[prevorfpos:nextorfpos+1], "prev:", prevpacbporf)
            logf(informant, "\t", orfid_list_informant[prevorfpos:nextorfpos+1], "next:", nextpacbporf)
            logf("gapsize:", intermediate_inwpcbg_gap,
                    "Q:", intermediate_query_gap,
                    "S:", intermediate_sbjct_gap,
                    "score:",intermediate_inwpcbg_score,
                    "minsr:",len(intermediate_inwpcbg_minsr),
                    "maxsr:",len(intermediate_inwpcbg_maxsr))

            if intermediate_sbjct_gap < intermediate_inwpcbg_minsr_len or\
            intermediate_query_gap < intermediate_inwpcbg_minsr_len and\
            intermediate_inwpcbg_score > min([prev_inwpcbg_score,next_inwpcbg_score]):
                # remove the least convincing of prev and next pacbporf
                if prev_inwpcbg_score < next_inwpcbg_score:
                    nodeQ = (OPTIONS.target, prevpacbporf.orfQ.id)
                    nodeS = (informant, prevpacbporf.orfS.id)
                    pacbporf_removal_list.append( ( prevpacbporf, nodeQ, nodeS ) )
                    logf("deleting:", informant, prevpacbporf )
                else:
                    nodeQ = (OPTIONS.target, nextpacbporf.orfQ.id)
                    nodeS = (informant, nextpacbporf.orfS.id)
                    pacbporf_removal_list.append( ( nextpacbporf, nodeQ, nodeS ) )
                    logf("deleting:", informant, nextpacbporf )

            elif intermediate_sbjct_gap < 0.9*intermediate_inwpcbg_maxsr_len and\
            intermediate_query_gap < 0.9*intermediate_inwpcbg_maxsr_len and\
            intermediate_inwpcbg_score > min([prev_inwpcbg_score,next_inwpcbg_score]):
                # remove the least convincing of prev and next pacbporf
                if prev_inwpcbg_score < next_inwpcbg_score:
                    nodeQ = (OPTIONS.target, prevpacbporf.orfQ.id)
                    nodeS = (informant, prevpacbporf.orfS.id)
                    pacbporf_removal_list.append( ( prevpacbporf, nodeQ, nodeS ) )
                    logf("deleting:", informant, prevpacbporf )
                else:
                    nodeQ = (OPTIONS.target, nextpacbporf.orfQ.id)
                    nodeS = (informant, nextpacbporf.orfS.id)
                    pacbporf_removal_list.append( ( nextpacbporf, nodeQ, nodeS ) )
                    logf("deleting:", informant, nextpacbporf )
           
            # continue with the next iteration
            offset = nextorfpos

    # if no pacbporfs are labeled for deletion -> break while loop
    if not pacbporf_removal_list: break

    # increase counter
    pacbporf_removal_cnt += len(pacbporf_removal_list)

    # remove pacbporfs from main PCG
    deletion_failed = False
    for pfObj,nodeQ,nodeS in pacbporf_removal_list:
        delete_key = pfObj.construct_unique_key(nodeQ,nodeS)
        for (key,n1,n2),pacbporfObj in PCG.pacbps.iteritems():
            if (key,n1,n2) == (delete_key,nodeQ,nodeS):
                _delete_pacbp(PCG,(delete_key,nodeQ,nodeS))
                break
            elif (n1,n2) == (nodeQ,nodeS):
                _key = pacbporfObj.construct_unique_key(nodeQ,nodeS)
                if _key == delete_key:
                    # delete by old key -> pacbporf has been edited somewhere
                    _delete_pacbp(PCG,(key,nodeQ,nodeS))
                    break
        else:
            # deletion failed! This is severe, because it can/will
            # result in an ethernal while loop.
            # Prevent this by hard-existing somewhat later
            deletion_failed = True
        logf("#", "DELETED::",nodeQ,nodeS,pfObj)
        
    # recreate inwpCBGs
    inwpcbgs = PCG2inwpCBGS(PCG)

    # check if the deletion was unsuccesfull -> if so, break out!
    if deletion_failed: break
    
# logf message
logf("#",stw.lap(),PCG,"pacbporfs removed by OrfID & gap analyses ",pacbporf_removal_cnt)


################################################################################
### Find gaps in the Orfid structure between informant and target,
### and fill these gaps by ClustalW obtained PabPORFS. This step saves
### tremendous amounts of time lateron in the algorithm (HMM, intronmapping)
################################################################################
pacbporf_replacement_cnt = 0
orfid_list_target = [ inwpCBG._get_target_node()[1] for inwpCBG in inwpcbgs ]
logf("TARGET orfids:", orfid_list_target)
for informant in GENE_INFORMANT_SET:
    orfid_list_informant = []
    for inwpCBG in inwpcbgs:
        if informant in inwpCBG.organism_set():
            orfid_list_informant.append( inwpCBG.get_organism_nodes(informant)[0][1] )
        else:
            orfid_list_informant.append(None)
    logf("INFORMANT orfids:",orfid_list_informant, informant)

    offset = 0
    while True:
        try:
            offset = orfid_list_informant.index(None,offset)
        except ValueError:
            break
        prevorfpos, nextorfpos = None, None
        for pos in range(offset-1,-1,-1):
            if orfid_list_informant[pos] != None:
                prevorfpos = int(pos)
                break
        if prevorfpos == None:
            offset+=1
            continue
        for pos in range(offset+1,len(orfid_list_informant),1):
            if orfid_list_informant[pos] != None:
                nextorfpos = int(pos)
                break
        if nextorfpos == None:
            offset+=1
            continue

        if orfid_list_informant[prevorfpos] != orfid_list_informant[nextorfpos]:
            offset = nextorfpos
            continue

        # print info
        logf(OPTIONS.target,"\t",orfid_list_target[prevorfpos:nextorfpos+1])
        logf(informant, "\t", orfid_list_informant[prevorfpos:nextorfpos+1])

        is_equal_final_orf = orfid_list_target[nextorfpos] == orfid_list_target[nextorfpos-1]
        is_equal_first_orf = orfid_list_target[prevorfpos] == orfid_list_target[prevorfpos+1]
        prevpacbporf = inwpcbgs[prevorfpos].get_pacbps_by_organisms(OPTIONS.target,informant)[0]
        nextpacbporf = inwpcbgs[nextorfpos].get_pacbps_by_organisms(OPTIONS.target,informant)[0]
        prevpacbp = pacb.conversion.pacbporf2pacbp(prevpacbporf)
        nextpacbp = pacb.conversion.pacbporf2pacbp(nextpacbporf)
        logf("==first",is_equal_first_orf,"==final",is_equal_final_orf)
        logf("prev:", prevpacbporf)
        logf("next:", nextpacbporf)
        
        if is_equal_first_orf and is_equal_final_orf:
            orfQobj = input[OPTIONS.target]['orfs'].get_orf_by_id( orfid_list_target[prevorfpos] )
            orfSobj = input[informant]['orfs'].get_orf_by_id( orfid_list_informant[prevorfpos] )
            replacements = [ prevpacbporf, nextpacbporf ]
        elif is_equal_first_orf:
            orfQobj = input[OPTIONS.target]['orfs'].get_orf_by_id( orfid_list_target[prevorfpos] )
            orfSobj = input[informant]['orfs'].get_orf_by_id( orfid_list_informant[prevorfpos] )
            replacements = [ prevpacbporf ]
        elif is_equal_final_orf:
            orfQobj = input[OPTIONS.target]['orfs'].get_orf_by_id( orfid_list_target[nextorfpos] )
            orfSobj = input[informant]['orfs'].get_orf_by_id( orfid_list_informant[nextorfpos] )
            replacements = [ nextpacbporf ]
        else:
            # none are equal!? possible, but unusual Orf pattern
            offset = nextorfpos
            continue

        # start building realignable fraction of Orf sequences
        orfQseq = orfQobj.protein_sequence.upper()
        orfSseq = orfSobj.protein_sequence.upper()
        
        if is_equal_first_orf and is_equal_final_orf:
            orfQseq = orfQobj.getaas(abs_pos_start=prevpacbp.query_start,
                                     abs_pos_end=nextpacbp.query_end)
            orfSseq = orfSobj.getaas(abs_pos_start=prevpacbp.sbjct_start,
                                     abs_pos_end=nextpacbp.sbjct_end)
            coords = [ prevpacbp.query_start, nextpacbp.query_end,
                       prevpacbp.sbjct_start, nextpacbp.sbjct_end]
        elif is_equal_first_orf:
            orfQseq = orfQobj.getaas(abs_pos_start=prevpacbp.query_start,
                                     abs_pos_end=orfQobj.protein_endPY)
            orfSseq = orfSobj.getaas(abs_pos_start=prevpacbp.sbjct_start,
                                     abs_pos_end=nextpacbp.sbjct_start)

            coords = [ prevpacbp.query_start, orfQobj.protein_endPY,
                       prevpacbp.sbjct_start, nextpacbp.sbjct_start]
        elif is_equal_final_orf:
            orfQseq = orfQobj.getaas(abs_pos_start=orfQobj.protein_startPY,
                                     abs_pos_end=nextpacbp.query_end)
            orfSseq = orfSobj.getaas(abs_pos_start=prevpacbp.sbjct_end,
                                     abs_pos_end=nextpacbp.sbjct_end)
            coords = [ orfQobj.protein_startPY, nextpacbp.query_end,
                       prevpacbp.sbjct_end, nextpacbp.sbjct_end]
        else:
            # not possible.
            pass

        # print some debugging info
        logf(OPTIONS.target, orfQobj, len(orfQseq), informant, orfSobj, len(orfSseq), coords)
        
        clwinput = {OPTIONS.target: orfQseq, informant: orfSseq }
        (alignedseqs,alignment) = clustalw( seqs=clwinput)

        # make pacbp from clustalw alignment
        newpacbp = pacbp_from_clustalw(
                    alignment=(
                            alignedseqs[OPTIONS.target],
                            alignment,
                            alignedseqs[informant]
                            ),
                    coords=coords
                    )

        # check if creation failed -> cannot be the case here, but safety first
        if not newpacbp: continue
        # strip unaligned fraction of this pacbp object, then check length
        newpacbp.strip_unmatched_ends()

        # split on (very) grove gaps
        (splittedpacbps, is_splitted) = pacb.splitting.split_pacb_on_gaps(newpacbp,gapsize=18)
        newpacbporfs = []
        for splittedpacbpObj in splittedpacbps:
            newpacbporf = pacbp2pacbporf(splittedpacbpObj,orfQobj,orfSobj)
            newpacbporf.extend_pacbporf_after_stops()
            newpacbporfs.append(newpacbporf)

        ########################################################################
        if VERBOSE:
            print prevpacbporf
            prevpacbp.print_protein(_linesize=150) 
            print nextpacbporf
            nextpacbp.print_protein(_linesize=150)
            print coords, (len(orfQseq), coords[1]-coords[0]),
            print (len(orfSseq), coords[3]-coords[2])
            for part in range(0,len(alignment),150):
                print alignedseqs[OPTIONS.target][part:part+150]
                print alignment[part:part+150]
                print alignedseqs[informant][part:part+150]
            for newpacbp in splittedpacbps:
                print "NEW::", newpacbp
                newpacbp.print_protein(_linesize=150)
        ########################################################################


        # now replace the old pacbporf(s) with the new pacbporf 
        nodeQ = (OPTIONS.target,orfQobj.id)
        nodeS = (informant,orfSobj.id)
        for pf in replacements:
            delete_key = None
            for (key,nQ,nS),pfObj in PCG.pacbps.iteritems():
                if nQ != nodeQ: continue
                if nS != nodeS: continue
                if pf.construct_unique_key(nodeQ,nodeS) ==\
                pfObj.construct_unique_key(nodeQ,nodeS):
                    delete_key = (key,nQ,nS)
                    break
            if delete_key:
                # delete this PacbPORF key from the PCG
                del( PCG.pacbps[delete_key] )
                # increase replacement counter
                pacbporf_replacement_cnt+=1

        # and store the new one(s)
        for newpacbporf in newpacbporfs:
            new_key = newpacbporf.construct_unique_key(nodeQ,nodeS)
            PCG.pacbps[(new_key,nodeQ,nodeS)] = newpacbporf

        # continue with the next iteration
        offset = nextorfpos

# logf message
logf("#",stw.lap(),PCG,"OrfID gaps discovered & filled:",pacbporf_replacement_cnt)


################################################################################
# ORFid analyses done. MileStone logging
################################################################################
if not VERBOSE: stdoutManager.reset()
logf("##", stwMS.lap(), "ORFid analyses: performed", PCG)

if VERBOSE and pacbporf_replacement_cnt:
    # recreate inwpcbgs
    inwpcbgs = PCG2inwpCBGS(PCG)
    print_inwpcbgstructure(inwpcbgs,ORGANISMS)
   

################################################################################
### starting next algorithm step: main HMM searches
################################################################################
if not VERBOSE: stdoutManager.buffer() # capture all STDOUT messages from deeper modules 
stwMS.stepname = "HMM"
stw.stepname = "HMM"

################################################################################
#### Start preparation for elaborate HMM searches of inwpCBGs
################################################################################

logf("#",stw.lap(),"Starting HMM searches")

if not OPTIONS.omithmm:
    # recreate inwpcbgs
    inwpcbgs = PCG2inwpCBGS(PCG)
    if VERBOSE: print_inwpcbgstructure(inwpcbgs,ORGANISMS)

HMM_MAXSR_MIN_AA_LENGTH                 = 20
HMM_MAXSR_HMM_MIN_BITSCORE              = -20
HMM_MAXSR_HMM_SPLIT_ON_GAPSIZE          = 4
HMM_MAXSR_PACBPORF_MIN_BITSCORE         = 10
HMM_MAXSR_PACBPORF_MIN_AA_LENGTH        = 15

HMM_LEFTSPRDIF_MIN_AA_LENGTH            = 20
HMM_LEFTSPRDIF_HMM_MIN_BITSCORE         = 0
HMM_LEFTSPRDIF_HMM_SPLIT_ON_GAPSIZE     = 100 # no split!
HMM_LEFTSPRDIF_PACBPORF_MIN_BITSCORE    = 10
HMM_LEFTSPRDIF_PACBPORF_MIN_AA_LENGTH   = 15

HMM_RIGTHSPRDIF_MIN_AA_LENGTH           = 20
HMM_RIGTHSPRDIF_HMM_MIN_BITSCORE        = 0
HMM_RIGTHSPRDIF_HMM_SPLIT_ON_GAPSIZE    = 100 # no split!
HMM_RIGTHSPRDIF_PACBPORF_MIN_BITSCORE   = 10
HMM_RIGTHSPRDIF_PACBPORF_MIN_AA_LENGTH  = 15


PACBPORF_HIGH_GAP_RATIO_THRESHOLD = 0.40

HMM_MAXSR_UNCORRECTED_MIN_AA_LENGTH          = 50 # n.a.
HMM_MAXSR_UNCORRECTED_HMM_MIN_BITSCORE       = -1000
HMM_MAXSR_UNCORRECTED_HMM_SPLIT_ON_GAPSIZE   = 20    # practically no splitting
HMM_MAXSR_UNCORRECTED_PACBPORF_MIN_BITSCORE  = -1000 # practically no filtering
HMM_MAXSR_UNCORRECTED_PACBPORF_MIN_AA_LENGTH = 100
HMM_MAXSR_UNCORRECTED_PACBPORF_MIN_GAP_RATIO_SCORE = PACBPORF_HIGH_GAP_RATIO_THRESHOLD 

#################################################################################
#### HMM MAXSR UNCORRECTED search for large poor indentity Orfs in inwpCBGs
#################################################################################
logf(stw.lap(),"HMM MAXSR UNCORRECTED searches start")

# TODO TODO protect this function by checking for long PacbPORF alignments
# having very *shitty* alignments with many gaps.
# TODO TODO this will prevent overprediction of False alignments

final_inwpcbg_is_reached = False
iteration_count = 0
while not final_inwpcbg_is_reached and PCG_HAS_HIGH_GAP_RATIO_PACBPORFS:
    # only allow this type of HMM search when PCG_HAS_HIGH_GAP_RATIO_PACBPORFS
    for pos in range(0,len(inwpcbgs)):
        # check for OPTIONS.omithmm -> no HMM profile searches!
        if OPTIONS.omithmm: continue
        inwpCBG = inwpcbgs[pos]

        if not GENE_IDENTIFIER_SET.difference(inwpCBG.organism_set()):
            # no organisms/genes/nodes missing -> continue
            continue

        # Check pacbporf.gap_ratio_score() in all PacbPORFs. In case of low
        # scores, this is a signal to break out of this function. Here, we
        # are looking for long, poor-idenity alignments of (in some species)
        # putatively single-exon genes with much coding sequence variation
        if max([ pacbporf.gap_ratio_score() for pacbporf in inwpCBG.pacbps.values() ]) <\
        HMM_MAXSR_UNCORRECTED_PACBPORF_MIN_GAP_RATIO_SCORE:
            # break out of this function
            continue

        # log message in VERBOSE mode
        logf("iteration_count:", iteration_count, inwpCBG)
    
        # get previous & next inwpCBG -> None in case of UNCORRECTED search!
        prev,next = None,None
    

        # create MAXSR HMM search profile
        fname_hmmbuild, hmmcoords = _create_hmm_profile(
                    inwpCBG,area="MAXSR",
                    prevcbg=prev,nextcbg=next,verbose=VERBOSE)
        # no MAXSR profile can be constructed
        if not fname_hmmbuild: continue
        if min([ end-sta for sta,end in hmmcoords.values()]) < HMM_MAXSR_UNCORRECTED_PACBPORF_MIN_AA_LENGTH:
            logf("NO HMM_MAXSR_UNCORRECTED_PACBPORF_MIN_AA_LENGTH:", inwpCBG)
            _file_cleanup([fname_hmmbuild])
            continue
        # check if OPTIONS.target still in hmmcoords; in exceptional cases,
        # it can have been removed from the HMMprofile
        if not hmmcoords.has_key(OPTIONS.target):
            _file_cleanup([fname_hmmbuild])
            continue
    

        # Variable which checks if any pacbporf is added.
        # If so, this is a trigger to break out of the for-loop over
        # the inwpcbgs, to recreate them and start all over
        is_any_pacbporf_added = False
    
        for informant in ORGANISMS:
            # continue if target or already seen informant species
            if informant == OPTIONS.target: continue
            if informant in inwpCBG.organism_set(): continue
            # Do NOT perform this on unigene informants.
            # It will likley result in arbitrary elongation of alignments
            if informant in UNIGENE_INFORMANT_SET: continue


            # create fasta search database with ORFs from this informant
            fname_fasta_db = _create_hmm_db(informant,input,inwpCBG,prev,next)
            # check for no elegiable ORFs at all
            if not fname_fasta_db: continue
    
            # run hmmsearch on this fasta database with hmmbuild file; maximum 2 hits
            results = hmmsearch_protein( fname_hmmbuild, fname_fasta_db, params= {'E':1,'A': 2})
    
            # logf preliminary HMM hits in VERBOSE mode
            if results: logf([ (res[0],res[-2]) for res in results ])
    
            # create pacbporfhmmlist list
            pacbporfhmmlist = hmmresults2splittedpacbps(results,hmmcoords,
                    OPTIONS.target,informant,inwpCBG,input,
                    gapsize=HMM_MAXSR_UNCORRECTED_HMM_SPLIT_ON_GAPSIZE,
                    min_bitscore=HMM_MAXSR_UNCORRECTED_HMM_MIN_BITSCORE)
    
            # check/add the pacbporfs to the PCG
            thepacbporfs = pacb.ordering.order_pacbporf_list(PCG.get_pacbps_by_organisms(OPTIONS.target,informant))
    
            # file cleanup of HMM fasta database
            _file_cleanup([fname_fasta_db])
    
            for hmmpacbporf in pacbporfhmmlist:
                # filter out to short ones
                if hmmpacbporf.get_unextended_length() < HMM_MAXSR_UNCORRECTED_PACBPORF_MIN_AA_LENGTH:
                    continue
    
                # filter out poor bitscores
                if hmmpacbporf.bitscore < HMM_MAXSR_UNCORRECTED_PACBPORF_MIN_BITSCORE:
                    continue
    
                rejected = is_hmmpacbporf_conflicting_with_pacbporflist(
                            hmmpacbporf, thepacbporfs)
    
                if not rejected:
                    # if here, add Node & PacbPORF to the PCG.
                    hmmpacbporf2PCG(hmmpacbporf,OPTIONS.target,informant,PCG)   
                    # set is_any_pacbporf_added variable to True
                    is_any_pacbporf_added = True
                    logf("# HMM pacbporf MAXSR UNCORRECTED ADDED::", iteration_count, informant, hmmpacbporf)
                
                    # re-get list of pacbporfs
                    thepacbporfs = pacb.ordering.order_pacbporf_list(
                            PCG.get_pacbps_by_organisms(OPTIONS.target,informant)
                            )
    
            # log message in VERBOSE mode        
            logf(inwpCBG, "MAXSR UNCORRECTED", informant, len(results), len(pacbporfhmmlist))
    
        # file cleanup of HMM profile
        _file_cleanup([fname_hmmbuild])

        if is_any_pacbporf_added:
            # **IMPORTANT** redo PCG.merge_high_gap_ratio_pacbporfs()
            PCG.merge_high_gap_ratio_pacbporfs(OPTIONS.target,verbose=True)
            # **IMPORTANT** recreate inwpcbgs and restart iteration
            inwpcbgs = PCG2inwpCBGS(PCG)
            if VERBOSE: print_inwpcbgstructure(inwpcbgs,ORGANISMS)
            logf(inwpCBG, "MAXSR UNCORRECTED", iteration_count, "ENDED")
            # increase iteration count for logging purposes
            iteration_count+=1
            # break out of for loop
            break

    else:
        # end of loop over inwpcbgs is reached -> break out of while loop!
        final_inwpcbg_is_reached = True

# HMM log message
logf("#",stw.lap(),"HMM MAXSR UNCORRECTED searches done")


#################################################################################
#### HMM MAXSR search for missing Orfs in inwpCBGs
#################################################################################
logf(stw.lap(),"HMM MAXSR searches start")

if not OPTIONS.omithmm:
    # recreate inwpcbgs
    inwpcbgs = PCG2inwpCBGS(PCG)
    if VERBOSE: print_inwpcbgstructure(inwpcbgs,ORGANISMS)

for pos in range(0,len(inwpcbgs)):
    # check for OPTIONS.omithmm -> no HMM profile searches!
    if OPTIONS.omithmm: continue
    inwpCBG = inwpcbgs[pos]

    if not GENE_IDENTIFIER_SET.difference(inwpCBG.organism_set()):
        # no organisms/genes/nodes missing -> continue
        continue

    # get previous & next inwpCBG
    prev,next = None,None
    if pos > 0: prev = inwpcbgs[pos-1]
    if pos < len(inwpcbgs)-1: next = inwpcbgs[pos+1]

    # log message in VERBOSE mode
    logf(inwpCBG,inwpCBG._get_target_node())
    # create MAXSR HMM search profile
    fname_hmmbuild, hmmcoords = _create_hmm_profile(
                inwpCBG,area="MAXSR",
                prevcbg=prev,nextcbg=next,verbose=VERBOSE)
    # no MAXSR profile can be constructed
    if not fname_hmmbuild: continue
    if min([ end-sta for sta,end in hmmcoords.values()]) < HMM_MAXSR_PACBPORF_MIN_AA_LENGTH:
        logf("NO HMM_MAXSR_PACBPORF_MIN_AA_LENGTH:", inwpCBG)
        _file_cleanup([fname_hmmbuild])
        continue
    # check if OPTIONS.target still in hmmcoords; in exceptional cases,
    # it can have been removed from the HMMprofile
    if not hmmcoords.has_key(OPTIONS.target):
        _file_cleanup([fname_hmmbuild])
        continue


    for informant in ORGANISMS:
        # continue if target or already seen informant species
        if informant == OPTIONS.target: continue
        if informant in inwpCBG.organism_set(): continue

        # create fasta search database with ORFs from this informant
        fname_fasta_db = _create_hmm_db(informant,input,inwpCBG,prev,next)
        # check for no elegiable ORFs at all
        if not fname_fasta_db: continue

        # run hmmsearch on this fasta database with hmmbuild file; maximum 2 hits
        results = hmmsearch_protein( fname_hmmbuild, fname_fasta_db, params= {'E':150,'A': 2})

        # logf preliminary HMM hits in VERBOSE mode
        if results: logf([ (res[0],res[-2]) for res in results ])

        # create pacbporfhmmlist list
        pacbporfhmmlist = hmmresults2splittedpacbps(results,hmmcoords,
                OPTIONS.target,informant,inwpCBG,input,
                gapsize=HMM_MAXSR_HMM_SPLIT_ON_GAPSIZE,
                min_bitscore=HMM_MAXSR_HMM_MIN_BITSCORE)

        # check/add the pacbporfs to the PCG
        thepacbporfs = pacb.ordering.order_pacbporf_list(PCG.get_pacbps_by_organisms(OPTIONS.target,informant))

        # file cleanup of HMM fasta database
        _file_cleanup([fname_fasta_db])

        for hmmpacbporf in pacbporfhmmlist:
            # filter out to short ones
            if hmmpacbporf.get_unextended_length() < HMM_MAXSR_PACBPORF_MIN_AA_LENGTH:
                continue

            # filter out poor bitscores
            if hmmpacbporf.bitscore < HMM_MAXSR_PACBPORF_MIN_BITSCORE:
                continue


            rejected = is_hmmpacbporf_conflicting_with_pacbporflist(
                        hmmpacbporf, thepacbporfs)

            if not rejected:
                # if here, add Node & PacbPORF to the PCG.
                hmmpacbporf2PCG(hmmpacbporf,OPTIONS.target,informant,PCG)   
                logf("# HMM pacbporf MAXSR ADDED::", informant, hmmpacbporf)
            
                # re-get list of pacbporfs
                thepacbporfs = pacb.ordering.order_pacbporf_list(
                        PCG.get_pacbps_by_organisms(OPTIONS.target,informant)
                        )

        # log message in VERBOSE mode        
        logf(inwpCBG, "MAXSR", informant, len(results), len(pacbporfhmmlist))

    # file cleanup of HMM profile
    _file_cleanup([fname_hmmbuild])

# HMM log message
logf("#",stw.lap(),"HMM MAXSR searches done")

if PCG_HAS_HIGH_GAP_RATIO_PACBPORFS:
    # perform again PCG.merge_high_gap_ratio_pacbporfs() after HMM searches
    PCG_HAS_HIGH_GAP_RATIO_PACBPORFS = PCG.merge_high_gap_ratio_pacbporfs(OPTIONS.target,verbose=True)
    logf(stw.lap(),PCG,"PCG.merge_high_gap_ratio_pacbporfs()",PCG_HAS_HIGH_GAP_RATIO_PACBPORFS)

#################################################################################
#### HMM LEFTSPRDIF search for missing Orfs in inwpCBGs
#################################################################################
logf(stw.lap(),"HMM LEFTSPRDIF searches start")

if not OPTIONS.omithmm:
    # recreate inwpcbgs
    inwpcbgs = PCG2inwpCBGS(PCG)
    if VERBOSE: print_inwpcbgstructure(inwpcbgs,ORGANISMS)

sprdif_hmmpacbporf_list = []
hmm_coord_check_list = []

for pos in range(0,len(inwpcbgs)):
    # check for OPTIONS.omithmm -> no HMM profile searches!
    if OPTIONS.omithmm: continue
    inwpCBG = inwpcbgs[pos]

    # check for left sprdif in inwpCBG
    if not inwpCBG.has_left_spanningrange_difference(
    sprdif_min_aa_length=HMM_LEFTSPRDIF_MIN_AA_LENGTH):
        continue

    # get previous & next inwpCBG
    prev,next = None,None
    if pos > 0: prev = inwpcbgs[pos-1]
    if pos < len(inwpcbgs)-1: next = inwpcbgs[pos+1]

    # create LEFTSPRDIF HMM search profile
    fname_hmmbuild, hmmcoords = _create_hmm_profile(
                inwpCBG,area="LEFTSPRDIF",
                sprdif_min_aa_length=HMM_LEFTSPRDIF_MIN_AA_LENGTH,
                prevcbg=prev,nextcbg=next,verbose=VERBOSE)

    # no LEFTSPRDIF profile can be constructed
    if not fname_hmmbuild: continue
    if min([ end-sta for sta,end in hmmcoords.values()]) < HMM_LEFTSPRDIF_PACBPORF_MIN_AA_LENGTH:
        _file_cleanup([fname_hmmbuild])
        continue

    # check if OPTIONS.target still in hmmcoords; in exceptional cases,
    # it can have been removed from the HMMprofile
    if not hmmcoords.has_key(OPTIONS.target):
        _file_cleanup([fname_hmmbuild])
        continue

    # check if this HMM profile is already evaluated
    ordered_coords = [ ( k,sta,end ) for k,(sta,end) in hmmcoords.iteritems() ]
    ordered_coords.sort()
    if ordered_coords in hmm_coord_check_list:
        _file_cleanup([fname_hmmbuild])
        continue

    # append hmmcoords to hmm_coord_check_list
    hmm_coord_check_list.append( ordered_coords )

    for informant in ORGANISMS:
        # continue if target
        if informant == OPTIONS.target: continue

        # check if, for this informant, a potential extention exists
        if hmmcoords.has_key(informant):
            thepacbporf = inwpCBG.get_pacbps_by_organisms(OPTIONS.target,informant)[0]
            hmmrange = Set(range(hmmcoords[informant][0],hmmcoords[informant][1]))
            diff = hmmrange.difference( thepacbporf.alignment_protein_range_sbjct())
            # no large enough AA difference to find HMM PacbPORF
            if len(diff) <= 10: continue


        # create fasta search database with ORFs from this informant
        fname_fasta_db = _create_hmm_db(informant,input,inwpCBG,prev,next)
        # check for no elegiable ORFs at all
        if not fname_fasta_db: continue

        # run hmmsearch on this fasta database with hmmbuild file; maximum 2 hits
        results = hmmsearch_protein( fname_hmmbuild, fname_fasta_db, params= {'E':150,'A': 2} )
        # create pacbporfhmmlist list
        pacbporfhmmlist = hmmresults2splittedpacbps(results,hmmcoords,
                OPTIONS.target,informant,inwpCBG,input,
                gapsize=HMM_LEFTSPRDIF_HMM_SPLIT_ON_GAPSIZE,
                min_bitscore=HMM_LEFTSPRDIF_HMM_MIN_BITSCORE)

        # file cleanup of HMM fasta database
        _file_cleanup([fname_fasta_db])

        for hmmpacbporf in pacbporfhmmlist:
            # filter out to short ones
            if hmmpacbporf.get_unextended_length() < HMM_LEFTSPRDIF_PACBPORF_MIN_AA_LENGTH:
                continue

            # filter out poor bitscores
            if hmmpacbporf.bitscore < HMM_LEFTSPRDIF_PACBPORF_MIN_BITSCORE:
                continue

            # get nodes of this hmmpacbporf
            nodeQ = (OPTIONS.target,hmmpacbporf.orfQ.id)
            nodeS = (informant,hmmpacbporf.orfS.id)
            
            if inwpCBG.has_edge(nodeQ,nodeS):
                # HMM-pacbp potentially overlaps with existing PacbPORF
                overlappacbp = inwpCBG.get_pacbps_by_nodes(nodeQ,nodeS)[0]
                pdata = overlappacbp.relatively_positioned_towards(hmmpacbporf)
                # check if HMM-pacbp can result in a long enough extention
                if pdata['Q2'][0] < HMM_LEFTSPRDIF_PACBPORF_MIN_AA_LENGTH:
                    continue
                if pdata['S2'][0] < HMM_LEFTSPRDIF_PACBPORF_MIN_AA_LENGTH:
                    continue

            if prev and prev.has_edge(nodeQ,nodeS):
                # HMM-pacbp potentially overlaps with existing PacbPORF
                overlappacbp = prev.get_pacbps_by_nodes(nodeQ,nodeS)[0]
                pdata = overlappacbp.relatively_positioned_towards(hmmpacbporf)
                # check if HMM-pacbp can result in a long enough extention
                if pdata['Q2'][2] < HMM_LEFTSPRDIF_PACBPORF_MIN_AA_LENGTH:
                    continue
                if pdata['S2'][2] < HMM_LEFTSPRDIF_PACBPORF_MIN_AA_LENGTH:
                    continue

            if prev and informant in prev.organism_set():
                # HMM-pacbp potentially overlaps with existing PacbPORF
                overlappacbp = prev.get_pacbps_by_organisms(OPTIONS.target,informant)[0]
                pdata = overlappacbp.relatively_positioned_towards(hmmpacbporf)
                # check if HMM-pacbp can result in a long enough extention
                if pdata['Q2'][2] < HMM_LEFTSPRDIF_PACBPORF_MIN_AA_LENGTH:
                    continue
                if pdata['S2'][2] < HMM_LEFTSPRDIF_PACBPORF_MIN_AA_LENGTH:
                    continue

            # if here -> store hmmpacbporf to sprdif_hmmpacbporf_list
            row = (hmmpacbporf.bitscore,hmmpacbporf,informant)
            sprdif_hmmpacbporf_list.append( row )
            logf("# HMM LEFTSPRDIF APPENDED TO hmmPCG::", informant, hmmpacbporf)

    # file cleanup of HMM profile
    _file_cleanup([fname_hmmbuild])

# HMM log message
logf("#",stw.lap(),"HMM LEFTSPRDIF searches done")


#################################################################################
#### HMM RIGTHSPRDIF search for missing Orfs in inwpCBGs
#################################################################################
# do NOT recreate inwpCBGs; LEFT & RIGTH SPRDIF searches done on same inwpCBGs
logf(stw.lap(),"HMM RIGTHSPRDIF searches done")

for pos in range(0,len(inwpcbgs)):
    # check for OPTIONS.omithmm -> no HMM profile searches!
    if OPTIONS.omithmm: continue
    inwpCBG = inwpcbgs[pos]

    # check for rigth sprdif in inwpCBG
    if not inwpCBG.has_rigth_spanningrange_difference(
    sprdif_min_aa_length=HMM_RIGTHSPRDIF_MIN_AA_LENGTH):
        continue

    # get previous & next inwpCBG
    prev,next = None,None
    if pos > 0: prev = inwpcbgs[pos-1]
    if pos < len(inwpcbgs)-1: next = inwpcbgs[pos+1]

    # create RIGTHSPRDIF HMM search profile
    fname_hmmbuild, hmmcoords = _create_hmm_profile(
                inwpCBG,area="RIGTHSPRDIF",
                sprdif_min_aa_length=HMM_RIGTHSPRDIF_MIN_AA_LENGTH,
                prevcbg=prev,nextcbg=next,verbose=VERBOSE)

    # no RIGTHSPRDIF profile can be constructed
    if not fname_hmmbuild: continue
    if min([ end-sta for sta,end in hmmcoords.values()]) < HMM_RIGTHSPRDIF_PACBPORF_MIN_AA_LENGTH:
        _file_cleanup([fname_hmmbuild])
        continue

    # check if OPTIONS.target still in hmmcoords; in exceptional cases,
    # it can have been removed from the HMMprofile
    if not hmmcoords.has_key(OPTIONS.target):
        _file_cleanup([fname_hmmbuild])
        continue

    # check if this HMM profile is already evaluated (as a LEFTSPRDIF)
    ordered_coords = [ ( k,sta,end ) for k,(sta,end) in hmmcoords.iteritems() ]
    ordered_coords.sort()
    if ordered_coords in hmm_coord_check_list:
        _file_cleanup([fname_hmmbuild])
        continue

    # append hmmcoords to hmm_coord_check_list
    hmm_coord_check_list.append( ordered_coords )

    for informant in ORGANISMS:
        # continue if target or already seen informant species
        if informant == OPTIONS.target: continue

        # check if, for this informant, a potential extention exists
        if hmmcoords.has_key(informant):
            thepacbporf = inwpCBG.get_pacbps_by_organisms(OPTIONS.target,informant)[0]
            hmmrange = Set(range(hmmcoords[informant][0],hmmcoords[informant][1]))
            diff = hmmrange.difference( thepacbporf.alignment_protein_range_sbjct())
            # no large enough AA difference to find HMM PacbPORF
            if len(diff) <= 10: continue

        # create fasta search database with ORFs from this informant
        fname_fasta_db = _create_hmm_db(informant,input,inwpCBG,prev,next)
        # check for no elegiable ORFs at all
        if not fname_fasta_db: continue

        # run hmmsearch on this fasta database with hmmbuild file; maximum 2 hits
        results = hmmsearch_protein( fname_hmmbuild, fname_fasta_db, params= {'E':150,'A': 2} )
        # create pacbporfhmmlist list
        pacbporfhmmlist = hmmresults2splittedpacbps(results,hmmcoords,
                OPTIONS.target,informant,inwpCBG,input,
                gapsize=HMM_RIGTHSPRDIF_HMM_SPLIT_ON_GAPSIZE,
                min_bitscore=HMM_RIGTHSPRDIF_HMM_MIN_BITSCORE)

        # file cleanup of HMM fasta database
        _file_cleanup([fname_fasta_db])

        for hmmpacbporf in pacbporfhmmlist:
            # filter out to short ones
            if hmmpacbporf.get_unextended_length() < HMM_RIGTHSPRDIF_PACBPORF_MIN_AA_LENGTH:
                continue

            # filter out poor bitscores
            if hmmpacbporf.bitscore < HMM_RIGTHSPRDIF_PACBPORF_MIN_BITSCORE:
                continue

            # get nodes of this hmmpacbporf
            nodeQ = (OPTIONS.target,hmmpacbporf.orfQ.id)
            nodeS = (informant,hmmpacbporf.orfS.id)
            
            if inwpCBG.has_edge(nodeQ,nodeS):
                # HMM-pacbp potentially overlaps with existing PacbPORF
                overlappacbp = inwpCBG.get_pacbps_by_nodes(nodeQ,nodeS)[0]
                pdata = overlappacbp.relatively_positioned_towards(hmmpacbporf)
                # check if HMM-pacbp can result in a long enough extention
                if pdata['Q2'][2] < HMM_RIGTHSPRDIF_PACBPORF_MIN_AA_LENGTH:
                    continue
                if pdata['S2'][2] < HMM_RIGTHSPRDIF_PACBPORF_MIN_AA_LENGTH:
                    continue

            if next and next.has_edge(nodeQ,nodeS):
                # HMM-pacbp potentially overlaps with existing PacbPORF
                overlappacbp = next.get_pacbps_by_nodes(nodeQ,nodeS)[0]
                pdata = overlappacbp.relatively_positioned_towards(hmmpacbporf)
                # check if HMM-pacbp can result in a long enough extention
                if pdata['Q2'][0] < HMM_RIGTHSPRDIF_PACBPORF_MIN_AA_LENGTH:
                    continue
                if pdata['S2'][0] < HMM_RIGTHSPRDIF_PACBPORF_MIN_AA_LENGTH:
                    continue

            if next and informant in next.organism_set():
                # HMM-pacbp potentially overlaps with existing PacbPORF
                overlappacbp = next.get_pacbps_by_organisms(OPTIONS.target,informant)[0]
                pdata = overlappacbp.relatively_positioned_towards(hmmpacbporf)
                # check if HMM-pacbp can result in a long enough extention
                if pdata['Q2'][0] < HMM_RIGTHSPRDIF_PACBPORF_MIN_AA_LENGTH:
                    continue
                if pdata['S2'][0] < HMM_RIGTHSPRDIF_PACBPORF_MIN_AA_LENGTH:
                    continue

            # if here -> store hmmpacbporf to sprdif_hmmpacbporf_list
            row = (hmmpacbporf.bitscore,hmmpacbporf,informant)
            sprdif_hmmpacbporf_list.append( row )
            logf("# HMM RIGTHSPRDIF APPENDED TO hmmPCG::", informant, hmmpacbporf)

    # file cleanup of HMM profile
    _file_cleanup([fname_hmmbuild])

# HMM log message
logf("#",stw.lap(),"HMM RIGTHSPRDIF searches done")

#################################################################################
#### filter HMM LEFTSPRDIF / RIGTHSPRDIF obtained PacbPORFS in a (hmm)PCG
#################################################################################
logf(stw.lap(),"HMM pacbporf -> PCG start")

if not OPTIONS.omithmm:
    # order the hmmpacbp list (on bitscore)
    sprdif_hmmpacbporf_list.sort()
    sprdif_hmmpacbporf_list.reverse()
    hmmPCG = PacbpCollectionGraph(crossdata={},blastmatrix=BLASTP_SIMILARITY_MATRIX)
    for informant in ORGANISMS:
        barcodes = []
        if informant == OPTIONS.target: continue
        crd_key = (OPTIONS.target,informant)
        for bitscore,hmmpacbporf,_informant in sprdif_hmmpacbporf_list:
            if _informant != informant: continue
            if hmmpacbporf.barcode() in barcodes: continue
            barcodes.append( hmmpacbporf.barcode() )
            logf(informant, "\t", hmmpacbporf)
            hmmpacbporf.print_protein_and_dna()
            # create crossdata key when not existing already
            if not hmmPCG._crossdata.has_key(crd_key):
                hmmPCG._crossdata[crd_key] = deepcopy(CROSSDATA_SUBDICT)

            # add to hmmPCG for further processing
            queryNode = (OPTIONS.target,hmmpacbporf.orfQ.id)
            sbjctNode = (informant,hmmpacbporf.orfS.id)
            hmmPCG.add_node(queryNode)
            hmmPCG.add_node(sbjctNode)
            hmmPCG.add_edge(queryNode,sbjctNode,wt=hmmpacbporf.bitscore)
            #pacbpkey = hmmpacbporf.construct_unique_key(queryNode,sbjctNode)
            #hmmPCG.pacbps[(pacbpkey,queryNode,sbjctNode)] = hmmpacbporf
            crd_pacbp_key = ( hmmpacbporf.bitscore, hmmpacbporf.length,
                                hmmpacbporf.orfQ.id,hmmpacbporf.orfS.id )
            hmmPCG._crossdata[crd_key]['accepted_pacbs'][crd_pacbp_key] =\
                hmmpacbporf
    
    hmmPCG.harvest_pacbps_from_crossdata()
    logf(hmmPCG, "hmmPCG.harvest_pacbps_from_crossdata()")
    hmmPCG.remove_nonlinear_pacbs()
    logf(hmmPCG, "hmmPCG.remove_nonlinear_pacbs()")
    hmmPCG.pacbps2pacbporfs(input)
    logf(hmmPCG, "hmmPCG.pacbps2pacbporfs(input)")
    hmmPCG.remove_noncoding_pacbpdnas(verbose=True)
    logf(hmmPCG, "hmmPCG.remove_noncoding_pacbpdnas()")
    hmmPCG.remove_inclusive_pacbps()
    logf(hmmPCG, "hmmPCG.remove_inclusive_pacbps()")
    hmmPCG.pacbporfs2pacbps()
    logf(hmmPCG, "hmmPCG.pacbporfs2pacbps()")
    hmmPCG.remove_alternative_pacbps(overlap_ratio=0.5) 
    logf(hmmPCG, "hmmPCG.remove_alternative_pacbps(overlap_ratio=0.5)")
    hmmPCG.correct_overlaps(GENECOMBIS)
    logf(hmmPCG, "hmmPCG.correct_overlaps(GENECOMBIS)")
    hmmPCG.remove_short_pacbps()
    logf(hmmPCG, "hmmPCG.remove_short_pacbps()")
    
    # store hmmpacbporfs from hmmPCG to main PCG
    for (key,queryNode,sbjctNode),hmmpacbporf in hmmPCG.pacbps.iteritems():
        new_key = hmmpacbporf.construct_unique_key(queryNode,sbjctNode)
        PCG.add_node(queryNode)
        PCG.add_node(sbjctNode)
        PCG.add_edge(queryNode,sbjctNode,wt=hmmpacbporf.bitscore)
        PCG.pacbps[(new_key,queryNode,sbjctNode)] = hmmpacbporf
        # logf message for addition to main PCG
        logf("# HMM pacbporf SPRDIF ADDED::", informant, hmmpacbporf)

    # first translate PCG pacbporfs back to pacbps
    PCG.pacbporfs2pacbps()
    # get rid of introduced overlap in the PCG
    PCG.correct_overlaps(GENECOMBIS)
    # remove (very) short pacbps
    PCG.remove_short_pacbps()
    # and restore all pacbps as extended pacbporf
    PCG.pacbps2pacbporfs(input)
    # restore all pacbps to EXTENDED pacbporfs
    PCG.extend_pacbporfs(input)

    if PCG_HAS_HIGH_GAP_RATIO_PACBPORFS:
        # perform again PCG.merge_high_gap_ratio_pacbporfs() after HMM searches
        PCG_HAS_HIGH_GAP_RATIO_PACBPORFS = PCG.merge_high_gap_ratio_pacbporfs(OPTIONS.target)

    # recreate inwpcbgs
    inwpcbgs = PCG2inwpCBGS(PCG)

    # some printing in VERBOSE mode
    if VERBOSE: print_inwpcbgstructure(inwpcbgs,ORGANISMS)

# HMM log message
logf("#",stw.lap(),"HMM SPRDIF PacbPORFs stored to PCG")


#################################################################################
#### HMM MINSR search for missing Orfs in **CENTRAL** inwpCBGs
#################################################################################
HMM_MINSR_MIN_PROFILE_AA_LENGTH         = 5
HMM_MINSR_MIN_AA_LENGTH                 = None # auto defined
HMM_MINSR_HMM_MIN_BITSCORE              = -100
HMM_MINSR_HMM_SPLIT_ON_GAPSIZE          = 2
HMM_MINSR_PACBPORF_MIN_BITSCORE         = 10
HMM_MINSR_PACBPORF_MIN_AA_LENGTH        = None # auto defined

logf(stw.lap(),"HMM CENTRALCBG MINSR searches start")

is_any_added = False
for pos in range(1,len(inwpcbgs)-1):
    # check for OPTIONS.omithmm -> no HMM profile searches!
    if OPTIONS.omithmm: continue
    inwpCBG = inwpcbgs[pos]

    if not GENE_IDENTIFIER_SET.difference(inwpCBG.organism_set()):
        # no organisms/genes/nodes missing -> continue
        continue

    # get previous & next inwpCBG
    prev = inwpcbgs[pos-1]
    next = inwpcbgs[pos+1]

    # log message in VERBOSE mode
    logf(inwpCBG,inwpCBG._get_target_node())
    # create MINSR HMM search profile
    # add additional (strict) filtering on strip_nonaligned_residues
    fname_hmmbuild, hmmcoords = _create_hmm_profile(
                inwpCBG,area="MINSR",
                strip_nonaligned_residues=True,
                prevcbg=prev,nextcbg=next,verbose=VERBOSE)
    
    # no MINSR profile can be constructed
    if not fname_hmmbuild: continue

    # check if OPTIONS.target still in hmmcoords; in exceptional cases,
    # it can have been removed from the HMMprofile
    if not hmmcoords.has_key(OPTIONS.target):
        _file_cleanup([fname_hmmbuild])
        continue

    # abandon the search for really small (<5AA) profiles
    # this corresponds to (initial) MIN_PABCP_AA_LENGTH threshold
    if max([b-a for (a,b) in hmmcoords.values()]) < HMM_MINSR_MIN_PROFILE_AA_LENGTH:
        _file_cleanup([fname_hmmbuild])
        continue

    # obtain minimal length of applied sequence(s) in the profile
    min_hmm_profile_size = min([b-a for (a,b) in hmmcoords.values()])
    HMM_MINSR_PACBPORF_MIN_AA_LENGTH = max([ 5, 0.75 * min_hmm_profile_size ])
       
    for informant in ORGANISMS:
        # continue if target or already seen informant species
        if informant == OPTIONS.target: continue
        if informant in inwpCBG.organism_set(): continue

        # create fasta search database with ORFs from this informant
        fname_fasta_db = _create_hmm_db(informant,input,inwpCBG,prev,next)

        # check for no elegiable ORFs at all
        if not fname_fasta_db: continue

        # run hmmsearch on this fasta database with hmmbuild file; maximum 2 hits
        results = hmmsearch_protein( fname_hmmbuild, fname_fasta_db, params= {'E':150,'A': 2})
        if results: logf([ (res[0],res[-2]) for res in results ])

        # create pacbporfhmmlist list
        pacbporfhmmlist = hmmresults2splittedpacbps(results,hmmcoords,
                OPTIONS.target,informant,inwpCBG,input,
                gapsize=HMM_MINSR_HMM_SPLIT_ON_GAPSIZE,
                min_bitscore=HMM_MINSR_HMM_MIN_BITSCORE)
        if pacbporfhmmlist: logf([ (pf.orfS.id,pf.get_unextended_length(),pf.bitscore) for pf in pacbporfhmmlist ])

        # check/add the pacbporfs to the PCG
        thepacbporfs = pacb.ordering.order_pacbporf_list(PCG.get_pacbps_by_organisms(OPTIONS.target,informant))

        # file cleanup of HMM fasta database
        _file_cleanup([fname_fasta_db])

        for hmmpacbporf in pacbporfhmmlist:
            # filter out to short ones
            if hmmpacbporf.get_unextended_length() < HMM_MINSR_PACBPORF_MIN_AA_LENGTH:
                continue

            # filter out poor bitscores
            if hmmpacbporf.bitscore < HMM_MINSR_PACBPORF_MIN_BITSCORE:
                continue

            rejected = is_hmmpacbporf_conflicting_with_pacbporflist(
                        hmmpacbporf, thepacbporfs)

            if not rejected:
                # if here, add Node & PacbPORF to the PCG.
                hmmpacbporf2PCG(hmmpacbporf,OPTIONS.target,informant,PCG)   
                is_any_added = True
                logf("# HMM pacbporf MINSR ADDED::", informant, hmmpacbporf)

                # re-get list of pacbporfs
                thepacbporfs = pacb.ordering.order_pacbporf_list(
                        PCG.get_pacbps_by_organisms(OPTIONS.target,informant)
                        )

    # file cleanup of HMM profile
    _file_cleanup([fname_hmmbuild])

# HMM log message
logf("#", stw.lap(),"HMM CENTRALCBG MINSR searches done")

# recreate inwpcbgs if any hmmpacbporf added to PCG
if is_any_added: inwpcbgs = PCG2inwpCBGS(PCG)

# logf messages in VERBOSE mode
logf(PCG,"after regular HMM searches")
if VERBOSE: print_inwpcbgstructure(inwpcbgs,ORGANISMS)


################################################################################
### check after HMM again if high_gap_ratio_pacbporfs exist
################################################################################
PCG_HAS_HIGH_GAP_RATIO_PACBPORFS = PCG.merge_high_gap_ratio_pacbporfs(
        OPTIONS.target,verbose=False)

if PCG_HAS_HIGH_GAP_RATIO_PACBPORFS:
    # we already know there are high gap ratio pacbporfs.
    # merge all PacbPORFs which have identical target and nodesGENE_IDENTIFIERs
    # in the inwpCBGs
    enforce_node_merges = Set()
    for pos in range(1,len(inwpcbgs)):
        prevInwpCBG = inwpcbgs[pos-1]
        nextInwpCBG = inwpcbgs[pos]
        prevNodeIdSet = Set([ prevInwpCBG.get_organism_nodes(org)[0] for org in prevInwpCBG.organism_set().intersection(GENE_IDENTIFIER_SET) ])
        nextNodeIdSet = Set([ nextInwpCBG.get_organism_nodes(org)[0] for org in nextInwpCBG.organism_set().intersection(GENE_IDENTIFIER_SET) ])
        if not prevNodeIdSet.symmetric_difference(nextNodeIdSet):
            enforce_node_merges.add(prevInwpCBG._get_target_node())

    print "enforce_node_merges::", enforce_node_merges
    PCG_HAS_HIGH_GAP_RATIO_PACBPORFS = PCG.merge_high_gap_ratio_pacbporfs(
            OPTIONS.target,enforce_nodes_to_be_merged = enforce_node_merges,
            verbose=True)

    # recreate inwpcbgs if PCG_HAS_HIGH_GAP_RATIO_PACBPORFS
    inwpcbgs = PCG2inwpCBGS(PCG)

logf("#", stw.lap(),PCG,"PCG.merge_high_gap_ratio_pacbporfs()",
        PCG_HAS_HIGH_GAP_RATIO_PACBPORFS)


################################################################################
# HMM search operations done. MileStone logging
################################################################################
if not VERBOSE: stdoutManager.reset()
logf("##", stwMS.lap(), "HMM searches performed")

################################################################################
### starting next algorithm step: tinyExon InwardsPointingCodingBlockGraph II
################################################################################
if not VERBOSE: stdoutManager.buffer() # capture all STDOUT messages from deeper modules 
stwMS.stepname = "TinyExon-II"
stw.stepname = "TinyExon-II"


################################################################################
#### Find TinyExons - complex
################################################################################

if not OPTIONS.omittinyexon and not PCG_HAS_HIGH_GAP_RATIO_PACBPORFS and\
not nt_identity_quarantined_input:
    # do tinyexonsearches; OMIT when PCG_HAS_HIGH_GAP_RATIO_PACBPORFS
    # or some informants have very high identity
    tinyexonpacbporfs = discover_tinyexons_complex(PCG,input,OPTIONS.target,GENE_INFORMANT_SET,UNIGENE_INFORMANT_SET)

    # translate tinyexonpacbporfs to tinyexonInwpCBGs
    tinyexonInwpCBGs = tinyexonpacbporfs2InwardsPointingCodingBlockGraphs(OPTIONS.target,tinyexonpacbporfs)
    
    is_any_added = False
    for tinyExonInwpCBG in tinyexonInwpCBGs:
        logf(tinyExonInwpCBG, tinyExonInwpCBG.get_bitscore())
        inconsistent_pacbporf_data = {}

        MINSR_tinyexon = tinyExonInwpCBG.minimal_spanning_range(organism=OPTIONS.target)
        is_novel_tinyexoncbg_accepted = True
        for inwp in inwpcbgs:
            analysesA = inwp.is_positioned_compatibly(tinyExonInwpCBG)
            analysesO = inwp.is_including(tinyExonInwpCBG)
            MINSR_current = inwp.minimal_spanning_range(organism=OPTIONS.target)
            logf(inwp, inwp.get_bitscore(), analysesA, analysesO)

            if analysesO.count(False) == 0:
                # tinyExonInwpCBG is COMPLETELY included in an
                # already existing inwpCBG
                logf("OVERLAP::", inwp.get_organism_nodes(OPTIONS.target),
                    tinyExonInwpCBG.get_organism_nodes(OPTIONS.target))
                # set is_novel_tinyexoncbg_accepted to False and break out
                is_novel_tinyexoncbg_accepted = False
                break

            if not MINSR_tinyexon.difference(MINSR_current) and\
            inwp.get_organism_nodes(OPTIONS.target) ==\
            tinyExonInwpCBG.get_organism_nodes(OPTIONS.target):
                # novel tinyexon CBG is already included as a
                # longer exon in the PCG. Omit it here
                logf( "OVERLAP::", inwp.get_organism_nodes(OPTIONS.target),
                        tinyExonInwpCBG.get_organism_nodes(OPTIONS.target) )
                is_novel_tinyexoncbg_accepted = False
                break


            # Check Orf transitions.
            # X - Y - X is NOT a likely Orf transition
            # This happens easily in existing inwpCBGs with high identity.
            # an out-of-frame (and PacbPORF-including) exon is proposed
            if analysesO.count(False) > 0:
                for pos in range(1,len(inwpcbgs)):
                    prevCBG,nextCBG = inwpcbgs[pos-1:pos+1]
                    prevCBGminsr = prevCBG.minimal_spanning_range(organism=OPTIONS.target)
                    nextCBGminsr = nextCBG.minimal_spanning_range(organism=OPTIONS.target)
                    if min(MINSR_tinyexon) >= max(prevCBGminsr) and\
                    max(MINSR_tinyexon) <= min(nextCBGminsr):
                        # tinyexonCBG is putatively placeable in between these
                        # check Orfids
                        identicalorfids = []
                        for org in prevCBG.organism_set().intersection(nextCBG.organism_set()).intersection(GENE_IDENTIFIER_SET):
                            nodePrev = prevCBG.node_by_organism(org)
                            nodeNext = nextCBG.node_by_organism(org)
                            if nodePrev==nodeNext:
                                if org in tinyExonInwpCBG.organism_set() and\
                                tinyExonInwpCBG.node_by_organism(org) in\
                                [ nodePrev, nodeNext ]:
                                    # tiny exon in another species, but hitting
                                    # to this longer (continious) exon in this
                                    # species. NOT an unlikely Orf transition
                                    identicalorfids.append(False)
                                else:
                                    # highly unlikely orf transition
                                    identicalorfids.append(True)
                            else:
                                identicalorfids.append(False)

                        # if there are only `a few` Trues, this is a strong
                        # signal for a false call
                        trues  = identicalorfids.count(True)
                        falses = identicalorfids.count(False)
                        if trues > falses/2:
                            logf( "UNLIKELY ORF TRANSITION::",identicalorfids )
                            # reject this tinyexon -> break out
                            is_novel_tinyexoncbg_accepted = False
                            break
                # check if this tinyExonInwpCBG is rejected
                if not is_novel_tinyexoncbg_accepted:
                    break

            # so some further filtering
            if not analysesA: continue
            if analysesA.count(False)==0: continue
    
        # check if we can add this tinyexon inwpCBG to the main PCG
        if is_novel_tinyexoncbg_accepted:
            logf("# TinyExon inwpCBG ADD TO PCG:", tinyExonInwpCBG,
                tinyExonInwpCBG.get_nodes(), PCG)
            is_any_added = True
            for key in inconsistent_pacbporf_data.keys(): _delete_pacbp(PCG,key)
            del( inconsistent_pacbporf_data )
            for (pacbpkey,nodeQ,nodeS),pacbporf in tinyExonInwpCBG.pacbps.iteritems():
                if not nodeQ in PCG.get_nodes(): PCG.add_node(nodeQ)
                if not nodeS in PCG.get_nodes(): PCG.add_node(nodeS)
                PCG.add_edge(nodeQ,nodeS,wt=pacbporf.bitscore)
                # label this tinyExonPacbpORF as non-deletable!
                pacbporf._IS_DELETE_PROTECTED = True
                # and store to the PCG
                PCG.pacbps[(pacbpkey,nodeQ,nodeS)] = pacbporf
    
    # recreate inwpCBGs after TinyExon inwpCBGs addition
    if is_any_added: inwpcbgs = PCG2inwpCBGS(PCG)

    if VERBOSE:
        # logf message in VERBOSE mode after addition
        logf(PCG)
        print_inwpcbgstructure(inwpcbgs,ORGANISMS)

logf("#", stw.lap(),"TinyExon - complex ready")

#################################################################################
### HMM MINSR search for missing Orfs in **CENTRAL** TINYEXON inwpCBGs
### here, we try to find additional alignments underneath just found
### TinyExon InwpCBGs only
#################################################################################
HMM_MINSR_MIN_PROFILE_AA_LENGTH         = 2
HMM_MINSR_MIN_AA_LENGTH                 = None # auto defined
HMM_MINSR_HMM_MIN_BITSCORE              = -100
HMM_MINSR_HMM_SPLIT_ON_GAPSIZE          = 1
HMM_MINSR_PACBPORF_MIN_BITSCORE         = 1
HMM_MINSR_PACBPORF_MIN_AA_LENGTH        = None # auto defined


is_any_added = False
for pos in range(1,len(inwpcbgs)-1):
    # check for OPTIONS.omithmm -> no TinyExon searches
    if OPTIONS.omittinyexon: continue
    # check for OPTIONS.omithmm -> no HMM profile searches!
    if OPTIONS.omithmm: continue

    inwpCBG = inwpcbgs[pos]

    # check if obtained by TINYEXON prediction
    if [ pf._IS_DELETE_PROTECTED for pf in inwpCBG.pacbps.values() ].count(True)==0:
        continue

    if not GENE_IDENTIFIER_SET.difference(inwpCBG.organism_set()):
        # no organisms/genes/nodes missing -> continue
        continue

    # get previous & next inwpCBG
    prev = inwpcbgs[pos-1]
    next = inwpcbgs[pos+1]

    # log message for inwpCBG
    logf(inwpCBG,inwpCBG._get_target_node())
    # create MINSR HMM search profile
    # **OMIT** additional (strict) filtering on strip_nonaligned_residues
    fname_hmmbuild, hmmcoords = _create_hmm_profile(
                inwpCBG,area="MINSR",
                prevcbg=prev,nextcbg=next,verbose=VERBOSE)
    
    # no MAXSR profile can be constructed
    if not fname_hmmbuild: continue

    # check if OPTIONS.target still in hmmcoords; in exceptional cases,
    # it can have been removed from the HMMprofile
    if not hmmcoords.has_key(OPTIONS.target):
        _file_cleanup([fname_hmmbuild])
        continue

    # abandon the search for really small (<3AA) profiles
    # this corresponds to (initial) MIN_PABCP_AA_LENGTH threshold
    if max([b-a for (a,b) in hmmcoords.values()]) < HMM_MINSR_MIN_PROFILE_AA_LENGTH:
        _file_cleanup([fname_hmmbuild])
        continue

    # obtain minimal length of applied sequence(s) in the profile
    min_hmm_profile_size = min([b-a for (a,b) in hmmcoords.values()])
    HMM_MINSR_PACBPORF_MIN_AA_LENGTH = max([HMM_MINSR_MIN_PROFILE_AA_LENGTH,min_hmm_profile_size-1])

    for informant in ORGANISMS:
        # continue if target or already seen informant species
        if informant == OPTIONS.target: continue
        if informant in inwpCBG.organism_set(): continue

        # create fasta search database with ORFs from this informant
        fname_fasta_db = _create_hmm_db(informant,input,inwpCBG,prev,next)

        # check for no elegiable ORFs at all
        if not fname_fasta_db: continue

        # run hmmsearch on this fasta database with hmmbuild file; maximum 2 hits
        results = hmmsearch_protein( fname_hmmbuild, fname_fasta_db, params= {'E':150,'A': 2})
        if results: logf([ (res[0],res[-2]) for res in results ])

        # create pacbporfhmmlist list
        pacbporfhmmlist = hmmresults2splittedpacbps(results,hmmcoords,
                OPTIONS.target,informant,inwpCBG,input,
                gapsize=HMM_MINSR_HMM_SPLIT_ON_GAPSIZE,
                min_bitscore=HMM_MINSR_HMM_MIN_BITSCORE)
        if pacbporfhmmlist: logf([ (pf.orfS.id,pf.get_unextended_length(),pf.bitscore) for pf in pacbporfhmmlist ])

        # check/add the pacbporfs to the PCG
        thepacbporfs = pacb.ordering.order_pacbporf_list(PCG.get_pacbps_by_organisms(OPTIONS.target,informant))

        # file cleanup of HMM fasta database
        _file_cleanup([fname_fasta_db])


        for hmmpacbporf in pacbporfhmmlist:
            # filter out to short ones
            if hmmpacbporf.get_unextended_length() < HMM_MINSR_PACBPORF_MIN_AA_LENGTH:
                continue

            # filter out poor bitscores
            if hmmpacbporf.bitscore < HMM_MINSR_PACBPORF_MIN_BITSCORE:
                continue

            rejected = is_hmmpacbporf_conflicting_with_pacbporflist(
                        hmmpacbporf, thepacbporfs)

            if not rejected:
                # if here, add Node & PacbPORF to the PCG.
                hmmpacbporf2PCG(hmmpacbporf,OPTIONS.target,informant,PCG)   
                is_any_added = True
                logf("# HMM pacbporf MINSR TinyExon ADDED::", informant, hmmpacbporf)
                # re-get list of pacbporfs
                thepacbporfs = pacb.ordering.order_pacbporf_list(
                        PCG.get_pacbps_by_organisms(OPTIONS.target,informant)
                        )

        # logf message in VERBOSE mode
        logf(inwpCBG, "CENTRALCBG MINSR TINYEXON", informant, len(results), len(pacbporfhmmlist))

    # file cleanup of HMM profile
    _file_cleanup([fname_hmmbuild])

    # recreate inwpCBGs after TinyExon inwpCBGs addition
    if is_any_added: inwpcbgs = PCG2inwpCBGS(PCG)

    if VERBOSE:
        # logf message in VERBOSE mode after addition
        logf(PCG)
        print_inwpcbgstructure(inwpcbgs,ORGANISMS)

# HMM log message
logf("#", stw.lap(),"HMM CENTRALCBG TINYEXON MINSR searches done")


################################################################################
# Find TinyExons - complex. MileStone logging
################################################################################
if not VERBOSE: stdoutManager.reset()
logf("##", stwMS.lap(), "Find TinyExons - complex: performed")


################################################################################
### starting next algorithm step: InwpCBG removal
################################################################################
if not VERBOSE: stdoutManager.buffer() # capture all STDOUT messages from deeper modules 
stwMS.stepname = "InwpCbgRemoval"
stw.stepname = "InwpCbgRemoval"
if not True in (QUIET,SILENT): print_inwpcbgstructure(inwpcbgs,ORGANISMS)


################################################################################
### Detect putative micro-syntheny: do Target and part of the Informants
### share neighbouring genes (in the same orientation?). In this case,
### removing poor informants will likely remove all informants which do
### NOT have this microsyntheny. This is bad. So, here we try to detect these
### putative microsynthenic gene orders and remove them. As a (very positive)
### consequence, the chance of predicting (False) gene fusions is highly
### decreased!
################################################################################
syntenicPCG = detect_and_remove_synteny(inwpcbgs,PCG,GENE_IDENTIFIER_SET)

if syntenicPCG:
    logf("# synteny in gene order detected; %s" % syntenicPCG)

    # check if an informant is completely deleted from the PCG
    for informant in Set(ORGANISMS).difference(PCG.organism_set()):
        status = removeinformantfrominput(informant,input)

    # recreate data structures
    ( ORGANISMS, GENECOMBIS,
      GENE_INFORMANT_SET,
      GENE_IDENTIFIER_SET,
      UNIGENE_INFORMANT_SET) = recreate_global_informantdata(PCG,OPTIONS.target)

    # check if an informant is not in ORGANISMS but still in input
    for informant in Set(input.keys()).difference(ORGANISMS):
        status = removeinformantfrominput(informant,input)

    # recreate inwpCBGs
    inwpcbgs = PCG2inwpCBGS(PCG)
    # print inwpCBG structure
    if VERBOSE: print_inwpcbgstructure(inwpcbgs,ORGANISMS)

# log message for syntenic gene order check
logf("#", stw.lap(), "syntenic gene order checked on gene loci (I):", syntenicPCG)


################################################################################
### remove gtg2GTG discrepancies
################################################################################
gtgdiscrepancyPCG = detect_and_remove_gtgdiscrepancy(inwpcbgs,PCG,GENE_IDENTIFIER_SET)

if gtgdiscrepancyPCG:
    logf("# gtgdiscrepancy in inwpCBGs detected; %s" % gtgdiscrepancyPCG)

    # check if an informant is completely deleted from the PCG
    for informant in Set(ORGANISMS).difference(PCG.organism_set()):
        status = removeinformantfrominput(informant,input)

    # recreate data structures
    ( ORGANISMS, GENECOMBIS,
      GENE_INFORMANT_SET,
      GENE_IDENTIFIER_SET,
      UNIGENE_INFORMANT_SET) = recreate_global_informantdata(PCG,OPTIONS.target)

    # check if an informant is not in ORGANISMS but still in input
    for informant in Set(input.keys()).difference(ORGANISMS):
        status = removeinformantfrominput(informant,input)

    # recreate inwpCBGs
    inwpcbgs = PCG2inwpCBGS(PCG)
    # print inwpCBG structure
    if VERBOSE: print_inwpcbgstructure(inwpcbgs,ORGANISMS)

# log message for gtgdiscrepancy check
logf("#", stw.lap(), "gtg2GTG discrepancy inwpCBGs detected (I):", gtgdiscrepancyPCG)


################################################################################
#### make (ordered list of) InwardsPointingCodingBlockGraph and try to
#### remove (noncoding) potention 5' and 3' UTR alignments
################################################################################
ncngPCG = detect_and_remove_utrornonegene_inwpcbgs(inwpcbgs,PCG)

if ncngPCG:
    logf("# non-coding / non-gene inwpCBGs detected; %s" % ncngPCG)

    # check if an informant is completely deleted from the PCG
    for informant in Set(ORGANISMS).difference(PCG.organism_set()):
        status = removeinformantfrominput(informant,input)

    # recreate data structures
    ( ORGANISMS, GENECOMBIS,
      GENE_INFORMANT_SET,
      GENE_IDENTIFIER_SET,
      UNIGENE_INFORMANT_SET) = recreate_global_informantdata(PCG,OPTIONS.target)

    # check if an informant is not in ORGANISMS but still in input
    for informant in Set(input.keys()).difference(ORGANISMS):
        status = removeinformantfrominput(informant,input)

    # recreate inwpCBGs
    inwpcbgs = PCG2inwpCBGS(PCG)
    # print inwpCBG structure
    if VERBOSE: print_inwpcbgstructure(inwpcbgs,ORGANISMS)

# log message for utrornonegene check
logf("#", stw.lap(), "utrornonegene inwpCBGs detected (I):", ncngPCG)


################################################################################
### place back quarantined informants
################################################################################

if nt_identity_quarantined_input:
    # place back the quarantined PacbPORFs and ndoes
    for (pacbpkey,nodeQ,nodeS),pacbporf in quarantinePCG.pacbps.iteritems():
        PCG.pacbps[(pacbpkey,nodeQ,nodeS)] = pacbporf
        PCG.add_node(nodeQ)
        PCG.add_node(nodeS)
        PCG.add_edge(nodeQ,nodeS,wt=pacbporf.bitscore)
    # update the main input data structure for the quarantined informants
    for informant in nt_identity_quarantined_input.keys():
        input.update(nt_identity_quarantined_input)
        logf("# quarantined informant (%s) placed back" % informant )

    # recreate data structures
    ( ORGANISMS, GENECOMBIS,
      GENE_INFORMANT_SET,
      GENE_IDENTIFIER_SET,
      UNIGENE_INFORMANT_SET) = recreate_global_informantdata(PCG,OPTIONS.target)

    # recreate inwpCBGs
    inwpcbgs = PCG2inwpCBGS(PCG)

    # remove inwpCBGs which are only covered by target & quarantined informants
    organism_checklist = nt_identity_quarantined_input.keys()
    organism_checklist.append( OPTIONS.target )
    is_any_deleted = False
    for inwpCBG in inwpcbgs:
        if not inwpCBG.organism_set().symmetric_difference(organism_checklist):
            is_any_deleted = True
            for pacbpkey in inwpCBG.pacbps.keys():
                _delete_pacbp(PCG,pacbpkey)
    if is_any_deleted:
        # check if an informant is completely deleted from the PCG
        for informant in Set(ORGANISMS).difference(PCG.organism_set()):
            status = removeinformantfrominput(informant,input)
            if informant in nt_identity_quarantined_input.keys():
                # yes, this can happen: high-nt identity informant is deleted
                # this can happen when identity has been found OUTSIDE ``main``
                # genemodel (putatively caused by False gene fusion prediction?)
                # Remove this informant from the nt_identity_quarantined_input dict
                del( nt_identity_quarantined_input[informant] )

        # recreate data structures
        ( ORGANISMS, GENECOMBIS,
          GENE_INFORMANT_SET,
          GENE_IDENTIFIER_SET,
          UNIGENE_INFORMANT_SET) = recreate_global_informantdata(PCG,OPTIONS.target)

        # recreate inwpCBGs
        inwpcbgs = PCG2inwpCBGS(PCG)

    # find blocks which seem NOT to be covered by high-identity informants,
    # and do an HMM search on these. This is (expected to be) a rare event.
    # Make shure there is not to strict filtering on surrounding inwpCBGs,
    # because it is likely caused by an out-of-frame PacbPORF which rejected
    # this Orf to be picked up (filtered out on overlap in PCG)
    is_any_added = False
    replace_pacbporf_key_list = []

    for inwpCBG in inwpcbgs:
        if not nt_identity_quarantined_input:
            # placed back informant is deleted on the moment that it was placed back
            continue

        if not Set(nt_identity_quarantined_input.keys()).difference(inwpCBG.organism_set()):
            # no organisms/genes/nodes missing -> continue
            continue
    
        # get previous & next inwpCBG
        prev,next = None,None
    
        # create MAXSR HMM search profile
        fname_hmmbuild, hmmcoords = _create_hmm_profile(
                    inwpCBG,area="MAXSR",
                    prevcbg=prev,nextcbg=next,verbose=VERBOSE)
        # no MAXSR profile can be constructed
        if not fname_hmmbuild: continue
        if min([ end-sta for sta,end in hmmcoords.values()]) < HMM_MAXSR_PACBPORF_MIN_AA_LENGTH:
            _file_cleanup([fname_hmmbuild])
            continue

        # check if OPTIONS.target still in hmmcoords; in exceptional cases,
        # it can have been removed from the HMMprofile
        if not hmmcoords.has_key(OPTIONS.target):
            _file_cleanup([fname_hmmbuild])
            continue
    
    
        for informant in nt_identity_quarantined_input.keys():
            # continue if target or already seen informant species
            if informant == OPTIONS.target: continue
            if informant in inwpCBG.organism_set(): continue
    
            # create fasta search database with ORFs from this informant
            fname_fasta_db = _create_hmm_db(informant,input,inwpCBG,prev,next)
            # check for no elegiable ORFs at all
            if not fname_fasta_db: continue
    
            # run hmmsearch on this fasta database with hmmbuild file; maximum 2 hits
            results = hmmsearch_protein( fname_hmmbuild, fname_fasta_db, params= {'E':150,'A': 2})
    
            # create pacbporfhmmlist list
            pacbporfhmmlist = hmmresults2splittedpacbps(results,hmmcoords,
                    OPTIONS.target,informant,inwpCBG,input,
                    gapsize=HMM_MAXSR_HMM_SPLIT_ON_GAPSIZE,
                    min_bitscore=HMM_MAXSR_HMM_MIN_BITSCORE)

            # check/add the pacbporfs to the PCG
            thepacbporfs = pacb.ordering.order_pacbporf_list(PCG.get_pacbps_by_organisms(OPTIONS.target,informant))
    
            # file cleanup of HMM fasta database
            _file_cleanup([fname_fasta_db])
    
            for hmmpacbporf in pacbporfhmmlist:
                # filter out to short ones
                if hmmpacbporf.get_unextended_length() < HMM_MAXSR_PACBPORF_MIN_AA_LENGTH:
                    continue

                # filter out high-nt-identity pacbporfs
                if hmmpacbporf.get_nt_identity() < MIN_QUARANTINE_HIGH_NT_IDENTITY:
                    continue
    
                # filter out poor bitscores
                if hmmpacbporf.bitscore < HMM_MAXSR_PACBPORF_MIN_BITSCORE:
                    continue
    
                rejected = is_hmmpacbporf_conflicting_with_pacbporflist(
                            hmmpacbporf, thepacbporfs)

                ################################################################
                if VERBOSE:
                    logf( hmmpacbporf, hmmpacbporf.get_nt_identity(),
                            hmmpacbporf.is_coding(), " REJECTED::", rejected )
                    hmmpacbporf.print_protein_and_dna()
                ################################################################
    
                if not rejected:
                    # if here, add Node & PacbPORF to the PCG.
                    hmmpacbporf2PCG(hmmpacbporf,OPTIONS.target,informant,PCG)   
                    is_any_added = True
                    logf("# HMM pacbporf MAXSR high-ntident ADDED::", informant, hmmpacbporf)
               
                    # re-get list of pacbporfs
                    thepacbporfs = pacb.ordering.order_pacbporf_list(
                            PCG.get_pacbps_by_organisms(OPTIONS.target,informant)
                            )
                else:
                    # see if we can replace this PacbPORF for another one
                    for currpacbporf in thepacbporfs:
                        if currpacbporf.overlap(hmmpacbporf) == 1.0 and\
                        currpacbporf.orfS.id != hmmpacbporf.orfS.id and\
                        hmmpacbporf.identityscore >= currpacbporf.identityscore:
                            the_replace_key = None
                            for (pacbpkey,nodeQ,nodeS) in PCG.pacbps.keys():
                                if nodeS == (informant,currpacbporf.orfS.id):
                                    the_replace_key = (pacbpkey,nodeQ,nodeS)
                                    replace_pacbporf_key_list.append(
                                        (the_replace_key,hmmpacbporf,informant)
                                        )
                                    break

            # logf message in verbose mode
            logf(inwpCBG, "MAXSR high nt-identity", informant, len(results), len(pacbporfhmmlist))
    
        # file cleanup of HMM profile
        _file_cleanup([fname_hmmbuild])

    # Check if data in replace_pacbporf_key_list
    # First element is pacbpkey which can be deleted
    # Second element is pacbporf to replace this one
    # Third element is the informant
    if replace_pacbporf_key_list:
        for (delete_key,new_hmm_pacbporf,informant) in replace_pacbporf_key_list:
            _delete_pacbp(PCG,delete_key)
            # if here, add Node & PacbPORF to the PCG.
            hmmpacbporf2PCG(hmmpacbporf,OPTIONS.target,informant,PCG) 

    # HMM log message
    logf(stw.lap(),"HMM MAXSR high-nt-identity searches done")    

    # recreate inwpcbgs in case a PacbPORF was added
    if is_any_added or replace_pacbporf_key_list: inwpcbgs = PCG2inwpCBGS(PCG)

    if VERBOSE:
        logf(PCG,"after high-ntidentity adding & HMM")
        print_inwpcbgstructure(inwpcbgs,ORGANISMS)

# log message high nt-identity quarantine placed back
logf("#", stw.lap(), "high nt-identity informants placed back:",
    quarantinePCG.organism_set())



################################################################################
# Check if OPTIONS.minimal_num_loci is still valid
################################################################################
if len(input) < OPTIONS.minimal_num_loci or\
len(GENE_IDENTIFIER_SET) < OPTIONS.minimal_num_loci or\
(OPTIONS.target and OPTIONS.target not in input.keys()):
    message = "## %s < %s (--minimal_num_loci) or target (%s) identifier removed" % (
        len(input), OPTIONS.minimal_num_loci, OPTIONS.target
        )
    # done -> abgpsysexit()
    logf(message)
    abgpsysexit(input,OPTIONS,message=message)



################################################################################
### remove gtg2GTG discrepancies
################################################################################
gtgdiscrepancyPCG = detect_and_remove_gtgdiscrepancy(inwpcbgs,PCG,GENE_IDENTIFIER_SET)

if gtgdiscrepancyPCG:
    logf("# gtgdiscrepancy in inwpCBGs detected; %s" % gtgdiscrepancyPCG)

    # check if an informant is completely deleted from the PCG
    for informant in Set(ORGANISMS).difference(PCG.organism_set()):
        status = removeinformantfrominput(informant,input)

    # recreate data structures
    ( ORGANISMS, GENECOMBIS,
      GENE_INFORMANT_SET,
      GENE_IDENTIFIER_SET,
      UNIGENE_INFORMANT_SET) = recreate_global_informantdata(PCG,OPTIONS.target)

    # check if an informant is not in ORGANISMS but still in input
    for informant in Set(input.keys()).difference(ORGANISMS):
        status = removeinformantfrominput(informant,input)

    # recreate inwpCBGs
    inwpcbgs = PCG2inwpCBGS(PCG)
    # print inwpCBG structure
    if VERBOSE: print_inwpcbgstructure(inwpcbgs,ORGANISMS)

# log message for gtgdiscrepancy check
logf("#", stw.lap(), "gtg2GTG discrepancy inwpCBGs detected (II):", gtgdiscrepancyPCG)


################################################################################
### Check AGAIN for putative micro-syntheny: do Target and part of the Informants
### share neighbouring genes (in the same orientation?). In this case,
### removing poor informants will likely remove all informants which do
### NOT have this microsyntheny. This is bad. So, here we try to detect these
### putative microsynthenic gene orders and remove them. As a (very positive)
### consequence, the chance of predicting (False) gene fusions is highly
### decreased!
################################################################################
syntenicPCG = detect_and_remove_synteny(inwpcbgs,PCG,GENE_IDENTIFIER_SET)

if syntenicPCG:
    logf("# synteny in gene order detected; %s" % syntenicPCG)

    # check if an informant is completely deleted from the PCG
    for informant in Set(ORGANISMS).difference(PCG.organism_set()):
        status = removeinformantfrominput(informant,input)

    # recreate data structures
    ( ORGANISMS, GENECOMBIS,
      GENE_INFORMANT_SET,
      GENE_IDENTIFIER_SET,
      UNIGENE_INFORMANT_SET) = recreate_global_informantdata(PCG,OPTIONS.target)

    # check if an informant is not in ORGANISMS but still in input
    for informant in Set(input.keys()).difference(ORGANISMS):
        status = removeinformantfrominput(informant,input)

    # recreate inwpCBGs
    inwpcbgs = PCG2inwpCBGS(PCG)
    # print inwpCBG structure
    if VERBOSE: print_inwpcbgstructure(inwpcbgs,ORGANISMS)

# log message for syntenic gene order check
logf("#", stw.lap(), "syntenic gene order checked on gene loci (II):", syntenicPCG)


################################################################################
### detect_and_remove_single_nonfinal_inwpcbg
################################################################################
nonfinalPCG = detect_and_remove_single_nonfinal_inwpcbg(inwpcbgs,PCG,
                GENE_IDENTIFIER_SET)

if nonfinalPCG:
    # check if an informant is completely deleted from the PCG
    for informant in Set(ORGANISMS).difference(PCG.organism_set()):
        status = removeinformantfrominput(informant,input)

    # recreate data structures
    ( ORGANISMS, GENECOMBIS,
      GENE_INFORMANT_SET,
      GENE_IDENTIFIER_SET,
      UNIGENE_INFORMANT_SET) = recreate_global_informantdata(PCG,OPTIONS.target)

    # check if an informant is not in ORGANISMS but still in input
    for informant in Set(input.keys()).difference(ORGANISMS):
        status = removeinformantfrominput(informant,input)

    # recreate inwpCBGs
    inwpcbgs = PCG2inwpCBGS(PCG)

# log message for detect_and_remove_single_nonfinal_inwpcbg
logf("#", stw.lap(), "detect_and_remove_single_nonfinal_inwpcbg:", nonfinalPCG)


################################################################################
### detect_and_remove_single_nonfirst_inwpcbg
################################################################################
nonfirstPCG = detect_and_remove_single_nonfirst_inwpcbg(inwpcbgs,PCG,
                GENE_IDENTIFIER_SET)

if nonfirstPCG:
    # check if an informant is completely deleted from the PCG
    for informant in Set(ORGANISMS).difference(PCG.organism_set()):
        status = removeinformantfrominput(informant,input)

    # recreate data structures
    ( ORGANISMS, GENECOMBIS,
      GENE_INFORMANT_SET,
      GENE_IDENTIFIER_SET,
      UNIGENE_INFORMANT_SET) = recreate_global_informantdata(PCG,OPTIONS.target)

    # check if an informant is not in ORGANISMS but still in input
    for informant in Set(input.keys()).difference(ORGANISMS):
        status = removeinformantfrominput(informant,input)

    # recreate inwpCBGs
    inwpcbgs = PCG2inwpCBGS(PCG)

# log message for detect_and_remove_single_nonfirst_inwpcbg
logf("#", stw.lap(), "detect_and_remove_single_nonfirst_inwpcbg:", nonfirstPCG)



################################################################################
### remove informants that are after microsynteny removal VERY POOR contributors....
################################################################################
if OPTIONS.disallow_informant_deletion == False:

    informant2summedbitscore = []
    for informant in ORGANISMS:
        if informant == OPTIONS.target: continue
        if not informant in PCG.organism_set():
            # apparently deleted in previous steps
            summedbits = 0
            summedlength = 0
            is_gene_informant = False
        else:
            pacbporfs = PCG.get_pacbps_by_organisms(OPTIONS.target,informant)
            if not pacbporfs:
                # apparently deleted in previous steps
                summedbits = 0
                summedlength = 0
                is_gene_informant = False
            else:
                summedbits = sum([pf.bitscore for pf in pacbporfs])
                summedlength = sum([ pf.get_unextended_length() for pf in pacbporfs ])
                if hasattr(pacbporfs[0].orfS,ORF_IS_UNIGENE_LABEL):
                    is_gene_informant = False
                else:
                    is_gene_informant = True
        # store to informant2summedbitscore
        informant2summedbitscore.append( ( summedbits, summedlength, informant, is_gene_informant) )

    informant2summedbitscore.sort()
    informant2summedbitscore.reverse()

    # logf message
    for pos in range(0,len(informant2summedbitscore)):
        summedbits, summedlength, informant, is_gene_informant = informant2summedbitscore[pos]
        logf("poorinformant analyses:", informant, summedbits, summedlength, is_gene_informant)

    # find the **best** non-high-nt-identity AND & non-unigene (so gene) informant
    best_gene_informant_pos = 0
    for pos in range(0,len(informant2summedbitscore)):
        summedbits, summedlength, informant, is_gene_informant = informant2summedbitscore[pos]
        if is_gene_informant and informant not in nt_identity_quarantined_input.keys():
            best_gene_informant_pos = pos
            break

    INFORMANT_REMOVAL_BITSCORE_RATIO = 0.15
    INFORMANT_REMOVAL_PACBPORFLENGTH_RATIO = 0.50
    INFORMANT_REMOVAL_BITSCORE_FOR_SHORT_LENGTH_RATIO = 0.40
    INFORMANT_REMOVAL_LENGTH_RATIO_PROTECTION = 0.80
    best_informant_bitscore = informant2summedbitscore[best_gene_informant_pos][0]
    best_informant_length   = informant2summedbitscore[best_gene_informant_pos][1]

    status = False
    for (score, length, informant, is_gene_informant) in informant2summedbitscore:
        score_ratio  = float(score)/best_informant_bitscore
        length_ratio = float(length)/best_informant_length
        remove_informant = False

        if is_gene_informant:
            if score_ratio < INFORMANT_REMOVAL_BITSCORE_RATIO and\
            length_ratio < INFORMANT_REMOVAL_LENGTH_RATIO_PROTECTION:
                remove_informant = True
        
            if length_ratio < INFORMANT_REMOVAL_PACBPORFLENGTH_RATIO and\
            score_ratio < INFORMANT_REMOVAL_BITSCORE_FOR_SHORT_LENGTH_RATIO and\
            length_ratio < INFORMANT_REMOVAL_LENGTH_RATIO_PROTECTION:
                remove_informant = True
        else:
            # informant is a (partial) unigene
            # only filter for poor similarity, NOT for short length coverage
            if score_ratio < INFORMANT_REMOVAL_BITSCORE_RATIO:
                remove_informant = True

            # but, do check how much of the raw unigene is covered
            # when this is very short -> remove it. This ratio can be obscured
            # by (long) UTRs in a partial unigene. But, it is safer to remove
            # it (because it can be as well a case of a very dissimilar gene locus)
            if float(length*3) / float(len(input[informant]['genomeseq'])) <\
            INFORMANT_REMOVAL_PACBPORFLENGTH_RATIO and\
            length_ratio < INFORMANT_REMOVAL_PACBPORFLENGTH_RATIO:
                remove_informant = True


        if remove_informant:
            # delete this informant from the input structure
            del( input[informant] )
            # remove its PacbPORFS
            delete_pacbp_keys = []
            for (pacbpkey,nodeQ,nodeS),pacbporf in PCG.pacbps.iteritems():
                if PCG.organism_by_node(nodeS) == informant:
                    delete_pacbp_keys.append((pacbpkey,nodeQ,nodeS))
            for pacbpkey in delete_pacbp_keys: _delete_pacbp(PCG,pacbpkey)
            # remove its nodes
            if informant in PCG.organism_set():
                for node in PCG.get_organism_nodes(informant): PCG.del_node(node)
            # some printing
            logf("# INFORMANT DELETED:", informant, is_gene_informant, "bits:",score_ratio, "ll:",length_ratio)
            status = True
        else:
            logf("# informant maintained:", informant, is_gene_informant, "bits:",score_ratio, "ll:",length_ratio)


    # recreate some global vars (if a block was deleted)
    if status:
        # recreate data structures
        ( ORGANISMS, GENECOMBIS,
          GENE_INFORMANT_SET,
          GENE_IDENTIFIER_SET,
          UNIGENE_INFORMANT_SET) = recreate_global_informantdata(PCG,OPTIONS.target)

        # check if an informant is not in ORGANISMS but still in input
        for informant in Set(input.keys()).difference(ORGANISMS):
            status = removeinformantfrominput(informant,input)
            
        # recreate inwpCBGs
        inwpcbgs = PCG2inwpCBGS(PCG)

    # log message for poor informant removal
    logf("#", stw.lap(), "poor informants removed")


################################################################################
# Check if OPTIONS.minimal_num_loci is still valid
################################################################################
if len(input) < OPTIONS.minimal_num_loci or\
(OPTIONS.target and OPTIONS.target not in input.keys()):
    message = "## %s < %s (--minimal_num_loci) or target (%s) identifier removed" % (
        len(input), OPTIONS.minimal_num_loci, OPTIONS.target
        )
    # done -> abgpsysexit()
    logf(message)
    abgpsysexit(input,OPTIONS,message=message)


################################################################################
### recreate blocks
################################################################################
blocks = PCG2blocks(PCG, verbose = False )


################################################################################
### remove_singleton_blocks_without_intron_signal
################################################################################

status = remove_singleton_blocks_without_intron_signal(blocks,PCG,ORGANISMS,OPTIONS)

# REDO translate to blocks (if a block was deleted)
if status:
    print "### remove_singleton_blocks_without_intron_signal"
    # check if an informant is completely deleted from the PCG
    for informant in Set(ORGANISMS).difference(PCG.organism_set()):
        del( input[informant] )
        logf("# INFORMANT DELETED:", informant )
    # recreate data structures
    ( ORGANISMS, GENECOMBIS,
      GENE_INFORMANT_SET,
      GENE_IDENTIFIER_SET,
      UNIGENE_INFORMANT_SET) = recreate_global_informantdata(PCG,OPTIONS.target)
    # recreate inwpCBGs
    inwpcbgs = PCG2inwpCBGS(PCG)
    # recreate blocks
    blocks = PCG2blocks(PCG)
    # print inwpCBG structure
    if VERBOSE: print_inwpcbgstructure(inwpcbgs,ORGANISMS)


################################################################################
### HARD remove_singleton_blocks_with_novel_orfs
################################################################################

status = remove_singleton_blocks_with_novel_orfs(blocks,PCG,ORGANISMS,OPTIONS,input)

# REDO translate to blocks (if a block was deleted)
if status:
    print "### remove_singleton_blocks_with_novel_orfs"
    # check if an informant is completely deleted from the PCG
    for informant in Set(ORGANISMS).difference(PCG.organism_set()):
        del( input[informant] )
        logf("# INFORMANT DELETED:", informant )
    # recreate data structures
    ( ORGANISMS, GENECOMBIS,
      GENE_INFORMANT_SET,
      GENE_IDENTIFIER_SET,
      UNIGENE_INFORMANT_SET) = recreate_global_informantdata(PCG,OPTIONS.target)
    # recreate inwpCBGs
    inwpcbgs = PCG2inwpCBGS(PCG)
    # recreate blocks
    blocks = PCG2blocks(PCG)
    # print inwpCBG structure
    if VERBOSE: print_inwpcbgstructure(inwpcbgs,ORGANISMS)


################################################################################
### second step: add (potential, strong) SignalP exon alignemnts
################################################################################
status = update_PCG_with_signalpexons(signalpexonseqs,PCG,OPTIONS)

if status:
    # signalp exon PacbpORFS added to PCG; recreate inwpCBGs
    inwpcbgs = PCG2inwpCBGS(PCG)

# logf message
logf("#", stw.lap(), PCG, "SignalP first exons predicted:", status)


################################################################################
### AGAIN:: detect_and_remove_single_nonfirst_inwpcbg
################################################################################
nonfirstPCG = detect_and_remove_single_nonfirst_inwpcbg(inwpcbgs,PCG,
                GENE_IDENTIFIER_SET)

if nonfirstPCG:
    # check if an informant is completely deleted from the PCG
    for informant in Set(ORGANISMS).difference(PCG.organism_set()):
        status = removeinformantfrominput(informant,input)

    # recreate data structures
    ( ORGANISMS, GENECOMBIS,
      GENE_INFORMANT_SET,
      GENE_IDENTIFIER_SET,
      UNIGENE_INFORMANT_SET) = recreate_global_informantdata(PCG,OPTIONS.target)

    # check if an informant is not in ORGANISMS but still in input
    for informant in Set(input.keys()).difference(ORGANISMS):
        status = removeinformantfrominput(informant,input)

    # recreate inwpCBGs
    inwpcbgs = PCG2inwpCBGS(PCG)


# log message for detect_and_remove_single_nonfirst_inwpcbg
logf("#", stw.lap(), "detect_and_remove_single_nonfirst_inwpcbg (II):", nonfirstPCG)


#################################################################################
#### HMM MAXSR search for missing Orfs in **FIRST** inwpCBGs
#################################################################################
HMM_FIRSTCBG_MAXSR_MIN_AA_LENGTH                 = 20
HMM_FIRSTCBG_MAXSR_HMM_MIN_BITSCORE              = -100
HMM_FIRSTCBG_MAXSR_HMM_SPLIT_ON_GAPSIZE          = 2
HMM_FIRSTCBG_MAXSR_PACBPORF_MIN_BITSCORE         = 5
HMM_FIRSTCBG_MAXSR_PACBPORF_MIN_AA_LENGTH        = 10
HMM_FIRSTCBG_MAXSR_MAX_INTRON_NT_LENGTH          = 175


# recreate inwpcbgs
inwpcbgs = PCG2inwpCBGS(PCG)
if VERBOSE: print_inwpcbgstructure(inwpcbgs,ORGANISMS)

logf(stw.lap(),"HMM FIRSTCBG MAXSR searches start")

# get info on (potential) first & final InwpCBG
posFirst,posFinal = get_first_and_final_inwpcbg_pos(inwpcbgs)

is_any_added = False
for pos in range(0,posFirst+1):
    # check for OPTIONS.omithmm -> no HMM profile searches!
    if OPTIONS.omithmm: continue
    inwpCBG = inwpcbgs[pos]

    # get next and next inwpCBG with reasonable coverage of informants
    next,next_for_fastadb = None, None
    if pos < len(inwpcbgs)-1: next = inwpcbgs[pos+1]
    nextpos = pos+1
    while nextpos < len(inwpcbgs)-1:
        if nextpos >= pos+3: break
        if not next_for_fastadb:
            next_for_fastadb = inwpcbgs[nextpos]
        elif inwpcbgs[nextpos].node_count() > (next_for_fastadb.node_count()+1):
            next_for_fastadb = inwpcbgs[nextpos]
        else:
            pass
        nextpos+=1

    # logf message for start HMM FIRST search for this inwpCBG
    logf(inwpCBG,inwpCBG._get_target_node())

    # create MAXSR HMM search profile
    fname_hmmbuild_maxsr, hmmcoords_maxsr = _create_hmm_profile(
                inwpCBG,area="MAXSR",
                nextcbg=next,verbose=VERBOSE)

    # create LEFTSPRDIF search profile
    fname_hmmbuild_leftsprdif, hmmcoords_leftsprdif = _create_hmm_profile(
                inwpCBG,area="LEFTSPRDIF",
                nextcbg=next,verbose=VERBOSE)

    # check if HMM profile(s) were constructed
    if not fname_hmmbuild_maxsr and not fname_hmmbuild_leftsprdif: continue

    # perform checks on maxsr HMM profile & coords
    if hmmcoords_maxsr:
        if min([ end-sta for sta,end in hmmcoords_maxsr.values()]) <\
        HMM_FIRSTCBG_MAXSR_PACBPORF_MIN_AA_LENGTH:
            logf("NO HMM_FIRSTCBG_MAXSR_PACBPORF_MIN_AA_LENGTH:", inwpCBG)
            _file_cleanup([fname_hmmbuild_maxsr])
            fname_hmmbuild_maxsr = None

        # check if OPTIONS.target still in hmmcoords; in exceptional cases,
        # it can have been removed from the HMMprofile
        if not hmmcoords_maxsr.has_key(OPTIONS.target):
            _file_cleanup([fname_hmmbuild_maxsr])
            fname_hmmbuild_maxsr = None

    # perform checks on leftsprdif HMM profile & coords
    if hmmcoords_leftsprdif:
        if min([ end-sta for sta,end in hmmcoords_leftsprdif.values()]) <\
            HMM_FIRSTCBG_MAXSR_PACBPORF_MIN_AA_LENGTH:
            logf("NO HMM_FIRSTCBG_LEFTSPRDIF_PACBPORF_MIN_AA_LENGTH:", inwpCBG)
            _file_cleanup([fname_hmmbuild_leftsprdif])
            fname_hmmbuild_leftsprdif = None

        if not hmmcoords_leftsprdif.has_key(OPTIONS.target):
            _file_cleanup([fname_hmmbuild_leftsprdif])
            fname_hmmbuild_leftsprdif = None

    # check if HMM profile(s) still exist
    if not fname_hmmbuild_maxsr and not fname_hmmbuild_leftsprdif: continue

    for informant in ORGANISMS:
        # continue if target or already seen informant species
        if informant == OPTIONS.target: continue

        # create fasta search database with ORFs from this informant
        fname_fasta_db = _create_hmm_db(informant,input,inwpCBG,None,
                next_for_fastadb,orf_must_have_start=True)

        # check for no elegiable ORFs at all
        if not fname_fasta_db: continue


        if fname_hmmbuild_maxsr and informant not in hmmcoords_maxsr.keys():
            # run hmmsearch on this fasta database with hmmbuild file; maximum 2 hits
            results = hmmsearch_protein( fname_hmmbuild_maxsr, fname_fasta_db, params= {'E':150,'A': 3})
            # logf raw HMM search results
            if results: logf([ (res[0],res[-2]) for res in results ])
            # create pacbporfhmmlist list
            pacbporfhmmlist = hmmresults2splittedpacbps(results,hmmcoords_maxsr,
                    OPTIONS.target,informant,inwpCBG,input,
                    gapsize=HMM_MAXSR_HMM_SPLIT_ON_GAPSIZE,
                    min_bitscore=HMM_FIRSTCBG_MAXSR_HMM_MIN_BITSCORE)

            # logf pacbporf HMM search results
            if VERBOSE and pacbporfhmmlist: logf([ (pf.orfS.id,pf.get_unextended_length(),pf.bitscore,pf.has_upstream_methionines()) for pf in pacbporfhmmlist ])

        elif fname_hmmbuild_leftsprdif and informant not in hmmcoords_leftsprdif.keys():
            # run hmmsearch on this fasta database with hmmbuild file; maximum 2 hits
            results = hmmsearch_protein( fname_hmmbuild_leftsprdif, fname_fasta_db, params= {'E':150,'A': 3})
            # logf raw HMM search results
            if results: logf([ (res[0],res[-2]) for res in results ])
            # create pacbporfhmmlist list
            pacbporfhmmlist = hmmresults2splittedpacbps(results,hmmcoords_leftsprdif,
                    OPTIONS.target,informant,inwpCBG,input,
                    gapsize=HMM_MAXSR_HMM_SPLIT_ON_GAPSIZE,
                    min_bitscore=HMM_FIRSTCBG_MAXSR_HMM_MIN_BITSCORE)

            # logf pacbporf HMM search results
            if VERBOSE and pacbporfhmmlist: logf([ (pf.orfS.id,pf.get_unextended_length(),pf.bitscore,pf.has_upstream_methionines()) for pf in pacbporfhmmlist ])

        else:
            # file cleanup of HMM fasta database
            _file_cleanup([fname_fasta_db])
            continue

        # file cleanup of HMM fasta database
        _file_cleanup([fname_fasta_db])

        # check/add the pacbporfs to the PCG
        thepacbporfs = pacb.ordering.order_pacbporf_list(PCG.get_pacbps_by_organisms(OPTIONS.target,informant))

        for hmmpacbporf in pacbporfhmmlist:
            # filter out to short ones
            if hmmpacbporf.get_unextended_length() < HMM_FIRSTCBG_MAXSR_PACBPORF_MIN_AA_LENGTH:
                continue

            # filter out poor bitscores
            if hmmpacbporf.bitscore < HMM_FIRSTCBG_MAXSR_PACBPORF_MIN_BITSCORE:
                continue

            # filter out to far remote alignments
            endQprev = hmmpacbporf._get_original_alignment_pos_end().query_dna_end
            staQthis = thepacbporfs[0]._get_original_alignment_pos_start().query_dna_start
            endSprev = hmmpacbporf._get_original_alignment_pos_end().sbjct_dna_end
            staSthis = thepacbporfs[0]._get_original_alignment_pos_start().sbjct_dna_start
            if (staQthis - endQprev) > HMM_FIRSTCBG_MAXSR_MAX_INTRON_NT_LENGTH:
                continue
            if (staSthis - endSprev) > HMM_FIRSTCBG_MAXSR_MAX_INTRON_NT_LENGTH:
                continue


            # filter out pacbporfs without has_upstream_methionines()
            if not hmmpacbporf.has_upstream_methionines():
                continue

            ####################################################################
            if VERBOSE:
                logf(informant, hmmpacbporf)
                hmmpacbporf.print_protein(_linesize=100)
            ####################################################################

            # calculate overlaps with current listed pacbporfs
            overlaps = [ pf.overlap(hmmpacbporf) for pf in thepacbporfs ]

            if overlaps[0] == 1.0:
                # already fully included in the first PacbPORF
                continue

            if overlaps[0] > 0.0 and hmmpacbporf.relatively_positioned_towards(thepacbporfs[0])['S1'][2] > 0:
                # Sbjct HMM position is located AFTER the first alignment -> no colinear alignments
                continue

            if len(overlaps) - overlaps.count(0.0) >= 2:
                # hmmpacbporf is overlapping with >=2 PacbPORFS -> ignore
                continue

            if overlaps.index(max(overlaps)) != 0:
                # hmmpacbporf is not overlapping with the first but another
                # pacbporf -> no extra alignment at 5'side of protein
                continue

            rejected = is_hmmpacbporf_conflicting_with_pacbporflist(
                        hmmpacbporf, thepacbporfs)

            if not rejected or informant in inwpCBG.organism_set():
                # if here, add Node & PacbPORF to the PCG.
                hmmpacbporf2PCG(hmmpacbporf,OPTIONS.target,informant,PCG)   
                is_any_added = True
                logf("# HMM pacbporf FIRST MAXSR ADDED::", informant, hmmpacbporf)

                if rejected and informant in inwpCBG.organism_set():
                    # correct overlap in PacbPORFs for this informant species only
                    status = PCG.correct_overlaps([(OPTIONS.target,informant)],verbose=False)
                    # recreate extended pacbporfs from pacbps;
                    # because PCG.correct_overlaps() has called PCG.pacbporfs2pacbps()
                    PCG.extend_pacbporfs(input)

                # re-get list of pacbporfs
                check_cnt = len(thepacbporfs)
                thepacbporfs = pacb.ordering.order_pacbporf_list(
                        PCG.get_pacbps_by_organisms(OPTIONS.target,informant)
                        )
                print "NEW PACBPORFS::", len(thepacbporfs), "was:", check_cnt

        # logf message for this HMM-searched inwpCBG
        logf(inwpCBG, "FIRSTCBG MAXSR", informant, len(results), len(pacbporfhmmlist))

    # file cleanup of HMM profile(s)
    _file_cleanup([fname_hmmbuild_maxsr])
    _file_cleanup([fname_hmmbuild_leftsprdif])

# HMM log message
logf("#", stw.lap(),"HMM FIRSTCBG MAXSR searches done")

if is_any_added and PCG_HAS_HIGH_GAP_RATIO_PACBPORFS:
    # perform again PCG.merge_high_gap_ratio_pacbporfs() after HMM searches
    PCG_HAS_HIGH_GAP_RATIO_PACBPORFS = PCG.merge_high_gap_ratio_pacbporfs(OPTIONS.target,verbose=True)
    logf(stw.lap(),PCG,"PCG.merge_high_gap_ratio_pacbporfs()",PCG_HAS_HIGH_GAP_RATIO_PACBPORFS)

# recreate inwpcbgs in case any pacbporf was added
if is_any_added: inwpcbgs = PCG2inwpCBGS(PCG)

if VERBOSE:
    logf(PCG,"after HMM FIRSTCBG MAXSR searches")
    print_inwpcbgstructure(inwpcbgs,ORGANISMS)


#################################################################################
#### HMM MAXSR search for missing Orfs in **FINAL** inwpCBGs
#################################################################################
HMM_FINALCBG_MAXSR_MIN_AA_LENGTH                 = 15
HMM_FINALCBG_MAXSR_HMM_MIN_BITSCORE              = -50
HMM_FINALCBG_MAXSR_HMM_SPLIT_ON_GAPSIZE          = 2
HMM_FINALCBG_MAXSR_PACBPORF_MIN_BITSCORE         = 5
HMM_FINALCBG_MAXSR_PACBPORF_MIN_AA_LENGTH        = 6

logf(stw.lap(),"HMM FINALCBG MAXSR searches start")

# get info on (potential) first & final InwpCBG
posFirst,posFinal = get_first_and_final_inwpcbg_pos(inwpcbgs)

is_any_added = False
for pos in range(posFinal,len(inwpcbgs)):
    # check for OPTIONS.omithmm -> no HMM profile searches!
    if OPTIONS.omithmm: continue
    inwpCBG = inwpcbgs[pos]

    if not GENE_IDENTIFIER_SET.difference(inwpCBG.organism_set()):
        # no organisms/genes/nodes missing -> continue
        continue

    # get prev and prev inwpCBG with reasonable coverage of informants
    prev,prev_for_fastadb = None, None
    if pos >= 1: prev = inwpcbgs[pos-1]
    prevpos = pos-1
    while prevpos >= 0:
        if prevpos <= pos-3: break
        if not prev_for_fastadb:
            prev_for_fastadb = inwpcbgs[prevpos]
        elif inwpcbgs[prevpos].node_count() > (prev_for_fastadb.node_count()+1):
            prev_for_fastadb = inwpcbgs[prevpos]
        else:
            pass
        prevpos-=1


    # logf message for start HMM FINALsearch for this inwpCBG
    logf(inwpCBG,inwpCBG._get_target_node())

    # create MAXSR HMM search profile
    fname_hmmbuild, hmmcoords = _create_hmm_profile(
                inwpCBG,area="MAXSR",
                prevcbg=prev,verbose=VERBOSE)

    # no MAXSR profile can be constructed
    if not fname_hmmbuild: continue
    if min([ end-sta for sta,end in hmmcoords.values()]) < HMM_FINALCBG_MAXSR_PACBPORF_MIN_AA_LENGTH:
        logf("NO HMM_FINALCBG_MAXSR_PACBPORF_MIN_AA_LENGTH:", inwpCBG)
        _file_cleanup([fname_hmmbuild])
        continue

    # check if OPTIONS.target still in hmmcoords; in exceptional cases,
    # it can have been removed from the HMMprofile
    if not hmmcoords.has_key(OPTIONS.target):
        _file_cleanup([fname_hmmbuild])
        continue

    for informant in ORGANISMS:
        # continue if target or already seen informant species
        if informant == OPTIONS.target: continue

        # do NOT omit informants already listed in the profile;
        # try to improve/elongate their final alignments too
        #if informant in inwpCBG.organism_set(): continue

        # create fasta search database with ORFs from this informant
        fname_fasta_db = _create_hmm_db(informant,input,inwpCBG,
                prev_for_fastadb,None)

        # check for no elegiable ORFs at all
        if not fname_fasta_db: continue

        # run hmmsearch on this fasta database with hmmbuild file; maximum 2 hits
        results = hmmsearch_protein( fname_hmmbuild, fname_fasta_db, params= {'E':150,'A': 3})
        if VERBOSE and results: logf([ (res[0],res[-2]) for res in results ])

        # create pacbporfhmmlist list
        pacbporfhmmlist = hmmresults2splittedpacbps(results,hmmcoords,
                OPTIONS.target,informant,inwpCBG,input,
                gapsize=HMM_FINALCBG_MAXSR_HMM_SPLIT_ON_GAPSIZE,
                min_bitscore=HMM_FINALCBG_MAXSR_HMM_MIN_BITSCORE)

        if VERBOSE and pacbporfhmmlist: logf([ (pf.orfS.id,pf.get_unextended_length(),pf.bitscore,pf.get_projected_tailing_stop_aa_difference()) for pf in pacbporfhmmlist ])

        # check/add the pacbporfs to the PCG
        thepacbporfs = pacb.ordering.order_pacbporf_list(PCG.get_pacbps_by_organisms(OPTIONS.target,informant))

        # file cleanup of HMM fasta database
        _file_cleanup([fname_fasta_db])

        # check if thepacbporfs exists; in freak cases, list can be empty
        # (informant removed, but not consistently)
        if not thepacbporfs: continue

        for hmmpacbporf in pacbporfhmmlist:
            # filter out to short ones
            if hmmpacbporf.get_unextended_length() < HMM_FINALCBG_MAXSR_PACBPORF_MIN_AA_LENGTH:
                continue

            # filter out poor bitscores
            if hmmpacbporf.bitscore < HMM_FINALCBG_MAXSR_PACBPORF_MIN_BITSCORE:
                continue

            if thepacbporfs[-1].overlap(hmmpacbporf) == 1.0:
                # already fully included in the final PacbPORF
                continue

            is_overlapping = False
            for pF in thepacbporfs:
                if pF.overlap(hmmpacbporf) >= 0.85 and\
                (pF.bitscore > hmmpacbporf.bitscore or pF.identityscore > hmmpacbporf.identityscore):
                    is_overlapping = True
                    break
            if is_overlapping:
                # already (fully) included in a better PacbPORF
                continue
                
            # filter out poor projected_tailing_stop_aa_difference
            # in comparison with current (final) pacbporf
            if informant in inwpCBG.organism_set():
                current_pacbporf = inwpCBG.get_pacbps_by_organisms(OPTIONS.target,informant)[0]
                if hmmpacbporf.get_projected_tailing_stop_aa_difference() >=\
                current_pacbporf.get_projected_tailing_stop_aa_difference():
                    continue

            # filter out poor projected_tailing_stop_aa_difference
            if hmmpacbporf.get_projected_tailing_stop_aa_difference() > 5:
                continue
 
            rejected = is_hmmpacbporf_conflicting_with_pacbporflist(
                        hmmpacbporf, thepacbporfs)

            if not rejected or informant in inwpCBG.organism_set():
                # if here, add Node & PacbPORF to the PCG.
                hmmpacbporf2PCG(hmmpacbporf,OPTIONS.target,informant,PCG)   
                is_any_added = True
                logf("# HMM pacbporf FINAL MAXSR ADDED::", informant, hmmpacbporf, rejected, informant in inwpCBG.organism_set())
                if rejected:
                    for pF in thepacbporfs:
                        print informant, pF
                
                if rejected and informant in inwpCBG.organism_set():
                    # correct overlap in PacbPORFs for this informant species only
                    status = PCG.correct_overlaps([(OPTIONS.target,informant)],verbose=False)
                    # recreate extended pacbporfs from pacbps;
                    # because PCG.correct_overlaps() has called PCG.pacbporfs2pacbps()
                    PCG.extend_pacbporfs(input)

                # re-get list of pacbporfs
                thepacbporfs = pacb.ordering.order_pacbporf_list(
                        PCG.get_pacbps_by_organisms(OPTIONS.target,informant)
                        )

    # file cleanup of HMM profile
    _file_cleanup([fname_hmmbuild])

# HMM log message
logf("#", stw.lap(),"HMM FINALCBG MAXSR searches done")

# recreate inwpcbgs in case any pacbporf was added
if is_any_added: inwpcbgs = PCG2inwpCBGS(PCG)

if VERBOSE:
    logf(PCG,"after HMM FINALCBG MAXSR searches")
    print_inwpcbgstructure(inwpcbgs,ORGANISMS)


#################################################################################
#### Find missing FinalExon of target species
#################################################################################

# get info on (potential) first & final InwpCBG
posFirst,posFinal = get_first_and_final_inwpcbg_pos(inwpcbgs)

is_any_added = False
for pos in range(posFinal,len(inwpcbgs)):

    # check for OPTIONS.omithmm -> no HMM profile searches!
    if OPTIONS.omithmm: continue
    inwpCBG = inwpcbgs[pos]
    prev,next = None,None
    if pos > 0: prev = inwpcbgs[pos-1]
    if pos < len(inwpcbgs)-1: next = inwpcbgs[pos+1]

    logf(inwpCBG,inwpCBG._get_target_node())

    for informant in GENE_INFORMANT_SET.intersection(inwpCBG.organism_set()):
        pacbporf = inwpCBG.get_pacbps_by_organisms(OPTIONS.target,informant)[0]
        endpos = pacbporf._get_original_alignment_pos_end()
        qstop  = pacbporf.orfQ.protein_endPY - endpos.query_pos
        sstop  = pacbporf.orfS.protein_endPY - endpos.sbjct_pos
        print informant, pacbporf,sstop,(qstop,)
        print pacbporf.get_projected_tailing_stop_aa_difference(),
        print pacbporf.get_projected_tailing_stop_nonaligned_aa_difference()

    # create RIGTHORFEND HMM search profile
    fname_hmmbuild, hmmcoords = _create_hmm_profile(
                inwpCBG,area="RIGTHORFEND",
                prevcbg=prev,verbose=True)

    # no RIGTHORFEND profile can be constructed
    if not fname_hmmbuild: continue

    # create fasta search database with ORFs from this informant
    fname_fasta_db = _create_hmm_db(OPTIONS.target,input,inwpCBG,prev,next)

    # check for no elegiable ORFs at all
    if not fname_fasta_db:
        _file_cleanup([fname_hmmbuild])
        continue

    # run hmmsearch on this fasta database with hmmbuild file; maximum 2 hits
    results = hmmsearch_protein( fname_hmmbuild, fname_fasta_db, params= {'E':150,'A': 3})
    if results: print [ (res[0],res[-2]) for res in results ]

    _file_cleanup([fname_hmmbuild])



################################################################################
# InwpCbgRemoval done. MileStone logging
################################################################################
if not VERBOSE: stdoutManager.reset()
logf("##", stwMS.lap(), "InwpCbgRemoval: performed")

if not True in (QUIET,SILENT):
    print PCG
    print_inwpcbgstructure(inwpcbgs,ORGANISMS)
    
################################################################################
### starting next algorithm step: intron predictions
################################################################################
if not VERBOSE: stdoutManager.buffer() # capture all STDOUT messages from deeper modules 
stwMS.stepname = "IntronPrediction"
stw.stepname = "IntronPrediction"


################################################################################
### ``REDO`` all the donor/acceptor site scanning.
### dunno why & where, but when sites are scanned a priori without
### allow non-canonical sites, and because the site scanning is cached, 
### non -canonical sites are not taken into account anymore.
### So, here HARD-RESET all splice site scanning with default settings,
### allowing for non-canonical sites
################################################################################

nodes_done = []
for (pacbpkey,nodeQ,nodeS), pacbporf in PCG.pacbps.iteritems():
    if not nodeQ in nodes_done:
        pacbporf.orfQ.scan_orf_for_pssm_donor_sites(forced=True)
        pacbporf.orfQ.scan_orf_for_pssm_acceptor_sites(forced=True)
        nodes_done.append(nodeQ)
    if not nodeS in nodes_done:
        pacbporf.orfS.scan_orf_for_pssm_donor_sites(forced=True)
        pacbporf.orfS.scan_orf_for_pssm_acceptor_sites(forced=True)
        nodes_done.append(nodeS)


logf("#", stw.lap(),"Donor & Acceptor sites (re)predicted")


################################################################################
### perform intron mapping/projecting in between the PacbPORFs
################################################################################

observed_introns = dict([ [ org, {} ] for org in ORGANISMS ])
assessed_interfaces = dict([ [ org, {} ] for org in ORGANISMS ])
intron2label = dict([ [ org, {} ] for org in ORGANISMS ])

for orgA,orgB in GENECOMBIS:

    if not orgA in PCG.organism_set(): continue
    if not orgB in PCG.organism_set(): continue
    pacbporfs = pacb.ordering.order_pacbporf_list(PCG.get_pacbps_by_organisms(orgA,orgB))

    print "COMBI:", orgA,orgB
    if len(pacbporfs) <= 1: continue

    for pos in range(1,len(pacbporfs)):
        prevPACBP = pacbporfs[pos-1]
        nextPACBP = pacbporfs[pos]

        overlapA = nextPACBP.overlap(prevPACBP)
        if pos < len(pacbporfs)-1:
            secondnextPACBP = pacbporfs[pos+1]
            overlapB = nextPACBP.overlap(secondnextPACBP)
        else:
            secondnextPACBP = None
            overlapB = None 

        # get OrfSet objects for query & sbjct
        qOrfs = input[orgA]['orfs']
        sOrfs = input[orgB]['orfs']

        # check if informant isa unigene -> then allow_sbjct_mapping == False
        if input[orgB].has_key('is_unigene') and input[orgB]['is_unigene']:
            allow_sbjct_mapping = False
        else:
            allow_sbjct_mapping = True

        # try to merge the pacbpsorf with an intron
        introns = pacb.connecting.connecting.merge_pacbporfs(
                prevPACBP,nextPACBP,
                qOrfs,sOrfs,
                allow_sbjct_mapping=allow_sbjct_mapping
                )

        if not len(introns['query']) and prevPACBP.orfQ.id != nextPACBP.orfQ.id\
        and ((overlapA and overlapA > 0.10) or (overlapB and overlapB > 0.10))\
        and secondnextPACBP:
            # there *SHOULD* be introns predicted for the query
            # and, nextPACBP is (largely) overlapping with its
            # neighbours. Maybe it is a non-sense alignment!?
            print "SECONDNEXT -> PacbPORF"
            nextPACBP = secondnextPACBP
            introns = pacb.connecting.connecting.merge_pacbporfs(prevPACBP,nextPACBP,qOrfs,sOrfs)


        # label tags the interface of PacbPORFs for the query
        label = (prevPACBP.orfQ.id,nextPACBP.orfQ.id)
        informant_label = (prevPACBP.orfS.id,nextPACBP.orfS.id)
        _update_observed_introns(introns['query'],label,informant_label,orgA,orgB,observed_introns,intron2label,assessed_interfaces)

        # label tags the interface of PacbPORFs for the sbjct
        label = (prevPACBP.orfS.id,nextPACBP.orfS.id)
        informant_label = (prevPACBP.orfQ.id,nextPACBP.orfQ.id)
        _update_observed_introns(introns['sbjct'],label,informant_label,orgB,orgA,observed_introns,intron2label,assessed_interfaces)


# logf message
logf("#", stw.lap(), "ABGP intron mapping & projecting done", PCG)

################################################################################
# IntronPrediction done. MileStone logging
################################################################################
if not VERBOSE: stdoutManager.reset()
logf("##", stwMS.lap(), "IntronPrediction: performed")

################################################################################
### starting next algorithm step: SequenceError and/or Disruptive mutation detection
################################################################################
if not VERBOSE: stdoutManager.buffer() # capture all STDOUT messages from deeper modules 
stwMS.stepname = "SequenceErrorsOrDisruptiveMutations"
stw.stepname = "SequenceErrorsOrDisruptiveMutations"


############################################################################
### apply insertion/deletion SequenceError coordinate corrections
### inserstion/deletion coordinates are calculated based on Orf overlap
### and *NOT* on the actual alignments (in order to have identical positions),
### Here, this is corrected for based on the PacbPORF alignments
############################################################################
for org,data in observed_introns.iteritems():
    if org != OPTIONS.target: continue
    if not data: continue
    replaced_seqerror_keys = []
    corrected_seqerror_data = {}
    for key,vlist in data.iteritems():
        if vlist[0].__class__.__name__ == "SequenceErrorConnectingOrfs" and\
        vlist[0].typeof in ['insertion','deletion']:
            seobjlist = apply_indel_position_correction(vlist)
            # create new key in introndata, intron2label, ....
            new_key = (seobjlist[0].start,seobjlist[0].end)
            corrected_seqerror_data[new_key] = seobjlist
            intron2label[OPTIONS.target][new_key] = intron2label[OPTIONS.target][key]
            # append new key to the list of replaced_seqerror_keys
            replaced_seqerror_keys.append( key )
            logf("#", "SequenceError coordinate correction:", key, "->", new_key)
           
    # update the corrected intron/seqerror keys
    data.update(corrected_seqerror_data)

    # remove all the replaced_seqerror_keys
    for old_key in replaced_seqerror_keys:
        del( data[old_key] )
        del( intron2label[OPTIONS.target][old_key] )
        logf("SequenceError old coord DELETION:", old_key)

# logf message EOF step
logf("#", stw.lap(), "SequenceError Coordinate correction done", PCG)


############################################################################
### Check all the predicted potential SequenceError interfaces for validity
############################################################################
GENE_MODEL_HAS_SEQUENCE_ERRORS = check_and_cleanup_sequenceerror_interfaces(
    observed_introns,inwpcbgs,PCG,OPTIONS,input,verbose=True)
       
# logf message EOF step
logf("#", stw.lap(), "SequenceError Assesment done; seqerrors present:",
    GENE_MODEL_HAS_SEQUENCE_ERRORS)

############################################################################
### If potential SequenceError interfaces: reasses all intron interfaces
### for hidden (in less conserved regions) SequenceErrors
############################################################################
if True or GENE_MODEL_HAS_SEQUENCE_ERRORS:
    GENE_MODEL_HAS_HIDDEN_SEQUENCE_ERRORS = check_intron_interfaces_for_hidden_sequenceerrors(
        observed_introns,inwpcbgs,PCG,OPTIONS,input,verbose=True)

    # cleanup intron2label dict after seqerror assessment & novel seqerror discovery
    intron_keys = intron2label[OPTIONS.target].keys()
    for key in intron_keys:
        if not observed_introns[OPTIONS.target].has_key(key):
            # intron is rejected because of SequenceError call
            del( intron2label[OPTIONS.target][key] )

    # update assessed_interfaces[OPTIONS.target] dict
    # TODO: here only the succesfull SeqErrors are added to
    # interface assesment, not the failed ones!
    for key in observed_introns[OPTIONS.target].keys():
        ilist =  observed_introns[OPTIONS.target][key]
        label = _get_main_interface( [ i._label for i in ilist ] )
        if not assessed_interfaces[OPTIONS.target].has_key(label):
            assessed_interfaces[OPTIONS.target][label] = [ i._reference for i in ilist ]
        if not intron2label[OPTIONS.target].has_key(key):
            intron2label[OPTIONS.target][key] = label
else:
    GENE_MODEL_HAS_HIDDEN_SEQUENCE_ERRORS = None

    
# logf message EOF step
logf("#", stw.lap(), "SequenceError Intron Assesment done; hidden seqerrors:",
    GENE_MODEL_HAS_HIDDEN_SEQUENCE_ERRORS)

    
################################################################################
# SequenceErrorsOrDisruptiveMutations. MileStone logging
################################################################################
if not VERBOSE: stdoutManager.reset()
logf("##", stwMS.lap(), "SequenceErrorsOrDisruptiveMutations: performed")


################################################################################
### starting next algorithm step: SmallTSSExon detection 
################################################################################
if not VERBOSE: stdoutManager.buffer() # capture all STDOUT messages from deeper modules 
stwMS.stepname = "SmallTSSExonDetection"
stw.stepname = "SmallTSSExonDetection"
    
    
################################################################################
# The next code works only in case OPTIONS.abinitio == False
# It detects small frontal exons in the provided gene structures
# and it it finds those, ABFGP tries harder to find small frontal TSS-exons
################################################################################
if OPTIONS.abinitio == False:
    small_first_exon_annotated = False
    for org in ORGANISMS:
        for warning in input[org]['warnings']:
            if warning.__class__.__name__ == 'SmallAnnotatedFirstExonWarning':
                small_first_exon_annotated = True
                break

    small_first_exons_added = 0
    if small_first_exon_annotated and input[OPTIONS.target]['orfid-genestructure']:
        firstexons = []
        for warning in input[OPTIONS.target]['warnings']:
            if warning.__class__.__name__ == 'SmallAnnotatedFirstExonWarning':
                firstorf   = input[OPTIONS.target]['orfs'].get_orf_by_id(input[OPTIONS.target]['orfid-genestructure'][0])
                firstexons = get_potention_first_exons_on_orf(firstorf,max_tinyexon_nt_length=75)
                break
        # check if there are firstexons
        if not firstexons:
            # the target species has not an annotated tiny first exon we can make use of....
            # so: try to harvest them from ... nearby Orfs
            informant_exon_objects = []
            for org in ORGANISMS:
                if org == OPTIONS.target: continue
                if not input[org]['orfid-genestructure']: continue
                for warning in input[org]['warnings']:
                    if warning.__class__.__name__ == 'SmallAnnotatedFirstExonWarning':
                        exonobjects = input[org]['gldobj'].as_exons(input=input[org])
                        informant_exon_objects.append( exonobjects[0] )
                        logf("#", org, exonobjects[0])
                        break
            # get some data on these tiny exons
            if informant_exon_objects:
                min_exon_nt_length = min([ te.length for te in informant_exon_objects ] ) - 3
                max_exon_nt_length = max([ te.length for te in informant_exon_objects ] ) + 3
                max_orf_start = min(inwpcbgs[0].minimal_spanning_range(organism=OPTIONS.target))*3
                min_orf_end   = max_orf_start - 150
                for orfX in input[org]['orfs'].get_elegiable_orfs(max_orf_start=max_orf_start,min_orf_end=min_orf_end):
                    # predict potential tinyexons
                    firstexons.extend( get_potention_first_exons_on_orf(orfX,
                        max_tinyexon_nt_length      = max_exon_nt_length,
                        min_tinyexon_nt_length      = min_exon_nt_length ) )
                    logf("#", "potential first orf:", orfX, len(firstexons) )

                # order firstexons on summed pssm score
                tmp = [ (exon.acceptor.pssm_score+exon.donor.pssm_score,exon) for exon in firstexons ]
                tmp.sort()
                tmp.reverse()
                firstexons = [ exon for score,exon in tmp ]
                # limit to best 5 scoring ones and hope for the best ....
                firstexons = firstexons[0:5]
                for fe in firstexons: logf("#", fe, fe.proteinsequence() )
     
        else:
            # get only a subset of the best X scoring firstexon
            firstexons = firstexons[0:3]
            for exon in firstexons:
                if (OPTIONS.target,exon.orf.id) not in PCG.get_organism_nodes(OPTIONS.target):
                    exon.orf._has_donor_sites_predicted = True
                    exon.orf._donor_sites = [ exon.donor ]
                    exon.orf._has_tss_predicted = True
                    exon.orf._tss_sites = [ exon.acceptor ]
                else:
                    # node already in PCG do *NOT* limit on donor sites,
                    # this can give grove errors in case firstexons[0] is not
                    # the correct one
                    pass
        
        for org in ORGANISMS:
            if org == OPTIONS.target: continue
            if not firstexons: continue
            if not input[org]['orfid-genestructure']: continue

            # get first orf with annotated tiny first exon
            firstorf = input[org]['orfs'].get_orf_by_id(input[org]['orfid-genestructure'][0])
            ###print org
            TSS_EXON_LINKED_WITH_INTRON = False
            for leadingQexon in firstexons:
                # check if TSS_EXON_LINKED_WITH_INTRON is True -> break
                if TSS_EXON_LINKED_WITH_INTRON == True:
                    break
                for leadingSexon in get_potention_first_exons_on_orf(firstorf,max_tinyexon_nt_length=100):
                    # check if TSS_EXON_LINKED_WITH_INTRON is True -> break
                    if TSS_EXON_LINKED_WITH_INTRON == True:
                        break

                    leadingSexon.orf = deepcopy(leadingSexon.orf)
                    leadingSexon.orf._has_donor_sites_predicted = True
                    leadingSexon.orf._donor_sites = [ leadingSexon.donor ]
                    leadingSexon.orf._has_tss_predicted = True
                    leadingSexon.orf._tss_sites = [ leadingSexon.acceptor ]

                    # align the sequences with clustalw
                    headerQ = str(OPTIONS.target) 
                    headerS = str(org)
                    protseqQ= leadingQexon.proteinsequence()
                    protseqS= leadingSexon.proteinsequence()
                    coords  = [ leadingQexon.protein_start(),
                                leadingQexon.protein_end(),
                                leadingSexon.protein_start(),
                                leadingSexon.protein_end(), ]
                    if leadingQexon.length == 2:
                        protseqQ = "M"
                        coords[1] = coords[0]+1
                    if leadingSexon.length == 2:
                        protseqS = "M"
                        coords[3] = coords[2]+1

                    (alignedseqs,alignment) =\
                    clustalw( seqs= { headerQ: protseqQ, headerS: protseqS } )
            
                    # make pacbp from clustalw alignment
                    pacbp = pacbp_from_clustalw(
                                alignment=(
                                        alignedseqs[headerQ],
                                        alignment,
                                        alignedseqs[headerS]
                                        ),
                                coords=coords
                                )
            

                    logf( org, "\t", leadingQexon.orf.id, leadingSexon.orf.id, "\t", "smallTSS pacbp:", pacbp)
                    
                    if not pacbp: continue
                    # strip unaligned fraction of this pacbp object, then check length
                    pacbp.strip_unmatched_ends()
                    if pacbp.bits <= 0: continue
                    if pacbp.identityscore < 0.40: continue

                    tsspacbporf = pacbp2pacbporf(pacbp,leadingQexon.orf,leadingSexon.orf)
                    tsspacbporf.extend_pacbporf_after_stops()

                    # get ordered list of pacbporf objects of target - informant
                    ordered_pacbporfs = pacb.ordering.order_pacbporf_list(
                        PCG.get_pacbps_by_organisms(OPTIONS.target,org))

                    # check if tsspacbporf is fully included in an already listed pacbporf
                    if max([ pf.overlap(tsspacbporf) for pf in ordered_pacbporfs ]) >= 1.0:
                        logf("tssexon included:",org,tsspacbporf)
                        continue

                    for pacbporf in ordered_pacbporfs[0:1]:
                        tssintrons = {'query':[],'sbjct':[]}

                        # try to merge the pacbpsorf with a convincing intron
                        introns = _filter_aligned_introns_on_pssm_entropy_combination(
                                pacb.connecting.connecting.merge_pacbporfs_with_introns(
                                    tsspacbporf,pacbporf,
                                    max_aa_offset=1,
                                    max_intron_nt_length=150,
                                ) )

                        logf("#", org, "\t", leadingQexon.orf.id, leadingSexon.orf.id, "\t", tsspacbporf, " introns::", len(introns) )

                        if introns:
                            max_apps_accep = max([ intronQ._apps_accep for intronQ,intronS in introns ])
                            for intronQ,intronS in introns:
                                # ignore poor positional alignment scores
                                if intronQ._apps_accep < 0.5*max_apps_accep: continue

                                # store tsspacbporf to the intron object for later addition to the PCG/inwpCBGS
                                intronQ._linked_to_pacbporfs = [ tsspacbporf ]
                                intronS._linked_to_pacbporfs = [ tsspacbporf ]
                                
                                # update tssintrons result dict
                                if intronQ.coords() not in [ intr.coords() for intr in tssintrons['query'] ]:
                                    tssintrons['query'].append(intronQ)
                                if intronS.coords() not in [ intr.coords() for intr in tssintrons['sbjct'] ]:
                                    tssintrons['sbjct'].append(intronS)

                        # try a more complex mapping of introns
                        introns = pacb.connecting.mapping.merge_pacbporfs_with_closeby_independant_introns(
                                    tsspacbporf,pacbporf,verbose=False,cig_max_aa_length=20
                                    )
                        for intronQ,intronS,cigpacbporf in introns:
                            intronQ._linked_to_pacbporfs.append(tsspacbporf)
                            intronS._linked_to_pacbporfs.append(tsspacbporf)
                            if intronQ.coords() not in [ intr.coords() for intr in tssintrons['query'] ]:
                                tssintrons['query'].append(intronQ)
                            if intronS.coords() not in [ intr.coords() for intr in tssintrons['sbjct'] ]:
                                tssintrons['sbjct'].append(intronS)

                        # label tags the interface of PacbPORFs for the query
                        label = (tsspacbporf.orfQ.id,pacbporf.orfQ.id)
                        informant_label = (tsspacbporf.orfS.id,pacbporf.orfS.id)
                        _update_observed_introns(tssintrons['query'],label,informant_label,OPTIONS.target,org,observed_introns,intron2label,assessed_interfaces)
                
                        # label tags the interface of PacbPORFs for the sbjct
                        label = (tsspacbporf.orfS.id,pacbporf.orfS.id)
                        informant_label = (tsspacbporf.orfQ.id,pacbporf.orfQ.id)
                        _update_observed_introns(tssintrons['sbjct'],label,informant_label,org,OPTIONS.target,observed_introns,intron2label,assessed_interfaces)

                        if tssintrons['query']:
                            # introns have been found! assume they are correct
                            # and do not allow next pacbporfs to contribute too
                            for intronQ in tssintrons['query']:
                                logf("#","SMALL TSS EXON predicted:", org, intronQ._linked_to_pacbporfs[0])
                            TSS_EXON_LINKED_WITH_INTRON = True
                            small_first_exons_added+=1
                            break

    # logf message
    logf("#", stw.lap(), small_first_exons_added,"small TSS exon(s) & introns predicted",PCG)


################################################################################
# SmallTSSExonDetection. MileStone logging
################################################################################
if not VERBOSE: stdoutManager.reset()
logf("##", stwMS.lap(), "SmallTSSExonDetection: performed")

################################################################################
### starting next algorithm step: stopless 3n intron polishing 
################################################################################
if not VERBOSE: stdoutManager.buffer() # capture all STDOUT messages from deeper modules 
stwMS.stepname = "Stopless3nIntronPolishing"
stw.stepname = "Stopless3nIntronPolishing"    
    
################################################################################
### detect (falsely) aligned stopless3n introns IN the PacbPORFs
################################################################################
DETECT_STOPLESS_3N_INTRONS_MIN_NT_INGAPPED_IDENTITY = 0.50
pacbp_keys_done = []
for inwpCBG in inwpcbgs:
    grScore = max([ pf.gap_ratio_score() for pf in inwpCBG.pacbps.values() ]) 
    if grScore >= PACBPORF_HIGH_GAP_RATIO_THRESHOLD:
        print "omit stopless3n:", inwpCBG
        continue

    # if here, search the pacbporfs of this inwpCBG for stopless3n introns
    for (k,n1,n2),pacbporf in inwpCBG.pacbps.iteritems():
        # check if this pacbporf is checked already
        if (k,n1,n2) in pacbp_keys_done:
            continue
        else:
            pacbp_keys_done.append( (k,n1,n2) )

        # omit pacbporfs which are (very) unlikely to contain aligned
        # stopless3n introns. Criteria are: presence of gaps & quit high
        # get_ungapped_nt_identity()
        if pacbporf.alignment_has_gaps() and pacbporf.get_ungapped_nt_identity() >\
        DETECT_STOPLESS_3N_INTRONS_MIN_NT_INGAPPED_IDENTITY:
            # find aligned stopless3n introns in the QUERY only
            stopless3nsQ = pacbporf.detect_aligned_stopless3n_intron_in_query()
        
            # set GFF fsource attribute for recognition of intron sources
            for intron in stopless3nsQ:
                intron._distance = 0  # !?!?!?!?!
                intron._linked_to_introns = []
                intron._linked_to_pacbporfs = []
                intron._gff['fsource'] = "ABGPstl3n"
                intron._apps_donor = 1.0
                intron._apps_accep = 1.0
        
            # label tags the interface of PacbPORFs for the stopless3n introns in query
            target_label    = (pacbporf.orfQ.id,pacbporf.orfQ.id)
            informant_label = (pacbporf.orfS.id,pacbporf.orfS.id)
            target          = OPTIONS.target
            informant       = n2[0]
            _update_observed_introns(stopless3nsQ,target_label,informant_label,
                    target,informant,observed_introns,intron2label,assessed_interfaces)

    
        # find aligned CONSERVED stopless3n introns
        conservedstopless = pacbporf.detect_conserved_stopless3n_intron(verbose=True)
    
        # set GFF fsource attribute for recognition of intron sources
        for intronQ,intronS in conservedstopless:
            # intron._distance = 0 # TAKEN CARE OF IN THE FUNCTION
            intronQ._linked_to_introns = []
            intronQ._linked_to_pacbporfs = []
            intronQ._gff['fsource'] = "ABGPstl3nCONS"
            intronQ._apps_donor = 1.0
            intronQ._apps_accep = 1.0
            intronS._linked_to_introns = []
            intronS._linked_to_pacbporfs = []
            intronS._gff['fsource'] = "ABGPstl3nCONS"
            intronS._apps_donor = 1.0
            intronS._apps_accep = 1.0
    
        # label tags the interface of PacbPORFs for the stopless3n introns in query
        target_label    = (pacbporf.orfQ.id,pacbporf.orfQ.id)
        informant_label = (pacbporf.orfS.id,pacbporf.orfS.id)
        target          = OPTIONS.target
        informant       = n2[0]

        # update conserved QUERY stopless3n introns
        _update_observed_introns([q for (q,s) in conservedstopless],
                target_label,informant_label,target,informant,
                observed_introns,intron2label,assessed_interfaces )
        # update conserved SBJCT stopless3n introns
        _update_observed_introns([s for (q,s) in conservedstopless],
                informant_label,target_label,informant,target,
                observed_introns,intron2label,assessed_interfaces )

################################################################################
### REMOVE -novel- stopless3n introns which have been found only once
################################################################################
MIN_SOLELY_OBSERVED_STOPLESS3N_COUNT = 2
for org,data in observed_introns.iteritems():
    intron_keys = data.keys()
    intron_keys.sort()
    for kk in intron_keys:
        if len(data[kk]) < MIN_SOLELY_OBSERVED_STOPLESS3N_COUNT:
            sources = Set([ intron._gff['fsource'] for intron in data[kk] ])
            if not sources.difference(['ABGPstl3n','ABGPstl3nCONS']):
                # delete this stopless3n intron
                del( observed_introns[org][kk] )
                logf("DELETING stopless3n::", org, kk)


# logf message
logf("#", stw.lap(), "stopless3n introns predicted",PCG)


################################################################################
### remove non-convincing stopless3n introns in terms of alignment presence
################################################################################
dna2protlength      = len(input[OPTIONS.target]['genomeseq'])/3
array_algpresence   = PCG2codingarray(PCG,OPTIONS.target,dna2protlength)
array_algsimilarity = PCG2similarityarray(PCG,OPTIONS.target,dna2protlength)

intron_keys = observed_introns[OPTIONS.target].keys()
intron_keys.sort()
removed_cnt = 0
for kk in intron_keys:
    # take any intron as example
    thisintron = observed_introns[OPTIONS.target][kk][0]
    if not thisintron.is_stopless_3n_intron(): continue
    # get presence & similarity score
    intron_aa_sta   = thisintron.start/3
    intron_aa_end   = thisintron.end/3
    intron_cnt      = len(observed_introns[OPTIONS.target][kk])
    intron_aa_length= thisintron.length/3
    intron_score    = intron_cnt * intron_aa_length
    presence_score  = sum(array_algpresence[intron_aa_sta:intron_aa_end+1])
    similarity_score= sum(array_algsimilarity[intron_aa_sta:intron_aa_end+1])
    informants      = [ intron._reference for intron in observed_introns[OPTIONS.target][kk] ]
    # check if any of the informants which did not deliver a stopless3n intron
    # has a pacbporf gap on the position of this locus
    informant_gap_scores = []
    for informant in GENE_INFORMANT_SET.difference(informants):
        infcodingarray = pacbporflist2codingarray(
                PCG.get_pacbps_by_organisms(OPTIONS.target,informant),
                'query',dna2protlength)

        if len(infcodingarray[intron_aa_sta:intron_aa_end+1]) == 0:
            # stopless3n intron out of range of this informant
            # This informant cannot deliver evidence for this area
            pass

        else:
            gap_score = float(sum(infcodingarray[intron_aa_sta:intron_aa_end+1])) /\
                        float(len(infcodingarray[intron_aa_sta:intron_aa_end+1]))
            if gap_score == 0.0 and sum(infcodingarray[0:intron_aa_sta]) == 0:
                # adjust if out of aligned range; in fact, do not store this
                # gap score. This informant cannot deliver evidence for this area
                pass
            else:
                informant_gap_scores.append( gap_score )

    if informant_gap_scores:
        gap_average_score = sum(informant_gap_scores) / float(len(informant_gap_scores))
    else:
        gap_average_score = 0.0

    # Check if there is any Orf transition anywhere in this interface
    # There is no Orf transition in the stopless3n introns ;-)
    # So, check on the relevant interface between the inwpCBGs for
    # the informants that did not deliver a stopless3n intron, and see if
    # any of those have a (potential) Orf transition
    potential_orf_transition = True
    for pos in range(1,len(inwpcbgs)):
        prevInwpCBG, nextInwpCBG = inwpcbgs[pos-1:pos+1]
        if (OPTIONS.target, thisintron.orfDonor.id) != prevInwpCBG._get_target_node():
            continue
        if (OPTIONS.target, thisintron.orfAcceptor.id) != nextInwpCBG._get_target_node():
            continue
        prevMINSR = prevInwpCBG.minimal_spanning_range(organism=OPTIONS.target)
        nextMINSR = nextInwpCBG.minimal_spanning_range(organism=OPTIONS.target)
        if intron_aa_sta in prevMINSR and intron_aa_end in prevMINSR:
            potential_orf_transition = False
            break
        elif intron_aa_sta in nextMINSR and intron_aa_end in nextMINSR:
            potential_orf_transition = False
            break
        elif intron_aa_sta >= min(prevMINSR) and intron_aa_end <= max(nextMINSR):
            prevNodes = [ prevInwpCBG.get_organism_nodes(org)[0] for org in GENE_INFORMANT_SET.difference(informants).intersection(prevInwpCBG.organism_set()) ]
            nextNodes = [ nextInwpCBG.get_organism_nodes(org)[0] for org in GENE_INFORMANT_SET.difference(informants).intersection(nextInwpCBG.organism_set()) ]
            prevOrgs  = [ org for (org,orfid) in prevNodes ]
            nextOrgs  = [ org for (org,orfid) in nextNodes ]
            # only get nodes of informants that are shared between both
            prevNodes = [ prevInwpCBG.get_organism_nodes(org)[0] for org in Set(prevOrgs).intersection(nextOrgs) ]
            nextNodes = [ nextInwpCBG.get_organism_nodes(org)[0] for org in Set(prevOrgs).intersection(nextOrgs) ]
            if prevNodes and nextNodes and not\
            Set(prevNodes).symmetric_difference(nextNodes):
                potential_orf_transition = False
                break

    # logf data on this intron
    logf("#", kk, (thisintron.orfDonor.id, thisintron.orfAcceptor.id),
            potential_orf_transition,
            len(observed_introns[OPTIONS.target][kk]),
            (intron_score,  similarity_score, presence_score, intron_aa_length),
            gap_average_score, informants, informant_gap_scores
            )
    # check if we can remove this stopless3n intron
    if similarity_score > (intron_score*5):
        del( observed_introns[OPTIONS.target][kk] )
        logf("# DELETING stopless3n::", OPTIONS.target, kk)
        removed_cnt+=1
    elif presence_score > intron_score and gap_average_score  == 1.0:
        del( observed_introns[OPTIONS.target][kk] )
        logf("# DELETING stopless3n::", OPTIONS.target, kk)
        removed_cnt+=1
    elif presence_score > intron_score and not potential_orf_transition and\
    gap_average_score >= 0.70:
        del( observed_introns[OPTIONS.target][kk] )
        logf("# DELETING stopless3n::", OPTIONS.target, kk)
        removed_cnt+=1
    else:
        pass
        
logf("#", stw.lap(), removed_cnt, "stopless3n introns removed",PCG)


################################################################################
# Stopless3nIntronPolishing. MileStone logging
################################################################################
if not VERBOSE: stdoutManager.reset()
logf("##", stwMS.lap(), "Stopless3nIntronPolishing: performed")


################################################################################
### starting next algorithm step: Finalization by consolidating introns
################################################################################
if not VERBOSE: stdoutManager.buffer() # capture all STDOUT messages from deeper modules 
stwMS.stepname = "Finalization"
stw.stepname = "Finalization"

################################################################################
### print all the found introns
################################################################################
print "\n\n"
for org,data in observed_introns.iteritems():
    if org != OPTIONS.target: continue
    if not data: continue
    exonobjects = input[org]['gldobj'].as_exons(input=input[org])
    print ">>>", org, "%s OBSERVED intron(s)" % len(data)
    intron_keys = data.keys()
    intron_keys.sort()
    for kk in intron_keys:
        introns = data[kk]
        example_intron = introns[0]
        if example_intron.__class__.__name__ == 'IntronConnectingOrfs':
            example_intron.assign_bp_and_ppts()
            has_bp      = example_intron.branchpoint != None
            has_ppt     = example_intron.ppt5p != None or example_intron.ppt3p != None
            has_bp_dist = example_intron.get_branchpoint_nt_distance()
        else:
            has_bp  = None
            has_ppt = None
            has_bp_dist = None

        interface = _get_main_interface( [ i._label for i in introns ] )
        is_existing_exon = kk in [ (exonobjects[i].start,exonobjects[i].end) for i in range(1,len(exonobjects),2) ]
        print kk, is_existing_exon, len(introns), interface,
        print has_bp, has_bp_dist, has_ppt,
        print assessed_interfaces[org][intron2label[org][kk]],
        print "%1.2f" % _calc_aligned_score(introns),
        print [ (i._reference, i._distance, i._gff['fsource']) for i in introns ]

logf("#", stw.lap(), "all ABFGP predicted introns listed")


################################################################################
### process the found introns
################################################################################
for org,data in observed_introns.iteritems():
    if org != OPTIONS.target: continue
    if not data: continue

    print ">>>", org, "%s OBSERVED intron(s)" % len(data)
    data = _filter_introns(data,input,PCG,ORGANISMS,OPTIONS)
    print ">>>", org, "%s FILTERED intron(s) introns" % len(data)
    data = _filter_introns_for_exon_msr(data,input,PCG,OPTIONS)
    print ">>>", org, "%s FILTERED intron(s) MSR" % len(data)

    intron_keys = data.keys()
    intron_keys.sort()
    for kk in intron_keys: print kk    

    ##print ">>>", org, "%s ANNOTATED intron(s)" % ((len(exonobjects)-1)/2)
    ##for i in range(0,len(exonobjects),2):     
    ##    exon = exonobjects[i]
    ##    print "\t", exon
    ##    print "\t", exon.orf.getaas(abs_pos_start=exon.acceptor.pos/3,abs_pos_end=exon.donor.pos/3)
    ##for i in range(1,len(exonobjects),2):     
    ##    existing_intron = exonobjects[i]
    ##    print existing_intron.start,existing_intron.end, existing_intron

    intron_keys = data.keys()
    intron_keys.sort()
    for kk in intron_keys:
        introns = data[kk]
        is_existing_exon = kk in [ (exonobjects[i].start,exonobjects[i].end) for i in range(1,len(exonobjects),2) ]
        print kk, is_existing_exon, len(introns), assessed_interfaces[org][intron2label[org][kk]],
        print "%1.2f" % _calc_aligned_score(introns), [ (i._reference, i._distance, i._gff['fsource']) for i in introns ],
        print introns[0]._label


    # add tinyexonPacbPORFS to the PCG if assigned by introns
    intron_keys = data.keys()
    intron_keys.sort()
    for kk in intron_keys:
        introns = data[kk]
        for intron in introns:
            if hasattr(intron,'_linked_to_pacbporfs'):
                for tinypacbporf in intron._linked_to_pacbporfs:
                    nodeQ = ( intron._organism,  tinypacbporf.orfQ.id )
                    nodeS = ( intron._reference, tinypacbporf.orfS.id )
                    if intron.__class__.__name__ == 'SequenceErrorConnectingOrfs':
                        # store ALL tinypacbporfs associated with SequenceError objects!
                        pass
                    elif max([ tinypacbporf.overlap(pF) for pF in PCG.get_pacbps_by_organisms(intron._organism,intron._reference) ]) == 1.0:
                        # fully included in an already listed PacbPORF -> ignore
                        continue

                    pkey  = tinypacbporf.construct_unique_key(nodeQ,nodeS)
                    if not nodeQ in PCG.get_nodes():
                        PCG.add_node(nodeQ)
                    if not nodeS in PCG.get_nodes():
                        PCG.add_node(nodeS)
                    if not PCG.has_edge(nodeQ,nodeS):
                        PCG.add_edge(nodeQ,nodeS,wt=tinypacbporf.bitscore)
                        PCG.pacbps[(pkey,nodeQ,nodeS)] = tinypacbporf
                        logf("#","tinypacbporf added:",nodeQ,nodeS,tinypacbporf)
                    else:
                        # edge already exists! Check if this (tiny) PacbPORF
                        # overlaps with an existing one
                        tinypacbporf_is_overlapping = False
                        for otherpacbp in PCG.get_pacbps_by_nodes(node1=nodeQ,node2=nodeS):
                            if pacb.ordering.issubsetorsuperset(otherpacbp,tinypacbporf):
                                tinypacbporf_is_overlapping = True
                                break
                        if not tinypacbporf_is_overlapping:
                            # add this tinypacbporf to the PCG
                            PCG.pacbps[(pkey,nodeQ,nodeS)] = tinypacbporf
                            logf("#","tinypacbporf added:",nodeQ,nodeS,tinypacbporf)


    # create gff lines for failed assessed interfaces
    gfflines = _create_failed_intron_gff(data,assessed_interfaces[org],intron2label[org])
    input[OPTIONS.target]['gffs'].extend( gfflines )


#################################################################################
### Now we have **the** PCG complete. Print it
################################################################################
# recreate inwpcbgs
inwpcbgs = PCG2inwpCBGS(PCG)
if VERBOSE:
    print PCG
    print_inwpcbgstructure(inwpcbgs,ORGANISMS)

logf("## gene-orfmodel:\t%s\t" % OPTIONS.target, input[OPTIONS.target]['orfid-genestructure'], input[OPTIONS.target]['proteinfref'] )
for org in GENE_INFORMANT_SET:
    if input[org]['orfid-genestructure']:
        logf("## gene-orfmodel:\t%s\t" % org, input[org]['orfid-genestructure'], input[org]['proteinfref'] )
for org in UNIGENE_INFORMANT_SET:
    if input[org]['orfid-unigenestructure']:
        logf("## UniGene-orfmodel:\t%s\t" % org, input[org]['orfid-unigenestructure'])

        
################################################################################
# Finalization. MileStone logging
################################################################################
if not VERBOSE: stdoutManager.reset()
logf("##", stwMS.lap(), "Finalization: performed")


################################################################################
### Prepare output files & database storage
################################################################################
if OPTIONS.output_creategff or OPTIONS.output_storetodb:
    # create GFF file
    if not VERBOSE: stdoutManager.buffer() # capture all STDOUT messages from deeper modules 
    introndata = observed_introns[OPTIONS.target]
    FREF = input[OPTIONS.target]['proteinfref']
    gff_fname   = create_intron_mapping_gff_file(input,PCG,introndata,
                assessed_interfaces,OPTIONS,FREF,organism=OPTIONS.target,
                verbose=True)
    if not VERBOSE: stdoutManager.reset()
    logf("##", stwMS.lap(),"result gff file:", gff_fname )

if OPTIONS.output_storetodb:
    # rewrite the fasta file -> header must equal used_fref variable
    fasta_fname = input[OPTIONS.target]['genomeseqfile']
    used_fref = open(gff_fname).readline().split("\t")[0]
    header, seq, descr = parseSingleFasta(open(fasta_fname).readlines())
    fasta_fname = osPathJoin(OPTIONS.outdir,used_fref+".fa")
    writeMultiFasta( {used_fref: seq }, fasta_fname )

    # if here: get fref, cleanup db & store2db
    dbcleanup(fref=used_fref)
    # grab status of storing succes; if one fails, all fails ...
    print "STORING GFF:  ", gff_fname
    print "STORING FASTA:", fasta_fname 
    status = gff2database(gff_fname,fasta_fname)
    # and remove the tmp-created fasta file
    _file_cleanup([fasta_fname])

    if not OPTIONS.output_creategff:
        _file_cleanup([gff_fname])

################################################################################
### Ready with ABFGP!
################################################################################

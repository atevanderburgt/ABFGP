"""
This library contains the optparse classes and functions for the ABGP suite of programs.
Command line arguments are - where applicable - standardized and therefor categorized.
In the ABGP programs themselves, a parser object is constructed by creating a parser
with the main function, and subsequently assemble the proper optiongroups by calling their function.

  def abgpoptparser()

Besides the main parser function, there are X pairs of functions with OptionGoup parser data and validation functions

  def generaloptions(parser)                    def validate_generaloptions(parser,options)
  def abgpinputoptions(parser)                  def validate_abgpinputoptions(parser,options)
  def abgpoptions(parser)                       def validate_abgpoptions(parser,options)
  def abgpinformantselectionoptions(parser)     def validate_abgpinformantselectionoptions(parser,options)
  def abgpoutputoptions(parser)                 def validate_abgpoutputoptions(parser,options)
  def dbwarehousesearchoptions(parser)          def validate_dbwarehousesearchoptions(parser,options)
  def dbwarehousesettingsoptions(parser)        def validate_dbwarehousesettingsoptions(parser,options)
  def dbwarehouseresultoptions(parser)          def validate_dbwarehouseresultoptions(parser,options)
  def dbwarehousefunctionoptions(parser)        def validate_dbwarehousefunctionoptions(parser,options)
  def genelocustatsoptions(parser)              def validate_genelocustatsoptions(parser,options)
  def abgpclustalwoptions(parser)

"""

# Python imports
import os, sys, time
from random import randint
# get command line option parsing classes and functions
from optparse import OptionParser, OptionGroup
from pythonlibs.optparsefunctions import *

# Abgp imports
from lib_fasta import parseSingleFastaHeaderFromFile, IsSingleFastaDna, IsSingleFastaProtein, parseDecoratedFasta
from lib_filevalidation import fastafilesequencelength
from abgpgenelocusdirectory import AbgpGeneLocusDirectory, IsAbgpGeneLocusDirectory
from abgpdbwarehouseminer import AbgpDbwarehouseMiner 

# import settings (for default option values)
from settings.dbwarehouse import *
from settings.genetreegraph import *
from settings.abgp import *


def abgpoptparser():
    """ """
    # initialize and return optparse object
    return OptionParser()

# end of function abgpoptparser


def genelocustatsoptions(parser):
    """
    Add ABGP genelocusdirectory statistics to parser object

    @type  parser: optparse parser object
    @param parser: optparse parser object
    """
    parser.add_option("-d", "--dir", dest="DIRECTORY", default="",
                type="str",
                action="callback",
                callback=isdir_callback,
                help="perform analyses on this AbgpGeneLocus Directory")
    parser.add_option("-v", "--verbose",
                dest="verbose",
                default=False,
                action="store_true",
                help="print as much intermediate data to STDOUT as possible")
    parser.add_option("-q", "--quiet",
                dest="quiet",
                default=False,
                action="store_true",
                help="print only major Warnings to stdout")
    parser.add_option("-e", "--explain",
                dest="EXPLAIN",
                default=False,
                action="store_true",
                help="explain this program & its output")

# end of function genelocustatsoptions


def validate_genelocustatsoptions(parser,options):
    """ 
    """
    # proces directory parameter 
    options.CURRENTWORKINGDIR = os.getcwd()
    if not options.EXPLAIN:
        #options.DIRECTORY = os.path.abspath(options.DIRECTORY)
        if IsAbgpGeneLocusDirectory(options.DIRECTORY):
            pass
        else:
            print "Exception.UnexpectedFileType.NoAbgpGeneLocusDirectory", "[-d|--dir]", options.DIRECTORY 
            sys.exit()

    # check stdout status message printing
    if options.quiet:
        options.verbose = False

# end of function validate_genelocustatsoptions 


def generaloptions(parser):
    """
    Add ABGP general options to parser object

    @type  parser: optparse parser object
    @param parser: optparse parser object
    """
    ####################################################################
    # general option group
    ####################################################################
    generalOptionGroup = OptionGroup(parser, "General Options")
    generalOptionGroup.add_option("-v", "--verbose",
                dest="verbose",
                default=False,
                action="store_true",
                help="print as much status messages to stdout as possible")
    generalOptionGroup.add_option("-q", "--quiet",
                dest="quiet",
                default=False,
                action="store_true",
                help="print only mayor milestone messages to stdout")
    generalOptionGroup.add_option("-s", "--silent",
                dest="silent",
                default=False,
                action="store_true",
                help="print not a single status messages to stdout, just the (final) outcome")
    generalOptionGroup.add_option("-o", "--outdir",
                dest="outdir",
                default=None,
                type="str",
                help="specify outdir where to create result file(s); "
                "auto-generated when not specified")
    generalOptionGroup.add_option("-f", "--force",
                dest="force",
                default=False,
                action="store_true",
                help="force file-overwrite of existing files and creation of new directories")
    parser.add_option_group(generalOptionGroup)

# end of function generaloptions


def abgpinputoptions(parser):
    """
    Add ABGP input options to parser object

    @type  parser: optparse parser object
    @param parser: optparse parser object
    """
    ####################################################################
    # input AbgpLocus, LocusDirectory and DNA files
    ####################################################################
    mainOptionGroup = OptionGroup(parser, "Main Options for assigning sequences for ABGP")
    mainOptionGroup.add_option("--loci", dest="loci", default=[],
                action="callback",
                callback=vararg_callback,
                help="perform gene prediction on these AbgpGeneLocus "
                "directorie(s): --loci <dirA> .. <dirX> ")
    mainOptionGroup.add_option("--dirwithloci", dest="dirwithloci", default=None,
                type="str",
                action="callback",
                callback=isdir_callback,
                help="perform gene prediction on this (preselected) LocusDirectory "
                "with AbgpGeneLocus directories")
    mainOptionGroup.add_option("--filewithloci", dest="filewithloci", default=None,
                type="str",
                action="callback",
                callback=isfile_callback,
                help="perform gene prediction on this file with absolute "
                "paths to AbgpGeneLocus directories")
    mainOptionGroup.add_option("--dna", dest="dnafiles", default=[],
                action="callback",
                callback=vararg_callback,
                help="perform gene prediction on single fasta dna file(s): --dna <fileA> ... <fileX> "
                "[length range %skb..%skb]" % (ABGP_MINIMAL_NT_LOCUS_LENGTH/1000,ABGP_MAXIMAL_NT_LOCUS_LENGTH/1000))
    mainOptionGroup.add_option("--multifasta", dest="multifasta", default=None,
                type="str",
                action="callback",
                callback=isfile_callback,
                help="perform gene prediction on multifasta DNA file: --multifasta <fileA>"
                "[length range %skb..%skb]" % (ABGP_MINIMAL_NT_LOCUS_LENGTH/1000,ABGP_MAXIMAL_NT_LOCUS_LENGTH/1000))
    parser.add_option_group(mainOptionGroup)

# end of function abgpinputoptions


def abgpclustalwoptions(parser):
    """
    Add ABGP clustalw options to parser object

    @type  parser: optparse parser object
    @param parser: optparse parser object
    """
    clustalwOptionGroup = OptionGroup(parser, "ABGP ClustalW specific Options")
    clustalwOptionGroup.add_option("--fasta",
                    dest="fasta",
                    default=[],
                    action="callback",
                    callback=vararg_callback,
                    help="include these genes/proteins for ClustalW"
                    " --fasta <proteinA> ... <proteinX>")
    clustalwOptionGroup.add_option("--outputproteins",
                    dest="outputproteins",
                    default=False,
                    action="store_true",
                    help="output protein multifasta input, not ClustalW output")
    clustalwOptionGroup.add_option("--omit_informants",
                    dest="omit_informants",
                    default=[],
                    action="callback",
                    callback=vararg_callback,
                    help="omit these genes/proteins for ClustalW"
                    " --omit_informants <infA> ... <infX>")
    clustalwOptionGroup.add_option("--use_informants",
                    dest="use_informants",
                    default=[],
                    action="callback",
                    callback=vararg_callback,
                    help="only use these genes/proteins for ClustalW"
                    " --use_informants <infA> ... <infX>")
    parser.add_option_group(clustalwOptionGroup)

# end of function abgpclustalwoptions


def abgpinformantselectionoptions(parser):
    """
    Add ABGP informant selection options to parser object

    @type  parser: optparse parser object
    @param parser: optparse parser object
    """
    ####################################################################
    # input AbgpLocus, LocusDirectory and DNA files
    ####################################################################
    infselOptionGroup = OptionGroup(parser, "Informant Selection Options")
    infselOptionGroup.add_option("--omit_informant_selection",
                dest="omit_informant_selection",
                default=False,
                action="store_true",
                help="omit selection of most usefull informant genes")
    infselOptionGroup.add_option("--target",
                dest="target",
                default=None, 
                type="str",
                help="use this Organism/Gene identifier asa target "
                "(and all others as informants)")
    # TODO: not implemented yet!
    infselOptionGroup.add_option("--required_informants",
                dest="required_informants",
                default=[],
                action="callback",
                callback=vararg_callback,
                help="force these informants to be selected"
                "(from a longer list of of informant genes):"
                " --required_informants <infA> ... <infX>")
    infselOptionGroup.add_option("--omit_informants",
                dest="omit_informants",
                default=[],
                action="callback",
                callback=vararg_callback,
                help="remove these informants prior to informant selection"
                "(from a longer list of of informant genes):"
                " --omit_informants <infA> ... <infX>")
    infselOptionGroup.add_option("--minimal_num_loci","--min_num_loci",
                dest="minimal_num_loci", type="int",
                default=ABGP_MINIMAL_NUM_LOCI,
                action="callback",
                callback_kwargs={'minval':ABGP_MINIMAL_NUM_LOCI_LOWER,
                                 'maxval':ABGP_MINIMAL_NUM_LOCI_UPPER},
                callback=valueinrange_callback,
                help="minimal number of gene loci to use as input"
                " [default: %default, range "+\
                str(ABGP_MINIMAL_NUM_LOCI_LOWER)+".."+\
                str(ABGP_MINIMAL_NUM_LOCI_UPPER)+"]")
    infselOptionGroup.add_option("--maximal_num_loci","--max_num_loci",
                dest="maximal_num_loci", type="int",
                default=ABGP_MAXIMAL_NUM_LOCI,
                action="callback",
                callback_kwargs={'minval':ABGP_MAXIMAL_NUM_LOCI_LOWER,
                                 'maxval':ABGP_MAXIMAL_NUM_LOCI_UPPER},
                callback=valueinrange_callback,
                help="maximal number of gene loci to use as input;"
                " if more applied, overflow of informants is discarded"
                " [default: %default, range "+\
                str(ABGP_MAXIMAL_NUM_LOCI_LOWER)+".."+\
                str(ABGP_MAXIMAL_NUM_LOCI_UPPER)+"]")
    parser.add_option_group(infselOptionGroup)

# end of function abgpinformantselectionoptions


def abgpoptions(parser):
    """
    Add ABGP algorithm parameter options to parser object

    @type  parser: optparse parser object
    @param parser: optparse parser object
    """
    ####################################################################
    # ABGP algorithm flow options
    ####################################################################
    abgpflowOptionGroup = OptionGroup(parser, "Options for ABGP algorithm flow")
    abgpflowOptionGroup.add_option("--minimal_gtgidentity",
                dest="minimal_gtgidentity",
                default=GTG_MINIMAL_IDENTITY,
                type="float",
                action="callback",
                callback_kwargs={
                    'minval': GTG_MINIMAL_IDENTITY_MINVAL,
                    'maxval': GTG_MINIMAL_IDENTITY_MAXVAL },
                callback=valueinrange_callback,
                help="minimal required GeneTreeGraph identityscore"
                " [default: %default, range "+\
                GTG_MINIMAL_IDENTITY_MINVAL_STRREPR+\
                ".."+\
                GTG_MINIMAL_IDENTITY_MAXVAL_STRREPR+"]")
    abgpflowOptionGroup.add_option("--minimal_gtgweakestnode",
                dest="minimal_gtgweakestnode",
                default=GTG_MINIMAL_WEAKESTNODE,
                type="float",
                action="callback",
                callback_kwargs={
                    'minval': GTG_MINIMAL_WEAKESTNODE_MINVAL,
                    'maxval': GTG_MINIMAL_WEAKESTNODE_MAXVAL },
                callback=valueinrange_callback,
                help="minimal required contribution of weakest gene Node "
                "to the GeneTreeGraph's identityscore"
                " [default: %default, range "+\
                GTG_MINIMAL_WEAKESTNODE_MINVAL_STRREPR+\
                ".."+\
                GTG_MINIMAL_WEAKESTNODE_MAXVAL_STRREPR+"]")
    abgpflowOptionGroup.add_option("--abinitio",
                dest="abinitio",
                default=False,
                action="store_true",
                help="do not use the annotated gene structure as prior knowledge")
    abgpflowOptionGroup.add_option("--gtganalyses",
                dest="gtganalyses",
                default=False,
                action="store_true",
                help="do GTG analyses only on this group of genes;"
                " print result to STDOUT")
    abgpflowOptionGroup.add_option("--cbganalyses",
                dest="cbganalyses",
                default=False,
                action="store_true",
                help="do CBG analyses only on this group of genes;"
                " print result to STDOUT")
    abgpflowOptionGroup.add_option("--stopcodonanalyses",
                dest="stopcodonanalyses",
                default=False,
                action="store_true",
                help="do stopcodon analyses only on this group of genes; "
                "print result to STDOUT. "
                "When final CBG not unambigiously found, nothing is printed!")
    abgpflowOptionGroup.add_option("--informantanalyses", 
                dest="informantanalyses",             
                default=False,
                action="store_true",            
                help="do informant analyses only on this group of genes;"
                " print result to STDOUT (requires --target)")
    parser.add_option_group(abgpflowOptionGroup)

# end of function abgpoptions


def abgpspecialoptions(parser):
    """
    Add ABGP special/debugging/speedup parameter options to parser object

    @type  parser: optparse parser object
    @param parser: optparse parser object
    """
    ####################################################################
    # ABGP algorithm special/debugging/speedup flow options
    ####################################################################
    abgpspecialOptionGroup = OptionGroup(parser,
        "Special/debugging options for ABGP algorithm flow")
    abgpspecialOptionGroup.add_option("--shortcut",
                dest="shortcut",
                default=False,
                action="store_true",
                help="USEFULL FOR NEAR-IDENTICAL GENES AND DEBUGGING PURPOSES:"
                " iterative blast is omitted after the first blast step")
    abgpspecialOptionGroup.add_option("--genomemask",
                dest="genomemask",
                default=[],
                action="callback",
                callback=vararg_callback,
                help="USE WITH CARE! Mask part of the genomic DNA sequences;"
                    "format org:sta..end")
    abgpspecialOptionGroup.add_option("--onlyusegenestructureorfs",
                dest="onlyusegenestructureorfs",
                default=[],
                action="callback",
                callback=vararg_callback,
                help="USE WITH CARE!! Use only those Orfs that are present in"
                    "the annotated gene structure -> large speedup but "
                    "possible performance drop")
    abgpspecialOptionGroup.add_option("--verbose_steps",
                type='choice',
                dest="verbose_steps",
                default=[],
                choices=['D1','D2','D3','D4','D5','D6','D7','D8','D9','MILESTONE'],
                action="callback",
                callback=vararg_callback_choices,
                help="Verbose mode only for these steps")
    abgpspecialOptionGroup.add_option("--omitksminxcbgs","--omit_ksminxcbgs",
                dest="omitksminxcbgs",
                default=False,
                action="store_true",
                help="USE WITH CARE! Do not search for K(s-x) CBGs ->"
                    "large speedup but possible performance drop!")
    parser.add_option_group(abgpspecialOptionGroup)

# end of function abgpspecialoptions



def abgppairwisespecialoptions(parser):
    """
    Add ABGP PAIRWISE special/debugging/speedup parameter options to parser object

    @type  parser: optparse parser object
    @param parser: optparse parser object
    """
    ####################################################################
    # ABGP algorithm special/debugging/speedup flow options
    ####################################################################
    abgpspecialOptionGroup = OptionGroup(parser,
        "Special/debugging options for ABGP algorithm flow")
    abgpspecialOptionGroup.add_option("--shortcut",
                dest="shortcut",
                default=False,
                action="store_true",
                help="USEFULL FOR TARGET GENE WITH CORRECT GENE STRUCTURE:"
                "2th blast is omitted after the first blast step and "
                "only --maximal_num_informants are used from the start")
    abgpspecialOptionGroup.add_option("--onlyusegenestructureorfs",
                dest="onlyusegenestructureorfs",
                default=[],
                action="callback",
                callback=vararg_callback,
                help="USE WITH CARE!! Use only those Orfs that are present in "
                    "the annotated gene structure(s) -> large speedup but "
                    "possible performance drop")
    abgpspecialOptionGroup.add_option("--omithmm",
                dest="omithmm",
                default=False,
                action="store_true",
                help="omit HMM search steps to obtain pairwise alignments " 
                    "-> large speedup but possible performance drop")
    abgpspecialOptionGroup.add_option("--omittinyexon",
                dest="omittinyexon",
                default=False,
                action="store_true",
                help="omit TinyExon search steps to obtain pairwise alignments "
                    "-> large speedup but possible performance drop")
    abgpspecialOptionGroup.add_option("--abinitio",
                dest="abinitio",
                default=False,
                action="store_true",
                help="do not use the TARGET annotated gene structure as "
                    "prior knowledge")
    parser.add_option_group(abgpspecialOptionGroup)

# end of function abgppairwisespecialoptions


def abgppairwiseinformantoptions(parser):
    """
    Add ABGP PAIRWISE informant options to parser object

    @type  parser: optparse parser object
    @param parser: optparse parser object
    """
    infOptionGroup = OptionGroup(parser, "Informant Options")
    infOptionGroup.add_option("--target",
                dest="target",
                default=None, 
                type="str",
                help="use this Organism/Gene identifier asa target "
                "(and all others as informants)")
    # TODO: not implemented yet!
    infOptionGroup.add_option("--required_informants",
                dest="required_informants",
                default=[],
                action="callback",
                callback=vararg_callback,
                help="force these informants to be selected "
                "(from a longer list of of informant genes):"
                " --required_informants <infA> ... <infX>")
    infOptionGroup.add_option("--omit_informants",
                dest="omit_informants",
                default=[],
                action="callback",
                callback=vararg_callback,
                help="remove these informants prior to informant selection"
                "(from a longer list of of informant genes):"
                " --omit_informants <infA> ... <infX>")
    infOptionGroup.add_option("--minimal_num_loci","--min_num_loci",
                dest="minimal_num_loci", type="int",
                default=5,
                action="callback",
                callback_kwargs={'minval':2,
                                 'maxval':50},
                callback=valueinrange_callback,
                help="minimal number of gene loci to use as input"
                " [default: %default, range "+\
                str(2)+".."+\
                str(50)+"]")
    infOptionGroup.add_option("--maximal_num_loci","--max_num_loci",
                dest="maximal_num_loci", type="int",
                default=15,
                action="callback",
                callback_kwargs={'minval':4,
                                 'maxval':50},
                callback=valueinrange_callback,
                help="maximal number of gene loci to use as input;"
                " if more applied, most usefull informants are selected"
                " [default: %default, range "+\
                str(4)+".."+\
                str(50)+"]")
    infOptionGroup.add_option("--disallow_informant_deletion",
                dest="disallow_informant_deletion",
                default=False,
                action="store_true",
                help="disallow automatic deletion of (very) poor informants")
    parser.add_option_group(infOptionGroup)

# end of function abgppairwiseinformantoptions



def abgpoutputoptions(parser):
    """
    Add ABGP algorithm parameter options to parser object

    @type  parser: optparse parser object
    @param parser: optparse parser object
    """
    ####################################################################
    # GGB & GFF options
    ####################################################################
    ggbandgffOptionGroup = OptionGroup(parser, "Options for GFF and GGB accessibility")
    ggbandgffOptionGroup.add_option("--storetoggbdb", dest="output_storetodb",
                default=False,
                action="store_true",
                help="store results in GGB (MYSQL) database [default: False]; "
                "specify DB_CONNECTION DB_NAME DB_HOST DB_USER DB_PASS in settings.ggb")
    ggbandgffOptionGroup.add_option("--nodetailedgff", dest="output_creategff",
                default=True,
                action="store_false",
                help="do NOT create detailed gff output file [ default: False ]")
    ggbandgffOptionGroup.add_option("--nofasta", dest="output_createfasta",
                default=True,
                action="store_false",
                help="do NOT create protein fasta output file [ default: False ]")
    ggbandgffOptionGroup.add_option("--storelocitoggbdb", dest="output_storelocitoggbdb",
                default=False,
                action="store_true",
                help="store raw AbgpGeneLoci data in GGB (MYSQL) database [default: False]; "
                "specify DB_CONNECTION DB_NAME DB_HOST DB_USER DB_PASS in settings.ggb")
    parser.add_option_group(ggbandgffOptionGroup)

# end of function abgpoutputoptions


def abgpunigeneoptions(parser):
    """
    Add ABGP unigene options to parser object

    @type  parser: optparse parser object
    @param parser: optparse parser object
    """
    ####################################################################
    # ABGP algorithm unigene benchmark options
    ####################################################################
    abgpUgOptionGroup = OptionGroup(parser,
        "Unigene options for ABGP")
    abgpUgOptionGroup.add_option("--unigenebenchmark",
                dest="unigenebenchmark",
                default=False,
                action="store_true",
                help="benchmark ABGP predicted GeneStructure on unigene")
    abgpUgOptionGroup.add_option("--omit_unigene_informants",
                dest="omit_unigene_informants",
                default=False,
                action="store_true",
                help="do NOT add unigenes as ABGP informants")
    parser.add_option_group(abgpUgOptionGroup)

# end of function abgpunigeneoptions


def dbwarehousefunctionoptions(parser):
    """
    Add DbWarehouse gzip/gunzip options for the DbWarehouse to parser object

    @type  parser: optparse parser object
    @param parser: optparse parser object
    """
    proceswarehouseOptionGroup = OptionGroup(parser, "Options for gzip/gunzip parts of the dbwarehouse for retrieval speedup")
    proceswarehouseOptionGroup.add_option("--gzip_organism",
                dest="gzip_organism",
                default=None,
                type="str",
                help="gzip all 'organism' flatfiles in dbwarehouse -> storage space --, retrieval speed --")
    proceswarehouseOptionGroup.add_option("--gunzip_organism",
                dest="gunzip_organism",
                default=None,
                type="str",
                help="gunzip all 'organism' flatfiles in dbwarehouse -> storage space ++, retrieval speed ++")
    proceswarehouseOptionGroup.add_option("--gzip_all",
                dest="gzip_all",
                default=False,
                action="store_true",
                help="gzip complete dbwarehouse -> storage space --, retrieval speed --")
    proceswarehouseOptionGroup.add_option("--gunzip_all",
                dest="gunzip_all",
                default=False,
                action="store_true",
                help="gunzip complete dbwarehouse -> storage space ++, retrieval speed ++")
    parser.add_option_group(proceswarehouseOptionGroup)

# end of function dbwarehousefunctionoptions 


def dbwarehouseresultoptions(parser):
    """
    Options for result of mining the DbWarehouse to parser object

    @type  parser: optparse parser object
    @param parser: optparse parser object
    """
    resultswarehouseOptionGroup = OptionGroup(parser, "Options for creating output directory for dbwarehouse search results")
    resultswarehouseOptionGroup.add_option("--createoutputdir",
                dest="createoutputdir",
                default=False,
                action="store_true",
                help="create output directory for DbwarehouseMiner results")
    resultswarehouseOptionGroup.add_option("--filewithloci",
                dest="filewithloci",
                default=False,
                action="store_true",
                help="print --filewithloci input file for ABFGP to STDOUT")
    parser.add_option_group(resultswarehouseOptionGroup)

# end of function dbwarehouseresultoptions 


def dbwarehousesearchoptions(parser):
    """
    Add search options for the DbWarehouse to parser object

    @type  parser: optparse parser object
    @param parser: optparse parser object
    """
    ####################################################################
    # search & find proteins/genes/identifiers options
    ####################################################################
    searchwarehouseOptionGroup = OptionGroup(parser, "Options for searching for a gene/protein identifier or sequence in dbwarehouse")
    searchwarehouseOptionGroup.add_option("--locus",
                dest="locus",
                default=None,
                type="str",
                action="callback",
                callback=isdir_callback,
                help="find matching gene loci to this AbgpGeneLocusDirectory"
                " -> then perform autoABGP")
    searchwarehouseOptionGroup.add_option("--identifier",
                dest="identifier",
                default=None,
                type="str",
                help="find identifier in ABGP sequence repository and find matching gene loci"
                " -> then perform autoABGP")
    searchwarehouseOptionGroup.add_option("--protein",
                dest="protein",
                default=None,
                type="str",
                action="callback",
                callback=isfile_callback,
                help="blastp single fasta protein to dbwarehouse to find matching gene loci"
                " -> then perform autoABGP")
    searchwarehouseOptionGroup.add_option("--dnasequence",
                dest="dnasequence",
                default=None,
                type="str",
                action="callback",
                callback=isfile_callback,
                help="blastx single fasta dna to dbwarehouse to find matching gene loci"
                " -> then perform autoABGP [range 2kb..15kb]")
    searchwarehouseOptionGroup.add_option("--genome",
                dest="genome",
                default=None,
                type="str",
                action="callback",
                callback=isfile_callback,
                help="gzipped multi-fasta genome to abstract a gene locus from. "
                "When combined with --identifier: -> protein -> tblastn -> "
                "create genelocus on genome -> then perform autoABGP. "
                "When combined with --protein: -> tblastn -> create genelocus "
                "on genome -> then perform autoABGP.")
    parser.add_option_group(searchwarehouseOptionGroup)

# end of function dbwarehousesearchoptions


def dbwarehousesettingsoptions(parser):
    """
    Add settings options for the DbWarehouse to parser object

    @type  parser: optparse parser object
    @param parser: optparse parser object
    """
    ####################################################################
    # dbwarehouse usage options
    ####################################################################
    dbwarehouseOptionGroup = OptionGroup(parser, "Options for adjusting dbwarehouse settings")
    dbwarehouseOptionGroup.add_option("--db",
                dest="dbwarehouse",
                default=DBWAREHOUSE_PATH,
                help="location of ABGP sequence repository;"
                " default in settings.dbwarehouse"
                " [default: '%default']")
    dbwarehouseOptionGroup.add_option("--search_mode",
                dest="search_mode",
                default=DBWAREHOUSE_SEARCH_METHOD,
                choices=DBWAREHOUSE_SEARCH_METHOD_LIST,
                help="search mode in warehouse [default: %default]. "
                "Searchmodes are SAFEORTHOLOGS, BDBH and SIMILARITY. "
                "Searchmode SAFEORTHOLOGS: BDBH approach combined with minimal distance towards closest hypothetical paralog (--safeorthologs_ratio). "
                "Default settings ensure selection of only 1:1:..:1:1 orthologous proteins. "
                "Searchmode BDBH: basic BiDirectionalBestHit approach. Paralogs are excluded, "
                "except when --allow_paralogs is given. "
                "Searchmode SIMILARITY: basic similarity approach. Paralogs are always used. "
                )
    dbwarehouseOptionGroup.add_option("--minimal_overlap_ratio",
                dest="minimal_overlap_ratio", type="float",
                default=DBWAREHOUSE_MINIMAL_OVERLAP_RATIO,
                action="callback",
                callback_kwargs={
                    'minval': DBWAREHOUSE_MINIMAL_OVERLAP_RATIO_MINVAL,
                    'maxval': DBWAREHOUSE_MINIMAL_OVERLAP_RATIO_MAXVAL },
                callback=valueinrange_callback,
                help="minimal required overlap ratio between 2 proteins"
                " [default: %default, range "+\
                DBWAREHOUSE_MINIMAL_OVERLAP_RATIO_MINVAL_STRREPR+\
                ".."+\
                DBWAREHOUSE_MINIMAL_OVERLAP_RATIO_MAXVAL_STRREPR+"]")
    dbwarehouseOptionGroup.add_option("--minimal_bitscore_ratio",
                dest="minimal_bitscore_ratio", type="float",
                default=DBWAREHOUSE_MINIMAL_BITSCORE_RATIO,
                action="callback",
                callback_kwargs={
                    'minval': DBWAREHOUSE_MINIMAL_BITSCORE_RATIO_MINVAL,
                    'maxval': DBWAREHOUSE_MINIMAL_BITSCORE_RATIO_MAXVAL },
                callback=valueinrange_callback,
                help="minimal required bitscore ratio between 2 proteins"
                " [default: %default, range "+\
                DBWAREHOUSE_MINIMAL_BITSCORE_RATIO_MINVAL_STRREPR+\
                ".."+\
                DBWAREHOUSE_MINIMAL_BITSCORE_RATIO_MAXVAL_STRREPR+"]")
    dbwarehouseOptionGroup.add_option("--maximal_length_ratio",
                dest="maximal_length_ratio", type="float",
                default=DBWAREHOUSE_MAXIMAL_LENGTH_RATIO,
                action="callback",
                callback_kwargs={
                    'minval': DBWAREHOUSE_MAXIMAL_LENGTH_RATIO_MINVAL,
                    'maxval': DBWAREHOUSE_MAXIMAL_LENGTH_RATIO_MAXVAL },
                callback=valueinrange_callback,
                help="maximal length ratio difference between 2 proteins"
                " [default: %default, range "+\
                DBWAREHOUSE_MAXIMAL_LENGTH_RATIO_MINVAL_STRREPR+\
                ".."+\
                DBWAREHOUSE_MAXIMAL_LENGTH_RATIO_MAXVAL_STRREPR+"]")
    dbwarehouseOptionGroup.add_option("--allow_paralogs",
                dest="allow_paralogs",
                default=False,
                action="store_true",
                help="allow selection of paralogous proteins")
    dbwarehouseOptionGroup.add_option("--safeorthologs_ratio",
                dest="safeorthologs_ratio", type="float",
                default=DBWAREHOUSE_SAFEORTHOLOGS_RATIO,
                action="callback",
                callback_kwargs={
                    'minval': DBWAREHOUSE_SAFEORTHOLOGS_RATIO_MINVAL,
                    'maxval': DBWAREHOUSE_SAFEORTHOLOGS_RATIO_MAXVAL },
                callback=valueinrange_callback,
                help="maximal ratio between 1th and 2th best hit to prevent putative inclusion of paralogs"
                " [default: %default, range "+\
                DBWAREHOUSE_SAFEORTHOLOGS_RATIO_MINVAL_STRREPR+\
                ".."+\
                DBWAREHOUSE_SAFEORTHOLOGS_RATIO_MAXVAL_STRREPR+"]"
                " (only used in combination with --allow_paralogs)")
    dbwarehouseOptionGroup.add_option("--maximal_num_loci",
                dest="maximal_num_loci", type="int",
                default=DBWAREHOUSE_DEFAULT_MAX_NUM_LOCI,
                action="callback",
                callback_kwargs={'minval':DBWAREHOUSE_MINIMAL_NUM_LOCI,'maxval':DBWAREHOUSE_MAXIMAL_NUM_LOCI},
                callback=valueinrange_callback,
                help="maximal number of most similar gene loci to start with"
                " -> then perform autoABGP"
                " [default: %default, range 2..25]")
    dbwarehouseOptionGroup.add_option("--minimal_num_loci",
                dest="minimal_num_loci", type="int",
                default=DBWAREHOUSE_MINIMAL_NUM_LOCI,
                action="callback",
                callback_kwargs={'minval':2,'maxval':7},
                callback=valueinrange_callback,
                help="minimal number of gene loci to use as input"
                " -> then perform autoABGP"
                " [default: %default, range 2..7]")
    dbwarehouseOptionGroup.add_option("--genomedirs_to_ignore",
                dest="genomedirs_to_ignore", default=[],
                action="callback",
                callback=vararg_callback,
                help="do not include proteins from these genome directories in the "
                "warehouse, represented by their short tag name(s): "
                "--genomedirs_to_ignore <dirtagA> .. <dirtagX> [ default: ignore none ]")
    dbwarehouseOptionGroup.add_option("--genomedirs_to_use",
                dest="genomedirs_to_use", default=[],
                action="callback",
                callback=vararg_callback,
                help="include only proteins from these genome directories in the "
                "warehouse, represented by their short tag name(s): "
                "--genomedirs_to_use <dirtagA> .. <dirtagX> [ default: use all ]")
    parser.add_option_group(dbwarehouseOptionGroup)

# end of function dbwarehousesettingsoptions


def validate_generaloptions(parser,options):
    """
    """
    if options.outdir:
        options.outdir = os.path.abspath(os.path.join(os.getcwd(),options.outdir))
        if not os.path.isdir(options.outdir):
            thisoption = parser.get_option('--outdir')
            print "Exception.OutputFolderDoesNotExist", "--outdir", options.outdir
            print "--outdir:", thisoption.help
            sys.exit()        
    else:
        # auto-define outdir!
        options.outdir = os.path.join(ABGP_OUTDIR_PATH,
            "%sR%s%s" % ( 
                time.strftime('abgpD%Y%m%dT%H%M%S'),
                str(time.time()).split('.')[1],
                randint(0,9) )
            )

    # check stdout status message printing
    if options.quiet:
        options.verbose = False
    if options.silent:
        options.quiet   = True
        options.verbose = False

# end of function validate_generaloptions


def validate_abgpinputoptions(parser,options):
    """
    """
    # are all applied options.loci indeed loci?
    for locusdir in options.loci:
        if IsAbgpGeneLocusDirectory(locusdir):
            pass
        else:
            print "Exception.UnexpectedFileType.NoAbgpGeneLocusDirectory", "[--locus|-d|--dir]", locusdir
            sys.exit()

    if options.filewithloci:
        # filename with absolute - recommended - full paths to loci directories is specified;
        # check if are indeed loci!
        IS_FIRST_FULLPATH = True
        for fullpath in open(options.filewithloci).readlines():
            fullpath = fullpath.strip()
            if not fullpath:
                continue
            elif fullpath and fullpath[0] == "#":
                continue
            if os.path.isdir(fullpath) and IsAbgpGeneLocusDirectory(fullpath):
                # check if already added (duplicate GLDs in file...
                if fullpath in options.loci: continue
                if IS_FIRST_FULLPATH and not options.target:
                    # assign this one as target gene
                    options.target = os.path.basename(fullpath).replace("/","")
                    IS_FIRST_FULLPATH = False
                # yes, an AbgpGeneLocusDirectory! append to options.loci
                options.loci.append( os.path.abspath(fullpath) )


    if options.dirwithloci:
        # directory with loci is specified, check if are indeed loci!
        for item in os.listdir(options.dirwithloci):
            fullpath = os.path.join(options.dirwithloci,item)
            if os.path.isdir(fullpath) and IsAbgpGeneLocusDirectory(fullpath):
                # yes, an AbgpGeneLocusDirectory! append to options.loci
                options.loci.append( fullpath )
                
    # options.multifasta; is it a multi-fasta DNA file, and are all loci of appropriate length?
    if options.multifasta:
        thisoption = parser.get_option('--multifasta')
        fname = options.multifasta
        # fname is an existing file?
        if not os.path.isfile(fname):
            print "Exception.FileDoesNotExist", "--multifasta", fname 
            print "--multifasta:", thisoption.help
            sys.exit()
        # fname IsSingleFastaDna?
        fname = os.path.abspath(fname)
        try:
            fastadict,hdrinfo = parseDecoratedFasta(open(fname))
        except:
            print "Exception.UnexpectedFileType.NoDnaMultiFastaFile", "--multifasta", fname
            print "--multifasta:", thisoption.help
            sys.exit()
        # check all the entries in the fasta file:
        # are they within the accepted length range?
        # are they DNA (ATGCN), contain IUPAC (URYMKWSBDHV)
        # in case `to much` (arbitrary threshold chosen here to be at most 5 characters) IUPAC -> Exception
        iupac_chars = list("URYMKWSBDHV")
        protein_chars = list("ARNDCQEGHILKMFPSTWYVBZX")
        headers = fastadict.keys()
        for hdr in headers:
            sequence = fastadict[hdr].upper()
            if not set(list(sequence)).difference(list("ATGCN")):
                seqlength = len(sequence)
                if not (seqlength >= ABGP_MINIMAL_NT_LOCUS_LENGTH and seqlength <= ABGP_MAXIMAL_NT_LOCUS_LENGTH):
                    # print warning message and continue but no sys.exit()
                    print "Warning.FastaLengthOutOfBounds: '%s' %snt, not %skb..%skb" % (hdr,seqlength,ABGP_MINIMAL_NT_LOCUS_LENGTH/1000,ABGP_MAXIMAL_NT_LOCUS_LENGTH/1000), "--multifasta", fname 
                    # discard this fastadict entry
                    del(fastadict[hdr])
                    continue
                else:
                    pass
            elif sum([ sequence.count(char) for char in iupac_chars ]) >= 5:
                # print warning message and continue but no sys.exit()
                print "Warning.ToMuchDegenerateDnaCharacters in '%s'" % hdr, "--multifasta", fname
                # discard this fastadict entry
                del(fastadict[hdr])
                continue
            else:
                print "Exception.Nonrecognized sequence type in '%s'. Might be RNA or PROTEIN, but not DNA" % hdr, "--multifasta", fname
                print "--multifasta:", thisoption.help
                sys.exit()                
            
        # Store the remaining fasta entries by overwriting options.multifasta
        options.multifasta = fastadict
                
    
    # options.dna are all IsSingleFastaDna?
    if options.dnafiles:
        thisoption = parser.get_option('--dna')
        for pos in range(0,len(options.dnafiles)):
            fname = options.dnafiles[pos]
            # fname is an existing file?
            if not os.path.isfile(fname):
                print "Exception.FileDoesNotExist", "--dna", fname 
                print "--dna:", thisoption.help
                sys.exit()
            # fname IsSingleFastaDna?
            fname = os.path.abspath(fname) 
            if not IsSingleFastaDna(fname):
                print "Exception.UnexpectedFileType.NoDnaSingleFastaFile", "--dna", fname
                print "--dna:", thisoption.help
                sys.exit()
            # sequence in fname in range (ABGP_MINIMAL_NT_LOCUS_LENGTH,ABGP_MAXIMAL_NT_LOCUS_LENGTH)?
            # default range is 2..20kb
            seqlength = fastafilesequencelength(fname)
            if not (seqlength >= ABGP_MINIMAL_NT_LOCUS_LENGTH and seqlength <= ABGP_MAXIMAL_NT_LOCUS_LENGTH):
                # print warning message and continue but no sys.exit()
                print "Warning.FastaLengthOutOfBounds: '%s' %snt, not %skb..%skb" % (hdr,seqlength,ABGP_MINIMAL_NT_LOCUS_LENGTH/1000,ABGP_MAXIMAL_NT_LOCUS_LENGTH/1000), "--multifasta", fname 
                continue
            # if here -> this is a correct dna locus file!
            options.dnafiles[pos] = fname


# end of function validate_abgpinputoptions


def validate_abgpoptions(parser,options):
    """
    """
    # check for presenence of both dirwithloci and filewithloci
    # although in fact there is no readon NOT to allow both,
    # it is deciced that only on of both options as input data is allowed
    if hasattr(options,"dirwithloci") and hasattr(options,"filewithloci"):
        if options.dirwithloci and options.filewithloci:
            message = " do not combine --dirwithloci and --filewithloci"
            print "Exception.WrongParameterCombination", message
            sys.exit()

    # check for presenence of none or either gtg/cbg/informant/stopcodon analyses
    analyses_attributes = [ 'cbganalyses', 'gtganalyses', 'stopcodonanalyses', 'informantanalyses' ]
    for attrA in analyses_attributes:
        if not hasattr(options,attrA): continue
        for attrB in analyses_attributes:
            if not hasattr(options,attrB): continue
            if attrA == attrB: continue
            # check if both are True
            if [ getattr(options,attrA), getattr(options,attrB) ] == [ True, True ]:
                message = " do not combine --%s and --%s" % (attrA,attrB) 
                print "Exception.WrongParameterCombination", message
                sys.exit()

    # any of the analyses attributes valid -> HARD-assign verbose/quiet/silent
    for attr in analyses_attributes:
        if not hasattr(options,attr): continue
        if getattr(options,attr) == True:
            # this requires silent mode!
            options.verbose = False
            options.quiet   = False
            options.silent  = True
            break

# end of function validate_abgpoptions


def validate_abgpspecialoptions(parser,options):
    """
    """
    pass

# end of function validate_abgpspecialoptions


def validate_abgpoutputoptions(parser,options):
    """
    """
    pass

# end of function validate_abgpoutputoptions


def validate_dbwarehousesearchoptions(parser,options):
    """
    """
    # check for presenence of none or either identifier / locus 
    if options.locus and options.identifier:
        print "Exception.WrongParameterCombination", " do not combine --locus and --identifier"
        sys.exit()

    # options.locus IsAbgpGeneLocusDirectory?
    if options.locus:
        if not IsAbgpGeneLocusDirectory(options.locus):
            # raise a friendly representation of NoAbgpGeneLocusDirectory
            print "Exception.UnexpectedFileType.NoAbgpGeneLocusDirectory", "--locus", options.locus
            sys.exit()
        else:
            locus = AbgpGeneLocusDirectory(options.locus)
            options.identifier = locus.protein_fref()
            print "UPDATED....NOOOOOOT"

    # options.protein IsSingleFastaProtein?
    if options.protein and not IsSingleFastaProtein(options.protein):
        print "Exception.UnexpectedFileType.NoProteinSingleFastaFile", "--protein", options.protein
        sys.exit()

    # options.dnasequence IsSingleFastaDna?
    if options.dnasequence:
        thisoption = parser.get_option('--dnasequence')
        if not IsSingleFastaDna(options.dnasequence):
            print "Exception.UnexpectedFileType.NoDnaSingleFastaFile", "--dnasequence", options.dnasequence
            print "--dnasequence:", thisoption.help
            sys.exit()
        # options.dnasequence in range 2kb..15kb?
        seqlength = fastafilesequencelength(options.dnasequence)
        if not (seqlength >= 2000 and seqlength <= 15000):
            print "incorrect length: %snt, not 2kb..15kb" % seqlength, "--dnasequence", options.dnasequence
            print "--dnasequence:", thisoption.help
            sys.exit()

    # options.identifier is in dbwarehouse?
    if options.identifier:
        dbminer = AbgpDbwarehouseMiner(optparseoptions=options)
        genomedir = dbminer.identifier2genomedir(options.identifier)
        if not genomedir:
            print "Exception.AbgpDbWarehouseMiner.IdentifierNotInWarehouse", "--identifier", options.identifier
            sys.exit()

        # find & set the locus of this identifier
        locusdir = dbminer.identifier2locusdir(options.identifier,genomedir=genomedir)
        options.locus = locusdir

    ## options.genome IsMultiFastaDna?
    #if options.genome and not IsMultiFastaDna(options.dnasequence,gzip=True):
    #    print "Exception.UnexpectedFileType.NoGzipDnaMultiFastaFile", "--dnasequence", options.dnasequence
    #    sys.exit()

# end of function validate_dbwarehousesearchoptions


def validate_abgpinformantselectionoptions(parser,options):
    """
    """
    if options.maximal_num_loci < options.minimal_num_loci:
        print "Exception.WrongParameterCombination",
        print "conflicting --maximal_num_loci (%s) < --minimal_num_loci (%s)" % (
                options.maximal_num_loci, options.minimal_num_loci )
        print ""
        print "--minimal_num_loci:",
        print parser.get_option('--minimal_num_loci').help
        print "--maximal_num_loci:",
        print parser.get_option('--maximal_num_loci').help
        sys.exit()

    if ( (options.informantanalyses and not options.target) or\
    ((options.omit_informant_selection or not options.target) ) and\
    ( len(options.loci) + len(options.dnafiles) ) >\
    options.maximal_num_loci):
        # to many GeneLoci, informant selection omitted...
        print "Exception.AbgpInputParameters.ToManyGeneLoci",
        print "--maximal_num_loci", options.maximal_num_loci,
        print "--target", options.target
        print ""
        print "Please:"
        print "    specify a target Gene\t(--target)"
        print "    increase maximal_num_loci\t(--maximal_num_loci)"
        if len(options.loci):
            print "    select less GeneLoci\t(--loci, --dirwithloci)"
        if len(options.dnafiles):
            print "    select less DNA sequences\t(--dna)" 
        print "Currently applied loci (%s) are:" % (
                len(options.loci) + len(options.dnafiles) )
        for locusabspath in options.loci:
            gld = AbgpGeneLocusDirectory(locusabspath)
            locuskey = gld._create_auto_key(
                identifier2organism=ABGP_LOCUS_IDENTIFIER2ORGANISM_MAPPING)
            print "    '%s'\t%s" % (locuskey,locusabspath)
        for dnaf in options.dnafiles: print "  ", dnaf
        print "--target:",
        print parser.get_option('--target').help
        print "--omit_informant_selection:",
        print parser.get_option('--omit_informant_selection').help
        print "--maximal_num_loci:",
        print parser.get_option('--maximal_num_loci').help
        sys.exit()

    # all informantselection parameters are okay
    return True

# end of function validate_abgpinformantselectionoptions


def validate_dbwarehouseresultoptions(parser,options):
    """
    """
    if not hasattr(options,"loci"):
        print "Exception.AbgpDbWarehouseMiner.NoMiningPerformed"
        return False
    elif len(options.loci) == 1:
        print "Exception.AbgpDbWarehouseMiner.NoSimilarSequencesFound", "--identifier", options.identifier
        return False
    elif len(options.loci) < options.minimal_num_loci:
        print "Exception.AbgpDbWarehouseMiner.ToFewSimilarSequencesFound", len(options.loci), "<", options.minimal_num_loci, "--identifier", options.identifier
        return False
    else:
        return True

# end of function validate_dbwarehouseresultoptions


def validate_dbwarehousesettingsoptions(parser,options):
    """
    """
    # check which genomedirs in dbwarehouse to include/exclude fromt he search
    if options.genomedirs_to_use and options.genomedirs_to_ignore:
        # hmm... dangerous usage -> only take genomedirs_to_use
        options.genomedirs_to_ignore = []
    
    # options.genomedirs_to_ignore are DbWareHouseGenomedirs? 
    if options.genomedirs_to_ignore:
        for i in range(0,len(options.genomedirs_to_ignore)):
            genomedir = options.genomedirs_to_ignore[i]
            if not os.path.isdir(os.path.join(options.dbwarehouse,genomedir)):
                thisoption = parser.get_option('--genomedirs_to_ignore')
                print "Exception.AbgpDbWarehouseMiner.NoGenomeDir", "--genomedirs_to_ignore", genomedir 
                print "--genomedirs_to_ignore:", thisoption.help
                sys.exit()
            else:
                options.genomedirs_to_ignore[i] = os.path.join(options.dbwarehouse,genomedir)
    
    # options.genomedirs_to_use are DbWareHouseGenomedirs?
    if options.genomedirs_to_use:
        for i in range(0,len(options.genomedirs_to_use)):
            genomedir = options.genomedirs_to_use[i]
            if not os.path.isdir(os.path.join(options.dbwarehouse,genomedir)):
                thisoption = parser.get_option('--genomedirs_to_use')
                print "Exception.AbgpDbWarehouseMiner.NoGenomeDir", "--genomedirs_to_use", genomedir
                print "--genomedirs_to_use:", thisoption.help
                sys.exit()
            else:
                options.genomedirs_to_use[i] = os.path.join(options.dbwarehouse,genomedir)

    if options.search_mode == 'SAFEORTHOLOGS':
        options.allow_paralogs = False
    elif options.search_mode == 'SIMILARITY':
        options.allow_paralogs = True 
    else:
        pass

# end of function validate_dbwarehousesettingsoptions


def validate_dbwarehousefunctionoptions(parser,options):
    """
    """
    if  [options.gzip_organism != None, options.gunzip_organism != None,
    options.gzip_all, options.gunzip_all].count(True) >= 2:
        print "Exception.WrongParameterCombination", " do not combine --g(un)zip_organism and --g(un)zip_all"
        sys.exit()

    # import the recreate_dbwarehouse_symlinks function
    from abgpdbwarehouseminer import (
        recreate_dbwarehouse_symlinks,
        remove_dbwarehouse_symlinks,
        )

    # options.gzip_organism
    if options.gzip_organism:
        if not os.path.isdir(os.path.join(options.dbwarehouse,options.gzip_organism)):
            thisoption = parser.get_option('--gzip_organism')
            print "Exception.AbgpDbWarehouseMiner.NoGenomeDir", "--gzip_organism", options.gzip_organism 
            print "--gzip_organism:", thisoption.help
            sys.exit()
        else:
            # remove all the secondary symmetrized files
            remove_dbwarehouse_symlinks(options.dbwarehouse,organism=options.gzip_organism)
            # gzip all the flatfiles of this organism 
            filepat1 = "blast.%s_x_%s.symmetrized" % (options.gzip_organism,options.gzip_organism) 
            filepat2 = "blast.%s_x_*.symmetrized" % (options.gzip_organism)
            filepat3 = "blast.*_x_%s.symmetrized" % (options.gzip_organism)
            for filepat in [ filepat1, filepat2, filepat3 ]:
                command = "gzip -f %s" % os.path.join(options.dbwarehouse,"_crossblastp",filepat)
                print command
                os.system(command)
            #### recreate the secondary symmetrized files
            ###recreate_dbwarehouse_symlinks(options.dbwarehouse,organism=options.gzip_organism)
            if options.verbose:
                print "--gzip_organism %s DONE" % options.gzip_organism
            sys.exit()

    # options.gunzip_organism
    if options.gunzip_organism:
        if not os.path.isdir(os.path.join(options.dbwarehouse,options.gunzip_organism)):
            thisoption = parser.get_option('--gunzip_organism')
            print "Exception.AbgpDbWarehouseMiner.NoGenomeDir", "--gunzip_organism", options.gunzip_organism
            print "--gunzip_organism:", thisoption.help
            sys.exit()
        else:
            # remove all the secondary symmetrized files
            remove_dbwarehouse_symlinks(options.dbwarehouse,organism=options.gzip_organism)
            # gunzip all the flatfiles of this organism
            filepat1 = "blast.%s_x_%s.symmetrized.gz" % (options.gunzip_organism,options.gunzip_organism)  
            filepat2 = "blast.%s_x_*.symmetrized.gz" % (options.gunzip_organism)
            filepat3 = "blast.*_x_%s.symmetrized.gz" % (options.gunzip_organism)
            for filepat in [ filepat1, filepat2, filepat3 ]:
                command = "gunzip %s" % os.path.join(options.dbwarehouse,"_crossblastp",filepat)
                print command
                os.system(command)
            #### recreate the secondary symmetrized files
            ###recreate_dbwarehouse_symlinks(options.dbwarehouse,organism=options.gzip_organism)
            if options.verbose:
                print "--gunzip_organism %s DONE" % options.gunzip_organism
            sys.exit()

    # options.gzip_all
    if options.gzip_all:
        # remove all the secondary symmetrized files
        remove_dbwarehouse_symlinks(options.dbwarehouse)
        # gzip all the flatfiles in the dbwarehouse
        command = "gzip -f %s" % os.path.join(options.dbwarehouse,"_crossblastp","blast.*_x_*.symmetrized")
        print command
        os.system(command)
        ##### recreate the secondary symmetrized files
        ### recreate_dbwarehouse_symlinks(options.dbwarehouse)
        if options.verbose:
            print "--gzip_all DONE"
        sys.exit()

    # options.gunzip_all
    if options.gunzip_all:
        # remove all the secondary symmetrized files
        remove_dbwarehouse_symlinks(options.dbwarehouse)
        # gunzip all the flatfiles in the dbwarehouse
        command = "gunzip %s" % os.path.join(options.dbwarehouse,"_crossblastp","blast.*_x_*.symmetrized.gz")
        print command
        os.system(command)
        #### recreate the secondary symmetrized files
        ###recreate_dbwarehouse_symlinks(options.dbwarehouse)
        if options.verbose:
            print "--gunzip_all DONE"
        sys.exit()

    # options.filewithloci and options.createoutputdir validation
    if options.filewithloci and options.createoutputdir:
        message = " do not combine --createoutputdir and --filewithloci"
        print "Exception.WrongParameterCombination", message
        sys.exit()

    if options.createoutputdir:
        if os.path.isdir(options.outdir) and options.force == False:
            print "OUTPUT DIRECTORY %s ALREADY EXIST; overwrite with --force" % options.outdir
            sys.exit()

# end of function validate_dbwarehousefunctionoptions


def validate_abgpunigeneoptions(parser,options,input):
    """
    @type  parser: optparse parser object
    @param parser: optparse parser object

    @type  options: parser.parse_args() options object
    @param options: parser.parse_args() options object

    @type  input: dict
    @param input: ABGP input data structure
    """
    if options.unigenebenchmark:
        if not options.target:
            message = " no --target specified"
            print "Exception.UnigeneBenchmarking", message
            sys.exit()
        elif not input.has_key(options.target):
            message = " --target not recognized"
            print "Exception.UnigeneBenchmarking", message
            sys.exit()
        elif not input[options.target]['orfid-unigenestructure']:
            message = " --target '%s' has no (full-length) unigene to benchmark" % ( options.target )
            print "Exception.UnigeneBenchmarking", message
            sys.exit()
        else:
            return True
    else:
        return True

# end of function validate_abgpunigeneoptions

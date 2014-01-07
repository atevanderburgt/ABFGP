#!/usr/bin/python
"""
Print informantion of given gene structures to STDOUT
"""

# ABGP imports
from abgpgenelocusdirectory import AbgpGeneLocusDirectory
from abgp_geneconfirmation import (
    readsequences,
    rungetorf,
    parseinputgff,
    annotatedgeneexonsizeevaluation,
    ntracksindnasequencecheck,
    confirmcanonicalsplicesites,
    )
from abgp_unigeneconfirmation import geneandunigeneconfirmation
from lib_sequencerepetitiveness import annotatedproteinsequencerepetitivenesscheck


# Import optparse functions and validators
from lib_optparse import (
    abgpoptparser,
    abgpinputoptions,
    validate_abgpinputoptions,
    validate_abgpoptions,
    )

# Imports (Python)
from sets import Set
from sys import exit as sysexit

# import settings
from settings.dbwarehouse import ABGP_LOCUS_IDENTIFIER2ORGANISM_MAPPING
from settings.gff.currentannotation import GFF_CDS_FMETHOD

################################################################################
# construct the command line option parser
# parse the command line & validate
################################################################################
parser = abgpoptparser()
abgpinputoptions(parser)
parser.remove_option('--dna')
parser.remove_option('--multifasta')
thisoption = parser.get_option('--loci')
thisoption.help = 'give statistics for these AbfgpGeneLocusDirectories: --loci <dirA> .. <dirX>'
thisoption = parser.get_option('--dirwithloci')
thisoption.help = 'give statistics for all the AbfgpGeneLocusDirectories in this --dirwithloci'
thisoption = parser.get_option('--filewithloci')
thisoption.help = 'give statistics for all the AbfgpGeneLocusDirectories listed in this --filewithloci'
(OPTIONS, args) = parser.parse_args()
OPTIONS.target = None # dummy option
validate_abgpinputoptions(parser,OPTIONS)
validate_abgpoptions(parser,OPTIONS)

################################################################################
################################################################################

if OPTIONS.loci:
    input = {} # input data structure, empty
    for locusdir in OPTIONS.loci:
        locus = AbgpGeneLocusDirectory(locusdir)
        locuskey = locus._create_auto_key(
                identifier2organism=ABGP_LOCUS_IDENTIFIER2ORGANISM_MAPPING
                )
        if input.has_key(locuskey):
            for suffix in list('abcdefghijklmnopqrstuvwxyz'):
                if not input.has_key(locuskey+suffix):
                    input.update( locus.toinputdict(key=locuskey+suffix) )
                    break
        else:
            input.update( locus.toinputdict(key=locuskey) )
else:
    print "Exception.NoAbgpGeneLocusDirectoriesProvided"
    print ""
    parser.print_help()
    sysexit()

################################################################################
# read input sequence files
################################################################################
input = readsequences(input)


################################################################################
# proces getorf data, input gff and confirm gene & unigene structure on Orfs
################################################################################
input = rungetorf(input)
input = parseinputgff(input)
input,gene_confirmation,unigene_confirmation = geneandunigeneconfirmation(input)


################################################################################
# do several sequence / annotation property checks
################################################################################
IS_REPETITIVE           = annotatedproteinsequencerepetitivenesscheck(input)
HAS_SMALL_OR_TINY_EXONS = annotatedgeneexonsizeevaluation(input)
HAS_N_SYMBOLS           = ntracksindnasequencecheck(input)
for org in input.keys():
    status,warnings = confirmcanonicalsplicesites(input[org]['genomeseq'],
        input[org]['gff-gene'],exon_fmethod=GFF_CDS_FMETHOD,verbose=False)
    input[org]['warnings'].extend(warnings)


################################################################################
# print Gene Structure information & warnings
################################################################################

for k in input.keys():
    print k, input[k]['proteinfref'], input[k]['gldobj'].gene_start_coordinate(),
    print "start coordinate in provided locus"
    exons = input[k]['gldobj'].as_exons(input=input[k])
    for exon in exons:
        try:
            seq = exon.orf.inputgenomicsequence[exon.start:exon.end].upper()
        except:
            # introns have no ORFS => take first elem, which is an exon...
            seq = exons[0].orf.inputgenomicsequence[exon.start:exon.end].upper()
        atperc = 100*(float(seq.count("A")+seq.count("T"))/len(seq))
        protseq = "" 
        if exon.__class__.__name__ != 'IntronConnectingOrfs' and exon.length <= 30:
            protseq = exon.proteinsequence()
        print "AT=%2.1f%s" % (atperc , "%"), exon, protseq
    if not exons:
        for track in input[k]['gff-gene']: print track
    for warn in input[k]['warnings']:
        print "\t", warn
    print ""

# check if the DNA exon structure equals the given protein sequence (None|True|False)
from dna2prot import dna2protein
print "# test if annotated gene model equals provided protein sequence in AbfgpGeneLocusDirectory"
for k in input.keys():
    exons = input[k]['gldobj'].as_exons()
    seqp = []
    for item in exons:
        if item.__class__.__name__.find('Exon') >= 0:
            seqp.append( item.dnasequence() )
    _trans = dna2protein("".join(seqp))
    _prot  = input[k]['proteinseq'].replace('*','')
    if not _prot:
        print None, "\t", k
        continue
    else:
        print _trans == _prot, "\t", k
    if dna2protein("".join(seqp))!=input[k]['proteinseq'].replace('*',''):
        translated = dna2protein("".join(seqp))
        for offset in range(0,len(input[k]['proteinseq']),100):
            _trans = translated[offset:offset+100]
            _prot  = input[k]['proteinseq'][offset:offset+100].replace("*","")
            print "prot:", input[k]['proteinseq'][offset:offset+100]
            print "gene:", translated[offset:offset+100],
            print _trans == _prot



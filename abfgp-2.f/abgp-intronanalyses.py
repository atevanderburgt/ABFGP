"""
Functions for analyses of introns in (annotated) genestructures 

"""

# import Python things
import os, re
from sets import Set
from sys import exit as sysExit


# YURAY branch point regular expression pattern
BRANCHPOINT_REPATTERN = re.compile("(T|C)T(G|A)A(T|C)",re.IGNORECASE)

# import ABGP functions
from abgp_unigeneconfirmation import geneandunigeneconfirmation
from lib_tcode import obtaintcodedata
from graphAbgp.ordering import order_list_by_attribute
from lib_abgpgff import exons2introns

# Gene Imports
from gene.stop import StopCodon 
from gene.splicesite import _score_splice_site

# Global varibale Import
from settings.gff import GFF_CDS_FMETHOD, GFF_UGEXON_FMETHOD
from settings.executables import TCODE_MAX_NONCODING, TCODE_MIN_CODING
from settings.splicesites import IC_DONOR_PATTERN_OFFSET, IC_ACCEPTOR_PATTERN_OFFSET

def explain():
    """ Explain the columns of the statistics output """
    txt = """# Please use -h or --help for command line options information\n""" +\
    """# Explanation of Intron statistics output columns, tab delimited:     \n""" +\
    """    ProteinIdentifier    str                                                \n""" +\
    """    OrganismIdentifier   str    possibly identical to ProteinIdentifier     \n""" +\
    """    Intron Count         int    number of introns in the annotated gene     \n""" +\
    """    Intron Rank          int    intron rank number                          \n""" +\
    """    UniGeneConfirmed     bool   True, False or None (no data)               \n""" +\
    """    Intron Phase         int    phase (0,1,2) or (F,?) when erroneous       \n""" +\
    """    Inframe Intron?      bool   True when % 3 == 0 and no * in this frame   \n""" +\
    """    Intron Start         int    nt start coordinate (GT....ag)              \n""" +\
    """    Intron End           int    nt end   coordinate (gt....AG)              \n""" +\
    """    Intron Length        int    nt intron length                            \n""" +\
    """    Donor PSSM           float  donor    site PSSM score ~(-4.0 .. 12.0)    \n""" +\
    """    Acceptor PSSM        float  acceptor site PSSM score ~(-4.0 .. 12.0)    \n""" +\
    """    Intron Tcode average float  average EMBOSS:TCODE score ( 0.35 .. 1.35 ) \n""" +\
    """    Intron Tcode symbol  str    0.35 .(N). 0.740 .(?). 0.950 .(C). 135      \n""" +\
    """    * Codon Count        int    number of stop codons in the intron         \n""" +\
    """    * Codon Phases       int    different phase observed among the * codons \n"""

    # return this txt string
    return txt

# end of function explain


def find_stop_codons(seq):
    """
    @type seq:  string
    @param seq: DNA sequence
    """
    pat = re.compile("(tag|taa|tga)", re.IGNORECASE)
    stopcodons = []
    for match in re.finditer(pat,seq):
        stop = StopCodon(match.start())
        stopcodons.append(stop)
    return stopcodons

# end of function find_stop_codons


def find_branch_point(seq, branchrepat = BRANCHPOINT_REPATTERN):
    """
    @type seq:  string
    @param seq: DNA sequence
    """
    branchpoints = []
    for match in re.finditer(branchrepat,seq):
        branchpoints.append( match.start() )
    branchpoints.reverse()
    return branchpoints

# end of function find_branch_point


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
    parser.remove_option('-q')

    # add a special parameter for only printing intron sequences
    parser.add_option("--sequenceonly",
                dest="sequenceonly",
                default=False,
                action="store_true",
                help="output multifasta file of introns sequences in stead of statistics")


    # parse the command line & validate
    (OPTIONS, args) = parser.parse_args()
    validate_genelocustatsoptions(parser,OPTIONS)

    # print explanation of this script & output
    if OPTIONS.EXPLAIN:
       print explain()
       sysExit()

    # change directory to the one specified by optparse option
    os.chdir(OPTIONS.DIRECTORY)

    try:
        locus = AbgpGeneLocusDirectory(OPTIONS.DIRECTORY)
        input = locus.toinputdict()
    except:
        print explain()
        sysExit()
        ## TEMPORARILY backwards-compatibility with old input_data_struct.txt file
        #input = eval(open('input_data_struct.txt').read().strip())

    # do geneandunigeneconfirmation
    input,gene_status,unigene_status = geneandunigeneconfirmation(input,verbose=False)

    if not OPTIONS.sequenceonly:
        # proces tcode data; no needed in `sequenceonly` output modus
        input = obtaintcodedata(input)

    for org in input.keys():

        if OPTIONS.verbose and not input[org]['gff-gene']:
            print "# No gene GFF data present in GeneLocusDirectory"
        if OPTIONS.verbose and not input[org]['gff-unigene']:
            print "# No unigene GFF data present in GeneLocusDirectory"

        # translate gene's exons to introns
        exons = []
        for gfftrack in input[org]['gff-gene']:
            if gfftrack[2] != GFF_CDS_FMETHOD: continue
            exons.append( gfftrack )
        exons.sort()
        introns = exons2introns(exons)

        # translate UNIgene's exons to introns
        unigeneintrons = []
        if input[org]['gff-unigene']:
            unigeneexons = []
            for gfftrack in input[org]['gff-unigene']:
                if gfftrack[2] != GFF_UGEXON_FMETHOD: continue
                unigeneexons.append( gfftrack )
            unigeneexons.sort()
            unigeneintrons = exons2introns(unigeneexons)

        cnt = 0
        anyOrf = input[org]['orfs'].orfs[0]                 
        for intron in introns:
            cnt+=1
            stagff, endgff = intron[3], intron[4]
            sta, end = stagff-1, endgff
            intronseq = input[org]['orfs'].sequence[sta:end]

            if OPTIONS.sequenceonly:
                print ">%s_%s_%s_%s\n%s" % (
                        OPTIONS.DIRECTORY.split("/")[-1], org,
                        cnt, len(introns), intronseq
                        )
                # continue -> no further stats!
                continue

            # score the splice sites
            staD, endD = sta-IC_DONOR_PATTERN_OFFSET[0], sta+2+IC_DONOR_PATTERN_OFFSET[1]
            staA, endA = end-2-IC_ACCEPTOR_PATTERN_OFFSET[0], end+IC_ACCEPTOR_PATTERN_OFFSET[1]
            donorseq = input[org]['orfs'].sequence[staD:endD]
            accepseq = input[org]['orfs'].sequence[staA:endA]
            donorPSSM = _score_splice_site(donorseq,splicetype='donor')
            accepPSSM = _score_splice_site(accepseq,splicetype='acceptor')

            if intronseq[0:2].upper() != "GT": print "Warning::NonCanonicalDonor", intronseq


            # find_branch_point
            branchpoints = find_branch_point(intronseq)
            if branchpoints:
                first_branchpoint = len(intronseq) - branchpoints[0]
            else:
                first_branchpoint = -1

            # calculate GC content
            gc = float(intronseq.upper().count("G") + intronseq.upper().count("C")) / len(intronseq)

            # get tcode data of this intron stretch
            tcode_score = anyOrf.average_tcode_score_of_range(
                    anyOrf._RAW_TCODE_DATA,sta,end)
            if tcode_score >= TCODE_MIN_CODING:
                tcode_symbol = "C"
            elif tcode_score <= TCODE_MAX_NONCODING:
                tcode_symbol = "N"
            else:
                tcode_symbol = "?"

            # is this intron confirmed by an unigene?
            if unigeneintrons and (stagff,endgff) in [ (ugint[3], ugint[4]) for ugint in unigeneintrons ]:
                ugconfirmedintron = True
            elif unigeneintrons:
                ugconfirmedintron = False
            else:
                ugconfirmedintron = None

            # get stop codon data in the intron
            stopcodons = find_stop_codons(intronseq)
            stopcodon_phases = list(Set([ stop.pos % 3 for stop in stopcodons ]))

            # get phase data of the intron
            try:
                orfDid = input[org]['orfid-genestructure'][cnt-1]
                orfAid = input[org]['orfid-genestructure'][cnt]
                orfD   = input[org]['orfs'].get_orf_by_id(orfDid)
                orfA   = input[org]['orfs'].get_orf_by_id(orfAid)
                phaseD = (sta - orfD.frame) % 3
                phaseA = (end - orfA.frame) % 3
                if phaseD == phaseA:
                    intronphase = phaseD
                    if (end-sta) % 3 != 0:
                        inframe = False
                    elif intronphase in stopcodon_phases:
                        inframe = False
                    else:
                        inframe = True
                else:
                    intronphase = "F"
                    if (end-sta) % 3 != 0: inframe = False
                    else:                  inframe = None
            except:
                intronphase = "?"
                if (end-sta) % 3 != 0: inframe = False
                else:                  inframe = None


            print "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%s\t%s\t%s\t%s\t%s" % (
                    OPTIONS.DIRECTORY.split("/")[-1], org, len(introns), cnt, 
                    ugconfirmedintron, intronphase, inframe, 
                    sta, end, (end-sta), gc, donorPSSM, accepPSSM,
                    tcode_score, tcode_symbol, 
                    len(stopcodons), len(stopcodon_phases), 
                    len(branchpoints), first_branchpoint
                    )

    # go back to OPTIONS.CURRENTWORKINGDIR
    os.chdir(OPTIONS.CURRENTWORKINGDIR)

    if OPTIONS.verbose:
        print "# EOF intronstatistics"

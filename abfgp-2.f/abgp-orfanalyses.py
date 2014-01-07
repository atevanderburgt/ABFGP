"""
Functions for analyses of Orfs in (annotated) genestructure vs. all Orfs 
"""

# import Python things
import os
from sys import exit as sysExit

# import ABGP functions
from abgp_geneconfirmation import geneconfirmation
from lib_tcode import obtaintcodedata
from graphAbgp.ordering import order_list_by_attribute
from lib_sequencerepetitiveness import proteinsequencelowcomplexityscore

def _orf_overlap_rank(theorf,orflist,orfpointer):
    """ """
    subset = [ theorf ]
    for orf in orflist[orfpointer+1:]:
        if orf.startPY > theorf.endPY:
            break
        subset.append(orf)
    for pointer in range(orfpointer-1,-1,-1):
        orf = orflist[pointer]
        if orf.endPY < theorf.startPY:
            break
        subset.insert(0,orf)
    # order the subset
    subset = order_list_by_attribute(subset,order_by='coding_propensity',reversed=True)
    rank=0 
    for orf in subset:
        rank+=1
        if orf.startPY == theorf.startPY:
            break
    # return the rank
    return rank

# end of _orf_overlap_rank 


def explain():
    """ Explain the columns of the statistics output """
    txt = """# Please use -h or --help for command line options information\n""" +\
    """# Explanation of Orf statistics output columns, tab delimited:        \n""" +\
    """    ProteinIdentifier    str                                                \n""" +\
    """    OrganismIdentifier	str    possibly identical to ProteinIdentifier     \n""" +\
    """    Rank Nr.             int    rank number of Orf Score                    \n""" +\
    """    Orf Start Coord      int    nt start coordinate (TGAsta...endTGA)       \n""" +\
    """    Orf End   Coord      int    nt end   coordinate (TGAsta...endTGA)       \n""" +\
    """    Orf Score            float  Orf AA Length * Orf Tcode average           \n""" +\
    """    Overlap Rank Nr.     int    Orf Score rank nr. of overlapping Orfs      \n""" +\
    """    Orf AA Length        int                                                \n""" +\
    """    Orf Tcode symbol     str    0.35 .(N). 0.740 .(?). 0.950 .(C). 135      \n""" +\
    """    Orf Tcode average    float  average EMBOSS:TCODE score ( 0.35 .. 1.35 ) \n""" +\
    """    IsCodingOrf          bool   TRUE if included in the annotated gene      \n""" +\
    """    LowComplexity        float  LowComplexity AA score ( 0.0 .. 1.0 )       \n"""
    # return this txt string
    return txt

# end of function explain

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

    # do geneconfirmation; this includes rungetorf() function call
    input,gene_status = geneconfirmation(input,verbose=True)

    # proces tcode data
    input = obtaintcodedata(input)

    for org in input.keys():

        # loop over all the Orf objects and score their coding propensity rank
        orflist = input[org]['orfs'].orfs
        orflist = order_list_by_attribute(orflist,'startPY')
        for orfpointer in range(0,len(orflist)):
            orf = orflist[orfpointer]
            # add data to tss object
            orf.coding_propensity = orf.tcode_score() * orf.length 

        for orfpointer in range(0,len(orflist)):
            orf = orflist[orfpointer]
            # define _orf_overlap_rank
            rank = _orf_overlap_rank(orf,orflist,orfpointer)
            orf._orf_overlap_rank = rank

        # order orflist by coding_propensity and score this rank 
        orflist = order_list_by_attribute(orflist,'coding_propensity',reversed=True)
        cnt = 1
        for orf in orflist:
            orf._rank_coding_propensity = cnt
            cnt+=1

        # print data to stdout
        cnt = 1

        for orf in orflist:
            lsc_score = proteinsequencelowcomplexityscore(orf.protein_sequence)
            print "%s\t%s\t%s\t%s\t%s\t%4.1f\t%s\t%s\t%s\t%1.3f\t%s\t%1.2f" % (
                     OPTIONS.DIRECTORY.split("/")[-1], org, cnt, orf.startPY, orf.endPY, orf.coding_propensity,
                     orf._orf_overlap_rank, orf.length, orf.tcode_symbolic(), orf.tcode_score(),
                     orf.id in input[org]['orfid-genestructure'], lsc_score,
                     )
            cnt+=1

    # go back to OPTIONS.CURRENTWORKINGDIR
    os.chdir(OPTIONS.CURRENTWORKINGDIR)

    if OPTIONS.verbose:
        print "# EOF orfstatistics"


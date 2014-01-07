"""
Functions for obtaining TCODE data for Alignment Based Gene Prediction
"""

# Module metadata
__authors__ = "Ate van der Burgt"
__license__ = "MIT"

# Imports
from parsers.tcode import *
import os

# import global variables
from settings.executables import (
    EXECUTABLE_TCODE,
    EXECUTABLE_TCODE_STEP,
    EXECUTABLE_TCODE_WINDOW,
    EXECUTABLE_TCODE_N_CORRECTION,
    ) 
from settings.abgp import USE_TCODE


def obtaintcodedata(input):
    """
    """
    for org in input.keys():
        if USE_TCODE:
            # make dict key for tcode output
            input[org]['tcode'] = []

            # check the input sequence for too long tracks of continious N's
            (check,correctedseq) = correct_sequence_for_n_tracks(
                        input[org]['genomeseq'],
                        window=EXECUTABLE_TCODE_N_CORRECTION
                        )
            if check == False:
                inputdnaseq = correctedseq
            else:
                inputdnaseq = input[org]['genomeseq']


            # check length of input sequence.
            # UniGenes/ESTs are often to short for TCODE!
            # DNA sequences <= EXECUTABLE_TCODE_WINDOW will
            # give hard errors -> avoid this. Better no TCODE data
            # then algorithm error. TCODE data absence is taken
            # into account furtheron in ABFGP
            if len(inputdnaseq) <= EXECUTABLE_TCODE_WINDOW:
                continue

            # run & parse tcode data 
            input[org]['tcode'] = parsetcodestdout( tcode( 
                        sequence=inputdnaseq,
                        step=EXECUTABLE_TCODE_STEP,
                        window=EXECUTABLE_TCODE_WINDOW,
                        EXECUTABLE_TCODE=EXECUTABLE_TCODE) )

            # TODO; solve this problem!!!!!
            if not input[org]['tcode']:
                print org, " TCODE OUTPUT FILE EMPTY -> TCODE STATUS N-tracks is FALSE!!!!"
                import sys
                sys.exit()

            # (F) store tcode data to orfs
            for orf in input[org]['orfs'].orfs: orf.add_tcode_data(input[org]['tcode'])

    # return the incremented input dict
    return input

# end of function obtaintcodedata



def intergenecity_tcode_analyses(prevcbg,nextcbg,input):
    """
    """
    prevomsr = prevcbg.overall_minimal_spanning_range()
    nextomsr = nextcbg.overall_minimal_spanning_range()
    tcode_nt_window = 100
    tcode_nt_step   = 20
    tcodedata = {}
    for org in prevcbg.organism_set().intersection(nextcbg.organism_set()):
        prevnode = prevcbg.node_by_organism(org)
        nextnode = nextcbg.node_by_organism(org)
        stacoord = max(prevomsr[prevnode])*3
        endcoord = min(nextomsr[nextnode])*3
        tclength = endcoord-stacoord
        # calculate average TCODE of complete window
        if tclength >= 1:
            tcode = sum([ elem[2] for elem in input[org]['tcode'][stacoord:endcoord] ])
            avtcode = tcode / tclength
        else:
            # no TCODE calcultation possible!
            continue

        # calculate lowest TCODE score of tcode_nt_window with tcode_nt_step in window
        tcodedata[org] = [ avtcode ]
        tmpscores = []
        for offset in range(0,tclength-tcode_nt_window+1,tcode_nt_step):
            sta = stacoord + offset
            end = stacoord + offset + tcode_nt_window
            tcode = sum([ elem[2] for elem in input[org]['tcode'][sta:end] ])
            tmpscores.append( tcode / (end-sta) )
        if not tmpscores:
            tcodedata[org] = [ avtcode ] * 4
        else:
            tcodedata[org].append( min(tmpscores) )
            tcodedata[org].append( sum(tmpscores)/len(tmpscores) )
            tcodedata[org].append( max(tmpscores) )

    # done!
    return tcodedata
            
# end of function intergenecity_tcode_analyses



def tcode_analyses_upstream_of_firstcbg(GSG,input,
    offsets = [700,600,500,400,300,200,100,0],
    stepsize=100, verbose=False
    ):
    """
    Obtain window-based TCODE data upstream of firstCBG

    @type  GSG: GenestructureOfCodingBlockGraphs
    @param GSG: GenestructureOfCodingBlockGraphs

    @type  input: dict
    @param input: input dictionary data structure

    @type  offsets: list
    @param offsets: (ordered) list of offset to calculate TCODE data for

    @type  stepsize: int
    @param stepsize: positive(!) integer; window/stepsize for TCODE calculation

    @type  verbose: Boolean
    @param verbose: print result table to STDOUT if True

    @rtype:  dict, list
    @return: dictionary with tcode data (keys==organisms), list with offsets used
    """
    firstcbg = GSG.get_first_cbg()
    # if no firstCBG available -> return empty data
    if not firstcbg: return {}, []
    firstomsr = firstcbg.overall_minimal_spanning_range()
    tcodedata = {}
    for node in firstcbg.get_ordered_nodes():
        tcodedata[firstcbg.organism_by_node(node)] = []
    for offset in offsets:
        for node in firstcbg.get_ordered_nodes():
            org   = firstcbg.organism_by_node(node)
            sta   = min(firstomsr[node])*3 - offset
            tcode = sum([ elem[2] for elem in input[org]['tcode'][sta:sta+stepsize] ]) / stepsize
            tcodedata[org].append(tcode)

    ############################################################################
    # in verbose node, print results to stdout in a nice table
    # example:
    # ## TCODE upstream of first CBG
    # offsets:   700     600     500     400     300     200     100     0 
    # AVERAGE    0.710   0.729   0.787   0.806   0.731   0.826   0.895   0.919 
    # fveg       0.600   0.806   0.836   0.809   0.765   0.792   0.938   0.812 
    # mgg        0.646   0.723   0.659   0.807   0.883   0.743   0.968   0.996 
    # ncu        0.788   0.721   0.721   0.777   0.648   0.870   0.807   0.870 
    # foxg       0.797   0.651   0.901   0.922   0.703   0.803   0.811   0.885 
    # fgsg       0.720   0.745   0.817   0.713   0.657   0.921   0.950   1.030 
    ############################################################################
    if verbose:
        # print results to stdout in a nice table
        print "## TCODE upstream of first CBG"
        print "offsets:",
        for offset in offsets:
            print "\t",offset,
        print ""
        print "AVERAGE",
        for i in range(0,len(offsets)):
            print "\t%1.2f" % (
                sum([vlist[i] for vlist in tcodedata.values()])/len(tcodedata)),
        print ""
        for org, vlist in tcodedata.iteritems():
            print org, "\t",
            for i in range(0,len(offsets)):
                print "\t%1.2f" % vlist[i],
            print ""
        print ""
    ############################################################################

    # return tcodedata and offsets used
    return tcodedata, offsets

# end of function tcode_analyses_upstream_of_firstcbg



def tcode_analyses_downstream_of_finalcbg(GSG,input,
    offsets = [0,100,200,300,400,500,600,700],
    stepsize=100, verbose=False
    ):
    """
    Obtain window-based TCODE data downstream of finalCBG

    @type  GSG: GenestructureOfCodingBlockGraphs
    @param GSG: GenestructureOfCodingBlockGraphs

    @type  input: dict
    @param input: input dictionary data structure

    @type  offsets: list
    @param offsets: (ordered) list of offset to calculate TCODE data for

    @type  stepsize: int
    @param stepsize: positive(!) integer; window/stepsize for TCODE calculation

    @type  verbose: Boolean
    @param verbose: print result table to STDOUT if True

    @rtype:  dict, list
    @return: dictionary with tcode data (keys==organisms), list with offsets used
    """
    finalcbg = GSG.get_final_cbg()
    # if no finalCBG available -> return empty data
    if not finalcbg: return {}, []
    finalomsr = finalcbg.overall_minimal_spanning_range()
    tcodedata = {}
    for node in finalcbg.get_ordered_nodes():
        tcodedata[finalcbg.organism_by_node(node)] = []
    for offset in offsets:
        for node in finalcbg.get_ordered_nodes():
            org   = finalcbg.organism_by_node(node)
            sta   = max(finalomsr[node])*3 + offset
            tcode = sum([ elem[2] for elem in input[org]['tcode'][sta:sta+stepsize] ]) / stepsize
            tcodedata[org].append(tcode)

    ############################################################################
    # in verbose node, print results to stdout in a nice table
    # example:
    # ## TCODE downstream of final CBG
    # offsets:   0       100     200     300     400     500     600     700 
    # AVERAGE    0.831   0.763   0.608   0.597   0.774   0.732   0.711   0.804 
    # fveg       0.961   0.885   0.533   0.537   0.810   0.816   0.746   0.783 
    # mgg        0.963   0.805   0.640   0.569   0.577   0.543   0.593   0.627 
    # ncu        0.815   0.636   0.595   0.670   0.749   0.617   0.667   0.896 
    # foxg       0.845   0.635   0.540   0.641   0.894   0.793   0.766   0.813 
    # fgsg       0.573   0.856   0.731   0.568   0.841   0.892   0.782   0.901 
    ############################################################################
    if verbose:
        # print results to stdout in a nice table
        print "## TCODE downstream of final CBG"
        print "offsets:",
        for offset in offsets:
            print "\t",offset,
        print ""
        print "AVERAGE",
        for i in range(0,len(offsets)):
            print "\t%1.2f" % (
                sum([vlist[i] for vlist in tcodedata.values()])/len(tcodedata)),
        print ""
        for org, vlist in tcodedata.iteritems():
            print org, "\t",
            for i in range(0,len(offsets)):
                print "\t%1.2f" % vlist[i],
            print ""
        print ""
    ############################################################################

    # return tcodedata and offsets used
    return tcodedata, offsets

# end of function tcode_analyses_downstream_of_finalcbg

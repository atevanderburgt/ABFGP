"""
Functions for printing analytic statistics of ABGP graphs
"""

from exceptions import (
    NoCodingBlockGraphApplied,
    NoGenestructureOfCodingBlockGraphsApplied,
    NoGeneTreeGraphApplied
    )


def average(values):
    """ """
    return float(sum(values))/float(len(values))
# end of function average


def printcbganalyses(cbg,GSG=None,GTG=None,prefix=None):
    """
    Print tab delimited line of CodingBlockGraph characteristics

    @type  cbg: CodingBlockGraph
    @param cbg: CodingBlockGraph
    
    @type  GSG: GenestructureOfCodingBlockGraphs
    @param GSG: GenestructureOfCodingBlockGraphs in which this CBG is from

    @type  GTG: GeneTreeGraph
    @param GTG: GeneTreeGraph of the GenestructureOfCodingBlockGraphs

    @type  prefix: * (string)
    @param prefix: first item to print (e.g. identifier, header, directory)

    @attention: elem1  = str prefix
    @attention: elem2  = TF  is this CBG accepted in the GSG or not?
    @attention: elem3  = int cbg.total_weight()
    @attention: elem4  = 0-1 gtg.identity() of this CBG
    @attention: elem5  = 0-1 GTG.identity() of the GTG of the GSG
    @attention: elem6  = 0>1 cbg.omsr_tcode_score(), TCODE score of this CBG
    @attention: elem7  = int cbg.omsr_length(), AA length of this CBG
    @attention: elem8  = 0>1 gtg / GTG identity ratio
    @attention: elem9  = 0>1 gtg / GTG relative graph alignment difference
    @attention: elem10 = 0>1 gtg / GTG absolute graph alignment difference
    @attention: elem11 = str (Organism identifier ordered) Orf identifier list
    """
    # input argument validation
    if not cbg or cbg.__class__.__name__ != 'CodingBlockGraph':
        raise NoCodingBlockGraphApplied
    if not GSG or GSG.__class__.__name__ != 'GenestructureOfCodingBlockGraphs':
        raise NoGenestructureOfCodingBlockGraphsApplied
    if not GTG or GTG.__class__.__name__ != 'GeneTreeGraph':
        raise NoGeneTreeGraphApplied

    # get the GeneTreeGraph of this CBG; can be identical to GTG !
    gtg = cbg.genetree()
    # get list (line) with statistics
    line = [
            prefix,
            GSG.is_codingblockgraph_already_in_genestructure(cbg),
            cbg.total_weight(),
            "%1.2f" % gtg.identity(),
            "%1.2f" % GTG.identity(),
            "%1.2f" % cbg.omsr_tcode_score(),
            cbg.omsr_length(),
            "%1.3f" % ( gtg.identity() / GTG.identity() ),
            "%1.3f" % GTG.graphalignmentdifference( gtg ),
            "%1.3f" % GTG.absolutegraphalignmentdifference( gtg ),
            # ordered Orf identifier list
            [ node[1] for node in cbg.get_ordered_nodes()],
            ]
    
    # print line (to STDOUT)
    print "\t".join( [ str(elem) for elem in line ] )

# end of function printcbganalyses


def printgtganalyses(gtg,prefix=None):
    """
    Print tab delimited line of GeneTreeGraph characteristics

    @type  gtg: GeneTreeGraph
    @param gtg: GeneTreeGraph of the GenestructureOfCodingBlockGraphs

    @type  prefix: * (string)
    @param prefix: first item to print (e.g. identifier, header, directory)

    @attention: elem1  = str prefix
    @attention: elem2  = 0-1 gtg.identity()   overall AA identity + similarity
    @attention: elem3  = 0-1 gtg.aaidentity() overall AA identity 
    @attention: elem4  = 0-1 gtg.ntidentity() overall NT identity
    @attention: elem5  = 0-1 gtg.bitscoreratio() overall bitscoreratio
    @attention: elem6  = 0-1 highest AA identityscore edge in the GTG
    @attention: elem7  = 0-1 lowest  AA identityscore edge in the GTG
    @attention: elem8  = str strongest connected node in the GTG
    @attention: elem9  = 0>1 strongest connected node, relative connectivity
    @attention: elem10 = str   weakest connected node in the GTG
    @attention: elem11 = 0-1   weakest connected node, relative connectivity
    @attention: elem12 = str string representation of the GTG
    """
    # input argument validation
    if not gtg or gtg.__class__.__name__ != 'GeneTreeGraph':
        raise NoGeneTreeGraphApplied
    
    weakestnode   = gtg.weakest_connected_node()
    strongestnode = gtg.strongest_connected_node()
    sc1 = gtg.get_node_weighted_connectivity_observed_vs_expected(strongestnode)
    sc2 = gtg.get_node_weighted_connectivity_observed_vs_expected(weakestnode)

    # get list (line) with statistics
    line = [
            prefix,
            "%1.2f" % gtg.identity(),
            "%1.2f" % gtg.aaidentity(),
            "%1.2f" % gtg.ntidentity(),
            "%1.2f" % gtg.bitscoreratio(),
            "%1.2f" % max(gtg.weights.values()),
            "%1.2f" % min(gtg.weights.values()),
            strongestnode,
            "%1.2f" % sc1,
            weakestnode,
            "%1.2f" % sc2,
            str(gtg),
           ]

    # print line (to STDOUT)
    print "\t".join( [ str(elem) for elem in line ] )
    
# end of function printgtganalyses


def printstopcodonanalyses(lastcbg,GTG=None,prefix=None):
    """
    Print tab delimited line of stopcodon(graph) of CodingBlockGraph

    @type  lastcbg: CodingBlockGraph
    @param lastcbg: CodingBlockGraph, presumably the correct final CBG
    
    @type  GTG: GeneTreeGraph
    @param GTG: GeneTreeGraph of the GenestructureOfCodingBlockGraphs

    @type  prefix: * (string)
    @param prefix: first item to print (e.g. identifier, header, directory)

    @attention: elem1  = str prefix
    @attention: elem2  = TFN ASCG is_optimal() overall outcome
    @attention: elem3  = TFN ASCG is_optimal_weight() stops not to far apart?
    @attention: elem4  = TFN ASCG is_optimal_tcode_ratio() is TCODE entropy oke?
    @attention: elem5  = TFN ASCG is_optimal_gtgweakestnode() as in the GTG?
    @attention: elem6  = 0-1 identity of GTG of GSG
    @attention: elem7  = 0-1 connectivity of weakestnode in the GTG
    @attention: elem8  = 0-1 identity of gtg of this finalCBG
    @attention: elem9  = 0-1 connectivity of weakestnode in this finalCBG
    @attention: elem10 = 0-1 ASCG average_weight() measure for stops distance
    @attention: elem11 = 0-1 ASCG weakest edge weight
    @attention: elem12 = 0-1 ASCG average edge weight
    @attention: elem13 = 0-1 ASCG strongest edge weight
    @attention: elem14 = 0>1 ASCG TCODE entropy 5p / 3p side
    @attention: elem15 = 0>1 ASCG TCODE 5p minimal value
    @attention: elem16 = 0>1 ASCG TCODE 5p average value
    @attention: elem17 = 0>1 ASCG TCODE 5p maximal value
    @attention: elem18 = 0>1 ASCG TCODE 3p minimal value
    @attention: elem19 = 0>1 ASCG TCODE 3p average value
    """
    # input argument validation
    if not lastcbg or lastcbg.__class__.__name__ != 'CodingBlockGraph':
        raise NoCodingBlockGraphApplied
    if not GTG or GTG.__class__.__name__ != 'GeneTreeGraph':
        raise NoGeneTreeGraphApplied

    # prepare data
    lastcbg.align_stop_codons()
    ascg = lastcbg._stopcodongraph
    gtg  = lastcbg.genetree()
    GTGweakest = GTG.weakest_connected_node()
    gtgweakest = gtg.weakest_connected_node()
    sc1 = GTG.get_node_weighted_connectivity_observed_vs_expected(GTGweakest)
    sc2 = gtg.get_node_weighted_connectivity_observed_vs_expected(gtgweakest)

    # get list (line) with statistics
    line = [
            prefix,
            "%s"    % ascg.is_optimal(),
            "%s"    % ascg.is_optimal_weight(),
            "%s"    % ascg.is_optimal_tcode_ratio(),
            "%s"    % ascg.is_optimal_gtgweakestnode(),
            "%1.2f" % GTG.identity(),
            "%1.2f" % sc1,
            "%1.2f" % gtg.identity(),
            "%1.2f" % sc2,
            ###"%2.1f" % ascg.total_weight(),
            "%1.2f" % ascg.average_weight(),
            "%1.2f" % min(ascg.weights.values()),
            "%1.2f" % average(ascg.weights.values()),
            "%1.2f" % max(ascg.weights.values()),
            "%1.2f" % ascg.average_tcode_entropy(),
            "%1.2f" % min(ascg._tcode5pscore.values()),
            "%1.2f" % average(ascg._tcode5pscore.values()),
            "%1.2f" % max(ascg._tcode5pscore.values()),
            "%1.2f" % min(ascg._tcode3pscore.values()),
            "%1.2f" % average(ascg._tcode3pscore.values()),
            "%1.2f" % max(ascg._tcode3pscore.values()),
            ]

    # print line (to STDOUT)
    print "\t".join( [ str(elem) for elem in line ] )

# end of function printstopcodonanalyses


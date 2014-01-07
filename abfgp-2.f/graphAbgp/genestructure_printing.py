"""
Class with functions for printing the statistics of the CBGs
in a GenestructureOfCodingBlockGraphs used in Aligment Based Gene Predictions
"""

# Module metadata
__authors__ = "Ate van der Burgt"
__license__ = "MIT"


# Python Imports
from copy import deepcopy

class GenestructureOfCodingBlockGraphsPrintingFunctions:
    """
    """
    def printGSGstructure(self,ignore_ignored=True):
        """
        Print CBG topology report for the GSG
        """
        print "########"*(len(self)+1)

        print "struc\t",
        for cbg in self.codingblockgraphs:
            if cbg.IS_FIRST and cbg.IS_LAST:
                print "FI/LA\t",
            elif cbg.IS_FIRST:
                print "FIRST\t",
            elif cbg.IS_LAST:
                print "LAST\t",
            elif ignore_ignored and cbg.IS_IGNORED:
                pass
            elif not ignore_ignored and cbg.IS_IGNORED:
                print "ignor\t",
            elif cbg.__class__.__name__ == 'LowSimilarityRegionCodingBlockGraph':
                print "-lsr-\t",
            else:
                print "- \t",
        print ""

        # print information for the aligned stop-codon graph
        print "*gra\t",
        for cbg in self.codingblockgraphs:
            if ignore_ignored and cbg.IS_IGNORED: continue
            if not cbg._stopcodongraph:
                cbg.align_stop_codons()
            if cbg._stopcodongraph:
                print "%1.2f\t" % cbg._stopcodongraph.average_weight(),
            else:
                print "n.a.\t",
        print ""


        # print information on have_all_starts_upstream_of_omsr
        print "TSS\t",
        for cbg in self.codingblockgraphs:
            if ignore_ignored and cbg.IS_IGNORED: continue
            print "%s\t" % cbg.have_all_starts_upstream_of_omsr(),
        print ""


        # print information on the edges in the CBGs
        print "edges\t",
        for cbg in self.codingblockgraphs:
            if ignore_ignored and cbg.IS_IGNORED: continue
            if cbg.__class__.__name__ == 'LowSimilarityRegionCodingBlockGraph':
                print "-lsr-\t",
            else:
                print "%1.1f\t" % cbg.connectivitysaturation(),
        print ""

        # print information on the PacbP(ORFs) in the CBGs
        print "PACBPS\t",
        for cbg in self.codingblockgraphs:
            if ignore_ignored and cbg.IS_IGNORED: continue
            if cbg.__class__.__name__ == 'LowSimilarityRegionCodingBlockGraph':
                print "-lsr-\t",
            else:
                print "%s-%s\t" % ( str(cbg.has_all_pacbps())[0], len(cbg.pacbps) ),
        print ""




        print "split\t",
        for cbg in self.codingblockgraphs:
            s1,s2,s3 = cbg.IS_SPLITTED, cbg.IS_5P_SPLITTED, cbg.IS_3P_SPLITTED
            if ignore_ignored and cbg.IS_IGNORED:
                continue
            elif (s1,s2,s3) == (True,True,True):
                print "3p 5p\t",
            elif (s1,s2,s3) == (True,True,False):
                print "5p\t",
            elif (s1,s2,s3) == (True,False,True):
                print "3p\t",
            elif (s1,s2,s3) == (False,False,False):
                print "- \t",
            else:
                print "FALSE\t",
        print ""

        print "cbgIF\t",
        for i in range(0,len(self)):
            printstring = ""
            if i==0: printstring += "na"
            else:
                cbg = self.codingblockgraphs[i]
                if ignore_ignored and cbg.IS_IGNORED: continue
                if self.has_acceptor_cbginterface(cbg):
                    if cbg.IS_5P_SPLITTED:
                        printstring += "<"
                    elif cbg._CBGinterface5p._optimal_aligned_acceptor:
                        phase = cbg._CBGinterface5p._optimal_aligned_acceptor.phase()
                        clnm  = cbg._CBGinterface5p._optimal_aligned_acceptor.__class__.__name__  
                        if phase in [0,1,2]:
                            printstring += str(phase)
                        elif clnm == 'AlignedAcceptorSiteWithPhaseShiftGraph':
                            printstring += "P"
                        else:
                            printstring += "?"
                    else:
                        printstring += "."
                else:
                    printstring += "."
                if self.cbginterface_is_optimal_acceptor(cbg):
                    printstring += "+"
                else:
                    printstring += "-"
            # append space
            printstring += " "
            if i==len(self)-1: printstring += "na"
            else:
                cbg = self.codingblockgraphs[i]
                if self.cbginterface_is_optimal_donor(cbg):
                    printstring += "+"
                else:
                    printstring += "-"
                if self.has_donor_cbginterface(cbg):
                    if cbg.IS_3P_SPLITTED:
                        printstring += ">"
                    elif cbg._CBGinterface3p._optimal_aligned_donor:
                        phase = cbg._CBGinterface3p._optimal_aligned_donor.phase()
                        clnm  = cbg._CBGinterface3p._optimal_aligned_donor.__class__.__name__
                        if phase in [0,1,2]:
                            printstring += str(phase)
                        elif clnm == 'AlignedDonorSiteWithPhaseShiftGraph':
                            printstring += "P"
                        else:
                            printstring += "?"
                    else:
                        printstring += "."
                else:
                    printstring += "."
            # print this generated string
            print printstring+"\t",
        print ""

        # add line for weather or not the CBG is optimal
        print "OPTIM\t",
        for cbg in self:
            statuslist = [ self._codingblock_prediction_status(cbg,org) for org in cbg.organism_set() ] 
            if False in statuslist:      print "False\t", 
            elif not True in statuslist: print "None\t",
            else:                        print "True\t",
        print ""

    # end of function printGSGstructure


    def printGTGanalyses(self,ignore_ignored=True):
        """
        Print CBG topology report for the GSG
        """
        GTG = self.genetree()
        print "########"*(len(self)+1)
        print "weight\t",
        for cbg in self.codingblockgraphs:
            if ignore_ignored and cbg.IS_IGNORED: continue
            print "%s\t" % cbg.total_weight(),
        print ""

        print "length\t",
        for cbg in self.codingblockgraphs:
            if ignore_ignored and cbg.IS_IGNORED: continue
            if cbg.__class__.__name__ == 'LowSimilarityRegionCodingBlockGraph':
                print "%s-%s\t" % (
                    min( cbg.minimal_spanning_range_sizes().values() ),
                    max( cbg.minimal_spanning_range_sizes().values() ),
                    ),
            else:
                print "%s \t" % min( cbg.minimal_spanning_range_sizes().values() ),
        print ""

        print "TCODE\t",
        for cbg in self:
            if ignore_ignored and cbg.IS_IGNORED: continue
            if cbg.__class__.__name__ == 'LowSimilarityRegionCodingBlockGraph':
                print "-lsr-\t",
            else:
                print "%1.3f\t" % cbg.msr_tcode_score(),
        print ""


        print "CXPDR\t",
        for cbg in self:
            if ignore_ignored and cbg.IS_IGNORED: continue
            if cbg.__class__.__name__ == 'LowSimilarityRegionCodingBlockGraph':
                print " --- \t",
            else:
                if not cbg._cexpander: cbg.cexpanderanalyses(projected_on=":::")
                cexpstring = cbg._cexpander.binarystring
                if len(cexpstring) == 0:
                    # no cexpander string !?
                    print "-----\t",
                elif cexpstring.count("1") == 0:
                    print "00000\t",
                elif cexpstring.count("0") == 0:
                    print "11111\t",
                elif cexpstring[0] == "0" and cexpstring[-1] == "0":
                    print "01110\t",
                elif cexpstring[0] == "0":
                    print "00111\t",
                elif cexpstring[-1] == "0":
                    print "11100\t",
                else:
                    print "11011\t",
        print ""

        print "CXPDR\t",
        for cbg in self:
            if ignore_ignored and cbg.IS_IGNORED: continue
            if cbg.__class__.__name__ == 'LowSimilarityRegionCodingBlockGraph':
                print " --- \t",
            else:
                if not cbg._cexpander: cbg.cexpanderanalyses(projected_on=":::")
                ratio = cbg._cexpander.uniformly_matched_ratio()
                if ratio == None:
                    # no cexpander binarystring !?
                    print "-----\t",
                else:
                    print "%1.2f\t" % ratio,
        print ""



        print "%ID\t",
        for cbg in self.codingblockgraphs:
            if ignore_ignored and cbg.IS_IGNORED: continue
            if cbg.__class__.__name__ == 'LowSimilarityRegionCodingBlockGraph':
                print "-lsr-\t",
            else:
                print "%1.3f\t" % cbg.genetree().identity(),
        print ""

        print "IDrat\t",
        for cbg in self.codingblockgraphs:
            if ignore_ignored and cbg.IS_IGNORED: continue
            if cbg.__class__.__name__ == 'LowSimilarityRegionCodingBlockGraph':
                print "-lsr-\t",
            else:
                print "%1.3f\t" % (cbg.genetree().identity() / GTG.identity()),
        print ""

        print "TOPdif\t",
        for cbg in self.codingblockgraphs:
            if ignore_ignored and cbg.IS_IGNORED: continue
            if cbg.__class__.__name__ == 'LowSimilarityRegionCodingBlockGraph':
                print "-lsr-\t",
            elif cbg.node_count() == self.EXACT_SG_NODE_COUNT:
                print "%1.3f\t" % ( GTG.graphalignmentdifference( cbg.genetree() ) ),
            else:
                GTGdelnode = deepcopy(GTG)
                for missingorg in GTG.organism_set().difference(cbg.organism_set()):
                    GTGdelnode.del_node(missingorg)
                print "%1.3f\t" % ( GTGdelnode.graphalignmentdifference( cbg.genetree() ) ),
        print ""

        print "ABSdif\t",
        for cbg in self.codingblockgraphs:
            if ignore_ignored and cbg.IS_IGNORED: continue
            if cbg.__class__.__name__ == 'LowSimilarityRegionCodingBlockGraph':
                print "-lsr-\t",
            elif cbg.node_count() == self.EXACT_SG_NODE_COUNT:
                print "%1.3f\t" % ( GTG.absolutegraphalignmentdifference( cbg.genetree() ) ),
            else:
                GTGdelnode = deepcopy(GTG)
                for missingorg in GTG.organism_set().difference(cbg.organism_set()):
                    GTGdelnode.del_node(missingorg)
                print "%1.3f\t" % ( GTGdelnode.absolutegraphalignmentdifference( cbg.genetree() ) ),
        print ""

    # end of function printGTGanalyses

# end of class GenestructureOfCodingBlockGraphsPrintingFunctions

 

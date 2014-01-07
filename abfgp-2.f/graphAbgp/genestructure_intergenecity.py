"""
Class with functions that can detect and resolve microsyntheny
in a GeneStructureOfCodingBlocks used in Alignment Based Gene Prediction
"""

# Module metadata
__authors__ = "Ate van der Burgt"
__license__ = "MIT"

# ABGP Imports
from lib_tcode import intergenecity_tcode_analyses
import recombination


# Global Variable Imports
from settings.executables import (
    TCODE_MAX_NONCODING,
    TCODE_MIN_CODING,
)
from settings.genestructure import (
    ABSOLUTE_MAX_INTRON_NT_LENGTH,
    MAX_INTRON_NT_LENGTH,
    MAX_INTERGENIC_MIN_NT_LENGTH,
    AVERAGE_INTERGENIC_MIN_NT_LENGTH,
    MIN_INTERGENIC_NT_LENGTH,
) 



INTERGENIC_TCODE_NT_WINDOW_SIZE     = 100
INTERGENIC_TCODE_NT_STEP_SIZE       = 20
INTERGENIC_TCODE_MAX_LOWEST_WINDOW  = TCODE_MAX_NONCODING - 0.05
INTERGENIC_TCODE_MAX_AVERAGE_WINDOW = TCODE_MIN_CODING
INTERGENIC_TCODE_MIN_STGRA_AV_WT    = 0.6

class IntergenecityFunctions:

    def resolve_microsyntheny(self,
        absolute_max_intron_nt_length=ABSOLUTE_MAX_INTRON_NT_LENGTH,
        max_intron_nt_length=MAX_INTRON_NT_LENGTH,
        max_intergenic_min_nt_length=MAX_INTERGENIC_MIN_NT_LENGTH,
        average_intergenic_min_nt_length=AVERAGE_INTERGENIC_MIN_NT_LENGTH,
        min_intergenic_nt_length=MIN_INTERGENIC_NT_LENGTH,
        apply_tcode_check=False,
        verbose=False):
        """ 
        Assign codingblocks outside of *the* genestructure that belong to another gene (in a synthenic region)
            
        @type  max_intron_nt_length: integer (positive)
        @param max_intron_nt_length: maximal intron size -> any larger distance MUST be CBGs of 2 distinct genes
            
        @type  min_intergenic_nt_length: integer (positive)
        @param min_intergenic_nt_length: minimum intergenic size -> any distance smaller MUST be CBGs of a single gene

        @type  max_intergenic_min_nt_distance: integer (positive)
        @param max_intergenic_min_nt_distance: represent the upper border of a small intergenic region
                
                
        @attention: *the* genestructure is supposed to be the most central one in case of doubt
        """     
                
        # List of list of position id's of CBGs together are (part of a) gene
        # When a distance is large enough for an intergenic region, the structure will
        # be something like [ [ 0,1,2], [3,4,5,6,7] ]. The first 3 CBGs are assigned to
        # another, frontal placed gene.
        structures = []
        for pos in range(0,len(self)-1):
            this = self.codingblockgraphs[pos]
            next = self.codingblockgraphs[pos+1]
            (mindist,maxdist,avdist) = _cbg_intergenic_distance_analyses(this,next)
            # check which nodes are mutual
            mutual    = this.mutual_nodes(next)

            ########################################################
            if verbose:
                print pos, pos+1,
                print "%s - %s - %s" % (mindist, avdist, maxdist),
                print mutual, "cbgIFopt?:",
                try:
                    print this._CBGinterface3p.optimalitycheck()
                except:
                    print "NO cbgIF!!"
            ########################################################

            if len(mutual) >= 1:
                # orf node(s) are shared -> this MUST be the same gene!
                pass
            elif maxdist > absolute_max_intron_nt_length:
                # largest distance is larger than the largest possible intron
                if not structures:  structures.append( range( 0, pos+1) )
                else:               structures.append( range( structures[-1][-1]+1, pos+1) )
            elif maxdist > max_intron_nt_length:
                # largest distance is larger than what is likely a very large intron
                # in this case: check the cbgIF
                if this._CBGinterface3p and this._CBGinterface3p.is_optimal() and\
                not next.have_all_starts_upstream_of_omsr():
                    # interface is optimal and 2th CBG is not a first TSS CBG.
                    # do not break the genestructure here.
                    pass
                else:
                    # break the genestructure in 2 parts
                    if not structures:  structures.append( range( 0, pos+1) )
                    else:               structures.append( range( structures[-1][-1]+1, pos+1) )
            elif mindist > min_intergenic_nt_length and maxdist > max_intergenic_min_nt_length:
                # the ``grey zone``, but very likely to be intergenic
                if this._CBGinterface3p and this._CBGinterface3p.is_compatible() and\
                this._CBGinterface3p.optimalitycheck().count(True) >= 2 and\
                not next.have_all_starts_upstream_of_omsr():
                    pass
                else:
                    if not structures:  structures.append( range( 0, pos+1) )
                    else:               structures.append( range( structures[-1][-1]+1, pos+1) )
            elif avdist > average_intergenic_min_nt_length:
                # very likely intergenic!
                if this._CBGinterface3p and this._CBGinterface3p.is_compatible() and\
                this._CBGinterface3p.optimalitycheck().count(True) >= 2 and\
                not next.have_all_starts_upstream_of_omsr():
                    pass
                else:
                    if not structures:  structures.append( range( 0, pos+1) )
                    else:               structures.append( range( structures[-1][-1]+1, pos+1) )
            elif maxdist <= min_intergenic_nt_length:
                # largest distance is smaller than ``smallest possible`` intergenic region
                pass
            else:
                # all other cases...many of these are spurious intergenic
                # but can as well be putative long introns. Do, if requested
                # for, the tcode_check. If not-> pass (and keep in the GSG)
                if apply_tcode_check:
                    # calculate TCODE averages
                    tcodescore = intergenecity_tcode_analyses(this,next,self.input)
                    averages = [ 0.0, 0.0, 0.0, 0.0 ]
                    for org,data in tcodescore.iteritems():
                        for i in range(0,len(data)): averages[i]+= data[i]
                    averages = [ item/len(tcodescore.keys()) for item in averages ]
                    omsrdists = this.omsr_distance_between_codingblocks(next)
                    lengthvar = [ float(a)/float(b) for a,b in recombination.pairwise(omsrdists.values()) ]
                    lengthvar = sum(lengthvar)/len(lengthvar)
                    ########################################################
                    if verbose:
                        for org,data in tcodescore.iteritems():
                            print org, "\t", data
                        print "AVER:\t", averages, "*graAV:", 
                        print this._stopcodongraph.average_weight(),
                        print "lenghtvar:", lengthvar
                        print this._CBGinterface3p
                    ########################################################

                    if averages[1] < INTERGENIC_TCODE_MAX_LOWEST_WINDOW and\
                    averages[2] < INTERGENIC_TCODE_MAX_AVERAGE_WINDOW and\
                    this._stopcodongraph.average_weight() >=\
                    INTERGENIC_TCODE_MIN_STGRA_AV_WT and\
                    (not this._CBGinterface3p or not this._CBGinterface3p.is_compatible()):
                        # yes, a somewhat less obvious intergenic boundary
                        if not structures:  structures.append( range( 0, pos+1) )
                        else:               structures.append( range( structures[-1][-1]+1, pos+1) )
                    else:
                       pass 
                else:
                    pass
        else:
            if not structures:  structures.append( range( 0, len(self) ) )
            else:               structures.append( range( structures[-1][-1]+1, len(self) ) )


        # now deal with the results
        if len(structures) > 1:
            # hmmmm... most likely a - rare(!?) - region of microsyntheny
            # in all species. Get the best structure out of these
            best_structure = self._define_best_structure(structures,verbose=verbose)
            self._update_best_structure_in_gsg(best_structure)
            ####################################################
            if verbose:
                print "microsyntheny", structures,
                print "MAX-intron:", max_intron_nt_length
                print "BEST:", best_structure
            ####################################################

    # end of function resolve_microsyntheny


    def _update_best_structure_in_gsg(self,best_structure):
        """ Update the results from the intergenecity analyses in the GSG """
        # set all CBGs outside of best_structures as IS_IGNORED
        for pos in range(0,len(self)):
            cbg = self.codingblockgraphs[pos]
            if pos not in best_structure:
                cbg.IS_IGNORED = True
                cbg.IS_FIRST   = False
                cbg.IS_LAST    = False

        if best_structure:
            # set the IS_FIRST and IS_LAST attribute
            self.codingblockgraphs[best_structure[0]].IS_FIRST = True
            self.codingblockgraphs[best_structure[-1]].IS_LAST = True

    # end of function _update_best_structure_in_gsg


    def _define_best_structure(self,structures,verbose=False):
        """
        Select **the** best structure from several intergenecity separated structures

        @type  structures: list
        @param structures: list with lists of integers, representing CBG positions in the GSG

        @rtype:  list
        @return: exactly ONE of the substructures in the list of structures
        """
        # return current single structure if there is only 1
        if len(structures) == 1: return structures[0]

        ########################################################
        if verbose:
            print "BEST STRUCTURE SELECTION:"
            print structures
            for cbg in self.codingblockgraphs: print cbg
        ########################################################

        # get input DNA sequence length to measure the offset on
        org2lengths = {}
        for org in self.input.keys():
            org2lengths[org] = len(self.input[org]['genomeseq'])
        # calculate average distances towards start & end of DNA sequence
        distances = []
        for i in range(0,len(structures)):
            struct = structures[i]
            dist5p, dist3p = {}, {} 
            for cbgpos in struct:
                cbg  = self.codingblockgraphs[cbgpos]
                for node, omsr in cbg.overall_minimal_spanning_range().iteritems():
                    org = cbg.organism_by_node(node)
                    if dist5p.has_key(org):
                        if min(omsr)*3 < dist5p[org]:
                            dist5p[org] = min(omsr)*3
                    else:
                        dist5p[org] = min(omsr)*3
                    if dist3p.has_key(org):
                        if ( org2lengths[org]  - ( max(omsr)*3 ) ) < dist3p[org]:
                            dist3p[org] = org2lengths[org]  - ( max(omsr)*3 )
                    else:
                        dist3p[org] = org2lengths[org]  - ( max(omsr)*3 )

            # obtain average distances
            dist5pav = sum(dist5p.values()) / len(dist5p)
            dist3pav = sum(dist3p.values()) / len(dist3p)
            distances.append( ( abs(dist5pav-dist3pav), i, dist5pav, dist3pav, dist5p.values(), dist3p.values() ) )

            #omsrMin = self.codingblockgraphs[struct[0]].overall_minimal_spanning_range()
            #omsrMax = self.codingblockgraphs[struct[-1]].overall_minimal_spanning_range()
            #dist5p, dist3p = 0, 0
            #for (orgid,orfid) in omsrMin.keys():
            #    dist5p += ( min( omsrMin[(orgid,orfid)] ) * 3 )
            #for (orgid,orfid) in omsrMax.keys():
            #    dist3p += ( org2lengths[orgid] - ( max( omsrMax[(orgid,orfid)] ) * 3 ) )
            #dist5p = float(dist5p) / float(len(omsrMin))
            #dist3p = float(dist3p) / float(len(omsrMax))
            #distances.append( ( abs(dist5p-dist3p), i, dist5p, dist3p ) )
            
        # order the distances => best struct is now on top!
        distances.sort()
        #########################################################
        if verbose: print "ORDERING:", distances
        #########################################################
        # position 2 in first tuple of distances is id of best struct
        best_structure = structures[ distances[0][1] ]

        # 09/03/2010: ordering based on subGSG partial weight is
        # DANGEROUS!! In case of microsyntheny, this still has a high risk
        # of falure. So -> stick to the original distance selection
        # for the *best* GSG 
        # obtain overall total_weight() of CBGs in the structures to
        # define the BEST structure
        #weights = []
        #for i in range(0,len(structures)):
        #    struct = structures[i]
        #    totwt = 0
        #    for pos in struct:
        #        totwt = self.codingblockgraphs[pos].total_weight()
        #    weights.append( ( totwt,i ) )
        ## reverse order the weights => best struct is now on top!
        #weights.sort()
        #weights.reverse()
        #########################################################
        #if verbose: print "ORDERING totwt:", weights 
        #########################################################
        ## position 2 in first tuple of weights is id of best struct
        ####best_structure = structures[ weights[0][1] ]

        # return this best structure
        return best_structure

     # end of function _define_best_structure 

# end of class IntergenecityFunctions


def _cbg_intergenic_distance_analyses(thisCBG,nextCBG):
    """
    Get distance data between 2 CodingBlockGraphs

    @type  thisCBG: CodingBlockGraph
    @param thisCBG: CodingBlockGraph

    @type  nextCBG: CodingBlockGraph
    @param nextCBG: CodingBlockGraph

    @attention: thisCBG MUST BE 5'/left/upstream of nextCBG

    @rtype:  tuple ( int(), int(), float() )
    @return: ( mindist, maxdist, avdist )
    """
    distances = thisCBG.distance_between_codingblocks( nextCBG ).values()
    # convert AA distances to nt distances
    distances = [ d*3 for d in distances ]
    mindist   = min(distances)
    maxdist   = max(distances)
    avdist    = sum(distances)/len(distances)
    # return obtained distance values
    return (mindist, maxdist, avdist)

# end of function _cbg_intergenic_distance_analyses

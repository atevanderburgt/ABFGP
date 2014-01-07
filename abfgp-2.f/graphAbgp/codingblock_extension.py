"""
Functions for extension of CodingBlockGraphs (CBGs)

Some definitions:

+ spanningrange difference
    The difference between the MAXSR and the OMSR, defined on one side of the CBG (left/5p or rigth/3p)
    Each node/organism in a CBG has its own spanningrange difference

+ overall minimal spanning range (OMSR)
    The AA coordinate range of a CBG that is covered by all pairwise alignments objects (PacbPORFS)
    Each node/organism in a CBG has its own OMSR

+ maximal spanning range (MAXSR)
    The largest AA coordinate range of a CBG, that is covered by at least a single pairwise alignments objects (PacbPORFS)
    Each node/organism in a CBG has its own MAXSR

"""

# graphAbgp Imports
from exceptions import *
import codingblock_splitting

# Python Imports
from sets import Set

# Global variables
from settings.codingblockgraph import (
    CBG_MIN_AA_LENGTH,
    CBG_EXTEND_PACBPORF_MIN_LENGTH_RATIO,
    CBG_EXTEND_PACBPORF_MIN_IDENTITYSCORE_RATIO,
    )

def extend_on_spanningrange_difference(self,side=None,
    sprdif_min_aa_length=CBG_MIN_AA_LENGTH,
    extended_pacbporf_min_length_ratio=CBG_EXTEND_PACBPORF_MIN_LENGTH_RATIO,
    extended_pacbporf_min_identityscore_ratio=CBG_EXTEND_PACBPORF_MIN_IDENTITYSCORE_RATIO,
    verbose=False):
    """
    """
    # input integrity check
    if side not in ('left','rigth'):
        message = "`side` must be 'left' or 'rigth', not '%s'" % side
        raise InproperlyAppliedArgument, message

    # the nodes that are represented in the spanningrange difference
    sprdif = codingblock_splitting.spanningrange_difference(self,side,
            sprdif_min_aa_length=sprdif_min_aa_length,
            sprdif_min_node_count=self.node_count(),
            correct_sprdif_all_nodes=False)

    if not sprdif:
        # nothing to be extended -> done
        return False

    ########################################################################
    if verbose:
        print "SPRDIF extended", side, "INIT"
        print self, len(self.pacbps), self.edge_count()
        for k,v in self.overall_minimal_spanning_range().iteritems():
            print k,len(v),"--",
        print ""
    ########################################################################

    # if this point is reached, this CBG is extendable
    # make Sets of the lists in sprdif values
    for node,vlist in sprdif.iteritems(): sprdif[node] = Set(vlist)

    # if here, start extending the codingblockgraph
    # IMPORTANT!!! DO NOT!!! clear_cache() of codingblock!!
    toberemoved = {}

    # loop over the pacbporfs and gather those that have
    # no/little overlap with the observed sprdif on that side
    for (k,n1,n2),pacbporf in self.pacbps.iteritems():
        qRange = pacbporf.alignment_protein_range_query()
        sRange = pacbporf.alignment_protein_range_sbjct()
        overlap = len(sprdif[n1].intersection(qRange)) + len(sprdif[n2].intersection(sRange))
        if overlap < (sprdif_min_aa_length/5)*2:
            # add this key and pacbporf to dict toberemoved
            toberemoved[(k,n1,n2)] = pacbporf
            ############################################################################
            if verbose:
                print pacbporf, "%s-%s" % ( min(sprdif[n1]), max(sprdif[n1]) ),
                print "%s-%s" % ( min(sprdif[n2]), max(sprdif[n2]) ),
                print len(sprdif[n1].intersection(qRange)),
                print len(sprdif[n2].intersection(sRange))
            ############################################################################

    # if no pacbps are registered for removal -> no extention
    if not toberemoved: return False

    # remove the pacbporfs that lack overlap with the sprdif
    for (k,n1,n2) in toberemoved.keys():
        del( self.pacbps[(k,n1,n2)] )
        self.del_edge(n1,n2)

    # make new pacbps for the missing edges and update edge weights
    self.make_pacbps_for_missing_edges(use_maxsr=True)

    IS_EXTENDED = False 
    for (k,n1,n2) in toberemoved.keys():
        if not self.has_edge(n1,n2):
            # clustalw creation of pacbp failed aparently
            # set back this one and continue
            pacbporf = toberemoved[(k,n1,n2)]
            self.add_edge(n1,n2,wt=pacbporf.bitscore) 
            self.pacbps[(k,n1,n2)] = pacbporf 
            if verbose: print "CLUSTALW FAILED!?"
            # now to the next one
            continue
        else:
            ####################################################
            if verbose:
                print "orig & new PacbPORF:", n1, n2
                print toberemoved[(k,n1,n2)]
                #toberemoved[(k,n1,n2)].print_protein() 
            ####################################################
            # this PacbPORF edge is recreated by ClustalW
            pass

        # get pacbporf of this cbg 
        pacbporf = self.get_pacbps_by_nodes(node1=n1,node2=n2)[0]
        #####################################################
        #if verbose:
        #    print pacbporf
        #    pacbporf.print_protein()
        #####################################################

        try:
            length_ratio        = len(pacbporf) / len(toberemoved[(k,n1,n2)])
            identityscore_ratio = pacbporf.identityscore / toberemoved[(k,n1,n2)].identityscore
        except:
            # except in case of ZeroDivisionError
            length_ratio        = 0.0
            identityscore_ratio = 0.0

        if length_ratio >= extended_pacbporf_min_length_ratio and\
        identityscore_ratio >= extended_pacbporf_min_identityscore_ratio:
            # set pacbporf.source to clustalw-EXTENDED
            pacbporf.source = 'clustalw-EXTENDED' 
            IS_EXTENDED = True
        else:
            # hmmm... it is even truncated! reset to the original one
            self.remove_pacbp(pacbporf,n1,n2)
            self.pacbps[(k,n1,n2)] = toberemoved[(k,n1,n2)]
            # log print messages in verbose mode
            ########################################################################
            if verbose:
                print "FAILED:  ", pacbporf
                print "FAILED:  ", length_ratio,
                print length_ratio >= extended_pacbporf_min_length_ratio,
                print identityscore_ratio,
                print identityscore_ratio >= extended_pacbporf_min_identityscore_ratio 
                print "RESTORED:", toberemoved[(k,n1,n2)]
            ########################################################################

    # final check for succes rate: does CBG still has OMSR?
    # do check this, first clear its cached OMSR
    self.clear_cache()
    if not self.has_overall_minimal_spanning_range():
        # hmm... exceptional case in which OMSR got lost due to
        # extended pacbporfs (no idea why this can happen...).
        # Restore one by one ( in random order ;-) the original
        # pacbporfs untill OMSR is restored
        keys = toberemoved.keys()
        while keys and not self.has_overall_minimal_spanning_range():
            (k,n1,n2) = keys.pop()
            # get current (updated) pacbporf and replace with original one
            pacbporf = self.get_pacbps_by_nodes(node1=n1,node2=n2)[0]
            self.remove_pacbp(pacbporf,n1,n2)
            self.pacbps[(k,n1,n2)] = toberemoved[(k,n1,n2)]
    # recreate the cache of the CBG
    self.create_cache()
    self.update_edge_weights_by_minimal_spanning_range()
    ########################################################################
    if verbose:
        print "SPRDIF extended", side, "AFTER"
        print self, len(self.pacbps), self.edge_count()
        for k,v in self.overall_minimal_spanning_range().iteritems():
            print k,len(v),"--",
        print ""
    ########################################################################
    self.create_cache()
    self.update_edge_weights_by_minimal_spanning_range()

    # return status of extension (True)
    return IS_EXTENDED

# end of function extend_on_spanningrange_difference


def extend_on_left_spanningrange_difference(self,
    sprdif_min_aa_length=CBG_MIN_AA_LENGTH):
    """
    Extend CBG on left spanningrange difference

    @attention: see extend_on_spanningrange_difference for documentation
    """
    return extend_on_spanningrange_difference(self,side='left',
            sprdif_min_aa_length=sprdif_min_aa_length)

# end of function extend_on_left_spanningrange_difference


def extend_on_rigth_spanningrange_difference(self,
    sprdif_min_aa_length=CBG_MIN_AA_LENGTH):
    """
    Extend CBG on rigth spanningrange difference

    @attention: see extend_on_spanningrange_difference for documentation
    """
    return self.extend_on_spanningrange_difference(self,side='rigth',
            sprdif_min_aa_length=sprdif_min_aa_length)

# end of function extend_on_rigth_spanningrange_difference

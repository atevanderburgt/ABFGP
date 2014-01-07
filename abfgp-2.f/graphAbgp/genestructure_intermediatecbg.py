"""
Class with functions for finding/optimizing intermediate CBGs
in a GeneStructureOfCodingBlocks used in Alignment Based Gene Prediction
"""

# Module metadata
__authors__ = "Ate van der Burgt"
__license__ = "MIT"

# graphAbgp Imports
from exceptions import *
import ordering
import conversion
from graph_exoncollection import *
from graph_pacbpcollection import PacbpCollectionGraph
from genestructure_cbgaddition import findmostlikelyCBG2GSGinsert 

# Gene Imports
from gene.donor import SpliceDonor
from gene.acceptor import SpliceAcceptor

# Abgp Imports
from codingblockgraphinterface import CodingBlockGraphInterface
from lib_tinyexononorf import bridge_two_pacbporfs_by_tinyexon
from lib_stopwatch import StopWatch
from lib_hmm import (
    cbghmmsearch2pacbpcollection,
    hmmbuild_protein, hmmsearch_protein, hmmhit2pacbp,
    create_hmmdb_for_neighbouring_cbgs,
    tinyexonhmmsearch
    )

from lib_clustalw import clustalwinput2cbg

from lib_fasta import writeMultiFasta, parseFasta
from abgp_warnings import UnexpectedEventWarning

# Python imports
from os import remove as osremove
from copy import deepcopy

# Global variables
from settings.genestructure import (
        # tinyexon settings
        TINYEXON_MAX_NT_LENGTH,
        TINYEXON_MIN_NT_LENGTH,
        TINYEXON_MAX_INTRON_NT_LENGTH, 
        # scaffold gap omsr offset
        SCAFFOLD_GAP_OMSR_OFFSET,
        # intermediate CBG removal
        INTERMEDIATE_CBG_REMOVAL_MAX_AA_LENGTH,
        INTERMEDIATE_CBG_REMOVAL_MAX_GTG_ID_RATIO,
        )
from settings.splicesites import (
        IC_DONOR_PATTERN_OFFSET,
        IC_ACCEPTOR_PATTERN_OFFSET,
        )
from settings.codingblockgraph import (
        MAX_SMALL_CBG_HMM_COMPLETION_GTG_TOPO_DIF,
        MAX_SMALL_CBG_HMM_COMPLETION_GTG_ABS_DIF,
        MAX_SMALL_CBG_HMM_COMPLETION_GTG_ID_RATIO,
        )

class IntermediateCodingBlockGraphFunctions:

    def non_scaffold_analysis_known_sites_intermediate_frameshiftcbg(self,
        min_exon_nt_length=TINYEXON_MIN_NT_LENGTH,
        max_exon_nt_length=TINYEXON_MAX_NT_LENGTH,
        max_intron_nt_length=TINYEXON_MAX_INTRON_NT_LENGTH,
        verbose=False):
        """
        Do non-scaffold analysis on the GSG with a potential frameshift CBG

        @type  min_exon_nt_length: integer
        @param min_exon_nt_length: positive minimum exon nt length

        @type  max_exon_nt_length: integer
        @param max_exon_nt_length: positive maximum exon nt length

        @type  max_intron_nt_length: integer
        @param max_intron_nt_length: positive maximum intron length to take into
                                     acount when selecting suitable ORFs

        @type  verbose: Boolean
        @param verbose: print extended debugging on STDOUT (True) or be quiet (False)
        """
        replacement_done = False
        for pos in range(len(self)-3,-1,-1):
            prev,inter,next = self.codingblockgraphs[pos:pos+3]
            # error-check CBG types (no lsrCBG)
            if prev.__class__.__name__ == 'LowSimilarityRegionCodingBlockGraph':
                continue
            if inter.__class__.__name__ == 'LowSimilarityRegionCodingBlockGraph':
                continue
            if next.__class__.__name__ == 'LowSimilarityRegionCodingBlockGraph':
                continue

            # error-check cbgIF of inter CBG; both cbgIFs must be incompatible
            if inter.omsrlength() > 10:
                continue
            if not inter._CBGinterface5p:
                continue
            if not inter._CBGinterface3p:
                continue
            if inter._CBGinterface5p.is_compatible():
                continue
            if inter._CBGinterface3p.is_compatible():
                continue

            # error-check cbgIF of prev CBG; cbgIF must be optimal donor
            if not prev._CBGinterface3p:
                continue
            if not prev._CBGinterface3p._optimal_aligned_donor:
                continue
            if prev._CBGinterface3p._optimal_aligned_donor.__class__.__name__ ==\
            'DonorSiteCollectionGraph':
                continue

            # error-check cbgIF of next CBG; cbgIF must be optimal acceptor
            if not next._CBGinterface5p:
                continue
            if not next._CBGinterface5p._optimal_aligned_acceptor:
                continue
            if next._CBGinterface5p._optimal_aligned_acceptor.__class__.__name__ ==\
            'AcceptorSiteCollectionGraph':
                continue

            # only run on all distinct nodes
            # TODO change function construct_intermediate_tinyexon_with_known_splicesites()
            if prev.mutual_nodes(next):
                continue

            # (temporarily) remove this inter(mediate) CBG
            inter = self.codingblockgraphs.pop(pos+1)
            ########################################################
            if verbose:
                print "TRYING TO REPLACE intermediateCBG:"
                print inter
            ########################################################

            # and run the function!
            status = self.construct_intermediate_tinyexon_with_known_splicesites(
                    prev,next,
                    min_exon_nt_length=min_exon_nt_length,
                    max_exon_nt_length=max_exon_nt_length,
                    max_intron_nt_length=max_intron_nt_length,
                    verbose=verbose)

            if not status:
                # not succesful -> reset the just removed inter CBG
                self.codingblockgraphs.insert(pos+1,inter)
                self.codingblockgraphs[pos]._CBGinterface3p = None
                self.codingblockgraphs[pos+1]._CBGinterface5p = None
                self.codingblockgraphs[pos+1]._CBGinterface3p = None
                self.codingblockgraphs[pos+2]._CBGinterface5p = None
                # and recreate the interface as it was...
                self.create_cbginterfaces()
                ############################################################
                if verbose: print "TRYING TO REPLACE intermediateCBG FAILED"
                ############################################################
            else:
                ############################################################
                if verbose: print "SUCCECFULL intermediateCBG REPLACED"
                ############################################################
                replacement_done = True

        # correct possible erroneous setted FIRST/LAST CBGs
        self.finalize_genestructure()

        # return Boolean value weather or not anything is replaced
        return replacement_done

    # end of function non_scaffold_analysis_known_sites_intermediate_frameshiftcbg


    def non_scaffold_analysis_known_sites(self,
        min_exon_nt_length=TINYEXON_MIN_NT_LENGTH,
        max_exon_nt_length=TINYEXON_MAX_NT_LENGTH,
        max_intron_nt_length=TINYEXON_MAX_INTRON_NT_LENGTH,
        verbose=False):
        """
        Do non-scaffold analysis on the GSG

        @type  min_exon_nt_length: integer
        @param min_exon_nt_length: positive minimum exon nt length

        @type  max_exon_nt_length: integer
        @param max_exon_nt_length: positive maximum exon nt length

        @type  max_intron_nt_length: integer
        @param max_intron_nt_length: positive maximum intron length to take into
                                     acount when selecting suitable ORFs

        @type  verbose: Boolean
        @param verbose: print extended debugging on STDOUT (True) or be quiet (False)
        """
        for pos in range(len(self)-2,0,-1):
            prev,next = self.codingblockgraphs[pos-1:pos+1]
            if prev.__class__.__name__ == 'LowSimilarityRegionCodingBlockGraph':
                continue
            if next.__class__.__name__ == 'LowSimilarityRegionCodingBlockGraph':
                continue
            # error-check cbgIF of prev CBG
            if not prev._CBGinterface3p:
                continue
            if not prev._CBGinterface3p._optimal_aligned_donor:
                continue
            if prev._CBGinterface3p._optimal_aligned_donor.__class__.__name__ ==\
            'DonorSiteCollectionGraph':
                continue
            # error-check cbgIF of next CBG
            if not next._CBGinterface5p:
                continue
            if not next._CBGinterface5p._optimal_aligned_acceptor:
                continue
            if next._CBGinterface5p._optimal_aligned_acceptor.__class__.__name__ ==\
            'AcceptorSiteCollectionGraph':
                continue
        
            # only run on all distinct nodes
            # TODO change function construct_intermediate_tinyexon_with_known_splicesites()
            if prev.mutual_nodes(next):
                continue
        
            # and run the function!
            self.construct_intermediate_tinyexon_with_known_splicesites(
                    prev,next,
                    min_exon_nt_length=min_exon_nt_length,
                    max_exon_nt_length=max_exon_nt_length,
                    max_intron_nt_length=max_intron_nt_length,
                    verbose=verbose)

        # correct possible erroneous setted FIRST/LAST CBGs
        self.finalize_genestructure()

    # end of function non_scaffold_analysis_known_sites


    def non_scaffold_analysis_random_sites(self,verbose=False,**kwargs):
        """
        Do non-scaffold analysis on the GSG

        @type  verbose: Boolean
        @param verbose: extended debugging on STDOUT (True) or quiet (False)

        @attention: for other **kwargs, see _get_single_highest_scoring_
                    intermediate_tinyexon_with_random_splicesites()

        @rtype:  list
        @return: list with Booleans for if novel CBG(s) are created
        """
        statusses = []
        for pos in range(len(self)-2,0,-1):
            prev,next = self.codingblockgraphs[pos-1:pos+1]
            if prev.__class__.__name__ == 'LowSimilarityRegionCodingBlockGraph':
                continue
            if next.__class__.__name__ == 'LowSimilarityRegionCodingBlockGraph':
                continue
            # error-check cbgIFs of prev and next CBG
            if not prev._CBGinterface3p:
                continue
            if not next._CBGinterface5p: continue
            if prev._CBGinterface3p.is_optimal():
                continue 

            # and run the function!
            is_created = self.construct_intermediate_tinyexon_with_random_splicesites(
                    prev,next,verbose=verbose,**kwargs)

            # append is_created to status list
            statusses.append( is_created ) 

        # correct possible erroneous setted FIRST/LAST CBGs
        self.finalize_genestructure()

        # return status list
        return statusses

    # end of function non_scaffold_analysis_random_sites


    def construct_intermediate_tinyexon_with_known_splicesites(self,prev,next,
        min_exon_nt_length=TINYEXON_MIN_NT_LENGTH,
        max_exon_nt_length=TINYEXON_MAX_NT_LENGTH,
        max_intron_nt_length=TINYEXON_MAX_INTRON_NT_LENGTH,
        verbose=False):
        """
        Construct a novel tiny CBG in between two (ordered) CodingBlockGraphs

        @type  self: GenestructureOfCodingBlockGraphs
        @param self: GenestructureOfCodingBlockGraphs instance
        
        @type  prev: CodingBlockGraph (or None)
        @param prev: CodingBlockGraph upstream/5p of potential novel tiny CBG
        
        @type  next: CodingBlockGraph (or None)
        @param next: CodingBlockGraph downstream/3p of potential novel tiny CBG

        @type  min_exon_nt_length: integer
        @param min_exon_nt_length: positive minimum exon nt length

        @type  max_exon_nt_length: integer
        @param max_exon_nt_length: positive maximum exon nt length
        
        @type  max_intron_nt_length: integer
        @param max_intron_nt_length: positive maximum intron length to take into
                                     acount when selecting suitable ORFs

        @type  verbose: Boolean
        @param verbose: print extended debugging on STDOUT (True) or be quiet (False)
        """
        # check if CBGinterfaces are present in both CBGs
        # if not -> this function is not possible!
        if not prev._CBGinterface3p:
            return False
        if not prev._CBGinterface3p._splicedonorgraph:
            return False
        if not prev._CBGinterface3p._splicedonorgraph.alignedsites:
            return False
        if not next._CBGinterface5p:
            return False
        if not next._CBGinterface5p._spliceacceptorgraph:
            return False
        if not next._CBGinterface5p._spliceacceptorgraph.alignedsites:
            return False

        # start the StopWatch for performance timing
        stw = StopWatch()
        stw.start() 

        # intialize empty ExonCollectionGraph
        ECG = ExonCollectionGraph()

        # Loop over the organism identifiers in prev (or next) CBG
        # and try to make intermediate tinyexon(s) by using the already
        # assigned optimal donor & acceptor site of prev and next CBG
        for organism in prev.organism_set():
            # get node identifier of organism identifier
            node = prev.node_by_organism(organism)

            # organism identifier must be present in the best aligned
            # splice site graph of prev and next CBG
            # if organism not in prev._CBGinterface3p._splicedonorgraph.alignedsites[0].organism_set(): continue
            # if organism not in next._CBGinterface5p._spliceacceptorgraph.alignedsites[0].organism_set(): continue
            if organism not in prev._CBGinterface3p._optimal_aligned_donor.organism_set():
                continue
            if organism not in next._CBGinterface5p._optimal_aligned_acceptor.organism_set():
                continue

            # get current optimal donor & acceptor object
            #donor = prev._CBGinterface3p._splicedonorgraph.alignedsites[0].get_organism_objects(organism)[0]
            #acceptor = next._CBGinterface5p._spliceacceptorgraph.alignedsites[0].get_organism_objects(organism)[0]
            donor = prev._CBGinterface3p._optimal_aligned_donor.get_organism_objects(organism)[0]
            acceptor = next._CBGinterface5p._optimal_aligned_acceptor.get_organism_objects(organism)[0]

            # get list of possible intermediate ORFs
            orflist  = self.input[organism]['orfs'].get_elegiable_orfs(
                    max_orf_start=acceptor.pos,
                    min_orf_end=donor.pos, 
                    )

            # get the Orf objects of prev and next CBG and current organism
            donorOrf = prev.get_orfs_of_graph(organism=organism)[0]
            accepOrf = next.get_orfs_of_graph(organism=organism)[0]

            # obtain list of possible intermediate tinyexons
            tinyexons= bridge_two_pacbporfs_by_tinyexon(
                    donorOrf,accepOrf,
                    preceding_donor_sites=[donor],
                    subsequent_acceptor_sites=[acceptor],
                    orflist=orflist,
                    min_tinyexon_nt_length=min_exon_nt_length,
                    max_tinyexon_nt_length=max_exon_nt_length,
                    max_tinyexon_intron_nt_length=max_intron_nt_length,
                    )

            for tinyexon in tinyexons:
                ################################################
                if verbose:
                    print organism, tinyexon.orf.id, tinyexon
                ################################################
                # add tinyexon and its node to the ExonCollectionGraph
                node = (organism,tinyexon.orf.id,tinyexon.start,tinyexon.end)
                ECG.add_node_and_object(node,tinyexon)

        # create edges in the ECG and find complete subgraphs
        ECG.create_edges()

        if len(ECG.organism_edge_set()) != self.EXACT_SG_EDGE_COUNT:
            ####################################################################
            if verbose:
                print "len(ECG.organism_edge_set()) != ", 
                print "GSG.EXACT_SG_EDGE_COUNT", len(ECG.organism_edge_set())
            ####################################################################
            # not at least each org2org edge present -> no Ks graph creatable
            return False

        # find complete graphs
        exon_graphs = ECG.find_fully_connected_subgraphs()

        ########################################################################
        if verbose:
            print "ECG.node_count():", ECG.node_count(),
            print "ECG.edge_count():", ECG.edge_count(),
            print "subgraphs:", len(exon_graphs)
        ########################################################################

        # when no complete exon_graphs -> return False
        if not exon_graphs: return False

        ########################################################################
        if verbose: print stw.lap(), len(exon_graphs), exon_graphs[0]
        ########################################################################

        # take the first graph as most likely
        exongraph = exon_graphs[0]

        # this graph must have at least 2 nodes i.o. to complete with HMMsearches
        if exongraph.organism_set_size() < 2: return False

        # convert to a CodingBlockGraph; when this returns a None -> not convertable
        newcbg = ExonCollectionGraph2CodingBlockGraph(exongraph)
        if newcbg == None: return False

        if newcbg.node_count() == self.EXACT_SG_NODE_COUNT:
            pass
        else:
            # complete with cbghmmsearch2pacbpcollection
            pacbpCollection = cbghmmsearch2pacbpcollection(newcbg,self.input,
                    next=next,prev=prev,hmmsearch_num_hits=2
                    )
            # get list of accepted CBGs
            accepted = conversion.pacbpCollection2AcceptedCodingBlockGraphs(
                    pacbpCollection,prev=prev,next=next
                    )
            if accepted and accepted[0].node_count() == prev.node_count():
                newcbg = accepted[0]
            else:
                # hmmsearch failed -> no compatible CBG -> return False
                return False

        # if this point is reached, see if we can make compatible cbgIFs
        # between (A) prev -- newcbg and (B) newcbg -- prev
        cbgIFa = CodingBlockGraphInterface(prev,newcbg)
        cbgIFa.harvest_splice_sites()
        cbgIFa.find_conserved_splice_sites()
        cbgIFb = CodingBlockGraphInterface(newcbg,next)
        cbgIFb.harvest_splice_sites()
        cbgIFb.find_conserved_splice_sites()

        ################################################################
        if verbose:
            print "# tiny intermediate non-scaffold exon:"
            print prev
            print cbgIFa
            print newcbg
            print cbgIFb
            print next
            if not cbgIFa.is_optimal():
                print cbgIFa
                cbgIFa.interfaceproperties()
            if not cbgIFb.is_optimal():
                print cbgIFb
                cbgIFb.interfaceproperties()
        ################################################################

        if cbgIFa.is_optimal() and cbgIFb.is_optimal():
            pass
        elif ( cbgIFa.is_optimal() and cbgIFb.is_compatible() ) or\
        ( cbgIFb.is_optimal() and cbgIFa.is_compatible() ):
            # left interface is perfect, rigth one has (minor) flaws
            # or
            # rigth interface is perfect, left one has (minor) flaws
            # Check if the currently existing interface is really a big mess.
            # If so -> accept. If not -> do not accept.
            phaseD = prev._CBGinterface3p._optimal_aligned_donor.phase()
            phaseA = prev._CBGinterface3p._optimal_aligned_acceptor.phase()
            if type(phaseD) == type(int()) and type(phaseA) == type(int()):
                if phaseD == phaseA: return False
                else:                pass # allow insertion of intermediate CBG!
            else:
                pass # allow insertion of intermediate CBG!
        else:
            # No logical is_optimal() interfaces. To risky to continue
            # trying to improve it, at the risk of high FP exon rate.
            return False

        # Hurray! we have 99% sure a new tiny intermediate exon discovered
        # 1)  Check if all other criteria are met
        # 2)  Prepare it for insertion into the GSG
        # 3)  Do actual insert
        checks = []
        newcbg.create_cache()
        cbggtg = newcbg.genetree()
        gtg = self.genetree()
        checks.append( gtg.graphalignmentdifference( cbggtg ) <=\
                MAX_SMALL_CBG_HMM_COMPLETION_GTG_TOPO_DIF )
        checks.append( gtg.absolutegraphalignmentdifference( cbggtg ) <=\
                MAX_SMALL_CBG_HMM_COMPLETION_GTG_ABS_DIF )
        checks.append( ( cbggtg.identity() / gtg.identity() ) >=\
                MAX_SMALL_CBG_HMM_COMPLETION_GTG_ID_RATIO )

        ################################################################
        if verbose: print checks
        ################################################################

        if not False in checks:
            # hurray! Try insert it into the GSG
            if self.add_codingblock(newcbg,only_try_adding=True):
                # insert will be succesfull! Set all cbgIF objects in place
                prev._CBGinterface3p   = cbgIFa
                newcbg._CBGinterface5p = cbgIFa
                newcbg._CBGinterface3p = cbgIFb
                next._CBGinterface5p   = cbgIFb
                # do actual insert and retrun status True!
                self.add_codingblock(newcbg)
                return True
            else:
                # insert will not be succesfull -> return False
                return False
        else:
            # checks failed -> return False
            return False

    # end of function construct_intermediate_tinyexon_with_known_splicesites


    def construct_intermediate_tinyexon_with_random_splicesites(self,prev,next,
        verbose=False,**kwargs):
        """
        Try to create a tinyCBG in between two adjacent CBGs with a poor CBGinterface
        
        @type  self: GenestructureOfCodingBlockGraphs
        @param self: GenestructureOfCodingBlockGraphs instance
        
        @type  prev: CodingBlockGraph (or None)
        @param prev: CodingBlockGraph upstream/5p of potential novel tiny CBG
        
        @type  next: CodingBlockGraph (or None)
        @param next: CodingBlockGraph downstream/3p of potential novel tiny CBG

        @attention: for other **kwargs, see _get_single_highest_scoring_intermediate_tinyexon_with_random_splicesites()

        @type  verbose: Boolean
        @param verbose: print extended debugging on STDOUT (True) or be quiet (False)

        @rtype:  Boolean
        @return: is a novel tiny CBG created or not
        """
        # check (again) if CBGinterfaces are present
        # if not -> this function is not possible!
        if not prev._CBGinterface3p:
            return False
        if not next._CBGinterface5p:
            return False
        if not next.node_count() == prev.node_count():
            return False

        # find the single highest scoring tinyexonCBG with
        # at least one is_optimal() cbgIF
        newtinycbg = self._get_single_highest_scoring_intermediate_tinyexon_with_random_splicesites(
                prev,next,verbose=verbose,**kwargs)

        # if newtinycbg == False -> no novel tinyexon creatable!
        if not newtinycbg: return False

        ############################################################
        if verbose:
            print "NOVEL tinyexonCBG with optimal interfaces:"
            print newtinycbg
            print newtinycbg._CBGinterface5p
            print newtinycbg._CBGinterface3p
        ############################################################

        # now check the interfaces; if only one is_optimal, do another
        # check between this newcbg and prev or next CBG in search
        # for yet ANOTHER tinyexonCBG. This situation DOES occur!!
        # An example is 2 tiny exons in NCU02366 ( MGG_03521 / abfgp_mgg0081 )
        cbgIFa = newtinycbg._CBGinterface5p
        cbgIFb = newtinycbg._CBGinterface3p

        if cbgIFa.is_optimal() and cbgIFb.is_optimal():
            # huray! novel tiny exon cbg found
            # save the cbgIFs and store into the GSG (self)
            NOVEL_TINYCBG_ADDED = self.add_codingblock(newtinycbg,omit_conditional_addition=True)
            if NOVEL_TINYCBG_ADDED:
                prev._CBGinterface3p = cbgIFa
                next._CBGinterface5p = cbgIFb
                # return True because a tinyexon CBG is added
                return True
            else:
                # not placeable in the GSG. That is VERY weird and
                # shouldnot be possible because it is verified in
                # the partialGSG in the _get_single_highest_scoring_
                # _intermediate_tinyexon_with_random_splicesites function
                message = "newtinyCBG placeable partialGSS but not in GSG"
                print UnexpectedEventWarning(message)
                return False
        elif cbgIFa.is_optimal():
            # The left/5p cbgIF is_optimal, the right/3p not.
            # That means that there is maybe another tinyexonCBG
            # in between newcbg and the CBG `next`
            ############################################################################
            if verbose: print "\nENTERING 2th iteration!!, prev - newcbg - ?? - next\n"
            ############################################################################
            newcbg2 = self._get_single_highest_scoring_intermediate_tinyexon_with_random_splicesites(
                    newtinycbg,next,verbose=verbose,**kwargs)

            if not newcbg2:
                # creation of another one failed -> do not store newcbg either
                return False

            # if here, check the interfaces
            cbgIFc = newcbg2._CBGinterface5p
            cbgIFd = newcbg2._CBGinterface3p
            # if one of both new interfaces is is_optimal,
            # store both of them. Actually, both of them SHOULD be
            # is_optimal, but now we have 2/3 cggIFs of 2 new tinyCBGs valid
            if cbgIFc.is_optimal or cbgIFd.is_optimal():
                # do not forget to update the cbgIF in between the 2 tiny new cbgs
                newtinycbg._CBGinterface3p = cbgIFc
                NOVEL_TINYCBG_1_ADDED = self.add_codingblock(
                        newtinycbg,omit_conditional_addition=True
                        )
                NOVEL_TINYCBG_2_ADDED = self.add_codingblock(
                        newcbg2,omit_conditional_addition=True
                        )
                if NOVEL_TINYCBG_1_ADDED and NOVEL_TINYCBG_2_ADDED:
                    prev._CBGinterface3p = cbgIFa
                    next._CBGinterface5p = cbgIFd # cbgIFd , not b!
                    # return True because a tinyexon CBG is added
                    return True
                else:
                    # not placeable in the GSG. That is VERY weird and
                    # shouldnot be possible because it is verified in
                    # the partialGSG in the _get_single_highest_scoring_
                    # _intermediate_tinyexon_with_random_splicesites function
                    message = "2 newtinyCBGs placeable partialGSS but not in GSG"
                    print UnexpectedEventWarning(message)
                    return False 
            else:
                # not a single cbgIF of the 2th tinycbg is_optimal
                return False

        elif cbgIFb.is_optimal():
            # The right/3p cbgIF is_optimal, the left/5p not.
            # That means that there is maybe another tinyexonCBG
            # in between `prev` CBG and newcbg
            ############################################################################
            if verbose: print "\nENTERING 2th iteration!!, prev - ?? - newcbg - next\n"
            ############################################################################
            newcbg2 = self._get_single_highest_scoring_intermediate_tinyexon_with_random_splicesites(
                    prev,newtinycbg,verbose=verbose,**kwargs)

            if not newcbg2:
                # creation of another one failed -> do not store newcbg either
                return False

            # if here, check the interfaces
            cbgIFc = newcbg2._CBGinterface5p
            cbgIFd = newcbg2._CBGinterface3p

            # if one of both new interfaces is is_optimal,
            # store both of them. Actually, both of them SHOULD be
            # is_optimal, but now we have 2/3 cggIFs of 2 new tinyCBGs valid
            if cbgIFc.is_optimal or cbgIFd.is_optimal():
                # do not forget to update the cbgIF in between the 2 tiny new cbgs
                newtinycbg._CBGinterface5p = cbgIFd
                NOVEL_TINYCBG_1_ADDED = self.add_codingblock(
                        newtinycbg,omit_conditional_addition=True
                        )
                NOVEL_TINYCBG_2_ADDED = self.add_codingblock(
                        newcbg2,omit_conditional_addition=True
                        )
                if NOVEL_TINYCBG_1_ADDED and NOVEL_TINYCBG_2_ADDED:
                    prev._CBGinterface3p = cbgIFa
                    next._CBGinterface5p = cbgIFd # cbgIFd , not b!
                    # return True because a tinyexon CBG is added
                    return True
                else:
                    # not placeable in the GSG. That is VERY weird and
                    # shouldnot be possible because it is verified in
                    # the partialGSG in the _get_single_highest_scoring_
                    # _intermediate_tinyexon_with_random_splicesites function
                    message = "2 newtinyCBGs placeable partialGSS but not in GSG"
                    print UnexpectedEventWarning(message)
                    return False
            else:
                # not a single cbgIF of the 2th tinycbg is_optimal
                return False

        else:
            # unexpected event; get_single_highest_scoring_
            # _intermediate_tinyexon_with_random_splicesites function
            # may only return newtinyCBGs with at least one cbgIF that
            # is_optimal(). Print warning and return False
            message = "newtinyCBG has no is_optimal() interface at all!"
            print UnexpectedEventWarning(message)
            return False

    # end of function construct_intermediate_tinyexon_with_random_splicesites


    def _get_single_highest_scoring_intermediate_tinyexon_with_random_splicesites(
        self,prev,next,
        min_exon_nt_length=TINYEXON_MIN_NT_LENGTH,
        max_exon_nt_length=TINYEXON_MAX_NT_LENGTH,
        max_intron_nt_length=TINYEXON_MAX_INTRON_NT_LENGTH,
        verbose=False):
        """
        Helper function for construct_intermediate_tinyexon_with_random_splicesites()
        
        @type  self: GenestructureOfCodingBlockGraphs
        @param self: GenestructureOfCodingBlockGraphs instance
        
        @type  prev: CodingBlockGraph (or None)
        @param prev: CodingBlockGraph upstream/5p of potential novel tiny CBG
        
        @type  next: CodingBlockGraph (or None)
        @param next: CodingBlockGraph downstream/3p of potential novel tiny CBG

        @type  min_exon_nt_length: integer
        @param min_exon_nt_length: positive minimum exon nt length

        @type  max_exon_nt_length: integer
        @param max_exon_nt_length: positive maximum intron exon nt length
        
        @type  max_intron_nt_length: integer
        @param max_intron_nt_length: positive maximum intron length to take into
                                     acount when selecting suitable ORFs

        @type  verbose: Boolean
        @param verbose: print extended debugging on STDOUT (True) or be quiet (False)

        @attention: only use this function from within
                    construct_intermediate_tinyexon_with_random_splicesites
        """
        # some variables that are not necessairy as function arguments
        # -- hmmsearch_num_hits; just enough compared to number of species/genes
        # -- take_best_per_phase_combi; only take the best (or best 2....?)
        #    highest scoring tiny exons per phase combi
        # -- take_best_tinyexons_per_organism; just several of the highest
        #    scoring tinyexons of current organism
        # -- omsr_2_mask_aa_length_correction; correct the max(prev) and
        #    min(next) OMSR for this value to limit the masked area a little bit
        hmmsearch_num_hits               = self.EXACT_SG_NODE_COUNT * 2
        take_best_per_phase_combi        = 1
        take_best_tinyexons_per_organism = 5
        omsr_2_mask_aa_length_correction = 6

        # get lengths of splice sites from settings.splicesite variables
        donor_pattern_length = sum(IC_DONOR_PATTERN_OFFSET) + 2
        accep_pattern_length = sum(IC_ACCEPTOR_PATTERN_OFFSET) + 2

        ############################################################
        if verbose:
            print prev
            print prev._CBGinterface3p._splicedonorgraph
            print next
            print next._CBGinterface5p._spliceacceptorgraph
        ############################################################

        # start the StopWatch for verbose logging
        stw = StopWatch()
        stw.start()

        # create a fasta database of ORFs in prev and next for hmm search
        fname_fasta_hmmdb = "tmp.tmp.tmp.mfa"

        # create a empty PCG; store the pacbporfs of both prev and next CBG
        PCG = PacbpCollectionGraph()
        for cbg in [prev,next]:
            # store the pacbps of this cbg to the new empty PCG
            for (k,qnode,snode), pacbporf in cbg.pacbps.iteritems():
                if qnode not in PCG.get_nodes(): PCG.add_node(qnode)
                if snode not in PCG.get_nodes(): PCG.add_node(snode)
                if not PCG.has_edge(qnode,snode):
                    PCG.add_edge(qnode,snode,pacbporf.bitscore)
                    PCG.pacbps[(k,qnode,snode)] = pacbporf

        # create fastaHMM database in the interface between prev and next CBG
        fasta_hmm_db_content = create_hmmdb_for_neighbouring_cbgs(self.input,prev,next)
        fh = open(fname_fasta_hmmdb,'w')
        fh.write(fasta_hmm_db_content)
        fh.close()

        ################################################################
        if verbose:
            dbseqs = parseFasta(open(fname_fasta_hmmdb).readlines())
            fastastats = {}
            for h,s in dbseqs.iteritems():
                org, txt, orfid = h.split("_")
                if not fastastats.has_key(org): fastastats[org] = [[],0]
                fastastats[org][0].append(orfid)
                fastastats[org][1]+= s.count("X")
            # print oneliner with content summary of the hmmdb
            print [ (org,data[0],data[1]) for org,data in fastastats.iteritems() ]
        ################################################################

        # get OMSRs of prev and next CBG to delimit elegiable orf set
        prevOmsr = prev.overall_minimal_spanning_range()
        nextOmsr = next.overall_minimal_spanning_range()

        # loop over the organisms in one of the CBGs,
        # and search the interface of prev and next CBGs
        # for PSSM-obtained tiny exons
        tinyexons_found = 0
        for organism in prev.organism_set():
            # nodes and ORFs of prev and next CBG
            prevNode = prev.node_by_organism(organism)
            nextNode = next.node_by_organism(organism)
            donorOrf = prev.get_orfs_of_graph(organism=organism)[0]
            accepOrf = next.get_orfs_of_graph(organism=organism)[0]

            # get list of possible intermediate ORFs
            orflist  = self.input[organism]['orfs'].get_elegiable_orfs(
                    max_orf_start=min(nextOmsr[nextNode])*3,
                    min_orf_end=max(prevOmsr[prevNode])*3,
                    )

            # remove the current donor & acceptor orf itself too
            # this means, we will NEVER find an inframe intron here,
            # but that chance is very small....
            TODO = True

            ############################################################
            if verbose: print organism, [ orf.id for orf in orflist ]
            ############################################################

            # now create dummyDonor & dummyAcceptor objects for each phase
            # and mine the orflist for tiny exons
            tinyexonlist = []
            for donorPhase in [0,1,2]:
                positionD = ( max(prevOmsr[prevNode]) + 1 ) * 3 + donorPhase
                dummyDonor = SpliceDonor(positionD," "*donor_pattern_length,donor="GT",pssm_score=10.0)
                dummyDonor.phase = donorPhase
                for accepPhase in [0,1,2]:
                    positionA = ( min(nextOmsr[nextNode]) + 1 ) * 3 + donorPhase
                    dummyAccep = SpliceAcceptor(positionA," "*accep_pattern_length,acceptor="AG",pssm_score=10.0)
                    dummyAccep.phase = accepPhase
                    tinyexons= bridge_two_pacbporfs_by_tinyexon(
                        donorOrf,accepOrf,
                        preceding_donor_sites=[dummyDonor],
                        subsequent_acceptor_sites=[dummyAccep],
                        orflist=orflist,
                        min_tinyexon_nt_length=min_exon_nt_length,
                        max_tinyexon_nt_length=max_exon_nt_length,
                        max_tinyexon_intron_nt_length=max_intron_nt_length,
                        )
                    if tinyexons:
                        # append only a single, highest PSSM scoring tiny exon for this phase combi
                        tinyexons = ordering.order_list_by_attribute(tinyexons,"pssm_score",reversed=True)
                        tinyexonlist.extend( tinyexons[0:take_best_per_phase_combi] )
                        tinyexons_found += len( tinyexons[0:take_best_per_phase_combi] )

            # order the tinyexons
            tinyexonlist =  ordering.order_list_by_attribute(tinyexonlist,"pssm_score",reversed=True)
            # apply only the best X to the ECG
            for tinyexon in tinyexonlist[0:take_best_tinyexons_per_organism]:
                ################################################################################
                if verbose: print organism, tinyexon.orf.id, tinyexon, tinyexon.proteinsequence()
                ################################################################################
                node = (organism,tinyexon.orf.id,tinyexon.start,tinyexon.end)

                # create single fasta hmm query for this tiny exon
                fname_fasta_hmm_query = "tinyexon_hmmquery_%s_%s.fa" % ( organism, tinyexon.orf.id )
                writeMultiFasta({'query': tinyexon.proteinsequence()},fname_fasta_hmm_query )

                # make hmmbuild file of the multiplealignment
                fname_hmmbuild = hmmbuild_protein( fname_fasta_hmm_query )
                # run hmmsearch on this fasta database with hmmbuild file
                results = hmmsearch_protein( fname_hmmbuild, fname_fasta_hmmdb, params= {'E':150,'A': hmmsearch_num_hits} )
                querynode = ( organism, tinyexon.orf.id )
                for hmmhit in results:
                    sbjctorfid = int(hmmhit[0].split('_')[-1])
                    sbjctorg   = hmmhit[0].replace("_orf_%s" % sbjctorfid,"")
                    if sbjctorg == organism: continue
                    sbjctnode = ( sbjctorg, sbjctorfid )
                    sbjctorf  = self.input[sbjctorg]['orfs'].get_orf_by_id(sbjctorfid)
                    querycoords = (tinyexon.protein_start(), tinyexon.protein_end() )
                    key, pacbporf = hmmhit2pacbp(tinyexon.orf, organism, querycoords, sbjctorf, sbjctorg, hmmhit)
                    #print pacbporf
                    if pacbporf:
                        if sbjctnode not in PCG.get_nodes(): PCG.add_node(sbjctnode)
                        if querynode not in PCG.get_nodes(): PCG.add_node(querynode)
                        PCG.add_edge( querynode, sbjctnode, pacbporf.bitscore )
                        PCG.pacbps[ key ] = pacbporf

                # file cleanup of hmm query and profile
                osremove(fname_fasta_hmm_query)
                osremove(fname_hmmbuild)

        # file cleanup of hmm search database
        osremove(fname_fasta_hmmdb)

        # check if any tiny exon was added to the PCG
        # if not -> no tinyexonCBG findable -> return False
        if tinyexons_found == 0:
            return False

        ################################################################################
        if verbose:
            print stw.lap()
            print PCG
            print len(PCG.pacbps), PCG.edge_count(), PCG.node_count()
        ################################################################################

        # create a deepcopy of this PCG to keep its pacbps
        dpcpPCG = deepcopy(PCG)
        if verbose: print stw.lap(), "deepcopied PCG"

        # TODO: next, we search for Ks graphs, no missing edges.
        # as soon as there are tinyexon(s) in 2 species/organisms/genes,
        # this edge will not be reported! Such a situation is solved
        # in the ExonCollectionGraph by creating missing edges. Maybe this
        # functionality can be copied to the PCG graph !? However, that might be
        # VERY dangerous is all other PCG graphs except this one. Because there
        # usually are many nodes in a PCG, htis function can run forever, and is likely
        # to over-predict nonsense edges (edges are created with ClustalW alignments!)
        TODO = True

        # split in complete graphs and order by total weight
        splitted_subgraphs = PCG.find_fully_connected_subgraphs(
            edges=self.EXACT_SG_NODE_COUNT-1 , max_missing_edges=0 )
        splitted_subgraphs = ordering.order_graphlist_by_total_weight(splitted_subgraphs)


        # create a partial GSG; import GenestructureOfCodingBlockGraphs here
        # to prevent circular import in the head of this library
        from graph_genestructure import GenestructureOfCodingBlockGraphs
        partialGSG = GenestructureOfCodingBlockGraphs(self.input)
        partialGSG.codingblockgraphs = [ prev, next ]


        ################################################################################
        if verbose:
            print "Ks graphs in the PCG:"
            for spl in splitted_subgraphs:
                print spl.total_weight(), "\t", spl.get_ordered_nodes()
            print "PartialGSG to test placeability in GSG:"
            print prev
            print next
        ################################################################################


        for newcbg in splitted_subgraphs:
            # in the list of splitted_subgraphs the prev and next CBG
            # itself will be present too! -> ignore these
            if len(prev.node_set().symmetric_difference(newcbg.get_nodes())) == 0:
                continue
            if len(next.node_set().symmetric_difference(newcbg.get_nodes())) == 0:
                continue
            if newcbg.__class__.__name__ == 'PacbpCollectionGraph':
                 continue
            # set minimal OMSR size to 1 for this tiny new CBG
            newcbg.MINIMAL_OVERAL_SPANNING_RANGE_SIZE = 1

            # recruite the pacbporf objects from parental deepcopied PCG
            try:
                newcbg._recrute_pacbporfs_from_parental_cbg(dpcpPCG)
            except:
                # pacbporf recruitement failed -> no OMSR!!
                ################################################################
                if verbose:
                    print False, newcbg.total_weight(), "\t",
                    print newcbg.get_ordered_nodes(), "NO OMSR!!"
                ################################################################
                continue


            # check for OMSR of this CBG -> no OMSR, no CBG!!
            if newcbg.has_overall_minimal_spanning_range():
                # update edge weights -> it is a tiny CBG!
                newcbg.update_edge_weights_by_minimal_spanning_range()

                # try if this thingy fits in the partialGSG
                is_placeable = partialGSG.add_codingblock(
                        newcbg,only_try_adding=True,omit_conditional_addition=True
                        )
                ################################################################
                if verbose: print is_placeable, newcbg
                ################################################################
                if is_placeable:
                    newcbg.create_cache()
                    ############################################################
                    if verbose:
                        newcbg.printmultiplealignment()
                        print newcbg._cexpander.binarystring,
                        print newcbg._cexpander.projected_on
                    ############################################################
                    # verify the cbgInterfaces of this newcbg;
                    # if they are optimal, we have a novel tinycbg discovered!

                    # create cbgIF between prev and newcbg CBG
                    cbgIFa = CodingBlockGraphInterface(prev,newcbg)
                    # overrule non-introns for mutual nodes
                    for node in prev.mutual_nodes(newcbg):
                        org = prev.organism_by_node(node)
                        cbgIFa._interface_is_intron[org] = False
                    cbgIFa.harvest_splice_sites()
                    cbgIFa.find_conserved_splice_sites()

                    # create cbgIF between newcbg and next CBG
                    cbgIFb = CodingBlockGraphInterface(newcbg,next)
                    # overrule non-introns for mutual nodes
                    for node in next.mutual_nodes(newcbg):
                        org = next.organism_by_node(node)
                        cbgIFb._interface_is_intron[org] = False
                    cbgIFb.harvest_splice_sites()
                    cbgIFb.find_conserved_splice_sites()

                    # try to (further) optimize the cbgIFs with the special
                    # optimizetinyexoninterface function if not optimal yet
                    if not cbgIFa.is_optimal():
                        cbgIFa.optimizetinyexoninterface(verbose=True)
                    if not cbgIFb.is_optimal():
                        cbgIFb.optimizetinyexoninterface(verbose=True)

                    ############################################################
                    if verbose: print cbgIFa, "\n", cbgIFb
                    ############################################################

                    if cbgIFa.is_optimal() or cbgIFb.is_optimal():
                        # this novel tinyexonCBG has an is_optimal()
                        # interface on at least one side.
                        newcbg._CBGinterface5p = cbgIFa
                        newcbg._CBGinterface3p = cbgIFb
                        # Return it after saving the cbgIFs into it
                        return newcbg
                    else:
                        # no is_optimal interface on either side
                        # this is very unlikely to be a bonafide CBG
                        # so, try the next new tinyexon CBG in the list
                        continue
                else:
                    # not placeable in the partialGSG -> continue
                    continue
            else:
                continue

        # if here, no new tinyCBG found
        return False

    # end of function _get_single_highest_scoring_intermediate_tinyexon_with_random_splicesites


    def cbgpair_blastp2pacbpcol_analyses(self,
        omit_non_cbg_orfs=False,
        omit_cbg_orfs=False,
        max_cbg_aa_distance=None,
        verbose=False):
        """
        """
        from cbgjunction2blastp import blastanalysescbgjunction
        novel_cbgs_created = 0
        for pos in self.cbgpositerator(reversed=True)[1:]:
            prevCBG,nextCBG = self.codingblockgraphs[pos:pos+2]
            # skip this interface in all these cases
            if prevCBG.IS_IGNORED:                              continue
            if prevCBG._short_name == "lsrCBG":                 continue
            if prevCBG.node_count() < self.EXACT_SG_NODE_COUNT: continue
            if nextCBG.IS_IGNORED:                              continue
            if nextCBG._short_name == "lsrCBG":                 continue
            if nextCBG.node_count() < self.EXACT_SG_NODE_COUNT: continue
            # skip identical node Sets; this function might be
            # called when lsrCBGs are not installed yet in the GSG!
            if not nextCBG.node_set().symmetric_difference(prevCBG.node_set()):
                continue
            if max_cbg_aa_distance and max_cbg_aa_distance <\
            min(prevCBG.distance_between_codingblocks(nextCBG).values()):
                continue 

            if prevCBG._CBGinterface3p and nextCBG._CBGinterface5p:
                # skip when interface is perfect or near-perfect;
                # prevCBG._CBGinterface3p == nextCBG._CBGinterface5p
                if prevCBG._CBGinterface3p.is_optimal():        continue
            else:
                # create the interface here on the fly!
                cbgIF = CodingBlockGraphInterface(prevCBG,nextCBG)
                cbgIF.harvest_splice_sites()
                cbgIF.find_conserved_splice_sites()
                # place the interface objects into the CBGs
                prevCBG._CBGinterface3p = cbgIF
                nextCBG._CBGinterface5p = cbgIF
                # skip when interface is perfect or near-perfect;
                # prevCBG._CBGinterface3p == nextCBG._CBGinterface5p
                if prevCBG._CBGinterface3p.is_optimal():        continue

            #if prevCBG._CBGinterface3p.is_optimal_donor():      continue
            #if prevCBG._CBGinterface3p.is_optimal_acceptor():   continue

            # create a partial GSG; import GenestructureOfCodingBlockGraphs here
            # to prevent circular import in the head of this library
            from graph_genestructure import GenestructureOfCodingBlockGraphs
            partGSG = GenestructureOfCodingBlockGraphs(self.input)
            partGSG.codingblockgraphs = [ prevCBG,nextCBG ]
            partGSG._GENETREE = self._GENETREE

            ####################################################################
            if verbose:
                print prevCBG
                print prevCBG._CBGinterface3p,
                print prevCBG._CBGinterface3p.is_optimal(),
                print prevCBG._CBGinterface3p.is_optimal_donor(),
                print prevCBG._CBGinterface3p.is_optimal_acceptor(),
                print prevCBG.distance_between_codingblocks(nextCBG)
                print nextCBG
            ####################################################################

            # get list of elegiable (tiny) CBGs in between both CBGs
            cbglist = blastanalysescbgjunction(self,prevCBG,nextCBG,
                            omit_non_cbg_orfs=omit_non_cbg_orfs,
                            omit_cbg_orfs=omit_cbg_orfs,
                            verbose=verbose)

            ####################################################################
            if verbose:
                if cbglist:
                    print "cbgpair_blastp2pacbpcol_analyses, cbglist:",
                    print len(cbglist)
                for cbg in cbglist:
                    print cbg
                    cbg.printmultiplealignment()
            ####################################################################

            partGSG = findmostlikelyCBG2GSGinsert(partGSG,cbglist,verbose=verbose)
            if len(partGSG) > 2:
                # update the slice into the main GSG
                self.codingblockgraphs.__setslice__(
                        pos, pos+2,
                        partGSG.codingblockgraphs )
                # increase counter
                novel_cbgs_created+=len(partGSG)-2
                # do NOT try other detection function; done here!
                continue

        # correct possible erroneous setted FIRST/LAST CBGs
        self.finalize_genestructure()

        # return counter of how much new CBGs inserted into the GSG
        return novel_cbgs_created 

    # end of function cbgpair_blastp2pacbpcol_analyses

    
    def discovered_independant_closeby_intron_gains(self,verbose=False):
        """
        """
        novel_cbgs_created = 0
        for pos in self.cbgpositerator(reversed=True)[1:]:
            prevCBG,nextCBG = self.codingblockgraphs[pos:pos+2]
            # skip this interface in all these cases
            if prevCBG.IS_IGNORED:                              continue
            if prevCBG._short_name == "lsrCBG":                 continue
            if prevCBG.node_count() < self.EXACT_SG_NODE_COUNT: continue
            if nextCBG.IS_IGNORED:                              continue
            if nextCBG._short_name == "lsrCBG":                 continue
            if nextCBG.node_count() < self.EXACT_SG_NODE_COUNT: continue
            # skip this interface when it is not set yet
            if not prevCBG._CBGinterface3p:                     continue
            if not nextCBG._CBGinterface5p:                     continue
            # skip when interface is perfect or near-perfect;
            # prevCBG._CBGinterface3p == nextCBG._CBGinterface5p
            if prevCBG._CBGinterface3p.is_optimal():            continue
            if prevCBG._CBGinterface3p.is_optimal_donor():      continue
            if prevCBG._CBGinterface3p.is_optimal_acceptor():   continue
    
    
            # if here, then start performing this function!
            ############################################################
            if verbose:
                print "\n\n>>", pos, "discovered_independant_closeby_intron_gains"
                print prevCBG
                print prevCBG._CBGinterface3p
                print nextCBG
            ############################################################


            # create a partial GSG; import GenestructureOfCodingBlockGraphs here
            # to prevent circular import in the head of this library
            from graph_genestructure import GenestructureOfCodingBlockGraphs
            partGSG = GenestructureOfCodingBlockGraphs(self.input)
            partGSG.codingblockgraphs = [ prevCBG,nextCBG ]
            partGSG._GENETREE = self._GENETREE
            bckp_cbgIF_prev = prevCBG._CBGinterface3p   
            bckp_cbgIF_next = nextCBG._CBGinterface5p

            # (try to) split on left tiny sprdif
            partGSG.split_cbgs_on_left_spanningrange_difference(
                    sprdif_min_aa_length=3,
                    sprdif_min_node_count=2,
                    ignore_first_cbg=True,
                    perform_cbgif_optimization=True,
                    verbose=verbose)

            ############################################################
            if verbose:
                print "discovered_independant_closeby_intron_gains:"
                for _cbg in partGSG.codingblockgraphs[1:-1]:
                    print "in partGSG:",_cbg
            ############################################################

            intermediate_cbg_is_added = False
            if len(partGSG) >= 3:
                # len() can be >3 in case of lsrCBG creation. Chance is very small...
                cbgIFaCheck = partGSG.codingblockgraphs[0]._CBGinterface3p.optimalitycheck()
                cbgIFbCheck = partGSG.codingblockgraphs[-1]._CBGinterface5p.optimalitycheck()
                if cbgIFaCheck.count(True) >= 2 and cbgIFbCheck.count(True) >= 2:
                    # yep, this perfectly fits; update the partGSG
                    # into the main GSG (self) with a slice update
                    self.codingblockgraphs.__setslice__(
                            pos, pos+2, partGSG.codingblockgraphs )
                    novel_cbgs_created += len(partGSG)-2
                    if verbose:
                        print ">> DISCOVERED independant_closeby_intron_gains"
                        for cbg in partGSG: print cbg
                    # do NOT try to do the split on the rigth side anymore
                    continue

            if not intermediate_cbg_is_added:
                # reset original cbgIF and remove potential IS_SPLITTED status
                prevCBG._CBGinterface3p = bckp_cbgIF_prev 
                nextCBG._CBGinterface5p = bckp_cbgIF_next 
                prevCBG.IS_3P_SPLITTED = False
                nextCBG.IS_5P_SPLITTED = False
                if not prevCBG.IS_5P_SPLITTED: prevCBG.IS_SPLITTED = False
                if not nextCBG.IS_3P_SPLITTED: nextCBG.IS_SPLITTED = False


            ####################################################################
            if verbose:
                print "\nno LEFT independant_closeby_intron_gains"
                print prevCBG
                print prevCBG._CBGinterface3p
                print nextCBG
                print ""
            ####################################################################
    
            # (re)create partial GSG
            partGSG = GenestructureOfCodingBlockGraphs(self.input)
            partGSG.codingblockgraphs = [ prevCBG,nextCBG ]
            partGSG._GENETREE = self._GENETREE
            bckp_cbgIF = prevCBG._CBGinterface3p


            # (try to) split on rigth tiny sprdif
            partGSG.split_cbgs_on_rigth_spanningrange_difference(
                    sprdif_min_aa_length=3,
                    sprdif_min_node_count=2,
                    ignore_final_cbg=True,
                    perform_cbgif_optimization=True,
                    verbose=verbose)
    
            intermediate_cbg_is_added = False
            if len(partGSG) >= 3:
                # len() can be >3 in case of lsrCBG creation. Chance is very small...
                cbgIFaCheck = partGSG.codingblockgraphs[0]._CBGinterface3p.optimalitycheck()
                cbgIFbCheck = partGSG.codingblockgraphs[-1]._CBGinterface5p.optimalitycheck()
                if cbgIFaCheck.count(True) >= 2 and cbgIFbCheck.count(True) >= 2:
                    # yep, this perfectly fits; update the partGSG
                    # into the main GSG (self) with a slice update
                    self.codingblockgraphs.__setslice__(
                            pos, pos+2, partGSG.codingblockgraphs )
                    novel_cbgs_created += len(partGSG)-2
                    intermediate_cbg_is_added = True
                    if verbose:
                        print ">> DISCOVERED independant_closeby_intron_gains"
                        for cbg in partGSG: print cbg
                    # continue to the next pair of CBGs
                    continue

            if not intermediate_cbg_is_added:
                # reset original cbgIF and remove potential IS_SPLITTED status
                # reset original cbgIF and remove potential IS_SPLITTED status
                prevCBG._CBGinterface3p = bckp_cbgIF_prev
                nextCBG._CBGinterface5p = bckp_cbgIF_next
                prevCBG.IS_3P_SPLITTED = False
                nextCBG.IS_5P_SPLITTED = False
                if not prevCBG.IS_5P_SPLITTED: prevCBG.IS_SPLITTED = False
                if not nextCBG.IS_3P_SPLITTED: nextCBG.IS_SPLITTED = False

            ####################################################################
            if verbose:
                print "\nno RIGTH independant_closeby_intron_gains"
                print prevCBG
                print prevCBG._CBGinterface3p
                print nextCBG
                print ""
            ####################################################################

        # correct possible erroneous setted FIRST/LAST CBGs
        self.finalize_genestructure()
    
        # done. Return the counter
        return novel_cbgs_created
    
    # end of function discovered_independant_closeby_intron_gains


    def scaffold_analysis(self,scaffold_gap_omsr_offset=SCAFFOLD_GAP_OMSR_OFFSET,
        verbose=False):
        """
        Create tinyexons in between CBGs where the scaffold assumes continuity
        but still contains a gap

        @type  scaffold_gap_omsr_offset: integer
        @param scaffold_gap_omsr_offset: additional AA positions around the
                                         scaffold gap to add to the sequence
                                         in search for potential tiny exons
                                         (default 1 AA)

        @attention: choose scaffold_gap_omsr_offset very small (default 1 AA)

        @type  verbose: Boolean
        @param verbose: print extended debugging on STDOUT (True) or be quiet (False)

        @explanation:   This requires some explanation. The scaffold of the GSG
                        is build up from the ordered Orf objects (labeled with
                        unique integer identifiers) of distinct organisms/genes
                        in the CBGs. A (single) intron gain or loss event will
                        result in 2 adjacent CBGs that have identical Orf
                        identifiers for some of the organism/gene indentifiers
                        but distinct Orf identifiers for some others. The first
                        case most likely resembles no intron but a continious
                        exon, whereas the second case definately resembles an
                        intron that bridges the CBGs for these organism
                        identifier(s). For the continious exon case(s), the
                        ordered CBGs provide a measurement of how many putative
                        Amino Acid positions should be available in/around the
                        interface between the CBGs. When the intron-case (distinct
                        Orf identifiers) seems to be (somewhat) shorter, this is
                        a sign that an additional tiny exon is present in this
                        interface. These tiny exons are picked up by this function
                        when the (short!!) Amino Acid sequence is modest to
                        strongly conserved.

        @rtype:  integer
        @return: integer count of how many novel tiny CBGs are created
        """
        # start the performance StopWatch timer
        stw = StopWatch(name='scaffoldAnalyses')
        stw.start()
        tinyexons_created = 0

        # take elegiable range (all interfaces, so omit first and last)
        cbgevalrange = range(len(self)-1,0,-1)
        for cbgpos in cbgevalrange:
            # get combinations of 2 neighbouring CBGs
            (first,second) = self.codingblockgraphs[cbgpos-1:cbgpos+1]
            poscombi = (cbgpos-1,cbgpos)

            # ignore if one of them IS_IGNORED
            if True in [ cbg.IS_IGNORED for cbg in (first,second) ]:
                continue
            # ignore if combination is created by a split
            if [ first.IS_3P_SPLITTED, second.IS_5P_SPLITTED ] == [True,True]:
                continue
            # ignore if all mutual nodes
            if not Set(first.get_nodes()).symmetric_difference(second.get_nodes()):
                continue
            # ignore if all distinct nodes 
            if not Set(first.get_nodes()).intersection(second.get_nodes()):
                continue
            # skip combinations that have an optimal CBGinterface
            if self.cbginterface_is_optimal_donor(first) and\
            self.cbginterface_is_optimal_acceptor(second):
                continue

            # gather distance per organism between codingblocks
            aadistomsr  = first.distance_between_codingblocks(second)
            aadistmaxsr = first.maxsr_distance_between_codingblocks(second)

            # gather mutual nodes & organisms
            mutual_nodes = Set(first.get_nodes()).intersection(second.get_nodes())
            mutual_orgs  = [ first.organism_by_node(node) for node in mutual_nodes ]    

            ############################################################
            if verbose:
                print first
                print mutual_nodes, "\t", mutual_orgs
                print second
                print stw.lap()
            ############################################################

            # continue if only a single mutual_node -> nothing to clustalw!
            # TODO: a single mutual_node is still maybe usefull. Make a workaround? 
            # TODO: as soon as such a case is identified, solve it. But, no such
            # TODO: case found yet for sure (02/11/09)
            if len(mutual_nodes) == 1: continue


            # check for bogus applied argument scaffold_gap_omsr_offset
            # it must be a value gte 0
            scaffold_gap_omsr_offset = max([0,scaffold_gap_omsr_offset])

            # gather sequence concerning the scaffold gap of the mutual nodes
            seqs, orfs, coords = {}, {}, {}
            for node in mutual_nodes:
                org = first.organism_by_node(node)
                sta = max( first.overall_minimal_spanning_range(node=node) ) -\
                        scaffold_gap_omsr_offset 
                end = min( second.overall_minimal_spanning_range(node=node) ) +\
                        scaffold_gap_omsr_offset 
                orf = first.get_orfs_of_graph(organism=org)[0]
                seq = orf.getaas(abs_pos_start=sta,abs_pos_end=end)
                seqs[org]   = seq
                orfs[org]   = orf
                coords[org] = [sta,end]

            # Check if sequences are succesfully obtained;
            # empty strings are possible for very small scaffold gaps
            if '' in seqs.values(): continue

            # get (a) protein similarity matrix of the CBG
            protsimmtrx = first.pacbps.values()[0].MATRIX

            # create novel CBG with ClustalW
            newcbg = clustalwinput2cbg(seqs,orfs,coords,mutual_nodes,
                    matrix = protsimmtrx, verbose=verbose)

            # if newcbg does not exist -> no proper K(s-x) CBG created
            if not newcbg:
                ################################################################
                if verbose: "clustalwinput2cbg did not produce CBG"
                ################################################################
                continue
            else:
                ################################################################
                if verbose:
                    print stw.lap(), "tinyHMM input:"
                    print newcbg
                    newcbg.printmultiplealignment()
                ################################################################

            # search the missing organisms for tiny exons 
            newcbglist = tinyexonhmmsearch(newcbg,self.input,
                    prev=first,next=second,verbose=False)

            # remove CBGs with a cexpander onlyzeros string
            acceptednewcbglist = []
            for newcbg in newcbglist:
                # remove non-CBGs (weird cases of 2-node-CBGs occur!)
                if newcbg.node_count() <= 2: continue
                # remove CBGs with unlikely nodes compared to its neighbours 
                if False == intermediateCBG_node_comparison(first,newcbg,second):
                    continue
                # update edge weight & create_cache() for cexpander data
                newcbg.update_edge_weights_by_minimal_spanning_range()
                newcbg.create_cache()
                if newcbg._cexpander.binarystring.find("1") >= 0:
                    acceptednewcbglist.append( newcbg )

            # order the CBG list on omsridentity()*total_weight()
            acceptednewcbglist = ordering.order_cbgraphlist(acceptednewcbglist)

            ####################################################################
            if verbose:
                print stw.lap(), "tinyHMM newcbglist:", len(acceptednewcbglist)
            ####################################################################

            # create a partialGSG of prev and next CBG ([ first, second ])
            # and try to add the novel tinyHMM CBGs in the interface
            from graph_genestructure import mutualorgsofcbgs2partialGSG
            partGSG = mutualorgsofcbgs2partialGSG(
                input=self.input,cbgs=[first,second]
                )
            # backup current cbgIF
            bckp_cbgIF_first  = first._CBGinterface3p
            bckp_cbgIF_second = second._CBGinterface5p 
            # copy genetree object from main GSG
            partGSG._GENETREE = self._GENETREE

            # store current partGSG length (==2) into variable
            curpartgsglen = len(partGSG)

            # (try to) insert the most likely tinyHMM CBGs in the partGSG
            partGSG = findmostlikelyCBG2GSGinsert(partGSG,acceptednewcbglist,
                    verbose=verbose)

            if len(partGSG) == curpartgsglen+1 and curpartgsglen == 2:
                nodesetA  = partGSG.codingblockgraphs[0].node_set()
                nodesetB  = partGSG.codingblockgraphs[1].node_set()
                nodesetC  = partGSG.codingblockgraphs[2].node_set()
                unionsetAC = nodesetA.union(nodesetC)
                if not nodesetB.difference(unionsetAC):
                    # All nodes in novel tiny exon CBG are covered
                    # by the neighboring CBGs. This represensts a case
                    # which was hopefully detected by the function
                    # discovered_independant_closeby_intron_gain().
                    # But, it didn't. Likely because it is not the case ;-)
                    # Here, detect and ignore these cases
                    print "DETECTED independant_closeby_intron_gain"
                    print "DETECTED in scaffold_analyses()"
                    print partGSG.codingblockgraphs[1]._CBGinterface5p
                    print partGSG.codingblockgraphs[1]
                    print partGSG.codingblockgraphs[1]._CBGinterface5p
                    print "DETECTED & omitted here..."

                    deleted = partGSG.codingblockgraphs.pop(1)
                    partGSG.codingblockgraphs[0]._CBGinterface3p = None
                    partGSG.codingblockgraphs[1]._CBGinterface5p = None
                    partGSG.create_cbginterfaces()

            if len(partGSG) > curpartgsglen:
                # yep, novel CBG added to partGSG -> update slice in main GSG
                stapos = poscombi[0]
                endpos = poscombi[1]+1
                self.codingblockgraphs.__setslice__(stapos,endpos,
                        partGSG.codingblockgraphs)
                # increate the tinyexons_created counter
                tinyexons_created += ( len(partGSG) - curpartgsglen )
                ################################################################
                if verbose:
                    print stw.lap(), "novel tinyHMM ADDED!!"
                    for tinyCBG in partGSG.codingblockgraphs[1:-1]:
                        print tinyCBG
                        tinyCBG.printmultiplealignment()
                        for key,pacbporf in tinyCBG.pacbps.iteritems():
                            if 'hca' in [key[0][0],key[1][0],key[2][0]]:
                                print pacbporf
                                pacbporf.print_protein()
                ################################################################
                continue
            else:
                # nope; make shure the CBGs in partGSG are intactly
                # left in place in the GSG (IS_SPLITTED attributes)
                # might have changed...
                first._CBGinterface3p  = bckp_cbgIF_first
                second._CBGinterface5p = bckp_cbgIF_second
                first.IS_3P_SPLITTED  = False
                second.IS_5P_SPLITTED = False
                if not first.IS_5P_SPLITTED:  first.IS_SPLITTED = False
                if not second.IS_3P_SPLITTED: second.IS_SPLITTED = False

        # correct possible erroneous setted FIRST/LAST CBGs
        self.finalize_genestructure()

        # return the number of tinyexons_created
        return tinyexons_created

    # end of function scaffold_analysis


    def skip_nonsense_intermediate_cbgs(self,verbose=False,**kwargs):
        """
        Remove short, low id% CBGs with crippled cbgIFs on both sides from the GSG
    
        @attention: see is_intermediate_cbg_elegiable_for_deletion for **kwargs
    
        @rtype:  Boolean
        @return: True or False, weather or not CBG(s) are removed
        """
    
        # loop in reversed order (in case of a deletion!)
        # over all intermediate CBGs (skip the first and last)
        IS_DELETED = False

        delete_positions = []

        for pos in range(len(self)-2,0,-1):
            # get current series of 3 CBGs: prev, this and next
            ( prevcbg, thiscbg, nextcbg ) = self.codingblockgraphs[pos-1:pos+2]
        
            # check if thisCBG is elegiable for deletion
            if not self.is_intermediate_cbg_elegiable_for_deletion(
                thiscbg,prevCBG=prevcbg,nextCBG=nextcbg,
                verbose=verbose,**kwargs):
                #################################################
                if verbose: print "NO DELETION!", pos 
                #################################################
                continue

            # append this position to putative delete positions
            delete_positions.append( pos )
    
            ############################################################
            if verbose:
                print "ELEGIABLE FOR DELETION!!!", pos
                print thiscbg._CBGinterface5p
                print thiscbg
                thiscbg.printmultiplealignment()
                print thiscbg._CBGinterface3p
            ############################################################

        # translate putative delete (integer) positions by merging
        # positions directly adjacent to each other
        grouped_delete_positions = []
        for pos in delete_positions:
            if not grouped_delete_positions:
                grouped_delete_positions = [ [ pos ] ]
            elif pos ==  grouped_delete_positions[-1][-1] - 1:
                grouped_delete_positions[-1].append( pos )
            else:
                grouped_delete_positions.append( [ pos ] )

        ############################################################
        if verbose:
            print "putative delete positions:", delete_positions
            print "grouped:", grouped_delete_positions
        ############################################################

        for positiongroup in grouped_delete_positions:
            # try to remove this (group of) CBGs that is selected
            # for putative deletion. The _try_remove_intermediate_cbg
            # function tries to create a nice cbgIF for the CBGs directly
            # neighboring the deletions
            is_deleted = _try_remove_intermediate_cbg(self,
                    positiongroup,verbose=verbose)

            if is_deleted:
                IS_DELETED = True
            elif not is_deleted and len(positiongroup) > 1:
                # okay, no succesfull deletion. If the effort contained
                # multiple positions, a second round of deletion is tried
                # but now on individual CBG positions
                for position in positiongroup:
                    # try to delete an individual position
                    is_deleted = _try_remove_intermediate_cbg(self,
                            [ position ], verbose=verbose ) 
                    if is_deleted:
                        # one of the >1 deletions was succesfull
                        # Do not try to delete the rest as well; this is
                        # a near-garanty for trouble ;-)
                        IS_DELETED = True
                        break
                    else:
                        pass    

            else:
                # no deletion ...
                pass                

        # correct possible erroneous setted FIRST/LAST CBGs
        self.finalize_genestructure()
   
        # return the IS_DELETED status variable
        return IS_DELETED
    
    # end of function skip_nonsense_intermediate_cbgs


    def is_intermediate_cbg_elegiable_for_deletion(self,cbg,
        prevCBG=None,nextCBG=None,
        cbg_max_aa_length=INTERMEDIATE_CBG_REMOVAL_MAX_AA_LENGTH,
        cbg_max_gtg_id_ratio=INTERMEDIATE_CBG_REMOVAL_MAX_GTG_ID_RATIO,
        verbose=False):
        """
        Is this intermediate CBG elegiable for deletion?
    
        @type  cbg: CodingBlockGraph
        @param cbg: CodingBlockGraph instance (with cbgIF objects!)
    
        @type  cbg_max_aa_length: positive integer
        @param cbg_max_aa_length: maximal AA OMSR length of the CBG to allow deletion
    
        @type  cbg_max_gtg_id_ratio: positive float
        @param cbg_max_gtg_id_ratio: maximal cbgGTG/theGTG ratio (0.0-1.0) to allow deletion
    
        @type  verbose: Boolean
        @param verbose: print extended debugging on STDOUT (True) or be quiet (False)
    
        @rtype:  Boolean
        @return: True or False
        """
        if cbg.__class__.__name__ == 'LowSimilarityCodingBlockGraph':
            if verbose: print "iicefd: lsrCBG"
            return False
        elif cbg.IS_SPLITTED:
            if verbose: print "iicefd: IS_SPLITTED"
            return False
        elif not cbg._CBGinterface5p:
            if verbose: print "iicefd: no cbgIF5p"
            return False
        elif not cbg._CBGinterface3p:
            if verbose: print "iicefd: no cbgIF3p"
            return False
        elif cbg.omsrlength() > cbg_max_aa_length:
            if verbose: print "iicefd: OMSR", cbg.omsrlength(), ">", cbg_max_aa_length
            return False
        elif cbg._CBGinterface5p.is_optimal():
            if verbose: print "iicefd: cbgIF5p optimal"
            return False
        elif cbg._CBGinterface3p.is_optimal():
            if verbose: print "iicefd: cbgIF3p optimal"
            return False
        elif cbg._CBGinterface5p.is_compatible() and\
        cbg._CBGinterface3p.is_compatible():
            if verbose: print "iicefd: cbgIF5p & cbgIF3p compatible"
            return False
        elif cbg.omsrlength() < cbg_max_aa_length/2 and prevCBG and nextCBG and\
        prevCBG._CBGinterface3p and nextCBG._CBGinterface5p and\
        prevCBG._CBGinterface3p.donor_phase() == nextCBG._CBGinterface5p.acceptor_phase():
            if verbose: print "iicefd: tinyCBG good cbgIF when skipped"
            return True
        elif ( cbg.genetree().identity() / self.genetree().identity() )\
        > cbg_max_gtg_id_ratio:
            if verbose: print "iicefd: gtgID% >", cbg_max_gtg_id_ratio 
            return False
        else:
            # all criteria are invalid: this CBG is elegiable
            # for deletion because it has all the characteristics
            # of a bogus multiple alignment
            if verbose: print "iicefd: DELETABLE!!"
            return True
    
    # end of function is_intermediate_cbg_elegiable_for_deletion


# end of class IntermediateCodingBlockGraphFunctions


def intermediateCBG_node_comparison(prevCBG,intermCBG,nextCBG,max_omsr_aa_gap=3):
    """
    Is this intermediate (tiny) CBG elegiable for the GSG given its neighbours? 

    @type  prevCBG: CodingBlockGraph
    @param prevCBG: CodingBlockGraph (with OMSR)
            
    @type  intermCBG: CodingBlockGraph
    @param intermCBG: CodingBlockGraph (with OMSR)

    @type  nextCBG: CodingBlockGraph
    @param nextCBG: CodingBlockGraph (with OMSR)

    @type  max_omsr_aa_gap: integer
    @param max_omsr_aa_gap: maximal AA distance between intermCBG and its
                            neighbour CBGs when all nodes are shared (default 3)

    @attention: prevCBG, itermCBG and nextCBG must be CBGs (with OMSR)
    @attention: prevCBG, itermCBG and nextCBG must be positionally ordered
    @attention: recommended not to adapt the max_omsr_aa_gap variable

    @rtype:  NoneBoolean
    @return: True (okay), False (not okay), None (not applicable)
    """
    if not prevCBG or not nextCBG:
        # not applicable -> return None
        return None
    if not intermCBG:
        raise NoCodingBlockGraphApplied

    # get the node Sets of these CBGs
    nodesP = prevCBG.node_set()
    nodesI = intermCBG.node_set()
    nodesN = nextCBG.node_set()

    if nodesI.difference( nodesP.union(nodesN) ):
        # Nodes in intermediateCBG that are absent in prevCBG and nextCBG
        # No grounds at all to discard this intermediateCBG (not applicable)
        return None
    elif not nodesP.difference(nodesN):
        # prevCBG and nextCBG have identical nodes AND intermCBG as well !?
        # Biological explanation(s):
        #   - A conserved inframe intron in all species
        #   - A (not dealth with) lsrCBG that is filled with a CBG
        #   - !?!?
        # In all cases: highly unlikely. Print some warning here
        return False
    elif not nodesP.difference(nodesI):
        # prevCBG and intermediateCBG have identical nodes.
        # This could be biologically True: a LowSimilarityRegion followed
        # by a stretch of higher similarity. Unlikely, but not impossible.
        # TODO: must we perform further checks here or not?
        # TODO: In the findmostlikelyCBG2GSGinsert() this case is dealth with
        # TODO: to some extend, but not perfectly yet.
        # TODO: True is returned, with it is a estimated guess, not a garantee
        return True
    elif not nodesN.difference(nodesI):
        # See the above explanation, but now for intermediateCBG and nextCBG
        return True
    else:
        # The case where this function was written for in the first place.
        # A (tiny)CBG that shares all nodes with the union of both parent nodes
        # Biological explanation(s):
        #   - independant cases of (closeby) intron gains
        #   - significantly lower sequence similarity in the precise intron
        #     region compared to the surrounding CBG regions
        #   - !?!?
        # The second case (low sequence similarity in the precise intron region)
        # will be a hard case to track and correctly recognize anyway.
        # The first case (independant closeby intron gain) is hard to detect,
        # and a special function is devoted to it. The intermediateCBG must
        # comply to some positioning criterion towards its neighbours to be
        # a proper case. Here, False cases are recognized.
        # The second case (low sequence similarity in the precise intron region)
        # will be a hard case to track anyway.
        PERFECTLY_FITTING_CLOSEBY_INTRON_GAIN = True
        intermOMSR = intermCBG.overall_minimal_spanning_range()
        for node in nodesI:
            if node in nodesN:
                if min(nextCBG.overall_minimal_spanning_range(node=node)) -\
                max(intermOMSR[node]) > max_omsr_aa_gap:
                    PERFECTLY_FITTING_CLOSEBY_INTRON_GAIN = False
                    break
            if node in nodesP:
                if max(prevCBG.overall_minimal_spanning_range(node=node)) -\
                min(intermOMSR[node]) < -max_omsr_aa_gap:
                    PERFECTLY_FITTING_CLOSEBY_INTRON_GAIN = False
                    break
        if not PERFECTLY_FITTING_CLOSEBY_INTRON_GAIN:
            return False
        else:
            # Likely a bona fide case of independant closeby intron gains
            return True

# end of function intermediateCBG_node_comparison


def _try_remove_intermediate_cbg(gsg,groupedpositions,verbose=False):
    """ """
    prevcbg = gsg.codingblockgraphs[ min(groupedpositions) - 1 ]
    nextcbg = gsg.codingblockgraphs[ max(groupedpositions) + 1 ]

    # check if the order nextcbg == prevcbg results in a perfect cbgIF    
    cbgIF = CodingBlockGraphInterface( prevcbg, nextcbg )
    cbgIF.harvest_splice_sites()
    cbgIF.find_conserved_splice_sites()
  
    ############################################################ 
    if verbose:
        print "NEW cbgIF IF DELETION WILL BE DONE:", groupedpositions
        print cbgIF
    ############################################################ 

    if cbgIF.is_optimal() or\
    ( cbgIF.is_optimal_donor() and cbgIF.is_compatible() ) or\
    ( cbgIF.is_optimal_acceptor() and cbgIF.is_compatible() ):
        # yes! The intermediate CBG is probably
        # a bogus multiple alignment. Remove it and
        # replace this cbgIF for the currently existing ones

        ########################################################
        if verbose:
            print "# deleting intermediate CBG", groupedpositions
            for position in groupedpositions:
                print gsg.codingblockgraphs[position]
        ########################################################

        # delete the CBGs from the GSG
        for position in groupedpositions:
            deletedcbg = gsg.codingblockgraphs.pop(position)

        # set the cbgIF in place
        prevcbg._CBGinterface3p = cbgIF
        nextcbg._CBGinterface5p = cbgIF

        # return status True for succesfull deletion
        return True

    else:
        # leave elegiable CBG in place; no strong interface
        # signal that pinpoints its deletion 
        return False

# end of function _try_remove_intermediate_cbg

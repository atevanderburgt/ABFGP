"""
Class with functions for finding/optimizing the final CBG
in a GeneStructureOfCodingBlocks used in Alignment Based Gene Prediction
"""

# Module metadata
__authors__ = "Ate van der Burgt"
__license__ = "MIT"

# grapgAbgp Imports
import graphAbgp
from exceptions import *
from graph_exoncollection import *
import ordering

# Gene Imports
from gene.codingblock import CodingBlockStart, CodingBlockEnd
from gene.exon import FinalExonOnOrf

# ABGP Imports
from lib_shorttailingexononorf import find_tailing_exon_on_orf
from lib_cexpander import cexpanderanalyses_omsr2orfend
from lib_hmm import cbghmmsearch2pacbpcollection
from lib_stopwatch import StopWatch
from codingblockgraphinterface import CodingBlockGraphInterface

# Python Imports
from copy import deepcopy
from sets import Set

# Global variables
from settings.genestructure import (
    SHORT_TAILINGEXON_MAX_NT_LENGTH, SHORT_TAILINGEXON_MAX_INTRON_NT_LENGTH,
    SHORT_TAILINGEXON_TAKE_MAX_BEST_ACCEPTORS, SHORT_TAILINGEXON_TAKE_MAX_BEST_ECGS,
    SHORT_TAILINGEXON_TAKE_MAX_BEST_CBGS
    )
from settings.codingblockgraph import (
    CBG_FINAL_SPRDIF_MIN_AA_LENGTH, CBG_FINAL_SPRDIF_MIN_NODE_COUNT,
    CBG_FINAL_SPRDIF_ONLY_IF_STOP_TW_RATIO_LTE, CBG_FINAL_SPRDIF_ONLY_IF_CBG_ID_GTE
    )
from settings.codingblockgraph import (
    MAX_CBG_REMOVE_NONSENSE_FINAL_AA_LENGTH,
    MAX_CBG_REMOVE_NONSENSE_FINAL_GTG_TOPO_DIF,
    MAX_CBG_REMOVE_NONSENSE_FINAL_GTG_ABS_DIF,
    MIN_CBG_REMOVE_NONSENSE_FINAL_GTG_ID_RATIO,
    MIN_CBG_REMOVE_NONSENSE_FINAL_TCODE_OMSR,
) 


class FinalCodingBlockGraphFunctions:
    """ """
    def get_final_cbg(self):
        """
        Get the CBG that is labeled as IS_LAST

        @rtype:  CodingBlockGraph
        @return: CodingBlockGraph instance
        """
        # get current final cbg
        for sg in self.codingblockgraphs:
            if sg.IS_LAST:
                return sg
        return None

    # end of function get_final_cbg


    def set_final_cbg(self,cbg):
        """
        Set this CBG in stead of the one that is labeled as IS_LAST

        @type  cbg:  CodingBlockGraph
        @param cbg: CodingBlockGraph instance
        """
        # get current final cbg
        for pos in range(0,len(self)):
            if self.codingblockgraphs[pos].IS_LAST: 
                self.codingblockgraphs[pos] = cbg
                break

    # end of function set_final_cbg


    def is_final_cbg_annotated_final_cbg(self):
        """
        Does the (current) last CBG correspond to all (annotated) final exons?

        @rtype:  Boolean
        @return: True only if all (annotated) final exons are in the final CBG
        """
        # no CBGs in the GSG yet!
        if len(self) == 0: return False

        # get set of nodes of last codingblock
        lastCBGnodeset = self.codingblockgraphs[-1].node_set()
        # generate a Set of last orf ids for the given gene structures
        annotatedlastnodeset = Set()
        for org in self.input.keys():
            if len(self.input[org]['orfid-genestructure']) >= 1:
                annotatedlastnodeset.add(
                    (org, self.input[org]['orfid-genestructure'][-1] )
                    )
        # check if all as final Exons annotated nodes are covered
        if len(annotatedlastnodeset) == len(lastCBGnodeset) and\
        len(lastCBGnodeset.difference(annotatedlastnodeset)) == 0:
            return True
        else:
            return False
        
    # end of function is_final_cbg_annotated_final_cbg


    def construct_final_tiny_cbg(self,
        max_exon_nt_length=SHORT_TAILINGEXON_MAX_NT_LENGTH,
        max_intron_nt_length=SHORT_TAILINGEXON_MAX_INTRON_NT_LENGTH,
        take_max_best_acceptors=SHORT_TAILINGEXON_TAKE_MAX_BEST_ACCEPTORS,
        take_max_best_ecgs=SHORT_TAILINGEXON_TAKE_MAX_BEST_ECGS,
        take_max_best_cbgs=SHORT_TAILINGEXON_TAKE_MAX_BEST_CBGS,
        maximal_current_stopcodongraph_average_weight=0.90,
        minimal_last_vs_new_identity_ratio=0.80,
        maximal_cexpander_cbg_tail_uniformity_aa_length=3,
        elegiable_donor_omsr_nt_offset=21,
        verbose=False):
        """
        Make a tiny final CBG by ``shooting tiny exons into the deep``
        """
        # get current last CBG
        last = self.get_final_cbg()

        # check if final tail of this CBG is uniformaly alignable
        cxpdrOutput = cexpanderanalyses_omsr2orfend(last)
        IS_UNIFORMLY_ALIGNED = True
        for trf in cxpdrOutput._transferblocks:
            if trf.binarystring[-maximal_cexpander_cbg_tail_uniformity_aa_length:].count("0"):
                IS_UNIFORMLY_ALIGNED = False
                break

        ############################################################
        if verbose:
            print "Cexpander uniformaly aligned:",
            print maximal_cexpander_cbg_tail_uniformity_aa_length,
            print "->", IS_UNIFORMLY_ALIGNED
            print "omsr:       ", last._cexpander.projected_on,
            print last._cexpander.binarystring
            trf = cxpdrOutput.get_transfer_of_projected_on(
                    last._cexpander.projected_on)
            if trf and trf != True:
                print "omsr2orfend:", last._cexpander.projected_on,
                print trf.binarystring
        ############################################################

        if IS_UNIFORMLY_ALIGNED:
            # break out of this function. Chance of overpredicting
            # a final tiny exon is bigger then finding a True one!
            return False

        # check if the stopcodongraph is not (very) good already
        if last._stopcodongraph.average_weight() >=\
        maximal_current_stopcodongraph_average_weight:
            # break out of this function. Chance of overpredicting
            # a final tiny exon is bigger then finding a True existing one
            return False

        # start the timer (performance benchmark in verbose mode)
        stw = StopWatch(name='stwFinalECG')
        stw.start()

        # get FinalExons on elegiable Orfs based on distance towards OMSR of
        # current last CBG and minimal acceptor site score
        omsr  = last.overall_minimal_spanning_range()
        maxsr = last.maximal_spanning_range()
        ECG = ExonCollectionGraph()

        ################################################################
        if verbose:
            print "currentLAST", last
            print last._stopcodongraph
            print last._stopcodongraph.is_optimal()
            for org in last.organism_set():
                print org, last._stopcodongraph.is_optimal(organism=org)
            for organism in last.organism_set():
                node = last.node_by_organism(organism)
                theorf = last.get_orfs_of_graph(organism=organism)[0]
                print organism, "\t", node, "\t", max(omsr[node]), "\t",
                print max(maxsr[node]), theorf.endPY/3
        ################################################################

        for organism in last.organism_set():
            node = last.node_by_organism(organism)
            # calculate an offset for the acceptor position
            # variable elegiable_acceptor_omsr_nt_offset is needed to
            # enlarge the OMSR definded offset. When the OMSR is by chance
            # a few nt or aa larger than the actual exon length, the true
            # acceptor position can be erroneously abandoned.
            offset = max(omsr[node]) * 3 - elegiable_donor_omsr_nt_offset 
            theorf = last.get_orfs_of_graph(organism=organism)[0]

            # check if this final orf is self can serve as a final extension
            remaining_orf_nt_length          = (theorf.protein_endPY - max(omsr[node])) * 3
            remaining_maxsr_nt_length        = (max(maxsr[node]) - max(omsr[node])) * 3
            remaining_maxsr_tostop_nt_length = (theorf.protein_endPY - max(maxsr[node])) * 3 


            FIND_NEW_FINAL_ORFS       = True
            STORE_CURRENT_ORF_AS_FIOO = False 
            if remaining_maxsr_nt_length >= max_exon_nt_length:
                # exceptionally large maxsr on rigth side of omsr
                # store as FIOO but to NOT search for an orf extension!
                ### FIND_NEW_FINAL_ORFS       = False # discarded 17/09/2009; when poos maxsr present, overruled!
                STORE_CURRENT_ORF_AS_FIOO = True
            elif remaining_maxsr_tostop_nt_length <= 18:
                # maxsr is less then 6 AA apart from stop on current orf
                #FIND_NEW_FINAL_ORFS       = False
                STORE_CURRENT_ORF_AS_FIOO = True
            elif remaining_orf_nt_length < max_exon_nt_length:
                # final piece of unaligned sequence is a perfect HMM seed
                STORE_CURRENT_ORF_AS_FIOO = True
            else:
                pass

            if STORE_CURRENT_ORF_AS_FIOO:
                cbs = CodingBlockStart( theorf.aapos2dnapos( max(omsr[node]) ) )
                # set pssm_score to (very) high; this rewards
                # using the current Orf as the last Orf
                cbs.pssm_score = 20.0
                fioo = FinalExonOnOrf(cbs,theorf.endPY,theorf)
                node = (organism,theorf.id,fioo.start,fioo.end)
                ECG.add_node_and_object(node,fioo)
                ################################################################
                if verbose:
                    print organism,theorf.id,"self==potential last exon", remaining_orf_nt_length
                    print organism, theorf.id, fioo, fioo.start,fioo.end, theorf.endPY
                ################################################################

            if not FIND_NEW_FINAL_ORFS:
                # quit here -> no orf extension of this CBG
                continue

            # get elegiable (new) final orfs
            orflist = self.input[organism]['orfs'].get_elegiable_orfs(
                    max_orf_start=offset+max_intron_nt_length,
                    min_orf_end=offset )
            ################################################################
            if verbose:
                print organism, [ orf.id for orf in orflist ], "offset:", offset, offset/3
            ################################################################
            for orf in orflist:
                results = find_tailing_exon_on_orf(
                        theorf,orf,
                        current_donor_pos=offset,
                        max_tailingexon_nt_length=max_exon_nt_length,
                        max_tailingexon_intron_nt_length=max_intron_nt_length,
                        )
                for exon,intron in results:
                    node = (organism,orf.id,exon.start,exon.end)
                    if node not in ECG.get_nodes():
                        ECG.add_node_and_object(node,exon)
                        if verbose: print organism, node, exon

        if verbose: print stw.lap(), "Exon objects gathered", ECG.node_count()

        # now take only the best `take_max_best_acceptors`
        # because there can be quite some of them!
        for organism in ECG.organism_set():
            objects = ordering.order_list_by_attribute( ECG.get_organism_objects(organism), order_by='pssm_score', reversed=True )
            for obj in objects[take_max_best_acceptors:]:
                node = (organism,obj.orf.id,obj.start,obj.end)
                ECG.del_node(node)
                if verbose: print "deleted:", node, obj.orf.id, obj.pssm_score

        ########################################################################
        if verbose:
            print stw.lap(), ">take_max_best_acceptors DELETED"
            for organism in ECG.organism_set():
                for obj in ordering.order_list_by_attribute(
                    ECG.get_organism_objects(organism),
                    order_by='pssm_score', reversed=True
                    ):
                    print "remaining", organism, obj.orf.id, obj.length, obj
        ######################################################################## 

        # only continue if all organisms are represented in the ECG
        if last.organism_set_size() > ECG.organism_set_size():
            if verbose: print "To few organisms/genes present -> return False"
            return False

        # create edges in the ECG between compatible phases and 
        # exon length, then make pacbps for these edges
        ECG.create_edges()
        ECG.make_pacbps_for_edges()
        if verbose:
            print stw.lap(), "edges + PACBPS created:", ECG.edge_count(), ECG.node_count(), len(ECG.pacbps)

        # search for complete graphs in this
        last_exon_graphs = ECG.find_fully_connected_subgraphs()

        ########################################################################
        if verbose: 
            print stw.lap(), "duration of ECG.find_fully_connected_subgraphs()",
            print len(last_exon_graphs)
        ########################################################################

        # only continue if there is an perfectly aligned last exon graph
        if not (last_exon_graphs and last_exon_graphs[0].connectivitysaturation() == 1.0):
            ####################################################################
            if verbose: print "no perfect aligned last exon graph -> return False"
            ####################################################################
            return False

        # convert to CodingBlockGraphs
        new_last_cbgs = []
        for leg in last_exon_graphs[0:take_max_best_ecgs]:
            cbg = ExonCollectionGraph2CodingBlockGraph(leg,is_last=True,lastCBG=last)
            if cbg != False and cbg != None and cbg.organism_set_size() == last.organism_set_size():
                # create cache of CBG and do final check on quality
                cbg.create_cache()
                if (cbg.total_weight() < 0 or cbg.omsrlength() <= 10) and\
                cbg._cexpander.binarystring.find("1") == -1:
                    # discard hardly alignable CBGs
                    continue
                # if here, then append this cbg as a possible novel final CBG
                new_last_cbgs.append( cbg )
                ################################################################
                if verbose: print "LEGcbg", cbg
                ################################################################

        ########################################################################
        if verbose: print stw.lap(), "ECGs converted to CBGs", len(new_last_cbgs)
        ########################################################################

        if not new_last_cbgs:
            ####################################################################
            if verbose: print "no ecgs convertable to CBGs -> return False"
            ####################################################################
            return False

        # order by total weight, get the optimal CBG and its corresponding ECG
        new_last_cbgs = ordering.order_graphlist_by_total_weight(new_last_cbgs)
        theNewLastCbg = None
        cbgIF = None


        # check all interfaces between the novel final CBGs and the previous
        # CBG. The best interface is added to the GSG!
        cbgif_accepted_new_last_cbgs = []
        already_checked_node_sets = []

        for newcbg in new_last_cbgs[0:take_max_best_cbgs]:
            lastExonGraph = newcbg._ExonCollectionGraph
            del( newcbg._ExonCollectionGraph )

            # check if it is not the extention of the current
            # last CBG (identical nodes)
            if len(last.node_set().symmetric_difference(newcbg.node_set())) == 0:
                if verbose: print "newCBG is the extention of current last CBG!!"
                continue

            # check if this combination of nodes (orfs) has not been tried already
            if newcbg.get_ordered_nodes() in already_checked_node_sets:
                ###############################################################
                if verbose: 
                    print "newCBG node set done earlier:", 
                    print newcbg.get_ordered_nodes()
                ###############################################################
                continue
            else:
                # append this set of nodes (as a list) to checklist
                already_checked_node_sets.append( newcbg.get_ordered_nodes() )

            # check if this new final tinyexon graph has a compatible interface
            # with the current last one
            cbgIF = CodingBlockGraphInterface(last,newcbg)
            cbgIF.harvest_splice_sites()
            distinct_orgs = []
            for node in lastExonGraph.get_nodes():
                exon = lastExonGraph.get_node_object(node)
                if exon.acceptor.__class__.__name__ == 'SpliceAcceptor':
                    distinct_orgs.append( lastExonGraph.organism_by_node(node) )
            cbgIF.allow_intron_in_organisms(distinct_orgs)
            cbgIF.find_conserved_splice_sites()
            # do NOT optimize -> consumes a lot of time and is helpfull
            # only in extreme cases...
            #cbgIF.optimize()

            if not cbgIF.is_compatible():
                ################################################################
                if verbose:
                    print "newCBG not a is_compatible() cbgIF"
                    print newcbg
                ################################################################
                continue

            # append to cbgif_accepted_new_last_cbgs
            newcbg._CBGinterface5p = cbgIF
            cbgif_accepted_new_last_cbgs.append(
                    (
                        cbgIF.optimalitycheck().count(True),
                        newcbg.total_weight(),
                        newcbg
                    )
                )

        ########################################################################
        if verbose:
            print stw.lap(), "cbgIFs checked %s/%s" % (
                len(cbgif_accepted_new_last_cbgs),
                len(new_last_cbgs[0:take_max_best_cbgs])
                )
        ########################################################################
        # now start by adding the highest scoring newcbg first
        cbgif_accepted_new_last_cbgs.sort()
        cbgif_accepted_new_last_cbgs.reverse()

        ########################################################################
        if verbose:
            print "candidate novel final CBGs:", len(cbgif_accepted_new_last_cbgs)
            for (true_cnt,totalwt,newcbg) in cbgif_accepted_new_last_cbgs:
                print true_cnt,totalwt,newcbg._CBGinterface5p
                print newcbg
        ########################################################################

        for (true_cnt,totalwt,newcbg) in cbgif_accepted_new_last_cbgs:
            # get the already created cbgIF from the newcbg graph
            cbgIF = newcbg._CBGinterface5p
    
            # now check 4 criteria:
            # (1) cbgIF.is_optimal() (2) >GTG.identity
            # (3) >STG.totalweight   (4) <STG.distance
            criteria = []
            criteria.append( cbgIF.is_optimal() )
            criteria.append( newcbg._stopcodongraph.total_weight() > last._stopcodongraph.total_weight() )
            criteria.append( newcbg.genetree().identity() > last.genetree().identity() )
            criteria.append( newcbg._stopcodongraph.stopcodon2omsrdistance() <= last._stopcodongraph.stopcodon2omsrdistance() )

            ####################################################################
            if verbose:
                print "TRYING ADDITION of final newcbg", criteria
                print true_cnt,totalwt,newcbg._CBGinterface5p
                print newcbg
            ####################################################################

            # check if there is only a single different node/orf changed in the newcbg
            # this is recognized by a symmetric_difference of size 2 
            # in this case, be very strict! This easily causes overprediction (FP) tiny exons 
            if len(last.node_set().symmetric_difference(newcbg.node_set())) == 2:
                # check if 4 criteria are valid;
                # a single False results in not accepting this new last tiny cbg
                if False in criteria:
                    if verbose: print "# NOVEL lastTinyExon discarded; single orf extension, criteria", criteria
                    # continue -> no new tiny CBG
                    continue

            # now start check the criteria.
            # if criteria[0] == True, means a fully is_optimal interface!
            # do not perform any additional check, just add!
            if criteria[0] == True:
                theNewLastCbg = newcbg
                break
            
            # total weight criterion -> new.tw() > last.tw()
            if criteria[1] == False:
                ##########################################################################
                if verbose:
                    print "# NOVEL lastTinyExon discarded; to low total weight"
                    print "#", newcbg._stopcodongraph
                ##########################################################################
                # continue -> no new tiny CBG
                continue

            # identity criterion -> allow a ratio i.s.o. new.id() > last.id()
            # this strict criterion (>) is applied for single-new-orf-CBGs
            if criteria[2] == False:
                ratio = newcbg.genetree().identity() / last.genetree().identity()
                if ratio < minimal_last_vs_new_identity_ratio:
                    ######################################################################
                    if verbose:
                        print "# NOVEL lastTinyExon discarded; to low identity"
                        print "#", newcbg._stopcodongraph, newcbg.genetree().identity()
                    ######################################################################
                    # continue -> no new tiny CBG
                    continue
 
            if criteria[3] == False:
                ##########################################################################
                if verbose:
                    print "# NOVEL lastTinyExon discarded; higher stopcodon2omsrdistance"
                    print "#", newcbg._stopcodongraph
                ##########################################################################
                # continue -> no new tiny CBG
                continue
 
            # if this point is reached, a new tiny last CBG has been found!
            theNewLastCbg = newcbg
            # break out of the for loop; store into the genestructure
            break



        # all okay -> ready for inserting the new CBG
        if theNewLastCbg and verbose:
            ################################################################################
            print "NEW FINAL TINY EXON FOUND!!"
            print theNewLastCbg
            print cbgIF, cbgIF.is_optimal(), cbgIF.is_acceptable()
            print cbgIF._optimal_aligned_donor, cbgIF.donor_phase()
            print cbgIF._optimal_aligned_acceptor, cbgIF.acceptor_phase()
            ################################################################################

        # hard-insert into the genestructure
        # using add_codingblock is likely to cause problems
        # because of the tinyness of the CBG
        if theNewLastCbg:
            for pos in range(0,len(self)):
                if self.codingblockgraphs[pos].IS_IGNORED: continue
                if self.codingblockgraphs[pos].IS_LAST:
                    thelast = self.codingblockgraphs[pos]
                    thelast.IS_LAST = False
                    newcbg.IS_LAST  = True
                    self.codingblockgraphs.insert(pos+1,theNewLastCbg)
                    # set the CBGInterface object in next and prev CBG
                    self.codingblockgraphs[pos]._CBGinterface3p = cbgIF
                    self.codingblockgraphs[pos+1]._CBGinterface5p = cbgIF
                    # break out; end of this function
                    break

            # done! return a True because newcbg is created & inserted
            return True
        else:
            # no newLastCbg found
            return False

    # end of function construct_final_tiny_cbg


    def split_final_cbg_on_spanningrange_difference(self,
        sprdif_min_aa_length=CBG_FINAL_SPRDIF_MIN_AA_LENGTH,
        sprdif_min_node_count=CBG_FINAL_SPRDIF_MIN_NODE_COUNT,
        sprdif_min_gtid_ratio=0.55,
        only_perform_if_stopcodon_tw_ratio_lte=CBG_FINAL_SPRDIF_ONLY_IF_STOP_TW_RATIO_LTE,
        only_preform_if_cbg_id_gte=CBG_FINAL_SPRDIF_ONLY_IF_CBG_ID_GTE ):
        """

        @type  sprdif_min_aa_length: integer
        @param sprdif_min_aa_length: minimal length of the sprdif in aa's
    
        @type  cbg_min_node_count: integer
        @param cbg_min_node_count: minimal number of nodes in a CBG to be elegiable for trying a split

        @type  sprdif_min_gtid_ratio: float
        @param sprdif_min_gtid_ratio:

        @type  only_perform_if_stopcodon_tw_ratio_lte: float 
        @param only_perform_if_stopcodon_tw_ratio_lte: run function only when lastCBG.stopcodongraph.totalweight <= threshold

        @type  only_preform_if_cbg_id_gte: float 
        @param only_preform_if_cbg_id_gte: run function only when lastCBG.genetree.identity() >= threshold

        """
        # get the CBG that is labelled as IS_LAST=True
        current_last = self.get_final_cbg()

        # check if we are alowed to peform this function
        # for groups of genes with very low identity, this function
        # is more likely to decrease the result then to improve the result

        # make AlignedStopCodonGraph
        current_last.align_stop_codons()
        tw_current = current_last._stopcodongraph.total_weight()
        ratio = tw_current / self.EXACT_SG_EDGE_COUNT

        # now check if it is alowed to enter the function: only_perform_if_...
        if only_perform_if_stopcodon_tw_ratio_lte and ratio > only_perform_if_stopcodon_tw_ratio_lte:
            return False
        if only_preform_if_cbg_id_gte and current_last.genetree().identity() < only_preform_if_cbg_id_gte:
            return False
            

        # check for rigth sprdif op requested size; if not => return False
        if not current_last.has_rigth_spanningrange_difference(
                sprdif_min_aa_length=sprdif_min_aa_length,
                sprdif_min_node_count=sprdif_min_node_count):
            # no rigth spanningrange difference -> done & return
            return False

        # make a deepcopy and clear cache of the one that will be processed
        last = deepcopy(current_last)
        last.clear_cache()

        # iteratively split
        splits = last.iteratively_split_codingblock_on_spanningrange_difference(
                side='rigth',
                sprdif_min_aa_length=sprdif_min_aa_length,
                sprdif_min_node_count=sprdif_min_node_count,
                )

        # was the split succesfull?
        if len(splits) == 1:
            # no splits => done here!
            return False

        # when here process the sprdif CBGs
        # 1) cbghmmsearch2pacbpcollection
        # 2) pacbpCollection2acceptedcbgs
        all_accepted_cbgs = []
        # loop over the splits; except for the most left one (the input `last` CBG)
        for splittedCBG in splits[1:]:
            if splittedCBG.__class__.__name__ == 'LowSimilarityRegionCodingBlockGraph': continue
            # complete with cbghmmsearch2pacbpcollection

            # get ratio of the GTG of this CBG
            ratio = splittedCBG.genetree().identity() / current_last.genetree().identity()

            # if ratio is bad -> do not perform!
            if sprdif_min_gtid_ratio and ratio < sprdif_min_gtid_ratio:
                continue

            pacbpCollection = cbghmmsearch2pacbpcollection(splittedCBG,self.input,
                    prev=last,
                    pacbp_min_length=sprdif_min_aa_length,
                    hmmsearch_num_hits=3
                    )

            # get list of accepted CBGs
            accepted =  conversion.pacbpCollection2AcceptedCodingBlockGraphs(pacbpCollection,prev=last)
            all_accepted_cbgs.extend( accepted )

        # if no accepted ones -> return False
        if not all_accepted_cbgs: return False

        # order graphs by total weight
        all_accepted_cbgs = ordering.order_graphlist_by_total_weight(all_accepted_cbgs)
        # and re-order on node occurrence: if a neighboring node is incorporated -> more likely!
        all_accepted_cbgs = ordering.reorder_cbgs_on_node_occurrence(all_accepted_cbgs,prev=last)

        # and now try to add the accepted cbgs into the genestructure
        # speedup the process by creating a tinyGSG object of only the last CBG
        # but, set the _GENETREE attribute to the genetree of the main GSG
        from graph_genestructure import GenestructureOfCodingBlockGraphs
        lastGSG = GenestructureOfCodingBlockGraphs(self.input)
        lastGSG.add_codingblock(current_last)
        lastGSG._GENETREE = self._GENETREE
        RETURN_STATUS_CBG_IS_ADDED = False

        for cbgL in all_accepted_cbgs:
            # only Ks CBG graphs are alowed here!
            if cbgL.node_count() != current_last.node_count(): continue

            if lastGSG.add_codingblock(cbgL,only_try_adding=True,
                max_cbg_gtg_topo_dif=self.MAX_CBG_FINAL_SPRDIF_GTG_TOPO_DIF,
                max_cbg_gtg_abs_dif=self.MAX_CBG_FINAL_SPRDIF_GTG_ABS_DIF,
                min_cbg_gtg_id_ratio=self.MIN_CBG_FINAL_SPRDIF_GTG_ID_RATIO,
                ):
                # it is addable; prepare for final addition to the genestructure
                lsrCBG = None
                cbgL.IS_SPLITTED     = False
                cbgL.IS_5P_SPLITTED  = False
                cbgL.IS_FIRST        = False
                cbgL.IS_LAST         = True
                current_last.IS_LAST = False
                # if identical nodes -> create a lsrCBG
                if not cbgL.node_set().difference(current_last.get_nodes()):
                    current_last.IS_SPLITTED    = True
                    current_last.IS_3P_SPLITTED = True
                    cbgL.IS_SPLITTED            = True
                    cbgL.IS_5P_SPLITTED         = True
                    lsrCBG = graphAbgp.codingblock_splitting.create_intermediate_lowsimilarity_region(
                            current_last, cbgL )
                    if not lsrCBG.node_count():
                        lsrCBG = None
                # now add the new last CBG
                status = self.add_codingblock(cbgL,
                        max_cbg_gtg_topo_dif=self.MAX_CBG_FINAL_SPRDIF_GTG_TOPO_DIF,
                        max_cbg_gtg_abs_dif=self.MAX_CBG_FINAL_SPRDIF_GTG_ABS_DIF,
                        min_cbg_gtg_id_ratio=self.MIN_CBG_FINAL_SPRDIF_GTG_ID_RATIO,
                        )
                status = lastGSG.add_codingblock(cbgL,
                        max_cbg_gtg_topo_dif=self.MAX_CBG_FINAL_SPRDIF_GTG_TOPO_DIF,
                        max_cbg_gtg_abs_dif=self.MAX_CBG_FINAL_SPRDIF_GTG_ABS_DIF,
                        min_cbg_gtg_id_ratio=self.MIN_CBG_FINAL_SPRDIF_GTG_ID_RATIO,
                        )
                # if added, update the return value (RETURN_STATUS_CBG_IS_ADDED)
                if status:
                    RETURN_STATUS_CBG_IS_ADDED = True
                    print cbgL
                    print cbgL.IS_5P_SPLITTED, cbgL.IS_SPLITTED, cbgL.IS_3P_SPLITTED
                # and add the intermediate lsrCBG when available
                if lsrCBG:
                    statusMainGSG = self.add_codingblock(lsrCBG)
                    statusLastGSG = lastGSG.add_codingblock(lsrCBG)
                    print "lsrCBG added:", statusMainGSG, statusLastGSG
            else:
                # not placeable in the genestructure
                pass

        # in exceptional cases, 2 CBGs can be added. In case the node_set() is identical,
        # yet another lsrCBG has to be created in between these 2 new CBGs
        # check this in the main GSG (NOT in the lastGSG; when a lsrCBG is added here,
        # splits are added to the surrounding CBGs. Because call-by-reference, these
        # splits are added to the main GSG (self) too, and adding the same lsrCBG
        # will fail (splitted CBGs are skipped!
        if RETURN_STATUS_CBG_IS_ADDED:
            self.finalize_genestructure()
            if self.join_false_inframe_introns():
               print "EXTRA lsrCBG added!!"
            # recreate interfaces if there is a new one created
            self.create_cbginterfaces()


        # return the return status True|False
        return RETURN_STATUS_CBG_IS_ADDED

    # end of function split_final_cbg_on_spanningrange_difference


    def final_codingblock_analyses(self,shift_max_cbg_number=2,
        ignore_cbg_max_aa_size=MAX_CBG_REMOVE_NONSENSE_FINAL_AA_LENGTH,
        ignore_cbg_min_identity_ratio=MIN_CBG_REMOVE_NONSENSE_FINAL_GTG_ID_RATIO,
        ignore_cbg_gtg_topo_dif=MAX_CBG_REMOVE_NONSENSE_FINAL_GTG_TOPO_DIF,
        ignore_cbg_gtg_abs_dif=MAX_CBG_REMOVE_NONSENSE_FINAL_GTG_ABS_DIF,
        ignore_cbg_gtg_tcode=MIN_CBG_REMOVE_NONSENSE_FINAL_TCODE_OMSR,
        verbose=False):
        """
        Ignore the final X CBGs that are not likely CBGs
        """
        # when only a single CBG -> no deletion possible
        if len(self) == 1: return False

        # list of (temporarily) deleted CBGs
        deleted = []

        ####################################################################
        if verbose: print "ENTERING final_codingblock_analyses FUNCTION !!"
        ####################################################################

        for iters in range(len(self)-1,max([0,len(self)-1-shift_max_cbg_number]),-1):
            # get the currently final CBG
            thiscbg = self.codingblockgraphs[-1]
            cbgif_is_optimal     = thiscbg._CBGinterface5p.is_optimal()
            cbgif_is_compatible  = thiscbg._CBGinterface5p.is_compatible()
            hmmsearch_in_sources = 'hmmsearch' in thiscbg.sources().keys()
            stopgra_is_optimal   = thiscbg._stopcodongraph.is_optimal()

            if cbgif_is_optimal:
                # as soon as an is_optimal() cbgIF is encountered -> break
                break
            elif stopgra_is_optimal and cbgif_is_compatible:
                # fair enough ;-)
                break
            elif not hmmsearch_in_sources and cbgif_is_compatible:
                # no hmmsearch component and a compatible cbgIF -> break
                break
            else:
                pass

            ####################################################################
            if verbose: print "elegiable for deletion:\n", thiscbg 
            ####################################################################

            # if here, then try to optimize the cbgIF of this CBG
            status = thiscbg._CBGinterface5p.optimize(verbose=verbose)
            cbgif_is_optimal     = thiscbg._CBGinterface5p.is_optimal()
            cbgif_is_compatible  = thiscbg._CBGinterface5p.is_compatible()

            ####################################################################
            if verbose: print "after cbgIF.optimize():\n", thiscbg._CBGinterface5p        
            ####################################################################

            # second try to keep this CBG
            if cbgif_is_optimal:
                # as soon as an is_optimal() cbgIF is encountered -> break
                break
            elif stopgra_is_optimal and cbgif_is_compatible:
                # fair enough ;-)
                break
            elif not hmmsearch_in_sources and cbgif_is_compatible:
                # no hmmsearch component and a compatible cbgIF -> break
                break
            else: 
                pass

            # check thiscbg for several properties; only remove it when
            # all criteria are met.
            # Get the genetree for comparison
            thegtg  = self.genetree()
            cbggtg  = thiscbg.genetree()
            prevcbg = self.codingblockgraphs[-2]
            checks = []
            checks.append( cbgif_is_compatible == False )
            checks.append( thiscbg.omsr_tcode_score() < ignore_cbg_gtg_tcode )
            checks.append( thiscbg.omsrlength() <= ignore_cbg_max_aa_size )
            checks.append( thegtg.graphalignmentdifference( cbggtg ) > ignore_cbg_gtg_topo_dif )
            checks.append( thegtg.absolutegraphalignmentdifference( cbggtg ) > ignore_cbg_gtg_abs_dif ) 
            checks.append( ( cbggtg.identity() / thegtg.identity() ) < ignore_cbg_min_identity_ratio )
            checks.append( len(thiscbg.mutual_nodes(prevcbg)) >= 1 )

            ####################################################################
            if verbose: print "ROUND", iters, checks
            ####################################################################

            # If all are False, no reason to believe this one is bogus
            # If all are True,  absolutely a bogus one.
            # Threshold for deletion when  at least 2 Trues.
            if checks.count(True) >= 2:
                pass
            else:
                # do not allow deletion of this CBG -> break
                break 

            # current IS_LAST not likely -> remove it!
            deleted.insert( 0, self.codingblockgraphs.pop() )
            # update the new IS_LAST cbg
            prevcbg.IS_LAST = True
            prevcbg._CBGinterface3p = None

            # try novel tiny tss construction
            if not prevcbg._stopcodongraph.is_optimal():
                IS_CREATED = self.construct_final_tiny_cbg(verbose=verbose)
            else:
                IS_CREATED = False

            ####################################################################
            if verbose:
                print "ROUND", iters, "of first_codingblock_analyses"
                print "deleted:", deleted[0]
                newlast = self.get_final_cbg()
                print "current:", newlast
                print "current:", newlast._stopcodongraph, 
                print newlast._stopcodongraph.is_optimal()
                print "created?", IS_CREATED
            ####################################################################

        # end of this effort. Is the IS_LAST now optimal?
        if deleted and self.get_final_cbg()._stopcodongraph.is_optimal():
            # CBGs removed and maybe novel tiny final CBG created
            return True 
        elif not deleted:
            # failed on the first one -> nothing to restore
            return False
        else:
            ####################################################################
            if verbose:
                print "DELETION FAILED ON HINDSIGHT"
                for item in deleted: print item
                print "RESTORING ORIGINAL LAST CBG"
            ####################################################################
            # failed; reset the removed CBGs to the que of CBGs
            self.get_final_cbg()._CBGinterface3p = deleted[0]._CBGinterface5p
            self.get_final_cbg().IS_LAST=False
            # actualy reset to the end of the que of the GSG
            for item in deleted:
                self.codingblockgraphs.append(item)
                # set to False, even the first one (in case >1 cbg in deleted)
                self.codingblockgraphs[-1].IS_LAST=False
            # and set IS_LAST back to the final one..
            self.codingblockgraphs[-1].IS_LAST = True
            # (Re)create cbgIFs
            self.create_cbginterfaces()
            # return status False
            return False

    # end of function final_codingblock_analyses



# end of class FinalCodingBlockGraphFunctions

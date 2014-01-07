"""
Class with functions for finding/optimizing the first CBG
in a GeneStructureOfCodingBlocks used in Alignment Based Gene Prediction
"""

# Module metadata
__authors__ = "Ate van der Burgt"
__license__ = "MIT"

# grapgAbgp Imports
import graphAbgp
from exceptions import *
from graph_exoncollection import *
import codingblock_splitting
import ordering
from subclass_sitealignment import sort_by_cumulative_score
import conversion

# Gene Imports
from gene.codingblock import CodingBlockStart, CodingBlockEnd
from gene.exon import FirstExonOnOrf

# ABGP Imports
from lib_shortleadingexononorf import find_leading_exon_on_orf
from lib_tinyexononorf import scan_orf_for_tiny_exon
from lib_hmm import cbgtinytssexonhmmsearch
from codingblockgraphinterface import CodingBlockGraphInterface
from lib_stopwatch import StopWatch

# Python Imports
from copy import deepcopy

# Global variables
from settings.genestructure import SHORT_LEADINGEXON_MAX_NT_LENGTH, SHORT_LEADINGEXON_MAX_NT_LENGTH, SHORT_LEADINGEXON_MAX_INTRON_NT_LENGTH  
from settings.codingblockgraph import CBG_MIN_AA_LENGTH, CBG_FIRST_SPRDIF_MIN_AA_LENGTH, CBG_FIRST_SPRDIF_MIN_NODE_COUNT
from settings.codingblockgraph import (
    MAX_CBG_REMOVE_NONSENSE_FIRST_AA_LENGTH,
    MAX_CBG_REMOVE_NONSENSE_FIRST_GTG_TOPO_DIF,
    MAX_CBG_REMOVE_NONSENSE_FIRST_GTG_ABS_DIF,
    MIN_CBG_REMOVE_NONSENSE_FIRST_GTG_ID_RATIO,
    MIN_CBG_REMOVE_NONSENSE_FIRST_TCODE_OMSR,
) 

class FirstCodingBlockGraphFunctions:
    """ """
    def get_first_cbg(self):
        """ 
        Get the CBG that is labeled as IS_FIRST
                
        @rtype:  CodingBlockGraph
        @return: CodingBlockGraph instance
        """                 
        # get current first cbg
        for cbg in self.codingblockgraphs:
            if cbg.IS_FIRST: 
                return cbg   
        return None         

    # end of function get_first_cbg


    def get_second_cbg(self):
        """ 
        Get the CBG that is directly after the one that is labeled as IS_FIRST
                
        @rtype:  CodingBlockGraph
        @return: CodingBlockGraph instance
        """                 
        # get current first cbg
        for pos in range(0,len(self)-1):
            if self.codingblockgraphs[pos].IS_FIRST:
                return self.codingblockgraphs[pos+1]
        return None         

    # end of function get_second_cbg


    def set_first_cbg(self,cbg):
        """
        Set this CBG in stead of the one that is labeled as IS_FIRST
        
        @type  cbg:  CodingBlockGraph
        @param cbg: CodingBlockGraph instance
        """
        # get current first cbg
        for pos in range(0,len(self)):
            if self.codingblockgraphs[pos].IS_FIRST: 
                self.codingblockgraphs[pos] = cbg
                break
           
    # end of function set_first_cbg


    def is_firt_cbg_annotated_first_cbg(self):
        """
        Does the (current) first CBG correspond to all (annotated) first exons?

        @rtype:  Boolean
        @return: True only if all (annotated) first exons are in the first CBG
        """
        # no CBGs in the GSG yet!
        if len(self) == 0: return False

        # get set of nodes of first codingblock in the list
        # do NOT use IS_FIRST attribute -> probably not set yet!
        firstCBGnodeset = self.codingblockgraphs[0].node_set()
        # generate a Set of first orf ids for the given gene structures
        annotatednodeset = Set()
        for org in self.input.keys():
            if len(self.input[org]['orfid-genestructure']) >= 1: 
                annotatednodeset.add(
                    (org, self.input[org]['orfid-genestructure'][0] )
                    )
        # check if all as first Exons annotated nodes are covered
        if len(annotatednodeset) == len(firstCBGnodeset) and\
        len(firstCBGnodeset.difference(annotatednodeset)) == 0:
            return True
        else:
            return False

    # end of function is_first_cbg_annotated_first_cbg


    def construct_first_tiny_cbg(self,
        max_exon_nt_length=SHORT_LEADINGEXON_MAX_NT_LENGTH,
        max_intron_nt_length=SHORT_LEADINGEXON_MAX_INTRON_NT_LENGTH,
        elegiable_acceptor_omsr_nt_offset=21,verbose=False):
        """
        Make a tiny first CBG by ``shooting tiny TSSexons into the deep``
        """
        first = self.get_first_cbg()
        ################################################################################
        if verbose:
            print "firstCBG:", first
            first.printmultiplealignment()
            if first._startcodongraph.alignedsites:
                print first._startcodongraph.alignedsites[0]
                if first._startcodongraph.alignedsites[0].__class__.__name__ == 'AlignedTranslationalStartSiteGraph':
                    print first._startcodongraph.alignedsites[0].is_optimal()
                else:
                    # not a aligned site but a collection
                    pass

            #### do acceptor alignment; just for printing
            ###first.harvest_elegiable_acceptor_sites()
            ###from subclass_sitealignment import sort_by_cumulative_score
            ###import conversion
            ###first._spliceacceptorgraph.alignedsites = sort_by_cumulative_score(
            ###    [ conversion.SpliceSiteCollectionGraph2AlignedSpliceSiteGraph(
            ###        algsite,max_node_count=first.organism_set_size()
            ###        ) for algsite in first._spliceacceptorgraph.find_conserved_sites() ] )
            ###print "ALIGNED ACCEPTORS:"
            ###print first._spliceacceptorgraph.alignedsites[0]

        ################################################################################

        omsr = first.overall_minimal_spanning_range()
        ECG = ExonCollectionGraph()
        for organism in first.organism_set():
            node = first.node_by_organism(organism)
            theorf = first.get_orfs_of_graph(organism=organism)[0]
        
            # check if this orf can serve as a first tiny exon part
            if theorf.has_start_upstream_of( min(omsr[node])+(elegiable_acceptor_omsr_nt_offset/3) ):
                if verbose: print organism,theorf.id,"has start upstream of omsr",min(omsr[node])
                cbe = CodingBlockEnd( theorf.aapos2dnapos( min(omsr[node]) ) )
                # set pssm_score to (very) high; this rewards
                # using the current Orf as the first Orf
                cbe.pssm_score = 20.0
                theorf.scan_orf_for_pssm_tss() 
                for tss in theorf._tss_sites:
                    # if TSS AA-pos < min(omsr)-AA-pos -> make a FIOO object!
                    if (tss.pos/3) < min(omsr[node]):
                        fioo = FirstExonOnOrf(tss,cbe,theorf)
                        fioonode = (organism,theorf.id,fioo.start,fioo.end)
                        ECG.add_node_and_object(fioonode,fioo)
                        if verbose: print "TSS in existing orf", theorf.id, tss.pos, min(omsr[node]), "omsrAA"
                    elif (tss.pos/3) < min(omsr[node])+(elegiable_acceptor_omsr_nt_offset/3):
                        # this can happen in cases like the following, examplified by the
                        # following (correct protein sequence) alignment
                        # FOXG_08492      MAPATQYELSPPPTDAVSAIAFAP... orfs 102 - 102 (! inframe intron!)
                        # FVEG_06184      MAPATQYELSPPPTDAVSAIAFAP... orfs  12 -  56
                        # FGSG_06219      MAPATQYELSPSPTDAVSSIAFAP... orfs  18 -  72
                        # NCU09744        MAPSTQFEVSQPPSDAISAVIFAP... orfs  59 -  23
                        # MGG_06111       ----kMFEAEPAPGDCPTAMKFAP... orfs  NA -  49
                        #                       :* . .* *. ::: ***
                        #                     ==================== OMSR of CBG (102-56-72-23-49)
                        # The pre-final CBG (102-56-72-23-49) has the Methionine of MGG in its OMSR
                        # Therefor, it is missed in the FEG and as a consequence in the `newcbg`
                        # when we do not look a few AAs into the OMSR region
                        altcbe = CodingBlockEnd( tss.pos+3 )
                        altcbe.pssm_score = 20.0
                        fioo = FirstExonOnOrf(tss,altcbe,theorf)
                        fioonode = (organism,theorf.id,fioo.start,fioo.end)
                        ECG.add_node_and_object(fioonode,fioo)
                        if verbose: print "TSS in existing orf AFTER OMSR!!", theorf.id, tss.pos, min(omsr[node]), "omsrAA"
                    else:
                        break

            # calculate an offset for the acceptor position
            # variable elegiable_acceptor_omsr_nt_offset is needed to
            # enlarge the OMSR definded offset. When the OMSR is by chance
            # a few nt or aa larger than the actual exon length, the true
            # acceptor position can be erroneously abandoned.
            offset = min(omsr[node]) * 3 + elegiable_acceptor_omsr_nt_offset

            # get elegiable orfs based on offset and intron length
            orflist = self.input[organism]['orfs'].get_elegiable_orfs(
                    max_orf_start=offset,
                    min_orf_end=offset-max_intron_nt_length,
                    has_starts=True )

            ####################################################################
            if verbose: print organism, "pot.orfs:", [orf.id for orf in orflist]
            ####################################################################

            # tmp store node and objects in a dict; we only allow a SINGLE,
            # highest scoring TSS exon of length <= 8. This because at this stage
            # quit easily a few of these tiny TSS options occur, which can - by chance -
            # result in wrong asignment of the actual TSS exon

            node2object = {}
            for orf in orflist:
                results = find_leading_exon_on_orf(
                        orf,theorf,
                        subsequent_acceptor_pos=offset,
                        max_leadingexon_nt_length=max_exon_nt_length,
                        max_leadingexon_intron_nt_length=max_intron_nt_length,
                        )
                # loop over the exons and store them; check for really short ones 
                # that only a single one is incorporated
                for exon,intron in results:
                    node = (organism,orf.id,exon.start,exon.end)
                    if not node2object:
                        node2object[node] = exon
                    elif exon.length > 8:
                        node2object[node] = exon
                    elif min([ ex.length for ex in node2object.values()]) > 8:
                        node2object[node] = exon
                    else:
                        # check if this one is higher scoring than already in the list
                        exonlist = ordering.order_list_by_attribute(node2object.values(),'length')
                        if exonlist[0].pssm_score >= exon.pssm_score:
                            pass  # do not store current <=8 nt exon
                        else:
                            # this new exon is higher! remove the old one!
                            replacethisexon = exonlist[0]
                            replacenode = None
                            for thisnode,storedexon in node2object.iteritems():
                                if str(storedexon) == str(replacethisexon):
                                    replacenode = thisnode
                                    ############################################
                                    if verbose:
                                        print "tiny TSS DISCARDED!", thisnode,
                                        print "len:%s, pssm:%2.2f" % (storedexon.length,
                                             storedexon.pssm_score)
                                    ############################################
                                    break
                            # remove the old, store the new
                            del( node2object[replacenode] )
                            node2object[node] = exon

            # store the small TSS exon nodes to the ECG
            for node,exon in node2object.iteritems():
                ECG.add_node_and_object(node,exon)
                ################################################################
                if verbose: print " ", node, "\t%s\t%1.2f\t%1.2f\t(%s)\t%s" % (
                        ("   %s" % exon.length)[-4:], exon.acceptor.pssm_score,
                        exon.donor.pssm_score, exon.donor.phase,
                        exon.proteinsequence() )
                ################################################################
        if verbose: print ECG

        # only continue if all organisms are represented in the ECG
        if first.organism_set_size() > ECG.organism_set_size():
            if verbose: print "%s > %s, not all organisms covered -> return False" % (
                    first.organism_set_size(), ECG.organism_set_size() )
            return False


        # create edges in the ECG
        ECG.create_edges()
        if verbose: print ECG
        if verbose: "edges created:", ECG.edge_count(), ECG.node_count()

        for edges in range(first.node_count()-1,1,-1):
            # search for complete graphs in this
            first_exon_graphs = ECG.find_fully_connected_subgraphs(edges=edges,verbose=verbose)
            if verbose:
                print "EDGES:", edges, len(first_exon_graphs), "first exon graphs, printing first 10"
                print "\n".join( [ "feg: %s" % feg for feg in first_exon_graphs[0:10] ] )
            # quit as soon as there are Aligned Exon Graphs found
            if first_exon_graphs: break
        else:
            # no first_exon_graphs found at all...
            first_exon_graphs = []

        # only continue if there is an perfectly aligned first exon graph
        if not first_exon_graphs:
            if verbose: print "no first exon graphs -> return False"
            return False

        # convert to CodingBlockGraphs
        new_first_cbgs = []
        _cbg2feg = {} # dict that checks which cbg originates from which feg
        for feg in first_exon_graphs[0:3]:
            if list(Set([ obj.donor.phase for obj in feg._node_object.values() ])) == [ None ]:
                # this is the extention of the current firstCBG itself (all parts
                # of existing Orfs, not a small exon joined with an intron).
                # Ignore this one -> because it is not extra information.
                # Moreover, ExonCollectionGraph2CodingBlockGraph will crash ;-)
                continue
            cbg = ExonCollectionGraph2CodingBlockGraph(feg,is_first=True,firstCBG=first)
            if cbg != None and cbg != False:
                new_first_cbgs.append( cbg )
                _cbg2feg[cbg] = feg

        if not new_first_cbgs:
            if verbose: print "no first exon graph convertable to new CBG -> return False"
            return False

        # order by total weight
        new_first_cbgs = ordering.order_graphlist_by_total_weight(new_first_cbgs)
        if verbose:
            print "FEGs (%s) converted to CBGs (%s)" % ( len(first_exon_graphs), len(new_first_cbgs) )
            for cbg in new_first_cbgs: print "cbg:", cbg

        # the first ECG is the best ECG; get this one and the FEG it originates from
        newcbg         = new_first_cbgs[0]
        firstExonGraph = _cbg2feg[newcbg]

        if newcbg.organism_set_size() != first.organism_set_size():
            if verbose: print "not all organisms available in feg CBG %s -> return False" % newcbg.organism_set_size()
            return False

        ################################################################################
        if verbose:
            print "NOVEL feg CBG found!!!"
            print newcbg
            print newcbg._splicedonorgraph.alignedsites[0]
            print newcbg._startcodongraph.alignedsites[0]
            newcbg.printmultiplealignment()
            print "haveallstarts()", newcbg.have_all_starts(),
            print "upstreamomsrstarts()", newcbg.have_all_starts_upstream_of_omsr()
        ################################################################################

        # final check; can we make an is_compatible CBGinterface for this novel tinyexon?
        cbgIF = CodingBlockGraphInterface(newcbg,first)

        # now harvest the spicesites
        cbgIF.harvest_splice_sites()

        # Check for which organisms there are introns allowed/needed
        # This is defined by the FEG donor boundary (donor or CodingBlockEnd) 
        distinct_orgs = []
        for node in firstExonGraph.get_nodes():
            _org = firstExonGraph.organism_by_node(node)
            exon = firstExonGraph.get_node_object(node)
            if exon.donor.__class__.__name__ == 'SpliceDonor':
                distinct_orgs.append( firstExonGraph.organism_by_node(node) )
                if cbgIF._interface_is_intron[_org] in [False,None]:
                    # interface_is_intron incorrectly interpreted by
                    # the CodingBlockGraphIntrerface code. Overrule it here!
                    print "OVERRULE _interface_is_intron", _org, "from: ", cbgIF._interface_is_intron[_org], "to:", True
                    cbgIF._interface_is_intron[_org] = True

        # apply the allowed introns per organism to the cbgIF
        cbgIF.allow_intron_in_organisms(distinct_orgs)
        cbgIF.find_conserved_splice_sites()
        ################################################################################
        if verbose:
            print "Distinct_orgs with donors, non-extensions:", distinct_orgs
            print cbgIF
            cbgIF.interfaceproperties()
        ################################################################################

        # okay, now this piece of code is putatively dangerous
        # the cbgIF is tried to *optimize*. This potentially
        # can cause FP first exons....
        if not cbgIF.is_optimal():
            cbgIF.optimize()
            if verbose: print cbgIF

        # now asses the predicted interface
        if not cbgIF.is_compatible():
            print "# NOVEL firstTinyExon discarded because of incompatible CBGinterface"
            print newcbg
            print cbgIF, cbgIF.is_optimal(), cbgIF.is_acceptable()
            print cbgIF._splicedonorgraph
            print cbgIF._optimal_aligned_donor, cbgIF.donor_phase()
            print cbgIF._optimal_aligned_acceptor, cbgIF.acceptor_phase()
            print cbgIF._spliceacceptorgraph
            # return False -> no new tiny CBG
            return False
        else:
            # allo okay -> ready for inserting the new CBG
            if verbose:
                ################################################################################
                print "NEW FIRST TINY TSS EXON FOUND!!"
                print newcbg
                print cbgIF, cbgIF.is_optimal(), cbgIF.is_acceptable()
                print cbgIF._optimal_aligned_donor, cbgIF.donor_phase()
                print cbgIF._optimal_aligned_acceptor, cbgIF.acceptor_phase()
                ################################################################################

        # and hard-insert into the genestructure
        # using add_codingblock is likely to cause problems
        # because of the tinyness of the CBG
        for pos in range(0,len(self)):
            if self.codingblockgraphs[pos].IS_IGNORED: continue
            if self.codingblockgraphs[pos].IS_FIRST:
                first = self.codingblockgraphs[pos]
                first.IS_FIRST = False
                # prepare inserting the new created cbg
                newcbg.IS_FIRST = True
                self.codingblockgraphs.insert(pos,newcbg)
                # set the CBGInterface object in next and prev CBG
                self.codingblockgraphs[pos]._CBGinterface3p = cbgIF
                self.codingblockgraphs[pos+1]._CBGinterface5p = cbgIF
                # break out; end of this function
                break

        # done! return a True because newcbg is created & inserted
        return True

    # end of function constuct_first_tiny_cbg


    def split_first_cbg_on_spanningrange_difference(self,
        sprdif_min_aa_length=CBG_FIRST_SPRDIF_MIN_AA_LENGTH,
        sprdif_min_node_count=CBG_FIRST_SPRDIF_MIN_NODE_COUNT,
        verbose=False):
        """                                   
        TODO: implement microGSG just as in split_final_cbg_on_... 
        """                                   
        # get the CBG that is labelled as IS_FIRST=True
        current_first = self.get_first_cbg()

        ################################################################
        if verbose:
           print "SplFirstCbgSprDif:",
           print sprdif_min_aa_length, sprdif_min_node_count
           print current_first._startcodongraph
           if current_first._startcodongraph.alignedsites:
               print current_first._startcodongraph.alignedsites[0]
           print current_first
           current_first.printmultiplealignment()
           for pacbporf in current_first.pacbps.values():
               print pacbporf
        ################################################################
            
        # check for left sprdif op requested size; if not => return False
        if not current_first.has_left_spanningrange_difference(
                sprdif_min_aa_length=sprdif_min_aa_length,
                sprdif_min_node_count=sprdif_min_node_count):
            # no left spanningrange difference -> done & return
            return False

        # make a deepcopy and clear cache of the one that will be processed
        if verbose: print current_first
        first = deepcopy(current_first)
        first.clear_cache()
        ################################################################
        if verbose: print "NoCache:", first
        ################################################################

        # iteratively split
        splits = first.iteratively_split_codingblock_on_spanningrange_difference(
                side='left',
                sprdif_min_aa_length=sprdif_min_aa_length,
                sprdif_min_node_count=sprdif_min_node_count,
                )

        ################################################################
        if verbose: print "/n".join(["sprdif:%s" % spl for spl in splits])
        ################################################################

        # was the split succesfull? If not return False and done
        if len(splits) == 1: return False

        # continue with the first non-lsrCBG split from the rigth
        # after the original (most rigth) one
        if len(splits) == 2:
            cbgMostLeft = splits[0]
        else:
            for splitpos in range(len(splits)-2,-1,-1):
                if splits[splitpos].__class__.__name__ == 'CodingBlockGraph':
                    cbgMostLeft = splits[splitpos]
                    break

        ################################################################
        if verbose:
            print "elegiable left split:",cbgMostLeft
            cbgMostLeft.printmultiplealignment()
        ################################################################

        # does this left splitted CBG have all starts upstream of OMSR
        # i.e. is it an elegiable first CBG?
        if not cbgMostLeft.have_all_starts(): return False

        # try to complete this splitted CBG asa tinyexonhmm
        cbglist = cbgtinytssexonhmmsearch(cbgMostLeft,self.input,next=first,verbose=verbose)

        TSSEXON_IS_ADDED = False
        for cbgL in cbglist:
            if cbgL.node_count() == self.EXACT_SG_NODE_COUNT and self.add_codingblock(cbgL,log_debug=True,
                max_cbg_gtg_topo_dif=self.MAX_CBG_FIRST_SPRDIF_GTG_TOPO_DIF,
                max_cbg_gtg_abs_dif=self.MAX_CBG_FIRST_SPRDIF_GTG_ABS_DIF,
                min_cbg_gtg_id_ratio=self.MIN_CBG_FIRST_SPRDIF_GTG_ID_RATIO,
                only_try_adding=True):

                ####################################################################
                if verbose:
                    # make TranslationalStartSiteCollectionGraph
                    cbgL.harvest_elegiable_tss_sites()
                    # do site alignment of TSS 
                    cbgL._startcodongraph.collection2alignedsites()
                    print "SPRDIFtssCBG:", cbgL
                    print "SPRDIFtssCBG:", cbgL._startcodongraph
                    if cbgL._startcodongraph and cbgL._startcodongraph.alignedsites:
                        print "SPRDIFtssCBG:", cbgL._startcodongraph.alignedsites[0]
                    else:
                        print "SPRDIFtssCBG:", None
                ####################################################################

                # it is addable! 
                # (1) create a CBGInterface; ignore if not is_compatible
                # (2) update IS_SPLITTED attributes / create lsrCBGs
                # (3) recreate cache
                # (4) actually add_codingblock()
                # (5) _update_is_first_cbg
                distinct_orgs = cbgL.organisms_with_different_orfs(current_first)
                if distinct_orgs:
                    # At least one organism/gene has a novel orf in cbgL 
                    # Create a cbgIF between cbgL -- current_first and
                    # check if it is a is_compatible interface  (NOT is_optimal!)
                    lsrCBG = None
                    cbgL.IS_SPLITTED    = False
                    cbgL.IS_3P_SPLITTED = False
                    ####################################################################
                    if verbose: print "potential new cbgL:", cbgL
                    ####################################################################
                    cbgIF = CodingBlockGraphInterface(cbgL,current_first)
                    cbgIF.harvest_splice_sites()
                    ####################################################################
                    if verbose: print "new cbgLIF:", cbgIF
                    ####################################################################
                    cbgIF.allow_intron_in_organisms(distinct_orgs)
                    ####################################################################
                    if verbose: print "new cbgLIF:", cbgIF, "non-introns removed"
                    ####################################################################
                    cbgIF.find_conserved_splice_sites()
                    ####################################################################
                    if verbose: print "new cbgLIF:", cbgIF, "sites aligned"
                    ####################################################################
                    # check the status of the cbgIF; -> continue if False
                    if not cbgIF.is_compatible(): continue
                    # set the interface in the cbgL object
                    cbgL._CBGinterface3p = cbgIF

                else:
                    # all nodes/orfs identical -> extention of the current_first CBG

                    # only accept this extension if the AlignedTSS is better as the current one
                    # this because otherwise this function is likely to overpredict things
                    # make TranslationalStartSiteCollectionGraph
                    cbgL.harvest_elegiable_tss_sites()
                    # do site alignment of TSS
                    cbgL._startcodongraph.collection2alignedsites()

                    currentTSS, newTSS = None, None
                    if current_first._startcodongraph and current_first._startcodongraph.alignedsites:
                        currentTSS = current_first._startcodongraph.alignedsites[0]
                    if cbgL._startcodongraph and cbgL._startcodongraph.alignedsites:
                        newTSS = cbgL._startcodongraph.alignedsites[0]
                    if not newTSS:
                        # no AlignedTSS in the cbg extention
                        continue
                    elif not currentTSS:
                        # no AlignedTSS in the starting situation, and now there is one!
                        pass
                    # both CBGs have an AlignedTSS. Which one is best?
                    elif currentTSS.node_count() == self.EXACT_SG_NODE_COUNT and\
                    newTSS.node_count() < self.EXACT_SG_NODE_COUNT:
                        # novelTSS lacks nodes -> discard
                        continue
                    elif newTSS.node_count() >= currentTSS.node_count() and\
                    newTSS.cumulative_score() > currentTSS.cumulative_score():
                        # better score for the novelTSS
                        pass
                    else:
                        # all other cases -> a deterioration in AlignedTSS. Discard
                        continue 

                    # create a lsrCBG, and do not check the cbgIF (lsrCBG IF always is_optimal!)
                    lsrCBG = None
                    lsrCBG = codingblock_splitting.create_intermediate_lowsimilarity_region(
                        cbgL, current_first )
                    if not lsrCBG or not lsrCBG.node_count():
                        # although this is an extention of the current_first CBG,
                        # lsrCBG creation failed. Do not accept this novel tssCBG
                        continue
                    else:
                        # lsrCBG creation succesfull
                        cbgL.IS_SPLITTED             = True
                        cbgL.IS_3P_SPLITTED          = True
                        current_first.IS_SPLITTED    = True
                        current_first.IS_3P_SPLITTED = True


                # (3) (re)create the cache of the cbgL and set some other properties
                cbgL.clear_cache()
                cbgL.create_cache()
                cbgL.IS_FIRST = True
                cbgL.IS_LAST  = False


                # (4) try to add to the genestructure graph
                status = self.add_codingblock(cbgL,
                        max_cbg_gtg_topo_dif=self.MAX_CBG_FIRST_SPRDIF_GTG_TOPO_DIF,
                        max_cbg_gtg_abs_dif=self.MAX_CBG_FIRST_SPRDIF_GTG_ABS_DIF,
                        min_cbg_gtg_id_ratio=self.MIN_CBG_FIRST_SPRDIF_GTG_ID_RATIO,
                        )

                if status:
                    # yes, storing the new IS_FIRST cbgL succesfull!
                    TSSEXON_IS_ADDED = True
                    # set current_first.IS_FIRST to False!
                    current_first.IS_FIRST = False

                    # run _update_is_first_cbg
                    self._update_is_first_cbg()

                    # add the lsrCBG if it exists & create interfaces
                    if lsrCBG:
                        from genestructure_lsrcbgs import prepare_lsrcbg_and_cbg_for_gsg_insertion
                        prepare_lsrcbg_and_cbg_for_gsg_insertion(cbgL,lsrCBG)
                        prepare_lsrcbg_and_cbg_for_gsg_insertion(lsrCBG,current_first)
                        # store the lsrCBG into the GSG
                        statusLsrCBG = self.add_codingblock(lsrCBG)
                        if not statusLsrCBG:
                            print "SeriousWarning: lsrCBG not insertable in GSG\n", lsrCBG
                    else:
                        # copy the new cbgIF to the currrent_first 5p side
                        current_first._CBGinterface5p = cbgL._CBGinterface3p

                    ############################################################
                    if verbose:
                        print "ADDED:", status, cbgL
                        print cbgL.printmultiplealignment()
                        if lsrCBG: print "ADDED:", statusLsrCBG, lsrCBG
                    ############################################################


                else:
                    # cbgL not alowed in the GSG -> continue to the next cbg
                    continue

            else:
                pass

        # return the status of tssexon addition
        return TSSEXON_IS_ADDED

    # end of function split_first_cbg_on_spanningrange_difference


    def _update_is_first_cbg(self):
        """
        """
        self.codingblockgraphs[0]._CBGinterface5p = None
        self.codingblockgraphs[0].IS_FIRST = True
        firstcbg = self.get_first_cbg()
        firstcbg.harvest_elegiable_tss_sites()
        firstcbg._startcodongraph.collection2alignedsites()
        for atssgra in firstcbg._startcodongraph.alignedsites:
            atssgra._codingblockgraph = firstcbg

    # end of function _update_is_first_cbg 


    def first_codingblock_analyses(self,shift_max_cbg_number=2,
        ignore_cbg_max_aa_size=MAX_CBG_REMOVE_NONSENSE_FIRST_AA_LENGTH,
        ignore_cbg_min_identity_ratio=MIN_CBG_REMOVE_NONSENSE_FIRST_GTG_ID_RATIO,
        ignore_cbg_gtg_topo_dif=MAX_CBG_REMOVE_NONSENSE_FIRST_GTG_TOPO_DIF,
        ignore_cbg_gtg_abs_dif=MAX_CBG_REMOVE_NONSENSE_FIRST_GTG_ABS_DIF,
        ignore_cbg_tcode=MIN_CBG_REMOVE_NONSENSE_FIRST_TCODE_OMSR,
        verbose=False):
        """
        Ignore the first X CBGs that are not likely CBGs
        """
        if self.get_first_cbg().have_all_starts_upstream_of_omsr():
            return False
        elif len(self) == 1:
            return False
        elif 'LowSimilarityRegionCodingBlockGraph' in\
        [ cbg.__class__.__name__ for cbg in self.codingblockgraphs[1:shift_max_cbg_number+1] ]:
            return False
        else:
            pass

        # print message that function is entered in verbode mode
        if verbose: print "ENTERING first_codingblock_analyses FUNCTION !!"

        # list of (temporarily) deleted CBGs
        deleted = []

        for iters in range(0,min([shift_max_cbg_number,len(self)])):
            thiscbg = self.codingblockgraphs[0]
            if thiscbg.have_all_starts_upstream_of_omsr():
                # done! this should be the new IS_FIRST cbg
                break
            elif not 'hmmsearch' in thiscbg.sources().keys():
                # if no hmmsearch component in this CBG -> do not remove it
                break
            elif thiscbg._CBGinterface3p.is_optimal():
                # if an optimal CBGinterface -> do not remove it
                break
            else:
                pass

            # check thiscbg for several properties; only remove it when
            # all criteria are met.
            # Get the genetree for comparison
            thegtg  = self.genetree()
            cbggtg  = thiscbg.genetree()
            checks = []
            checks.append( thiscbg.omsr_tcode_score() < ignore_cbg_tcode )
            checks.append( thiscbg.omsrlength() <= ignore_cbg_max_aa_size )
            checks.append( thegtg.graphalignmentdifference( cbggtg ) > ignore_cbg_gtg_topo_dif )
            checks.append( thegtg.absolutegraphalignmentdifference( cbggtg ) > ignore_cbg_gtg_abs_dif ) 
            checks.append( ( cbggtg.identity() / thegtg.identity() ) < ignore_cbg_min_identity_ratio )
            if len(self) <= 1:
                # final CBG to check -> no mutual nodes analyses
                checks.append( False )
            else:
                # get nextcbg and check for nodes occuring in both CBGs
                nextcbg = self.codingblockgraphs[1]
                checks.append( len(thiscbg.mutual_nodes(nextcbg)) >= 1 )

            if True or verbose:
                print "ROUND", iters, checks

            # If all are **False**, no reason to believe this one is bogus
            # If all are **True **,  absolutely a bogus one.
            # Threshold for deletion when  at least 2 Trues.
            if checks.count(True) >= 2:
                pass
            else:
                # do not allow deletion of this CBG -> break
                break 

            # current IS_FIRST != have_all_starts_upstream_of_omsr()
            # so -> remove it!
            deleted.append( self.codingblockgraphs.pop(0) )

            # update current firstCBG to the new IS_FIRST
            self._update_is_first_cbg()

            # try novel tiny tss construction
            if not self.get_first_cbg().have_all_starts_upstream_of_omsr():
                IS_CREATED = self.construct_first_tiny_cbg(verbose=True)
            else:
                IS_CREATED = False

            ####################################################################
            if verbose:
                print "ROUND", iters, "of first_codingblock_analyses"
                print "deleted:", deleted[-1]
                print "current:", self.get_first_cbg()
                print "current:", self.get_first_cbg().have_all_starts_upstream_of_omsr()
                print "created?", IS_CREATED
            ####################################################################

        # end of this effort. Is the IS_FIRST now optimal?
        if self.get_first_cbg().have_all_starts_upstream_of_omsr():
            # CBGs removed and maybe novel tiny CBG created
            return True 
        elif not deleted:
            # failed on the first one -> nothing to restore
            return False
        else:
            # failed; reset the removed CBGs to the que of CBGs
            self.get_first_cbg()._startcodongraph = None
            self.get_first_cbg()._CBGinterface5p = deleted[-1]._CBGinterface3p
            self.get_first_cbg().IS_FIRST=False
            deleted.reverse()
            # actualy insert to the front of the que of the GSG
            for item in deleted:
                self.codingblockgraphs.insert(0,item)
                # set to False, even the first one (in case >1 cbg in deleted)
                self.codingblockgraphs[0].IS_FIRST=False
            # and set IS_FIRST back to the frontal one..
            self.codingblockgraphs[0].IS_FIRST = True
            # return status False
            return False

    # end of function first_codingblock_analyses


    def first_codingblock_analyses_remove_falsecbgs(self,shift_max_cbg_number=2,
        ignore_cbg_max_aa_size=MAX_CBG_REMOVE_NONSENSE_FIRST_AA_LENGTH,
        ignore_cbg_min_identity_ratio=MIN_CBG_REMOVE_NONSENSE_FIRST_GTG_ID_RATIO,
        ignore_cbg_gtg_topo_dif=MAX_CBG_REMOVE_NONSENSE_FIRST_GTG_TOPO_DIF,
        ignore_cbg_gtg_abs_dif=MAX_CBG_REMOVE_NONSENSE_FIRST_GTG_ABS_DIF,
        ignore_cbg_tcode=MIN_CBG_REMOVE_NONSENSE_FIRST_TCODE_OMSR,
        verbose=False):
        """
        Ignore the first X CBGs that are not likely CBGs
        """
        # correct shift_max_cbg_number for the length of the GSG
        shift_max_cbg_number = min([shift_max_cbg_number,len(self)-1])

        if len(self) == 1:
            return False
        elif 'LowSimilarityRegionCodingBlockGraph' in\
        [ cbg.__class__.__name__ for cbg in self.codingblockgraphs[1:shift_max_cbg_number+1] ]:
            return False
        else:
            pass

        # make list of booleans for weather or not there might be an
        # alternative firstCBG -> True must occur in this list
        cbgshavestarts = [ cbg.have_all_starts_upstream_of_omsr() for cbg in\
                            self.codingblockgraphs[0:shift_max_cbg_number+1] ]

        # check if CBGs with all Methiones exist; ignore current first
        if cbgshavestarts[1:].count(True) == 0:
            # no alternative more downstream CBG with methiones
            return False

        # Check this same frontal lists of CBGs for `correctness`.
        cbgscorrectness= []
        cbglengths     = []
        cbgomsrtcode   = []
        for cbg in self.codingblockgraphs[0:shift_max_cbg_number+1]:
            tmplist = []
            for org in cbg.organism_set():
                tmplist.append( self._codingblock_prediction_status(cbg,org) )
            if False in tmplist:
                cbgscorrectness.append( False )
            else:
                cbgscorrectness.append( True )
            cbglengths.append( cbg.omsrlength() )
            cbgomsrtcode.append( cbg.omsr_tcode_score() )

        # Only continue if there is an incorrect CBG among these
        if cbgscorrectness[0:shift_max_cbg_number].count(False) == 0:
            # no CBGs are listed as erroneous
            return False

        ########################################################################
        if verbose:
            firstCBG = self.get_first_cbg()
            print "first_codingblock_analyses_remove_falsecbgs:"
            print firstCBG
            print firstCBG._startcodongraph
            if firstCBG._startcodongraph:
                if firstCBG._startcodongraph.alignedsites:
                    print firstCBG._startcodongraph.alignedsites[0]
                else:
                    print None
            else:
               print None
            print cbgshavestarts
            print cbgscorrectness
            print cbglengths
            print cbgomsrtcode
            print "THRESHOLDS:"
            print "tcode <", ignore_cbg_tcode
            print "omsrlength <=", ignore_cbg_max_aa_size
        ########################################################################

        # now do a more elaborate analyses on this subset of frontal CBGs
        # We are in search for small, creapy CBGs in front of a strong TSS CBG
        cbganalyses = []
        for pos in range(0,shift_max_cbg_number+1):
            thiscbg = self.codingblockgraphs[pos]
            # check thiscbg for several properties; only remove it when
            # all criteria are met.
            # Get the genetree for comparison
            thegtg  = self.genetree()
            cbggtg  = thiscbg.genetree()
            checks = []
            checks.append( thiscbg.omsr_tcode_score() < ignore_cbg_tcode )
            checks.append( thiscbg.omsrlength() <= ignore_cbg_max_aa_size )
            checks.append( thegtg.graphalignmentdifference( cbggtg ) > ignore_cbg_gtg_topo_dif )
            checks.append( thegtg.absolutegraphalignmentdifference( cbggtg ) > ignore_cbg_gtg_abs_dif )
            checks.append( ( cbggtg.identity() / thegtg.identity() ) < ignore_cbg_min_identity_ratio )
            if pos == shift_max_cbg_number:
                # final CBG to check -> no mutual nodes analyses
                checks.append( False )
            else:
                # get nextcbg and check for nodes occuring in both CBGs
                nextcbg = self.codingblockgraphs[pos+1]
                checks.append( len(thiscbg.mutual_nodes(nextcbg)) == 0 )

            ####################################################################
            if verbose: print "ROUND", pos, checks
            ####################################################################

            # If all are **False**, no reason to believe this one is bogus
            # If all are **True **,  absolutely a bogus one.
            # Threshold for deletion when  at least 2 Trues.
            if checks.count(True) >= 2:
                # this is a shitty CBG
                cbganalyses.append( False )
            else:
                # this is a CBG that seems fine
                cbganalyses.append( True )

        ####################################################################
        if verbose: print "SUMMARIZED CBG ANALYSES:", cbganalyses 
        ####################################################################

        if not cbganalyses:
            # no CBG position remaining that can be omitted
            return False
        elif not True in cbganalyses:
            # Not a single 100% correct CBG.
            # Hopeless to improve this GSG
            return False
        elif [cbganalyses[0],cbgscorrectness[0],cbgshavestarts[0]].count(False) == 0:
            # current first CBG remaines suitable....
            return False
        elif not False in cbgscorrectness:
            # All CBGs are ~correct. Leave GSG intact
            return False
        else:
            pass

        # Judgement day, by example, shift_max_cbg_number == 2, GSG >= 3
        # cbgshavestarts    [True, False, True]
        # cbgscorrectness   [True, False, True]
        # cbganalyses       [False, False, True]
        # delete just before the first all Trues in the 3 lists
        deleteuntillpos = 0
        for pos in range(0,shift_max_cbg_number+1):
            if [cbganalyses[pos],cbgscorrectness[pos],cbgshavestarts[pos]].count(False) == 0:
                deleteuntillpos = pos
                break
        # check deleteuntillpos; must be > 0
        if deleteuntillpos:
            deleted = self.codingblockgraphs[0:deleteuntillpos]
            ###################################################################
            if verbose:
                for cbg in deleted: print "DEL:",cbg
            ###################################################################
            self.codingblockgraphs = self.codingblockgraphs[deleteuntillpos:]
            self.codingblockgraphs[0].IS_FIRST = True
            self.codingblockgraphs[0]._CBGinterface5p = None
            # update this novel IS_FIRST CBG
            firstCBG = self.get_first_cbg()
            # make TranslationalStartSiteCollectionGraph
            firstCBG.harvest_elegiable_tss_sites()
            # do site alignment of TSS
            firstCBG._startcodongraph.collection2alignedsites()
            # set cbg attribute to all aligned sites -> required for is_optimal check
            for atssgra in firstCBG._startcodongraph.alignedsites:
                atssgra._codingblockgraph = firstCBG
            # return True -> we changed the GeneStructure
            return True
        else:
            # no deletion/change...
            return False

    # end of function first_codingblock_analyses_remove_falsecbgs


    def construct_first_tiny_cbg_by_most_likely_acceptor_phase(self,
        max_exon_nt_length              = 200,   # not to large, but quit large
        max_intron_nt_length            = 200,   # keep short; to long -> overpredictions
        allow_non_canonical_donor       = True,  # use it; construct_first_tiny_cbg likely
                                                 # failed because of this aspect 
        allow_non_canonical_acceptor    = False, # not implemented at all
        enlarge_acceptor_5p_boundary_by = 0,     # no enlargement needed, but do not shrink either
        enlarge_acceptor_3p_boundary_by = -5,    # descrease the area
        min_tss_pssm_score              = 5.0,   # keep it quit high
        min_donor_pssm_score            = 5.0,   # keep it quit high
        min_acceptor_pssm_score         = 3.5,   # keep quit high but not to high;
                                                 # acceptors in fungi not so strong signal
        verbose=False):
        """
        Make a tiny first CBG by ``shooting tiny TSSexons into the deep``
        """
        first = self.get_first_cbg()

        # only allow this function to run when IS_FIRST cbg is still unlikely
        if first.have_all_starts_upstream_of_omsr(): return False

        # do acceptor alignment of first CBG.
        first.harvest_elegiable_acceptor_sites(
                    enlarge_5p_boundary_by = enlarge_acceptor_5p_boundary_by,
                    enlarge_3p_boundary_by = enlarge_acceptor_3p_boundary_by,
                    min_acceptor_pssm_score= min_acceptor_pssm_score,
                    )
        first._spliceacceptorgraph.alignedsites = sort_by_cumulative_score(
            [ conversion.SpliceSiteCollectionGraph2AlignedSpliceSiteGraph(
              algsite,max_node_count=first.organism_set_size()
            ) for algsite in first._spliceacceptorgraph.find_conserved_sites()])
        # check if there is a strong acceptor signal

        ########################################################################
        if verbose:
            print "firstCBG:", first
            print "ALIGNED ACCEPTORS:"
            if first._spliceacceptorgraph:
                print first._spliceacceptorgraph
                if first._spliceacceptorgraph.alignedsites:
                    print first._spliceacceptorgraph.alignedsites[0]
        ########################################################################

        if not first._spliceacceptorgraph.alignedsites:
            # no proper aligned sites
            return False

        # get the best aligned acceptor site and do check 
        # on its object type and its pssm objects (acceptors)
        bestalgaccepgra = first._spliceacceptorgraph.alignedsites[0]

        if bestalgaccepgra.__class__.__name__ == 'AcceptorSiteCollectionGraph':
            # no poper aligned site
            return False
        if bestalgaccepgra.node_count() != first.organism_set_size():
            # not all species/genes covered
            return False

        if min([ accobj.pssm_score for accobj in bestalgaccepgra ._node_object.values() ]) <\
        min_acceptor_pssm_score:
            # to low acceptor site PSSM score
            for accobj in bestalgaccepgra ._node_object.values(): print accobj
            return False

        omsr = first.overall_minimal_spanning_range()
        ECG = ExonCollectionGraph()
        for organism in first.organism_set():
            node     = first.node_by_organism(organism)
            theorf   = first.get_orfs_of_graph(organism=organism)[0]
            acceptor = bestalgaccepgra.get_organism_objects(organism)[0]
            ####################################################################
            if verbose: print organism, acceptor
            ####################################################################

            # get elegiable orfs based on acceptor and intron length
            orflist = self.input[organism]['orfs'].get_elegiable_orfs(
                    max_orf_start= acceptor.pos,
                    min_orf_end= acceptor.pos - max_intron_nt_length,
                    has_starts=True )

            ####################################################################
            if verbose: print organism, "pot.orfs:", [orf.id for orf in orflist]
            ####################################################################

            # tmp store node and objects in a dict; we only allow a SINGLE,
            # highest scoring TSS exon of length <= 8. This because at this stage
            # quit easily a few of these tiny TSS options occur, which can - by chance -
            # result in wrong asignment of the actual TSS exon

            node2object = {}
            for orf in orflist:
                results = find_leading_exon_on_orf(
                        orf,theorf,
                        allow_non_canonical_donor=allow_non_canonical_donor,
                        subsequent_acceptor=acceptor,
                        max_leadingexon_nt_length=max_exon_nt_length,
                        max_leadingexon_intron_nt_length=max_intron_nt_length,
                        min_donor_pssm_score = min_donor_pssm_score,
                        min_tss_pssm_score   = min_tss_pssm_score, 
                        )

                # loop over the exons and store them; check for really short ones 
                # that only a single one is incorporated
                for exon,intron in results:
                    node = (organism,orf.id,exon.start,exon.end)
                    node2object[node] = exon

            # order the objects in node2object on total pssm score
            orderedobjects = ordering.order_list_by_attribute(node2object.values(),'pssm_score',reversed=True)

            # store the small TSS exon nodes to the ECG
            cnt = 0
            for exonobj in orderedobjects:
                cnt += 1
                # only allow exons with high enough PSSM scores
                if exonobj.donor.pssm_score < min_donor_pssm_score or\
                exonobj.acceptor.pssm_score < min_tss_pssm_score:
                    continue
                node = None
                for n,obj in node2object.iteritems():
                    if obj == exonobj:
                        node = n
                        break
                ECG.add_node_and_object(node,exonobj)
                ################################################################
                if verbose: print " ", cnt, node, "\t%s\t%1.2f\t%1.2f\t(%s)\t%s" % (
                        ("   %s" % exonobj.length)[-4:], exonobj.acceptor.pssm_score,
                        exonobj.donor.pssm_score, exonobj.donor.phase,
                        exonobj.proteinsequence() )
                ################################################################

        # only continue if all organisms are represented in the ECG
        if first.organism_set_size() > ECG.organism_set_size():
            ####################################################################
            if verbose: print "%s > %s, not all organisms covered -> return False" % (
                    first.organism_set_size(), ECG.organism_set_size() )
            ####################################################################
            return False

        # create edges in the ECG; create edges AMONG ALL NODES!
        # thant means, not only by exons that are +/- equally sized
        ECG.create_edges(max_length_nt_difference=max_exon_nt_length,max_length_ratio=0.0)
        ########################################################################
        if verbose: "edges created:", ECG.edge_count(), ECG.node_count()
        ########################################################################

        # find complete grapghs in the ECG
        first_exon_graphs = ECG.find_fully_connected_subgraphs()
        ########################################################################
        if verbose:
            print len(first_exon_graphs)
            for feg in first_exon_graphs[0:10]:
                print feg
        ########################################################################

        # convert to CodingBlockGraphs
        new_first_cbgs = []
        _cbg2feg = {} # dict that checks which cbg originates from which feg
        for feg in first_exon_graphs[0:3]:
            cbg = ExonCollectionGraph2CodingBlockGraph(feg,is_first=True,firstCBG=first)
            if cbg != None and cbg != False:
                new_first_cbgs.append( cbg )
                _cbg2feg[cbg] = feg
                if verbose: print cbg

        if not new_first_cbgs:
            ####################################################################
            if verbose: print "no first exon graph convertable to new CBG -> return False"
            ####################################################################
            return False

        # order by total weight; first ECG is the best ECG!
        new_first_cbgs = ordering.order_graphlist_by_total_weight(new_first_cbgs)
        newcbg         = new_first_cbgs[0]
        firstExonGraph = _cbg2feg[newcbg]

        ########################################################################
        if verbose:
            print "FEGs (%s) converted to CBGs (%s)" % ( 
                len(first_exon_graphs), len(new_first_cbgs) )
            for cbg in new_first_cbgs: print "cbg:", cbg
        ########################################################################

        if newcbg.organism_set_size() != first.organism_set_size():
            ####################################################################
            if verbose:
                print "not all orgs in feg CBG %s -> return False" % (
                    newcbg.organism_set_size() )
            ####################################################################
            return False

        ################################################################################
        if verbose:
            print "NOVEL feg CBG found!!!"
            print newcbg
            print newcbg._splicedonorgraph.alignedsites[0]
            print newcbg._startcodongraph.alignedsites[0]
            newcbg.printmultiplealignment()
            print "haveallstarts()", newcbg.have_all_starts(),
            print "upstreamomsrstarts()", newcbg.have_all_starts_upstream_of_omsr()
        ################################################################################

        # save current best AlignedSpliceSiteGraphs
        algdsgra = newcbg._splicedonorgraph.alignedsites[0]

        # make cbgIF by allowing non-canonical donors 
        cbgIF = CodingBlockGraphInterface(newcbg,first)
        cbgIF.harvest_acceptor_sites(
                enlarge_5p_boundary_by = +5,
                enlarge_3p_boundary_by = -5 )
        cbgIF.harvest_donor_sites(allow_non_canonical=True)
        # override _USE_ENTROPY attribure in cbgIF for is_optimal_site() function
        cbgIF._USE_ENTROPY = False

        # hard-set the AlignedDonorSiteGraph of the fegCBG graph
        cbgIF._splicedonorgraph.alignedsites.insert( 0, algdsgra )
        cbgIF._optimal_aligned_donor = algdsgra 
        # hard-set the AlignedAcceptorSiteGraph of the (current) first CBG
        cbgIF._spliceacceptorgraph.alignedsites.insert( 0, bestalgaccepgra )
        cbgIF._optimal_aligned_acceptor = bestalgaccepgra 

        ########################################################################
        if verbose:
            print cbgIF
            cbgIF.interfaceproperties()
        ########################################################################

        if cbgIF.optimalitycheck().count(True) >= 2:
            # yes, this will become the new first CBG!
            first.IS_FIRST = False
            first._CBGinterface5p = cbgIF
            newcbg.IS_FIRST = True
            newcbg._CBGinterface3p = cbgIF
            self.codingblockgraphs.insert(0,newcbg)
            return True 
        else:
            # nope, not a compatible interface...
            return False 

    # end of function construct_first_tiny_cbg_by_most_likely_acceptor_phase


    def DEPRECATED_first_codingblock_analyses(self,sprdif_min_node_count=2,sprdif_min_aa_length=CBG_MIN_AA_LENGTH):
        """
        Assign the optimal first CBG by removing non-favourable ones in front of a more favorable one

        """

        # Step 1. Split the second CBG on left spanningrange difference
        for pos in range(0,len(self)):
            if self.codingblockgraphs[pos].IS_FIRST:
                first = self.codingblockgraphs[pos]
                if not first.IS_SPLITTED and first.node_count() < len(self.input.keys()) and len(self)-1 > pos:
                    thepos = pos+1
                    second = self.codingblockgraphs[thepos]

                    # iteratively split on spanningrange difference
                    splits = second.iteratively_split_codingblock_on_spanningrange_difference(
                            side='left',
                            sprdif_min_aa_length=sprdif_min_aa_length,
                            sprdif_min_node_count=sprdif_min_node_count
                            )
                    for split in splits:
                        print "SPLITTED", split.node_count(), split

                    # set the last split back into the genestructure
                    self.codingblockgraphs[thepos] = splits[-1]
                    # and insert novel splitted CBGs if available
                    for splitpos in range(len(splits)-2,-1,-1):
                        newsplittedcbg = splits[splitpos]
                        self.codingblockgraphs.insert(thepos,newsplittedcbg)

                # break out after this analyses
                break

        # Step 2. Now assign the optimal CBG by confirming have_all_starts and
        for pos in range(0,len(self)):
            if self.codingblockgraphs[pos].IS_FIRST:
                first = self.codingblockgraphs[pos]
                compatible_first_cbg = True
                if not first.have_all_starts() and first.node_count() < len(self.input):
                    # no, impossible first CBG
                    compatible_first_cbg = False
                elif not first.have_all_starts_upstream_of_omsr() and first.node_count() < len(self.input):
                    # no, impossible first CBG
                    compatible_first_cbg = False
                else:
                    pass
                # check the compatibility of this CBG as first CBG
                if compatible_first_cbg: break

                # now check if it is possible to shift the IS_FIRST attribute
                compatible_next_is_first_cbg = True
                if pos == len(self)-1:
                    # this is the last CBG in the complete list!
                    compatible_next_is_first_cbg = False
                elif first.IS_LAST:
                    # CBG is assigned as IS_LAST -> no shift possible!
                    compatible_next_is_first_cbg = False
                else:
                    pass

                if compatible_next_is_first_cbg:
                    # okay, set the IS_FIRST attribute to the next CBG
                    first.IS_FIRST = False
                    first.IS_IGNORED = True
                    first.IS_3P_SPLITTED = False
                    if not first.IS_5P_SPLITTED:
                        first.IS_SPLITTED = False
                    next = self.codingblockgraphs[pos+1]
                    next.IS_FIRST = True
                    next.IS_5P_SPLITTED = False
                    if not next.IS_3P_SPLITTED:
                        next.IS_SPLITTED = False

    # end of function first_codingblock_analyses

    
    def fix_nonsense_first_cbg_with_tss5p_noncoding_alignment_part_in_second_cbg(self,
        ignore_cbg_max_aa_size=MAX_CBG_REMOVE_NONSENSE_FIRST_AA_LENGTH,
        ignore_cbg_min_identity_ratio=MIN_CBG_REMOVE_NONSENSE_FIRST_GTG_ID_RATIO,
        ignore_cbg_gtg_topo_dif=MAX_CBG_REMOVE_NONSENSE_FIRST_GTG_TOPO_DIF,
        ignore_cbg_gtg_abs_dif=MAX_CBG_REMOVE_NONSENSE_FIRST_GTG_ABS_DIF,
        ignore_cbg_gtg_tcode=MIN_CBG_REMOVE_NONSENSE_FIRST_TCODE_OMSR,
        verbose=False):
        """
        Remove a nonsense firs CBG in case the 2th CBG likely has a utr5p alignment
        """
        # do some function entrance checks. Not oke -> no action, return False
        if not self.codingblockgraphs[0].IS_FIRST == True:
            # first CBG in que must be IS_FIRST
            return False
        elif self.codingblockgraphs[0].have_all_starts_upstream_of_omsr():
            # first CBG must be unacceptable as firstcbg
            return False
        elif self.codingblockgraphs[0].IS_LAST:
            # there is only a single CBG in GSG; function not applicable
            return False
        elif self.codingblockgraphs[1].__class__.__name__ ==\
        'LowSimilarityRegionCodingBlockGraph':
            # 2th CBG must not be a lsrCBG
            return False
        elif self.codingblockgraphs[1]._CBGinterface5p and\
        self.codingblockgraphs[1]._CBGinterface5p.is_optimal():
            # CBG interface between 1th and 2th CBG must be bogus
            return False
        elif len(self) == 1:
            # there must be a 2th CBG in the GSG
            return False
        elif self.codingblockgraphs[1].have_all_starts_upstream_of_omsr():
            # second must be seemingly unacceptable as well
            return False
        else:
            # Prerequisites for this function seem okay.
            # Now, do a more elaborate check on the current first CBG
            # to see if it is a potentially non-coding bogus alignment.
            # Get the genetree for comparison
            thiscbg = self.codingblockgraphs[0]
            thegtg  = self.genetree()
            cbggtg  = thiscbg.genetree()
            checks  = []
            checks.append( thiscbg.omsr_tcode_score() < ignore_cbg_gtg_tcode )
            checks.append( thiscbg.omsrlength() <= ignore_cbg_max_aa_size )
            checks.append( thegtg.graphalignmentdifference( cbggtg ) > ignore_cbg_gtg_topo_dif )
            checks.append( thegtg.absolutegraphalignmentdifference( cbggtg ) > ignore_cbg_gtg_abs_dif ) 
            checks.append( ( cbggtg.identity() / thegtg.identity() ) < ignore_cbg_min_identity_ratio )
            if verbose: print "CHECKS:", checks
            # If all are False, no reason to believe this one is bogus
            # If all are True,  absolutely a bogus one.
            # Threshold for deletion when  at least 2 Trues.
            if checks.count(True) >= 2:
                pass
            else:
                # do not allow deletion of this CBG -> break
                return False
    
        # check if this second CBG has a 'hidden' TSSgra due to utr5p noncoding area
        if has_tss5p_noncoding_alignment_part(self.codingblockgraphs[1],verbose=verbose):
            # hard-place the Aligned TSS gra that separates the utr5p from the TSS
            _assign_tssgra_with_tss5p_noncoding_alignment_part(self.codingblockgraphs[1])
            # and fix the GSG by removing the current first CBG
            deleted = self.codingblockgraphs.pop(0)
            self.codingblockgraphs[0].IS_FIRST = True
            self.codingblockgraphs[0]._CBGinterface5p = None
            # ready! return status True
            return True
        else:
            return False
    
    # end of function fix_nonsense_first_cbg_with_tss5p_noncoding_alignment_part_in_second_cbg


    def optimize_tssgraph(self,minimal_pssm_ratio=0.70,verbose=False):
        """
        Optimize the alignedTSSgraph of the first CBG

        @type  minimal_pssm_ratio: float
        @param minimal_pssm_ratio: minimal ratio to replace a more upstream TSS.

        @type  verbose: Boolean
        @param verbose: print status/debugging messages to STDOUT

        @rtype:  NoneBoolean
        @return: None (not applicable), True (changed/optimized) or False (not changed)
        """
        # get the first CBG / labelled as IS_FIRST
        firstcbg = self.get_first_cbg()

        # if not _startcodongraph (no TranslationalStartSiteCollectionGraph),
        # then this function is not applicable -> return None
        if not firstcbg._startcodongraph: return None

        # if no _startcodongraph.alignedsites (no AlignedTranslationalStartSiteGraphs),
        # then this function is not applicable -> return None
        if not firstcbg._startcodongraph.alignedsites: return None

        # get the best/highest scoring aligned TSS graph
        bestalgtssgra = firstcbg._startcodongraph.alignedsites[0]

        # check if all organisms/genes are covered in this algTSSgraph
        if not bestalgtssgra.node_count() == firstcbg.node_count():
            return None
        if not bestalgtssgra.node_count() == self.EXACT_SG_NODE_COUNT:
            return None

        # okay, start optimizing!
        ####################################################################
        if verbose:
            print bestalgtssgra
            print bestalgtssgra.get_ordered_nodes()
            for organism in firstcbg._startcodongraph.organism_set():
                for tss in firstcbg._startcodongraph.get_organism_objects(organism=organism):
                    print organism, tss
            print ""
        ####################################################################
        
        improve_tss_node2newnode = {}
        alternative_nodes_found = False
        for organism in firstcbg._startcodongraph.organism_set():
            currenttss  = bestalgtssgra.get_organism_objects(organism=organism)[0]
            cbgnode     = firstcbg.node_by_organism(organism)
            currentnode = ( cbgnode[0], cbgnode[1], currenttss.pos/3, currenttss.pos)
            for tss in firstcbg._startcodongraph.get_organism_objects(organism=organism):
                if tss.pos == currenttss.pos: continue
                if tss.pos < currenttss.pos and\
                tss.pssm_score / currenttss.pssm_score >= minimal_pssm_ratio:
                    ########################################################
                    if verbose:
                        print "CUR:", organism, currenttss
                        print "NEW:", organism, tss
                    ########################################################
                    tssnode = ( cbgnode[0], cbgnode[1], tss.pos/3, tss.pos)
                    improve_tss_node2newnode[currentnode] = tssnode
                    alternative_nodes_found = True
                    break
            else:
                improve_tss_node2newnode[currentnode] = None

        # if no alternative nodes found -> return False, not optimizable
        if not alternative_nodes_found: return False
            
        # create list of ordeded nodes we need to search for
        search_nodes = []
        for keynode, valuenode in improve_tss_node2newnode.iteritems():
            # if valuenode != None -> a better TSSnode, if not, not alternative
            if valuenode: search_nodes.append( valuenode )
            else:         search_nodes.append( keynode )
        search_nodes.sort()

        for algtssgra in firstcbg._startcodongraph.alignedsites:
            if algtssgra.get_ordered_nodes() == search_nodes:
                # yes, this aligned TSS graph exists!
                # append it to the front of the list of
                # alignedsites and exit function!
                firstcbg._startcodongraph.alignedsites.insert(0, algtssgra)
                ############################################################
                if verbose:
                    print "# NEW algTSSgra found and set as optimal:"
                    print firstcbg._startcodongraph.alignedsites[0]
                ############################################################
                return True
        else:
            return False

    # end of function optimize_tssgraph


    def construct_upstream_exon_by_most_likely_acceptor_phase(self,
        max_exon_nt_length              = 75,    # not to large, but quit large
        min_exon_nt_length              = 10,    # small.
        max_intron_nt_length            = 110,   # keep short; to long -> overpredictions
        min_intron_nt_length            = 45,    # small MIN_INTRON_NT_LENGTH
        allow_non_canonical_donor       = True,  # use it; construct_first_tiny_cbg likely
                                                 # failed because of this aspect 
        allow_non_canonical_acceptor    = False, # not implemented at all
        enlarge_acceptor_5p_boundary_by = 0,     # no enlargement needed, but do not shrink either
        enlarge_acceptor_3p_boundary_by = -5,    # descrease the area
        min_donor_pssm_score            = 5.0,   # keep it quit high
        min_acceptor_pssm_score         = 3.5,   # keep quit high but not to high;
                                                 # acceptors in fungi not so strong signal
        verbose=False):
        """
        Make a tiny first CBG by ``shooting tiny TSSexons into the deep``
        """
        first = self.get_first_cbg()

        # only allow this function to run when IS_FIRST cbg is still unlikely
        if first.have_all_starts_upstream_of_omsr(): return False

        # do acceptor alignment of first CBG.
        first.harvest_elegiable_acceptor_sites(
                    enlarge_5p_boundary_by = enlarge_acceptor_5p_boundary_by,
                    enlarge_3p_boundary_by = enlarge_acceptor_3p_boundary_by,
                    min_acceptor_pssm_score= min_acceptor_pssm_score,
                    )
        first._spliceacceptorgraph.alignedsites = sort_by_cumulative_score(
            [ conversion.SpliceSiteCollectionGraph2AlignedSpliceSiteGraph(
              algsite,max_node_count=first.organism_set_size()
            ) for algsite in first._spliceacceptorgraph.find_conserved_sites()])
        # check if there is a strong acceptor signal

        ########################################################################
        if verbose:
            print "firstCBG:", first
            print "ALIGNED ACCEPTORS:"
            if first._spliceacceptorgraph:
                print first._spliceacceptorgraph
                if first._spliceacceptorgraph.alignedsites:
                    print first._spliceacceptorgraph.alignedsites[0]
        ########################################################################

        if not first._spliceacceptorgraph.alignedsites:
            # no proper aligned sites
            return False

        # get the best aligned acceptor site and do check 
        # on its object type and its pssm objects (acceptors)
        bestalgaccepgra = first._spliceacceptorgraph.alignedsites[0]

        if bestalgaccepgra.__class__.__name__ == 'AcceptorSiteCollectionGraph':
            # no poper aligned site
            return False
        if bestalgaccepgra.node_count() != first.organism_set_size():
            # not all species/genes covered
            return False
        if min([ accobj.pssm_score for accobj in bestalgaccepgra ._node_object.values() ]) <\
        min_acceptor_pssm_score:
            # to low acceptor site PSSM score
            for accobj in bestalgaccepgra ._node_object.values(): print accobj
            return False

        omsr = first.overall_minimal_spanning_range()
        ECG = ExonCollectionGraph()
        for organism in first.organism_set():
            node     = first.node_by_organism(organism)
            theorf   = first.get_orfs_of_graph(organism=organism)[0]
            acceptor = bestalgaccepgra.get_organism_objects(organism)[0]
            ####################################################################
            if verbose: print organism, acceptor
            ####################################################################

            # get elegiable orfs based on acceptor and intron length
            orflist = self.input[organism]['orfs'].get_elegiable_orfs(
                    max_orf_start= acceptor.pos,
                    min_orf_end= acceptor.pos - max_intron_nt_length )

            ####################################################################
            if verbose: print organism, "pot.orfs:", [orf.id for orf in orflist]
            ####################################################################

            # tmp store node and objects in a dict; we only allow a SINGLE,
            # highest scoring TSS exon of length <= 8. This because at this stage
            # quit easily a few of these tiny TSS options occur, which can - by chance -
            # result in wrong asignment of the actual TSS exon


            min_donor_pssm_score = 2.0
            max_donor_pssm_score = 2.0

            node2object = {}
            for orf in orflist:
                results = scan_orf_for_tiny_exon(orf,
                        max_intron_nt_length = max_intron_nt_length,
                        min_intron_nt_length = min_intron_nt_length,
                        min_tinyexon_nt_length   = min_exon_nt_length,
                        max_tinyexon_nt_length   = max_exon_nt_length,
                        allow_non_canonical_donor=allow_non_canonical_donor,
                        subsequent_acceptor_site=acceptor,
                        min_donor_pssm_score = min_donor_pssm_score,
                        min_acceptor_pssm_score   = min_acceptor_pssm_score, 
                        )

                # loop over the exons and store them; check for really short ones 
                # that only a single one is incorporated
                for exon in results:
                    node = (organism,orf.id,exon.start,exon.end)
                    node2object[node] = exon

            # order the objects in node2object on total pssm score
            orderedobjects = ordering.order_list_by_attribute(node2object.values(),'pssm_score',reversed=True)

            # store the small TSS exon nodes to the ECG
            cnt = 0
            for exonobj in orderedobjects:
                cnt += 1
                # only allow exons with high enough PSSM scores
                if exonobj.donor.pssm_score < min_donor_pssm_score or\
                exonobj.acceptor.pssm_score < min_acceptor_pssm_score:
                    continue
                node = None
                for n,obj in node2object.iteritems():
                    if obj == exonobj:
                        node = n
                        break
                ECG.add_node_and_object(node,exonobj)
                ################################################################
                if verbose: print " ", cnt, node, "\t%s\t%1.2f\t%1.2f\t(%s)\t%s" % (
                        ("   %s" % exonobj.length)[-4:], exonobj.acceptor.pssm_score,
                        exonobj.donor.pssm_score, exonobj.donor.phase,
                        exonobj.proteinsequence() )
                ################################################################

        ## only continue if all organisms are represented in the ECG
        #if first.organism_set_size() > ECG.organism_set_size():
        #    ####################################################################
        #    if verbose: print "%s > %s, not all organisms covered -> return False" % (
        #            first.organism_set_size(), ECG.organism_set_size() )
        #    ####################################################################
        #    return False

        # create edges in the ECG; create edges AMONG ALL NODES!
        # thant means, not only by exons that are +/- equally sized
        ECG.create_edges(max_length_nt_difference=15,max_length_ratio=0.6)
        ########################################################################
        if verbose: "edges created:", ECG.edge_count(), ECG.node_count()
        ########################################################################

        # find complete grapghs in the ECG
        first_exon_graphs = ECG.find_fully_connected_subgraphs()
        ########################################################################
        if verbose:
            print len(first_exon_graphs)
            for feg in first_exon_graphs[0:10]:
                print feg
        ########################################################################

        # convert to CodingBlockGraphs
        new_first_cbgs = []
        _cbg2feg = {} # dict that checks which cbg originates from which feg
        for feg in first_exon_graphs[0:3]:
            cbg = ExonCollectionGraph2CodingBlockGraph(feg,is_first=True,firstCBG=first)
            if cbg != None and cbg != False:
                new_first_cbgs.append( cbg )
                _cbg2feg[cbg] = feg
                if verbose: print cbg

        if not new_first_cbgs:
            ####################################################################
            if verbose: print "no exon graph convertable to new CBG -> return False"
            ####################################################################
            return False


        new_first_cbg = new_first_cbgs[0]
        new_first_cbg.IS_FIRST = True
        from graph_genestructure import mutualorgsofcbgs2partialGSG
        partGSG = mutualorgsofcbgs2partialGSG(self.input,cbgs=[ new_first_cbg ])
        partGSG.construct_first_tiny_cbg()
        print partGSG
        for cbg in partGSG:
            print cbg
            cbg.printmultiplealignment()
            print cbg._CBGinterface3p

    # end of function construct_upstream_exon_by_most_likely_acceptor_phase

# end of class FirstCodingBlockGraphFunctions


def has_tss5p_noncoding_alignment_part(cbg,max_5p_3p_ratio=0.45,verbose=False):
    """
    Does this CBG have a (long) 5putr upstream of the highest aligned TSS graph?

    @type  max_5p_3p_ratio: float
    @param max_5p_3p_ratio: ratio (0.0-1.0) 5p/3p from best TSS to create a CBG split

    @type  verbose: Boolean
    @param verbose: print status/debugging messages to STDOUT

    @rtype:  Boolean
    @return: True or False
    """

    # get the optimal K(s) Aligned TSS graph for this CBG
    bestTSSgra = _get_optimal_ks_tssgra(cbg)

    ####################################################################
    if verbose:
        print cbg
        cbg.printmultiplealignment()
        print bestTSSgra
    ####################################################################

    # no optimal aligned TSSgraph available...
    if not bestTSSgra: return False

    # translate TSS graph to AA postions per organism
    tssaacoords = {}
    for k, obj in bestTSSgra._node_object.iteritems():
        cbgnode = ( k[0], k[1] )
        aapos   = k[2]
        tssaacoords[cbgnode] = aapos

    # obtain bitscore normalized by aa length for the fraction
    # 5p and 3p of the highest scoring TSS per organism
    omsr = cbg.overall_minimal_spanning_range()
    ratios = []
    for (key,nodeQ,nodeS), pacbporf in cbg.pacbps.iteritems():
        orgQ = cbg.organism_by_node(nodeQ)
        orgS = cbg.organism_by_node(nodeS)
        # calculate 5p side bitscore part
        bitscore5porgQ = pacbporf.bitscore_slice_by_abs_protein_query(
            min(omsr[nodeQ]), tssaacoords[nodeQ]+1 )
        bitscore5porgS = pacbporf.bitscore_slice_by_abs_protein_sbjct(
            min(omsr[nodeS]), tssaacoords[nodeS]+1 )

        ratio5porgQ = float(bitscore5porgQ) / ( tssaacoords[nodeQ] - min(omsr[nodeQ]) )
        ratio5porgS = float(bitscore5porgS) / ( tssaacoords[nodeS] - min(omsr[nodeS]) )
        ratio5p = min([ ratio5porgQ, ratio5porgS ])
        # calculate 3p side bitscore part
        bitscore3porgQ = pacbporf.bitscore_slice_by_abs_protein_query(
            tssaacoords[nodeQ], max(omsr[nodeQ])+1 )
        bitscore3porgS = pacbporf.bitscore_slice_by_abs_protein_sbjct(
            tssaacoords[nodeS], max(omsr[nodeS])+1 )
        ratio3porgQ = float(bitscore3porgQ) / ( max(omsr[nodeQ]) - tssaacoords[nodeQ] )
        ratio3porgS = float(bitscore3porgS) / ( max(omsr[nodeS]) - tssaacoords[nodeS] )
        ratio3p = min([ ratio3porgQ, ratio3porgS ])
        # append to ratios
        ratios.append( ( ratio5p, ratio3p ) ) 
        ####################################################################
        if verbose:
            print pacbporf
            print nodeQ, nodeS, "\t", (bitscore5porgQ,bitscore5porgS),
            print (bitscore3porgQ,bitscore3porgS),
            print "%1.1f %1.1f" % (ratio5p,ratio3p), 
            #print min(omsr[nodeQ]), tssaacoords[nodeQ] , max(omsr[nodeQ])
            print "Q5p:", tssaacoords[nodeQ]+1 - min(omsr[nodeQ]), "S5p:", tssaacoords[nodeS]+1 - min(omsr[nodeS]),
            print "Q3p:", max(omsr[nodeQ])+1 - tssaacoords[nodeQ], "S3p:", max(omsr[nodeS])+1 - tssaacoords[nodeS]  
        ####################################################################


    if ratios:
        ratio5p = sum([ item[0] for item in ratios]) / len(ratios)
        ratio3p = sum([ item[1] for item in ratios]) / len(ratios)
        if ratio3p == 0.0:
            overall_5p_3p_ratio = 0.0 # ZeroDivision Error avoided!
        else:
            overall_5p_3p_ratio = ratio5p / ratio3p      
    else:
        # very unusual / impossible case of no gaterhed ratios
        ratio5p = 0.0
        ratio3p = 0.0
        overall_5p_3p_ratio = 0.0
    ####################################################################
    if verbose:
        print "5p of TSS:", ratio5p
        print "3p of TSS:", ratio3p
        print "5p/3p:", overall_5p_3p_ratio
    ####################################################################

    if overall_5p_3p_ratio == 0.0:
        return False
    elif overall_5p_3p_ratio > max_5p_3p_ratio:
        return False
    else:
        # ratio is low enough to represent a (putative) utr5p alignment part
        return True

# end of function has_tss5p_noncoding_alignment_part


def _get_optimal_ks_tssgra(cbg):
    """
    Get the highest scoring K(s) aligned TSS graph from this CBG

    @attention: USE WITH CARE! -> cbg._startcodongraph is OVERWRITTEN!
    @attention: use only from within has_tss5p_noncoding_alignment_part()
    @attention: use only from within _assign_tssgra_with_tss5p_noncoding_alignment_part()

    @type  cbg: CodingBlockGraph
    @param cbg: CodingBlockGraph instance (with cbgIF objects!)

    @rtype:  AlignedTranslationalStartSiteGraph or False
    @return: AlignedTranslationalStartSiteGraph or False
    """
    cbg.harvest_elegiable_tss_sites(skip_nonelegiable_sites=False)
    if not cbg._startcodongraph: return False
    cbg._startcodongraph.collection2alignedsites()
    if not cbg._startcodongraph.alignedsites: return False

    # get the highest scoring TSS graph that is present in ALL organisms
    # agree with the nodes of the cbg
    bestTSSgra = None
    for gra in cbg._startcodongraph.alignedsites:
        if gra.__class__.__name__ ==\
        'AlignedTranslationalStartSiteGraph' and\
        gra.node_count() == cbg.node_count():
            bestTSSgra = gra
            break
    if not bestTSSgra:
        return False
    else:
        return bestTSSgra

# end of function _get_optimal_ks_tssgra


def _assign_tssgra_with_tss5p_noncoding_alignment_part(cbg,verbose=False):
    """
    Helper function for fix_nonsense_first_cbg_with_tss5p_noncoding_alignment_part_in_second_cbg
    """
    # get the optimal K(s) Aligned TSS graph for this CBG
    bestTSSgra = _get_optimal_ks_tssgra(cbg)

    # and insert this graph in front of the list of alignedsites
    cbg._startcodongraph.alignedsites.insert(0,bestTSSgra)

# end of function _assign_tssgra_with_tss5p_noncoding_alignment_part


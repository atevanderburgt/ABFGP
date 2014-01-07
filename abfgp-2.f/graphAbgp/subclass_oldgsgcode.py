"""
Old code in GSG that is not wise to delete yet
"""

# Module metadata
__authors__ = "Ate van der Burgt"
__license__ = "MIT"


# Python Imports
from sets import Set
from copy import deepcopy

# Global variables
from settings.codingblockgraph import *
from settings.genestructure import *
from settings.inframeintron import *


class DeprecatedGSGFunctions:
    """
    """
    def constuct_leading_exon_codingblockgraph(self,organism=None):
        """
        """
        SHORT_LEADINGEXON_MAXIMAL_NT_TSS_OFFSET = 75

        print "#"*10, "constuct_leading_exon_codingblockgraph", organism
        scaffold_organism = 'ncu'
        scaffoldStartNode = ('ncu',500,1500)
        sg = self.codingblockgraphs[0]
        for organism in self.codingblockgraphs[0].organism_set():
            if organism == scaffold_organism: continue
            for ( (a,b,c,d),(g1,o1),(g2,o2) ), pacbporf in sg.pacbps.iteritems():
                if organism not in (g1,g2):   continue
                if scaffold_organism == g1:   pass
                elif scaffold_organism == g2: pacbporf = swap_query_and_sbjct(pacbporf)
                else:                         continue
                algpos = pacbporf.dnaposition_query( scaffoldStartNode[2], forced_return=True )
                algstartpos = pacbporf._positions[pacbporf._original_alignment_pos_start]
                print organism, algpos, pacbporf._positions[0], "ORIG:", pacbporf._positions[pacbporf._original_alignment_pos_start]
                print "ALGSTARTPOS: ORIG:", algstartpos
                print pacbporf
                # take the first (and only) orf of this organism
                orf_of_org = self.codingblockgraphs[0].get_orfs_of_graph(organism=organism)[0]
                print "SPLICE SITE RANGE:", sg._spliceacceptorgraph.get_consideredsplicesiterange(organism), len(orf_of_org._acceptor_sites), len(self.codingblockgraphs[0]._spliceacceptorgraph.get_organism_objects(organism))
                print orf_of_org, "ID:", orf_of_org.id
                #print [ a.pos for a in orf_of_org._acceptor_sites ]
                print [ a.pos for a in self.codingblockgraphs[0]._spliceacceptorgraph.get_organism_objects(organism) ]

                for orf in self.input[organism]['orfs'].orfs:
                    for acceptor in self.codingblockgraphs[0]._spliceacceptorgraph.get_organism_objects(organism):
                        # recalculate the length we are looking for!
                        lengthLookingFor = algstartpos.query_dna_start - scaffoldStartNode[2]
                        if acceptor.pos > algstartpos.sbjct_dna_start:
                            lengthLookingFor = lengthLookingFor + ( acceptor.pos - algstartpos.sbjct_dna_start )
                        else:
                            lengthLookingFor = lengthLookingFor - ( algstartpos.sbjct_dna_start - acceptor.pos )
                        res = find_leading_exon_on_orf(orf,orf_of_org,subsequent_acceptor=acceptor)
                        if res:
                            for exon,intron in res:
                                dist = lengthLookingFor-(exon.donor.pos-exon.acceptor.pos)
                                if abs(dist) < SHORT_LEADINGEXON_MAXIMAL_NT_TSS_OFFSET:
                                    print exon.acceptor.pos, exon.donor.pos, " .. ", intron.acceptor.pos, "[", exon.length, intron.length,
                                    print intron.acceptor.phase, "]", "(", orf.id, organism, ")",
                                    print exon.acceptor.pssm_score + exon.donor.pssm_score + intron.acceptor.pssm_score,
                                    print "\t", lengthLookingFor, dist
            print ""

    # end of function find_leading_exon


    def final_codingblock_analyses(self,sprdif_min_node_count=2,sprdif_min_aa_length=CBG_MIN_AA_LENGTH):
        """
        Check for which codingblock will be assigned as the STOP coding block
        """
        MIN_STOPCODON_WT_DISTANCE = 3.0     # TODO TODO: based on 5 organisms, but value depends on the number of organisms in the CBG (4-5-6-7)
        OPT_STOPCODON_WT_DISTANCE = 7.0     # TODO TODO: based on 5 organisms, but value depends on the number of organisms in the CBG (4-5-6-7)

        # position of where the hypothetical last CBG is now
        last = 0
        # make the aligned stop codon graphs
        for pos in range(0,len(self)):
            sg = self.codingblockgraphs[pos]
            if sg.IS_IGNORED: continue    # skip IGNORED codingblocks
            if sg.__class__.__name__ == 'LowSimilarityRegionCodingBlockGraph':
                continue

            # make AlignedStopCodonGraph
            print sg
            print sg.overall_minimal_spanning_range_sizes()
            sg.align_stop_codons()
            last = pos

        for pos in range(len(self)-1,0,-1):
            this = self.codingblockgraphs[pos]
            prev = self.codingblockgraphs[pos-1]

            # ignore a whole bunch of CBG combinations
            if this.IS_IGNORED: continue     # skip IGNORED codingblocks
            if prev.IS_IGNORED: continue     # skip IGNORED codingblocks
            if prev.IS_3P_SPLITTED: continue # skip 5P-3P SPLITTED codingblock combis
            if this.__class__.__name__ == 'LowSimilarityRegionCodingBlockGraph':
                continue
            if prev.__class__.__name__ == 'LowSimilarityRegionCodingBlockGraph':
                continue

            if this._stopcodongraph.total_weight() < MIN_STOPCODON_WT_DISTANCE and\
            prev._stopcodongraph.total_weight() >= OPT_STOPCODON_WT_DISTANCE and\
            len(Set(this.get_nodes()).intersection(prev.get_nodes())) == 0:
                # Discrepancy in alignments of orf lenghts, and not a single Orf
                # is shared by both CBGs -> intron in all species or intergenic region??
                # H0 -> intergenic region, prev.IS_LAST
                # H1 -> intron, strong intron signals (in all species) are needed!

                # some tmp printing
                distances = prev.distance_between_codingblocks( this ).values()
                mindist   = min(distances)
                maxdist   = max(distances)
                print "Discrepancy", pos-1, pos, prev._stopcodongraph.total_weight(), this._stopcodongraph.total_weight(), "AAdist:", mindist, maxdist
                print pos-1, pos
                print prev
                print this

                # check for aligned splicesites that contradict this
                # make SpliceSiteCollectionGraphs

                # no `projected_acceptors` and `forced_codingblock_ends` needed here!
                this.harvest_elegiable_acceptor_sites(projected_acceptors={},forced_codingblock_ends={},prev=prev)
                this._spliceacceptorgraph.collection2alignedsites(edges=self.EXACT_SG_NODE_COUNT-1,minimal_edges=2)

                # no `projected_donors` and `forced_codingblock_ends` needed here!
                prev.harvest_elegiable_donor_sites(projected_donors={},forced_codingblock_ends={},next=this)
                prev._splicedonorgraph.collection2alignedsites(edges=self.EXACT_SG_NODE_COUNT-1,minimal_edges=2)

                #if not this._spliceacceptorgraph:
                #    # make SpliceSiteCollectionGraphs
                #    # no `projected_acceptors` and `forced_codingblock_ends` needed here!
                #    this.harvest_elegiable_acceptor_sites(projected_acceptors={},forced_codingblock_ends={},prev=prev)
                #    this._spliceacceptorgraph.collection2alignedsites(edges=self.EXACT_SG_NODE_COUNT-1,minimal_edges=2)
                #if not prev._splicedonorgraph:
                #    # make SpliceSiteCollectionGraphs
                #    # no `projected_acceptors` and `forced_codingblock_ends` needed here!
                #    prev.harvest_elegiable_donor_sites(projected_donors={},forced_codingblock_ends={},next=this)
                #    prev._splicedonorgraph.collection2alignedsites(edges=self.EXACT_SG_NODE_COUNT-1,minimal_edges=2)

                if this._spliceacceptorgraph.alignedsites and prev._splicedonorgraph.alignedsites:
                    hypo_donor = prev._splicedonorgraph.alignedsites[0]
                    hypo_accep = this._spliceacceptorgraph.alignedsites[0]
                    if hypo_donor.organism_set_size() == prev.organism_set_size() and\
                    hypo_accep.organism_set_size() == this.organism_set_size() and\
                    hypo_donor.phase() == hypo_accep.phase():
                        # Unexpected, but H1 seems True. Although a dubious
                        # stop site total weight, there is a strong intron signal here
                        continue

                # no, not a strong signal; shift the STOP site to here!
                print "YAHOOOOO -> CHANGING LAST!"
                prev.IS_LAST = True
                this.IS_LAST = False
                this.IS_IGNORED = True
                for p in range(pos+1,len(self)):
                    self.codingblockgraphs[p].IS_IGNORED = True
                    self.codingblockgraphs[p].IS_LAST = False
                # and return the notification of the change
                return (pos, last)
        else:
            # if we reach this point, no changes done!
            return (last, last)

    # end of function final_codingblock_analyses


    def WORKING_construct_intermediate_tinyexon_with_known_splicesites(self,prev,next,
        min_exon_nt_length=5,
        max_exon_nt_length=21,
        max_intron_nt_length=SHORT_LEADINGEXON_MAX_INTRON_NT_LENGTH):
        """
        """
        print prev
        print prev._splicedonorgraph.alignedsites[0]
        print next
        print next._spliceacceptorgraph.alignedsites[0]
        ECG = ExonCollectionGraph()
        for organism in prev.organism_set():
            node = prev.node_by_organism(organism)
            donor_cnt = prev._splicedonorgraph.alignedsites[0].organism_set_size()
            accep_cnt = next._spliceacceptorgraph.alignedsites[0].organism_set_size()
            if organism in prev._splicedonorgraph.alignedsites[0].organism_set():
                preceding_donor_sites = [
                        prev._splicedonorgraph.alignedsites[0].get_organism_objects(organism)[0]
                        ]
                donor = preceding_donor_sites[0]
            else:
                preceding_donor_sites = prev._splicedonorgraph.get_organism_objects(organism)
                donor = None
            if organism in next._spliceacceptorgraph.alignedsites[0].organism_set():
                subsequent_acceptor_sites = [
                        next._spliceacceptorgraph.alignedsites[0].get_organism_objects(organism)[0]
                        ]
                acceptor = subsequent_acceptor_sites[0]
            else:
                subsequent_acceptor_sites = next._spliceacceptorgraph.get_organism_objects(organism)
                acceptor = None
            if not donor or not acceptor:
                preceding_donor_sites = prev._splicedonorgraph.get_organism_objects(organism)
                subsequent_acceptor_sites = next._spliceacceptorgraph.get_organism_objects(organism)
                tmp_donor = preceding_donor_sites[0]
                tmp_accep = subsequent_acceptor_sites[0]
                for acc in subsequent_acceptor_sites:
                    if acc.pos > tmp_accep.pos: tmp_accep = acc
                for don in preceding_donor_sites:
                    if don.pos < tmp_donor.pos: tmp_donor = don
                orflist  = self.input[organism]['orfs'].get_elegiable_orfs(
                    max_orf_start=tmp_accep.pos,
                    min_orf_end=tmp_donor.pos,
                    )
                print organism, len(orflist), [ orf.id for orf in orflist ]
                print tmp_donor
                print tmp_accep

            else:
                orflist  = self.input[organism]['orfs'].get_elegiable_orfs(
                    max_orf_start=acceptor.pos,
                    min_orf_end=donor.pos,
                    )
            donorOrf = prev.get_orfs_of_graph(organism=organism)[0]
            accepOrf = next.get_orfs_of_graph(organism=organism)[0]

            if not donor or not acceptor:
                results = []
                for orf in orflist:
                    exons = scan_orf_for_tiny_exon(orf,
                        min_acceptor_pos=tmp_donor.pos,
                        max_donor_pos=tmp_accep.pos,
                        min_tinyexon_nt_length=min_exon_nt_length,
                        max_tinyexon_nt_length=max_exon_nt_length,
                        )
                    results.extend(exons)

            else:
                results = bridge_two_pacbporfs_by_tinyexon(
                    donorOrf,accepOrf,
                    preceding_donor_sites=preceding_donor_sites,
                    subsequent_acceptor_sites=subsequent_acceptor_sites,
                    orflist=orflist,
                    min_tinyexon_nt_length=min_exon_nt_length,
                    max_tinyexon_nt_length=max_exon_nt_length,
                    max_tinyexon_intron_nt_length=max_intron_nt_length,
                    min_donor_pssm_score = -2.0,
                    min_acceptor_pssm_score = -2.0,
                    )
            for tinyexon in results:
                print organism, tinyexon.orf.id, tinyexon
                node = (organism,tinyexon.orf.id,tinyexon.start,tinyexon.end)
                ECG.add_node_and_object(node,tinyexon)

        # create edges in the ECG
        ECG.create_edges()
        # search for complete graphs in this
        print ECG.node_count(), ECG.edge_count()
        exon_graphs = ECG.find_fully_connected_subgraphs()

        print exon_graphs[0]

        ## convert to a CodingBlockGraph
        #newcbg = ExonCollectionGraph2CodingBlockGraph(firstExonGraph,firstCBG=first)
        #print firstExonGraph
        #print newcbg
        #print newcbg._splicedonorgraph.alignedsites[0]
        #print newcbg._startcodongraph.alignedsites[0]

    # end of function construct_intermediate_tinyexon_with_known_splicesites



# end of class DeprecatedGSGFunctions 

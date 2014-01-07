"""
Functions for blastp analyses of a CBG junction
"""

# Module metadata
__authors__ = "Ate van der Burgt"
__license__ = "MIT"

# graphAbgp Imports
from graphAbgp import PacbpCollectionGraph

# Abgp imports
from lib_fasta import parseFasta,writeMultiFasta
from lib_blastp import blastall_seq2db, formatdb
from lib_hmm import create_hmmdb_for_neighbouring_cbgs
from lib_stopwatch import StopWatch
from abgp_etc import _file_cleanup
import pacb
from listofcodingblockgraphs import ListOfCodingBlockGraphs

# Python Imports
from sets import Set

# Global variable Imports
CBG_JUNCTION_BLAST2PACBPCOL_EXTRA_BLASTP_PARAMS = {'F': 'F', 'e': '20','W': 2 } # Wordsize=2!!
CBG_JUNCTION_BLAST2PACBPCOL_OMSR_2_AA_MASK = 4

def blastanalysescbgjunction(gsg,prevCBG,nextCBG,
    omit_cbg_orfs = False,
    omit_non_cbg_orfs = False,
    extra_blastp_params=CBG_JUNCTION_BLAST2PACBPCOL_EXTRA_BLASTP_PARAMS,
    omsr_2_mask_aa_length_correction=CBG_JUNCTION_BLAST2PACBPCOL_OMSR_2_AA_MASK,
    verbose=False):
    """
    """
    ############################################################
    if verbose:
        stw = StopWatch('blastanalysescbgjunction')
        stw.start()
    ############################################################
    orfs = {}
    if not omit_cbg_orfs:
        # gather Orfs from prevCBG and nextCBG
        for org,orflist, in prevCBG.get_orfs_of_graph().iteritems():
            orf = orflist[0]
            orfs[(org,orf.id)] = orf
        for org,orflist, in nextCBG.get_orfs_of_graph().iteritems():
            orf = orflist[0]
            orfs[(org,orf.id)] = orf

    ############################################################
    if verbose:
        print stw.lap(), "orfs (1):",len(orfs)
        print _format_orf_nodes_to_string(orfs.keys())
    ############################################################

    # create masked fasta database in a dict
    fastadbmfa = parseFasta(
        create_hmmdb_for_neighbouring_cbgs(
            gsg.input,prevCBG,nextCBG,
            omsr_2_mask_aa_length_correction=omsr_2_mask_aa_length_correction,
            ).split("\n")
        )

    ############################################################
    if verbose: print stw.lap(), "fasta db (1):",len(fastadbmfa)
    ############################################################

    # remove ORFs that do not belong to prevCBG and nextCBG,
    # or that DO belong to prevCBG and nextCBG, or neither
    fastaheaders = fastadbmfa.keys()
    for header in fastaheaders:
        org,orfid = header.split("_orf_")
        orfid = int(orfid)
        node = (org,orfid)

        # check for the omit_non_cbg_orfs criterion
        add_orf = False
        if omit_non_cbg_orfs:
            if node not in orfs:
               del(fastadbmfa[header])
        else:
            add_orf = True

        # check for the omit_cbg_orfs criterion
        if omit_cbg_orfs and node in orfs:
            del(fastadbmfa[header])

        if add_orf:
            # get this Orf and add to orfs
            orfs[node] = gsg.input[org]['orfs'].get_orf_by_id(orfid)

    ############################################################
    if verbose:
        print stw.lap(), "fasta db (2):",len(fastadbmfa)
        print _format_fastadbmfa_nodes_to_string(fastadbmfa.keys())
    ############################################################

    ############################################################
    if verbose:
        print stw.lap(), "orfs (2):",len(orfs)
        print _format_orf_nodes_to_string(orfs.keys())
    ############################################################

    # no query/sbjct range left at all
    if not fastadbmfa: return [] 

    # check if all organisms are still covered
    orgSet = Set([ k.split("_orf_")[0] for k in fastadbmfa.keys()])
    if orgSet.symmetric_difference(gsg.organism_set()):
        return [] 

    # create !single! fasta database
    fastadbname = prevCBG.barcode()+"_"+nextCBG.barcode()+".mfa"
    writeMultiFasta(fastadbmfa,fastadbname)
    formatdb(fname=fastadbname)

    # remap the identifiers of the orf objects i.o.t....
    multifastas = {}
    blastdbs = {}
    pacbpcol    = PacbpCollectionGraph()
    dpcpacbpcol = PacbpCollectionGraph() # ``deepcopied`` variant for pacbps

    ############################################################
    if verbose: print stw.lap(), "blastp starting"
    ############################################################

    for orgQ,orgS in prevCBG.pairwisecrosscombinations_organism():

        for nodeQ,orfQ in orfs.iteritems():
            # only blast the (masked) Orfs of orgQ
            if prevCBG.organism_by_node(nodeQ) != orgQ: continue
            # get the masked protein sequence of this orfObj
            header = orgQ+"_orf_"+str(orfQ.id)
            # check if key exists in fastadbmfa. In a case where
            # an Orf is masked out completely, it is absent here!
            if not fastadbmfa.has_key(header): continue
            protseq = fastadbmfa[orgQ+"_orf_"+str(orfQ.id)]
            # run blast_seqs2db
            blastrec = blastall_seq2db(orfQ.id,protseq,fastadbname,
                    extra_blastp_params=extra_blastp_params)
            # omit empty blast records
            if len(blastrec.alignments) == 0: continue

            for alignment in blastrec.alignments:
                # get sbjct Org and Orf identifiers
                _orgS,_orfSid = alignment.title.replace(">","").split("_orf_")
                if _orgS != orgS: continue
                nodeS = (_orgS,int(_orfSid))
                orfS  = orfs[nodeS]
               
                # take only the *best* HSP (highest scoring first one)
                hsp = alignment.hsps[0]

                # correct to absolute positions
                hsp.query_start = hsp.query_start + orfQ.protein_startPY
                hsp.sbjct_start = hsp.sbjct_start + orfS.protein_startPY

                # initialize the PacbP
                pacbporf = pacb.conversion.pacbp2pacbporf(
                        pacb.PacbP(blastp_hsp=hsp),orfQ,orfS)

                ################################################################
                if verbose:
                    print pacbporf, orgQ,orgS, orfQ
                    print pacbporf.query
                    print pacbporf.match
                    print pacbporf.sbjct
                ################################################################

                # create nodes; ( Organism Identifier, Orf Identifier )
                nodeQ = ( orgQ, orfQ.id )
                nodeS = ( orgS, orfS.id )
                uqkey = pacbporf.construct_unique_key(nodeQ,nodeS)
                if not nodeQ in pacbpcol.get_nodes(): pacbpcol.add_node(nodeQ)
                if not nodeS in pacbpcol.get_nodes(): pacbpcol.add_node(nodeS)
                pacbpcol.add_edge(nodeQ,nodeS,wt=pacbporf.bitscore)
                # store to dpcpacbpcol -> pacbpcol is broken in pieces lateron!
                dpcpacbpcol.pacbps[(uqkey,nodeQ,nodeS)] = pacbporf

    ############################################################
    if verbose: print stw.lap(), "blastp done"
    ############################################################

    # file cleanup
    _file_cleanup(multifastas.values())
    _file_cleanup(["formatdb.log"])
    _file_cleanup([ fname+".*" for fname in blastdbs.values()])

    # check if all Organism/Gene identifiers are covered in PacbPs
    if not pacbpcol.organism_set_size() == gsg.organism_set_size():
        return [] 

    # ``deepcopy`` PacbPcollection pacbpcol to dpcpacbpcol
    # In dpcpacbpcol the actual PacbPORFs are stores & kept,
    # whereas pacbpcol itself is splitted in CBGs (which
    # function does not yet (!?) take the actual pacbps into account)
    dpcpacbpcol.add_nodes( pacbpcol.get_nodes() )
    for (uqkey,nodeQ,nodeS) in dpcpacbpcol.pacbps.keys():
        (bitscore,length,orfQid,orfSid) = uqkey
        dpcpacbpcol.add_edge(nodeQ,nodeS,wt=bitscore)

    ################################################################
    if verbose:
        print pacbpcol
        print "PCG bitscores:",
        print [ p.bitscore for p in dpcpacbpcol.pacbps.values() ]
        print "PCG nodes:", dpcpacbpcol.get_ordered_nodes()
    ################################################################

    #### do some transformations on the pacbpcol
    ####pacbpcol.remove_low_connectivity_nodes(min_connectivity=gsg.EXACT_SG_NODE_COUNT-1)
    ####splittedCBGs = pacbpcol.find_fully_connected_subgraphs(
    ####        edges=gsg.node_count()-1 , max_missing_edges=0 )
    ##### convert to list of CBGs and do some transformations
    ####cbgList = ListOfCodingBlockGraphs(splittedCBGs,input={},crossdata={})
    ####cbgList.remove_all_but_complete_cbgs()
    ####cbgList.remove_cbgs_with_lt_nodes(gsg.EXACT_SG_NODE_COUNT)
    ####cbgList.harvest_pacbps_from_pacbpcollection(dpcpacbpcol)
    ####cbgList.remove_cbgs_without_omsr()
    ####cbgList.update_edge_weights_by_minimal_spanning_range()
    ####cbgList.order_list_by_attribute(order_by='total_weight',reversed=True)

    min_connectivity = max([1,gsg.EXACT_SG_NODE_COUNT-1-2])
    pacbpcol.remove_low_connectivity_nodes(min_connectivity=min_connectivity)
    max_missing_edges = gsg.EXACT_SG_NODE_COUNT - 3
    splittedCBGs = pacbpcol.find_fully_connected_subgraphs(
            edges=gsg.node_count()-1 , max_missing_edges=max_missing_edges )
    # convert to list of CBGs and do some transformations
    cbgList = ListOfCodingBlockGraphs(splittedCBGs,input={},crossdata={})
    cbgList.remove_all_but_cbgs()
    cbgList.harvest_pacbps_from_pacbpcollection(dpcpacbpcol)
    cbgList.make_pacbps_for_missing_edges()
    cbgList.remove_all_but_complete_cbgs()
    cbgList.remove_cbgs_with_lt_nodes(gsg.EXACT_SG_NODE_COUNT)
    cbgList.remove_cbgs_without_omsr()
    cbgList.update_edge_weights_by_minimal_spanning_range()
    cbgList.order_list_by_attribute(order_by='total_weight',reversed=True)

    # and create_cache() for these CBGs
    for cbg in cbgList: cbg.create_cache()

    ####################################################################
    if verbose:
        print stw.lap(), "CBGs created", len(cbgList)
        for newcbg in cbgList: print "new:", newcbg
    ####################################################################

    # return list with CBGs
    return cbgList.codingblockgraphs

# end of function blastanalysescbgjunction


def _format_orf_nodes_to_string(nodes):
    """ orf node is of formant (org,orfid) """
    _data = []
    orgs = list(Set([node[0] for node in nodes]))
    orgs.sort()
    for org in orgs:
        orfids = []
        for node in nodes:
            if node[0] == org: orfids.append(node[1])
        orfids.sort()
        if orfids: _data.append( "%s: %s" % (
            org, " ".join([ str(id) for id in orfids ]) ) )
    return "\t"+" ".join(_data)

# end of function _format_orf_nodes_to_string


def _format_fastadbmfa_nodes_to_string(nodes):
    """ fastadbmfa node is of formant "<org>_orf_<orfid>" """
    orfnodes = []
    for node in nodes:
        org,orfid = node.split("_orf_")
        orfnodes.append( ( org, int(orfid) ) )
    return _format_orf_nodes_to_string(orfnodes)

# end of function _format_fastadbmfa_nodes_to_string


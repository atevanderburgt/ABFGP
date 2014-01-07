"""
Functions for obtaining frameshift (data) of a CodingBlockGraph
"""

# Module metadata
__authors__ = "Ate van der Burgt"
__license__ = "MIT"

# graphAbgp Imports
from graphAbgp import PacbpCollectionGraph
from graphAbgp.codingblock_ordering import relatively_positioned_towards

# Abgp imports
from listofcodingblockgraphs import ListOfCodingBlockGraphs
from lib_fasta import writeMultiFasta
from lib_blastp import blastall_seq2db, formatdb
from lib_orfset import GetorfOrfSet
from abgp_etc import _file_cleanup
import pacb

# Python Imports


def get_frameshifted_cbg(cbg,input,verbose=True):
    """
    Get a CBG with frameshifts (in some of if Orfs) compared to this CBG

    @type  cbg: CodingBlockGraph
    @param cbg: CodingBlockGraph to check for frameshifts

    @type  input: dict
    @param input: input <dict data structure> with lists of Orfs

    @type  verbose: Boolean
    @param verbose: print intermediate info to STDOUT for debugging purposes

    @rtype:  CodingBlockGraph or None
    @return: CodingBlockGraph (when existing) or None
    """

    # get elegiable lists of Orfs
    orfs = _get_elegiable_frameshift_orfsets(cbg,input)

    # check how many Orfs are elgiable...
    if sum([ len(l.orfs) for l in orfs.values() ]) == cbg.node_count():
        # no frameshift possible here...
        return None

    # remap the identifiers of the orf objects i.o.t....
    multifastas = {}
    blastdbs = {}
    pacbpcol    = PacbpCollectionGraph()
    dpcpacbpcol = PacbpCollectionGraph() # ``deepcopied`` variant for pacbps

    for org in orfs.keys():
        # REMAP fastaheaders as ids to retrieve the Orfs after blast..
        for orf in orfs[org].orfs: orf.fastaheader = str(orf.id)
        fname = "%s_frameshiftcbg_%s.mfa" % (org,cbg.barcode())
        writeMultiFasta(orfs[org].tofastadict(),fname)
        multifastas[org] = fname
        ########################################################################
        if verbose:
            print "ORFS:", org, len(orfs[org].orfs),
            print [ orf.id for orf in orfs[org].orfs ],
            print [ str(orf) for orf in orfs[org].orfs ]
        ########################################################################

    for orgQ,orgS in cbg.pairwisecrosscombinations_organism():
        # create blastdb if it does not exist yet
        if not blastdbs.has_key(orgS):
            formatdb(fname=multifastas[orgS])
            blastdbs[orgS] = multifastas[orgS]

        for orfQ in orfs[orgQ].orfs:
            # run blast_seqs2db
            blastrec = blastall_seq2db(orfQ.id,orfQ.protein_sequence,
                        dbname="./"+blastdbs[orgS])
            if len(blastrec.alignments) == 0: continue

            for alignment in blastrec.alignments:
                # obtain coordinates from sbjct orf identifier
                orfid = alignment.title.replace(">","").split(" ")[0].replace("_","")
                orfS = orfs[orgS].get_orf_by_id(int(orfid))

                nodeQ = ( orgQ, orfQ.id )
                nodeS = ( orgS, orfS.id )
                if nodeQ in cbg.get_nodes() and nodeS in cbg.get_nodes():
                    pacbporf = cbg.get_pacbps_by_nodes(node1=nodeQ,node2=nodeS)[0]

                else:
                    # take only the *best* HSP (highest scoring first one)
                    hsp = alignment.hsps[0]
    
                    # correct to absolute positions
                    hsp.query_start = hsp.query_start + orfQ.protein_startPY
                    hsp.sbjct_start = hsp.sbjct_start + orfS.protein_startPY
    
                    # initialize the PacbP
                    pacbporf = pacb.conversion.pacbp2pacbporf(
                            pacb.PacbP(blastp_hsp=hsp),orfQ,orfS)
                    ############################################################
                    if verbose: print "NEW:", pacbporf
                    ############################################################

                uqkey = pacbporf.construct_unique_key(nodeQ,nodeS)
                if not nodeQ in pacbpcol.get_nodes(): pacbpcol.add_node(nodeQ)
                if not nodeS in pacbpcol.get_nodes(): pacbpcol.add_node(nodeS)
                pacbpcol.add_edge(nodeQ,nodeS,wt=pacbporf.bitscore)
                # store to dpcpacbpcol -> pacbpcol is broken in pieces lateron!
                dpcpacbpcol.pacbps[(uqkey,nodeQ,nodeS)] = pacbporf

    # file cleanup
    _file_cleanup(multifastas.values())
    _file_cleanup(["formatdb.log"])
    _file_cleanup([ fname+".*" for fname in blastdbs.values()])

    if not pacbpcol.organism_set_size() == cbg.organism_set_size():
        ############################################################
        if verbose: print "org_set_size() PCG < CBG" 
        ############################################################
        # no CBG on the reverse strand
        return None

    # ``deepcopy`` PacbPcollection
    dpcpacbpcol.add_nodes( pacbpcol.get_nodes() )
    for (uqkey,nodeQ,nodeS) in dpcpacbpcol.pacbps.keys():
        (bitscore,length,orfQid,orfSid) = uqkey
        dpcpacbpcol.add_edge(nodeQ,nodeS,wt=bitscore)

    ############################################################################
    if verbose:
        print pacbpcol, "bitscores:",
        print [ pacbporf.bitscore for pacbporf in dpcpacbpcol.pacbps.values() ]
    ############################################################################

    # do some transformations on the pacbpcol
    pacbpcol.remove_low_connectivity_nodes(min_connectivity=cbg.node_count()-1)
    splittedCBGs = pacbpcol.find_fully_connected_subgraphs(
            edges=cbg.node_count()-1 , max_missing_edges=0 )
    # convert to list of CBGs and do some transformations
    cbgList = ListOfCodingBlockGraphs(splittedCBGs,input={},crossdata={})
    cbgList.remove_all_but_cbgs()
    cbgList.remove_cbgs_with_lt_nodes(cbg.node_count())
    cbgList.harvest_pacbps_from_pacbpcollection(dpcpacbpcol)
    cbgList.remove_cbgs_without_omsr()
    cbgList.update_edge_weights_by_minimal_spanning_range()
    cbgList.order_graphlist_by_total_weight_and_identity()

    ############################################################################
    if verbose:
        print "FScbgs (%s)" % len(cbgList)
        for fscbg in cbgList: print fscbg 
    ############################################################################

    if not cbgList:
        # no (better) frameshifted CBG
        return None
    elif cbgList and not cbgList[0].node_set().symmetric_difference(cbg.node_set()):
        # best CBG is not frameshifted, but CBG itself
        return None
    else:
        # score the difference between the frameshifted and current CBG
        score_cbg   = cbg.total_weight() * cbg.omsr_identityscore()
        score_fscbg = cbgList[0].total_weight() * cbgList[0].omsr_identityscore()
        # check overlap between the frameshifted and current CBG 
        a,b,c,d,e,f,g = relatively_positioned_towards(cbgList[0],cbg)

        ########################################################################
        if verbose:
            print "CBG", cbg
            cbg.printmultiplealignment()
            for fscbg in cbgList:
                print "fsCBG:", fscbg
                fscbg.printmultiplealignment()
        ########################################################################

        if (c,d) == ( (0, 0, 1), (1, 0, 0) ) or (c,d) == ( (0, 0, 1), (1, 0, 0) ): 
            # CBG and frameshifted CBG do not share a single AA overlap...
            # This does not represent a frameshifted CBG as we searched for
            return False
        elif score_fscbg > score_cbg: 
            # return the highest scoring, frameshifted CBG
            return cbgList[0]
        else:
            # no, still not convinced that this is a frameshifted CBG
            return False

# end of function get_frameshifted_cbg


def _get_elegiable_frameshift_orfsets(cbg,input,extra_nt_offset=3):
    """
    Get elegiable lists of Orfs that might represent a frameshifted CBG

    @type  cbg: CodingBlockGraph
    @param cbg: CodingBlockGraph to reversecomplement

    @type  input: dict
    @param input: input <dict data structure> with lists of Orfs

    @rtype:  dictionary
    @return: CBG's Organism Identifiers as keys, GetorfOrfSet object as values
    """
    # initialize empty OrfLists
    retdict = {}
    # get current Orfs of graph
    orfs = cbg.get_orfs_of_graph()
    for node,omsr in cbg.overall_minimal_spanning_range().iteritems():
        org = cbg.organism_by_node(node)
        retdict[org] = GetorfOrfSet()
        thisorf = orfs[org][0]
        #orfsubset = input[org]['orfs'].get_elegiable_orfs(
        #        max_orf_start = min(omsr)*3+thisorf.frame+extra_nt_offset,
        #        min_orf_end= (max(omsr)+1)*3+thisorf.frame-extra_nt_offset
        #        )
        orfsubset = input[org]['orfs'].get_elegiable_orfs(
                max_orf_start = (max(omsr)+1)*3,
                min_orf_end= min(omsr)*3
                )
        retdict[org].orfs.extend( orfsubset )
    # return dict with orfList objects
    return retdict

# end of function _get_elegiable_frameshift_orfsets

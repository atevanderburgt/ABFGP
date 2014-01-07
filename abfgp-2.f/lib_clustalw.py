from parsers.clustalw import *
from copy import deepcopy
from lib_proteinsimilaritymatrix import make_clustalw_alignment_match


def clustalwinput2cbg(seqs,orfs,coords,nodes,
    matrix = None,
    minimal_overall_spanning_range_size = 3,
    verbose=False):
    """

    @type  seqs: dict
    @param seqs: dict with ORGANISM IDENTIFIER as keys, sequences as values

    @type  orfs: dict
    @param orfs: dict with ORGANISM IDENTIFIER as keys, Orf objects as values

    @type  coords: dict
    @param coords: dict with ORGANISM IDENTIFIER as keys, [ sta, end ] as values

    @type  nodes: list
    @param nodes: list with nodes corresponding to the ORGANISM IDENTIFIER in the dictionaries

    @attention: coordinates in coords should correspond to the sequneces in seqs!

    """
    # do clustalw and strip_alignment_for_exterior_gaps
    (algseqs,algm) = clustalw(seqs=seqs)
    ####################################################################
    if verbose: print seqs, "\n", algseqs, "\n", algm, "\n", coords
    ####################################################################
    _testalgseqs,_testalgm,_testcoords = strip_alignment_for_exterior_gaps(
        deepcopy(algseqs),deepcopy(algm),deepcopy(coords))
    if not _testalgm:
        ####################################################################
        if verbose: print "NO ALGM\n", seqs, "\n", _testalgseqs, "\n", _testalgm
        ####################################################################
        # alignment completely vanished by `strip_alignment_for_exterior_gaps`
        return None

    # do required import here (prevent circular imports)
    from graphAbgp.graph_codingblock import CodingBlockGraph
    from graphAbgp.exceptions import NoOverallMinimalSpanningRange
    from pacb import conversion as pacbconversion

    if not matrix:
        raise "No ProteinSimilarityMatrix applied!"

    # translate the clustalw alignment into an artificial CBG
    newcbg = CodingBlockGraph()
    newcbg.add_nodes(nodes)
    pacbp_is_none = False
    for nodeA,nodeB in newcbg.pairwisecrosscombinations_node():
        orgA = newcbg.organism_by_node(nodeA)
        orgB = newcbg.organism_by_node(nodeB)

        # create stripped alignments for this pair of sequences
        # do not forget to make deepcopies of the data structures!
        subcoords  = { orgA: coords[orgA], orgB: coords[orgB] }
        subalgseqs = { orgA: algseqs[orgA], orgB: algseqs[orgB] }
        _algseqs,_algm,_coords = strip_alignment_for_exterior_gaps(
            deepcopy(subalgseqs),deepcopy(algm),deepcopy(subcoords) )

        # recreate a pairwise ClustalW alignment string
        _algm = make_clustalw_alignment_match(
                _algseqs[orgA],_algseqs[orgB],
                matrix = matrix.matrix )

        # _algseqs keys are organisms, not nodes!
        alignment  = ( _algseqs[orgA], _algm, _algseqs[orgB] )
        paircoords = ( _coords[orgA][0], _coords[orgA][1],
                       _coords[orgB][0], _coords[orgB][1] )
        pacbp = pacbconversion.pacbp_from_clustalw(
                alignment=alignment,coords=paircoords)
        if pacbp == None:
            # pacbp is not creatable -> break i.o.t. return None
            pacbp_is_none = True
            break
        pacbporf = pacbconversion.pacbp2pacbporf(pacbp,orfs[orgA],orfs[orgB])
        ####################################################################
        if verbose:
            print orgA, orgB, pacbporf
            for item in alignment: print item
            print paircoords
        ####################################################################
        wt = pacbporf.bitscore
        pacbpkey = pacbporf.construct_unique_key(nodeA,nodeB)
        newcbg.add_edge(nodeA,nodeB,wt=wt)
        newcbg.pacbps[(pacbpkey,nodeA,nodeB)] = pacbporf

    # check if all pacbporfs are created succesfully
    if pacbp_is_none: return None

    # update edge weight by OMSR and return
    newcbg.MINIMAL_OVERAL_SPANNING_RANGE_SIZE =\
        minimal_overall_spanning_range_size

    if newcbg.has_overall_minimal_spanning_range():
        newcbg.update_edge_weights_by_minimal_spanning_range()
        try:
            newcbg.correct_pacbpgaps_nearby_omsr()
            return newcbg
        except NoOverallMinimalSpanningRange:
            return None
    else:
        return None

# end of function clustalwinput2cbg


def sprdif2clustalw2cbg(cbg,sprdif,SCAFFOLD_GAP_OMSR_OFFSET=1,verbose=False):
    """ """
    # gather sequence concerning the scaffold gap of the mutual nodes
    seqs, orfs, coords = {}, {}, {}
    nodes = sprdif.keys()
    for node in sprdif.keys():
        org = cbg.organism_by_node(node)
        sta = min( sprdif[node] ) - SCAFFOLD_GAP_OMSR_OFFSET
        end = max( sprdif[node] ) + SCAFFOLD_GAP_OMSR_OFFSET
        orf = cbg.get_orfs_of_graph(organism=org)[0]
        # correct a priori for out-of-range exceptions
        # due to SCAFFOLD_GAP_OMSR_OFFSET
        sta = max([ sta, orf.protein_startPY ])
        end = min([ end, orf.protein_endPY ])
        seq = orf.getaas(abs_pos_start=sta,abs_pos_end=end)
        seqs[org]   = seq
        orfs[org]   = orf
        coords[org] = [sta,end]

    # get (a) protein similarity matrix of the CBG
    protsimmtrx = cbg.pacbps.values()[0].MATRIX

    # call the clustalwinput2cbg function
    return clustalwinput2cbg(seqs,orfs,coords,nodes,matrix=protsimmtrx)

# end of function sprdif2clustalw2cbg


def WORKING_sprdif2clustalw2cbg(cbg,sprdif,SCAFFOLD_GAP_OMSR_OFFSET=1,verbose=False):
    """ """
    # gather sequence concerning the scaffold gap of the mutual nodes
    seqs, orfs, coords = {}, {}, {}
    for node in sprdif.keys():
        org = cbg.organism_by_node(node)
        sta = min( sprdif[node] ) - SCAFFOLD_GAP_OMSR_OFFSET
        end = max( sprdif[node] ) + SCAFFOLD_GAP_OMSR_OFFSET
        orf = cbg.get_orfs_of_graph(organism=org)[0]
        # correct a priori for out-of-range exceptions
        # due to SCAFFOLD_GAP_OMSR_OFFSET
        sta = max([ sta, orf.protein_startPY ])
        end = min([ end, orf.protein_endPY ])
        seq = orf.getaas(abs_pos_start=sta,abs_pos_end=end)
        seqs[org]   = seq
        orfs[org]   = orf
        coords[org] = [sta,end]

    # do clustalw and strip_alignment_for_exterior_gaps
    (algseqs,algm) = clustalw(seqs=seqs)
    ####################################################################
    if verbose: print seqs, "\n", algseqs, "\n", algm, "\n", coords
    ####################################################################
    _testalgseqs,_testalgm,_testcoords = strip_alignment_for_exterior_gaps(
        deepcopy(algseqs),deepcopy(algm),deepcopy(coords))
    if not _testalgm:
        ####################################################################
        if verbose: print "NO ALGM\n", seqs, "\n", _testalgseqs, "\n", _testalgm
        ####################################################################
        # alignment completely vanished by `strip_alignment_for_exterior_gaps`
        return None

    # do required import here (prevent circular imports)
    from graphAbgp.graph_codingblock import CodingBlockGraph
    from graphAbgp.exceptions import NoOverallMinimalSpanningRange
    from pacb import conversion as pacbconversion
    from lib_cexpander import cexpander_checkCBG4omsrbordergaps, ZeroUniformlyAlignedPositions

    # translate the clustalw alignment into an artificial CBG
    newcbg = CodingBlockGraph()
    newcbg.add_nodes(sprdif.keys())
    pacbp_is_none = False
    for nodeA,nodeB in newcbg.pairwisecrosscombinations_node():
        orgA = cbg.organism_by_node(nodeA)
        orgB = cbg.organism_by_node(nodeB)

        # create stripped alignments for this pair of sequences
        # do not forget to make deepcopies of the data structures!
        subcoords  = { orgA: coords[orgA], orgB: coords[orgB] }
        subalgseqs = { orgA: algseqs[orgA], orgB: algseqs[orgB] }
        _algseqs,_algm,_coords = strip_alignment_for_exterior_gaps(
            deepcopy(subalgseqs),deepcopy(algm),deepcopy(subcoords) )

        # get a/the ProteinSimilarityMatrix from the original PacbP(ORF)
        # and then recreate a pairwise ClustalW alignment string
        protsimmtrx = cbg.get_pacbps_by_nodes(node1=nodeA,node2=nodeB)[0].MATRIX
        _algm = make_clustalw_alignment_match(
                _algseqs[orgA],_algseqs[orgB],
                matrix = protsimmtrx.matrix )

        # _algseqs keys are organisms, not nodes!
        alignment  = ( _algseqs[orgA], _algm, _algseqs[orgB] )
        paircoords = ( _coords[orgA][0], _coords[orgA][1],
                       _coords[orgB][0], _coords[orgB][1] )
        pacbp = pacbconversion.pacbp_from_clustalw(
                alignment=alignment,coords=paircoords)
        if pacbp == None:
            # pacbp is not creatable -> break i.o.t. return None
            pacbp_is_none = True
            break
        pacbporf = pacbconversion.pacbp2pacbporf(pacbp,orfs[orgA],orfs[orgB])
        ####################################################################
        if verbose:
            print orgA, orgB, pacbporf
            for item in alignment: print item
            print paircoords
        ####################################################################
        wt = pacbporf.bitscore
        pacbpkey = pacbporf.construct_unique_key(nodeA,nodeB)
        newcbg.add_edge(nodeA,nodeB,wt=wt)
        newcbg.pacbps[(pacbpkey,nodeA,nodeB)] = pacbporf

    # check if all pacbporfs are created succesfully
    if pacbp_is_none: return None

    # update edge weight by OMSR and return
    newcbg.MINIMAL_OVERAL_SPANNING_RANGE_SIZE = 3
    if newcbg.has_overall_minimal_spanning_range():
        newcbg.update_edge_weights_by_minimal_spanning_range()
        try:
            newcbg.correct_pacbpgaps_nearby_omsr()
            return newcbg
        except NoOverallMinimalSpanningRange:
            return None
    else:
        return None

# end of function sprdif2clustalw2cbg



def WORKING_sprdif2clustalw2cbg(cbg,sprdif,SCAFFOLD_GAP_OMSR_OFFSET=0,verbose=False):
    """ """
    # gather sequence concerning the scaffold gap of the mutual nodes
    seqs, orfs, coords = {}, {}, {}
    for node in sprdif.keys():
        org = cbg.organism_by_node(node)
        sta = min( sprdif[node] ) - SCAFFOLD_GAP_OMSR_OFFSET
        end = max( sprdif[node] ) + SCAFFOLD_GAP_OMSR_OFFSET
        orf = cbg.get_orfs_of_graph(organism=org)[0]
        seq = orf.getaas(abs_pos_start=sta,abs_pos_end=end)
        seqs[org]   = seq
        orfs[org]   = orf
        coords[org] = [sta,end]

    # do clustalw and strip_alignment_for_exterior_gaps
    (_algseqs,_algm) = clustalw(seqs=seqs)
    ####################################################################
    if verbose: print seqs, "\n", _algseqs, "\n", _algm
    ####################################################################
    _algseqs,_algm,coords = strip_alignment_for_exterior_gaps(_algseqs,_algm,coords)
    if not _algm:
        ####################################################################
        if verbose: print "NO ALGM.??\n", seqs, "\n", _algseqs, "\n", _algm
        ####################################################################
        # alignment completely vanished by `strip_alignment_for_exterior_gaps`
        return None

    # do required import here (prevent circular imports)
    from graphAbgp.graph_codingblock import CodingBlockGraph
    from graphAbgp.exceptions import NoOverallMinimalSpanningRange
    from pacb import conversion as pacbconversion
    from lib_cexpander import cexpander_checkCBG4omsrbordergaps, ZeroUniformlyAlignedPositions

    # translate the clustalw alignment into an artificial CBG
    newcbg = CodingBlockGraph()
    newcbg.add_nodes(sprdif.keys())
    pacbp_is_none = False
    for nodeA,nodeB in newcbg.pairwisecrosscombinations_node():
        orgA       = cbg.organism_by_node(nodeA)
        orgB       = cbg.organism_by_node(nodeB)
        # _algseqs keys are organisms, not nodes!
        alignment  = ( _algseqs[orgA], _algm, _algseqs[orgB] )
        paircoords = ( coords[orgA][0], coords[org][1], coords[orgB][0], coords[orgB][1] )
        pacbp = pacbconversion.pacbp_from_clustalw(alignment=alignment,coords=paircoords)
        if pacbp == None:
            # pacbp is not creatable -> break i.o.t. return None
            pacbp_is_none = True
            break
        pacbporf = pacbconversion.pacbp2pacbporf(pacbp,orfs[orgA],orfs[orgB])
        wt = pacbporf.bitscore
        pacbpkey = pacbporf.construct_unique_key(nodeA,nodeB)
        newcbg.add_edge(nodeA,nodeB,wt=wt)
        newcbg.pacbps[(pacbpkey,nodeA,nodeB)] = pacbporf

    # check if all pacbporfs are created succesfully
    if pacbp_is_none: return None

    # update edge weight by OMSR and return
    newcbg.MINIMAL_OVERAL_SPANNING_RANGE_SIZE = 3
    if newcbg.has_overall_minimal_spanning_range():
        newcbg.update_edge_weights_by_minimal_spanning_range()
        try:
            newcbg.correct_pacbpgaps_nearby_omsr()
            return newcbg
        except NoOverallMinimalSpanningRange:
            return None
        #try:
        #    status = cexpander_checkCBG4omsrbordergaps(newcbg)
        #    return newcbg 
        #except NoOverallMinimalSpanningRange:
        #    return None
        #except ZeroUniformlyAlignedPositions:
        #    return None
        #except:
        #    return None
    else:
        return None

# end of function sprdif2clustalw2cbg

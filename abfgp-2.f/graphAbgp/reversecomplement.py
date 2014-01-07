"""
Functions for obtaining the reversecomplement (data) of a CodingBlockGraph
"""

# Module metadata
__authors__ = "Ate van der Burgt"
__license__ = "MIT"

# graphAbgp Imports
from graphAbgp import PacbpCollectionGraph
from graphAbgp.basic import BasicCodingBlockGraph

# Abgp imports
from lib_fasta import _reversecomplement, writeMultiFasta
from dna2prot import dna2protein,dna2proteinbyframe
from lib_proteinsimilaritymatrix import alignment2bitscore
from lib_blastp import blastall_seq2db, formatdb
from gene.orf import TcodeCodingOrf
from lib_orfset import OrfSet
from abgp_etc import _file_cleanup
import pacb
from listofcodingblockgraphs import ListOfCodingBlockGraphs

# Python Imports
from sets import Set

# Global Variable Import
from settings.executables import EXECUTABLE_GETORF_MINSIZE
STOP_CODONS     = ['TGA','TAA','TAG']
REV_STOP_CODONS = ['TCA','TTA','CTA']
CBG_REVERSECOMPLEMENT_MIN_REMOVAL_RATIO   = 1.15
CBG_REVERSECOMPLEMENT_MIN_DOUBTFULL_RATIO = 1.01
CBG_REVERSECOMPLEMENT_MIN_OVERLAP_RATIO   = 0.80 # or higher?

from settings.gff import (
    GFF_REVCBG_SIMILARITY_FSOURCE,
    GFF_REVCBG_SIMILARITY_FMETHOD, 
    GFF_REVCBG_SIMILARITY_GCLASS,
    )

class ReversecomlementCodingBlockGraph(BasicCodingBlockGraph):
    """
    """
    def __init__(self,**kwargs):
        """
        Initialize a ReversecomlementCodingBlockGraph
        """
        # Initialize from PacbpCollectionGraph
        BasicCodingBlockGraph.__init__(self,**kwargs)

        # short name tag for printing
        self._short_name = "revCBG"

        # defaults for togff() function
        self._gff = {
            'fref':     None,
            'fstrand':  '-',
            'fsource':  GFF_REVCBG_SIMILARITY_FSOURCE,
            'fmethod':  GFF_REVCBG_SIMILARITY_FMETHOD,
            'gclass':   GFF_REVCBG_SIMILARITY_GCLASS,
            'fphase':   '.',
            }

    # end of function __init__


    def __str__(self):
        """ """
        organdorfs = self.get_ordered_nodes()
        identity = " " # single space!
        if self.has_overall_minimal_spanning_range():
            allmsr = self.reversecomplement_overall_minimal_spanning_range()
            themsr = [ "OMSR" ]
            for node in self.get_ordered_nodes():
                themsr.append( "%s(%s):%s..%s" % (
                        self.organism_by_node(node),
                        node[1],
                        min(allmsr[node]),
                        max(allmsr[node]) ) )
            msrlen   = min(self.overall_minimal_spanning_range_sizes().values())
            identity = " id=%1.2f " % self.omsr_identityscore()
        else:
            msrlen = 0
            themsr = ["NO OMSR"]
            identity = " " # single space!
            for node in self.get_ordered_nodes():
                themsr.append( "%s(%s)" % (
                        self.organism_by_node(node),
                        node[1] ) )

        # obtain pacbp overshoot string
        pacbp_overshoot = " "
        if len(self.pacbps) - self.edge_count() :
            pacbp_overshoot = " pacbp:+%s " % (len(self.pacbps)-self.edge_count())

        # return the string representation of this revCBG
        return "<%s wt=%s%slen=%s cs=%1.1f %s[%s]>" % (
                self._short_name,
                self.total_weight(),
                identity,
                msrlen,
                self.connectivitysaturation(),
                pacbp_overshoot,
                " ".join(themsr)
                )
    
    # end of function __str__


    def reversecomplement_overall_minimal_spanning_range(self,organism=None,node=None):
        """
        Get the reversecomplement OMSR, so in forward (+) strand AA coords.

        @type  organism: *
        @param organism: organism identifier (or None)

        @type  node: *
        @param node: Node identifier (or None)

        @rtype:  dictionary (or Set if organism is specified)
        @return: dictionary with nodes (keys) and spanning range Sets
                 (values), or only the spanning range Set if an organism
                 or node identifier was specified
        """
        # first obtain regular OMSR data
        omsrdata = self.overall_minimal_spanning_range(
                    organism=organism,node=node)
        # now convert this to reversecomplement data based on Orf start-end
        if node or organism:
            if not organism: organism = self.organism_by_node(node)
            orf = self.get_orfs_of_graph(organism=organism)[0]
            minomsr = orf.protein_endPY-max(omsrdata) + orf.protein_startPY
            maxomsr = orf.protein_endPY-min(omsrdata) + orf.protein_startPY
            return Set(range(minomsr,maxomsr+1))
        else:
            for node in omsrdata.keys():
                organism = self.organism_by_node(node)
                orf = self.get_orfs_of_graph(organism=organism)[0]
                minomsr = orf.protein_endPY-max(omsrdata[node]) + orf.protein_startPY
                maxomsr = orf.protein_endPY-min(omsrdata[node]) + orf.protein_startPY
                omsrdata[node] = Set(range(minomsr,maxomsr+1))
            # return updated omsrdata range
            return omsrdata

    # end of function reversecomplement_overall_minimal_spanning_range


    def togff(self,organism,gff={}):
        """ """
        if organism not in self.organism_set():
            raise OrganismNotPresentInGraph

        # update gff data as stored in object
        self._gff.update(gff)

        # get orf of graph
        theorf = self.get_orfs_of_graph(organism=organism)[0]
        gffdata = []

        for (k,(org1,orf1),(org2,orf2)), pacbporf in self.pacbps.iteritems():
            if organism == org1:
                pacbporf.unextend_pacbporf()
                start = pacbporf._positions[0].query_dna_start+1
                end   = pacbporf._positions[-1].query_dna_end
                group       = "%s-%s-orfid-%s" % (gff['fref'],org2,orf2)
            elif organism == org2:
                pacbporf.unextend_pacbporf()
                start = pacbporf._positions[0].sbjct_dna_start+1
                end   = pacbporf._positions[-1].sbjct_dna_end
                group       = "%s-%s-orfid-%s" % (gff['fref'],org1,orf1)
            else:
                continue
    
            # initialize column9data dict
            column9data = {
                            'Note'      : "Similarity on the reverse strand",
                            'Length'    : pacbporf.length,
                            'Bitscore'  : pacbporf.bitscore,
                            'Similarity': pacbporf.similarity,
                            'Identity'  : pacbporf.identity,
                            }
            # column9data dict -> column9data string
            tmp = []
            for key,value in column9data.iteritems():
                if str(value).find(" ") > -1: value = "'%s'" % value
                tmp.append( "; %s %s" % (key,value) )
            column9string = "".join(tmp)

            # append GFF tuple
            gffdata.append( (
                self._gff['fref'],
                self._gff['fsource'],
                self._gff['fmethod'],
                start,
                end,
                "%1.3f" % pacbporf.identityscore,
                self._gff['fstrand'],
                self._gff['fphase'],
                "%s %s%s" % ( self._gff['gclass'], group, column9string )
                ) )

        # return gffdata list
        return gffdata

    # end of function togff

# end of class ReversecomlementCodingBlockGraph


def CodingBlockGraph2ReversecomlementCodingBlockGraph(cbg):
    """ """
    revcbg = ReversecomlementCodingBlockGraph()
    revcbg.nodes   = cbg.nodes
    revcbg.weights = cbg.weights
    revcbg.pacbps  = cbg.pacbps
    return revcbg

# end of function CodingBlockGraph2ReversecomlementCodingBlockGraph


def omsr_overlap_ratio(cbg,revcbg):
    """
    """
    # first check: are all Organism Identifiers from revCBG covered by CBG?
    # This can happen when K(s-x) revCBGs are taken into account
    if cbg.organism_set().difference(revcbg.organism_set()):
        # return 0.0 overlap ratio
        return 0.0

    # check if CBG and revCBG are overlapping
    omsrCBG = cbg.overall_minimal_spanning_range()
    omsrREV = revcbg.reversecomplement_overall_minimal_spanning_range()

    # create data on OMSR overlaps
    omsrsummary = []
    for nodeCBG, omsrdataCBG in omsrCBG.iteritems():
        org         = cbg.organism_by_node(nodeCBG)
        nodeREV     = revcbg.node_by_organism(org)
        omsrdataREV = omsrREV[nodeREV]
        # compare omsrdataCBG & omsrdataREV
        omsrsummary.append( [
            len(omsrdataCBG.difference(omsrdataREV)),
            len(omsrdataCBG.intersection(omsrdataREV)) ] )
    # create OMSR overlap summary
    omsr_intersection = sum([ i for (d,i) in omsrsummary ] )
    omsr_difference   = sum([ d for (d,i) in omsrsummary ] )
    omsr_ratio = omsr_intersection / float(omsr_difference + omsr_intersection)

    # return this ratio
    return omsr_ratio

# end of function omsr_overlap_ratio


def is_reversecomplement(cbg,revcbg=None,verbose=False):
    """
    Is this CBG more likely a ReversecomplementCodingBlockGraph ?

    @type  cbg: CodingBlockGraph
    @param cbg: CodingBlockGraph to reversecomplement

    @type  verbose: Boolean
    @param verbose: print intermediate info to STDOUT for debugging purposes

    @rtype:  NoneBoolean
    @return: False, None or True
    """
    if not revcbg: revcbg = get_reversecomplement(cbg)
    if not revcbg: return False

    # obtain omsr_overlap_ratio between CBG and revCBG
    omsr_ratio = omsr_overlap_ratio(cbg,revcbg)

    # check if the OMSR overlap is large enough to represent a revCBG
    if omsr_ratio < CBG_REVERSECOMPLEMENT_MIN_OVERLAP_RATIO:
        return False

    # compare total_weight of both CBG and revCBG
    tw_revcbg = float(revcbg.total_weight())
    tw_cbg    = float(cbg.total_weight())

    if tw_revcbg / tw_cbg < CBG_REVERSECOMPLEMENT_MIN_DOUBTFULL_RATIO:
        return False
    elif tw_revcbg / tw_cbg > CBG_REVERSECOMPLEMENT_MIN_REMOVAL_RATIO:
        return True
    else:
        return None

# end of function is_reversecomplement


def get_reversecomplement(cbg,verbose=False):
    """
    Get *the* ReversecomplementCodingBlockGraph of this CBG

    @type  cbg: CodingBlockGraph
    @param cbg: CodingBlockGraph to reversecomplement

    @type  verbose: Boolean
    @param verbose: print intermediate info to STDOUT for debugging purposes

    @rtype:  ReversecomplementCodingBlockGraph or None
    @return: ReversecomplementCodingBlockGraph (when existing) or None
    """
    # get ReversecomplementCBGs in frame 0 and 2
    revcbg_frame0 = get_reverse_cbg(cbg,0,verbose=verbose)
    revcbg_frame2 = get_reverse_cbg(cbg,2,verbose=verbose)
    # return the ReversecomplementCBG with highest total_weight
    if not revcbg_frame0 and not revcbg_frame2:
        return None
    elif revcbg_frame0 and revcbg_frame2:
        if revcbg_frame0.total_weight() > revcbg_frame2.total_weight():
            return revcbg_frame0
        else:
            return revcbg_frame2
    elif revcbg_frame0:
        return revcbg_frame0
    else:
        return revcbg_frame2

# end of function get_reversecomplement


def get_reverse_cbg(cbg,frame,verbose=False):
    """
    Get the ReversecomplementCodingBlockGraph in requested frame of this CBG

    @type  cbg: CodingBlockGraph
    @param cbg: CodingBlockGraph to reversecomplement

    @type  frame: integer
    @param frame: 0,1 or 2

    @type  verbose: Boolean
    @param verbose: print intermediate info to STDOUT for debugging purposes

    @rtype:  ReversecomplementCodingBlockGraph or None
    @return: ReversecomplementCodingBlockGraph (when existing) or None
    """
    min_orf_length = (cbg.omsrlength()/2)*3
    orfs = get_reverse_strand_orfsets(cbg,frame,min_orf_length=min_orf_length)

    # remap the identifiers of the orf objects i.o.t....
    multifastas = {}
    blastdbs = {}
    pacbpcol    = PacbpCollectionGraph()
    dpcpacbpcol = PacbpCollectionGraph() # ``deepcopied`` variant for pacbps

    for org in orfs.keys():
        fname = "%s_reversecbg_%s.mfa" % (org,cbg.barcode())
        writeMultiFasta(orfs[org].tofastadict(),fname)
        multifastas[org] = fname
        ########################################################################
        if verbose:
            print "ORFS:", org, len(orfs[org].orfs),
            print [len(o.protein_sequence) for o in orfs[org].orfs ]
        ########################################################################

    revpacbps = {}
    for orgQ,orgS in cbg.pairwisecrosscombinations_organism():
        # create blastdb if it does not exist yet
        if not blastdbs.has_key(orgS):
            formatdb(fname=multifastas[orgS])
            blastdbs[orgS] = multifastas[orgS]

        revpacbporfs = {}
        for orfQ in orfs[orgQ].orfs:
            # run blast_seqs2db
            blastrec = blastall_seq2db(orfQ.id,orfQ.protein_sequence,
                        dbname="./"+blastdbs[orgS])
            if len(blastrec.alignments) == 0: continue

            for alignment in blastrec.alignments:
                # obtain coordinates from sbjct orf identifier
                orfS = orfs[orgS].get_orf_by_id(alignment.title.replace(">",""))
                # take only the *best* HSP (highest scoring first one)
                hsp = alignment.hsps[0]
                # skip if hsp is very short
                if len(hsp.query) < cbg.omsrlength()/2: continue

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
                    ###pacbporf.print_protein_and_dna()
                ################################################################

                nodeQ = ( orgQ, orfQ.protein_startPY )
                nodeS = ( orgS, orfS.protein_startPY )
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
    cbgList.remove_all_but_complete_cbgs()
    cbgList.harvest_pacbps_from_pacbpcollection(dpcpacbpcol)
    cbgList.remove_cbgs_without_omsr()
    cbgList.update_edge_weights_by_minimal_spanning_range()
    cbgList.order_list_by_attribute(order_by='total_weight',reversed=True)

    ############################################################################
    if verbose:
        for revcbg in cbgList:
            print "revCBG:", revcbg
    ############################################################################

    if not cbgList:
        # no CBG on the reverse strand
        return None
    else:
        # return the highest scoring CBG as a ReversecomlementCodingBlockGraph
        return CodingBlockGraph2ReversecomlementCodingBlockGraph(
                cbgList.codingblockgraphs[0])

# end of function get_reverse_cbg


def get_reverse_strand_orfsets(cbg,frame,
    min_orf_length=EXECUTABLE_GETORF_MINSIZE):
    """
    TODO....

    @type  cbg: CodingBlockGraph
    @param cbg: CodingBlockGraph to reversecomplement

    @type  frame: integer
    @param frame: 0,1 or 2

    @type  min_orf_length: integer
    @param min_orf_length: minimal Orf length (in nt's, stop-to-stop)

    @rtype:  dictionary
    @return: CBG's Organism Identifiers as keys, OrfSet object as values
    """
    # initialize empty OrfLists
    retdict = {}
    for org in cbg.organism_set(): retdict[org] = OrfSet()

    for (pacbpkey,nodeQ,nodeS),pacbporf in cbg.pacbps.iteritems():
        orgQ = cbg.organism_by_node(nodeQ)
        orgS = cbg.organism_by_node(nodeS)
        orfQlist, orfSlist = get_overlapping_orfs_on_negative_strand(
                pacbporf,frame,min_orf_length=min_orf_length)

        # increment list of orfs from the Query organism
        if not retdict[orgQ].orfs:
            retdict[orgQ].orfs = orfQlist.orfs
        else:
            for orfObj in orfQlist.orfs:
                if not retdict[orgQ].get_orf_by_id(orfObj.id):
                    retdict[orgQ].orfs.append(orfObj)

        # increment list of orfs from the Sbjct organism
        if not retdict[orgS].orfs:
            retdict[orgS].orfs = orfSlist.orfs
        else:
            for orfObj in orfSlist.orfs:
                if not retdict[orgS].get_orf_by_id(orfObj.id):
                    retdict[orgS].orfs.append(orfObj)

    # return dict with orfList objects
    return retdict

# end of function get_reverse_strand_orfsets


def get_overlapping_orfs_on_negative_strand(pacbporf,frame,
    min_orf_length=EXECUTABLE_GETORF_MINSIZE):
    """
    TODO ...

    @type  pacbporf: PacbPORF
    @param pacbporf: PacbPORF object

    @type  frame: integer
    @param frame: 0,1 or 2

    @type  min_orf_length: integer
    @param min_orf_length: minimal Orf length (in nt's, stop-to-stop)

    @rtype:  tuple()
    @return: ( OrfSet, OrfSet ) from Query and Sbjct of PacbPORF
    """
    # initialise coordinates and sequences (for current frame)
    dnaQsta = frame + pacbporf.query_dna_start
    dnaQend = frame + pacbporf.query_dna_end
    dnaSsta = frame + pacbporf.sbjct_dna_start
    dnaSend = frame + pacbporf.sbjct_dna_end
    dnaQseq = pacbporf.orfQ.inputgenomicsequence
    dnaSseq = pacbporf.orfS.inputgenomicsequence

    # get sequences & coordinates
    dnaQsta,dnaQend,dnaQrev,protQrev =\
        _get_reverse_strand_coords_and_seqs(dnaQseq,dnaQsta,dnaQend)
    dnaSsta,dnaSend,dnaSrev,protSrev =\
        _get_reverse_strand_coords_and_seqs(dnaSseq,dnaSsta,dnaSend)

    # get orfs on this reverse strand
    revQorfs = _split_sequences_in_orfs(dnaQsta,dnaQend,dnaQrev,protQrev,
                                        min_orf_length=min_orf_length)
    revSorfs = _split_sequences_in_orfs(dnaSsta,dnaSend,dnaSrev,protSrev,
                                        min_orf_length=min_orf_length)
    # return orfslists
    return (revQorfs,revSorfs)

# end of function get_overlapping_orfs_on_negative_strand


def _split_sequences_in_orfs(start,end,dnaseq,protseq,
    min_orf_length=EXECUTABLE_GETORF_MINSIZE):
    """
    TODO ...

    @type  start: integer
    @param start:

    @type  end: integer
    @param end:

    @type  dnaseq: string
    @param dnaseq: DNA sequence string

    @type  protseq: string
    @param protseq: Protein sequence string

    @attention: protseq *must* be the exact translation from dnaseq

    @type  min_orf_length: integer
    @param min_orf_length: minimal Orf length (in nt's, stop-to-stop)

    @rtype:  OrfSet 
    @return: OrfSet object containing list of Orfs
    """
    orfListObj = OrfSet(sequence=dnaseq)
    orfListObj.strand = '-'
    offset,prevoffset = 0, 0
    protseqlen = len(protseq)
    while True:
        offset = protseq.find("*",offset)
        if offset == -1 and prevoffset == protseqlen:
            # EOF protein sequence reached
            break
        elif offset == -1:
            # novel ORF at end of protseq found (EOF dnaseq, no *)
            orfDNAsta = 0 
            orfDNAend = end - (prevoffset*3)
            if orfDNAend-orfDNAsta >= min_orf_length:
                orfObj = _create_new_rev_orf(orfDNAsta,orfDNAend,
                                offset,prevoffset,dnaseq,protseq)
                # append orfObj to orfListObj
                orfListObj.orfs.append(orfObj)
            # EOF protein sequence reached
            break
        elif offset == 0:
            pass
        else:
            # novel ORF found!
            orfDNAsta = end - (offset*3)
            orfDNAend = end - (prevoffset*3)
            if orfDNAend-orfDNAsta >= min_orf_length:
                orfObj = _create_new_rev_orf(orfDNAsta,orfDNAend,
                                offset,prevoffset,dnaseq,protseq)
                # append orfObj to orfListObj
                orfListObj.orfs.append(orfObj)

        # set prevoffset, increase offset
        prevoffset=offset+1
        offset+=1

    # return list with (reverse) Orfs
    return orfListObj

# end of function _split_sequences_in_orfs


def _create_new_rev_orf(orfDNAsta,orfDNAend,offset,prevoffset,dnaseq,protseq):
    """
    @attention: Helper function of _split_sequences_in_orfs
    """
    id = "%s_%s" % (orfDNAsta,orfDNAend)
    # make NEW dymmy genome sequence for each RevOrf ...
    dummyinputgenomicsequence = "n"*(orfDNAsta-(prevoffset*3))+dnaseq
    return TcodeCodingOrf(
            id,protseq[prevoffset:offset],
            inputgenomicsequence=dummyinputgenomicsequence,
            force_id=id,
            start=orfDNAsta+1, # python->GGB/GFF coord
            end=orfDNAend )

# end of function _create_new_rev_orf


def _mirrir_coords(s,e,mains,maine):
    """ """
    ms = maine-e
    me = maine-s
    return ms, me

# end of function _mirrir_coords

def _get_reverse_strand_coords_and_seqs(dnaseq,start,end):
    """
    TODO ....

    @type  dnaseq: string
    @param dnaseq: DNA sequence string

    @type  start: integer
    @param start:

    @type  end: integer
    @param end:

    @rtype:  tuple
    @return: (start,end,dnarevseq,protseq)
    """
    # obtain position of the exterior stop codons on the other strand
    # take care with EOF sequence; after while loop +3; so stop at -6 end
    if end >= len(dnaseq)-6:
        # end offset already in the danger zone of reaching EOF sequence
        # shift untill in reach (3nt offset) of EOF sequence
        while end >= len(dnaseq)-3: end-=3
    else:
        while end < len(dnaseq)-6:
            if dnaseq[end:end+3].upper() in REV_STOP_CODONS: break
            end += 3
    # take care with EOF sequence; after while loop -3; so stop at +6 end
    if start < 6:
        # start offset already in the danger zone of reaching EOF sequence
        # shift untill in reach (3nt offset) of EOF sequence
        while start < 3: start += 3
    else:
        while start >= 6:
            if dnaseq[start-3:start].upper() in REV_STOP_CODONS: break
            start -= 3
    # get DNA & Protein sequences on the reverse strand
    start-=3
    end+=3
    dnarevseq = _reversecomplement(dnaseq[start:end])
    protseq   = dna2protein(dnarevseq)
    # return coords & sequences
    return (start,end,dnarevseq,protseq)

# end of function _get_reverse_strand_coords_and_seqs




"""

for cbg in GSG:

for cbg in GSG:
    bitsum = sum([ pacbporf.bits for pacbporf in cbg.pacbps.values() ])
    print cbg
    orfdict = get_reverse_strand_orfsets(cbg,0)
    

from lib_fasta import _reversecomplement
from dna2prot import dna2protein,dna2proteinbyframe
from lib_proteinsimilaritymatrix import alignment2bitscore
from lib_blastp import blastall_seq2seq
from graphAbgp.reversecomplement import *
for cbg in GSG:
    bitsum = sum([ pacbporf.bits for pacbporf in cbg.pacbps.values() ])
    framebitsums = [ 0, 0, 0 ]
    for pacbporf in cbg.pacbps.values():
        STOP_CODONS = ['TGA','TAA','TAG']
        REV_STOP_CODONS = ['TCA','TTA','CTA']
        for frame in [0,1,2]:
            dnaQsta = frame + pacbporf.query_dna_start
            dnaQend = frame + pacbporf.query_dna_end
            dnaSsta = frame + pacbporf.sbjct_dna_start
            dnaSend = frame + pacbporf.sbjct_dna_end
            print frame, dnaQsta, dnaQend, dnaSsta, dnaSend
            while True:
                if pacbporf.orfQ.inputgenomicsequence[dnaQend:dnaQend+3].upper() in REV_STOP_CODONS:
                    break
                dnaQend += 3
            while True:
                if pacbporf.orfQ.inputgenomicsequence[dnaQsta-3:dnaQsta].upper() in REV_STOP_CODONS:
                    break
                dnaQsta -= 3
            while True:
                if pacbporf.orfS.inputgenomicsequence[dnaSend:dnaSend+3].upper() in REV_STOP_CODONS:
                    break
                dnaSend += 3
            while True:
                if pacbporf.orfS.inputgenomicsequence[dnaSsta-3:dnaSsta].upper() in REV_STOP_CODONS:
                    break
                dnaSsta -= 3

            dnaQrev = _reversecomplement(pacbporf.orfQ.inputgenomicsequence[dnaQsta-3:dnaQend+3])
            dnaSrev = _reversecomplement(pacbporf.orfS.inputgenomicsequence[dnaSsta-3:dnaSend+3])
            protQ = dna2protein(dnaQrev)
            protS = dna2protein(dnaSrev)
            print frame, dnaQsta, dnaQend, dnaSsta, dnaSend,
            print _reversecomplement(pacbporf.orfQ.inputgenomicsequence[dnaQsta-3:dnaQsta]),
            print _reversecomplement(pacbporf.orfQ.inputgenomicsequence[dnaQend:dnaQend+3])
            print protQ
            print protS


            revorfdata = []
            protseq = protQ
            dnasta  = dnaQsta-3
            dnaend  = dnaQend+3
            offset,prevoffset = 0, 0
            while True:
                offset = protQ.find("*",offset)
                if offset == -1: break
                if offset == 0:
                    pass
                else:
                    orfDNAsta = dnasta+(prevoffset*3)
                    orfDNAend = dnasta+(offset*3)
                    if orfDNAend-orfDNAsta > (len(pacbporf)/2)*3:
                        print prevoffset, offset, protseq[prevoffset:offset], orfDNAsta, orfDNAend
                # set prevoffset, increase offset
                prevoffset=offset+1
                offset+=1


        q,s = pacbporf.get_unextended_aligned_dna_sequences()
        revq = _reversecomplement(q)
        revs = _reversecomplement(s)
        revQaa0 = dna2proteinbyframe(revq,0).replace("X","-")
        revQaa1 = dna2proteinbyframe(revq,1).replace("X","-")
        revQaa2 = dna2proteinbyframe(revq,2).replace("X","-")
        revSaa0 = dna2proteinbyframe(revs,0).replace("X","-")
        revSaa1 = dna2proteinbyframe(revs,1).replace("X","-")
        revSaa2 = dna2proteinbyframe(revs,2).replace("X","-")
        framebitsums[0] += alignment2bitscore(revQaa0,"",revSaa0,matrix=DEFAULT_MATRIX.matrix)
        framebitsums[1] += alignment2bitscore(revQaa1,"",revSaa1,matrix=DEFAULT_MATRIX.matrix)
        framebitsums[2] += alignment2bitscore(revQaa2,"",revSaa2,matrix=DEFAULT_MATRIX.matrix)
        print revQaa0
        print revSaa0
        blastoutput = blastall_seq2seq(fastadata=(('Q',revQaa0.replace("-","")),('S',revSaa0.replace("-",""))))
        print dir(blastoutput)
        try:
            besthsp = blastoutput.alignments[0].hsps[0]
            print besthsp.query, besthsp.bits
            print besthsp.match
            print besthsp.sbjct
        except IndexError:
            pass
    print cbg.total_weight(),"\t", bitsum, "\t", bitsum > max(framebitsums), "\t%1.2f\t" % (float(max(framebitsums))/float(bitsum)), framebitsums


"""

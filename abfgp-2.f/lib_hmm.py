"""

"""

# Module metadata
__authors__ = "Ate van der Burgt"
__license__ = "MIT"

# import hmmparser items
from parsers.hmm import hmmbuild_protein, hmmsearch_protein

# Imports
import os
from sets import Set
from copy import deepcopy
from re import ( 
    compile as re_compile,
    match as re_match,
    )

# import ABGP items
from pythonlibs.uniqueness import get_random_string_tag
from lib_tinyexononorf import scan_orf_for_tiny_exon
from lib_fasta import parseFasta, writeMultiFasta
from lib_clustalw import clustalw, strip_alignment_for_exterior_gaps
from lib_stopwatch import StopWatch
from pacb import PacbP, PacbPORF
import pacb
from graphAbgp import conversion

def _get_cbg_hmmsearch_results(cbg, inputdict, org, prev, next,
    orf_must_have_start=False,max_intron_nt_length=500,
    return_hmm_query_coord_dict=False,
    hmmsearch_num_hits=4,verbose=False):
    """
    Get HMM hits on ORFs for this organism, based in a HMM of a (partial/incomplete) CodingBlockGraph

    @type  cbg: CodingBlockGraph 
    @param cbg: CodingBlockGraph instance for which a HMM will be constructed

    @type  inputdict: dict 
    @param inputdict: <input data structure> 

    @type  org: * (presumably string)
    @param org: Organism identifier that is recognizable in the <input data structure>

    @type  prev: CodingBlockGraph (or None)
    @param prev: CodingBlockGraph upstream/5p of the cbg that must be completed

    @type  next: CodingBlockGraph (or None)
    @param next: CodingBlockGraph downstream/3p of the cbg that must be completed

    @attention: `prev` and `next` CodingBlockGraphs reduce the search space of ORFs
                 to scan with hmm. This speeds up & improves the quality of the outcome!

    @type  orf_must_have_start: Boolean
    @param orf_must_have_start: must the (new/completed) CodingBlockGraph cbg have alignable Methionines?

    @type  max_intron_nt_length: integer
    @param max_intron_nt_length: positive maximum intron length to take into acount when selecting suitable ORFs

    @type  hmmsearch_num_hits: integer
    @param hmmsearch_num_hits: (small) positive number of best HMM hits to take into account

    @type  verbose: Boolean
    @param verbose: report an extended debugging-report to STDOUT (True) or be quiet (False)
    """

    # try to limit searchspace by prev and next CBG
    prevNode, nextNode = None, None
    prevMin,  nextMax  = None, None
    maskcoords = []
    # check if (informant) organism is in the prev CBG AND if this CBG
    # has an OMSR -> not per se the case!
    if prev and org in prev.organism_set() and prev.has_overall_minimal_spanning_range():
        prevNode = prev.node_by_organism(org)
        omsr = prev.overall_minimal_spanning_range(organism=org)
        prevMin = (max(omsr)+1)*3
        maskcoords.append( ( 0, max(omsr) ) )
    # check if (informant) organism is in the next CBG AND if this CBG
    # has an OMSR -> not per se the case!
    if next and org in next.organism_set() and next.has_overall_minimal_spanning_range():
        nextNode = next.node_by_organism(org)
        omsr = next.overall_minimal_spanning_range(organism=org)
        nextMax = min(omsr)*3
        aaseqlen = len(inputdict[org]['genomeseq'])/3
        maskcoords.append( ( min(omsr), aaseqlen ) )
    if not prev and next and nextMax:
        prevMin = nextMax - max_intron_nt_length
    if not next and prev and prevMin:
        nextMax = prevMin + max_intron_nt_length 


    # make unique filename for hmm database and hmmbuild file
    fname_base = "%s_%s_range_%s_%s" % (
            org,
            get_random_string_tag(),
            prevMin, nextMax
            )
    fname_fasta_db = "%s_database.fa" % fname_base
    fname_fasta_alignment = "%s.fa" % fname_base


    # get elegiable sets of orfs from prev and next
    if not orf_must_have_start:
        elegiable_orfs = inputdict[org]['orfs'].get_elegiable_orfs(
                min_orf_end = prevMin,
                max_orf_start = nextMax
                )
    else:
        # ORFs *must* have starts => searching for a TSS exon/CBG
        elegiable_orfs = inputdict[org]['orfs'].get_elegiable_orfs(
                min_orf_end = prevMin,
                max_orf_start = nextMax,
                has_starts=True
                )

    ####################################################################
    if verbose:
        if len(elegiable_orfs) > 10:
            print "hmm-elegibable orfs:", org, len(elegiable_orfs), "/",
            print len(inputdict[org]['orfs'].orfs), "prevMin:", prevMin,
            print "nextMax:", nextMax
        else:
            print "hmm-elegibable orfs:", org, [ orf.id for orf in elegiable_orfs ], "/",
            print len(inputdict[org]['orfs'].orfs) 
    ####################################################################

    # check elegiable orfs; can be zero in case of a very tiny region to check
    if not elegiable_orfs:
        if return_hmm_query_coord_dict:
            # return the hmm hits & the basename AND the hmm profile coordinates!
            return [], fname_base, {}
        else:
            # return the hmm hits & the basename of the hmm files for later cleanup
            return [], fname_base


    # write multiple alignment to file
    if return_hmm_query_coord_dict:
        # retrieve the hmm profile multiplealignment coordinates too
        fastastring, hmmcoords = cbg.multiplealignment_and_coords()
        fastaseqs = parseFasta( fastastring.split("\n") )
    else:
        fastaseqs = parseFasta( cbg.multiplealignment().split("\n") )

    if not fastaseqs:
        # no fasta multiple alignment left -> return no results!
        if return_hmm_query_coord_dict:
            # return the hmm hits & the basename AND the hmm profile coordinates!
            return [], fname_base, hmmcoords
        else:
            # return the hmm hits & the basename of the hmm files for later cleanup
            return [], fname_base


    # write multiple alignment input file
    writeMultiFasta(fastaseqs,fname_fasta_alignment)

    # write masked orfs to fasta database
    db_fasta = inputdict[org]['orfs'].tomaskedfasta(coords=maskcoords,orflist=elegiable_orfs,header_prefix=org) 
    if orf_must_have_start:
        if len(db_fasta.strip()) == 0:
            # no UNmasked suitable ORFs remaining!
            # This is recognized lateron in this function 
            pass
        else:
            # mask out all AAs before the first start
            lines = db_fasta.split("\n")
            for linenr in range(0,len(lines)):
                line = lines[linenr]
                if line[0] != ">":
                    mpos = line.find("M")
                    if mpos > 0:
                        line = "X"*mpos+line[mpos:]
                    lines[linenr] = line
            # recreate db_fasta string
            db_fasta = "\n".join(lines)

    # write masked orfs to fasta database
    fh = open(fname_fasta_db,'w')
    fh.write( db_fasta )
    fh.close()

    # make shure that there where orfs written to file;
    # in case very little orfs are selected and all are masked -> no files!
    seqs_in_db = parseFasta(open(fname_fasta_db).readlines())
    if not seqs_in_db:
        # return no results!
        if return_hmm_query_coord_dict:
            # return the hmm hits & the basename AND the hmm profile coordinates!
            return [], fname_base, hmmcoords
        else:
            # return the hmm hits & the basename of the hmm files for later cleanup
            return [], fname_base

    # make hmmbuild file of the multiplealignment
    fname_hmmbuild = hmmbuild_protein( fname_fasta_alignment )

    if fname_hmmbuild:
        # run hmmsearch on this fasta database with hmmbuild file
        results = hmmsearch_protein( fname_hmmbuild, fname_fasta_db, params= {'E':150,'A': hmmsearch_num_hits} )
    else:
        # no fname_hmmbuild file => creation failed for some reason
        results = []

    if return_hmm_query_coord_dict:
        # return the hmm hits & the basename AND the hmm profile coordinates!
        return results, fname_base, hmmcoords
    else:
        # return the hmm hits & the basename of the hmm files for later cleanup
        return results, fname_base

# end of function _get_cbg_hmmsearch_results


def _hmmcleanup(fname_base):
    """
    """
    for extention in ['.fa','.hmmbuild','.output','_database.fa']:
        if os.path.exists(fname_base+extention):
            os.remove(fname_base+extention)

# end of function _hmmcleanup


def cbgtinytssexonhmmsearch(cbg,inputdict,prev=None,next=None,verbose=False):
    """
    Print data about potential hits!

    @type  cbg: CodingBlockGraph
    @param cbg: CodingBlockGraph to complete by hmmsearch

    @type  inputdict: dict
    @param inputdict: ``input data structure`` that contains lists of orfs

    @type  prev: CodingBlockGraph or None
    @param prev: previous CodingBlockGraph in an ordered genestructure

    @type  next: CodingBlockGraph or None
    @param next: next CodingBlockGraph in an ordered genestructure

    @type  verbose: Boolean
    @param verbose: report an extended debugging-report to STDOUT (True) or be quiet (False)

    @rtype:  list 
    @return: list of completed possible tss exons or the single input cbg 
    """
    # first some checks on the cbg; can it serve asa tiny start CBG?
    if not cbg.have_all_starts():
        return [ cbg ]

    # make a backup if the original one must be returned
    # and then clear_cache of the CBG
    bckp_cbg = deepcopy(cbg)
    cbg.clear_cache()

    # decrease minimal OMSR size
    cbg.MINIMAL_OVERAL_SPANNING_RANGE_SIZE = 2

    # trim the pacbps to the first occurence of a Methionine
    for pacbporf in cbg.pacbps.values():

        posObj = pacbporf._get_original_alignment_pos_start()
        while 'M' not in [ posObj.query, posObj.sbjct ] and pacbporf._original_alignment_pos_start < pacbporf._original_alignment_pos_end-1: 
            pacbporf._original_alignment_pos_start+=1
            posObj = pacbporf._get_original_alignment_pos_start()

        # TODO: rewrite this truncation of pacbps by first checking if there is a Methionine in the OMSR region!!

        if pacbporf._original_alignment_pos_end - pacbporf._original_alignment_pos_start < 2:
            # pacbporf is gone! Return the input
            return [ bckp_cbg ]

        # no problems; pacbporf truncated, so rescore
        pacbporf._score_alignment()

    # cgeck - double-check: does cbg still have an OMSR?
    if not cbg.has_overall_minimal_spanning_range():
        return [ bckp_cbg ]

    ############################################################
    if verbose:
        print "tss-hmm-cbg:", cbg
        print cbg.printmultiplealignment()
        for pacbp in cbg.pacbps.values():
            print pacbp
    ############################################################

    # get a collection of suitable pacbps based on the hmm from this cbg
    pacbpcol = cbghmmsearch2pacbpcollection(cbg,inputdict,prev=prev,next=next,
        orf_must_have_start=True,
        max_intron_nt_length=150,
        pacbp_min_length=None,
        pacbp_min_bitscore=None,
        hmmsearch_num_hits=3,verbose=verbose,genetreegraph=None)

    ############################################################
    if verbose: print "PacbpColHMM nodes:", pacbpcol.get_nodes()
    ############################################################

    # get list of accepted CBGs
    accepted = []
    for newcbg in conversion.pacbpCollection2AcceptedCodingBlockGraphs(pacbpcol,next=next):
        if newcbg.have_all_starts_upstream_of_omsr():
            accepted.append(newcbg)
            if verbose: print "newtsshmmCBG", newcbg

    # check what to return; are the accepted cbgs?
    if accepted:
        return accepted
    else:
        return [ bckp_cbg ]

# end of function cbgtinytssexonhmmsearch 


def cbghmmsearch(cbg,inputdict,prev=None,next=None,max_intron_nt_length=500,verbose=False):
    """
    Print data about potential hits!

    @type  cbg: CodingBlockGraph
    @param cbg: CodingBlockGraph to complete by hmmsearch

    @type  inputdict: dict
    @param inputdict: ``input data structure`` that contains lists of orfs

    @type  prev: CodingBlockGraph or None
    @param prev: previous CodingBlockGraph in an ordered genestructure

    @type  next: CodingBlockGraph or None
    @param next: next CodingBlockGraph in an ordered genestructure

    @type  verbose: Boolean
    @param verbose: report an extended debugging-report to STDOUT (True) or be quiet (False)

    @rtype:  None
    @return: None
    """

    # make a backup if the original one must be returned
    # and then clear_cache of the CBG
    bckp_cbg = deepcopy(cbg)
    cbg.clear_cache()
    for org in inputdict.keys():
        # ignore orgs/nodes that are found already
        if org in cbg.organism_set(): continue

        # get hmm results for this organism
        results, fname_base = _get_cbg_hmmsearch_results(cbg,inputdict,org,prev,next,max_intron_nt_length=max_intron_nt_length)

        # reorder hits based on surrounding cbgs
        results = _reorder_hits_on_occurrence_in_surrounding_cbgs(results,org,prev,next)

        for hit in results:
            sbjct_header, sbjct_start, sbjct_end, query_start, query_end, query, match, sbjct, score, expect = hit
            orfid = int(sbjct_header.split("_")[-1])
            hmmnode = (org,orfid)
            ##thisorf = inputdict[org]['orfs'].orfs[orfid]
            thisorf = inputdict[org]['orfs'].get_orf_by_id(orfid)


            ##sbjct_aa_start  = sbjct_start - 1 + thisorf.protein_startPY  # correct for -1 offset
            ##sbjct_aa_end    = sbjct_end + thisorf.protein_startPY
            ##sbjct_dna_start = thisorf.aapos2dnapos(sbjct_aa_start)
            ##sbjct_dna_end   = thisorf.aapos2dnapos(sbjct_aa_end)
            ##print org, score, expect
            ##print org, score, expect, thisorf, thisorf.id, "startAA:", thisorf.protein_startPY
            ##print query, query_start, query_end
            ##print match
            ##print sbjct, sbjct_dna_start, sbjct_dna_end, "\t\t", (sbjct_aa_start, sbjct_aa_end), "\t\t", prevMin,  nextMax

            # expand the CBG with this orf of an additional organism
            #print cbg
            #cbg.printmultiplealignment()

            # expand the CBG with this orf of an additional organism
            cbg = hmmfoundorf2cbg(cbg,thisorf,org,hit,verbose=verbose)
            if verbose: print "HMMincrease:", cbg
            break

        # remove the created files!
        _hmmcleanup(fname_base)

        if cbg.has_overall_minimal_spanning_range():
            # create a new bckp
            bckp_cbg = deepcopy(cbg)
        else:
            return bckp_cbg


    # and return this (increased) cbg
    if cbg.has_overall_minimal_spanning_range():
        if verbose: print "cbghmmsearchOUTPUT", cbg
        return cbg
    else:
        if verbose: print "UNCHANGED OUTPUT", bckp_cbg
        return bckp_cbg

# end of function cbghmmsearch


def _reorder_hits_on_occurrence_in_surrounding_cbgs(results,org,prev,next):
    """
    """
    reorder = []
    for pos in range(0,len(results)):
        orfid = int(results[pos][0].split("_")[-1])
        node = (org,orfid)
        if prev and next and node in prev.get_nodes() and node in next.get_nodes():
            # present in both CBGs! Unlikely, but can happen
            reorder.append(pos)
        elif prev and node in prev.get_nodes():
            reorder.append(pos)
        elif next and node in next.get_nodes():
            reorder.append(pos)
        else:
            pass
    # check if reordering is needed
    if reorder:
        reordered_results = []
        for pos in reorder:
            reordered_results.append( results[pos] )
        for pos in range(0,len(results)):
            if pos not in reorder:
                reordered_results.append( results[pos] )
        # and set reordered_results to results!
        results = reordered_results

    # return (reordered) results
    return results

# end of function _reorder_hits_on_occurrence_in_surrounding_cbgs



def tinyexonhmmsearch(cbg,inputdict,prev=None,next=None,max_intron_nt_length=500,verbose=False):
    """
    Complete a CodingBlockGraph that misses a (tiny) alignment from on or more organisms

    @type  cbg: CodingBlockGraph
    @param cbg: CodingBlockGraph to complete by hmmsearch

    @type  inputdict: dict
    @param inputdict: ``input data structure`` that contains lists of orfs

    @type  prev: CodingBlockGraph or None
    @param prev: previous CodingBlockGraph in an ordered genestructure

    @type  next: CodingBlockGraph or None
    @param next: next CodingBlockGraph in an ordered genestructure

    @type  verbose: Boolean
    @param verbose: report an extended debugging-report to STDOUT (True) or be quiet (False)

    @rtype:  []
    @return: list of CodingBlockGraph instance
    """

    # make a backup if the original one must be returned
    # and then clear_cache of the CBG
    bckp_cbg = deepcopy(cbg)
    cbg.clear_cache()

    pacbpcol = cbghmmsearch2pacbpcollection(cbg,inputdict,prev=prev,next=next,
        max_intron_nt_length=max_intron_nt_length,
        pacbp_min_length=None,
        pacbp_min_bitscore=None,
        hmmsearch_num_hits=3,verbose=verbose,genetreegraph=None)

    if verbose: print "PacbpColHMM nodes:", pacbpcol.get_nodes()

    # get list of accepted CBGs
    accepted = []
    for newcbg in conversion.pacbpCollection2AcceptedCodingBlockGraphs(pacbpcol,next=next):
        accepted.append(newcbg)
        if verbose: print "newTINYexonCBG:", newcbg

    # check what to return; are the accepted cbgs?
    if accepted:
        return accepted
    else:
        return [ bckp_cbg ]

# end of function tinyexonhmmsearch



def OLD_WORKING_tinyexonhmmsearch(cbg,inputdict,prev=None,next=None,max_intron_nt_length=500,verbose=False):
    """
    Complete a CodingBlockGraph that misses a (tiny) alignment from on or more organisms

    @type  cbg: CodingBlockGraph
    @param cbg: CodingBlockGraph to complete by hmmsearch

    @type  inputdict: dict
    @param inputdict: ``input data structure`` that contains lists of orfs

    @type  prev: CodingBlockGraph or None
    @param prev: previous CodingBlockGraph in an ordered genestructure

    @type  next: CodingBlockGraph or None
    @param next: next CodingBlockGraph in an ordered genestructure

    @rtype:  CodingBlockGraph
    @return: CodingBlockGraph instance
    """

    # make a backup if the original one must be returned
    # and then clear_cache of the CBG
    bckp_cbg = deepcopy(cbg)
    cbg.clear_cache()

    for org in inputdict.keys():
        # ignore orgs/nodes that are found already
        if org in cbg.organism_set(): continue

        # get hmm results for this organism
        results, fname_base = _get_cbg_hmmsearch_results(cbg,inputdict,org,prev,next,max_intron_nt_length=max_intron_nt_length)

        # reorder hits based on surrounding cbgs
        results = _reorder_hits_on_occurrence_in_surrounding_cbgs(results,org,prev,next)

        for hit in results:
            print 'HMM-HIT', hit[0], hit[-2], hit[-1]


        cbg_extended_with_tinyexon = False
        # get coordinates of previous and next codingblock
        # this in order to compare with the intron/exon prediction
        prevNode, nextNode = None, None
        prevMin,  nextMax  = None, None
        if prev and org in prev.organism_set():
            prevNode = prev.node_by_organism(org)
            omsr = prev.overall_minimal_spanning_range(organism=org)
            prevMin = (max(omsr)+1)*3
        if next and org in next.organism_set():
            nextNode = next.node_by_organism(org)
            omsr = next.overall_minimal_spanning_range(organism=org)
            nextMax = min(omsr)*3

        for hit in results:
            if cbg_extended_with_tinyexon: break
            # get orf, node and AA and DNA coordinates of this sbjct hit
            sbjct_header, sbjct_start, sbjct_end, query_start, query_end, query, match, sbjct, score, expect = hit
            orfid = int(sbjct_header.split("_")[-1])
            hmmnode = (org,orfid)
            ##thisorf = inputdict[org]['orfs'].orfs[orfid]
            thisorf = inputdict[org]['orfs'].get_orf_by_id(orfid)

            sbjct_aa_start  = sbjct_start - 1 + thisorf.protein_startPY  # correct for -1 offset
            sbjct_aa_end    = sbjct_end + thisorf.protein_startPY
            sbjct_dna_start = thisorf.aapos2dnapos(sbjct_aa_start)
            sbjct_dna_end   = thisorf.aapos2dnapos(sbjct_aa_end)

            if prev and next and hmmnode in prev.get_nodes() and hmmnode in next.get_nodes():
                # tinyexon bridges both CBGs! A very weird case, possible in case of
                # a non-detected lowsimilarityregion around an intron position
                # no tinyexon searches with introns needed -> it is not an tinyexon!
                tinyexons = [ True ] # Hack; actual data IN tinyexons is not needed
                pass

            elif prev and hmmnode in prev.get_nodes():
                # hypothetical tinyexon is an extension of the previous CBG
                # find elegiable tinyexons by limiting to size range and approximate position
                max_tinyexon_nt_length = (sbjct_dna_end - sbjct_dna_start) + 12
                tinyexons = scan_orf_for_tiny_exon(thisorf,max_tinyexon_nt_length=max_tinyexon_nt_length,
                        min_acceptor_pos=None,max_donor_pos=nextMax)
                # if no tinyexons -> no introns possible -> continue
                if not tinyexons: continue
 
            elif next and hmmnode in next.get_nodes():
                # hypothetical tinyexon is an extension of the next CBG
                # hypothetical tinyexon is an extension of the previous CBG
                # find elegiable tinyexons by limiting to size range and approximate position
                max_tinyexon_nt_length = (sbjct_dna_end - sbjct_dna_start) + 12
                tinyexons = scan_orf_for_tiny_exon(thisorf,max_tinyexon_nt_length=max_tinyexon_nt_length,
                        min_acceptor_pos=prevMin,max_donor_pos=None)
                # if no tinyexons -> no introns possible -> continue
                if not tinyexons: continue

            else:
                # find elegiable tinyexons by limiting to size range and approximate position
                max_tinyexon_nt_length = (sbjct_dna_end - sbjct_dna_start) + 12
                tinyexons = scan_orf_for_tiny_exon(thisorf,max_tinyexon_nt_length=max_tinyexon_nt_length,
                        min_acceptor_pos=prevMin,max_donor_pos=nextMax)
                # if no tinyexons -> no introns possible -> continue 
                if not tinyexons: continue 

            #print hit
            #print thisorf, thisorf.protein_startPY
            #print thisorf.protein_sequence, "tinyexons:", len(tinyexons)
            #print sbjct_dna_start, sbjct_dna_end, "\t\t", (sbjct_aa_start, sbjct_aa_end), "\t\t", prevMin,  nextMax
            #for te in tinyexons:
            #    print te.acceptor.pos, te.donor.pos, te

            # if there are tinyexons -> we have found a CBG with tinyexon
            if tinyexons:
                # expand the CBG with this tinyexon of an additional organism
                cbg = hmmfoundorf2cbg(cbg,thisorf,org,hit,verbose=verbose)
                cbg_extended_with_tinyexon = True
                print cbg
                break


        # remove the created files!
        _hmmcleanup(fname_base)

        if cbg.has_overall_minimal_spanning_range():
            bckp_cbg = deepcopy(cbg)
        else:
            # no OMSR after this step -> done here!
            return bckp_cbg 

    # and return this (increased) cbg
    if cbg.has_overall_minimal_spanning_range():
        cbg.create_cache()
        print "cbghmmsearchOUTPUT", cbg
        print cbg.genetree()
        return cbg
    else:
        print "UNCHANGED OUTPUT", bckp_cbg
        return bckp_cbg

# end of function tinyexonhmmsearch



def hmmfoundorf2cbg(cbg,neworf,neworg,hmmhit,verbose=False):
    """

    @type  verbose: Boolean
    @param verbose: report an extended debugging-report to STDOUT (True) or be quiet (False)
    """
    # trim hmmhit for unmatched characters
    sbjct_header, sbjct_start, sbjct_end, query_start, query_end, query, match, sbjct, score, expect = hmmhit
    while match and match[0] == ' ':
        query = query[1:]
        match = match[1:]
        sbjct = sbjct[1:]
        sbjct_start+=1
        query_start+=1
    while match and match[-1] == ' ':
        query = query[0:-1]
        match = match[0:-1]
        sbjct = sbjct[0:-1]
        sbjct_end-=1
        query_end-=1
    
    # get orf, node and AA and DNA coordinates of this sbjct hit
    sbjct_aa_start  = sbjct_start - 1 + neworf.protein_startPY  # correct for -1 offset
    sbjct_aa_end    = sbjct_end + neworf.protein_startPY
    sbjct_dna_start = neworf.aapos2dnapos(sbjct_aa_start)
    sbjct_dna_end   = neworf.aapos2dnapos(sbjct_aa_end)
    sbjctNode       = (neworg,neworf.id)
    query = query.replace(".","-")
    sbjct = sbjct.replace(".","-").upper()
    # get the input hmm sequences
    fastaseqs = parseFasta( cbg.multiplealignment().split("\n") )
    org2seq = {}
    for h,seq in fastaseqs.iteritems():
        org,etc = h.split('_orf_')
        org2seq[org] = seq

    # find gap positions in the hmm query
    gappos = []
    if query.find("-") >= 0:
        gpos = 0
        while True:
            gpos = query.find("-",gpos)
            if gpos >= 0:
                gappos.append(gpos)
                gpos+=1
            else:
                break
    if verbose:
        print "hmmQ:", query, gappos, query_start, query_end, "gaps:", query.count('-'), len(query)
        print "hmmM:", match
        print "hmmS:", sbjct, sbjctNode, sbjct_aa_start, sbjct_aa_end, "len:", sbjct_aa_end-sbjct_aa_start , len(sbjct)
    newpacbporfs = {}
    allomsr = cbg.overall_minimal_spanning_range() 
    pacbps_are_created = True
    for org in cbg.organism_set():
        orf  = cbg.get_orfs_of_graph(organism=org)[0]
        queryNode = (org,orf.id)
        omsr = allomsr[queryNode]

        # get querysequence from multiplealignment fasta
        queryseq = org2seq[org]
        queryseq = queryseq[ query_start-1 : query_end ]
        # project gaps with the hmmbuild into the sequence
        for gpos in gappos:
            queryseq = queryseq[0:gpos]+"-"+queryseq[gpos:]

        # find query sequence position on Orf
        query_aa_start = orf.protein_startPY + orf.protein_sequence.find(queryseq.replace('-',''))
        query_aa_end   = query_aa_start + len(queryseq) - queryseq.count('-')

        combi = [ queryNode, sbjctNode ]
        combi.sort()
        if verbose:
            print "hmmq:", queryseq, queryNode, query_aa_start, query_aa_end, "len:", query_aa_end-query_aa_start, len(queryseq), 
        # make a deepcopy; sbjct is needed unchanged for the next iteration
        # in the for loop, but here we want to trim of gap sequences
        sbjctseq = deepcopy(sbjct)
        sbjctaastart = deepcopy(sbjct_aa_start)
        sbjctaaend   = deepcopy(sbjct_aa_end)
        while queryseq and queryseq[0] == '-':
            queryseq = queryseq[1:]
            sbjctseq = sbjctseq[1:]
            sbjctaastart+=1
        while sbjctseq and sbjctseq[0] == '-':
            queryseq = queryseq[1:]
            sbjctseq = sbjctseq[1:]
            query_aa_start+=1
        while queryseq and queryseq[-1] == '-':
            queryseq = queryseq[0:-1]
            sbjctseq = sbjctseq[0:-1]
            sbjctaaend-=1
        while sbjctseq and sbjctseq[-1] == '-':
            queryseq = queryseq[0:-1]
            sbjctseq = sbjctseq[0:-1]
            query_aa_end-=1
        if verbose:
            print queryseq, sbjctseq, (query_aa_start, query_aa_end), (sbjctaastart,sbjctaaend)
        if queryseq and sbjctseq:
            if combi == [ queryNode, sbjctNode ]:
                pacbpinput = (queryseq,sbjctseq,query_aa_start,sbjctaastart)
                pacbp      = PacbP(input=pacbpinput)
                pacbp.source = 'hmmsearch'
                pacbporf   = PacbPORF(pacbp,orf,neworf)
                pacbporf.extend_pacbporf_after_stops()
                newpacbporfs[tuple(combi)] = pacbporf
            else:
                pacbpinput = (sbjctseq,queryseq,sbjctaastart,query_aa_start)
                pacbp      = PacbP(input=pacbpinput)
                pacbp.source = 'hmmsearch'
                pacbporf   = PacbPORF(pacbp,neworf,orf)
                pacbporf.extend_pacbporf_after_stops()
                newpacbporfs[tuple(combi)] = pacbporf
        else:
            pacbps_are_created = False
            break

    if pacbps_are_created:
        # now add the newpacbporfs to the graph
        cbg.add_node(sbjctNode)
        for (n1,n2), pacbporf in newpacbporfs.iteritems():
            wt = pacbporf.bitscore
            pacbpkey = pacbporf.construct_unique_key(n1,n2)
            cbg.add_edge(n1,n2,wt=wt)
            cbg.pacbps[(pacbpkey,n1,n2)] = pacbporf
            if verbose:
                print pacbporf
    # and return this increased cbg
    return cbg

# end of function hmmfoundorf2cbg 




def cbghmmsearch2pacbpcollection(cbg,inputdict,prev=None,next=None,
    orf_must_have_start=False,
    max_intron_nt_length=500,
    pacbp_min_length=None,
    pacbp_min_bitscore=None,
    hmmsearch_num_hits=3,verbose=False,genetreegraph=None):
    """
    Print data about potential hits!

    @type  cbg: CodingBlockGraph
    @param cbg: CodingBlockGraph to complete by hmmsearch

    @type  inputdict: dict
    @param inputdict: ``input data structure`` (that contains lists of orfs)

    @type  prev: CodingBlockGraph or None
    @param prev: previous CodingBlockGraph in an ordered genestructure

    @type  next: CodingBlockGraph or None
    @param next: next CodingBlockGraph in an ordered genestructure

    @type  max_intron_nt_length: (positive) integer
    @param max_intron_nt_length: TODO

    @type  pacbp_min_length: (positive) integer or None
    @param pacbp_min_length: TODO

    @type  pacbp_min_bitscore: integer or None
    @param pacbp_min_bitscore: TODO

    @type  hmmsearch_num_hits: (small positive) integer 
    @param hmmsearch_num_hits: TODO

    @type  genetreegraph: GeneTreeGraph object
    @param genetreegraph: TODO

    @type  verbose: Boolean 
    @param verbose: report an extended debugging-report to STDOUT (True) or be quiet (False) 

    @rtype:  PacbpCollectionGraph
    @return: PacbpCollectionGraph instance
    """

    # make a PacbpCollectionGraph of the CBG
    pcg = conversion.CodingBlockGraph2PacbpCollectionGraph(cbg)

    hmmsearched_orgs = []
    hmmfound_nodes   = {}
    for org in inputdict.keys():
        # ignore orgs/nodes that are found already
        if org in cbg.organism_set(): continue
        hmmsearched_orgs.append(org)

        ############################################################################
        if verbose:
            fn = "cbghmmsearch2pacbpcollection"
            print fn+" about to do HMM", org
            print fn+" pacbp_min_length   :", pacbp_min_length
            print fn+" pacbp_min_bitscore :", pacbp_min_bitscore
            print cbg
        ############################################################################

        # get hmm results for this organism
        results, fname_base = _get_cbg_hmmsearch_results(cbg,inputdict,org,prev,next,
                orf_must_have_start=orf_must_have_start,
                max_intron_nt_length=max_intron_nt_length,
                hmmsearch_num_hits=hmmsearch_num_hits,
                verbose=verbose
                )
        # reordering of hits is not needed, because all hits are proccessed!

        if verbose and orf_must_have_start:
            for hit in results: print 'HMM-HIT', hit[0], hit[-2], hit[-1],
            print ''


        if verbose and results:
            for hit in results: print 'HMM-HIT', hit[0], hit[-2], hit[-1],
            print ''
        elif verbose:
            print "NO HMM HITS!!"
        else:
            pass 

        for hit in results:
            # translate tuple `hit` into separate variables
            sbjct_header, sbjct_start, sbjct_end, query_start, query_end, query, match, sbjct, score, expect = hit
            # get orfid from sbjct_header and get the Orf object
            orfid   = int(sbjct_header.split("_")[-1])
            ##thisorf = inputdict[org]['orfs'].orfs[orfid]
            thisorf = inputdict[org]['orfs'].get_orf_by_id(orfid)

            # get fasta-like representation of this sequence
            aastart = thisorf.protein_startPY + sbjct_start - 1
            aaend   = thisorf.protein_startPY + sbjct_end
            header  = "%s_orf_%s_%s_%s" % (org,orfid,aastart,aaend)
            seq     = thisorf.getaas(abs_pos_start=aastart,abs_pos_end=aaend)
            # store fasta representation in hmmfound_nodes
            hmmfound_nodes[header] = seq

            # create pacbporf of this CBG and the hmm hits
            pacbporfs = hmmhit2pacbps(cbg,thisorf,org,hit,verbose=verbose)

            # add new nodes & edges of pacbporfs to pacbpcollection
            for ((wt,length,orfid1,orfid2), n1, n2),pacbporf in pacbporfs.iteritems():
                if pacbp_min_length and len(pacbporf) < pacbp_min_length: continue
                if pacbp_min_bitscore and pacbporf.bitscore < pacbp_min_bitscore: continue
                if n1 not in pcg.get_nodes(): pcg.add_node(n1)
                if n2 not in pcg.get_nodes(): pcg.add_node(n2)
                pcg.add_edge(n1,n2,wt=wt)
                if verbose:
                    print "HMM:", n1, n2, pacbporf 
                # add the pacbporf themselve to pcg.pacbps
                pcg.pacbps[( (wt,length,orfid1,orfid2), n1, n2 )] = pacbporf

        # remove the created files!
        _hmmcleanup(fname_base)

    # create regular expression pattern for hmmfound nodes
    hmmfound_nodes_regex_pattern = re_compile("^(\w+)_orf_(\d+)_(\d+)_(\d+)$")

    # make pacbporfs between the hmmhits of distinct organisms
    if len(hmmsearched_orgs) > 1:
        for header1,seq1 in hmmfound_nodes.iteritems():
            org1,orfid1,aa1start,aa1end = re_match(hmmfound_nodes_regex_pattern,header1).groups()
            aa1start, aa1end, orfid1 = int(aa1start), int(aa1end), int(orfid1)
            node1 = (org1,orfid1)
            # check if node1 is indeed incorporated in the PacbpCollection
            if node1 not in pcg.get_nodes(): continue
            # get orf object of first node
            orf1 = inputdict[org1]['orfs'].get_orf_by_id(orfid1)

            for header2,seq2 in hmmfound_nodes.iteritems():
                # ignore identical headers
                if header1 == header2: continue 
                org2,orfid2,aa2start,aa2end = re_match(hmmfound_nodes_regex_pattern,header2).groups()
                aa2start, aa2end, orfid2 = int(aa2start), int(aa2end), int(orfid2)
                node2 = (org2,orfid2)
                # check if node2 is indeed incorporated in the PacbpCollection
                if node2 not in pcg.get_nodes(): continue
                # ignore identical organisms
                if org1 == org2: continue

                # do only the diagonal of the cross
                edge = [node1,node2]
                edge.sort()
                if edge != [node1,node2]: continue

                # get orf object of second node
                ##orf2 = inputdict[org2]['orfs'].orfs[orfid2]
                orf2 = inputdict[org2]['orfs'].get_orf_by_id(orfid2)

                # align the sequences with clustalw
                seqs = { header1: seq1, header2: seq2 }
                (alignedseqs,alignment) = clustalw(seqs=seqs)
                # make pacbp from clustalw alignment
                pacbp = pacb.conversion.pacbp_from_clustalw(
                        alignment=(
                                alignedseqs[header1],
                                alignment,
                                alignedseqs[header2]
                                ),
                        coords=(aa1start,aa1end,aa2start,aa2end)
                        )

                # only accept pacbps that forfill the requirements
                if not pacbp: continue
                if pacbp_min_length and len(pacbp) < pacbp_min_length: continue
                if pacbp_min_bitscore and pacbp.bitscore < pacbp_min_bitscore: continue

                # make pacbporf and extend
                pacbporf   = pacb.PacbPORF(pacbp,orf1,orf2)
                pacbporf.extend_pacbporf_after_stops()
                if verbose:
                    print node1, node2, "crossconnection after HMM"
                    print pacbporf

                # and add edge & pacbporf to pacbpcollection
                pcg.add_edge(node1,node2,wt=pacbporf.bitscore)
                key = pacbporf.construct_unique_key(node1,node2)
                pcg.pacbps[(key,node1,node2)] = pacbporf

    # return the pacbpcollection
    return pcg

# end of function cbghmmsearch2pacbpcollection



def hmmhit2pacbps(cbg,neworf,neworg,hmmhit,verbose=False):
    """
    """
    ####################################################################################
    if verbose:
        # start timing this complete operation
        stw = StopWatch('hmmhit2pacbps')
        stw.start()
    ####################################################################################

    # trim hmmhit for unmatched characters
    sbjct_header, sbjct_start, sbjct_end, query_start, query_end, query, match, sbjct, score, expect = hmmhit
    while match and match[0] == ' ':
        query = query[1:]
        match = match[1:]
        sbjct = sbjct[1:]
        sbjct_start+=1
        query_start+=1
    while match and match[-1] == ' ':
        query = query[0:-1]
        match = match[0:-1]
        sbjct = sbjct[0:-1]
        sbjct_end-=1
        query_end-=1
    
    # get orf, node and AA and DNA coordinates of this sbjct hit
    sbjct_aa_start  = sbjct_start - 1 + neworf.protein_startPY  # correct for -1 offset
    sbjct_aa_end    = sbjct_end + neworf.protein_startPY
    sbjctNode       = (neworg,neworf.id)
    query           = query.replace(".","-")
    sbjct           = sbjct.replace(".","-").upper()

    ####################################################################################
    if verbose:
        print "hmmhit2pacbps CREATING pacbps for organism/orf: (%s,%s)" % (
                neworg,neworf.id)
        print "hmmhit2pacbps Q '%s'" % query
        print "hmmhit2pacbps m '%s'" % match
        print "hmmhit2pacbps S '%s'" % sbjct
    ####################################################################################

    # get the input hmm sequences
    fastaseqs = parseFasta( cbg.multiplealignment().split("\n") )
    org2seq = {}
    for h,seq in fastaseqs.iteritems():
        org,etc = h.split('_orf_')
        org2seq[org] = seq
    if verbose: print "hmmhit2pacbps; org2seq.keys() %s lengths: %s" % (
            org2seq.keys(),
            [ len(_seq) for _seq in org2seq.values() ],
            )

    # find gap positions in the hmm query
    gappos = []
    if query.find("-") >= 0:
        gpos = 0
        while True:
            gpos = query.find("-",gpos)
            if gpos >= 0:
                gappos.append(gpos)
                gpos+=1
            else:
                break

    ####################################################################################
    if verbose:
        print "hmmhit2pacbps", stw.lap(), "preparations done"
        print "hmmQ:", query, gappos, query_start, query_end,
        print "gaps:", query.count('-'), len(query)
        print "hmmM:", match
        print "hmmS:", sbjct, sbjctNode, sbjct_aa_start, sbjct_aa_end,
        print "len:", sbjct_aa_end-sbjct_aa_start , len(sbjct)
    ####################################################################################

    # return dictionary with edges (keys) and new pacbporfs (pacbporfs)
    newpacbporfs = {}

    for org in cbg.organism_set():
        # when organism not part of the HMM search profile -> continue
        if org not in org2seq.keys(): continue

        # get Orf object of this organism identifier
        orf  = cbg.get_orfs_of_graph(organism=org)[0]

        # get querysequence from multiplealignment fasta
        queryNode = (org,orf.id)
        queryseq  = org2seq[org]
        queryseq  = queryseq[ query_start-1 : query_end ]
        # project gaps with the hmmbuild into the sequence
        for gpos in gappos:
            queryseq = queryseq[0:gpos]+"-"+queryseq[gpos:]

        # find query sequence position on Orf
        query_aa_start = orf.protein_startPY + orf.protein_sequence.find(queryseq.replace('-',''))
        query_aa_end   = query_aa_start + len(queryseq) - queryseq.count('-')

        # make node combi (the key in the return dictionary)
        combi = [ queryNode, sbjctNode ]
        combi.sort()

        ################################################################################
        if verbose:
            print "hmmq:", queryseq, queryNode, query_aa_start, query_aa_end,
            print "len:", query_aa_end-query_aa_start, len(queryseq), 
        ################################################################################

        # make a deepcopy; sbjct is needed unchanged for the next iteration
        # in the for loop, but here we want to trim of gap sequences
        sbjctseq = deepcopy(sbjct)
        sbjctaastart = deepcopy(sbjct_aa_start)
        sbjctaaend   = deepcopy(sbjct_aa_end)
        while queryseq and queryseq[0] == '-':
            queryseq = queryseq[1:]
            sbjctseq = sbjctseq[1:]
            sbjctaastart+=1
        while sbjctseq and sbjctseq[0] == '-':
            queryseq = queryseq[1:]
            sbjctseq = sbjctseq[1:]
            query_aa_start+=1
        while queryseq and queryseq[-1] == '-':
            queryseq = queryseq[0:-1]
            sbjctseq = sbjctseq[0:-1]
            sbjctaaend-=1
        while sbjctseq and sbjctseq[-1] == '-':
            queryseq = queryseq[0:-1]
            sbjctseq = sbjctseq[0:-1]
            query_aa_end-=1

        # check opposing gap positions which can be
        # the result of very gapped HMMSEARCH profiles
        for pos in range(len(sbjctseq)-1,-1,-1):
            if queryseq[pos] == '-' and sbjctseq[pos] == '-':
                queryseq = queryseq[0:pos]+queryseq[pos+1:]
                sbjctseq = sbjctseq[0:pos]+sbjctseq[pos+1:]

        if queryseq and sbjctseq:
            if combi == [ queryNode, sbjctNode ]:
                pacbpinput = (queryseq,sbjctseq,query_aa_start,sbjctaastart)
                pacbp      = PacbP(input=pacbpinput)
                pacbp.source = 'hmmsearch'
                pacbporf   = PacbPORF(pacbp,orf,neworf)
                pacbporf.extend_pacbporf_after_stops()
                pacbpkey = pacbporf.construct_unique_key(queryNode,sbjctNode)
                newpacbporfs[(pacbpkey,queryNode,sbjctNode)] = pacbporf
            else:
                pacbpinput = (sbjctseq,queryseq,sbjctaastart,query_aa_start)
                pacbp      = PacbP(input=pacbpinput)
                pacbp.source = 'hmmsearch'
                pacbporf   = PacbPORF(pacbp,neworf,orf)
                pacbporf.extend_pacbporf_after_stops()
                pacbpkey = pacbporf.construct_unique_key(sbjctNode,queryNode)
                newpacbporfs[(pacbpkey,sbjctNode,queryNode)] = pacbporf
        else:
            # not a query and a sbjct sequence. This will result lateron 
            # in non-complete CBGs, but here it is not a problem or an error
            pacbps_are_created = False
            pass

        ################################################################################
        if verbose: print stw.lap()
        ################################################################################

    ################################################################################
    if verbose:
        # function hmmhit2pacbps DONE!
        print "hmmhit2pacbps CREATED (%s)\n" % len(newpacbporfs)
    ################################################################################

    # return the created pacbporfs
    return newpacbporfs

# end of function hmmhit2pacbps 


def hmmhit2pacbp(queryorf,queryorg,querycoords,sbjctorf,sbjctorg,hmmhit,verbose=False):
    """
    """
    # trim hmmhit for unmatched characters
    ( sbjct_header, sbjct_start, sbjct_end,
      query_start, query_end,
      query, match, sbjct, score, expect ) = hmmhit

    while match and match[0] == ' ':
        query = query[1:]
        match = match[1:]
        sbjct = sbjct[1:]
        sbjct_start+=1
        query_start+=1
    while match and match[-1] == ' ':
        query = query[0:-1]
        match = match[0:-1]
        sbjct = sbjct[0:-1]
        sbjct_end-=1
        query_end-=1

    # get orf, node and AA and DNA coordinates of this sbjct hit;
    # correct for -1 offset in start coordinate!!
    sbjct_aa_start  = sbjct_start - 1 + sbjctorf.protein_startPY
    sbjct_aa_end    = sbjct_end + sbjctorf.protein_startPY
    sbjctNode       = (sbjctorg,sbjctorf.id)
    query           = query.replace(".","-").upper()
    sbjct           = sbjct.replace(".","-").upper()

    ############################################################################
    if verbose:
        print "hmmhit2pacbp CREATING pacbps for organism/orf: (%s,%s)" % (
                sbjctorg,sbjctorf.id)
        print "hmmhit2pacbp Q '%s'" % query
        print "hmmhit2pacbp m '%s'" % match
        print "hmmhit2pacbp S '%s'" % sbjct
        print "hmmQ:", query, query_start, query_end, "gaps:",
        print query.count('-'), len(query)
        print "hmmM:", match
        print "hmmS:", sbjct, sbjctNode, sbjct_aa_start, sbjct_aa_end,
        print "len:", sbjct_aa_end-sbjct_aa_start , len(sbjct)
    ############################################################################

    # get Node and sequence of the query
    queryNode = (queryorg,queryorf.id)
    queryseq  = deepcopy(query)

    # calculate query sequence position on queryorf
    query_aa_start = querycoords[0] + query_start - 1
    query_aa_end   = query_aa_start + len(queryseq) - queryseq.count('-')

    # make node combi (the key in the return dictionary)
    combi = [ queryNode, sbjctNode ]
    combi.sort()

    ############################################################################
    if verbose:
        print "hmmq:", queryseq, queryNode, query_aa_start, query_aa_end,
        print "len:", query_aa_end-query_aa_start, len(queryseq),
    ############################################################################

    # make a deepcopy; sbjct is needed unchanged for the next iteration
    # in the for loop, but here we want to trim of gap sequences
    sbjctseq = deepcopy(sbjct)
    sbjctaastart = deepcopy(sbjct_aa_start)
    sbjctaaend   = deepcopy(sbjct_aa_end)
    while queryseq and queryseq[0] == '-':
        queryseq = queryseq[1:]
        sbjctseq = sbjctseq[1:]
        sbjctaastart+=1
    while sbjctseq and sbjctseq[0] == '-':
        queryseq = queryseq[1:]
        sbjctseq = sbjctseq[1:]
        query_aa_start+=1
    while queryseq and queryseq[-1] == '-':
        queryseq = queryseq[0:-1]
        sbjctseq = sbjctseq[0:-1]
        sbjctaaend-=1
    while sbjctseq and sbjctseq[-1] == '-':
        queryseq = queryseq[0:-1]
        sbjctseq = sbjctseq[0:-1]
        query_aa_end-=1

    if queryseq and sbjctseq:
        ################################################################
        if len(queryseq) != len(sbjctseq):
            # this will result in a exception to be raised:
            # pacb.exceptions.InproperlyAppliedArgument
            # print data here about what went wrong, then
            # just let the error be raised
            print queryseq, len(queryseq), sbjctseq, len(sbjctseq)
            print hmmhit
            print "Q:", query_aa_start, query_aa_end,
            print query_aa_end - query_aa_start,
            print "S:", sbjctaastart, sbjctaaend,
            print sbjctaaend - sbjctaastart
        ################################################################
        if combi == [ queryNode, sbjctNode ]:
            pacbpinput = (queryseq,sbjctseq,query_aa_start,sbjctaastart)
            pacbp      = PacbP(input=pacbpinput)
            pacbp.source = 'hmmsearch'
            pacbporf   = PacbPORF(pacbp,queryorf,sbjctorf)
            pacbporf.strip_unmatched_ends()
            if pacbporf.length==0:
                # Pacbp creation failed!
                return False, None
            else:
                pacbporf.extend_pacbporf_after_stops()
                pacbpkey = pacbporf.construct_unique_key(queryNode,sbjctNode)
                # return unique key and pacbporf
                return (pacbpkey,queryNode,sbjctNode), pacbporf
        else:
            pacbpinput = (sbjctseq,queryseq,sbjctaastart,query_aa_start)
            pacbp      = PacbP(input=pacbpinput)
            pacbp.source = 'hmmsearch'
            pacbporf   = PacbPORF(pacbp,sbjctorf,queryorf)
            pacbporf.strip_unmatched_ends()
            if not pacbporf:
                # Pacbp creation failed!
                return False, None
            else:
                pacbporf.extend_pacbporf_after_stops()
                pacbpkey = pacbporf.construct_unique_key(sbjctNode,queryNode)
                # return unique key and pacbporf
                return (pacbpkey,sbjctNode,queryNode), pacbporf
    else:
        # Pacbp creation failed!
        return False, None

# end of function hmmhit2pacbp



def create_hmmdb_for_neighbouring_cbgs(
    inputdict,prev, next,
    organism=None,
    orf_must_have_start=False,
    omsr_2_mask_aa_length_correction=6,
    max_intron_nt_length=500,
    verbose=False):
    """

    TODO: just written; test this function AND replace it by removing it in old code!!

    """
    # list variable with fasta db strings per organism
    fasta_db_strings = []

    for org in prev.organism_set():
        # check if this organism is requested for or not
        if organism and org != organism: continue

        # try to limit searchspace by prev and next CBG
        prevNode, nextNode = None, None
        prevMin,  nextMax  = None, None
        maskcoords = []
        if prev and org in prev.organism_set():
            prevNode = prev.node_by_organism(org)
            omsr = prev.overall_minimal_spanning_range(organism=org)
            prevMin = (max(omsr)+1)*3
            maskcoords.append( ( 0, max(omsr)-omsr_2_mask_aa_length_correction ) )
        if next and org in next.organism_set():
            nextNode = next.node_by_organism(org)
            omsr = next.overall_minimal_spanning_range(organism=org)
            nextMax = min(omsr)*3
            aaseqlen = len(inputdict[org]['genomeseq'])/3
            maskcoords.append( ( min(omsr)+omsr_2_mask_aa_length_correction, aaseqlen ) )
        if not prev and next and nextMax:
            prevMin = nextMax - max_intron_nt_length
        if not next and prev and prevMin:
            nextNode = next.node_by_organism(org)
            omsr = next.overall_minimal_spanning_range(organism=org)
            nextMax = min(omsr)*3
            aaseqlen = len(inputdict[org]['genomeseq'])/3
            maskcoords.append( ( min(omsr), aaseqlen ) )
        if not prev and next and nextMax:
            prevMin = nextMax - max_intron_nt_length
        if not next and prev and prevMin:
            nextMax = prevMin + max_intron_nt_length


        # get elegiable sets of orfs from prev and next
        if not orf_must_have_start:
            elegiable_orfs = inputdict[org]['orfs'].get_elegiable_orfs(
                min_orf_end = prevMin,
                max_orf_start = nextMax
                )
        else:
            # ORFs *must* have starts => searching for a TSS exon/CBG
            elegiable_orfs = inputdict[org]['orfs'].get_elegiable_orfs(
                min_orf_end = prevMin,
                max_orf_start = nextMax,
                has_starts=True
                )

        ####################################################################
        if verbose:
            if len(elegiable_orfs) > 10:
                print "hmm-elegibable orfs:", org, len(elegiable_orfs), "/",
                print len(inputdict[org]['orfs'].orfs), "prevMin:", prevMin,
                print "nextMax:", nextMax
            else:
                print "hmm-elegibable orfs:", org, [ orf.id for orf in elegiable_orfs ], "/",
                print len(inputdict[org]['orfs'].orfs)
        ####################################################################

        # check elegiable orfs; can be zero in case of a very tiny region to check
        if not elegiable_orfs:
            continue

        # write masked orfs to fasta database
        # variable `db_fasta` isa string!
        db_fasta = inputdict[org]['orfs'].tomaskedfasta(coords=maskcoords,orflist=elegiable_orfs,header_prefix=org)
        if orf_must_have_start:
            if len(db_fasta.strip()) == 0:
                # no UNmasked suitable ORFs remaining!
                # This is recognized lateron in this function
                pass
            else:
                # mask out all AAs before the first start
                lines = db_fasta.split("\n")
                for linenr in range(0,len(lines)):
                    line = lines[linenr]
                    if line[0] != ">":
                        mpos = line.find("M")
                        if mpos > 0:
                            line = "X"*mpos+line[mpos:]
                        lines[linenr] = line
                # recreate db_fasta string
                db_fasta = "\n".join(lines)

        # append to fasta_db_strings list variable
        fasta_db_strings.append( db_fasta )

    # return the overall gathered database of (masked) ORFs
    return "\n".join( fasta_db_strings )

# end of function create_hmmdb_for_neighbouring_cbgs




def hmmsearch_results_for_inwpcbg(inwpcbg, inputdict, org, prev, next,
    orf_must_have_start=False,max_intron_nt_length=500,
    return_hmm_query_coord_dict=False,
    hmmsearch_num_hits=4,verbose=False):
    """
    Get HMM hits on ORFs for this organism, based in a HMM of a (partial/incomplete) CodingBlockGraph

    @type  inwpcbg: CodingBlockGraph 
    @param inwpcbg: CodingBlockGraph instance for which a HMM will be constructed

    @type  inputdict: dict 
    @param inputdict: <input data structure> 

    @type  org: * (presumably string)
    @param org: Organism identifier that is recognizable in the <input data structure>

    @type  prev: CodingBlockGraph (or None)
    @param prev: CodingBlockGraph upstream/5p of the inwpcbg that must be completed

    @type  next: CodingBlockGraph (or None)
    @param next: CodingBlockGraph downstream/3p of the inwpcbg that must be completed

    @attention: `prev` and `next` CodingBlockGraphs reduce the search space of ORFs
                 to scan with hmm. This speeds up & improves the quality of the outcome!

    @type  orf_must_have_start: Boolean
    @param orf_must_have_start: must the (new/completed) CodingBlockGraph inwpcbg have alignable Methionines?

    @type  max_intron_nt_length: integer
    @param max_intron_nt_length: positive maximum intron length to take into acount when selecting suitable ORFs

    @type  hmmsearch_num_hits: integer
    @param hmmsearch_num_hits: (small) positive number of best HMM hits to take into account

    @type  verbose: Boolean
    @param verbose: report an extended debugging-report to STDOUT (True) or be quiet (False)
    """

    # try to limit searchspace by prev and next CBG
    prevNode, nextNode = None, None
    prevMin,  nextMax  = None, None
    maskcoords = []
    # check if (informant) organism is in the prev CBG AND if this CBG
    # has an OMSR -> not per se the case!
    if prev and org in prev.organism_set() and prev.has_overall_minimal_spanning_range():
        prevNode = prev.node_by_organism(org)
        omsr = prev.overall_minimal_spanning_range(organism=org)
        prevMin = (max(omsr)+1)*3
        maskcoords.append( ( 0, max(omsr) ) )
    # check if (informant) organism is in the next CBG AND if this CBG
    # has an OMSR -> not per se the case!
    if next and org in next.organism_set() and next.has_overall_minimal_spanning_range():
        nextNode = next.node_by_organism(org)
        omsr = next.overall_minimal_spanning_range(organism=org)
        nextMax = min(omsr)*3
        aaseqlen = len(inputdict[org]['genomeseq'])/3
        maskcoords.append( ( min(omsr), aaseqlen ) )
    if not prev and next and nextMax:
        prevMin = nextMax - max_intron_nt_length
    if not next and prev and prevMin:
        nextMax = prevMin + max_intron_nt_length 


    # make unique filename for hmm database and hmmbuild file
    fname_base = get_random_string_tag()
    fname_fasta_db = "hmm_%s_database.fa" % fname_base
    fname_fasta_alignment = "%s.fa" % fname_base

    # get elegiable sets of orfs from prev and next
    if not orf_must_have_start:
        elegiable_orfs = inputdict[org]['orfs'].get_elegiable_orfs(
                min_orf_end = prevMin,
                max_orf_start = nextMax
                )
    else:
        # ORFs *must* have starts => searching for a TSS exon/CBG
        elegiable_orfs = inputdict[org]['orfs'].get_elegiable_orfs(
                min_orf_end = prevMin,
                max_orf_start = nextMax,
                has_starts=True
                )

    ####################################################################
    if verbose:
        if len(elegiable_orfs) > 10:
            print "hmm-elegibable orfs:", org, len(elegiable_orfs), "/",
            print len(inputdict[org]['orfs'].orfs), "prevMin:", prevMin,
            print "nextMax:", nextMax
        else:
            print "hmm-elegibable orfs:", org, [ orf.id for orf in elegiable_orfs ], "/",
            print len(inputdict[org]['orfs'].orfs) 
    ####################################################################

    # check elegiable orfs; can be zero in case of a very tiny region to check
    if not elegiable_orfs:
        if return_hmm_query_coord_dict:
            # return the hmm hits & the basename AND the hmm profile coordinates!
            return [], fname_base, {}
        else:
            # return the hmm hits & the basename of the hmm files for later cleanup
            return [], fname_base

    min_sprdif_left_aa_size = 20
    if inwpcbg.has_left_spanningrange_difference(sprdif_min_aa_length=min_sprdif_left_aa_size):
        sprdL = inwpcbg.left_spanningrange_difference(sprdif_min_aa_length=min_sprdif_left_aa_size)
        fastaseqs, coords = inwpcbg.get_left_sprdif_proteinsequences_and_coords(sprdif_min_aa_length=min_sprdif_left_aa_size)
        # rename key to node string representations
        nodes = fastaseqs.keys()
        for node in nodes:
            #str_node = "%s_%s" % node
            str_node = "%s" % node[0]   # organism identifier expected, not Node

            fastaseqs[str_node]   = fastaseqs[node]
            coords[str_node] = coords[node]
            del( fastaseqs[node] )
            del( coords[node] )

        # remove to short ones
        seqdata = []
        for k,seq in fastaseqs.iteritems():
            seqdata.append( ( len(seq), k ) )
        seqdata.sort()
        seqdata.reverse()
        while ( float(seqdata[-1][0]) / seqdata[0][0] ) < 0.70:
            ll,node = seqdata.pop(-1)
            del(fastaseqs[node])
            del(coords[node])
            print "DELETED::", node, ll, float(ll) / seqdata[0][0]

    else:
        if return_hmm_query_coord_dict:
            # return the hmm hits & the basename AND the hmm profile coordinates!
            return [], fname_base, {}
        else:
            # return the hmm hits & the basename of the hmm files for later cleanup
            return [], fname_base

    if not fastaseqs or len(fastaseqs)==1:
        # no fasta multiple alignment left -> return no results!
        if return_hmm_query_coord_dict:
            # return the hmm hits & the basename AND the hmm profile coordinates!
            return [], fname_base, {}
        else:
            # return the hmm hits & the basename of the hmm files for later cleanup
            return [], fname_base


    (alignedseqs,alignment) = clustalw( seqs= fastaseqs )
    #for k,algseq in alignedseqs.iteritems():
    #    print algseq,k, min(coords[k]),max(coords[k])
    #print alignment


    #alignedseqs,alignment,coords = strip_alignment_for_exterior_gaps(
    #    deepcopy(alignedseqs),deepcopy(alignment),deepcopy(coords) )
    #
    #print "AFTER STRIPPING:"
    #for k,algseq in alignedseqs.iteritems():
    #    print algseq,k, min(coords[k]),max(coords[k])
    #print alignment


    # write multiple alignment input file
    writeMultiFasta(alignedseqs,fname_fasta_alignment)

    # write masked orfs to fasta database
    db_fasta = inputdict[org]['orfs'].tomaskedfasta(coords=maskcoords,orflist=elegiable_orfs,header_prefix=org) 
    if orf_must_have_start:
        if len(db_fasta.strip()) == 0:
            # no UNmasked suitable ORFs remaining!
            # This is recognized lateron in this function 
            pass
        else:
            # mask out all AAs before the first start
            lines = db_fasta.split("\n")
            for linenr in range(0,len(lines)):
                line = lines[linenr]
                if line[0] != ">":
                    mpos = line.find("M")
                    if mpos > 0:
                        line = "X"*mpos+line[mpos:]
                    lines[linenr] = line
            # recreate db_fasta string
            db_fasta = "\n".join(lines)

    # write masked orfs to fasta database
    fh = open(fname_fasta_db,'w')
    fh.write( db_fasta )
    fh.close()

    # make shure that there where orfs written to file;
    # in case very little orfs are selected and all are masked -> no files!
    seqs_in_db = parseFasta(open(fname_fasta_db).readlines())
    if not seqs_in_db:
        # return no results!
        if return_hmm_query_coord_dict:
            # return the hmm hits & the basename AND the hmm profile coordinates!
            return [], fname_base, coords
        else:
            # return the hmm hits & the basename of the hmm files for later cleanup
            return [], fname_base

    # make hmmbuild file of the multiplealignment
    fname_hmmbuild = hmmbuild_protein( fname_fasta_alignment )

    if fname_hmmbuild:
        # run hmmsearch on this fasta database with hmmbuild file
        results = hmmsearch_protein( fname_hmmbuild, fname_fasta_db, params= {'E':150,'A': hmmsearch_num_hits} )
    else:
        # no fname_hmmbuild file => creation failed for some reason
        results = []

    if return_hmm_query_coord_dict:
        # return the hmm hits & the basename AND the hmm profile coordinates!
        return results, fname_base, coords
    else:
        # return the hmm hits & the basename of the hmm files for later cleanup
        return results, fname_base

# end of function hmmsearch_results_for_inwpcbg

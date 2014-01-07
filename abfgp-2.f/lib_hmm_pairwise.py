"""
"""

# Module metadata
__authors__ = "Ate van der Burgt"
__license__ = "MIT"

# import hmmparser items
from parsers.hmm import hmmbuild_protein, hmmsearch_protein

# Imports
from os import remove as osRemove
from sets import Set
from copy import deepcopy

# Pacb Imports
from pacb import PacbP, PacbPORF
from pacb.conversion import pacbporf2pacbp, pacbp2pacbporf
from pacb.ordering import order_pacbp_list
from pacb.overlap import correct_overlap_for_sbjct, correct_overlap_for_query
from pacb.splitting import split_pacb_on_gaps
from pacb import swap_query_and_sbjct

from graphAbgp.graph_pacbpcollection import pacbporf2PCG

# ABGP items
from pythonlibs.uniqueness import get_random_string_tag
from lib_fasta import parseFasta, writeMultiFasta
from lib_clustalw import (
    clustalw,
    strip_alignment_for_exterior_gaps,
    strip_overall_nonaligned_residues,
    strip_poorly_supported_tails,
    )
from lib_crossblast import _order_list_by_attribute
from lib_stopwatch import StopWatch

def _remove_short_sprdif_contributors(coords,lengthratio=0.40,verbose=False):
    """ Helper function for _create_hmm_profile """
    # remove short contributors to SPRDIF
    lengths = []
    for node,vlist in coords.iteritems():
        lengths.append( ( max(vlist)-min(vlist), node ) )
    lengths.sort()
    lengths.reverse()
    while True:
        if lengths[0][0] == 0:
            coords = {}
            break
        if float(lengths[-1][0]) / lengths[0][0] < lengthratio:
            ll,node = lengths.pop(-1)
            del(coords[node])
            ####################################################################
            if verbose:
                print "SPRDIF DELETED::",
                print node, ll, float(ll) / lengths[0][0]
            ####################################################################
        else:
            break
    # return coord list
    return coords

# end of function _remove_short_sprdif_contributors


def _rename_dict_keys_to_strings(*args):
    """ Helper function for _create_hmm_profile """
    # rename key to node string representations
    nodes = args[0].keys()
    for node in nodes:
        str_node = "%s" % node[0]   # organism identifier expected, not Node
        for pos in range(0,len(args)):
            args[pos][str_node] = args[pos][node]
            del( args[pos][node] )

    # return the (list of) dict arguments
    return args

# _rename_dict_keys_to_strings

def _create_hmm_profile(cbg,area="OMSR",prevcbg=None,nextcbg=None,
    strip_nonaligned_residues=False,
    verbose=False,**kwargs):
    """
    """
    # area must be one of 
    # OMSR MINSR MAXSR
    # LEFTSPRDIF RIGTHSPRDIF
    # OMSRANDLEFTSPRDIF OMSRANDRIGTHSPRDIF
    # RIGTHORFEND

    # update to default value
    if not kwargs.has_key('sprdif_min_aa_length'):
        kwargs['sprdif_min_aa_length'] = 20

    if area == "OMSR":
        if cbg.has_overall_minimal_spanning_range():
            coords = cbg.overall_minimal_spanning_range()
        else:
            return None, {}
    elif area == "MINSR":
        if cbg.has_minimal_spanning_range():
            coords = cbg.minimal_spanning_range()
        else:
            return None, {}
    elif area == "MAXSR":
        if cbg.has_maximal_spanning_range():
            coords = cbg.maximal_spanning_range()
        else:
            return None, {}
    elif area == "LEFTSPRDIF":
        if cbg.has_left_spanningrange_difference(**kwargs):
            coords = cbg.left_spanningrange_difference(**kwargs)
        else:
            return None, {}
    elif area == "RIGTHSPRDIF":
        if cbg.has_rigth_spanningrange_difference(**kwargs):
            coords = cbg.rigth_spanningrange_difference(**kwargs)
        else:
            return None, {}
    elif area == "OMSRANDLEFTSPRDIF":
        kwargs['sprdif_min_aa_length'] = 20
        if not cbg.has_overall_minimal_spanning_range() or\
        not cbg.has_left_spanningrange_difference(**kwargs):
            return None, {}
        # if here, start preparing coords
        coords = cbg.left_spanningrange_difference(**kwargs)
        # remove short contributors to left SPRDIF
        coords = _remove_short_sprdif_contributors(coords,verbose=verbose)
        # increase coord range by OMSR area
        omsr = cbg.overall_minimal_spanning_range()
        for node,coordrange in coords.iteritems():
            coords[node] = Set( range( min(coordrange), max(omsr[node])+1 ) )
    elif area == "OMSRANDRIGTHSPRDIF":
        kwargs['sprdif_min_aa_length'] = 20
        if not cbg.has_overall_minimal_spanning_range() or\
        not cbg.has_rigth_spanningrange_difference(**kwargs):
            return None, {}
        # if here, start preparing coords
        coords = cbg.rigth_spanningrange_difference(**kwargs)
        # remove short contributors to left SPRDIF
        coords = _remove_short_sprdif_contributors(coords,verbose=verbose)
        # increase coord range by OMSR area
        omsr = cbg.overall_minimal_spanning_range()
        for node,coordrange in coords.iteritems():
            coords[node] = Set( range( min(omsr[node]), max(coordrange)+1 ) )
    elif area == "RIGTHORFEND":
        # area in between MAXSR and orfend
        if not cbg.has_maximal_spanning_range(): return None, {}
        # get coords & obtain Orf ends
        coords = cbg.maximal_spanning_range()
        nodes = coords.keys()
        for node in nodes:
            organism = cbg.organism_by_node(node)
            theorf = cbg.get_orfs_of_graph(organism=organism)[0]
            coords[node] = range(max(coords[node])+1,theorf.protein_endPY)
            # remove zero-length ranges
            if len(coords[node]) == 0: del(coords[node])
    else:
        raise "WHAT ELSE!?"

    ############################################################################
    if verbose: print area, sum([(max(v)-min(v)) for k,v in coords.iteritems()]),len(coords)
    ############################################################################

    # decrease coord range by prevcbg if applicable
    if area in ["MAXSR","LEFTSPRDIF","OMSRANDLEFTSPRDIF"] and prevcbg:
        omsr = prevcbg.overall_minimal_spanning_range()
        for org in cbg.organism_set().intersection( prevcbg.organism_set() ):
            # omsr/coords have Node keys -> translate to Organism keys
            nodeCbg  = cbg.get_organism_nodes(org)[0]
            nodePrev = prevcbg.get_organism_nodes(org)[0]
            # check if node not deleted earlier in coords dict
            if not coords.has_key(nodeCbg): continue
            if not omsr.has_key(nodePrev): continue
            sta = max( [ max(omsr[nodePrev])+1, min(coords[nodeCbg]) ] )
            end = max(coords[nodeCbg])+1
            coords[nodeCbg] = Set(range(sta,end))
            if not coords[nodeCbg]: del( coords[nodeCbg] )

    # decrease coord range by nextcbg if applicable
    if area in ["MAXSR","RIGTHSPRDIF","OMSRANDRIGTHSPRDIF"] and nextcbg:
        omsr = nextcbg.overall_minimal_spanning_range()
        for org in cbg.organism_set().intersection( nextcbg.organism_set() ):
            # omsr/coords have Node keys -> translate to Organism keys
            nodeCbg  = cbg.get_organism_nodes(org)[0]
            nodeNext = nextcbg.get_organism_nodes(org)[0]
            # check if node not deleted earlier in coords dict
            if not coords.has_key(nodeCbg): continue
            if not omsr.has_key(nodeNext): continue
            sta = min(coords[nodeCbg])
            end = min( [ min(omsr[nodeNext]), max(coords[nodeCbg])+1 ] )
            coords[nodeCbg] = Set(range(sta,end))
            if not coords[nodeCbg]: del( coords[nodeCbg] )

    # check if coords still present
    if not coords: return None, {}

    ############################################################################
    if verbose: print area, sum([(max(v)-min(v)) for k,v in coords.iteritems()]),len(coords)
    ############################################################################

    # do/redo _remove_short_sprdif_contributors id required
    if area in ["MAXSR","LEFTSPRDIF","RIGTHSPRDIF",
    "OMSRANDLEFTSPRDIF","OMSRANDRIGTHSPRDIF","RIGTHORFEND"]:
        coords = _remove_short_sprdif_contributors(coords)

    ############################################################################
    if verbose: print area, sum([(max(v)-min(v)) for k,v in coords.iteritems()]),len(coords)
    ############################################################################

    # check if at least 2 sequences/nodes are remaining
    if len(coords) <= 1: return None, {}

    # check sprdif_min_aa_length if applicable
    if area in ["RIGTHSPRDIF","LEFTSPRDIF","OMSRANDRIGTHSPRDIF",
    "OMSRANDLEFTSPRDIF"]:
        maxlength = max([ len(vlist) for vlist in coords.values() ])
        if maxlength < kwargs['sprdif_min_aa_length']:
            return None, {}

    # if here, obtain sequences and build HMM search profile

    # get fasta sequences and 
    fastaseqs = cbg._get_sequences_by_coords(coords)

    # rewrite dict (node) keys to string keys
    fastaseqs, coords = _rename_dict_keys_to_strings(fastaseqs, coords)

    # remove empty sequence strings from fastaseqs dict
    empty_seq_keys = []
    for k,seq in fastaseqs.iteritems():
        if seq == "" or len(seq) == 1:
            empty_seq_keys.append(k)
    for k in empty_seq_keys:
        del(coords[k])
        del(fastaseqs[k])

    # check (again) if at least 2 sequences/nodes are remaining
    if len(coords) <= 1: return None, {}

    # rewrite coords to (min,max) tuple
    coords = dict([ (key,[min(vlist),max(vlist)+1]) for key,vlist in coords.iteritems() ])

    # perform clustalw multiple alignment
    (alignedseqs,alignment) = clustalw( seqs= fastaseqs )


    # strip exterior gaps in case of OMSR/MINSR area
    if area in ["OMSR","MINSR"]:
        alignedseqs,alignment,coords = strip_alignment_for_exterior_gaps(
                deepcopy(alignedseqs),deepcopy(alignment),deepcopy(coords) )


    # strip poorly conserved residues in case of RIGTHORFEND
    if area in ["RIGTHORFEND"]:
        alignedseqs,alignment,coords = strip_poorly_supported_tails(
            deepcopy(alignedseqs),deepcopy(alignment),deepcopy(coords),0.20 )


    # strip_overall_nonaligned_residues if requested for: THIS IS VERY RIGID!
    if strip_nonaligned_residues:
        alignedseqs,alignment,coords = strip_overall_nonaligned_residues(
                deepcopy(alignedseqs),deepcopy(alignment),deepcopy(coords) )
        # check if alignment was completely consumed or not
        if not alignment or len(alignment) <= 1:
            return None, {}


    ############################################################################
    if verbose:
        print "## HMM clustalw input profile:",prevcbg!=None,area,nextcbg!=None
        for node,algseq in alignedseqs.iteritems():
            print algseq, node, coords[node]
        print alignment
    ############################################################################

    # make unique filename for hmm profile file
    fname_hmm_profile = "hmmbuild_profile_%s.hmmprof" % get_random_string_tag()

    # write multiple alignment input file
    writeMultiFasta(alignedseqs,fname_hmm_profile)

    # make hmmbuild file of the multiplealignment
    fname_hmmbuild_file = hmmbuild_protein( fname_hmm_profile )

    # remove hmm profile multiple alignment file
    osRemove(fname_hmm_profile)

    # return HMM serach profile filename
    return fname_hmmbuild_file, coords

# end of function _create_hmm_profile



def _create_hmm_db(organism,inputdict,cbg,prev,next,
    orf_must_have_start=False,max_intron_nt_length=200,
    verbose=False):
    """
    Create fasta ORF database for a organism in a CBG and its viscinity

    @type  organism: * (presumably string)
    @param organism: Organism identifier recognizable in <input data structure>

    @type  inputdict: dict 
    @param inputdict: <input data structure> 

    @type  cbg: CodingBlockGraph or related object
    @param cbg: CodingBlockGraph upstream/5p of the cbg that must be completed

    @type  prev: CodingBlockGraph or related object (or None)
    @param prev: CodingBlockGraph upstream/5p of cbg that must be completed

    @type  next: CodingBlockGraph or related object (or None)
    @param next: CodingBlockGraph downstream/3p of cbg that must be completed

    @attention: `prev` and `next` CodingBlockGraphs reduce the search space of
                ORFs to scan with the HMM profile. This Speeds up and
                improves the quality of results.

    @type  orf_must_have_start: Boolean
    @param orf_must_have_start: only allow ORFs with methionines as sbjct ORFs

    @type  max_intron_nt_length: integer
    @param max_intron_nt_length: positive maximum intron length to take
                                 into acount when selecting suitable ORFs

    @type  verbose: Boolean
    @param verbose: report debugging-report on STDOUT (True) or be quiet (False)
    """

    # fullpath filename of result hmm multi fasta database
    fname_hmm_db_mfa = None
    if not cbg: return fname_hmm_db_mfa

    # (1) try to limit searchspace by prev and next CBG
    prevNode, nextNode = None, None
    prevMin,  nextMax  = None, None
    maskcoords = []

    # (1a) check if (informant) organism is in the prev CBG AND if this CBG
    # has an OMSR -> not per se the case!
    if prev and organism in prev.organism_set() and\
    prev.has_overall_minimal_spanning_range():
        prevNode = prev.node_by_organism(organism)
        try:
            omsr = prev.overall_minimal_spanning_range(organism=organism)
            prevMin = (max(omsr)+1)*3
            maskcoords.append( ( 0, max(omsr) ) )
        except KeyError:
            # hmmm.... block has an OMSR, but not for this organism!??!!?
            pass


    # (1b) check if (informant) organism is in the next CBG AND if this CBG
    # has an OMSR -> not per se the case!
    if next and organism in next.organism_set() and\
    next.has_overall_minimal_spanning_range():
        nextNode = next.node_by_organism(organism)
        try:
            omsr = next.overall_minimal_spanning_range(organism=organism)
            nextMax = min(omsr)*3
            aaseqlen = len(inputdict[organism]['genomeseq'])/3
            maskcoords.append( ( min(omsr), aaseqlen ) )
        except KeyError:
            # hmmm.... block has an OMSR, but not for this organism!??!!?
            pass

    # (1c) limit search space if only prev or next was specified
    if not prev and next and nextMax:
        prevMin = nextMax - max_intron_nt_length
    if not next and prev and prevMin:
        nextMax = prevMin + max_intron_nt_length 

    # (2a) get elegiable sets of orfs from prev and next
    if not orf_must_have_start:
        elegiable_orfs = inputdict[organism]['orfs'].get_elegiable_orfs(
                min_orf_end = prevMin, max_orf_start = nextMax
                )
    else:
        # ORFs *must* have starts => searching for a TSS exon/CBG
        elegiable_orfs = inputdict[organism]['orfs'].get_elegiable_orfs(
                min_orf_end = prevMin, max_orf_start = nextMax,
                has_starts=True
                )

    # (2b) check orf count; can be zero in case of a very tiny region to check
    if not elegiable_orfs: return fname_hmm_db_mfa

    # (3) write masked orfs to fasta database multi line string
    db_fasta = inputdict[organism]['orfs'].tomaskedfasta(
            coords=maskcoords,
            orflist=elegiable_orfs,
            header_prefix=organism) 
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

    ############################################################################
    if verbose:
        if len(elegiable_orfs) > 10:
            orfidlist = len(elegiable_orfs)
        else:
            orfidlist = [ orf.id for orf in elegiable_orfs ]
        print "hmm-elegibable orfs:", organism, orfidlist, "/",
        print len(inputdict[organism]['orfs'].orfs), "prevMin:", prevMin,
        if prev:
            print prev.has_overall_minimal_spanning_range(),
        else:
            print None,
        print "nextMax:", nextMax,
        if next:
            print next.has_overall_minimal_spanning_range()
        else:
            print None
    ############################################################################

    # (4) make unique filename for hmm database file
    fname_base = get_random_string_tag()
    fname_hmm_db_mfa = "hmm_database_%s_%s.fa" % (fname_base,organism)

    # (5) write masked orfs to fasta database
    fh = open(fname_hmm_db_mfa,'w')
    fh.write( db_fasta )
    fh.close()

    # (6) make shure that there where orfs written to file;
    # in case very little orfs are selected and all are masked -> no files!
    seqs_in_db = parseFasta(open(fname_hmm_db_mfa).readlines())
    if not seqs_in_db:
        # delete this (empty) file
        osRemove( fname_hmm_db_mfa )
        return None

    # (7) return hmm search database filename
    return fname_hmm_db_mfa

# end of function _create_hmm_db



def is_hmmpacbporf_conflicting_with_pacbporflist(hmmpacbporf,pacbporflist):
    """ """
    IS_HMMPACBP_CONFLICTING = False
    for pacbporf in pacbporflist:
        # check if positioned compatibly
        if not pacbporf.is_postioned_compatibly(hmmpacbporf):
            overlap = False # init printing variable
            IS_HMMPACBP_CONFLICTING = True
            break
        # check if not overlapping
        overlap = pacbporf.overlap(hmmpacbporf)
        if overlap == 0.0:
            pass
        elif overlap <= 0.25:
            # correct for slightly overlapping PacbPORFS
            # Lazy... not willing to check orientation of
            # PacbPs here; let the overlap function handle it
            thispacbp = pacbporf2pacbp(pacbporf)
            hmmpacbp  = pacbporf2pacbp(hmmpacbporf)

            _prev,_next = order_pacbp_list([thispacbp,hmmpacbp])
            _prev, _next, status1 = correct_overlap_for_sbjct(
                        _prev, _next , verbose=False )
            _prev, _next, status2 = correct_overlap_for_query(
                        _prev, _next , verbose=False)

            if hmmpacbp.length == 0:
                IS_HMMPACBP_CONFLICTING = True
                break
            if thispacbp.length == 0:
                print "FatalWarning: HMM overlap caused PacbPORF to dissapear"
                IS_HMMPACBP_CONFLICTING = True
                break

            # Okay! Convert back to the pacbporf & the hmmpacbporf
            hmmpacbporf = pacbp2pacbporf(hmmpacbp,
                    hmmpacbporf.orfQ,hmmpacbporf.orfS)

        else:
            IS_HMMPACBP_CONFLICTING = True
            break

    # return binary outcome of overlap conflict
    return IS_HMMPACBP_CONFLICTING

# end of function is_hmmpacbporf_conflicting_with_pacbporflist



def hmmresults2splittedpacbps(results,hmmcoords,target,informant,inwpCBG,input,gapsize=2,min_bitscore=0):
    """ """
    hmm_pacbporf_list = []
    for hmmhit in results:
        ( sbjct_header, sbjct_start, sbjct_end, query_start, query_end,
          query, match, sbjct, score, expect ) = hmmhit
        if score < min_bitscore: continue
        _org,orfid = hmmhit[0].split('_orf_')
        orfSbjct = input[informant]['orfs'].get_orf_by_id( int(orfid) )
        orfQuery = inwpCBG.get_orfs_of_graph(target)[0]
        querycoords = ( min(hmmcoords[target]),max(hmmcoords[target]) )
        key_data, hmmpacbporf = hmmhit2pacbp(
                orfQuery,target,querycoords,
                orfSbjct,informant,hmmhit)

        # check if hmmpacbporf creation was succesfull
        if not hmmpacbporf: continue

        # if here, unextend and split on gapsize
        (pacbpkey,qNode,sNode) = key_data 
        hmmpacbporf.unextend_pacbporf()
        splittedhmmpacbporfs, splittedstatus =\
            split_pacb_on_gaps(hmmpacbporf,gapsize=gapsize)

        # loop over the splitted ones and store high(er) scoring fractions
        for splittedhmmpf in splittedhmmpacbporfs:
            # added code to strip unmatched ends. Should not
            # be neccesarily anymore, but just to be certain
            # no leading/trailing gaps are there
#            if '-' in [ splittedhmmpf.query[0], splittedhmmpf.sbjct[0],
#            splittedhmmpf.query[-1], splittedhmmpf.sbjct[-1] ]:
#                hmmpacbp = pacbporf2pacbp(splittedhmmpf)
#                hmmpacbp.strip_unmatched_ends()
#                if not hmmpacbp: continue
#                if len(hmmpacbp) <= 1: continue
#                # if here, make again a pacbporf of the pacbp
#                splittedhmmpf = pacbp2pacbporf(hmmpacbp,splittedhmmpf.orfQ,splittedhmmpf.orfS)

            if splittedhmmpf.bitscore < min_bitscore: continue
            # check if query/sbjct must be swapped
            queryNode = (target,orfQuery.id)
            sbjctNode = (informant,orfSbjct.id)
            if qNode == queryNode:
                pass
            elif qNode == sbjctNode:
                # swap query and sbjct!
                splittedhmmpf = swap_query_and_sbjct(splittedhmmpf)
            else:
                # whaaaat else !?
                raise "UNEXPECTED EVENT"

            # append to hmm_pacbporf_list
            hmm_pacbporf_list.append( splittedhmmpf )

    # return bitscore ordered list of hmmpacbporfs
    return _order_list_by_attribute(hmm_pacbporf_list,
            order_by='bitscore',reversed=True)

# end of function hmmresults2splittedpacbps


def hmmpacbporf2PCG(hmmpacbporf,target,informant,PCG,source=''):
    """ """

    hmmpacbporf.source = 'HMM'
    hmmpacbporf._gff['fsource'] = 'HMM'
    pacbporf2PCG(hmmpacbporf,target,informant,PCG,source=source)

    # MOVED TO graphAbgp.graph_pacbpcollection
    #queryNode = (target,hmmpacbporf.orfQ.id)
    #sbjctNode = (informant,hmmpacbporf.orfS.id)
    #hmmpacbporf.extend_pacbporf_after_stops()
    #hmmpacbporf.source = 'HMM'
    #hmmpacbporf._gff['fsource'] = 'HMM'
    #pacbpkey = hmmpacbporf.construct_unique_key(queryNode,sbjctNode)
    #PCG.add_node(sbjctNode)
    #PCG.add_edge(queryNode,sbjctNode,wt=pacbpkey[0])
    #PCG.pacbps[(pacbpkey,queryNode,sbjctNode)] = hmmpacbporf

# end of function hmmpacbporf2PCG


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

    ############################################################################
    if verbose:
        print "hmmq:", queryseq, queryNode, query_aa_start, query_aa_end,
        print "len:", query_aa_end-query_aa_start, len(queryseq)
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

    # NEW NEW code in december 2010. Since inwpCBGs are implemented, HMM
    # profiles are build from clustalw alignments which have loosely aligned
    # tails (SPRDIF sequences). Problem with HMM is, that in the result file
    # no information is written on where in teh constructed HMM this hit
    # starts. This **sucks** because special care was taken in ABFGP code to
    # make shure the exact aa-coordinates of the applied sequences to ClustalW
    # are known. Hmmbuild here nullifies this effort by not giving start
    # coordinates. Therefore, we have to check the exact start position
    # of the HMM match on the queryorf.
    if queryseq.replace("-","") != queryorf.getaas(query_aa_start,query_aa_end):
        # obtain (search) query sequence, replace gaps by X symbol
        searchqueryseq = queryseq.upper().replace("-","X")
        # count length of the query sequence; here IGNORE THE GAPS!!
        seqlen = len(queryseq.upper().replace("-",""))

        # make fasta sequence dictionary
        seqdict = {
            'query_hmm': searchqueryseq,
            'query_orf': queryorf.protein_sequence,
            }

        # make coords dictionary for remapping
        coords = {
            'query_hmm':[0,seqlen],
            'query_orf':[queryorf.protein_startPY,queryorf.protein_endPY],
            }

        # perform clustalw multiple alignment
        (alignedseqs,alignment) = clustalw( seqs= seqdict )
        # strip exterior gaps
        alignedseqs,alignment,coords = strip_alignment_for_exterior_gaps(
            deepcopy(alignedseqs),deepcopy(alignment),deepcopy(coords) )

        if alignedseqs['query_hmm'].count("-") > 0:
            # in (very) exceptional cases, gaps can be introduced in the
            # clustalw alignment in the HMM seq. This normally does not
            # occur! Fix this here by placing gaps in sbjctseq too.
            sbjctseq_as_list = list(sbjctseq)
            for pos in range(0,len(alignedseqs['query_hmm'])):
                if alignedseqs['query_hmm'][pos] == "-":
                    sbjctseq_as_list.insert(pos,"-")
                if alignedseqs['query_hmm'].find("-",pos) == -1:
                    break
            sbjctseq = "".join(sbjctseq_as_list)

        ########################################################################
        if verbose:
            print "\t", "FALSE::", sbjctseq, "[ WITH GAPS,SBJCT ]" 
            print "\t", "FALSE::", queryseq, "[ WITH GAPS ]" 
            for k,algseq in alignedseqs.iteritems():
                print "\t", "FALSE::", algseq, k, coords[k], len(algseq)
            print "\t", "FALSE::", sbjctseq, "SBJCT", len(sbjctseq)
            print "\t", "FALSE::", alignment, "ALMNT", len(alignment)
            print "\t", "SOLVED:", len(alignedseqs['query_orf']) == len(sbjctseq)
        ########################################################################
    
        # update query sequence & coordinates
        if len(alignedseqs['query_orf']) == len(sbjctseq):
            queryseq       = alignedseqs['query_orf']
            query_aa_start = coords['query_orf'][0]
            query_aa_end   = coords['query_orf'][1]
        else:
            # still not identical lengths. ClustalW recovery of HMM hit
            # failed miserably. For now: omit
            # TODO: resolve this case!!
            # example: --filewithloci examples/bilal/CFU_830450.bothss.csv
            # ## HMM clustalw input profile: False MAXSR True
            # FPKGCESGKFINWKTFKANGVNLGAWLAKEKTHDPVW foxga [561, 598]
            # FQRACR--KFID-ETLSAHAL---EWESKEIVPPEVW CFU [357, 388]
            # hmmhit2pacbp CREATING pacbps for organism/orf: (NP1064101[anid],1)
            # hmmhit2pacbp Q 'FQKACRSGKFIDWKTLKANALNLGEWLAKEKVHD'
            # hmmhit2pacbp m '+ ka +   F  W   k  + nLG Wl  E   d'
            # hmmhit2pacbp S 'YTKAFQ--PF-SWSSAKVRGANLGGWLVQEASID'
            # hmmQ: FQKACRSGKFIDWKTLKANALNLGEWLAKEKVHD 1 34 gaps: 0 34
            # hmmM: + ka +   F  W   k  + nLG Wl  E   d
            # hmmS: YTKAFQ--PF-SWSSAKVRGANLGGWLVQEASID ('NP1064101[anid]', 1) 33 64 len: 31 34
            # hmmq: FQKACRSGKFIDWKTLKANALNLGEWLAKEKVHD ('CFU', 91) 357 391 len: 34 34
            #         FALSE:: YTKAFQ---------PF-SWSS-----------------AKVR----------GANLGG--W-LVQEASID [ WITH GAPS,SBJCT ]
            #         FALSE:: FQKACRSGKFIDWKTLKANALNLGEWLAKEKVHD [ WITH GAPS ]
            #         FALSE:: FQKACR-------SGKFIDWKT-----------------LKAN----------ALNLGE--W-LAKEKVH query_hmm [0, 33] 70
            #         FALSE:: FQRACRKFIDETLSAHALEWESKEIVPPEVWQRFAEANMLIPNLAALASRMVGEIGIGNAFWRLSVQGLR query_orf [357, 427] 70
            #         FALSE:: YTKAFQ---------PF-SWSS-----------------AKVR----------GANLGG--W-LVQEASID SBJCT 71
            #         FALSE:: **:***       *.: ::*::                 * .*           :.:*:  * *: : :: ALMNT 70
            #         SOLVED: False
            # Pacbp creation failed!
            return False, None

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
            print query_aa_end - query_aa_start, "len:", len(queryseq)
            print "S:", sbjctaastart, sbjctaaend,
            print sbjctaaend - sbjctaastart, "len:",len(sbjctseq)
        ################################################################
        pacbpinput = (queryseq,sbjctseq,query_aa_start,sbjctaastart)
        pacbp      = PacbP(input=pacbpinput)
        # remove consistent internal gaps caused hy HMM profile search
        pacbp.strip_consistent_internal_gaps()
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
        # Pacbp creation failed!
        return False, None

# end of function hmmhit2pacbp

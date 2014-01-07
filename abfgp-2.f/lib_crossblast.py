"""
Functions in use for iterative crossblastp in ABFGP

# Main functions
  def createblastdbs(input,GSG,OPTIONS,dbfraction=None,organism=None,acceptorfids=[],rejectorfids=[]):
  def iterativecrossblastp(input,crossdata,GTgraph=None,CBgraph=None,GSgraph=None, ... ):

# Helper functions
  def _blastorfset2blastdb(geneQ,geneS,blastdbfname,input,crossdata,GSgraph, ... ):
  def _get_orfids_from_crossdata(geneid,crossdata):
  def _order_list_by_attribute(inputlist,order_by='',reversed=False):


"""
# Abfgp imports
from lib_blastp import formatdb, blastall_seq2db 
from pythonlibs.uniqueness import get_random_string_tag
from abgp_etc import _blastdb_cleanup
import pacb

# Python imports
from sets import Set
from os.path import join as osPathJoin

# import linearization settings
from settings.pacbp import LINEARIZATION_STARTWITHBEST, LINEARIZATION_ACCEPTANCE_MARGIN, LINEARIZATION_WEIGTHED

def createblastdbs(input,GSG,OPTIONS,dbfraction=None,organism=None,acceptorfids=[],rejectorfids=[]):
    """
    (Re)create blast-db's by masking the areas thar are incorporated in the GSG

    @type  input: dict
    @param input: `input` data structure dictionary

    @type  GSG: GenestructureOfCodingBlockGraphs
    @param GSG: GenestructureOfCodingBlockGraphs instance

    @type  OPTIONS: optparse options instance
    @param OPTIONS: optparse options instance (with attribute 'abinitio')

    @type  dbfraction: string
    @param dbfraction: None, 'all', 'GSGupstream', 'GSGcentral', 'GSGdownstream', 'annotation'

    @type  organism: organism identifier
    @param organism: only recreate blastdb for this organism/gene identifier

    @type  acceptorfids: list with integers
    @param acceptorfids: list of orf ids to accept

    @type  rejectorfids: list with integers
    @param rejectorfids: list of orf ids to reject

    @attention: acceptorfids and rejectorfids are only used when organism is specified!
    """
    seqsindb = {}
    for org in input.keys():
        # if organism is given, do only this one
        if organism and org!=organism: continue
        # acceptorfids anc rejectorfids only valid in combi with `organism`
        if not organism: acceptorfids, rejectorfids = [], [] 
  
        # assign blast database name / multi fasta file and open filehandle
        uniquetag = get_random_string_tag()
        fname = '%s-blastdb-%s.fa' % (uniquetag,org)
        fullpath = osPathJoin(OPTIONS.outdir,fname)
        fh = open(fullpath,'w')
        seqsindb[org] = 0

        # distinct cases possible:
        if len(GSG):
            # there is already a GSG, so this is not the first blast iteration
            # do not apply a shortcut when OPTIONS.abinitio == False
            coords = GSG.omsr2mask(org)
            if dbfraction == 'GSGupstream':
                # take only orfs LEFT of the first CBG in GSG
                max_orf_nt_start = max(GSG[0].overall_minimal_spanning_range(organism=org)) * 3
                orflist = input[org]['orfs'].get_elegiable_orfs(max_orf_start=max_orf_nt_start,
                        acceptorfids=acceptorfids,rejectorfids=rejectorfids)
            elif dbfraction == 'GSGdownstream':
                # take only orfs RIGTH of the last CBG in GSG
                min_orf_nt_end = min(GSG[-1].overall_minimal_spanning_range(organism=org)) * 3
                orflist = input[org]['orfs'].get_elegiable_orfs(min_orf_end=min_orf_nt_end,
                        acceptorfids=acceptorfids,rejectorfids=rejectorfids)
            elif dbfraction == 'GSGcentral':
                # take only orfs in between FIRST and LAST CBG in GSG (can be only one CBG!)
                max_orf_nt_start = max(GSG[-1].overall_minimal_spanning_range(organism=org)) * 3
                min_orf_nt_end   = min(GSG[0].overall_minimal_spanning_range(organism=org)) * 3
                orflist = input[org]['orfs'].get_elegiable_orfs(min_orf_end=min_orf_nt_end,
                        max_orf_start=max_orf_nt_start,
                        acceptorfids=acceptorfids,rejectorfids=rejectorfids)
            else:
                # dbfraction equals 'all' or None -> no limitation, just take all orfs!
                # do only the general limitation on sublists of orfids
                orflist = input[org]['orfs'].get_elegiable_orfs(
                        acceptorfids=acceptorfids,rejectorfids=rejectorfids)

            # create masked fasta of this sequence part only
            newfasta = input[org]['orfs'].tomaskedfasta(coords=coords,orflist=orflist,header_prefix=org)
            # write to file and count accessions in this file -> seqsindb[org]
            fh.write(newfasta)
            seqsindb[org] = newfasta.count(">")

        else:
            # No filled GSG objects -> no a priori knowledge yet
            # When dbfraction=='annotated' and !OPTIONS.abinitio -> take annotated orfs only
            # TODO: dbfraction is not checked/used here -> just OPTIONS.abinitio
            for orf in input[org]['orfs'].orfs:
                # in case not abinitio, make only a db of orfs in teh current annotation!
                if OPTIONS.abinitio == False and orf.id not in input[org]['orfid-genestructure']:
                    continue
                if orf.id in rejectorfids:
                    # ignore Orfs that are listed as to-be-ignored
                    continue
                if acceptorfids and orf.id not in acceptorfids:
                    # ignore Orfs that are not listed as to-be-accepted
                    continue
                # write fasta of orf to file
                fh.write(orf.tofasta(header="%s_orf_%s" % (org,orf.id))+"\n")
                # increase seqsindb[org] counter
                seqsindb[org]+=1

        # close the filehandle
        fh.close()
        # run formatdb
        formatdb(fname=fullpath)
        # set name of blastdb in infodict
        input[org]['blastdb'] = fullpath

    # return the counter of how much orf sequences are stored in the blast database
    return seqsindb

# end of function createblastdbs


def _blastorfset2blastdb(geneQ,geneS,blastdbfname,input,crossdata,GSgraph,
    blastoptions = None,
    elegiable_orfsQ=[],
    elegiable_orfsS_ids=[],
    logging=False):
    """
    """
    hitcnt = 0
    for orfQ in elegiable_orfsQ:
        # check if protein sequence present in Orf object
        # in obscure cases of unigenes, no protein sequence present!
        if not orfQ.protein_sequence: continue

        # make unique node identifier and blast header
        nodeQ = (geneQ,orfQ.id)
        header = "%s_orf_%s" % (geneQ,orfQ.id)

        # do the blastp!
        blastrec = blastall_seq2db(header, orfQ.protein_sequence,
                dbname=blastdbfname,
                extra_blastp_params=blastoptions.extra_blastp_params )

        # check if blast failed (then, blastrec == False)
        if not blastrec: continue
    
        # check if there are any hits/hsps!
        if len(blastrec.alignments) == 0:
            # no hits; continue
            continue
    
        for alignment in blastrec.alignments:
            # get back orfpointerB from the SBJCT and create nodeS
            _parts = alignment.title.split("_")
            geneS = "_".join(_parts[0:-2]).replace('>','')
            _orfpointerS = int(_parts[-1])
            nodeS = (geneS,_orfpointerS)

            # ignore hit if nodeS orfid not occurring in the NON-empty list elegiable_orfsS_ids
            if elegiable_orfsS_ids and _orfpointerS not in elegiable_orfsS_ids:
                continue

            # get the Orf object of this sbjct sequence
            orfS = input[geneS]['orfs'].get_orf_by_id(_orfpointerS)

            # loop over the HSPs
            for hsp in alignment.hsps:

                # If hits are really tiny (happens in case of BLOSUM45 matrix),
                # discard them directly before precious time is lost...
                if len(hsp.query) <= blastoptions.BLASTP_DIRECTLY_IGNORE_TINY_HITS:
                    continue

                # correct to absolute positions
                hsp.query_start = hsp.query_start + orfQ.protein_startPY
                hsp.sbjct_start = hsp.sbjct_start + orfS.protein_startPY
                hsp.query_end = hsp.query_end + orfQ.protein_startPY
                hsp.sbjct_end = hsp.sbjct_end + orfS.protein_startPY

                # VERY exceptional case: HSP starts or ends with a gap
                # I expect this is an error in Blastp .... 
                strip_exterior_gaps(hsp)

                if hsp.query.find(" ") > 0:
                    # VERY exceptional case: erroneously NCBI parsed HSP:
                    # Score 8 (7 bits), expectation 1.7e+01, alignment length 41
                    # Query:    1622 STHTYDAC                                                TRCI----PFVDTGHKHENPTEALLDSTA 1654
                    #                TR      P  +  H +  P+++   S++
                    # Sbjct:     489 TEHIYLHT                                                TRSTWPPKPPTNASHANTKPSKSHHRSSS 525
                    if len(hsp.query.split(" ")[-1]) == len(hsp.match):
                        hsp.query = hsp.query.split(" ")[-1]
                        hsp.sbjct = hsp.sbjct.split(" ")[-1]
                    elif len(hsp.query.split(" ")[0]) == len(hsp.match):
                        hsp.query = hsp.query.split(" ")[0]
                        hsp.sbjct = hsp.sbjct.split(" ")[0]
                    elif len(hsp.query) == len(hsp.match):
                        # spaces in both query/match/sbjct
                        while hsp.query.find(" ") > 0:
                            pos =  hsp.query.find(" ")
                            if hsp.sbjct[pos] == " " and hsp.match[pos] == " ":
                                hsp.query = hsp.query[0:pos] + hsp.query[pos+1:]
                                hsp.match = hsp.match[0:pos] + hsp.match[pos+1:]
                                hsp.sbjct = hsp.sbjct[0:pos] + hsp.sbjct[pos+1:]
                            else:
                                # HPS is not repairable -> quit trying
                                break
                    elif len(hsp.query) == len(hsp.sbjct):
                        while hsp.query.find(" ") > 0:
                            pos =  hsp.query.find(" ")
                            hsp.query = hsp.query[0:pos] + hsp.query[pos+1:]
                            hsp.sbjct = hsp.sbjct[0:pos] + hsp.sbjct[pos+1:]
                        # recreate alignment match string is done upon
                        # creation of PacbP object

                    else:
                        pass

                # VERY exceptional case: HSP starts or ends with a gap
                # I expect this is an error in Blastp .... 
                strip_exterior_gaps(hsp)
                try:
                    pacbp = pacb.PacbP(blastp_hsp=hsp,MATRIX=blastoptions.MATRIX)
                except:
                    # VERY exceptional miscelaneous cases: erroneously NCBI parsed HSP:
                    print hsp
                    print "'%s' X" % hsp.query, len(hsp.query), hsp.query_start, hsp.query_end
                    print "'%s' X" % hsp.match, len(hsp.match)
                    print "'%s' X" % hsp.sbjct, len(hsp.sbjct), hsp.sbjct_start, hsp.sbjct_end
                    pacbp = pacb.PacbP(blastp_hsp=hsp,MATRIX=blastoptions.MATRIX)


                # make pacbp of this hsp
                pacbp = pacb.PacbP(blastp_hsp=hsp,MATRIX=blastoptions.MATRIX)

                # if logging is requested for, print this pacbp to STDOUT
                if logging:
                    print ">>> Q", nodeQ, orfQ.tcode_symbolic(), "S", nodeS, orfS.tcode_symbolic(), pacbp, blastoptions.MATRIX.name, hsp.expect, hsp.bits
                    print ">>>", blastoptions.extra_blastp_params
                    if pacbp.length > 100:
                        print pacbp.query[0:40]+'.'*7+str(pacbp.length-80)+'.'*7+pacbp.query[-40:]
                        print pacbp.match[0:40]+'.'*7+str(pacbp.length-80)+'.'*7+pacbp.match[-40:]
                        print pacbp.sbjct[0:40]+'.'*7+str(pacbp.length-80)+'.'*7+pacbp.sbjct[-40:]
                    else:
                        print pacbp.query
                        print pacbp.match
                        print pacbp.sbjct

                # blastoptions.BLASTP_HSP_MINIMAL_LENGTH represents the minimal
                # length of the aligned part. (To) short pacbp's are abandoned
                if pacbp.length < blastoptions.BLASTP_HSP_MINIMAL_LENGTH:
                    if pacbp.identityscore == float(pacbp.length):
                        # escape for 100% identical tiny pacbps
                        pass
                    elif pacbp.identity + pacbp.similarity == pacbp.length:
                        # escape for 100% similar tiny pacbps
                        pass
                    else:
                        #  pacbp is to small. Discard!
                        if logging: print "to small..."
                        continue

                # check if the pacbp is not conflicting with the currect GSG graph
                # if so, ignore now because it will not yield a proper edge in an (accepted) CBG!
                if GSgraph and len(GSgraph) and GSgraph.is_pacbp_conflicting_with_genestructure(pacbp,orgQ=geneQ,orgS=geneS):
                    ###print "GSGconflict!", nodeQ,nodeS, GSgraph.is_pacbp_conflicting_with_genestructure(pacbp,orgQ=geneQ,orgS=geneS), len(pacbp)
                    continue

                # here we have a potentially accepted pacbp.
                # make a/the unique key of this pacbp
                key = (pacbp.bits, pacbp.length, orfQ.id,_orfpointerS)

                # check for evalue criterion
                if (blastoptions.BLASTP_HSP_MAXIMAL_EXPECT or blastoptions.BLASTP_HSP_MAXIMAL_EXPECT==0.0) and hsp.expect > blastoptions.BLASTP_HSP_MAXIMAL_EXPECT:
                    # pacbp is long enough but has a to high evalue
                    crossdata[(geneQ,geneS)]['lowscoring_pacbs'][key] = pacbp
                    if logging: print "to low bitscore or expect"
                    continue

                # check for bitscore criterion
                if (blastoptions.BLASTP_HSP_MINIMAL_BITS or blastoptions.BLASTP_HSP_MINIMAL_BITS==0) and pacbp.bits < blastoptions.BLASTP_HSP_MINIMAL_BITS:
                    # pacbp is long enough but has a to low bitscore
                    crossdata[(geneQ,geneS)]['lowscoring_pacbs'][key] = pacbp
                    if logging: print "to low bitscore or expect"
                    continue

                # !!Hurray!! an accepted pacbp. Store to crossdata
                # store it to the 'accepted_pacbs' dict of crossdata
                crossdata[(geneQ,geneS)]['accepted_pacbs'][key] = pacbp
                hitcnt+=1
                if logging: print "ACCEPTED"
                # done -> check next orf!

    # return the filled crossdata structure
    return crossdata, hitcnt

# end of function _blastorfset2blastdb


def iterativecrossblastp(input,crossdata,OPTIONS,GTgraph=None,CBgraph=None,GSgraph=None,
    dbfraction = None, blastoptions=None, logging=False ):
    """
    Itrative crossblast function; all is taken care for here

    @type  input: dict
    @param input: `input` dictionary with data, e.g. the list of Orf objects

    @type  crossdata: dict
    @param crossdata: `crossdata` dictionary with the gathered Pacbps in previous iterations

    @type  GTgraph: GeneTreeGraph
    @param GTgraph: GeneTreeGraph instance, in use for ... ??

    @type  GSgraph: GenestructureOfCodingblocksGraph
    @param GSgraph: GenestructureOfCodingblocksGraph instance, with fully maturated
                    CodingBlockGraphs of this genestructure. When applied, any
                    new Pacbp/hsp that conflicts with the genestructure, is rejected. 

    @type  CBgraph: PacbpCollectionGraph
    @param CBgraph: PacbpCollectionGraph instance, with nodes (orfs) and edges
                    (PacbPs/hsps) discovered in a previous iteration

    @type  dbfraction: string
    @param dbfraction: None, 'all', 'GSGupstream', 'GSGcentral', 'GSGdownstream', 'annotated', 'HCP'

    @type  blastoptions: CrossBlastOptions object
    @param blastoptions: CrossBlastOptions object
    """
    # initialize the return counter dict
    crossblasthitcnt = {}

    # get organism edges from GeneTreeGraph (GTG, GTgraph) edge strength
    # If no GTG applied yet, take (unordered) from crossdata
    if GTgraph and GTgraph.node_count():
        ordered_edges = GTgraph.pairwisecrosscombinations_organism(order_by='identity')
    else:
        ordered_edges = crossdata.keys()

    # dict that incrementally gathers the ids of encountered orf nodes
    # these lists will serve as limitations for next pair wise blasts
    encounteredOrfNodeIds = {}

    # loop over the ordered_edges, gather elegiable orf Q list and do the blast!
    for (geneQ,geneS) in ordered_edges:
        # initialize empty counter element
        crossblasthitcnt[(geneQ,geneS)] = {'orfsQ':0,'hits':0}

        # get elegiable_orf_set based on `dbfraction` parameter
        if GSgraph and len(GSgraph)>=1:
            # if this is not the first pairwise cross, the nodes (orfs) encountered
            # in previous pairs are a limitation for this pair wise blast!
            encounteredQorfids = []
            if encounteredOrfNodeIds.has_key(geneQ):
                encounteredQorfids = encounteredOrfNodeIds[geneQ]

            # check if the blast database is already recreated
            if not encounteredOrfNodeIds.has_key(geneS):
                # remove current existing database of this informant 
                _blastdb_cleanup( { geneS: input[geneS] })
                # and make new ones
                sbjctcnts = createblastdbs(input,GSgraph,OPTIONS,dbfraction=dbfraction,organism=geneS)
                if logging: print "new ABINITIO blastdb: ", geneS, sbjctcnts

            # determine action based on `dbfraction` parameter
            if dbfraction == 'GSGupstream':
                # take only orfs LEFT of the first CBG in GSG
                max_orf_nt_start = max(GSgraph[0].overall_minimal_spanning_range(organism=geneQ)) * 3
                orfQlist = input[geneQ]['orfs'].get_elegiable_orfs(
                        acceptorfids=encounteredQorfids,max_orf_start=max_orf_nt_start)
            elif dbfraction == 'GSGdownstream':
                # take only orfs RIGTH of the last CBG in GSG
                min_orf_nt_end = min(GSgraph[-1].overall_minimal_spanning_range(organism=geneQ)) * 3
                orfQlist = input[geneQ]['orfs'].get_elegiable_orfs(
                        acceptorfids=encounteredQorfids,min_orf_end=min_orf_nt_end)
            elif dbfraction == 'GSGcentral':
                # take only orfs in between FIRST and LAST CBG in GSG (can be only one CBG!)
                max_orf_nt_start = max(GSgraph[-1].overall_minimal_spanning_range(organism=geneQ)) * 3
                min_orf_nt_end   = min(GSgraph[0].overall_minimal_spanning_range(organism=geneQ)) * 3
                orfQlist = input[geneQ]['orfs'].get_elegiable_orfs(
                        acceptorfids=encounteredQorfids,min_orf_end=min_orf_nt_end,
                        max_orf_start=max_orf_nt_start)
            else:
                # dbfraction equals 'all' or None -> no limitation, just take all orfs!
                # only to the general limimation on previously encountered orf nodes
                # this is in fact an exception -> wrong use of parameters!
                orfQlist = input[geneQ]['orfs'].get_elegiable_orfs(acceptorfids=encounteredQorfids)

        elif dbfraction == 'annotation':
            # take only the subset of orfs in the current gene's annotation
            orfQlist = [ input[geneQ]['orfs'].get_orf_by_id(id) for id in input[geneQ]['orfid-genestructure'] ]
            # if no orfid-genestructure known -> just pretend a kind of ab initio
            # and thus take BLASTP_LONGEST_ORFS
            if not orfQlist:
                orfQlist = _order_list_by_attribute(input[geneQ]['orfs'].orfs,order_by='length',reversed=True)[0:blastoptions.BLASTP_LONGEST_ORFS]

        elif dbfraction == 'iterative':
            # take subset of Query Orfs found in previous blasts (in crossdata)
            acceptorfids = Set()
            for _geneQ,_geneS in ordered_edges:
                for (bits,length,orfQid,orfSid) in crossdata[(_geneQ,_geneS)]['accepted_pacbs'].keys():
                    acceptorfids.add(orfQid)
            _cnt = len(acceptorfids)
            for (bits,length,orfQid,orfSid) in crossdata[(geneQ,geneS)]['accepted_pacbs'].keys():
                if orfQid in acceptorfids:
                    acceptorfids.remove(orfQid)
            orfQlist = input[geneQ]['orfs'].get_elegiable_orfs(acceptorfids=list(acceptorfids))

        elif dbfraction == 'HCP':
            # Take High Coding Potential Orfs:
            # take *not* TCODE non-coding and at least of a certain length
            orfQlist = input[geneQ]['orfs'].get_elegiable_orfs(
                    tcode_is_noncoding = False,
                    min_orf_length=250,
                    rejectorfids=input[geneQ]['orfid-genestructure'] )

        elif dbfraction == 'longest':
            # take the subset of BLAST_LONGEST_ORFS; somewhat ~20 orfs !?
            orfQlist = _order_list_by_attribute(input[geneQ]['orfs'].orfs,order_by='length',reversed=True)[0:blastoptions.BLASTP_LONGEST_ORFS]

        else:
            # no GSgraph yet -> limitation on GSG makes no sense
            # no (recognized) dbfraction -> just take all orfs!
            orfQlist = input[geneQ]['orfs'].orfs

        # perform the actual crossblast
        crossdata, hitcnt = _blastorfset2blastdb(geneQ,geneS,input[geneS]['blastdb'],input,crossdata,GSgraph,
            blastoptions=blastoptions,
            elegiable_orfsQ=orfQlist,
            logging=False)
        # set the counter values
        crossblasthitcnt[(geneQ,geneS)]['orfsQ'] = len(orfQlist)
        crossblasthitcnt[(geneQ,geneS)]['hits']  = hitcnt

        if GTgraph:
            # Do here already a first linearization. Orf nodes after this
            # pair wise blast (started with the highest scoring one!) 
            # are the only orfs that are used to search for in newer blasts
            pacbps = crossdata[(geneQ,geneS)]['accepted_pacbs']
            prev = len(pacbps)
            ( accepted, nonlinear ) = pacb.linearization.linearise_pacbs(pacbps,
                    ACCEPTANCE_MARGIN=LINEARIZATION_ACCEPTANCE_MARGIN,
                    start_with_bests=LINEARIZATION_STARTWITHBEST,
                    is_weigthed=LINEARIZATION_WEIGTHED)
            # update crossdata dictionaries
            crossdata[(geneQ,geneS)]['accepted_pacbs'] = accepted
            crossdata[(geneQ,geneS)]['rejected_pacbs_nl'].update(nonlinear)

            # Get lists of orfids from crossdata
            # When this list is requested for the first time,
            # make a new (smaller) blast database!
            if not encounteredOrfNodeIds.has_key(geneQ):
                encounteredOrfNodeIds[geneQ] = _get_orfids_from_crossdata(geneQ,crossdata)
                # remove current existing database of this informant 
                _blastdb_cleanup( { geneS: input[geneQ] })
                # and make new ones
                sbjctcnts = createblastdbs(input,GSgraph,OPTIONS,dbfraction=dbfraction,
                        organism=geneQ,acceptorfids=encounteredOrfNodeIds[geneQ])
                if logging: print "new blastdb: ", geneQ, sbjctcnts
            if not encounteredOrfNodeIds.has_key(geneS):
                encounteredOrfNodeIds[geneS] = _get_orfids_from_crossdata(geneS,crossdata)
                # remove current existing database of this informant 
                _blastdb_cleanup( { geneS: input[geneS] })
                # and make new ones
                sbjctcnts = createblastdbs(input,GSgraph,OPTIONS,dbfraction=dbfraction,
                        organism=geneS,acceptorfids=encounteredOrfNodeIds[geneS])
                if logging: print "new blastdb: ", geneS, sbjctcnts

            # print message for how big the limitation is/was
            if logging: print (geneQ,geneS), prev, len(accepted), len(nonlinear), "orfids:", len(encounteredOrfNodeIds[geneQ]), len(encounteredOrfNodeIds[geneS])
 
    # All crossblastp steps are done.
    # Return the input data structure.
    return crossdata, crossblasthitcnt

# end of function iterativecrossblastp


def _get_orfids_from_crossdata(geneid,crossdata):
    """
    Get all orf ids from a gene/organism identifier from crossdata dict structure

    @type  geneid: gene/organism identifier
    @param geneid: gene/organism identifier

    @type  crossdata: dict
    @param crossdata: `crossdata` dictionary with the gathered Pacbps in previous iterations

    @return: list
    @return: list with orf ids (integers)
    """
    orfids = Set() 
    for g1,g2 in crossdata.keys():
        if geneid == g1:   pointer = 2
        elif geneid == g2: pointer = 3
        else:
            continue 
        # get orfids from the keys by pointer
        for key in crossdata[(g1,g2)]['accepted_pacbs'].keys(): orfids.add( key[pointer] )
    return list(orfids)    

# end of function _get_orfids_from_crossdata


def _order_list_by_attribute(inputlist,order_by='',reversed=False):
    """
    Helper function to order objects in a list by any given attribute
    """
    if not inputlist: return inputlist
    if hasattr(inputlist[0],str(order_by)):
        orderedlist = [ ( getattr(o,str(order_by)), o ) for o in inputlist ]
        orderedlist.sort()
        if reversed: orderedlist.reverse()
        return [ o for (attr,o) in orderedlist ]
    else:
        return inputlist
    
# end of function _order_list_by_attribute


def strip_exterior_gaps(hsp):
    """ When SEVERELY ERRONEOUS (!!) BlastHSP is encountered, it has exterior gaps & spaces!"""
    while hsp.query[0] == "-":
        hsp.query = hsp.query[1:].strip()
        hsp.match = hsp.match[1:].strip()
        hsp.sbjct = hsp.sbjct[1:].strip()
        hsp.sbjct_start+=1
    while hsp.sbjct[0] == "-":
        hsp.sbjct = hsp.sbjct[1:].strip()
        hsp.match = hsp.match[1:].strip()
        hsp.query = hsp.query[1:].strip()
        hsp.query_start+=1
    #while self.query[-1] == "-":
    #    self.query_end+=1
    #    self.query = self.query[0:-1]
    #while self.sbjct[-1] == "-":
    #    self.sbjct_start+=1
    #    self.sbjct = self.sbjct[1:]

# end of function strip_exterior_gaps


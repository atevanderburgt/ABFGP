"""
"""
from sets import Set
import pacb
from graphAbgp import (
    InwardsPointingCodingBlockGraph,
    GeneTreeGraph,
    )
from graphAbgp.graph_pacbpcollection import _delete_pacbp

from abgp_exceptions import InproperlyAppliedArgument

# Global variable Imports
from settings.genestructure import (
    ORF_IS_UNIGENE_LABEL,
    INFORMANT_HIGH_NT_IDENTITY,
    )
SINGLETON_BLOCK_MINIMAL_CNT   = 4
SINGLETON_BLOCK_MINIMAL_RATIO = 0.20
from settings.executables import TCODE_MAX_NONCODING


def _pacbporfs2wt_by_bitscore(pacbporflist):
    """ """
    if not pacbporflist: return 0
    return sum([pf.bitscore for pf in pacbporflist])

# end of function _pacbporfs2wt_by_bitscore


def _pacbporfs2wt_by_ntidentity(pacbporflist):
    """ """
    if not pacbporflist: return 0.0
    ntidentities = []
    lengths = []
    for pf in pacbporflist:
        idperc = pf.get_nt_identity()
        length = pf.get_unextended_length()
        ntidentities.append( idperc * length )
        lengths.append( length )
    ntidentity = sum(ntidentities) / sum(lengths)
    return ntidentity

# end of function _pacbporfs2wt_by_ntidentity


def _pcg2gtg(pcg,target,identifier_list=[],by='bitscore'):
    """ """
    # make a (artificially fully connected) GeneTreeGraph
    gtg = GeneTreeGraph()
    gtg.add_node(target)
    for informant in pcg.organism_set():
        if informant == target: continue
        if identifier_list and informant not in identifier_list: continue
        pacbporfs = pcg.get_pacbps_by_organisms(target,informant)
        if by == 'bitscore':
            weight = _pacbporfs2wt_by_bitscore(pacbporfs)
        elif by == 'ntidentity':
            weight = _pacbporfs2wt_by_ntidentity(pacbporfs)
        else:
            message = "by not in 'bitscore,ntidentity' but '%s'" % by
            raise InproperlyAppliedArgument,message
        # add node and edge with weight
        gtg.add_node(informant)
        gtg.add_edge( target,informant, wt = weight)

    # make artificially missed edges between the informants
    for informant in pcg.organism_set():
        if informant == target: continue
        if identifier_list and informant not in identifier_list: continue
        # get weight of target-informant combination
        weight = gtg.weights[(target,informant)]
        for org in pcg.organism_set():
            if org == target: continue
            if org == informant: continue
            if identifier_list and org not in identifier_list: continue
            if not gtg.has_edge( informant, org ):
                gtg.add_edge(informant,org,wt=weight)
            elif gtg.has_edge( informant, org ) and\
            gtg.weights[(informant, org)] > weight:
                # update edge weight by the LOWEST weight
                gtg.set_edge_weight(informant,org,wt=weight)
            else:
                pass

    #return the gtg
    return gtg

# end of function _pcg2gtg


def pcg2gtg_by_bitscore(pcg,target,identifier_list=[]):
    """ """
    return _pcg2gtg(pcg,target,identifier_list=identifier_list,by='bitscore')

# end of function pcg2gtg_by_bitscore


def pcg2gtg_by_ntidentity(pcg,target,identifier_list=[]):
    """ """
    return _pcg2gtg(pcg,target,identifier_list=identifier_list,by='ntidentity')

# end of function pcg2gtg_by_ntidentity


def pcg2gtg_by_identity(pcg,target,identifier_list=[]):
    """ """
    return _pcg2gtg(pcg,target,identifier_list=identifier_list,by='ntidentity')

# end of function pcg2gtg_by_identity


def DEPRE_pcg2gtg_by_identity(pcg,target,identifier_list=[]):
    """ """
    # make a (artificially fully connected) GeneTreeGraph
    gtg = GeneTreeGraph()
    gtg.add_node(target)
    for informant in pcg.organism_set():
        if informant == target: continue
        if identifier_list and informant not in identifier_list: continue
        pacbporfs = pcg.get_pacbps_by_organisms(target,informant)
        ntidentities = []
        lengths = []
        for pf in pacbporfs:
            idperc = pf.get_nt_identity()
            length = pf.get_unextended_length()
            ntidentities.append( idperc * length )
            lengths.append( length )
        ntidentity = sum(ntidentities) / sum(lengths)
        gtg.add_node(informant)
        gtg.add_edge( target,informant, wt = ntidentity )

    # make artificially missed edges between the informants
    for informant in pcg.organism_set():
        if informant == target: continue
        if identifier_list and informant not in identifier_list: continue
        for org in pcg.organism_set():
            if org == target: continue
            if org == informant: continue
            if identifier_list and org not in identifier_list: continue
            if not gtg.has_edge( informant, org ):
                gtg.add_edge(informant,org,wt=ntidentity)
            elif gtg.has_edge( informant, org ) and\
            gtg.weights[(informant, org)] > ntidentity:
                # update edge weight by the LOWEST nt-identity
                gtg.set_edge_weight(informant,org,wt=ntidentity)
            else:
                pass

    #return the gtg
    return gtg

# end of function pcg2gtg_by_identity


def pacbplist2queryMSR(pacbplist):
    """ """
    # first, create an query MSR for organisms that are covered only once
    spos = pacbplist[0]._get_original_alignment_pos_start()
    epos = pacbplist[0]._get_original_alignment_pos_end()
    startQ, endQ = spos.query_pos, epos.query_pos + 1
    MSR = Set(range(startQ,endQ))
    for otherpacbp in pacbplist[1:]:
        spos = otherpacbp._get_original_alignment_pos_start()
        epos = otherpacbp._get_original_alignment_pos_end()
        startQ, endQ = spos.query_pos, epos.query_pos + 1
        MSR = MSR.intersection(range(startQ,endQ))
    return MSR

# end of function pacbplist2queryMSR


def PCG2inwpCBGS(PCG,omit_singletons=False):
    """ Translate the PCG to list of InwardsPointingCodingBlockGraph """
    blocks = PCG2blocks(PCG)    
    omit_block_positions = []
    if omit_singletons:
        # Find singletons PacbPORFs -> match to blocks
        # and store the list index of these in omit_block_positions
        # This excludes very poor / 5' / 3' blocks from the PCG
        organisms = list(PCG.organism_set())
        omit_block_positions = select_singleton_block_index(blocks,organisms)

    inwpcbgs = []
    for pos in range(0,len(blocks)):
        if pos in omit_block_positions: continue
        block = blocks[pos]
        iwpCBG = InwardsPointingCodingBlockGraph()
        (pacbpkey,nodeQ,nodeS) = block[0]
        iwpCBG.add_node(nodeQ)
        for key in block:
            (pacbpkey,nodeQ,nodeS) = key
            pacbporf = PCG.pacbps[(pacbpkey,nodeQ,nodeS)]
            iwpCBG.add_node(nodeS)
            iwpCBG.add_edge(nodeQ,nodeS,wt=pacbpkey[0])
            iwpCBG.pacbps[key] = PCG.pacbps[key]
        inwpcbgs.append( iwpCBG )
    # return list of inwpcbgs
    return inwpcbgs

# end of function PCG2inwpCBGS



def _select_overlapping_pacbps_based_on_pacbp_seed(seedlist,pacbp_keys,PCG):
    """ Helper function for PCG2blocks """

    # return list with overlapping keys
    overlapping_keys = []

    # select all pacbps that share query overlap with the seed
    for k2 in pacbp_keys:
        if k2 in seedlist: continue
        otherpacbp = PCG.pacbps[k2]
        # get minimal overlap with (one of the) seed pacbps
        overlap = min([ PCG.pacbps[k1].query_overlap(otherpacbp) for k1 in seedlist ])
        if overlap > 0.10: # at least 10% overlap; was 0 percent!
            overlapping_keys.append( k2 )

    ##########################################################
    # get rid of pacbps/edges not belonging to the
    # query or sbjct node (other Orf)
    ##########################################################
    del_keys = []
    for (blastkey,nodeQ,nodeS) in overlapping_keys:
        if nodeQ not in [ seedNodeQ for (k,seedNodeQ,nS) in seedlist ]:
            del_keys.append((blastkey,nodeQ,nodeS))
            continue
        (orgSid,orfSid) = nodeS
        for (k,seedNodeQ,seedNodeS) in seedlist:
            if orgSid == seedNodeS[0] and orfSid != seedNodeS[1]:
                del_keys.append((blastkey,nodeQ,nodeS))
                break
    # delete pacbpkeys labeled for deletion from overlapping_keys
    for del_key in del_keys: overlapping_keys.remove(del_key)

    # insert the seed(s) at the begin of overlapping_keys
    for k1 in seedlist: overlapping_keys.insert(0,k1)

    # return list of overlap
    return overlapping_keys

# end of function _select_overlapping_pacbps_based_on_pacbp_seed


def PCG2blocks(PCG, verbose = False):
    """ """
    ############################################################################
    if verbose:
        print PCG
        target_org = PCG.pacbps.keys()[0][1][0]
        for nodeQ,nodeS in PCG.weights.keys():
            if nodeQ[0] != target_org: continue
            pacbporfs = PCG.get_pacbps_by_nodes(node1=nodeQ,node2=nodeS)
            if len(pacbporfs) >= 2:
                print nodeQ,nodeS
                for pacbporf in pacbporfs: print pacbporf
    ############################################################################
    pacbp_keys = PCG.pacbps.keys()
    pacbp_keys.sort()
    blocks = {}
    for k1 in pacbp_keys:
        theblastkey,thenodeQ,thenodeS = k1
        seen = False
        for bl in blocks.values():
            if k1 in bl:
                seen = True
                break
        if seen: continue

        # seed for a new block / inwpCBG. Initialize list with single pacbpkey
        thepacbp = PCG.pacbps[k1]
        overlapping_keys = [ k1 ]

        # get list of overlapping pacbp keys based on this seed
        seedlist = [ k1 ]
        overlapping_keys = _select_overlapping_pacbps_based_on_pacbp_seed(seedlist,pacbp_keys,PCG)

        # check informant organism coverage in list of overlapping_keys
        sbjctOrgs = [ orgS for k,(orgQ,orfQ),(orgS,orfS) in overlapping_keys ]

        if len(Set(sbjctOrgs)) == len(sbjctOrgs):
            # print the content of this block
            blockMSR = pacbplist2queryMSR([PCG.pacbps[k] for k in overlapping_keys])
            if not blockMSR:
                continue

            # store to blocks
            block_key = ( min(blockMSR), max(blockMSR) )
            blocks[block_key] = overlapping_keys

            # break out by continueing
            continue


        ###############################################################
        # If here, then there are duplicate organism keys in the list;
        # throw out all except the one that is contributing to to MSR
        ###############################################################
        # first, create an query MSR for organisms that are covered only once
        spos = thepacbp._get_original_alignment_pos_start()
        epos = thepacbp._get_original_alignment_pos_end()
        startQ, endQ = spos.query_pos, epos.query_pos + 1
        theMSR = Set(range(startQ,endQ))
        for pacbpkey in overlapping_keys[1:]:
            blastkey,nodeQ,(orgS,orfS) = pacbpkey
            if sbjctOrgs.count(orgS) == 1:
                otherpacbp = PCG.pacbps[pacbpkey]
                spos = otherpacbp._get_original_alignment_pos_start()
                epos = otherpacbp._get_original_alignment_pos_end()
                startQ, endQ = spos.query_pos, epos.query_pos + 1
                theMSR = theMSR.intersection(range(startQ,endQ))

        if not theMSR:
            ################################################################
            if verbose:
                print " WARNING!!!! NO MSR REMAINING HERE.....",
                print thepacbp, thenodeQ,thenodeS
                for pacbpkey in overlapping_keys[1:]:
                    blastkey,nodeQ,(orgS,orfS) = pacbpkey
                    otherpacbp = PCG.pacbps[pacbpkey]
                    if sbjctOrgs.count(orgS) == 1:
                        print (orgS,orfS), True, otherpacbp
                for pacbpkey in overlapping_keys[1:]:
                    blastkey,nodeQ,(orgS,orfS) = pacbpkey
                    otherpacbp = PCG.pacbps[pacbpkey]
                    if sbjctOrgs.count(orgS) != 1:
                        print (orgS,orfS), False, otherpacbp
            ################################################################
            # solve no MSR here. It is caused because 2 informants do
            # not overlap on target coordinates. Find these cases
            # by interating over the overlap selection by taking these
            # non-overlapping informant pacbps together with the
            # initial (other organism) informant seed pacbp
            msr_breaking_keys = Set()
            for p1 in range(1,len(overlapping_keys)):
                pacbpkey1 = overlapping_keys[p1]
                blastkey1,nodeQ1,(orgS1,orfS1) = pacbpkey1
                if sbjctOrgs.count(orgS1) != 1: continue
                otherpacbp1 = PCG.pacbps[pacbpkey1]
                for p2 in range(1,len(overlapping_keys)):
                    if p2 <= p1: continue
                    pacbpkey2 = overlapping_keys[p2]
                    if pacbpkey1 == pacbpkey2: continue
                    blastkey2,nodeQ2,(orgS2,orfS2) = pacbpkey2
                    if sbjctOrgs.count(orgS2) != 1: continue
                    otherpacbp2 = PCG.pacbps[pacbpkey2]
                    overlap = otherpacbp1.query_overlap(otherpacbp2)
                    if not overlap:
                        msr_breaking_keys.add(pacbpkey1)
                        msr_breaking_keys.add(pacbpkey2)
                        ########################################################
                        if verbose:
                            print "NO::", overlap, (orgS1,orfS1), (orgS2,orfS2)
                        ########################################################

            for second_seed_key in msr_breaking_keys:
                # check if this msr_breaking_key is in a previous iteration
                # through this exact code resolved. In that case, it would
                # create duplicate blocks here -> HELL!!
                # Check this here...
                seen = False
                for bl in blocks.values():
                    if k1 in bl and second_seed_key in bl:
                        seen = True
                        break
                if seen: continue

                new_seedlist = [ k1, second_seed_key ]
                resolved_overlapping_keys = _select_overlapping_pacbps_based_on_pacbp_seed(new_seedlist,pacbp_keys,PCG)

                # check informant organism coverage in list of overlapping_keys
                resolved_sbjctOrgs = [ orgS for k,(orgQ,orfQ),(orgS,orfS) in resolved_overlapping_keys ]
                resolved_blockMSR = pacbplist2queryMSR([PCG.pacbps[k] for k in resolved_overlapping_keys])

                if len(Set(resolved_sbjctOrgs)) == len(resolved_sbjctOrgs)\
                and len(resolved_blockMSR):
                    # yep, resolved! Store to blocks
                    block_key = ( min(resolved_blockMSR), max(resolved_blockMSR) )
                    blocks[block_key] = resolved_overlapping_keys
                    ############################################################
                    if verbose:
                        print "SECOND::", second_seed_key,
                        print len(Set(resolved_sbjctOrgs)) == len(resolved_sbjctOrgs),
                        print len(resolved_blockMSR)
                    ############################################################
                else:
                    # hmm... this did not help....
                    # We expect that this is not possible to happen.
                    # And, if it does, it is very likely not something that
                    # will contribute to proper gene structure prediction.
                    # So, leave it as it is....
                    print "NON-OVERLAP SEARCH DID NOT HELP FOR::", second_seed_key

        else:
            # yes, this resulted in a clear-cut MSR.
            # select pacbp from organisms that are present twice or more
            # by checking overlap with the MSR
            org2overlap = {}
            for pacbpkey in overlapping_keys[1:]:
                blastkey,nodeQ,(orgS,orfS) = pacbpkey
                if sbjctOrgs.count(orgS) >= 2:
                    if not org2overlap.has_key(orgS): org2overlap[orgS] = []
                    otherpacbp = PCG.pacbps[pacbpkey]
                    spos = otherpacbp._get_original_alignment_pos_start()
                    epos = otherpacbp._get_original_alignment_pos_end()
                    startQ, endQ = spos.query_pos, epos.query_pos + 1
                    msroverlap = theMSR.intersection(range(startQ,endQ))
                    if msroverlap:
                        overlap_ratio = float(len(msroverlap)) / len(theMSR)
                    else:
                        overlap_ratio = 0.0
                    org2overlap[orgS].append((overlap_ratio,pacbpkey))
            # remove all except the highest overlap ratio per organism
            for orgS,overlap_list in org2overlap.iteritems():
                overlap_list.sort()
                # remove all except the final (==highest overlap)
                for overlap_ratio,del_key in overlap_list[0:-1]:
                    overlapping_keys.remove(del_key)

            # print the content of this block
            blockMSR = pacbplist2queryMSR([PCG.pacbps[k] for k in overlapping_keys])
            if not blockMSR:
                # WARNING!!!! NO MSR REMAINING HERE.....
                print "NO MSR REMAINING!!!!!"
                for _key in overlapping_keys:
                    print "NOMSR:", _key
                    print "NOMSR:", PCG.pacbps[_key]
                continue

            else:
                # store to blocks
                block_key = ( min(blockMSR), max(blockMSR) )
                blocks[block_key] = overlapping_keys

    ########################################################################
    if verbose:
        print "# RESULT OF PCG2blocks():"
        for msr,overlapping_keys in blocks.iteritems():
            print "BLOCK::", len(overlapping_keys), overlapping_keys[0],
            print len(msr), msr
            for k in overlapping_keys:
                otherpacbp = PCG.pacbps[k]
                print k, thepacbp.query_overlap(otherpacbp), thepacbp.is_coding()
            if len(overlapping_keys) == 1:
                print "SINGLETON::", PCG.pacbps[overlapping_keys[0]]
            print ""
    ########################################################################

    # return ORDERED blocks
    block_keys = blocks.keys()
    block_keys.sort()
    return [ blocks[k] for k in block_keys ]

# end of function PCG2blocks


def blocks2orfidstruct(blocks,PCG,ORGANISMS,OPTIONS):
    """ """
    orfidstruct = {}
    # blockMSRs contains the MSR regions of the blocks
    blockMSRs = []
    for block in blocks:
        blockMSR = pacbplist2queryMSR([PCG.pacbps[k] for k in block ])
        orfidlist = []
        for org in ORGANISMS:
            if org == OPTIONS.target:
                orfidlist.append( block[0][1][1] )
            else:
                for blastkey,nodeQ,(orgS,orfS) in block:
                    if orgS == org:
                        orfidlist.append( orfS )
                        break
                else:
                    orfidlist.append("-")
        key = min(blockMSR),max(blockMSR)
        orfidstruct[key] = orfidlist
    if float(len(block))/len(ORGANISMS) > 0.50:
        # convincing enough to be part of the Exon structure
        blockMSRs.append(key)
    return orfidstruct, blockMSRs

# end of function blocks2orfidstruct


def print_inwpcbgstructure(inwpcbglist,ORGANISMS,GTG=None):
    """
    """
    if not inwpcbglist: return False
    target = inwpcbglist[0]._get_target_organism()
    separator = "########"*(len(inwpcbglist)+2)

    # check which organisms are unigenes
    UNIGENE_INFORMANT_LIST = []
    GENE_INFORMANT_LIST = []
    for org in ORGANISMS:
        if org == target: continue
        for inwpCBG in inwpcbglist:
            if org in inwpCBG.organism_set():
                orfobj = inwpCBG.get_orfs_of_graph(organism=org)[0]
                if hasattr(orfobj,ORF_IS_UNIGENE_LABEL):
                    UNIGENE_INFORMANT_LIST.append(org)
                    break
                else:
                    GENE_INFORMANT_LIST.append(org)
                    break

    print separator
    data = []
    for inwpCBG in inwpcbglist:
        if inwpCBG.has_minimal_spanning_range():
            minsr = inwpCBG.minimal_spanning_range(organism=target)
            maxsr = inwpCBG.maximal_spanning_range(organism=target)
            data.append((min(minsr),max(minsr),min(maxsr),max(maxsr)))
        else:
            data.append(("-","-","-","-"))
    for elemA,elemB,elemC,elemD in data:
        print elemA,"\t",
    print "MINSR (MIN)"
    for elemA,elemB,elemC,elemD in data:
        print elemB,"\t",
    print "MINSR (MAX)"
    print separator
    for pos in range(0,len(data)):
        elemA,elemB,elemC,elemD = data[pos]
        if pos > 0 and pos < len(data)-1 and\
        (elemC,elemD) == data[pos-1][2:] and (elemC,elemD) == data[pos+1][2:]:
            print "-\t",
        else:
            print elemC,"\t",
    print "MAXSR (MIN)"
    for pos in range(0,len(data)):
        elemA,elemB,elemC,elemD = data[pos]
        if pos > 0 and pos < len(data)-1 and\
        (elemC,elemD) == data[pos-1][2:] and (elemC,elemD) == data[pos+1][2:]:
            print "-\t",
        else:
            print elemD,"\t",
    print "MAXSR (MAX)"

    for inwpCBG in inwpcbglist:
        score = inwpCBG.get_exon_uniformity_score()
        if score == None: print "None\t",
        else:             print "%1.3f\t" % score, 
    print "ExonUniformity"
    print separator
    for inwpCBG in inwpcbglist:
        print inwpCBG.has_target_orf_methionine(),"\t",
    print "methionine?"
    for inwpCBG in inwpcbglist:
        print inwpCBG.have_informant_orfs_methionines(),"\t",
    print "alignedMs?"
    for inwpCBG in inwpcbglist:
        print inwpCBG.count_orfs_with_methionines(),"\t",
    print "MethionineCnt"    
    for inwpCBG in inwpcbglist:
        print [ pf.has_upstream_methionines() for pf in inwpCBG.pacbps.values() ].count(True),"\t",
    print "UpstrMethionineCnt"
    for inwpCBG in inwpcbglist:
        print [ pf.has_upstream_tss() for pf in inwpCBG.pacbps.values() ].count(True),"\t",
    print "UpstrTSSCnt"
    for inwpCBG in inwpcbglist:
        print "%2.1f\t" % inwpCBG.get_average_upstream_methionine_pssm_score(),
    print "UpstrMethPSSM"
    for inwpCBG in inwpcbglist:
        print inwpCBG.count_orfs_with_signalpeptides(),"\t",
    print "SignalpCnt"
    for inwpCBG in inwpcbglist:
        score = inwpCBG.get_signalp_score()
        if score: print "%1.2f\t" % (score),
        else:     print "None\t",
    print "SignalpScore"
    for inwpCBG in inwpcbglist:
        print inwpCBG.count_orfs_labeled_as_first_exon(),"\t",
    print "IsFirstOrf?"
    for inwpCBG in inwpcbglist:
        print inwpCBG.count_orfs_labeled_as_annotated_exon(),"\t",
    print "IsAnnotatedOrf?"
    for inwpCBG in inwpcbglist:
        print inwpCBG.count_orfs_labeled_as_final_exon(),"\t",
    print "IsFinalOrf?"
    for inwpCBG in inwpcbglist:
        d1 = inwpCBG.get_projected_tailing_stop_aa_difference()
        d2 = inwpCBG.get_projected_tailing_stop_nonaligned_aa_difference()
        if d1 >= 1000: d1 = "%se3" % int(round(float(d1)/1000))
        if d2 >= 1000: d2 = "%se3" % int(round(float(d2)/1000))
        print "%s,%s\t" % ( d1, d2 ),
    print "TailingStopDist"
    print separator


    if GTG:
        total_wt = 0.0
        for informant in GTG.organism_set():
            if informant == target: continue
            total_wt+=GTG.weights[(target,informant)]
        for inwpCBG in inwpcbglist:
            gtg = pcg2gtg_by_identity(inwpCBG,target,identifier_list=GTG.organism_set())
            difference = 0.0
            for informant in GTG.organism_set().intersection(gtg.organism_set()):
                if informant == target: continue
                wtGTG = GTG.weights[(target,informant)]
                wtgtg = gtg.weights[(target,informant)]
                difference += ( wtgtg - wtGTG )
            for informant in GTG.organism_set().difference(gtg.organism_set()):
                wtGTG = GTG.weights[(target,informant)]
                difference -= wtGTG
            print "%1.2f\t" % ( difference / total_wt ),
        print "GTGdiff"

        for inwpCBG in inwpcbglist:
            print "%s\t" % inwpCBG.get_cexpander_uniformly_aligned_count(),
        print "CEXPANDER"

        print separator

    # print sumed bitscore based on MSXSR; measure for trustability
    for inwpCBG in inwpcbglist:
        print "%s\t" % inwpCBG.get_bitscore(),
    print "SUM(bits)"


    for inwpCBG in inwpcbglist:
        print "%1.2f\t" % inwpCBG.get_bits_per_aa(),
    print "BITS/AA)"


    # print !!ESTIMATED!! overall identity (measured on MAXSR, not OMSR!)
    for inwpCBG in inwpcbglist:
        print "%1.2f\t" % inwpCBG.get_identityscore(),
    print "ID%estimate"

    for inwpCBG in inwpcbglist:
        try:    print "%1.2f\t" % inwpCBG.tcode(),
        except: print "n.a.\t",
    print "TCODE"

    print separator

    # min/avr/max nt distance between CBGs
    _data = []
    for pos in range(1,len(inwpcbglist)):
        this = inwpcbglist[pos-1]
        nextlist = inwpcbglist[pos:]
        dist = this.nt_spacing_between_codingblocks(nextlist).values()
        _data.append( ( min(dist), sum(dist)/len(dist), max(dist) ) )

    for value in [ a for a,b,c in _data ]:
        print "%s\t" % value, 
    print "-\tnt-DIST (min)"
    for value in [ b for a,b,c in _data ]:
        print "%s\t" % value, 
    print "-\tnt-DIST (avr)"
    for value in [ c for a,b,c in _data ]:
        print "%s\t" % value, 
    print "-\tnt-DIST (max)"

    print separator

    _data = []
    for pos in range(1,len(inwpcbglist)):
        this = inwpcbglist[pos-1]
        nextlist = inwpcbglist[pos:]
        tcodedist = this.tcode_spacing_between_codingblocks(nextlist)
        _data.append(tcodedist)
    for tcodedist in _data:
        if tcodedist:
            print "%s/%s\t" % (
                [ score <= TCODE_MAX_NONCODING for score in tcodedist.values()].count(True),
                len(tcodedist) ),
        else:
            print "0\t",
    print "-\ttcode-DIST (cnt)"
    for tcodedist in _data:
        if tcodedist: print "%1.2f\t" % (sum(tcodedist.values())/len(tcodedist)),
        else:         print "None\t",
    print "-\ttcode-DIST (avr)"
    for tcodedist in _data:
        if tcodedist.has_key(target): print "%1.2f\t" % (tcodedist[target]),
        else:                         print "None\t",
    print "-\ttcode-DIST (target)"

    print separator


    # count annotated Orf count transitions

    _data = []
    for inwpCBG in inwpcbglist:
        _data.append( set(inwpCBG.get_orfs_labeled_as_annotated_exon().keys()) )
        
    cnts = []
    for offset in range(0,len(inwpcbglist)):
        cnt5p = set()
        cnt3p = set()
        for pos in range(0,offset+1): cnt5p.update(_data[pos])
        for pos in range(offset+1,len(inwpcbglist)): cnt3p.update(_data[pos])
        cnts.append( ( len(cnt5p), len(cnt3p) ) ) 

    for elem in cnts:  print "%s,%s\t" % elem,
    print "KnownOrfTransition"

    print separator


    for ugid in UNIGENE_INFORMANT_LIST:
        for inwpCBG in inwpcbglist:
            if ugid in inwpCBG.organism_set():
                orfid = inwpCBG.get_orfs_of_graph(organism=ugid)[0].id
                print orfid,"\t",
            else:
                print "-\t",
        print ugid
    print separator
    for inwpCBG in inwpcbglist:
        orfid = inwpCBG.get_orfs_of_graph(organism=target)[0].id
        print orfid,"\t",
    print target
    print separator
    for org in ORGANISMS:
        if org == target: continue
        if org in UNIGENE_INFORMANT_LIST: continue
        for inwpCBG in inwpcbglist:
            if org in inwpCBG.organism_set():
                orfid = inwpCBG.get_orfs_of_graph(organism=org)[0].id
                print orfid,"\t",
            else:
                print "-\t",
        try:
            print org, "\t%1.2f" % GTG.weights[(target,org)]
        except:
            print org
    print separator

# end of function print_inwpcbgstructure


def print_orfidstruct(orfidstruct,input,ORGANISMS,OPTIONS):
    """ """
    keys = orfidstruct.keys()
    keys.sort()
    print "########"*(len(keys)+2)
    for k in keys: print k[0],"\t",
    print ""
    for k in keys: print k[1],"\t",
    print ""
    print "########"*(len(keys)+2)
    for k in keys:
        print len(orfidstruct[k])-orfidstruct[k].count("-")-1,"\t",
    print ""
    print "########"*(len(keys)+2)
    for pos in range(0,len(ORGANISMS)):
        org  = ORGANISMS[pos]
        if input[org].has_key('is_unigene') and input[org]['is_unigene']:
            for k in keys:
                print orfidstruct[k][pos],"\t",
            print org
    print "########"*(len(keys)+2)

    # is there a Methionine in the TARGET gene on this Orf?
    for pos in range(0,len(ORGANISMS)):
        org  = ORGANISMS[pos]
        if org == OPTIONS.target:
            for k in keys:
                thisorf = input[org]['orfs'].get_orf_by_id( orfidstruct[k][pos] )
                print thisorf.has_methionine(),"\t",
            print "methionine?"

    # are there a Methionines in the INFORMANT genes on this Orf?
    counts = []
    for k in keys:
        methionine_check = []
        for pos in range(0,len(ORGANISMS)):
            org  = ORGANISMS[pos]
            if input[org].has_key('is_unigene') and input[org]['is_unigene']:
                continue
            elif org == OPTIONS.target:
                continue
            elif orfidstruct[k][pos] != "-":
                thisorf = input[org]['orfs'].get_orf_by_id( orfidstruct[k][pos] )
                methionine_check.append( thisorf.has_methionine() )
            else:
                pass
        counts.append( methionine_check.count(True) )
        if methionine_check.count(True) == len(methionine_check):
            print True,"\t",
        elif methionine_check.count(False) == len(methionine_check):
            print None,"\t",
        else:
            print False,"\t",
    print "alignedMs?"
    for pos in range(0,len(keys)):
        print counts[pos],"\t",
    print "MethionineCnt"


    # check if this Orf is annotated as the FirstExon.
    for k in keys:
        cnt = 0
        for pos in range(0,len(ORGANISMS)):
            if orfidstruct[k][pos] == "-":
                pass
            else:
                org = ORGANISMS[pos]
                if input[org]['orfid-genestructure']:
                    if orfidstruct[k][pos] == input[org]['orfid-genestructure'][0]:
                        cnt+=1
        print cnt,"\t",
    print "IsFirstOrf?"


    # check if this Orf is annotated as the FinalExon.
    for k in keys:
        cnt = 0
        for pos in range(0,len(ORGANISMS)):
            if orfidstruct[k][pos] == "-":
                pass
            else:
                org = ORGANISMS[pos]
                if input[org]['orfid-genestructure']:
                    if orfidstruct[k][pos] == input[org]['orfid-genestructure'][-1]:
                        cnt+=1
        print cnt,"\t",
    print "IsFinalOrf?"


    print "########"*(len(keys)+2)
    for pos in range(0,len(ORGANISMS)):
        org  = ORGANISMS[pos]
        if org == OPTIONS.target:
            for k in keys:
                print orfidstruct[k][pos],"\t",
            print org
    print "########"*(len(keys)+2)
    for pos in range(0,len(ORGANISMS)):
        org  = ORGANISMS[pos]
        if input[org].has_key('is_unigene') and input[org]['is_unigene']:
            continue
        elif org == OPTIONS.target:
            continue
        else:
            for k in keys:
                print orfidstruct[k][pos],"\t",
            print org
    print "########"*(len(keys)+2)
    print "## Current annotated gene's orfid-structure"
    for org in ORGANISMS:
        if input[org]['orfid-genestructure']:
            print "##%s\tgene-orfmodel:" % org,
            print input[org]['orfid-genestructure'], input[org]['proteinfref']

# end of function print_orfidstruct


def remove_poorly_covered_unannotated_inwpcbgs(inwpcbgs,PCG,GENE_INFORMANT_SET,verbose=False):
    """ """
    reject_inwpcbg_list = []
    correct_inwpcbg_list = []
    for pos in range(0,len(inwpcbgs)):
        inwpCBG = inwpcbgs[pos]
        if inwpCBG.count_orfs_labeled_as_annotated_exon() <= 1 and\
        float(inwpCBG.organism_set_size())/len(GENE_INFORMANT_SET) <= 0.33:
            reject_inwpcbg_list.append( inwpCBG )
        else:
            correct_inwpcbg_list.append( inwpCBG )

    # get all pacbp keys belonging to the reject_inwpcbg_list inwpcbgs ONLY
    reject_pacbpkeys = []
    for rejectInwpCBG in reject_inwpcbg_list:
        for pacbpkey in rejectInwpCBG.pacbps.keys():
            # check if this pacbpkey is occuring in a non-removed inwpCBG
            is_occurring_in_correct_inwpcbg = False
            for inwp in correct_inwpcbg_list:
                if pacbpkey in inwp.pacbps.keys():
                    is_occurring_in_correct_inwpcbg = True
                    break
            # if is_occurring_in_correct_inwpcbg, continue and do not delete
            if is_occurring_in_correct_inwpcbg:
                continue
            # store to reject_pacbpkeys when not stored already
            if pacbpkey not in reject_pacbpkeys:
                reject_pacbpkeys.append(pacbpkey)

    # delete all PacbPORFs in reject_pacbpkeys from the PCG
    for key in reject_pacbpkeys:
        # remove from main PCG
        _delete_pacbp(PCG,key)


    # return counters of how much inwpcbgs/pacbporfs have been deleted
    return ( len(reject_inwpcbg_list), len(reject_pacbpkeys) )

# end of function remove_poorly_covered_unannotated_inwpcbgs



def remove_noncoding_blocks(blocks,PCG):
    """ Remove blocks which PacbPORFs are ALL non-coding """
    is_any_block_removed = False
    for pos in range(len(blocks)-1,-1,-1):
        block = blocks[pos]
        ARE_ALL_NONCODING = True
        for pacbpkey in block:
            if not PCG.pacbps.has_key(pacbpkey):
                # already deleted in a previous block deletion!
                continue
            pacbporf = PCG.pacbps[pacbpkey]
            # check for high nt identity alignments
            if pacbporf.get_nt_identity() >= INFORMANT_HIGH_NT_IDENTITY:
                ARE_ALL_NONCODING = False
                break
            # do coding_sequence_periodicity_test    
            ( is_lowest_conservation_of_3_positions,
            is_highest_ratio_of_3_frames ) = pacbporf.is_coding()
            if is_highest_ratio_of_3_frames == True:
                ARE_ALL_NONCODING = False
                break
        # check if all are non-coding; if so -> delete complete block!
        if ARE_ALL_NONCODING:
            # A block is labeled for deletion.
            # Update return value to True
            is_any_block_removed = True
            for pacbpkey in block:
                # delete this pacbporf
                _delete_pacbp(PCG,pacbpkey,None)

    # return status weather or not a block was removed
    return is_any_block_removed

# end of function remove_noncoding_blocks


def remove_noncoding_blockorigins(blocks,PCG):
    """ Remove blocks for which its parental PacbPORF is non-coding AND
        removal of this PacbPORF yields a sub-block (of another existing block)
    """
    is_any_block_removed = False
    for pos in range(len(blocks)-1,-1,-1):
        block = blocks[pos]
        parent_pacbporf_key = block[0] 
        if not PCG.pacbps.has_key(parent_pacbporf_key):
            # pacbp/blockorigin already removed in previous iteration
            continue
        # get the pacbporf which is the origin of this block / inwpCBG
        parent_pacbporf = PCG.pacbps[ parent_pacbporf_key ]
        # do coding_sequence_periodicity_test    
        ( is_lowest_conservation_of_3_positions,
        is_highest_ratio_of_3_frames ) = parent_pacbporf.is_coding()
        if is_highest_ratio_of_3_frames or len(block) == 1:
            # ignore this block; parent PacbPORF is probably coding
            continue
    
        # check if the rest of the block could be contained in another block
        REMOVE_PARENT_PACBP = False
        for _pos in range(len(blocks)-1,-1,-1):
            if pos == _pos: continue
            _block = blocks[_pos]
            if not Set(block[1:]).difference(_block):
                # all other PacbPORFs fully included in another block
                REMOVE_PARENT_PACBP = True
                break
        # remove this parent pacbp if labeled as such
        if REMOVE_PARENT_PACBP:
            # A block is labeled for deletion.
            pacbporf = PCG.pacbps[parent_pacbporf_key]
            _delete_pacbp(PCG,parent_pacbporf_key,pacbporf)
            # Update return value to True
            is_any_block_removed = True

    # return status weather or not a block was removed
    return is_any_block_removed

# end of function remove_noncoding_blockorigins


def remove_unlikely_start_blocks(blocks,PCG,organism):
    """
    """
    ############################################################################
    #454     478     481     496     522     529     536     586     649 
    #459     504     488     501     561     536     564     616     678 
    ############################################################################
    #1       1       3       2       1       1       5       5       1 
    ############################################################################
    ############################################################################
    #True    False   True    True    False   False   True    False   True meth..
    #True    True    True    True    None    True    False   False   True alig..
    #1       1       3       2       0       1       4       3       1    MeCnt.
    ############################################################################
    #66      20      66      66      27      56      57      28      93   CFU
    ############################################################################
    #39      -       39      -       -       -       32      -       -    cfua
    #-       7       -       -       -       24      7       83      -    cheC5
    #-       -       -       -       -       -       76      76      -    mycfi
    #-       -       68      -       -       -       87      29      -    ptr
    #-       -       28      28      -       -       77      77      -    sno
    #-       -       -       81      12      -       -       12      81   ss1
    ############################################################################
    ### Current annotated gene's orfid-structure
    ###CFU   gene-orfmodel: [66, 57, 28] CFU_838770
    ###cfua  gene-orfmodel: [39, 32] CFU_839553
    ###cheC5 gene-orfmodel: [24, 7, 83] CHEC5_080265
    ###mycfi gene-orfmodel: [76, 76, 111] MYCFI_71435
    ###ptr   gene-orfmodel: [68, 87, 29] PTRT_05805
    ###sno   gene-orfmodel: [28, 77] SNOT_15614
    ###ss1   gene-orfmodel: [81, 12, 41] SS1T_05007
    #
    # GeneModels are 100% certain correct on the 5' side (ClustalW confirmed)
    # This function will remove these nodes - blocks
    # ('CFU',20),('cheC5',7)     block 2
    # ('CFU',56),('cheC5',24)    block 6

    # return status variable
    is_any_block_removed = False

    for informant in PCG.organism_set():
        if informant == organism: continue
        pacbporfs = pacb.ordering.order_pacbporf_list(
            PCG.get_pacbps_by_organisms(organism,informant))
        looping = True
        while pacbporfs and looping:
            first = pacbporfs[0]
            if first.orfQ.has_methionine() and first.orfS.has_methionine():
                looping = False
                break
            nodeQ = (organism, first.orfQ.id)
            nodeS = (informant,first.orfS.id)
            # check in which block this PacbPORF is in
            for block in blocks:
                if len(block) == 1 and block[0][1:] == (nodeQ,nodeS):
                    # This one can be deleted!
                    key = block[0]
                    if not PCG.pacbps.has_key(key):
                        # already deleted in a previous block deletion!
                        continue
                    status = _delete_pacbp(PCG,key,PCG.pacbps[key])
                    if status:
                        is_any_block_removed = True
                        # pop first pacbporf from list
                        pacbporfs.pop(0)
                        # break out of the for loop -> skip EOF blocks reached
                        break
            else:
                # EOF blocks reached -> break out
                looping = False

    # return status weather or not a block was removed
    return is_any_block_removed

# end of function remove_unlikely_start_blocks


def remove_unlikely_orf_transition_blocks(blocks,PCG):
    """
    """
    #################################################################
    #2       1       3       3       3       2 
    #################################################################
    #################################################################
    #58      5       58      74      74      30      CFU
    #################################################################
    #40      40      40      39      73      73      cfua
    #-       -       16      1       1       89      mycfi
    #100     -       100     100     100     -       mycfia
    #################################################################
    #
    # Example (CFU_840834)
    # block 1-2-3 proses Orf transition 58-5-58,
    # with block 1 and 3 properly supported and block 2 poorly supported

    MAX_BLOCK_SIZE = 2

    is_any_block_removed = False
    # loop reversed over the blocks, IGNORE 1th and last block
    for pos in range(len(blocks)-2,0,-1):
        block = blocks[pos]
        # WARNING!! current implementation expects 1 here
        # (a block supported by a single PACBP)
        if len(block) > MAX_BLOCK_SIZE:
            continue

        prevblock = blocks[pos-1]
        nextblock = blocks[pos+1]
        block_edges     = [ (nodeQ,nodeS) for (pk,nodeQ,nodeS) in block ]
        prevblock_edges = [ (nodeQ,nodeS) for (pk,nodeQ,nodeS) in prevblock ]
        nextblock_edges = [ (nodeQ,nodeS) for (pk,nodeQ,nodeS) in nextblock ]
        bothblocks = []
        bothblocks.extend(prevblock)
        for elem in nextblock:
            if elem not in bothblocks:
                bothblocks.append( elem )
        bothblockBitScore = sum([ key[0] for (key,n1,n2) in bothblocks ])
        blockBitScore = sum([ key[0] for (key,n1,n2) in block ])

        criterionA = len(Set(prevblock_edges).intersection(nextblock_edges)) >= MAX_BLOCK_SIZE
        criterionB = len(Set(block_edges).intersection(nextblock_edges)) >= 1
        criterionC = len(Set(block_edges).intersection(prevblock_edges)) >= 1
        # TODO: criterionB and criterionC can be used as criteria to
        # disallow the deletion. TODO: implement this !?

        if len(Set(prevblock_edges).intersection(nextblock_edges)) >= 2 or\
        ( len(Set(prevblock_edges).intersection(nextblock_edges))==1 and\
        not Set(prevblock_edges).intersection(block_edges) and\
        not Set(nextblock_edges).intersection(block_edges) and\
        bothblockBitScore/blockBitScore > 20):
            # PREV and NEXT block are graphs from identical Orfs,
            # with this singleton block breaking this apart
            for pacbpkey in block:
                if not PCG.pacbps.has_key(pacbpkey):
                    # already deleted in a previous block deletion!
                    continue
                # delete this pacbporf
                pacbporf = PCG.pacbps[pacbpkey]
                _delete_pacbp(PCG,pacbpkey,pacbporf)

            # Update return value to True
            is_any_block_removed = True

    # return status weather or not a block was removed
    return is_any_block_removed

# end of function remove_unlikely_orf_transition_blocks



def _pacbplist2queryMSR_check_if_PCG_has_key(PCG,block):
    """ """
    _pacbps = []
    for k in block:
        if PCG.pacbps.has_key(k):
            _pacbps.append( PCG.pacbps[k] )
        else:
            # PacbP deletion in previous iteration
            pass
    # return pacbplist2queryMSR 
    return pacbplist2queryMSR(_pacbps)

# end of function _pacbplist2queryMSR_check_if_PCG_has_key


def remove_overlapping_blocks_with_conflicting_orfs(blocks,PCG):
    """
    """
    ############################################################################
    #454     481     496     522     536     586     649 
    #459     488     501     561     564     616     678 
    ############################################################################
    #1       3       2       1       5       5       1 
    ############################################################################
    ############################################################################
    #True    True    True    False   True    False   True    methionine?
    #True    True    True    None    False   False   True    alignedMs?
    #1       3       2       0       4       3       1       MethionineCnt
    ############################################################################
    #66      66      66      27      57      28      93      CFU
    ############################################################################
    #39      39      -       -       32      -       -       cfua
    #-       -       -       -       7       83      -       cheC5
    #-       -       -       -       76      76      -       mycfi
    #-       68      -       -       87      29      -       ptr
    #-       28      28      -       77      77      -       sno
    #-       -       81      12      -       12      81      ss1
    ############################################################################
    ## Current annotated gene's orfid-structure
    ##CFU   gene-orfmodel: [66, 57, 28] CFU_838770
    ##cfua  gene-orfmodel: [39, 32] CFU_839553
    ##cheC5 gene-orfmodel: [24, 7, 83] CHEC5_080265
    ##mycfi gene-orfmodel: [76, 76, 111] MYCFI_71435
    ##ptr   gene-orfmodel: [68, 87, 29] PTRT_05805
    ##sno   gene-orfmodel: [28, 77] SNOT_15614
    ##ss1   gene-orfmodel: [81, 12, 41] SS1T_05007
    ##
    #
    # GeneModels are 100% certain correct on the 5' side (ClustalW confirmed)
    # This function will remove these nodes - blocks
    # ('CFU',27),('ss1',12)      block 4
    #
    # The query MSR overlaps with block 5, and block 5 is better supported

    MAX_BLOCK_SIZE = 2
    MAX_AA_OVERLAP = 8

    is_any_block_removed = False
    # loop reversed over the blocks
    for pos in range(len(blocks)-1,-1,-1):
        block = blocks[pos]
        if len(block) > MAX_BLOCK_SIZE: continue
        prevblock = None
        nextblock = None
        if pos > 0:
            prevblock = blocks[pos-1]
            prevblockScore = sum([ key[0] for (key,n1,n2) in prevblock ])
        if pos < len(blocks)-1:
            nextblock = blocks[pos+1]
            nextblockScore = sum([ key[0] for (key,n1,n2) in nextblock ])

        blockMSR = _pacbplist2queryMSR_check_if_PCG_has_key(PCG,block)
        blockScore = sum([ key[0] for (key,n1,n2) in block ])
        if not blockMSR: continue
        if prevblock and blockScore < prevblockScore:
            prevMSR = _pacbplist2queryMSR_check_if_PCG_has_key(PCG,prevblock)
            if len(blockMSR.intersection(prevMSR)) > MAX_AA_OVERLAP:
                # delete this block!
                for pacbpkey in block:
                    if not PCG.pacbps.has_key(pacbpkey):
                        # already deleted in a previous block deletion!
                        continue
                    # delete this pacbporf
                    pacbporf = PCG.pacbps[pacbpkey]
                    _delete_pacbp(PCG,pacbpkey,pacbporf)
                # pop out the block from the list of blocks
                blocks.pop(pos)
                # Update return value to True
                is_any_block_removed = True
                # continue to next block in blocks
                continue
        if nextblock and blockScore < nextblockScore:
            nextMSR = _pacbplist2queryMSR_check_if_PCG_has_key(PCG,nextblock)
            if len(blockMSR.intersection(nextMSR)) > MAX_AA_OVERLAP:
                # delete this block!
                for pacbpkey in block:
                    if not PCG.pacbps.has_key(pacbpkey):
                        # already deleted in a previous block deletion!
                        continue
                    # delete this pacbporf
                    pacbporf = PCG.pacbps[pacbpkey]
                    _delete_pacbp(PCG,pacbpkey,pacbporf)
                # pop out the block from the list of blocks
                blocks.pop(pos)
                # Update return value to True
                is_any_block_removed = True
                # continue to next block in blocks
                continue

    # return status weather or not a block was removed
    return is_any_block_removed

# end of function remove_overlapping_blocks_with_conflicting_orfs



################################################################################
### asses the singleton PacbPORFs for mappable introns if not -> remove
################################################################################

def select_singleton_block_index(blocks,ORGANISMS):
    """ """
    indexes = []
    for pos in range(0,len(blocks)):
        block = blocks[pos]
        if len(block) == 1 and len(ORGANISMS) >= SINGLETON_BLOCK_MINIMAL_CNT or\
        float(len(block)) / len(ORGANISMS) <= SINGLETON_BLOCK_MINIMAL_RATIO:
            indexes.append(pos)
    # return list of positions
    return indexes

# end of function select_singleton_block_index


def select_singleton_block_keys(blocks,PCG,ORGANISMS):
    """ """
    singleton_pacbporfs = {}
    for pos in select_singleton_block_index(blocks,ORGANISMS):
        block = blocks[pos]
        for (blastkey,nodeQ,nodeS) in block:
            (orgQ,orfQ) = nodeQ
            (orgS,orfS) = nodeS
            if not PCG.has_edge(nodeQ,nodeS): continue
            try:
                pacbporfs = PCG.get_pacbps_by_nodes(node1=nodeQ,node2=nodeS)
            except:
                # no pacbporfs of this edge anymore in the PCG
                continue
            # if >1 returned -> which one !?!? TODO: solve this
            if len(pacbporfs) > 1: continue
            singletonpacbporf = pacbporfs[0]
            # store to singleton_pacbporfs dict
            if singleton_pacbporfs.has_key(orgS):
                singleton_pacbporfs[orgS].append( (blastkey,nodeQ,nodeS) )
            else:
                singleton_pacbporfs[orgS] = [ (blastkey,nodeQ,nodeS) ]
    # return dict of singleton_pacbporfs
    return singleton_pacbporfs

# end of function select_singleton_block_keys


def remove_singleton_blocks_without_intron_signal(blocks,PCG,ORGANISMS,OPTIONS):
    """ """
    # select singleton pacbporfs
    singleton_block_keys = select_singleton_block_keys(blocks,PCG,ORGANISMS)
    is_any_block_removed = False
    print "SBC::", singleton_block_keys
    for org in singleton_block_keys.keys():
        pacbporfs = pacb.ordering.order_pacbporf_list(PCG.get_pacbps_by_organisms(OPTIONS.target,org))
        for singleton_key in singleton_block_keys[org]:
            print "### singleton:", org, singleton_key, "PACBPORFS:", len(pacbporfs), [pf.bitscore for pf in pacbporfs]
            (blastkey,nodeQ,nodeS) = singleton_key
            # check if not removed in the previous block
            if not PCG.has_edge(nodeQ,nodeS):
                continue
            # in select_singleton_block_keys() function UNIQUE edges are
            # selected, so just take the first list element
            singleton = PCG.get_pacbps_by_nodes(node1=nodeQ,node2=nodeS)[0]
            SINGLETON_IS_INVALID = False
            for pos in range(1,len(pacbporfs)):
                prevPACBP = pacbporfs[pos-1]
                nextPACBP = pacbporfs[pos]
                if not (pacb.comparison.IsIdenticalPacbPORF(prevPACBP,singleton)\
                or pacb.comparison.IsIdenticalPacbPORF(nextPACBP,singleton)):
                    continue
                # if here, this is the interface of this singleton
                identQorf = prevPACBP.orfQ.id == nextPACBP.orfQ.id
                identSorf = prevPACBP.orfS.id == nextPACBP.orfS.id
                if identQorf and identSorf:
                    # no criterion to test invalidity here !?
                    pass
                elif identQorf:
                    introns = pacb.connecting.projecting.merge_pacbporfs_by_intron_in_sbjct(
                        prevPACBP,nextPACBP,
                        projected_intron_max_nt_offset=0,
                        projected_intron_max_aa_offset=0)
                    if not introns:
                        SINGLETON_IS_INVALID = True
                        break
                elif identSorf:
                    introns = pacb.connecting.projecting.merge_pacbporfs_by_intron_in_query(
                        prevPACBP,nextPACBP,
                        projected_intron_max_nt_offset=0,
                        projected_intron_max_aa_offset=0)
                    if not introns:
                        seqerror = pacb.connecting.sequenceerror.merge_pacbporf_with_sequenceerror_in_query(prevPACBP,nextPACBP)
                        if not seqerror:
                            SINGLETON_IS_INVALID = True
                            # do NOT break ; check the potential other interface too
                        else:
                            print "SEQUENCE-ERROR!!", org, seqerror, singleton, SINGLETON_IS_INVALID
                            SINGLETON_IS_INVALID = False
                            break
                else:
                    # no criterion to test invalidity here !?
                    pass

            # test invalidity. If so ->  delete!
            if SINGLETON_IS_INVALID:
                # A block is labeled for deletion.
                # Update return value to True
                is_any_block_removed = True
                # delete this singleton pacbporf
                _delete_pacbp(PCG,singleton_key,singleton)
                print "SINGLETON::", singleton, len(blocks),
                # pop out the singleton_key from the list of blocks
                for block in blocks:
                    if singleton_key in block:
                        block.remove(singleton_key)
                        break
                # remove empty blocks
                while [] in blocks: blocks.remove([])
                print len(blocks)

    # return status weather or not a block was removed
    return is_any_block_removed

# end of function remove_singleton_blocks_without_intron_signal



def remove_singleton_blocks_with_novel_orfs(blocks,PCG,ORGANISMS,OPTIONS,input):
    """ """
    # select singleton pacbporfs
    singleton_block_keys = select_singleton_block_keys(blocks,PCG,ORGANISMS)
    is_any_block_removed = False
    for org,pacbpkeys in singleton_block_keys.iteritems():
        for singleton_key in pacbpkeys:
            (pacbpkey,nodeQ,nodeS) = singleton_key
            (orgQ,orfQid) = nodeQ
            print "SingletonNovelOrf:", org, singleton_key,
            print orfQid in input[OPTIONS.target]['orfid-genestructure']
            if orfQid not in input[OPTIONS.target]['orfid-genestructure']:
                _delete_pacbp(PCG,singleton_key,None)
                # Update return value to True
                is_any_block_removed = True

    # return status weather or not a block was removed
    return is_any_block_removed

# end of function remove_singleton_blocks_with_novel_orfs



def remove_singleton_blocks(blocks,PCG,ORGANISMS,OPTIONS,input):
    """ """
    # select singleton pacbporfs
    singleton_block_keys = select_singleton_block_keys(blocks,PCG,ORGANISMS)
    is_any_block_removed = False
    for org,pacbpkeys in singleton_block_keys.iteritems():
        for singleton_key in pacbpkeys:
            (pacbpkey,nodeQ,nodeS) = singleton_key
            singleton_pacbporf = PCG.pacbps[singleton_key]
            _delete_pacbp(PCG,singleton_key,singleton_pacbporf)
            # Update return value to True
            is_any_block_removed = True

    # return status weather or not a block was removed
    return is_any_block_removed

# end of function remove_singleton_blocks

"""
Tinyexon detection function for ABFGP pairwise mode
"""

# Module metadata
__authors__ = "Ate van der Burgt"
__license__ = "MIT"

# Python Imports
from sets import Set

# PACB Imports
from pacb.pacbp import PacbP
from pacb.connecting.functions import _update_kwargs
from pacb.connecting.orfs import get_potential_tiny_exons_on_orf
from pacb.connecting.projecting import merge_pacbporfs_by_intron_in_query
from pacb.connecting.mapping import merge_pacbporfs_with_introns
from pacb.conversion import (
    exononorfs2pacbporf,
    pacbp2pacbporf,
    )
from pacb.ordering import order_pacbporf_list


# ABFGP Imports
from graphAbgp import InwardsPointingCodingBlockGraph
from pythonlibs.ordering import order_list_by_attribute
from lib_proteinsimilaritymatrix import (
    ProteinSimilarityMatrix,
    make_alignment_match,
    )


from lib_hmm_pairwise import is_hmmpacbporf_conflicting_with_pacbporflist

# Global Variable Imports
from settings.splicesites import KWARGS_TINYEXON_PAIRWISE
from settings.emboss import BLOSUM62_PATH
TINYEXON_MATRIX_PATH=BLOSUM62_PATH
TINYEXON_MATRIX = ProteinSimilarityMatrix(fname=TINYEXON_MATRIX_PATH)


def discover_tinyexons_simple(*args,**kwargs):
    """
    @attention: passing function for _discover_tinyexons_pairwise()
    @attention: see _discover_tinyexons_pairwise for *args and **kwargs
    """
    kwargs['final_discovery_ratio'] = 0.7
    return _discover_tinyexons_pairwise(['PP'],*args,**kwargs)

# end of function discover_tinyexons_simple



def discover_tinyexons_complex(*args,**kwargs):
    """
    @attention: passing function for _discover_tinyexons_pairwise()
    @attention: see _discover_tinyexons_pairwise for *args and **kwargs
    """
    kwargs['final_discovery_ratio'] = 0.6
    return _discover_tinyexons_pairwise(['PP','QP','QQ','QQug'],*args,**kwargs)

# end of function discover_tinyexons_complex


def _discover_tinyexons_pairwise(perform_steps,PCG,input,target,GENE_INFORMANT_SET,UNIGENE_INFORMANT_SET,
    min_discovery_count=2,final_discovery_ratio=0.6,verbose=False):
    """
    @attention: PCG *must* contain PacbPORFs, not PacbPs objects
    """

    # predict tiny exons in target and informants; omit unigene informants
    tinyexondata = _get_tinyexon_dict(input,
            omit_identifier_list=UNIGENE_INFORMANT_SET)

    ############################################################################
    if verbose:
        print "tinyexons:", KWARGS_TINYEXON_PAIRWISE
        print "tinyexons OMIT:", UNIGENE_INFORMANT_SET
        print "tinyexons:",
        print [ (k,len(tinyexondata[k])) for k in tinyexondata.keys() ]
        for tinyexon in tinyexondata[target]:
            print target, tinyexon, tinyexon.proteinsequence()
        for informant in tinyexondata.keys():
            if informant == target: continue
            for tinyexon in tinyexondata[informant]:
                print informant, tinyexon, tinyexon.proteinsequence()
            
    ############################################################################

    # PP step; 1 conserved tinyexon in both Query and Informant(s)
    if 'PP' in perform_steps:
        # find putatively conserved tiny exons in target and informants
        protmatches = _find_pp_tinyexons(target,tinyexondata,
                min_discovery_count=min_discovery_count)
    
        # get PP tinyexon PacbpORF data
        tinyexons = _convert_tinyexon_proteinmatches_to_pacbporfs(
                target,protmatches,tinyexondata,PCG,
                min_discovery_count=min_discovery_count)
    else:
        tinyexons = {}

    # QP step; 1 intron additionally gained in the Query / target gene
    if 'QP' in perform_steps:
        # get QP tinyexon PacbpORF data
        qp_tinyexons = _find_qp_and_pq_tinyexons_as_pacbporfs(target,tinyexondata,PCG,
                min_discovery_count=min_discovery_count)
    else:
        qp_tinyexons = {}

    # QQ step; 2 introns gained in the Query / target gene
    if 'QQ' in perform_steps:
        # get QQ tinyexon PacbpORF data
        qq_tinyexons = _find_qq_tinyexons_as_pacbporfs(target,tinyexondata,PCG,
                min_discovery_count=min_discovery_count)
    else:
        qq_tinyexons = {}

    # QQ step; 2 introns gained in the Query / target gene; find in UNIGENES
    if 'QQug' in perform_steps:
        ugtinyexondata = {}
        ugtinyexondata[target] = tinyexondata[target]
        for informant in UNIGENE_INFORMANT_SET:
            ugtinyexondata[informant] = [] # empty list
        # get QQ tinyexon PacbpORF data of unigenes
        qq_ug_tinyexons = _find_qq_tinyexons_as_pacbporfs(target,ugtinyexondata,PCG,
                min_discovery_count=1)
    else:
        qq_ug_tinyexons = {}

    # merge QP, QQ and QQ_ug tinyexons into main PP tinyexon PacbPORF data
    for datadict in [qp_tinyexons,qq_tinyexons,qq_ug_tinyexons]:
        for key,data in datadict.iteritems():
            if not tinyexons.has_key(key):
                tinyexons[key] = data
            else:
                # update by merging
                for (pacbporf,rejected,informant) in data:
                    _update_tinyexon_pacbporf_dict(tinyexons,
                        key,pacbporf,rejected,informant)


    ############################################################################
    if verbose:
        for k,data in tinyexons.iteritems():
            print k, len(data), [ b for a,b,c in data ],
            print [ c for a,b,c in data ],
            print [ a._tinyexon_label for a,b,c in data ],
            print [ a.bitscore for a,b,c in data ],
            print [ "%1.2f" % a.identityscore for a,b,c in data ],
            print "INIT"

    ############################################################################

    # cleanup tinyexon pacbporfs that have been observed only once
    _remove_dict_elements_with_short_value_list(
            tinyexons,min_value_list_size=2)

    # cleanup tinyexon pacbporfs that are only conflicting
    del_keys = []
    for key,data in tinyexons.iteritems():
        conflicts = [ b for a,b,c in data ]
        if conflicts.count(False) == 0:
            del_keys.append(key)
    # actually delete these keys
    for key in del_keys: del(tinyexons[key])


    # cleanup tinyexon pacbporfs for which, when decreased for conflicting
    # PacbPORFs, are observed to little
    del_keys = []
    for key,data in tinyexons.iteritems():
        conflicts = [ b for a,b,c in data ]
        if len(data) - conflicts.count(True) < min_discovery_count:
            del_keys.append(key)
    # actually delete these keys
    for key in del_keys: del(tinyexons[key])
    

    # cleanup tinyexon pacbporfs with low coverage ratio on GENE_INFORMANT_SET
    del_keys = []
    for key,data in tinyexons.iteritems():
        perfectppintrons = [ _has_pp_tinyexonpacbporf_perfect_introns(a,target,c,PCG) for a,b,c in data ]
        encountered_score = float( len(data) + perfectppintrons.count(True) )
        if encountered_score / len(GENE_INFORMANT_SET) >= final_discovery_ratio:
            pass
        elif ( len(data) - [a._tinyexon_label for a,b,c in data].count("PP")) >= 5:
            # nicely covered or at least a certain number of convincing tinyexons
            pass
        else:
            # mark for deletion!
            del_keys.append(key)
    # actually delete these keys
    for key in del_keys: del(tinyexons[key])

    ############################################################################
    if verbose:
        for k,data in tinyexons.iteritems():
            print k, len(data), [ b for a,b,c in data ],
            print [ c for a,b,c in data ],
            print [ a._tinyexon_label for a,b,c in data ],
            print [ a.bitscore for a,b,c in data ],
            print [ "%1.2f" % a.identityscore for a,b,c in data ]
    ############################################################################

    return tinyexons

# end of function _discover_tinyexons_pairwise


def _has_pp_tinyexonpacbporf_perfect_introns(tinyexonPF,target,informant,PCG):
    """ """
    # check if a (perfect) introns can be mapped
    is_confirmed_with_introns = False

    if tinyexonPF._tinyexon_label != 'PP': return False

    # get ordered PacbPORFS for this informant
    thepacbporfs = order_pacbporf_list(PCG.get_pacbps_by_organisms(target,informant))

    for pos in range(1,len(thepacbporfs)):
        prevPF,nextPF = thepacbporfs[pos-1],thepacbporfs[pos]
        if prevPF.distance_towards(tinyexonPF) > 0 and\
        tinyexonPF.distance_towards(nextPF) > 0:
            intronsPREV = merge_pacbporfs_with_introns(
                    prevPF,tinyexonPF,max_aa_offset=0,
                    max_intron_nt_length=None)
            intronsNEXT = merge_pacbporfs_with_introns(
                    tinyexonPF,nextPF,max_aa_offset=0,
                    max_intron_nt_length=None)
            if len(intronsPREV) >= 1 and len(intronsNEXT) >= 1:
                perfect_prev_intron = False
                perfect_next_intron = False
                for intronQ,intronS in intronsPREV:
                    intronQ.assign_bp_and_ppts()
                    intronS.assign_bp_and_ppts()
                    if intronQ.branchpoint and intronS.branchpoint:
                        perfect_prev_intron = True
                        break
                for intronQ,intronS in intronsNEXT:
                    intronQ.assign_bp_and_ppts()
                    intronS.assign_bp_and_ppts()
                    if intronQ.branchpoint and intronS.branchpoint:
                        perfect_next_intron = True
                        break
                # check if both intron options have a perfect candidate
                if perfect_prev_intron and perfect_next_intron:
                    is_confirmed_with_introns = True
            # break out
            break

    # return is_confirmed_with_introns status
    return is_confirmed_with_introns

# end of function _has_pp_tinyexonpacbporf_perfect_introns
                

def _convert_tinyexon_proteinmatches_to_pacbporfs(target,protmatches,
    tinyexondata,PCG,min_discovery_count=2):
    """  """
    target_tinyexon_pacbporf_data = {}
    # fish these protein matches from the tinyexons and convert to PacbPORFs
    for informant in tinyexondata.keys():
        if informant == target: continue
        thepacbporfs = order_pacbporf_list(
                PCG.get_pacbps_by_organisms(target,informant))
        for exonQ in tinyexondata[target]:
            if exonQ.orf.id in [ pf.orfQ.id for pf in thepacbporfs ]: continue
            if exonQ.proteinsequence() not in protmatches.keys(): continue
            for exonS in tinyexondata[informant]:
                if exonS.length > exonQ.length: break
                if exonS.proteinsequence() not in protmatches.keys(): continue
                # omit non-identical exons
                if not _are_tinyexons_similar(exonQ,exonS): continue
                # if here: similar exons. make PacbPORF
                pacbporf = exononorfs2pacbporf(exonQ,exonS,matrix=TINYEXON_MATRIX)
                if not pacbporf: continue
                # check if placeable in PCG/pacbporflist
                rejected = [ pf.is_postioned_compatibly(pacbporf) for pf in thepacbporfs ].count(False) > 0

                # label pacbporf as found by tinyexon PP
                pacbporf._tinyexon_label = "PP"

                # store to target_tinyexon_pacbporf_data
                key = (exonQ.proteinsequence(),exonQ.start)
                _update_tinyexon_pacbporf_dict(target_tinyexon_pacbporf_data,
                    key,pacbporf,rejected,informant)

    # cleanup tinyexon protein matches that have been observed to litte
    _remove_dict_elements_with_short_value_list(
            target_tinyexon_pacbporf_data,
            min_value_list_size=min_discovery_count)

    # return target_tinyexon_pacbporf_data
    return target_tinyexon_pacbporf_data

# end of function _convert_tinyexon_proteinmatches_to_pacbporfs


def _get_tinyexon_dict(input,omit_identifier_list=[],**kwargs):
    """ """
    _update_kwargs(kwargs,KWARGS_TINYEXON_PAIRWISE)
    tinyexondata = {}
    for orgid in input.keys():
        if orgid in omit_identifier_list: continue
        tinyexondata[orgid] = []
        for orfObj in input[orgid]['orfs'].orfs:
            tinyexondata[orgid].extend( get_potential_tiny_exons_on_orf(
                        orfObj,**kwargs ) )
        tinyexondata[orgid] = order_list_by_attribute(
                        tinyexondata[orgid],order_by='length')
    # return dict with predicted tinyexons
    return tinyexondata

# end of function _get_tinyexon_dict


def _are_tinyexons_identical(exonQ,exonS):
    """ """
    if exonS.length != exonQ.length: return False
    if exonQ.donor.phase != exonS.donor.phase: return False
    if exonQ.proteinsequence() != exonS.proteinsequence(): return False
    # if here: identical!
    return True

# end of function _are_tinyexons_identical


def _are_tinyexons_similar(exonQ,exonS,min_identity_count=1,max_dissimilar_count=0):
    """ """
    if exonS.length != exonQ.length: return False
    if exonQ.donor.phase != exonS.donor.phase: return False
    if exonQ.proteinsequence() == exonS.proteinsequence(): return True
    # if here: test similarity
    match = make_alignment_match(
            exonQ.proteinsequence(),
            exonS.proteinsequence(),
            matrix=TINYEXON_MATRIX.matrix )
    if match.count(" ") <= max_dissimilar_count and\
    match.count("*") >= min_identity_count:
        return True
    else:
        return False

# end of function _are_tinyexons_similar


def _find_qp_or_pq_match_on_orfobj(exon,orfObj,min_identity_count=1,max_dissimilar_count=0):
    """ """
    tinyexonmatches = []
    protseq = exon.proteinsequence()
    protlen = len(protseq)
    for offset in range(0,orfObj.protein_length-protlen):
        seqpart = orfObj.protein_sequence[offset:offset+protlen]
        match = make_alignment_match(protseq,seqpart,matrix=TINYEXON_MATRIX.matrix)
        if match.count(" ") <= max_dissimilar_count and\
        match.count("*") >= min_identity_count:
            aapos  = orfObj.protein_startPY + offset
            dnapos = orfObj.aapos2dnapos(aapos)
            if exon.acceptor.phase == 2: dnapos-=1
            if exon.acceptor.phase == 1: dnapos-=2
            dnaseq = orfObj.inputgenomicsequence[dnapos-2:dnapos+exon.length+2].upper()
            if dnaseq[0:2] == 'AG' or dnaseq[-2:] in ['GT','GC']:
                tinyexonmatches.append( (seqpart,aapos) )

    # return list of tinyexon match tuples
    return tinyexonmatches

# end of function _find_qp_or_pq_match_on_orfobj


def _find_match_on_orfobj(exon,orfObj,min_identity_count=1,max_dissimilar_count=0):
    """ """
    tinyexonmatches = []
    protseq = exon.proteinsequence()
    protlen = len(protseq)
    for offset in range(0,orfObj.protein_length-protlen):
        seqpart = orfObj.protein_sequence[offset:offset+protlen]
        match = make_alignment_match(protseq,seqpart,matrix=TINYEXON_MATRIX.matrix)

        if protseq == "SGWNAA" and seqpart in ['SGFNSA','SGWNAA','GLFNSV','SGFTSA','GGFTSA','GDFNAV','GKFNTI','SGFNSA','GNFTTI','GGGSTN','GDFSAV','GKFNTI','GAFTSA']:
            maxQ = TINYEXON_MATRIX.scorealignment(protseq,protseq)
            maxS = TINYEXON_MATRIX.scorealignment(seqpart,seqpart)
            print True, protseq, "'%s'" % match, seqpart, (TINYEXON_MATRIX.scorealignment(protseq,seqpart),maxQ,maxS), orfObj

        #if protseq == "LSPSM":
        #    maxQ = TINYEXON_MATRIX.scorealignment(protseq,protseq)
        #    maxS = TINYEXON_MATRIX.scorealignment(seqpart,seqpart)
        #    print False, protseq, "'%s'" % match, seqpart, (TINYEXON_MATRIX.scorealignment(protseq,seqpart),maxQ,maxS), orfObj


        if match.count(" ") <= max_dissimilar_count and\
        match.count("*") >= min_identity_count:
            aapos  = orfObj.protein_startPY + offset
            tinyexonmatches.append( (seqpart,aapos) )

        elif len(seqpart) >= 5 and match.count("*") >= min_identity_count and\
        TINYEXON_MATRIX.scorealignment(protseq,seqpart) > 0 and\
        match.count(" ") <= max_dissimilar_count+1:
            # escape for longer tinyexons; relax constrain a little bit
            aapos  = orfObj.protein_startPY + offset
            tinyexonmatches.append( (seqpart,aapos) )

        else:
            pass


    # return list of tinyexon match tuples
    return tinyexonmatches

# end of function _find_match_on_orfobj

def _find_pp_tinyexons(target,tinyexondata,min_discovery_count=2):
    """ """
    protmatches = {}
    for informant in tinyexondata.keys():
        if informant == target: continue
        for exonQ in tinyexondata[target]:
            for exonS in tinyexondata[informant]:
                if exonS.length > exonQ.length: break
                is_identical = _are_tinyexons_identical(exonQ,exonS)
                is_similar = False
                if not is_identical:
                    is_similar = _are_tinyexons_similar(exonQ,exonS)

                # abandon non-similar tinexon protein sequences
                if not is_identical and not is_similar: continue

                # if here, store this match
                if protmatches.has_key(exonQ.proteinsequence()):
                    protmatches[exonQ.proteinsequence()]+=1
                else:
                    protmatches[exonQ.proteinsequence()]=1
                if is_similar:
                    if protmatches.has_key(exonS.proteinsequence()):
                        protmatches[exonS.proteinsequence()]+=1
                    else:
                        protmatches[exonS.proteinsequence()]=1

    # cleanup tinyexon protein matches that have been observed to litte
    del_keys = []
    for seq,cnt in protmatches.iteritems():
        if cnt < min_discovery_count:
            del_keys.append(seq)
    # actual delete here
    for seq in del_keys: del( protmatches[seq] )

    # return protein matches
    return protmatches

# end of function _find_pp_tinyexons


def _find_qp_and_pq_tinyexons_as_pacbporfs(target,tinyexondata,PCG,min_discovery_count=2):
    """ """
    target_tinyexon_pacbporf_data = {}
    for informant in tinyexondata.keys():
        if informant == target: continue
        thepacbporfs = order_pacbporf_list(
                PCG.get_pacbps_by_organisms(target,informant))
        for exonQ in tinyexondata[target]:
            if exonQ.orf.id in [ pf.orfQ.id for pf in thepacbporfs ]: continue
            for orfObj in PCG.get_orfs_of_graph(organism=informant):
                tinyexonmatches = _find_qp_or_pq_match_on_orfobj(exonQ,orfObj)
                for (aaseq,aapos) in tinyexonmatches:
                    # make pacbporf object
                    pacbpobj = PacbP(input=(
                            exonQ.proteinsequence(), aaseq,
                            exonQ.orf.dnapos2aapos(exonQ.start), aapos ) )
                    pacbporfobj = pacbp2pacbporf(pacbpobj,exonQ.orf,orfObj)
                    pacbporfobj.extend_pacbporf_after_stops()
    
                    # remove included pacbporfs
                    is_suborsuperset = False
                    for accepted_pacbporf in thepacbporfs:
                        if pacbporfobj.issubsetorsuperset(accepted_pacbporf):
                            is_suborsuperset = True
                            break
                    if is_suborsuperset:
                        continue

                    # check if a (perfect) intron can be projected
                    is_confirmed_by_intron_projection = False
                    for accepted_pacbporf in thepacbporfs:
                        if accepted_pacbporf.orfS.id == pacbporfobj.orfS.id:
                            if min(accepted_pacbporf.alignment_dna_range_query()) > min(pacbporfobj.alignment_dna_range_query()):
                                try:
                                    introns = merge_pacbporfs_by_intron_in_query(
                                        pacbporfobj,accepted_pacbporf,
                                        max_aa_offset=0,
                                        max_intron_nt_length=None)
                                        #max_intron_nt_length=140)
                                except IndexError:
                                    # unexpected event: TODO: solve in merge_pacbporfs_by_intron_in_query
                                    introns = []

                            else:
                                try:
                                    introns = merge_pacbporfs_by_intron_in_query(
                                        accepted_pacbporf,pacbporfobj,
                                        max_aa_offset=0,
                                        max_intron_nt_length=None)
                                        #max_intron_nt_length=140)
                                except IndexError:
                                    # unexpected event: TODO: solve in merge_pacbporfs_by_intron_in_query
                                    introns = []

                            if len(introns) >= 1:
                                is_confirmed_by_intron_projection = True
                                break

                    # continue if not is_confirmed_by_intron_projection
                    if not is_confirmed_by_intron_projection: continue

                    # check if placeable in PCG/pacbporflist
                    rejected = [ pf.is_postioned_compatibly(pacbporfobj) for pf in thepacbporfs ].count(False) > 0

                    # label pacbporf as found by tinyexon QP
                    pacbporfobj._tinyexon_label = "QP"

                    # store to target_tinyexon_pacbporf_data
                    key = (exonQ.proteinsequence(),exonQ.start)
                    _update_tinyexon_pacbporf_dict(
                            target_tinyexon_pacbporf_data,
                            key,pacbporfobj,rejected,informant)


    # cleanup tinyexon protein matches that have been observed to litte
    _remove_dict_elements_with_short_value_list(
            target_tinyexon_pacbporf_data,
            min_value_list_size=min_discovery_count)

    # return target_tinyexon_pacbporf_data
    return target_tinyexon_pacbporf_data

# end of function _find_qp_and_pq_tinyexons_as_pacbporfs



def _find_qq_tinyexons_as_pacbporfs(target,tinyexondata,PCG,min_discovery_count=2):
    """ """
    target_tinyexon_pacbporf_data = {}
    for informant in tinyexondata.keys():
        if informant == target: continue
        thepacbporfs = order_pacbporf_list(
                PCG.get_pacbps_by_organisms(target,informant))
        for exonQ in tinyexondata[target]:
            if exonQ.orf.id in [ pf.orfQ.id for pf in thepacbporfs ]: continue
            for (prevpos,nextpos) in [ (pos-1,pos) for pos in range(1,len(thepacbporfs)) ]:
                prevPF = thepacbporfs[prevpos]
                nextPF = thepacbporfs[nextpos]
                if prevPF.orfS.id == nextPF.orfS.id:

                    # check if PacbPORFs are positioned more or less okay
                    if prevPF.distance_towards(nextPF) > 20: continue

                    # check if exonQ is positioned ~between these PacbPORFs
                    if exonQ.orf.dnapos2aapos(exonQ.end) < max(prevPF.alignment_protein_range_query())-12:
                        continue
                    if exonQ.orf.dnapos2aapos(exonQ.start) > min(nextPF.alignment_protein_range_query())+12:
                        continue

                    # check if gap can be projected already by a perfect intron
                    introns = merge_pacbporfs_by_intron_in_query(
                                prevPF,nextPF,max_aa_offset=1)
                    # if introns found => continue
                    if introns: continue

                    # orfObj is the orfS of prevPF or nextPF (just take any)
                    orfObj = prevPF.orfS
                    # assign elegiable range of tinyexon match on SBJCT
                    aapos_sbjct_range = range(
                            max(prevPF.alignment_protein_range_sbjct())-12,
                            min(nextPF.alignment_protein_range_sbjct())+12
                            )

                    tinyexonmatches = _find_match_on_orfobj(exonQ,orfObj)
                    for (aaseq,aapos) in tinyexonmatches:
                        # check if the match is obtained in the expected
                        # sbjct AA range; if not, ignore the match
                        if aapos not in aapos_sbjct_range: continue

                        # make pacbporf object
                        pacbpobj = PacbP(input=(
                                exonQ.proteinsequence(), aaseq,
                                exonQ.orf.dnapos2aapos(exonQ.start), aapos ) )
                        pacbporfobj = pacbp2pacbporf(pacbpobj,exonQ.orf,orfObj)
                        pacbporfobj.extend_pacbporf_after_stops()
        
                        # remove included pacbporfs
                        is_suborsuperset = False
                        for accepted_pacbporf in thepacbporfs:
                            if pacbporfobj.issubsetorsuperset(accepted_pacbporf):
                                is_suborsuperset = True
                                break
                        if is_suborsuperset:
                            continue
    

                        # check if 2 (perfect) introns can be projected
                        introns5p = merge_pacbporfs_by_intron_in_query(
                                prevPF,pacbporfobj,
                                max_aa_offset=1,
                                max_intron_nt_length=None)
                                #max_intron_nt_length=140)
                        introns3p = merge_pacbporfs_by_intron_in_query(
                                pacbporfobj,nextPF,
                                max_aa_offset=1,
                                max_intron_nt_length=None)
                                #max_intron_nt_length=140)

                        # continue if not is_confirmed_by_intron_projection
                        if not introns5p or not introns3p: continue
    
                        # check if placeable in PCG/pacbporflist
                        distPrev = prevPF.distance_towards(pacbporfobj)
                        distNext = pacbporfobj.distance_towards(nextPF)
                        ovrlPrev = pacbporfobj.overlap(prevPF)
                        ovrlNext = pacbporfobj.overlap(nextPF)
                        if distPrev and distNext:
                            rejected = False
                        elif not distPrev and ovrlPrev:
                            rejected = False
                        elif not distNext and ovrlNext:
                            rejected = False
                        elif ovrlPrev and ovrlNext:
                            rejected = False
                        else:
                            rejected = True

                        print "OKAY", exonQ.proteinsequence(), aaseq, rejected, informant, (distPrev,distNext,ovrlPrev,ovrlNext)

                        # label pacbporf as found by tinyexon QQ
                        pacbporfobj._tinyexon_label = "QQ"

                        # store to target_tinyexon_pacbporf_data
                        key = (exonQ.proteinsequence(),exonQ.start)
                        _update_tinyexon_pacbporf_dict(
                                target_tinyexon_pacbporf_data,
                                key,pacbporfobj,rejected,informant)


    # cleanup tinyexon protein matches that have been observed to litte
    _remove_dict_elements_with_short_value_list(
            target_tinyexon_pacbporf_data,
            min_value_list_size=min_discovery_count)

    # return target_tinyexon_pacbporf_data
    return target_tinyexon_pacbporf_data

# end of function _find_qq_tinyexons_as_pacbporfs



def _remove_dict_elements_with_short_value_list(d,min_value_list_size=2):
    """ """
    del_keys = []
    for key,data in d.iteritems():
        if len(data) < min_value_list_size:
            del_keys.append(key)
    # actual delete here
    for key in del_keys: del( d[key] )

# end of function _remove_dict_elements_with_short_value_list


def _update_tinyexon_pacbporf_dict(tedict,key,pacbporf,rejected,informant):
    """ """
    if tedict.has_key(key) and informant in [ c for a,b,c in tedict[key] ]:
        for pos in range(0,len(tedict[key])):
            _pacbporf,_rejected,_informant = tedict[key][pos]
            if _informant == informant:
                if pacbporf.bitscore > _pacbporf.bitscore or (
                pacbporf.bitscore == _pacbporf.bitscore and\
                rejected == False and _rejected == True):
                    # replace; this pacbporf is more logical than current one
                    tedict[key][pos] = (pacbporf,rejected,informant)
                else:
                    # do not replace
                    pass
                # break out forloop
                break
    elif tedict.has_key(key):
        # append to tedict[key]
        tedict[key].append(
            (pacbporf,rejected,informant) )
    else:
        # create new element in tedict
        tedict[key] =\
            [ (pacbporf,rejected,informant) ]

# end of function _update_tinyexon_pacbporf_dict



def tinyexonpacbporfs2InwardsPointingCodingBlockGraphs(target,tinyexonpacbporfs):
    """ translate tinyexonpacbporfs to tinyexonInwpCBGs """
    tinyexonInwpCBGs = []
    for key,datalist in tinyexonpacbporfs.iteritems():
        inwpCBG = InwardsPointingCodingBlockGraph()
        (pacbporf,rejected,informant) = datalist[0]
        nodeQ = (target,pacbporf.orfQ.id)
        inwpCBG.add_node(nodeQ)
        for (pacbporf,rejected,informant) in datalist:
            nodeS = (informant,pacbporf.orfS.id)
            inwpCBG.add_node(nodeS)
            inwpCBG.add_edge(nodeQ,nodeS,wt=pacbporf.bitscore)
            pacbpkey = pacbporf.construct_unique_key(nodeQ,nodeS)
            inwpCBG.pacbps[(pacbpkey,nodeQ,nodeS)] = pacbporf
    
        tinyexonInwpCBGs.append( ( inwpCBG.get_bitscore(), inwpCBG ) )

    # order by bitscore
    tinyexonInwpCBGs.sort()
    tinyexonInwpCBGs.reverse()
    tinyexonInwpCBGs = [ inwp for bitscore,inwp in tinyexonInwpCBGs ]
    # return bitscore ordered list of tinyexonInwpCBGs
    return tinyexonInwpCBGs

# end of function tinyexonpacbporfs2InwardsPointingCodingBlockGraphs

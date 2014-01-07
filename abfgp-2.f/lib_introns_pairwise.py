"""
Functions for mapped/projected/etc introns in Pairwise AB(F)GP
"""

# Module metadata
__authors__ = "Ate van der Burgt"
__license__ = "MIT"


# ABGP Imports
from lib_codingarray import (
    PCG2codingarray,
    PCG2similarityarray,
    )
from lib_abgpgff import get_raw_abfgp_genestruture
from settings.dbwarehouse import _get_organism_full_name
from settings.genestructure import OPTIMAL_BRACNHPOINT_TO_ACCEPTOR_DISTANCE

from gene.intron.comparison import (
    _branchpoint_comparison,
    _polypyrimidinetract_comparison,
    _algsimilarity_comparison,
    _finetune_splicesite_comparison,
    )

# Python Imports
from sets import Set


def _calc_aligned_score(introns):
    """ """
    fsource2score = {
        'ABGPprojectingTE': 2.0,
        'ABGPprojecting':   2.0,
        'ABGPmapping':      1.0,
        'ABGPcig':          1.0,
        'ABGPphs':          1.0,
        'ABGPbridgeing':    0.5,
        'ABGPstl3n':        1.0,
        'ABGPstl3nCONS':    1.0,
        }
    score = 0.0
    for intron in introns:
        if intron._gff.has_key('fsource') and fsource2score.has_key(intron._gff['fsource']):
            max_score = fsource2score[intron._gff['fsource']]
        else:
            max_score = 1.0
        if hasattr(intron,"_distance"):
            dist_ratio = 1.0 / ( 1+( intron._distance/3 ) )
        else:
            dist_ratio = 0.0
        score += (max_score * dist_ratio)
    return score

# end of function _calc_aligned_score


def _get_distance(intron):
    """ """
    for attr in ["_distance","distance"]:
        if hasattr(intron,attr):
            return abs(getattr(intron,attr))
# end of function _get_distance


def _filter_identical_interface_introns_on_best_pssm(intronsubdata,verbose=False):
    """ """
    # if nothing passed -> EOF function
    if not intronsubdata: return {}

    # (orderable) list to store the scores & associated intron key into
    score2intronkeys = []

    for intronkey in intronsubdata.keys():
        intronexmpl = intronsubdata[intronkey][0]
        if intronexmpl.__class__.__name__ == "SequenceErrorConnectingOrfs":
            # this effectively removes the SequenceError here!
            pssmscore = 0.0
            overallscore = pssmscore * intronexmpl._apps_donor * intronexmpl._apps_accep
            score2intronkeys.append( ( overallscore, intronkey ) )
        else:
            #pssmscore = intronexmpl.donor.pssm_score + intronexmpl.acceptor.pssm_score
            #overallscore = pssmscore * intronexmpl._apps_donor * intronexmpl._apps_accep
            #score2intronkeys.append( ( overallscore, intronkey ) )
            pssmscore = intronexmpl.donor.pssm_score + intronexmpl.acceptor.pssm_score
            apps_donors = [ intron._apps_donor for intron in intronsubdata[intronkey] ]
            apps_acceps = [ intron._apps_accep for intron in intronsubdata[intronkey] ]
            overallscore = pssmscore * (sum(apps_donors)/len(apps_donors)) *\
                           (sum(apps_acceps)/len(apps_acceps))
            score2intronkeys.append( ( overallscore, intronkey ) )

    # now return only the highest scoring intron based on PSSM/entropy
    score2intronkeys.sort()
    score2intronkeys.reverse()
    best_intron_key = score2intronkeys[0][1]
    # return subdict of only this intron
    return { best_intron_key: intronsubdata[ best_intron_key ] }

# end of function _filter_identical_interface_introns_on_best_pssm


def _filter_identical_interface_introns_on_best_msr(intronsubdata,
    abfgp_raw_exons,PCG,input,OPTIONS,verbose=False):
    """ """
    # if nothing passed -> EOF function
    if not intronsubdata: return {}

    ############################################################################
    # TODO: this must be moved to a single position in the code!!
    ############################################################################
    # CREATE alignment_presence/alignment_similarity tracks
    dna2protlength      = len(input[OPTIONS.target]['genomeseq'])/3
    array_algpresence   = PCG2codingarray(PCG,OPTIONS.target,dna2protlength)
    array_algsimilarity = PCG2similarityarray(PCG,OPTIONS.target,dna2protlength)
    ############################################################################
    # TODO: this must be moved to a single position in the code!!
    ############################################################################

    # check on which interface we are going to filter
    labels = []
    for intronlist in intronsubdata.values():
        for i in intronlist:
            labels.append( i._label )
    interface = _get_main_interface(labels)

    prevexon,nextexon = None,None
    for (_prev,_next) in [ (_pos-1,_pos) for _pos in range(1,len(abfgp_raw_exons)) ]:
        prevexon = abfgp_raw_exons[_prev]
        nextexon = abfgp_raw_exons[_next]
        if (prevexon.orf.id,nextexon.orf.id) == interface:
            break
        elif prevexon.orf.id == interface[1]:
            # first exon NOT in abfgp_raw_exons
            # this can happen when final exon is assigned by special function ->
            # intron in introndata, pacbporf in intron._linked_to_pacbporfs
            prevexon,nextexon = None,prevexon
            break           
        elif nextexon.orf.id == interface[0]:
            # final exon NOT in abfgp_raw_exons
            # this can happen when final exon is assigned by special function ->
            # intron in introndata, pacbporf in intron._linked_to_pacbporfs
            prevexon,nextexon = nextexon, None
            break
        else:
            pass
    else:
        if len(abfgp_raw_exons) == 1:
            if abfgp_raw_exons[0].orf.id == interface[1]:
                prevexon,nextexon = None,abfgp_raw_exons[0]
            elif abfgp_raw_exons[0].orf.id == interface[0]:
                prevexon,nextexon = abfgp_raw_exons[0],None
            else:
                if verbose: print "INTERFACE not present:", interface
                # interface not present in (RAW!) genestructure
                # function is not applicable any further
                return intronsubdata
        else:
            if verbose: print "INTERFACE not present:", interface
            # interface not present in (RAW!) genestructure
            # function is not applicable any further
            return intronsubdata

    # (orderable) list to store the scores & associated intron key into
    score2intronkeys = []

    keys = intronsubdata.keys()
    keys.sort()
    for intronkey in keys:
        # create Intron signal scores
        (intron_nt_start,intron_nt_stop) = intronkey
        intron_aa_start = intron_nt_start/3
        intron_aa_stop  = intron_nt_stop/3
        intron_length   = intron_aa_stop-intron_aa_start
        if intron_nt_stop-intron_nt_start in [2,3,4]:
            # avoid ZeroDivisionError; a SequenceErrorConnectingOrfs!
            algprs_score = 0
            intron_cnt   = len(intronsubdata[intronkey])
            intron_score = 0
            intron_ratio = 0.0
        else:
            algprs_score = sum(array_algpresence[intron_aa_start:intron_aa_stop+1])
            intron_cnt   = len(intronsubdata[intronkey])
            intron_score = (intron_length)*intron_cnt
            intron_ratio = float(algprs_score)/intron_score

        if prevexon:
            # calculate prev EXON signal score
            sta = prevexon.protein_start()
            end = intron_aa_start+1
            len5p  = end-sta-1
            # avoid ZeroDivisionError
            if len5p == 0: end,len5p = end+1,len5p+1
            algprs_prevexon_score = sum(array_algpresence[sta:end])
            exon5p_ratio = float(algprs_prevexon_score) / len5p
        else:
            algprs_prevexon_score = 0
            exon5p_ratio = 0.0
            len5p = 0

        if nextexon:
            # calculate next EXON signal score
            sta = intron_aa_stop
            end = nextexon.protein_end()+1
            len3p  = end-sta-1
            # avoid ZeroDivisionError
            if len3p == 0: end,len3p = end+1,len3p+1
            algprs_nextexon_score = sum(array_algpresence[sta:end])
            exon3p_ratio = float(algprs_nextexon_score) / len3p
        else:
            algprs_nextexon_score = 0
            exon3p_ratio = 0.0
            len3p = 0

        # avoid ZeroDivisionError; arbitrarily low intron_ratio
        # Do not set much lower; score ratio will become very high!
        if intron_ratio == 0.0: intron_ratio = 0.5

        score = (exon5p_ratio+exon3p_ratio)/intron_ratio
        score2intronkeys.append( (score, intronkey) )

        ########################################################################
        if verbose:
            example = intronsubdata[intronkey][0]
            print ">>SCORE:", score, intronkey,
            print (intron_score, algprs_score, intron_ratio),
            try:
                print "pssm: %1.2f - %1.2f" % (example.donor.pssm_score,example.acceptor.pssm_score),
            except:
                # SequenceError..
                print "None - None",
            print "5pExon: (%s,%s,%1.2f)" % ( algprs_prevexon_score, len5p, exon5p_ratio/intron_ratio ),
            print "3pExon: (%s,%s,%1.2f)" % ( algprs_nextexon_score, len3p, exon3p_ratio/intron_ratio ),
            print score
        ########################################################################

    # now return only the highest scoring intron based on AlgPres/AlgSim data
    score2intronkeys.sort()
    score2intronkeys.reverse()
    best_intron_key = score2intronkeys[0][1]
    # calculate relative signal strength by score division
    if score2intronkeys[-1][0] == 0.0:
        # avoid ZeroDivisionError
        algsim_signal_ratio = score2intronkeys[0][0] / 1.0
    else:
        algsim_signal_ratio = score2intronkeys[0][0] / score2intronkeys[-1][0]

    # return subdict of only this intron
    ###return { best_intron_key: intronsubdata[ best_intron_key ] }
    return { best_intron_key: algsim_signal_ratio }


# end of function _filter_identical_interface_introns_on_best_msr





def _filter_introns_for_exon_msr(introndata,input,PCG,OPTIONS):
    """ """

    ############################################################################
    # TODO: this must be moved to a single position in the code!!
    ############################################################################
    # CREATE alignment_presence/alignment_similarity tracks
    dna2protlength      = len(input[OPTIONS.target]['genomeseq'])/3
    array_algpresence   = PCG2codingarray(PCG,OPTIONS.target,dna2protlength)
    array_algsimilarity = PCG2similarityarray(PCG,OPTIONS.target,dna2protlength)
    ############################################################################
    # TODO: this must be moved to a single position in the code!!
    ############################################################################

    keys = introndata.keys()
    keys.sort()

    # what is the largest stack of introns delivered by informants?
    max_intron_coverage = float(max([ len(introndata[key]) for key in keys ]))
    poor_intron_coverage_ratio = 0.25

    for key in keys:

        if key not in introndata.keys():
            # intron cascade deleted...
            continue

        # create Intron and Exon signal scores
        (intron_nt_start,intron_nt_stop) = key
        intron_aa_start = intron_nt_start/3
        intron_aa_stop  = intron_nt_stop/3
        intron_coverage = len(introndata[key])
        algprs_score    = sum(array_algpresence[intron_aa_start:intron_aa_stop+1])
        algsim_score    = sum(array_algsimilarity[intron_aa_start:intron_aa_stop+1])
        algprs_zeros    = [ p==0 for p in array_algpresence[intron_aa_start:intron_aa_stop+1] ].count(True)
        intron_cnt      = len(introndata[key])
        intron_expl     = introndata[key][0]
        intron_score    = (intron_aa_stop-intron_aa_start)*intron_cnt

        print "@@@@", key, intron_score, algprs_score, algsim_score,
        print (intron_coverage,max_intron_coverage), "orfs::",
        print (introndata[key][0].orfDonor.id,introndata[key][0].orfAcceptor.id)

        delete_intron = None

        # first main check for very shitty stopless3n introns
        if intron_coverage/max_intron_coverage <= poor_intron_coverage_ratio and\
        intron_expl.is_stopless_3n_intron():
            intron_expl.assign_bp_and_ppts()
            if not intron_expl.branchpoint:
                delete_intron = True
            elif intron_expl.get_branchpoint_nt_distance() >\
            max(OPTIMAL_BRACNHPOINT_TO_ACCEPTOR_DISTANCE) + 5:
                # even with an additional extra offset of 5, still
                # a very unlikely branchpoint for a stopless3n intron
                delete_intron = True
            else:
                delete_intron = False

        # now do more elaborate check if this intron must be deleted
        if delete_intron:
            pass
        elif intron_score > algprs_score:
            # perfect intron signal, leave intact
            delete_intron = False
        elif intron_score > algsim_score:
            if intron_expl.is_stopless_3n_intron():
                # check intron evidence: is there an Orf transition in informants
                check = []
                for actualintron in introndata[key]:
                    print "@@@@", actualintron, actualintron._gff['fsource'], actualintron.is_tcode_coding(), actualintron.is_tcode_noncoding()
                    if actualintron._gff['fsource'] == "ABGPprojecting" and\
                    actualintron.is_tcode_noncoding() == False and\
                    actualintron._informant_label[0] ==\
                    actualintron._informant_label[1]:
                        # intron is obtained by projection on continious Orf
                        # in the informant; this very likely represent a False
                        # intron call in case of a poor similarity stretch
                        # in this protein alignment
                        check.append( False )
                    else:
                        check.append( True )
                if check.count(True) == 0:
                    delete_intron = True
                else:
                    # still accept as a valid intron
                    delete_intron = False

            else:
                # still accept as a valid intron
                delete_intron = False
        elif intron_expl.__class__.__name__ == 'SequenceErrorConnectingOrfs':
            # leave as it is; will be filtered upon lateron
            delete_intron = False
        elif intron_expl.is_stopless_3n_intron():
            # discard this intron!
            delete_intron = True
        else:
            if intron_score*4 < algprs_score and intron_cnt <= 2:
                # very poor intron, presence score >> intron_score
                delete_intron = True
            elif intron_score*2 < algprs_score and intron_cnt <= 2 and not intron_expl.branchpoint:
                delete_intron = True
            else:
                # although the exon signal is stronger as the intron signal,
                # there *must* be an intron here. So, leave intact
                delete_intron = False

        ########################################################################
        summedscore = _calc_aligned_score(introndata[key])
        interface   = _get_main_interface([i._label for i in introndata[key]])
        print "delete:", delete_intron, " %1.2f" % summedscore,
        print interface, key, intron_cnt,
        print intron_expl.is_stopless_3n_intron(),
        print intron_score,(algprs_score,algsim_score,algprs_zeros),
        try:
            print "%1.2f - %1.2f" % (intron_expl.donor.pssm_score,intron_expl.acceptor.pssm_score)
        except:
            print "None - None"
        ########################################################################

        # continue if intron is labeled to be retained
        if not delete_intron: continue

        # if here, this intron(s) can be deleted
        # when it represents a tinyexon interface, cascade removal!
        for del_intron in introndata[key]:
            _cascade_delete_intron(del_intron,introndata,forced_delete=True)

        # delete this intron            
        del( introndata[key] )

    # return filtered intron data
    return introndata

# end of function _filter_introns_for_exon_msr


def _filter_introns(introndata,input,PCG,ORGANISMS,OPTIONS,verbose=False):
    """ """
    # calculate score per intron position and make 2 data structures:
    # 1) dict key2score which links intron keys to ontained score
    # 2) list with tuples (score,intronkey) which will be ordered on score
    keyscores = []
    key2score = {}
    for key in introndata.keys():
        score = _calc_aligned_score(introndata[key])
        keyscores.append( ( score, key ) )
        key2score[key] = score
    keyscores.sort()
    keyscores.reverse()
    # now loop over each of the introns, highest scoring one first
    # inside this loop, loop again over each of the introns and
    # remove all introns that are obtained on the same interface
    # BUT have a lower overall score

    encountered_interfaces = []
    
    # threshold for different threatment of introns
    # >  SCORE_THRESHOLD: introns with (some) confidence
    # <= SCORE_THRESHOLD: singletons, poor aligned introns etc
    if len(ORGANISMS)-1 >= 8:
        # 8 or more informants -> increase threshold a little bit here
        SCORE_THRESHOLD = 3.0
    else:
        SCORE_THRESHOLD = 2.0


    ############################################################################
    # TODO: this must be moved to a single position in the code!!
    ############################################################################

    # data about the annotated gene structure
    geneobjects = input[OPTIONS.target]['gldobj'].as_exons()
    known_introns = [ geneobjects[pos] for pos in range(1,len(geneobjects),2) ]
    known_exons   = [ geneobjects[pos] for pos in range(0,len(geneobjects),2) ]
    known_orf_ids = [ exon.orf.id for exon in known_exons ] 
    abfgp_raw_exons = get_raw_abfgp_genestruture(OPTIONS.target,PCG,
                        introndata,known_exons,verbose=False)

    if verbose:
        for exon in abfgp_raw_exons: print "RAW EXON::", exon

    ############################################################################
    # TODO: this must be moved to a single position in the code!!
    ############################################################################


    for score,key in keyscores:
        if not introndata.has_key(key):
            # intron has been filtered out in the meantime
            continue

        # make unique list of encountered interfaces
        intfA = _get_main_interface([i._label for i in introndata[key]])
        if intfA not in encountered_interfaces:
            encountered_interfaces.append( intfA )

        if score <= SCORE_THRESHOLD:
            # no not process the `garbage` with this function
            continue

        # if here, this intron will be most likely maintained
        # store its interface to the list of encountered interfaces

        for _key,_score in key2score.iteritems():
            if _key == key:
                # this intron itself -> omit
                continue
            if not introndata.has_key(key):
                # HIGHER scoring intron has been filtered out in the meantime
                # this can happen in case MSR-PSSM overrule intron score filter
                continue
            if not introndata.has_key(_key):
                # intron has been filtered out in the meantime
                continue
            if  _score >= score:
                # higher score -> can not be deleted!
                continue


            if score > SCORE_THRESHOLD:
                pass
            elif score <= SCORE_THRESHOLD and (len(introndata[key]) >= 3 and len(introndata[_key]) == 1):
                # allow this secondary intron to be deleted
                pass
            elif score <= SCORE_THRESHOLD:
                key_distances  = Set([intron._distance for intron in introndata[key] ])
                _key_distances = Set([intron._distance for intron in introndata[_key] ])
                key_sources    = Set([intron._gff['fsource'] for intron in introndata[key] ])
                _key_sources   = Set([intron._gff['fsource'] for intron in introndata[_key] ])
                if list(_key_sources) == ['ABGPmapping'] and 0 not in _key_distances:
                    TODO = True

            else:
                # do NOT allow this secondary intron to be deleted
                continue

            # now check on which interface (orfD,orfA) the intron was obtained
            # intfA is defined above
            intfB = _get_main_interface( [ i._label for i in introndata[_key] ] )
            if intfA[0] != intfB[0] and intfA[1] != intfB[1]:
                # not obtained on the same interface -> continue
                continue


            # if interfaces are not (Exact) identical, check for
            # overlap between the predicted introns. In case they do not
            # overlap at all, check if they are the result of tinyexon efforts.
            # If so, introns can be `protected` here from deletion
            if not Set(range(key[0],key[1])).intersection(range(_key[0],_key[1])):
                # if here, then both introns do NOT overlap at all
                if intfA != intfB:
                    # not on the overall indentical interface
                    continue
                elif intfA == intfB and intfA[0] == intfA[1]:
                    # stopless3n introns which do not overlap at all
                    # 2 stopless3n introns in the same Orf!?
                    # unlikely, but yes. MYCGR_105888
                    continue
                else:
                    pass

                secondary_intron_is_protected = False
                for intron in introndata[key]:
                    if intron._gff['fsource'] in ['ABGPprojectingTE','ABGPmappingTE']:
                        if hasattr(intron,'_linked_to_introns'):
                           for cdi in intron._linked_to_introns:
                               protected_key = ( cdi.donor.pos, cdi.acceptor.pos )
                               if protected_key == _key:
                                   secondary_intron_is_protected = True
                                   #############################################
                                   if verbose: print protected_key, "PROTECTED BY:", key
                                   #############################################
                                   break
                # check if secondary_intron_is_protected == True
                # continue (and not delete this intron)
                if secondary_intron_is_protected == True:
                    continue
                else:
                    ############################################################
                    if verbose:
                        print "NOT OVERLAPPING REMOVE: =>", "main",
                        print key, "second.", _key, intfA, intfB
                    ############################################################
                    pass

            else:
                # Overlapping introns. When both introns have a score higher
                # then SCORE_THRESHOLD, check filter on best MSR & best pssm.
                # When the (marginally) lower scoring intron is appointed as
                # **the** intron by both MSR & best pssm analyses, overrrule
                # the score criterion and delete the higher scoring intron
                # (based on mapping/projecting distance)

                intronObjmax = introndata[key][0]
                intronObjsub = introndata[_key][0]


                if _score > SCORE_THRESHOLD or _score/score >= 0.70:
                    # compare these 2 introns detailedly
                    intronObjmax = introndata[key][0]
                    intronObjsub = introndata[_key][0]

                    # calculate length differences on donor/acceptor position
                    dpos_diff=abs(intronObjmax.donor.pos-intronObjsub.donor.pos)
                    apos_diff=abs(intronObjmax.acceptor.pos-intronObjsub.acceptor.pos)
                    intron_length_difference = dpos_diff + apos_diff

                    # compare these 2 introns together with any other intron(s)
                    # predicted on the Orf interface
                    tmp_dict = {
                            key:  introndata[key],
                            _key: introndata[_key] }
                    best_intron_msr =\
                            _filter_identical_interface_introns_on_best_msr(
                            tmp_dict,abfgp_raw_exons,PCG,input,OPTIONS)
                    best_intron_pssm =\
                            _filter_identical_interface_introns_on_best_pssm(
                            tmp_dict)


                    if intron_length_difference > 100:
                        # Very large length difference! Here, we likely have
                        # by-chance aligned coding nnGTnn nnAGnn residues.
                        # No balanced filtering here, but very rigidly based
                        # on overal exon signal
                        if [ _key ] == best_intron_msr.keys() and\
                        [ _key ] == best_intron_pssm.keys() and\
                        best_intron_msr.values()[0] >= 10.0:
                            ####################################################
                            if verbose:
                                # some printing
                                print "MSR-PSSM OVERRULED score::",
                                print (key,score), (_key,_score),
                                print "msr:", best_intron_msr.keys(),
                                print "pssm:", best_intron_pssm.keys(),
                                print "algExonRatio::",
                                print best_intron_msr.values()[0]
                            ####################################################
                            # delete the higher scoring intron
                            for del_intron in introndata[key]:
                                _cascade_delete_intron(del_intron,introndata)
                                        # delete this intron            
                            del( introndata[key] )
                            # break out by continue
                            continue

                        else:
                            # do not delete here, just go on in function
                            pass
                    elif intronObjmax.__class__.__name__ == 'SequenceErrorConnectingOrfs':
                        # allow the alternative, lower scoring intron to be deleted 
                        pass
                    elif intronObjsub.__class__.__name__ == 'SequenceErrorConnectingOrfs':
                        # allow the alternative, lower scoring intron to be deleted 
                        pass
                    else:
                        # detailed comparison of introns
                        intronObjmax.assign_bp_and_ppts()
                        intronObjsub.assign_bp_and_ppts()
                        status1 = _finetune_splicesite_comparison(intronObjmax.donor,intronObjsub.donor)
                        status2 = _finetune_splicesite_comparison(intronObjmax.acceptor,intronObjsub.acceptor)
                        status3 = _branchpoint_comparison(intronObjmax,intronObjsub)
                        status4 = _polypyrimidinetract_comparison(intronObjmax,intronObjsub)
                        check_statusses = [ status1, status2, status3 ]
    
                        # now evaluate. If [ _key ] == best_intron_msr.keys()
                        # and [ _key ] == best_intron_pssm.keys() and
                        # donor/acceptor/branchpoint is a clear improvement,
                        # overrule intron filtering on score.
                        if [ _key ] == best_intron_msr.keys() and\
                        check_statusses.count(False) == 0 and\
                        check_statusses.count(True) >= 1:
                            ########################################################
                            if verbose:
                                # some printing
                                print "MSR-PSSM OVERRULED score::",
                                print (key,score), (_key,_score),
                                print "msr:", best_intron_msr.keys(),
                                print "pssm:", best_intron_pssm.keys(),
                                print [status1,status2,status3,status4]
                            ########################################################
                            # delete the higher scoring intron
                            for del_intron in introndata[key]:
                                _cascade_delete_intron(del_intron,introndata)
                                        # delete this intron            
                            del( introndata[key] )
                            # break out by continue
                            continue
                        else:
                            ########################################################
                            if verbose:
                                # some printing
                                print "**NO** MSR-PSSM OVERRULED score::",
                                print (key,score), (_key,_score),
                                print "msr:", best_intron_msr.keys(),
                                print "pssm:", best_intron_pssm.keys(),
                                print [status1,status2,status3,status4],
                                print intronObjmax.get_branchpoint_nt_distance(),
                                print intronObjsub.get_branchpoint_nt_distance()
                            ########################################################
                            pass



            # if here, this intron(s) can be deleted
            # when it represents a tinyexon interface, cascade removal!
            for del_intron in introndata[_key]:
                _cascade_delete_intron(del_intron,introndata)

            # delete this intron            
            del( introndata[_key] )


    # second step. Filter out the garbage
    for interface in encountered_interfaces:
        this_interface_introns = {}
        for key in introndata.keys():
            score = _calc_aligned_score(introndata[key])
            intfA = _get_main_interface([i._label for i in introndata[key]])
            if score > SCORE_THRESHOLD:
                continue
            elif intfA != interface:
                continue
            else:
                pass
            this_interface_introns[ key ] = introndata[key]
        # check if still introns exist on this interface
        if not this_interface_introns:
            continue
        # calculate AlgPres / AlgSim score now...
        print "Garbage Filtering on:", interface, len(this_interface_introns)
        best_intron = _filter_identical_interface_introns_on_best_msr(
            this_interface_introns,abfgp_raw_exons,PCG,input,OPTIONS)

        if len(best_intron) == len(this_interface_introns):
            # garbage filtering on interface failed => interface not in genestructure?
            # do filtering on summed pssm score
            best_intron = _filter_identical_interface_introns_on_best_pssm(
                this_interface_introns)

        # remove all suboptimal introns
        for intronkey in this_interface_introns.keys():
            if not introndata.has_key(intronkey):
                # already deleted in previous iteration
                continue

            if intronkey != best_intron.keys()[0]:
                # if here, this intron(s) can be deleted
                # when it represents a tinyexon interface, cascade removal!
                for del_intron in introndata[intronkey]:
                    _cascade_delete_intron(del_intron,introndata)
        
                # delete this intron            
                del( introndata[intronkey] )

    # return filtered introndata
    return introndata

# end of function _filter_introns


def _cascade_delete_intron(intron,introndata,forced_delete=False):
    """ """
    if intron._gff['fsource'] in ['ABGPprojectingTE','ABGPmappingTE']:
        if hasattr(intron,'_linked_to_introns'):
            # get current intron key & score
            cur_key   = ( intron.donor.pos, intron.acceptor.pos )
            cur_score = _calc_aligned_score(introndata[cur_key])
            for cdi in intron._linked_to_introns:
                # Cascade Delete Intron -> originating from a
                # tinyexon effort. So, delete all introns which
                # originated from the same effort!
                _casc_key = ( cdi.donor.pos, cdi.acceptor.pos )
                if _casc_key == cur_key:
                    # do not delete this intron itself ;-)
                    continue
                elif introndata.has_key(_casc_key):
                    # get score of this linked intron
                    _score = _calc_aligned_score(introndata[_casc_key])
                    # only delete if _score < cur_score
                    if _score < cur_score or forced_delete:
                        print "cascade->", _casc_key, _score
                        del( introndata[_casc_key] )
                else:
                    # already deleted in this/previous iteration
                    pass

# end of function _cascade_delete_intron


def _get_main_interface(labels):
    """ Simple majority voting of list elements """
    ordered = [ ( labels.count(label), label ) for label in Set(labels) ]
    ordered.sort()
    return ordered[-1][1]
# end of function _get_main_interface


def _create_failed_intron_gff(introndata,assessed_interfaces,intron2label):
    """ """
    failed_introns = {}
    for kk in introndata.keys():
        start,end = kk
        introns = introndata[kk]
        label = intron2label[kk]
        value = assessed_interfaces[label]
        if start in [ s for s,e in failed_introns.keys() ]:
            if end < min([ e for s,e in failed_introns.keys() ]):
                # replace this key
                for s,e in failed_introns.keys():
                    if start == s and e == min([ _e for _s,_e in failed_introns.keys() ]):
                        del( failed_introns[(s,e)] )
                        failed_introns[kk] = value
                        break
        elif end in [ e for s,e in failed_introns.keys() ]:
            if start > max([ s for s,e in failed_introns.keys() ]):
                # replace this key
                for s,e in failed_introns.keys():
                    if end == e and s == max([ _s for _s,_e in failed_introns.keys() ]):
                        del( failed_introns[(s,e)] )
                        failed_introns[kk] = value
                        break
        else:
            # new key -> store
            failed_introns[kk] = value

    # create `failed` intron gff lines
    gfflines = []
    for ( start,end ) in failed_introns.keys():
        assesed_org_list = failed_introns[( start,end )]
        for orgid in assesed_org_list:
            if orgid in [ intron._reference for intron in introndata[(start,end)] ]:
                continue
            # if here, a failed organism
            gclass  = "failedintron" # "%s_%s" % (start, end)
            gname = "%s_%s_%s" % (orgid, start, end)

            orgSfullname = _get_organism_full_name(orgid)
            #if ABGP_ORGANISM2FULLSPECIESNAME_MAPPING.has_key(orgid):
            #    orgSfullname = ABGP_ORGANISM2FULLSPECIESNAME_MAPPING[orgid]
            #elif ABGP_ORGANISM2FULLSPECIESNAME_MAPPING.has_key(orgid[0:-1]):
            #    orgSfullname = ABGP_ORGANISM2FULLSPECIESNAME_MAPPING[orgid[0:-1]]
            #else:
            #    orgSfullname = orgid
            # create gclass/gname and GFF decoration string
            lastcol = """%s %s; Reference %s; Note '%s'""" % (
                gclass, gname, orgid, orgSfullname
            )
            # get fmethod/class from 1 of the intron(s)
            example = introndata[( start,end )][0]
            fmethod = example.__class__.__name__
            gffline = ( None, 'ABGPfailed', fmethod, start+1, end, ".", "+", ".", lastcol )
            gfflines.append( gffline )

    # return gfflines list
    return gfflines

# end of function _create_failed_intron_gff


def _update_observed_introns(introns,label,informant_label,organism,reference,
    observed_introns,intron2label,assessed_interfaces):
    """ """
    for intron in introns:
        intron._organism =  organism
        intron._reference = reference
        intron._label = label
        intron._informant_label = informant_label
        kk = intron.coords()
        if observed_introns[organism].has_key(kk):
            # Check if this reference organism - intron combi already exists.
            # This can occur for tinyTSSexons/introns
            if reference in [ _intron._reference for _intron in observed_introns[organism][kk] ]:
                ###print "REFERENCE ALREADY EXISTS!", label, reference, intron
                pass
            else:
                observed_introns[organism][kk].append( intron )
        else:
            observed_introns[organism][kk] = [ intron ]
        if not intron2label[organism].has_key(kk):
            intron2label[organism][kk] = label

    if assessed_interfaces[organism].has_key(label):
        # Check if reference already exists for this label.
        # This can occur for tinyTSSexons, repetitive/small PabcPORFs etc
        if reference not in assessed_interfaces[organism][label]:
            assessed_interfaces[organism][label].append( reference )
    else:
        assessed_interfaces[organism][label] = [ reference ]

# end of function _update_observed_introns

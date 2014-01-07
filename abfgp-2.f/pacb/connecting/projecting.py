"""
PacbPORF connection by missing intron projections
"""

# Module metadata
__authors__ = "Ate van der Burgt"
__license__ = "MIT"

# Pacb Imports
from pacb.validators import IsPacbPORF, IsIdenticalOrfs
from pacb.conversion import pacbp2pacbporf
from pacb.exceptions import InproperlyAppliedArgument
from pacb.connecting.orfs import (
    merge_orfs_with_intron,
    merge_orfs_with_two_tinyexons
    )
from pacb.connecting.functions import (
    _update_kwargs,
    set_apps_intron_query,
    set_apps_intron_sbjct,
    )

# Gene structure Imports
from gene.splicesite import (
    _score_splice_site,
    get_shared_nucleotides_at_splicesite,
    )
from gene.intron import (
    IntronConnectingOrfs,
    ProjectedIntronConnectingOrfs,
    )

# Other Imports
from parsers.fasta import parseFasta

# Python Imports
from os import popen as osPopen, remove as osRemove
from sets import Set

# Global variable Imports
from settings.executables import EXECUTABLE_SFM

from settings.splicesites import (
    IC_ACCEPTOR_PATTERN_OFFSET,
    IC_DONOR_PATTERN_OFFSET,
    KWARGS_PROJECTED_INTRON,
    KWARGS_PROJECTED_TINYEXON,
    )

AUSO = IC_ACCEPTOR_PATTERN_OFFSET[0] # 6nts; Acceptor UpStream Offset in IC PSSM
DDSO = IC_DONOR_PATTERN_OFFSET[1]    # 4nts; Donor DownStream Offset in IC PSSM


def merge_pacbporfs_by_intron_tinyexon_intron_in_query(pacbporfD,pacbporfA,
    orfSetObject,**kwargs):
    """ """
    return _merge_pacbporfs_by_tinyexon_and_two_introns(
        pacbporfD,pacbporfA,orfSetObject,"query",**kwargs)

# end of function merge_pacbporfs_by_intron_tinyexon_intron_in_query


def merge_pacbporfs_by_intron_tinyexon_intron_in_sbjct(pacbporfD,pacbporfA,
    orfSetObject,**kwargs):
    """ """
    return _merge_pacbporfs_by_tinyexon_and_two_introns(
        pacbporfD,pacbporfA,orfSetObject,"sbjct",**kwargs)

# end of function merge_pacbporfs_by_intron_tinyexon_intron_in_sbjct


def merge_pacbporfs_by_intron_in_sbjct(pfD,pfA,**kwargs):
    """ """
    # Orfs of QUERY must be identical
    IsIdenticalOrfs(pfD.orfQ,pfA.orfQ)
    return _merge_pacbporfs_by_intron(pfD,pfA,"sbjct",**kwargs)

# end of function merge_pacbporfs_by_intron_in_sbjct


def merge_pacbporfs_by_intron_in_query(pfD,pfA,**kwargs):
    """ """
    # Orfs of SBJCT must be identical
    IsIdenticalOrfs(pfD.orfS,pfA.orfS)
    return _merge_pacbporfs_by_intron(pfD,pfA,"query",**kwargs)

# end of function merge_pacbporfs_by_intron_in_query


def merge_pacbporfs_by_inframe_intron_in_sbjct(pfD,pfA,**kwargs):
    """ """
    # Inframe intron in Sbjct -> Sbjct Orfs must be identical too!!
    IsIdenticalOrfs(pfD.orfS,pfA.orfS)
    return merge_pacbporfs_by_intron_in_sbjct(pfD,pfA,**kwargs)

# end of function merge_pacbporfs_by_inframe_intron_in_sbjct


def merge_pacbporfs_by_inframe_intron_in_query(pfD,pfA,**kwargs):
    """ """
    # Inframe intron in Query -> Query Orfs must be identical too!!
    IsIdenticalOrfs(pfD.orfQ,pfA.orfQ)
    return merge_pacbporfs_by_intron_in_query(pfD,pfA,**kwargs)

# end of function merge_pacbporfs_by_inframe_intron_in_query


def project_splicesites_on_pacbporfs_with_lacking_intron_in_query(*args,**kwargs):
    """
    @attention: DEPRECATED NAME; alias for merge_pacbporfs_by_intron_in_sbjct
    """
    return merge_pacbporfs_by_intron_in_sbjct(*args,**kwargs)

# end of function project_splicesites_on_pacbporfs_with_lacking_intron_in_query


def project_splicesites_on_pacbporfs_with_lacking_intron_in_sbjct(*args,**kwargs):
    """
    @attention: DEPRECATED NAME; alias for merge_pacbporfs_by_intron_in_query
    """
    return merge_pacbporfs_by_intron_in_query(*args,**kwargs)

# end of function project_splicesites_on_pacbporfs_with_lacking_intron_in_sbjct


def merge_pacbporfs_by_two_tinyexons_in_query(pacbporfD,pacbporfA,
    orfSetObject,**kwargs):
    """ """
    return _merge_pacbporfs_by_two_tinyexons(
        pacbporfD,pacbporfA,orfSetObject,"query",**kwargs)

# end of function merge_pacbporfs_by_two_tinyexons_in_query


def merge_pacbporfs_by_two_tinyexons_in_sbjct(pacbporfD,pacbporfA,
    orfSetObject,**kwargs):
    """ """
    return _merge_pacbporfs_by_two_tinyexons(
        pacbporfD,pacbporfA,orfSetObject,"sbjct",**kwargs)

# end of function merge_pacbporfs_by_two_tinyexons_in_sbjct


################################################################################
#### Below, the actual functions for projection of introns                  ####
################################################################################


def _merge_pacbporfs_by_intron(pfD,pfA,queryorsbjct,verbose=False,**kwargs):
    """
    Project splicesites from SBJCT intron on continious QUERY PacbPORFs

    @type  pfD: PacbPORF object
    @param pfD: PacbPORF object that has to deliver (aligned) donor sites

    @type  pfA: PacbPORF object
    @param pfA: PacbPORF object that has to deliver (aligned) acceptor sites

    @type  queryorsbjct: string
    @param queryorsbjct: literal string 'query' or 'sbjct'

    @type  verbose: Boolean
    @param verbose: print debugging info to STDOUT when True

    @attention: see pacb.connecting.merge_orfs_with_intron for **kwargs)

    @rtype:  list
    @return: list with ProjectedIntrons (from Sbjct on Query)
    """
    # input validation
    IsPacbPORF(pfD)
    IsPacbPORF(pfA)

    # edit **kwargs dictionary for some forced attributes
    _update_kwargs(kwargs,KWARGS_PROJECTED_INTRON)

    ### if not kwargs.has_key('projected_intron_max_nt_offset'):
    ###    kwargs['projected_intron_max_nt_offset'] = PROJECTED_INTRON_MAX_NT_OFFSET
    ### if not kwargs.has_key('projected_intron_max_aa_offset'):
    ###    kwargs['projected_intron_max_aa_offset'] = PROJECTED_INTRON_MAX_AA_OFFSET

    # settings for minimal alignment entropy score
    min_donor_site_alignment_entropy = 0.0
    min_acceptor_site_alignment_entropy = 0.0


    ELEGIABLE_SPLICE_SITE_AA_RANGE = 75

    sposD = pfD._get_original_alignment_pos_start()
    eposD = pfD._get_original_alignment_pos_end()
    sposA = pfA._get_original_alignment_pos_start()
    eposA = pfA._get_original_alignment_pos_end()
    if queryorsbjct == "query":
        # Orfs of SBJCT must be identical
        IsIdenticalOrfs(pfD.orfS,pfA.orfS)
        donorOrf = pfD.orfQ
        accepOrf = pfA.orfQ
        prjctOrf = pfD.orfS # pfD.orfS == pfA.orfS
        dStart = sposD.query_dna_start  # ALIGNED start of donorPacbPORF
        dEnd   = pfD.query_dna_end      # ABSOLUTE end of donorPacbPORF
        aStart = pfA.query_dna_start    # ABSOLUTE start of acceptorPacbPORF
        aEnd   = eposA.query_dna_end    # ALIGNED end of acceptorPacbPORF
        outOfAlignmentAttribute = "sbjct_dna_start"
        # calculate elegiable splice site range
        qdr = pfD.alignment_dna_range_query()
        qar = pfA.alignment_dna_range_query()
        min_donor_pos = max([ min(qdr), max(qdr)-(ELEGIABLE_SPLICE_SITE_AA_RANGE*3) ])
        max_accep_pos = min([ max(qar), min(qar)+(ELEGIABLE_SPLICE_SITE_AA_RANGE*3) ])

    elif queryorsbjct == "sbjct":
        # Orfs of QUERY  must be identical
        IsIdenticalOrfs(pfD.orfQ,pfA.orfQ)
        donorOrf = pfD.orfS
        accepOrf = pfA.orfS
        prjctOrf = pfD.orfQ # pfD.orfQ == pfA.orfQ
        dStart = sposD.sbjct_dna_start  # ALIGNED start of donorPacbPORF
        dEnd   = pfD.sbjct_dna_end      # ABSOLUTE end of donorPacbPORF
        aStart = pfA.sbjct_dna_start    # ABSOLUTE start of acceptorPacbPORF
        aEnd   = eposA.sbjct_dna_end    # ALIGNED end of acceptorPacbPORF
        outOfAlignmentAttribute = "query_dna_start"
        # calculate elegiable splice site range
        sdr = pfD.alignment_dna_range_sbjct()
        sar = pfA.alignment_dna_range_sbjct()
        min_donor_pos = max([ min(sdr), max(sdr)-(ELEGIABLE_SPLICE_SITE_AA_RANGE*3) ])
        max_accep_pos = min([ max(sar), min(sar)+(ELEGIABLE_SPLICE_SITE_AA_RANGE*3) ])

    else:
        message = "'queryorsbjct' (%s), not 'query' or 'sbjct'" % queryorsbjct
        raise InproperlyAppliedArgument, message

    # predict introns only in `queryorsbjct` Orfs
    # introns is a list of IntronConnectingOrfs objects
    introns = merge_orfs_with_intron(donorOrf,accepOrf,
            min_donor_pos=min_donor_pos,
            max_acceptor_pos=max_accep_pos,
            order_by='length',**kwargs)

    # return list with projected introns
    projected_introns = []

    # gather unique donor and acceptor positions from list
    # of IntronConnectingOrfs
    for intron in introns:
        # break if intron is to large
        if kwargs['max_intron_nt_length'] and intron.length > kwargs['max_intron_nt_length']: break
        # continue if intron is to small
        if kwargs['min_intron_nt_length'] and intron.length < kwargs['min_intron_nt_length']: continue
        # continue if intron has non-canonical features


        # check if intron.start is on pfD;
        # inframe-introns can be projected outside of pfD/pfA area
        if intron.start <= dStart: continue
        if intron.start >= dEnd:   continue

        # check if intron.end is on pfA;
        # inframe-introns can be projected outside of pfD/pfA area
        if intron.end <= aStart: continue
        if intron.end >= aEnd:   continue

        if queryorsbjct == "sbjct":
            # get positions of donor & acceptor in the PacbPORF alignment
            donorPositionPos, phaseD = pfD.dnaposition_sbjct(intron.donor.pos,forced_return=True)
            accepPositionPos, phaseA = pfA.dnaposition_sbjct(intron.acceptor.pos,forced_return=True)
            # calculate projected distance on QUERY
            posDposQuery = pfD._positions[donorPositionPos].query_pos
            posAposQuery = pfA._positions[accepPositionPos].query_pos
            aaDistance   = posAposQuery - posDposQuery
        else:
            # get positions of donor & acceptor in the PacbPORF alignment
            donorPositionPos, phaseD = pfD.dnaposition_query(intron.donor.pos,forced_return=True)
            accepPositionPos, phaseA = pfA.dnaposition_query(intron.acceptor.pos,forced_return=True)
            # calculate binary entropy from projected position on SBJCT
            posDposSbjct = pfD._positions[donorPositionPos].sbjct_pos
            posAposSbjct = pfA._positions[accepPositionPos].sbjct_pos
            aaDistance   = posAposSbjct - posDposSbjct

        # calculate binary entropy score
        entropyDonorSbjct   = pfD.alignment_entropy(donorPositionPos,method='donor')
        entropyAcceptorSbjct= pfA.alignment_entropy(accepPositionPos,method='acceptor')

        # do distance check upon (projected) intron acceptance
        if abs(aaDistance) <= kwargs['max_aa_offset']:

            # check if we've runned out of the aligned part
            outofalignedpacbporf = False

            # get the projected donor position; mind the gap on this spot ;-)
            while pfD._positions[donorPositionPos].isa_gap and donorPositionPos > 0 :
                donorPositionPos -= 1
            else:
                projected_donor_position = getattr(pfD._positions[donorPositionPos],outOfAlignmentAttribute) + phaseD
                if donorPositionPos == 0 and pfD._positions[donorPositionPos].isa_gap:
                    print "WarningThatIsTackled::outofalignedpacbporf::donor"
                    outofalignedpacbporf = True

            # get the projected acceptor position; mind the gap on this spot ;-)
            while pfA._positions[accepPositionPos].isa_gap and len(pfA._positions) > accepPositionPos+1:
                accepPositionPos += 1
            else:
                projected_accep_position = getattr(pfA._positions[accepPositionPos],outOfAlignmentAttribute) + phaseA
                if accepPositionPos == len(pfA._positions)-1 and pfA._positions[accepPositionPos].isa_gap:
                    print "WarningThatIsTackled::outofalignedpacbporf::acceptor"
                    outofalignedpacbporf = True

            if not outofalignedpacbporf:
                ################################################################
                # set some meta-data properties to the intron object
                ################################################################
                # add distance score to intron
                intron._distance = abs(aaDistance)*3

                # add Alignment Positional Periphery Score into objects
                if queryorsbjct == "query":
                    succes = set_apps_intron_query(intron,pfD,pfA)
                else:
                    succes = set_apps_intron_sbjct(intron,pfD,pfA)
        
                # set GFF fsource attribute for recognition of intron sources
                intron._gff['fsource'] = "ABGPprojecting"

                # make a ProjectedIntronConnectingOrfs object
                pico = ProjectedIntronConnectingOrfs(prjctOrf,
                        projected_donor_position,
                        projected_accep_position)
                intron.binary_entropy_donor = entropyDonorSbjct
                intron.binary_entropy_acceptor = entropyAcceptorSbjct
                pico.add_projected_intron( intron )
                pico.phase = intron.phase
                projected_introns.append( pico )

                ################################################################
                if verbose:
                    print "PROJ::", intron._distance,
                    print (pfD.orfQ.id, pfA.orfQ.id),
                    print (pfD.orfS.id, pfA.orfS.id),
                    print "%s-%snt" % (intron.donor.pos, intron.acceptor.pos),
                    print "%2.1f,%2.1f" % (intron.donor.pssm_score, intron.acceptor.pssm_score),
                    print "%2.1f,%2.1f" % (intron.binary_entropy_donor,intron.binary_entropy_acceptor)
                ################################################################

        if aaDistance > kwargs['max_aa_offset']:
            # break out; ordered by length can never result in
            # a proper projected intron
            break


    # filter out less relevant ones compared to complete set of results
    projected_introns = _filter_projected_introns(projected_introns)

    # and return a list of ProjectedIntronConnectingOrfs
    return projected_introns

# end of function _merge_pacbporfs_by_intron


def _filter_projected_introns(picolist):
    """
    Helper function for _merge_pacbporfs_by_intron;
    Filter projection list for most likely one(s)
    """
    if not picolist: return picolist


    # settings for minimal alignment entropy score
    min_donor_site_alignment_entropy = 0.0
    min_acceptor_site_alignment_entropy = 0.0

    return_list = []
    best_intron_pos = []
    ignore_intron_pos = []
    introns = [ pico.projected_introns[0] for pico in picolist ]
    min_dist = min([ intron._distance for intron in introns ])
    max_entr = max([ intron.binary_entropy_donor + intron.binary_entropy_acceptor for intron in introns ])
    max_pssm = max([ intron.donor.pssm_score + intron.acceptor.pssm_score for intron in introns ])

    # if there are `perfect` entropy positioned examples,
    # remove all the garbage examples (but leave poor examples untouched)
    has_perfect_entropy_intron = False
    for intron in introns:
        if intron.binary_entropy_donor > min_donor_site_alignment_entropy and\
        intron.binary_entropy_acceptor > min_acceptor_site_alignment_entropy:
            has_perfect_entropy_intron = True
            break
    if has_perfect_entropy_intron:
        for pos in range(0,len(introns)):
            intron = introns[pos]
            # filter for SUB-zero, not for the threshold!
            if intron.binary_entropy_donor <= 0.0 or\
            intron.binary_entropy_acceptor <= 0.0:
                ignore_intron_pos.append( pos )

    for pos in range(0,len(introns)):
        intron = introns[pos]
        if intron._distance == min_dist:
            best_intron_pos.append( pos )
        if intron.binary_entropy_donor + intron.binary_entropy_acceptor == max_entr:
            best_intron_pos.append( pos )
        if intron.donor.pssm_score + intron.acceptor.pssm_score == max_pssm:
            best_intron_pos.append( pos )
    for unique_pos in Set(best_intron_pos):
        if unique_pos not in ignore_intron_pos or not\
        Set(best_intron_pos).symmetric_difference(ignore_intron_pos):
            return_list.append( picolist[unique_pos] )
    # return the list of filtered projected introns
    return return_list

# end of function _filter_projected_introns


def _merge_pacbporfs_by_tinyexon_and_two_introns(pacbporfD,pacbporfA,
    orfSetObject,queryorsbjct,verbose = False, **kwargs):
    """
    Merge 2 PacbPORF objects by introns

    @attention: see pacb.connecting.merge_orfs_with_intron for **kwargs)

    @type  pacbporfD: PacbPORF object
    @param pacbporfD: PacbPORF object that has to deliver PSSM donor objects

    @type  pacbporfA: PacbPORF object
    @param pacbporfA: PacbPORF object that has to deliver PSSM acceptor objects

    @type  orfSetObject: object with elegiable Orfs
    @param orfSetObject: object with elegiable Orfs

    @type  queryorsbjct: string
    @param queryorsbjct: literal string 'query' or 'sbjct'

    @type  verbose: Boolean
    @param verbose: print debugging info to STDOUT when True

    @rtype:  list
    @return: list with ( intron, ExonOnOrf, intron ) on the query sequence
    """
    # input validation
    IsPacbPORF(pacbporfD)
    IsPacbPORF(pacbporfA)

    # edit **kwargs dictionary for some forced attributes
    _update_kwargs(kwargs,KWARGS_PROJECTED_TINYEXON)

    MAX_TINYEXON_NT_LENGTH = 33
    MIN_TINYEXON_NT_LENGTH = 6

    tinyexons = []
    if queryorsbjct == "query":
        donorOrf = pacbporfD.orfQ
        accepOrf = pacbporfA.orfQ
        prjctOrf = pacbporfD.orfS
        alignedDonorRange = pacbporfD.alignment_dna_range_query()
        alignedAccepRange = pacbporfA.alignment_dna_range_query()
    elif queryorsbjct == "sbjct":
        donorOrf = pacbporfD.orfS
        accepOrf = pacbporfA.orfS
        prjctOrf = pacbporfD.orfQ
        alignedDonorRange = pacbporfD.alignment_dna_range_sbjct()
        alignedAccepRange = pacbporfA.alignment_dna_range_sbjct()
    else:
        message = "'queryorsbjct' (%s), not 'query' or 'sbjct'" % queryorsbjct
        raise InproperlyAppliedArgument, message

    for dObj in donorOrf._donor_sites:
        # do not make a projection OVER the aligned area
        if dObj.pos < min(alignedDonorRange): continue
        if queryorsbjct == "query":
            (dPos,dPhase) = pacbporfD.dnaposition_query(dObj.pos,forced_return=True)
        else:
            (dPos,dPhase) = pacbporfD.dnaposition_sbjct(dObj.pos,forced_return=True)
        try:
            algDobj = pacbporfD._positions[dPos]
        except IndexError:
            # site out of range of PacbPORF -> break
            break
        for aObj in accepOrf._acceptor_sites:
            # do not make a projection OVER the aligned area
            if aObj.pos > max(alignedAccepRange): continue
            if queryorsbjct == "query":
                (aPos,aPhase) = pacbporfA.dnaposition_query(aObj.pos,forced_return=True)
            else:
                (aPos,aPhase) = pacbporfA.dnaposition_sbjct(aObj.pos,forced_return=True)
            try:
                algAobj = pacbporfA._positions[aPos]
            except IndexError:
                # site out of range of PacbPORF -> break
                break
            if queryorsbjct == "query":
                posDsbjct = algDobj.sbjct_dna_start + dPhase
                posAsbjct = algAobj.sbjct_dna_start + aPhase
            else:
                posDsbjct = algDobj.query_dna_start + dPhase
                posAsbjct = algAobj.query_dna_start + aPhase
            distance = posAsbjct - posDsbjct
            if distance >= MAX_TINYEXON_NT_LENGTH:
                break
            if distance < MIN_TINYEXON_NT_LENGTH:
                continue

            ####################################################
            # generate a ScanForMatches pattern file
            ####################################################
            # example pattern: 6...6 AG NNGNNANNANNGN[2,0,0] GT 3...3
            query = list(prjctOrf.inputgenomicsequence[posDsbjct:posAsbjct])
            # mask all non-phase0 nucleotides to N residues;
            # this represents the regularexpression for a specific
            # peptide sequence
            firstphasepositions = range( 3-dPhase % 3, len(query), 3)
            for pos in range(0,len(query)):
                if pos not in firstphasepositions:
                    query[pos] = "N"
            # calculate a ~50% mismatch number
            mismatches =  max([ 0, (len(query) - query.count("N"))/2 ])
            # write the pattern to string and subsequently to file
            # example pattern: 6...6 AG NNGNNANNANNGN[2,0,0] GT 3...3
            if kwargs['allow_non_canonical_donor']:
                sfmpat = "%s...%s AG %s[%s,0,0] G (T | C) %s...%s" % (
                    AUSO,AUSO,"".join(query),mismatches,DDSO,DDSO)
            else:
                sfmpat = "%s...%s AG %s[%s,0,0] GT %s...%s" % (
                    AUSO,AUSO,"".join(query),mismatches,DDSO,DDSO)

            ####################################################
            if verbose:
                print (pacbporfD.orfQ.id,pacbporfA.orfQ.id),
                print distance, dObj, aObj
                print sfmpat
            ####################################################

            fname = "sfmpat_tinyexon_%s_%s_%s_%s" % (
                        donorOrf.id,
                        accepOrf.id,
                        posDsbjct,
                        posAsbjct,
                        )
            fh = open(fname,'w')
            fh.write(sfmpat+"\n")
            fh.close()

            ####################################################
            # run ScanForMatches
            ####################################################
            command = """echo ">myseq\n%s" | %s %s | tr "[,]" "\t\t#" | """ +\
                      """tr -d "\n " | sed "s/>/\\n>/g" | tr "#" "\t" | """ +\
                      """awk -F'\t' '{ if (NF==4 && $2>%s && $3<%s) """ +\
                      """{ print $1"["$2","$3"]\\n"$4 } }' """
            command = command % (
                        donorOrf.inputgenomicsequence,
                        EXECUTABLE_SFM,fname,
                        dObj.pos+(kwargs['min_intron_nt_length']-3),
                        aObj.pos-(kwargs['min_intron_nt_length']-3) )
            co = osPopen(command)
            matches = parseFasta(co.readlines())
            co.close()

            # filter matches for:
            # (1) correct donor & acceptor phase
            # (2) high enough donor & acceptor site scores
            for hdr,seqmatch in matches.iteritems():
                startQ,stopQ = [ int(item) for item in hdr.split(":")[1][1:-1].split(",") ]
                exonQstart   = startQ + AUSO + 2 - 1
                exonQstop    = stopQ  - DDSO - 2

                ####################################
                # get Orf object of tinyexon
                ####################################
                tinyexonorf = None
                # select the Orf on which the tinyexon is located
                for orfObj in orfSetObject.get_elegiable_orfs(
                max_orf_start=exonQstart,min_orf_end=exonQstop):
                    orfPhase = (exonQstart - orfObj.startPY) % 3
                    if orfPhase == dPhase:               
                        tinyexonorf = orfObj
                        break
                else:
                    # No tinyexonorf assigned!! Iin case a regex matched
                    # over a STOP-codon or the regex length is smaller
                    # then the smallest Orf, no Orf can be assigned
                    continue

                # filter for donor & acceptor score            
                dScore = _score_splice_site(seqmatch[-9:],splicetype='donor')
                aScore = _score_splice_site(seqmatch[0:11],splicetype='acceptor')
                if dScore < kwargs['min_donor_pssm_score']:
                    continue
                if aScore < kwargs['min_acceptor_pssm_score']:
                    continue

                # scan Orf for splicesites
                tinyexonorf.scan_orf_for_pssm_splice_sites(
                        splicetype="donor",
                        min_pssm_score=kwargs['min_donor_pssm_score'],
                        allow_non_canonical=kwargs['allow_non_canonical_donor'],
                        non_canonical_min_pssm_score=kwargs['non_canonical_min_donor_pssm_score'])
                tinyexonorf.scan_orf_for_pssm_splice_sites(
                        splicetype="acceptor",
                        min_pssm_score=kwargs['min_acceptor_pssm_score'],
                        allow_non_canonical=kwargs['allow_non_canonical_acceptor'],
                        non_canonical_min_pssm_score=kwargs['non_canonical_min_acceptor_pssm_score'])

                # get 1th intron donor object
                intron1_aObj = None
                for a in tinyexonorf._acceptor_sites:
                    if a.pos == exonQstart:
                        intron1_aObj = a
                        break
                else:
                    # pseudo-acceptorsite as found be SFM regex
                    # is not a valid acceptor site of high enough score
                    # continue to next iteration of (hdr,seqmatch) pair
                    continue

                # get 2th intron donor object
                intron2_dObj = None
                for d in tinyexonorf._donor_sites:
                    if d.pos == exonQstop:
                        intron2_dObj = d
                        break
                else:
                    # pseudo-donorsite as found be SFM regex
                    # is not a valid acceptor site of high enough score
                    # continue to next iteration of (hdr,seqmatch) pair
                    continue


                # check if introns are of elegiable lengths
                if (intron1_aObj.pos-dObj.pos) > kwargs['max_intron_nt_length']:
                    continue
                if (aObj.pos-intron2_dObj.pos) > kwargs['max_intron_nt_length']:
                    continue

                ####################################################
                if True or verbose:
                    # if here, a candidate!!!
                    print (pacbporfD.orfQ.id,tinyexonorf.id,pacbporfA.orfQ.id),
                    print hdr, dScore, aScore
                    print seqmatch
                ####################################################

                # append to found tinyexons
                query_data      = ( tinyexonorf, exonQstart, exonQstop )
                sbjct_data      = ( prjctOrf, posDsbjct, posAsbjct )
                splicesite_data = ( dObj,intron1_aObj, intron2_dObj, aObj )
                tinyexons.append( ( query_data, sbjct_data, splicesite_data ) )


            # file cleanup
            osRemove(fname)

    # return - End Of Function - if no tinyexons are found
    if not tinyexons:
        return []

    ####################################
    # select the **best** tinyexon
    ####################################
    (query_data,sbjct_data,splicesite_data) = tinyexons[0]
    orfQ,query_dna_start,query_dna_end = query_data
    orfS,sbjct_dna_start,sbjct_dna_end = sbjct_data
    (intron1_dObj,intron1_aObj,intron2_dObj,intron2_aObj) = splicesite_data

    ####################################################
    if verbose:
        print "tinyexon orf:", orfQ
        print "tinyexon orf:", intron1_aObj
        print "tinyexon orf:", intron2_dObj
    ####################################################

    ####################################
    # make tinyexon PacbPORF
    ####################################
    startQaa = orfQ.dnapos2aapos(query_dna_start) -1
    startSaa = orfS.dnapos2aapos(sbjct_dna_start) -1
    stopQaa  = orfQ.dnapos2aapos(query_dna_end) +1
    stopSaa  = orfS.dnapos2aapos(sbjct_dna_end) +1
    # check for directly leading stop codon on tinyexon
    while startQaa <= orfQ.protein_startPY:
        startQaa+=1
        startSaa+=1
        query_dna_start+=3
        sbjct_dna_start+=3
    while startSaa <= orfS.protein_startPY:
        startQaa+=1
        startSaa+=1
        query_dna_start+=3
        sbjct_dna_start+=3
    # check for directly tailing stop codon on tinyexon
    while stopQaa > orfQ.protein_endPY:
        stopQaa-=1
        stopSaa-=1
        query_dna_end-=3
        sbjct_dna_end-=3
    while stopSaa > orfS.protein_endPY:
        stopQaa-=1
        stopSaa-=1
        query_dna_end-=3
        sbjct_dna_end-=3
    # get sequences
    qAAseq = orfQ.getaas(abs_pos_start=startQaa,abs_pos_end=stopQaa)
    sAAseq = orfS.getaas(abs_pos_start=startSaa,abs_pos_end=stopSaa)

    ####################################################
    if verbose or len(qAAseq) != len(sAAseq):
        # if unequal lengths, error will be raised upon PacbP.__init__()
        print orfQ, qAAseq, startQaa, stopQaa, (stopQaa-startQaa),
        print (query_dna_start,query_dna_end)
        print orfS, sAAseq, startSaa, stopSaa, (stopSaa-startSaa),
        print (sbjct_dna_start,sbjct_dna_end)
        print orfQ.inputgenomicsequence[query_dna_start-2:query_dna_end+2]
        print orfS.inputgenomicsequence[sbjct_dna_start-2:sbjct_dna_end+2]
    ####################################################

    # initialize extended tinyexon PacbPORF
    from pacb import PacbP
    pacbp = PacbP(input=( qAAseq, sAAseq, startQaa, startSaa ) )
    pacbp.strip_unmatched_ends()
    pacbporf = pacbp2pacbporf(pacbp,orfQ,orfS)
    pacbporf.extend_pacbporf_after_stops()
    pacbporf.source = 'ABGPprojectingTE'

    ####################################
    # make introns
    ####################################
    intron1 = IntronConnectingOrfs(
                intron1_dObj, intron1_aObj, None,
                donorOrf,pacbporf.orfQ )
    intron2 = IntronConnectingOrfs(
                intron2_dObj, intron2_aObj, None,
                pacbporf.orfQ, accepOrf )


    ################################################################
    # set some meta-data properties to the intron objects
    ################################################################
    # add distance score to intron
    intron1._distance = 0
    intron2._distance = 0

    # add Alignment Positional Periphery Score into objects
    if queryorsbjct == "query":
        succes = set_apps_intron_query(intron1,pacbporfD,pacbporf)
        succes = set_apps_intron_query(intron2,pacbporf,pacbporfA)
    else:
        succes = set_apps_intron_sbjct(intron1,pacbporfD,pacbporf)
        succes = set_apps_intron_sbjct(intron2,pacbporf,pacbporfA)

    # set GFF fsource attribute for recognition of intron sources
    intron1._gff['fsource'] = "ABGPprojectingTE"
    intron2._gff['fsource'] = "ABGPprojectingTE"

    # create _linked_to_xxx attributes
    intron1._linked_to_pacbporfs = [ pacbporf ]
    intron2._linked_to_pacbporfs = [ pacbporf ]
    intron1._linked_to_introns   = [ intron2 ]
    intron2._linked_to_introns   = [ intron1 ]

    ####################################################
    if verbose:
        print pacbporf
        pacbporf.print_protein_and_dna()
        print intron1
        print intron2
        if False:
            # printing data when this function needs to be debugged:
            print ""
            print intron1
            print intron2
            print ""
            print pacbporfD
            pacbporfD.print_protein_and_dna()
            print ""
            print pacbporf
            pacbporf.print_protein_and_dna()
            print ""
            print pacbporfA
            pacbporfA.print_protein_and_dna()
            import sys
            sys.exit()
    ####################################################

    # return introns and intermediate tinyexon PacbPORF
    return [(intron1,intron2,pacbporf)]

# end of function _merge_pacbporfs_by_tinyexon_and_two_introns


def _merge_pacbporfs_by_two_tinyexons(pacbporfD,pacbporfA,
    orfSetObject,queryorsbjct,verbose = False, **kwargs):
    """ """
    # edit **kwargs dictionary for some forced attributes
    _update_kwargs(kwargs,KWARGS_PROJECTED_TINYEXON)

    tinyexons = []
    sposD = pacbporfD._get_original_alignment_pos_start()
    eposD = pacbporfD._get_original_alignment_pos_end()
    sposA = pacbporfA._get_original_alignment_pos_start()
    eposA = pacbporfA._get_original_alignment_pos_end()
    if queryorsbjct == "query":
        donorOrf = pacbporfD.orfQ
        accepOrf = pacbporfA.orfQ
        prjctOrf = pacbporfD.orfS
        dStart,dEnd = sposD.query_dna_start, eposD.query_dna_end
        aStart,aEnd = sposA.query_dna_start, eposA.query_dna_end
    elif queryorsbjct == "sbjct":
        donorOrf = pacbporfD.orfS
        accepOrf = pacbporfA.orfS
        prjctOrf = pacbporfD.orfQ
        dStart,dEnd = sposD.sbjct_dna_start, eposD.sbjct_dna_end
        aStart,aEnd = sposA.sbjct_dna_start, eposA.sbjct_dna_end
    else:
        message = "'queryorsbjct' (%s), not 'query' or 'sbjct'" % queryorsbjct
        raise InproperlyAppliedArgument, message

    # get all potential combinations of two tinyexons
    tinyexoncombis = merge_orfs_with_two_tinyexons(
                donorOrf, accepOrf,
                donorOrf._donor_sites,
                accepOrf._acceptor_sites,
                orfSetObject.orfs,
                )

    results = []

    for dObj in donorOrf._donor_sites:
        if queryorsbjct == "query":
            (dPos,dPhase) = pacbporfD.dnaposition_query(dObj.pos,forced_return=True)
        else:
            (dPos,dPhase) = pacbporfD.dnaposition_sbjct(dObj.pos,forced_return=True)
        try:
            algDobj = pacbporfD._positions[dPos]
        except IndexError:
            # site out of range of PacbPORF -> break
            break

        # check if dObj is on pfD;
        # introns of tinyexons can be projected outside of pfD/pfA area
        if dObj.pos < dStart: continue

        for aObj in accepOrf._acceptor_sites:
            if queryorsbjct == "query":
                (aPos,aPhase) = pacbporfA.dnaposition_query(aObj.pos,forced_return=True)
            else:
                (aPos,aPhase) = pacbporfA.dnaposition_sbjct(aObj.pos,forced_return=True)
            try:
                algAobj = pacbporfA._positions[aPos]
            except IndexError:
                # site out of range of PacbPORF -> break
                break

            # check if aObj is on pfA;
            # introns of tinyexons can be projected outside of pfD/pfA area
            if aObj.pos > aEnd: continue

            if queryorsbjct == "query":
                posDsbjct = algDobj.sbjct_dna_start + dPhase
                posAsbjct = algAobj.sbjct_dna_start + aPhase
            else:
                posDsbjct = algDobj.query_dna_start + dPhase
                posAsbjct = algAobj.query_dna_start + aPhase
            distance = posAsbjct - posDsbjct
            if distance >= (kwargs['max_tinyexon_nt_length']*2):
                break
            if distance < (kwargs['min_tinyexon_nt_length']*2):
                continue

            filtered_tinyexoncombis = _filter_tinyexoncombis(tinyexoncombis,
                    min_length = distance,
                    max_length = distance,
                    min_first_acceptor_pos = dObj.pos + kwargs['min_tinyexon_intron_nt_length'],
                    max_final_donor_pos = aObj.pos - kwargs['min_tinyexon_intron_nt_length'],
                    phase_final_donor = aObj.phase,
                    phase_first_acceptor= dObj.phase,
                    )

            if not filtered_tinyexoncombis: continue

            ####################################################################
            if verbose:
                print distance, dObj, aObj, len(tinyexoncombis),
                print len(filtered_tinyexoncombis)
            ####################################################################

            for exon1,intron,exon2 in filtered_tinyexoncombis:
                # make preceding intron
                preceding_intron = IntronConnectingOrfs(
                    dObj,exon1.acceptor,
                    None,donorOrf,exon1.orf )

                # make subsequent intron
                subsequent_intron = IntronConnectingOrfs(
                    exon2.donor, aObj,
                    None,exon2.orf,accepOrf)

                ################################################################
                if verbose:
                    print "\t", exon1, exon1.proteinsequence(),
                    print preceding_intron.phase, exon1.donor.phase,
                    print subsequent_intron.phase, preceding_intron.shared_aa,
                    print intron.shared_aa, subsequent_intron.shared_aa 
                    print "\t", exon2, exon2.proteinsequence()
                ################################################################

                # get prjctOrf sequence for comparison
                correctionA = 0
                if aObj.phase != 0:
                    # INCLUDE the final AA which is broken by the splicesite
                    correctionA=1
                if queryorsbjct == "query":
                    startPos,_phase = pacbporfD.dnaposition_query(dObj.pos,forced_return=True)
                    stopPos,_phase  = pacbporfA.dnaposition_query(aObj.pos,forced_return=True)
                    start = pacbporfD._positions[startPos].sbjct_pos
                    stop  = pacbporfA._positions[stopPos].sbjct_pos + correctionA
                else:
                    startPos,_phase = pacbporfD.dnaposition_sbjct(dObj.pos,forced_return=True)
                    stopPos,_phase  = pacbporfA.dnaposition_sbjct(aObj.pos,forced_return=True)
                    start = pacbporfD._positions[startPos].query_pos
                    stop  = pacbporfA._positions[stopPos].query_pos + correctionA

                if stop <= start:
                    # tinyexon is so tiny that is does not have a single
                    # full aligned AA -> discard here
                    continue

                # actually get the prjctOrf sequence
                aaseq = prjctOrf.getaas(abs_pos_start=start,abs_pos_end=stop)

                # initialize a PacbP for the combination of both tinyexons
                # afterwards, check if the indentityscore is > 0.XX
                from pacb import PacbP
                seqparts = [ preceding_intron.shared_aa,
                             exon1.proteinsequence(),
                             intron.shared_aa,
                             exon2.proteinsequence(),
                             subsequent_intron.shared_aa ]

                ################################################################
                if verbose or len("".join(seqparts)) != len(aaseq):
                    print pacbporfD
                    print exon1.orf, exon2.orf, prjctOrf
                    print pacbporfA
                    print seqparts
                    print aaseq, len(aaseq), len("".join(seqparts)), (start,stop)
                    print "'%s'" % queryorsbjct,
                    print "Q", (algDobj.query_pos, algAobj.query_pos),
                    print "S", (algDobj.sbjct_pos, algAobj.sbjct_pos)
                    print "distance:", distance, kwargs['max_tinyexon_nt_length'],
                    print (posDsbjct, posAsbjct),
                    print "Q-dna:", ( algDobj.query_dna_start, dPhase, algAobj.query_dna_start, aPhase ),
                    print "S-dna:", ( algDobj.sbjct_dna_start, dPhase, algAobj.sbjct_dna_start, aPhase )
                ################################################################

                # ignore by continue when sequences not identical in length
                if len("".join(seqparts)) != len(aaseq): continue

                testpacbp = PacbP(input=( "".join(seqparts), aaseq, 0, 0) )
                testpacbp.strip_unmatched_ends()

                if not ( testpacbp.identityscore > 0.60 and\
                (float(testpacbp.length) / len(aaseq)) > 0.70 ):
                    # not a very convincing alignment
                    continue

                ################################################################
                if verbose:
                    print testpacbp
                    testpacbp.print_protein()
                ################################################################

                # if here, succesfully mapped 2 tiny exons!!
                # get all sequences/coordinates in place for
                # pacbporf formation
                orfQ1   = exon1.orf
                orfS1   = prjctOrf
                orfQ2   = exon2.orf
                orfS2   = prjctOrf
                seqQ1   = exon1.proteinsequence()
                seqQ2   = exon2.proteinsequence()
                coordQ1 = exon1.acceptor.pos / 3
                coordS1 = start
                coordQ2 = exon2.acceptor.pos / 3
                coordS2 = start + len(seqparts[0]) + len(seqparts[1]) + len(seqparts[2])
                seqS1   = aaseq[0:(len(seqparts[0])+len(seqparts[1]))]
                seqS2   = aaseq[-(len(seqparts[3])+len(seqparts[4])):]
                if len(seqparts[0]):
                    seqS1 = seqS1[1:]
                    coordS1 += 1
                if len(seqparts[4]):
                    seqS2 = seqS2[:-1]

                if queryorsbjct == "sbjct": 
                    # swap query <-> sbjct
                    orfQ1,orfS1 = orfS1,orfQ1 
                    orfQ2,orfS2 = orfS2,orfQ2
                    seqQ1,seqS1 = seqS1,seqQ1
                    seqQ2,seqS2 = seqS2,seqQ2
                    coordQ1,coordS1 = coordS1,coordQ1
                    coordQ2,coordS2 = coordS2,coordQ2

                ################################################################
                if verbose:
                    print "tinypacbporf1:", seqQ1, seqQ2, coordQ1, coordQ2
                    print "tinypacbporf2:", seqS1, seqS2, coordS1, coordS2
                ################################################################


                # make pacbporfs
                pacbp1 = PacbP(input=( seqQ1, seqS1, coordQ1, coordS1) )
                pacbp1.strip_unmatched_ends()
                tinypacbporf1 = pacbp2pacbporf(pacbp1,orfQ1,orfS1)
                tinypacbporf1.extend_pacbporf_after_stops()
                pacbp2 = PacbP(input=( seqQ2, seqS2, coordQ2, coordS2) )
                pacbp2.strip_unmatched_ends()
                tinypacbporf2 = pacbp2pacbporf(pacbp2,orfQ2,orfS2)
                tinypacbporf2.extend_pacbporf_after_stops()

                ################################################################
                if verbose:
                    print tinypacbporf1
                    tinypacbporf1.print_protein_and_dna()
                    print tinypacbporf2
                    tinypacbporf2.print_protein_and_dna()
                ################################################################


                ################################################################
                # set some meta-data properties to the intron objects
                ################################################################
                # add distance score to intron
                preceding_intron._distance  = 0
                intron._distance            = 0
                subsequent_intron._distance = 0
            
                # add Alignment Positional Periphery Score into objects
                if queryorsbjct == "query":
                    succes = set_apps_intron_query(preceding_intron,pacbporfD,tinypacbporf1)
                    succes = set_apps_intron_query(intron,tinypacbporf1,tinypacbporf2)
                    succes = set_apps_intron_query(subsequent_intron,tinypacbporf2,pacbporfA)
                else:
                    succes = set_apps_intron_sbjct(preceding_intron,pacbporfD,tinypacbporf1)
                    succes = set_apps_intron_sbjct(intron,tinypacbporf1,tinypacbporf2)
                    succes = set_apps_intron_sbjct(subsequent_intron,tinypacbporf2,pacbporfA)
            
                # set GFF fsource attribute for recognition of intron sources
                preceding_intron._gff['fsource']  = "ABGPprojectingTE"
                intron._gff['fsource']            = "ABGPprojectingTE"
                subsequent_intron._gff['fsource'] = "ABGPprojectingTE"


                # create _linked_to_xxx attributes
                preceding_intron._linked_to_pacbporfs = [ tinypacbporf1, tinypacbporf2 ]
                intron._linked_to_pacbporfs = [ tinypacbporf1, tinypacbporf2 ]
                subsequent_intron._linked_to_pacbporfs = [ tinypacbporf1, tinypacbporf2 ]
                preceding_intron._linked_to_introns   = [ intron,subsequent_intron ]
                intron._linked_to_introns             = [ preceding_intron,subsequent_intron ]
                subsequent_intron._linked_to_introns  = [ intron,preceding_intron ]

                ################################################################
                # append to results
                ################################################################
                results.append( (
                    preceding_intron,
                    intron,
                    subsequent_intron,
                    tinypacbporf1,
                    tinypacbporf2,
                    ) )


    # return 3 introns and 2 intermediate tinyexon PacbPORFs (per row)
    return results

# end of function _merge_pacbporfs_by_two_tinyexons


def _filter_tinyexoncombis(tinyexoncombis,
    min_length=None,
    max_length=None,
    min_first_acceptor_pos=None,
    max_final_donor_pos=None,
    phase_final_donor=None,
    phase_first_acceptor=None):
    """
    Helper function for _merge_pacbporfs_by_two_tinyexons
    """
    elegiable = []
    for (exon1,intron,exon2) in tinyexoncombis:
        cur_length = exon1.length + exon2.length
        if min_length and cur_length < min_length:
            continue
        if max_length and cur_length > max_length:
            continue
        if (phase_final_donor or phase_final_donor==0) and\
        exon2.donor.phase != phase_final_donor:
            continue
        if (phase_first_acceptor or phase_first_acceptor==0) and\
        exon1.acceptor.phase != phase_first_acceptor:
            continue
        if max_final_donor_pos and exon2.donor.pos > max_final_donor_pos:
            continue
        if min_first_acceptor_pos and exon1.acceptor.pos < min_first_acceptor_pos:
            continue

        # if here: elegiable combi!
        elegiable.append( (exon1,intron,exon2) )

    # return list of elegiable tinyexon combis
    return elegiable
        
# end of function _filter_tinyexoncombis

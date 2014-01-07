"""
Functions for obtaining SignalP data for Alignment Based Gene Prediction
"""

# Module metadata
__authors__ = "Ate van der Burgt"
__license__ = "MIT"

# Pacb Imports
from pacb.connecting.orfs import get_signalp_first_exon_on_orf
from pacb.exceptions import CoordinateOutOfRange
from pacb.conversion import pacbp_from_clustalw, pacbp2pacbporf
from pacb.ordering import order_pacbporf_list
from pacb import PacbPCOORDS

# graphAbgp Imports
from graphAbgp.graph_pacbpcollection import pacbporf2PCG


# Gene Imports
from gene.signalpeptide import SignalPSignalPeptide
from gene.start import TranslationalStartSite

# Other ABGP Imports
from lib_clustalw import clustalw
from pythonlibs.ordering import order_list_by_attribute

# Python Imports
from sets import Set

# Global Variable Imports
from settings.dbwarehouse import _get_organism_full_name
USE_SIGNALP = True
SIGNALP_FIRSTEXON_MAX_INTRON_NT_LENGTH = 250


def obtainlocussignalpdata(input):
    """ """
    cnt = 0
    for org in input.keys():
        if USE_SIGNALP:
            # store the locus SignalP predictions to the Orf objects
            cnt += input[org]['gldobj'].load_locus_signalp_to_orfs()
            input[org]['gldobj'].load_protein_signalp_to_orf()

    # return number of loaded signalpeptides
    return cnt

# end of function obtainlocussignalpdata



class ProjectedSignalPSignalPeptide(SignalPSignalPeptide):
    """ """
    pass


def create_projected_signalp_for_informants(organism,PCG,input,
    minimal_aa_overlap=None,verbose=False):
    """ """
    # return list with ProjectedSignalPeptides
    return_projsignalpeptides = []

    for orgS in PCG.organism_set():
        if organism == orgS: continue
        if verbose: print orgS, organism
        # dict in order to create unique projections
        projections = {}
        for pacbporf in PCG.get_pacbps_by_organisms(organism,orgS):

            if not (pacbporf.orfS._has_signalp_sites_predicted and\
            pacbporf.orfS._signalp_sites):
                # no SignalPeptides on this Orf of the informant
                continue

            if minimal_aa_overlap and not pacbporf.orfQ._signalp_sites:
                # overlap required AND no SignalPeptides on the target Orf
                continue
            else:
                # calculate DNA range of target gene's SignalPeptides
                qdnarange = []
                for spQ in pacbporf.orfQ._signalp_sites:
                    qdnarange.extend( range(spQ.start,spQ.end) )


            # get full name tag of informant
            orgSfullname =  "%s [%s]" % (
                    _get_organism_full_name(orgS,truncate=True),
                    input[orgS]['proteinfref'] )

            # loop over the informant's SignalPeptides
            for spS in pacbporf.orfS._signalp_sites:
                prjQSPstart = pacbporf.dnapos_sbjct2query(spS.start)
                prjQSPend   = pacbporf.dnapos_sbjct2query(spS.end)
                prjQTSSstart= pacbporf.dnapos_sbjct2query(spS.tss.start)
                prjQTSSend = pacbporf.dnapos_sbjct2query(spS.tss.end)
                coords = [ prjQSPstart, prjQSPend, prjQTSSstart,prjQTSSend ]
                if CoordinateOutOfRange in coords:
                    ############################################################
                    if verbose: print "OUT-OF-RANGE:", spS
                    ############################################################
                    # projection on query NOT possible (out of range of Orf)                    
                    continue

                ################################################################
                if verbose:
                    print orgS, pacbporf, len(pacbporf.orfQ._signalp_sites), 
                    print len(pacbporf.orfS._signalp_sites)
                ################################################################

                # calculate AA overlap with target organism SignalPeptides    
                sdnarange = Set(range(prjQSPstart,prjQSPend))
                dna_overlap = sdnarange.intersection(qdnarange)
                aa_overlap  = len(dna_overlap)/3

                if minimal_aa_overlap and aa_overlap < minimal_aa_overlap:
                    # minimal overlap required and not achieved
                    continue
    
                # create a ProjectedSignalPSignalPeptide
                prjTSS = TranslationalStartSite(
                    prjQTSSstart,"n"*19,
                    pssm_score=spS.tss.pssm_score )
                prjTSS._gff['fmethod'] = 'projectedTSSpssm'
                prjTSS._gff['column9data'] = {'Informant': orgSfullname}
                prjSignalP_gff = {
                    'fmethod'    : 'projectedSignalPeptide',
                    'column9data': {'Informant': orgSfullname},
                    'gname'      : "%s_%s_%s" % (orgS,prjQSPstart,prjQSPend) }
                prjSignalP = ProjectedSignalPSignalPeptide(
                    prjQSPstart,prjQSPend,
                    spS.pssm_score,tss=prjTSS,
                    gff=prjSignalP_gff )
                prjSignalP._gff['fmethod'] = 'projectedSignalPeptide'
                prjSignalP._gff['column9data'] = {'Informant': orgSfullname}

                # store to prjections dict
                projections[(prjSignalP.start,prjSignalP.end)] = prjSignalP

                ################################################################
                if verbose: print prjSignalP, aa_overlap
                ################################################################

        # store unique projections to return_projsignalpeptides
        for sgp in projections.values():
            return_projsignalpeptides.append(sgp)

    # return list of ProjectedSignalPeptides
    return return_projsignalpeptides

# end of function create_projected_signalp_for_informants


def predict_signalp_exons(input,GENE_IDENTIFIER_SET,OPTIONS,verbose=False):
    """ """
    signalpexon_in_organisms = {}
    for organism in GENE_IDENTIFIER_SET:
        for orf in input[organism]['orfs'].orfs:
            signalpexon = get_signalp_first_exon_on_orf(orf)
            if signalpexon:
                if signalpexon_in_organisms.has_key(organism):
                    signalpexon_in_organisms[organism].append( signalpexon )
                else:
                    signalpexon_in_organisms[organism] = [ signalpexon ]
                ################################################################
                if verbose:
                    print ">>",organism, signalpexon
                    print signalpexon.proteinsequence()
                ################################################################

    # return dict with signalpexon_in_organisms
    return signalpexon_in_organisms

# end of function predict_signalp_exons


def update_PCG_with_signalpexons(signalpexonseqs,PCG,OPTIONS,
    min_pacbporf_identityscore=0.20,verbose=True):
    """ """
    if not signalpexonseqs.has_key(OPTIONS.target): return False
    is_any_pacbporf_added = False
    for targetSPexon in signalpexonseqs[OPTIONS.target]:
        target = OPTIONS.target
        for informant,infSPlist in signalpexonseqs.iteritems():
            if informant == OPTIONS.target: continue
            # check if informant has been deleted in the meanwhile
            if informant not in PCG.organism_set(): continue
            # list to store signalp exons into
            signalpexon_pacbp_list = []
            # get ordered pacbporfs fromt he PCG
            thepacbporfs = order_pacbporf_list(PCG.get_pacbps_by_organisms(OPTIONS.target,informant))
            if not thepacbporfs:
                # no alignments present for this organism (can happen!)
                continue
            for informantSPexon in infSPlist:
                coords  = [ targetSPexon.protein_start(),
                            targetSPexon.protein_end(),
                            informantSPexon.protein_start(),
                            informantSPexon.protein_end(), ]

                # prior to making ClustalW-PacbP, check PacbPCOORD placeability
                # into the list of pacbporfs
                pacbpCoordsObj = PacbPCOORDS(input=(
                        targetSPexon.proteinsequence(),
                        informantSPexon.proteinsequence(),
                        targetSPexon.protein_start(),
                        informantSPexon.protein_start(),
                        ) )

                if False in [ pacbpCoordsObj.is_positioned_compatibly(pacbporf) for pacbporf in thepacbporfs ]:
                    # *NOT* placable in current ordered list of PacbPORFS
                    continue

                dist = pacbpCoordsObj.distance_towards(thepacbporfs[0])
                if dist > SIGNALP_FIRSTEXON_MAX_INTRON_NT_LENGTH/3:
                    # WAY TO FAR in front of current gene structure parts.
                    # Do not allow (pooras a *NOT* placable in current ordered list of PacbPORFS
                    continue
                elif dist == 0:
                    # NOT placeable in front of the rest of the PacbPORFS.
                    continue
                else:
                    pass

                # perform ClustalW alignment on the SP exons
                    (alignedseqs,alignment) =\
                clustalw( seqs= { 
                    OPTIONS.target: targetSPexon.proteinsequence(),
                    informant: informantSPexon.proteinsequence() } )

                # make pacbp from clustalw alignment
                pacbp = pacbp_from_clustalw(
                            alignment=(
                                    alignedseqs[OPTIONS.target],
                                    alignment,
                                    alignedseqs[informant]
                                    ),
                            coords=coords
                            )

                # is there any alignment constructed?
                if not pacbp: continue

                # ignore (very) poor identyscore alignments
                if pacbp.identityscore < min_pacbporf_identityscore: continue

                # if here make extended pacbpORF
                signalpexonPacbpORF = pacbp2pacbporf(pacbp,
                        targetSPexon.orf,informantSPexon.orf)
                signalpexonPacbpORF.extend_pacbporf_after_stops()
                # and store in signalpexon_pacbp_list
                signalpexon_pacbp_list.append( signalpexonPacbpORF )

                ################################################################
                if verbose:
                    print alignedseqs[OPTIONS.target], OPTIONS.target
                    print alignment
                    print alignedseqs[informant], informant
                    if pacbp:
                        print pacbp, (OPTIONS.target, targetSPexon.orf.id),
                        print (informant, informantSPexon.orf.id),
                        print "DISTANCE::", dist
                        pacbp.print_protein()
                        print ""
                ################################################################

            # If there are signalpexon-guided pacbporfs found, store the one
            # with the highest bitscore
            if signalpexon_pacbp_list:
                signalpexon_pacbp_list = order_list_by_attribute(
                        signalpexon_pacbp_list,order_by='bits',reversed=True)
                # store best bitscoring pacbporf to PCG
                signalp_pacbporf = signalpexon_pacbp_list[0]
                pacbporf2PCG(signalp_pacbporf,OPTIONS.target,informant,PCG,source='SignalP-ClustalW') 
                is_any_pacbporf_added = True
                ####################################################################
                if verbose:
                    print "SignalP Exon added to PCG:", signalp_pacbporf, informant
                ####################################################################
            else:
                pass

    # return pointer is_any_pacbporf_added
    return is_any_pacbporf_added

# end of function update_PCG_with_signalpexons

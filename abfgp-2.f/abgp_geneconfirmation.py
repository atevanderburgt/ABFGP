"""
Functions for checking the given (annotated) gene structure on a DNA sequence of Orfs
"""

# Module metadata
__authors__ = "Ate van der Burgt"
__license__ = "MIT"

# Python Imports 
import os
from re import finditer
from sets import Set

# import ABGP functions
import dna2prot
from gff import (
        parsegfffile, mirrorgfflist,
        order_gff_list, filtergffs4fmethod,
        gffs2coordset, coordset2gfftracks,
        gffs2txt,
        )
from lib_orfset import GetorfOrfSet
from graphAbgp.ordering import order_list_by_attribute
from abgp_warnings import (
        CodingExonNotFoundAsOrfWarning,
        PotentialSequenceErrorWarning,
        SmallAnnotatedExonWarning,
        SmallAnnotatedFirstExonWarning,
        SmallAnnotatedFinalExonWarning,
        SmallAnnotatedIntronWarning,
        DnaSequenceContainsNsymbolsWarning,
        NonCanonicalSpliceSiteWarning,
        ImplausibleSpliceSiteWarning,
        GeneStructureIsNotMappableOnOrfsWarning,
        IncompleteGeneStructureWarning,
        )

from settings.genestructure import (
        FIRST_ANNOTATED_EXON_LABEL,
        FINAL_ANNOTATED_EXON_LABEL,
        IS_ANNOTATED_EXON_LABEL,
        ORF_IS_UNIGENE_LABEL,
        )

from settings.gff import (
        # Unigene tracks
        GFF_UNIGENE_GCLASS,
        GFF_UGEXON_FMETHOD,
        GFF_UNIGENE_FMETHOD,
        # Gene tracks
        GFF_GENE_GCLASS,
        GFF_CDS_FMETHOD, GFF_EXON_FMETHOD,
        GFF_GENE_ELEMENT_FMETHOD_LIST,
        GFF_TRANSCRIPT_GCLASS_LIST,
        GFF_GENESTART_FMETHOD, GFF_GENESTOP_FMETHOD,
        # UTR tracks
        GFF_UTR5_FSOURCE, GFF_UTR5_FMETHOD,
        GFF_UTR3_FSOURCE, GFF_UTR3_FMETHOD,
        )
from settings.splicesites import (
        CANONICAL_DONOR_SITES,
        NON_CANONICAL_DONOR_SITES,
        CANONICAL_ACCEPTOR_SITES,
        )
from settings.executables import (
        EXECUTABLE_UNIGENEANNOTATION,
        EXECUTABLE_GFF2FASTA,
        PYTHON_PATH,
        )


def dna2protein(_sequence,_frame):
    """
    Translate a DNA sequence into AAs for a given frame

    @type  _sequence: string
    @param _sequence: DNA seqeunce string

    @type  _frame: integer
    @param _frame: desired reading frame (0,1,2)

    @rtype:  string
    @return: Protein sequence
    """
    _untill = (len(_sequence)-_frame) % 3
    if not _untill:  _untill = len(_sequence)
    else:            _untill = -_untill
    # and do the actual translation
    return dna2prot.dna2protein(_sequence[_frame:_untill])

# end of function dna2protein


def annotatedgeneexonsizeevaluation(input,small_exon_nt_threshold=25,small_intron_nt_threshold=42):
    """
    Check the size of the annotated coding exons
    
    @type  input: dict
    @param input: input dictionary data structure

    @type  verbose: Boolean
    @param verbose: print discrepancies to STDOUT (True) or not (False, default)

    @rtype:  Boolean
    @return: are small coding exons present in the available annotations?
    """
    HAS_TINY_EXONS = False
    keysandtype = [
            ('gff-gene', 'orfid-genestructure', 'annotated gene'),
            ('gff-unigene', 'orfid-unigenestructure', 'unigene')
            ]
    for org in input.keys(): 
        for (keyA,keyB,typeofevidence) in keysandtype:
            if input[org][keyA]:
                ###for gfftrack in input[org][keyA]:
                ###    if gfftrack[2] != GFF_CDS_FMETHOD: continue
                codingexons = filtergffs4fmethod(input[org][keyA],GFF_CDS_FMETHOD)
                for pos in range(0,len(codingexons)):
                    gfftrack = codingexons[pos]
                    start, end = int(gfftrack[3])-1, int(gfftrack[4])
                    if ( end - start ) <= small_exon_nt_threshold:
                        HAS_TINY_EXONS = True

                        # a small exon! create a warning message
                        message = "organism: '%s' %s (gff) start: %s end: %s length: %s" % (
                                    org, typeofevidence, (start+1), end,  (end - start) )
                        # a Try to find the id of the orf
                        if input[org][keyB]:
                            for id in input[org][keyB]:
                                orf = input[org]['orfs'].get_orf_by_id(id)
                                if orf.start <= start+1 and orf.end >= end:
                                    message+=" on %s" % (str(orf))
                                    break
                        else:
                            pass

                        # append the Warning
                        input[org]['warnings'].append( SmallAnnotatedExonWarning(message) )

                        # check if it is a small First/Final exon
                        if pos == 0:
                            input[org]['warnings'].append( SmallAnnotatedFirstExonWarning(message) )
                        elif pos == len(codingexons)-1:
                            input[org]['warnings'].append( SmallAnnotatedFinalExonWarning(message) )

                # find small annotated introns
                for pos in range(1,len(codingexons)):
                    gfftrackP = codingexons[pos-1]
                    gfftrackN = codingexons[pos]
                    sta = int(gfftrackP[4])
                    end = int(gfftrackN[3])-1
                    if end-sta < small_intron_nt_threshold:
                        # a small intron! create a warning message
                        message = "organism: '%s' %s (gff) start: %s end: %s length: %s" % (
                                    org, typeofevidence, (sta+1), end,  (end - sta) )
                        # append the Warning
                        input[org]['warnings'].append( SmallAnnotatedIntronWarning(message) )



    # return the Boolean status
    return HAS_TINY_EXONS
 
# end of function annotatedgeneexonsizeevaluation


def ntracksindnasequencecheck(input,verbose=False):
    """
    Check the size of the annotated coding exons
    
    @type  input: dict
    @param input: input dictionary data structure

    @type  verbose: Boolean
    @param verbose: print discrepancies to STDOUT (True) or not (False, default)

    @rtype:  Boolean
    @return: are N-nucleotides present in the DNA slices? 
    """
    HAS_N_TRACKS = False
    for org in input.keys(): 
        if input[org]['genomeseq'].lower().find("n") > -1:
            HAS_N_TRACKS = True
            # create and print the Warning in verbose mode
            message = "organism/gene: '%s', %s N-nucleotides" % (
                org, input[org]['genomeseq'].lower().count("n") )
            input[org]['warnings'].append( 
                DnaSequenceContainsNsymbolsWarning(message)
                )
            ####################################################
            if verbose: print input[org]['warnings'][-1]
            ####################################################
     
    # return the Boolean status
    return HAS_N_TRACKS 
 
# end of function ntracksindnasequencecheck 


def genestructuregff2orfs(gfftracks,orfset,typeofunigene='n.a.',verbose=False):
    """
    Obtain the Orfs on which a (Uni)Gene structure of gffs is located

    @type  gfftracks: list
    @param gfftracks: list of gff tuple tracks 

    @type  orfset: OrfSet object 
    @param orfset: OrfSet object instance with Orfs of a DNA sequence 

    @type  typeofunigene: string 
    @param typeofunigene: given by unigeneannotation.py (or n.a.)

    @type  verbose: Boolean
    @param verbose: print discrepancies to STDOUT (True) or not (False, default)

    @rtype:  tuple 
    @return: tuple of ( list of orfids, TRACKS_ARE_PROPERLY_MATCHED boolean )
    """
    orfids = []
    TRACKS_ARE_PROPERLY_MATCHED = True
    readingframe = 0

    for pos in range(0,len(gfftracks)):
        gfftrack = gfftracks[pos]
        start, end = int(gfftrack[3])-1, int(gfftrack[4])
        # find the orf to which this thing belongs
        seen = False
        orfsubset = orfset.get_elegiable_orfs(max_orf_start=start+2,min_orf_end=end-3)

        for orf in orfsubset:
            # end offset+4 is needed for exon on the final frontier
            # of the Orfs; nnntaGT means +3 for the TAG and +1 for the splicedonor
            # start offset-2 likewise;  tAGnnnn for the 2nt of the spliceacceptor
            if orf.startPY-2 <= start and orf.endPY+4 >= end:
                dnaseq  = orf.inputgenomicsequence[start:end]
                readingframes = [ readingframe ]

                # do not take a fixed reading frame -> one of previous exons not matched!
                if not TRACKS_ARE_PROPERLY_MATCHED or (typeofunigene in ['utr3p','fragm'] and pos==0):
                    # unigene-specific correction; unigenes of types
                    # 'utr3p','fragm' are protentially NOT starting in frame
                    # to correct this, assume any frame possible and check which
                    # one is the correct one.
                    readingframes = [0,1,2,3]
                for frame in readingframes:
                    protseq = dna2protein(dnaseq,frame)
                    if len(protseq) > 0 and protseq[-1] == '*':
                        # unigenes sometimes have TGA included
                        dnaseq = dnaseq[0:-3]
                        protseq = protseq[0:-1]

                    # find protseq in the ORF's protein sequence
                    inorf = orf.protein_sequence.find(protseq)

                    if inorf > -1:
                        for pos in finditer(protseq,orf.protein_sequence):
                            startDNA = orf.startPY + pos.start()*3 - frame
                            if start == startDNA:
                                orfids.append( orf.id )
                                seen = True
                                #print readingframe, len(dnaseq), len(dnaseq) % 3, "new:", (len(dnaseq) % 3 - readingframe) % 3
                                # update the readingframe variable
                                readingframe = 3 - ( (len(dnaseq) % 3 - frame) % 3 )
                                # break the for-loop
                                break

                        # check if this was the orf of this exon
                        # if so, break the readingframe for-loop
                        if seen: break

                # break if seen
                if seen: break

        if not seen:
            # hmmm....(uni)gene exon not retrievable on an Orf. Why??
            # print orfs that should fit based on coordinates 
            #for orf in orfsubset:
            #    print "MATCHING ORF", orf, orf.frame, start, end,
            #    print orf.inputgenomicsequence[start-2:start+3],
            #    print "..", orf.inputgenomicsequence[end-3:end+2]


            # if here, then this (uni)gene exon is definately not mappable
            if verbose:
                # print the CodingExonNotFoundAsOrfWarning
                print CodingExonNotFoundAsOrfWarning("GFF Track not found as Orf", gfftrack)
            TRACKS_ARE_PROPERLY_MATCHED = False
            # update the readingframe anyway (on UGexon length) 
            start, end = int(gfftrack[3])-1, int(gfftrack[4])
            readingframe = 3 - ( ((end-start) % 3 - readingframe) % 3 )

    # return the result
    return (orfids, TRACKS_ARE_PROPERLY_MATCHED)

# end of function genestructuregff2orfs
 

def potentialsequenceerror(gfftracks,orfset,verbose=False):
    """
    Obtain the Orfs on which a (Uni)Gene structure of gffs is located

    @type  gfftracks: list
    @param gfftracks: list of gff tuple tracks

    @type  orfset: OrfSet object
    @param orfset: OrfSet object instance with Orfs of a DNA sequence

    @type  verbose: Boolean
    @param verbose: print discrepancies to STDOUT (True) or not (False, default)

    @rtype:  tuple
    @return: tuple of ( list of orfids, TRACKS_ARE_PROPERLY_MATCHED boolean )
    """
    readingframe = 0
    orfids = []
    for gfftrack in gfftracks:
        start, end = int(gfftrack[3])-1, int(gfftrack[4])

        if orfset.get_elegiable_orfs(max_orf_start=start+2,min_orf_end=end-3):
            continue

        # get list of orfs that can build up this missing orf
        orfsublist = orfset.get_elegiable_orfs(max_orf_start=end,min_orf_end=start)

        # print orf assembly message
        assembleorf(orfsublist,start,end)

# end of function potentialsequenceerror 


def assembleorf(orflist,start,end,max_errors=2,minimum_errors_only=True):
    """ """
    orflist = order_list_by_attribute(orflist,order_by='length',reversed=True)
    orfidcombis = [] 
    for startpos in range(0,len(orflist)):
        centralorf = orflist[startpos]
        orfset = [ ( centralorf.startPY, centralorf.endPY, centralorf.frame, centralorf ) ]
        # correct 'end' orf gff track position with STOP codon length!
        minpos = min([ tup[0] for tup in orfset ])
        maxpos = max([ tup[1] for tup in orfset ])
        while minpos > start or maxpos < end-3:
            for pos in range(0,len(orflist)):
                #if minpos <= start and maxpos >= end-3: break
                orf = orflist[pos]
                if orf.id in [ tup[-1].id for tup in orfset ]: continue

                if minpos > start and orf.startPY < minpos and orf.endPY >= minpos-3:
                    orfset.insert(0, ( orf.startPY, orf.endPY, orf.frame, orf ) )
                    minpos = min([ tup[0] for tup in orfset ])
                    maxpos = max([ tup[1] for tup in orfset ])
                    # check if this orf even OVERLAPS ths current group of orfs
                    if orf.endPY == maxpos:
                        orfset.append( ( orf.startPY, orf.endPY, orf.frame, orf ) )
                    # goto next iteration
                    break
                if maxpos < end-3 and orf.endPY > maxpos and orf.startPY <= maxpos+3:
                    orfset.append( ( orf.startPY, orf.endPY, orf.frame, orf ) )
                    minpos = min([ tup[0] for tup in orfset ])
                    maxpos = max([ tup[1] for tup in orfset ])
                    # check if this orf even OVERLAPS ths current group of orfs
                    if orf.startPY == minpos:
                        orfset.insert(0, ( orf.startPY, orf.endPY, orf.frame, orf ) )
                    # goto next iteration
                    break
        # calculate minimal amount of sequence changes needed
        frames = [ tup[2] for tup in orfset ]
        distance = _getorfgroupframedistance(frames)
        # only store if not to much errors
        if not minimum_errors_only and distance[0] > max_errors:
            pass
        else:
            if ( sum(distance), distance, orfset ) not in orfidcombis: orfidcombis.append( ( sum(distance), distance, orfset ) )

    # now find the most likely explanation
    orfidcombis.sort()
    minimum_error_count = None    
    for (summed,distance,orfset) in orfidcombis:
        # define minimum number of sequence errors
        if minimum_error_count == None: minimum_error_count=distance[0]
        # break when more than minimum sequence errors
        if minimum_errors_only and distance[0] > minimum_error_count: break
        # if NO sequence errors -> break as well (not an error, but a missed tiny Orf!)
        orfids = [ tup[-1].id for tup in orfset ]
        frames = [ tup[2] for tup in orfset ]

        # oldprinting style, replaced by PotentialSequenceErrorWarning 
        #print orfids, frames, distance
        #for tup in orfset: print tup[-1], tup[-1].frame 

        # print PotentialSequenceErrorWarning
        message = [ str(distance[1:]), "(insertions,deletions,errors)" ]
        for tup in orfset: message.append( str(tup[-1]) )
        print PotentialSequenceErrorWarning( message )

# end of function assembleorf 


def _getorfgroupframedistance(frames):
    """
    @rtype:  tuple
    @return: (total,insertions,deletions,errors)
    """
    insertions, deletions, errors = 0, 0, 0
    for pos in range(1,len(frames)):
        if frames[pos-1] == frames[pos]: errors+=1
        elif frames[pos]-frames[pos-1] == 1: deletions+=1
        else: insertions+=1
    return (len(frames)-1,insertions,deletions,errors)

# end of function _getorfgroupframedistance  


def readsequences(input):
    """
    Read sequence files into input dict

    @type  input: dict
    @param input: input dictionary data structure

    @rtype:  dict
    @return: incremented input dictionary data structure
    """
    for org in input.keys():

        # (A) open protein sequence file
        if input[org]['proteinseqfile']:
            input[org]['proteinseq'] = "".join(open(input[org]['proteinseqfile']
                ).readlines()[1:]).replace("\n","")
        else:
            input[org]['proteinseq'] = ""
        # (B) open DNA sequence file
        input[org]['genomeseq']  = "".join(open(input[org]['genomeseqfile']
                ).readlines()[1:]).replace("\n","")

    # return the incremented input dict
    return input

# end of function readsequences


def rungetorf(input):
    """
    Run getorf and parse the output file in Orf objects for the input DNA sequences

    @type  input: dict
    @param input: input dictionary data structure

    @rtype:  dict
    @return: incremented input dictionary data structure
    """
    for org in input.keys():

        # check if sequence files are alrady read
        if not input[org]['proteinseq']:
            if input[org]['proteinseqfile']: 
                input[org]['proteinseq'] = "".join(open(input[org]['proteinseqfile']
                   ).readlines()[1:]).replace("\n","")
        if not input[org]['genomeseq']:
            input[org]['genomeseq']  = "".join(open(input[org]['genomeseqfile']
                   ).readlines()[1:]).replace("\n","")

        ## (C) run& parse getorf
        if not input[org].has_key('orfs') or not input[org]['orfs']:
            orfsetobject = GetorfOrfSet(sequence=input[org]['genomeseq'])
            orfsetobject.accession = input[org]['FREF']
            orfsetobject.getorfs()
            input[org]['orfs'] = orfsetobject 

    # return the incremented input dict
    return input

# end of function rungetorf


def parseinputgff(input):
    """
    Parse all the input gff file(s) into the tracks of the genes and the unigenes

    @type  input: dict
    @param input: input dictionary data structure

    @rtype:  dict
    @return: incremented input dictionary data structure
    """
    for org in input.keys():

        # (A) set dict keys for gene, unigene and etcetera
        for key in ['gffs','gff-gene','gff-unigene']:
            if not input[org].has_key(key):
                input[org][key] = []

        # (B) parse input gff file
        if input[org].has_key('gff') and input[org]['gff']:
            # NEW since 01/07/2009; when a GeneLocusDirectoy is applied,
            # coordinates must be recalculated from absolute to 1-based
            # furthermore, negative strand coords must be mirrored
            offset = 0
            strand = '+'

            if input[org].has_key('locusfile'):
                has_locusfile = True
                locusGfft = open(input[org]['locusfile']).readlines()[0].split("\t")
                offset    = int(locusGfft[3]) - 1 # GGB/GFF 1 coord == Python 0 coord
                length    = int(locusGfft[4]) - int(locusGfft[3]) + 1
                strand    = locusGfft[6]
            else:
                has_locusfile = False


            if input[org]['gff-gene'] or input[org]['gff-unigene']:
                pass
            else:
                # get gfflines of both gene and unigene structure
                gfflines = parsegfffile( input[org]['gff'], offset=offset )
                if input[org].has_key('unigenefile') and input[org]['unigenefile']:
                    gfflines.extend( parsegfffile( input[org]['unigenefile'], offset=offset ) )

                # loop over the lines and place all where they belong
                for gff in gfflines:
                    column9keys = [ c9part.split(" ")[0] for c9part in gff[8].split("; ") ]
                    if GFF_UNIGENE_GCLASS == column9keys[0]:
                        # track belonging to an unigene
                        input[org]['gff-unigene'].append( gff )
                    elif GFF_GENE_GCLASS == column9keys[0]:
                        # track belonging to current gene's annotation
                        input[org]['gff-gene'].append( gff )
                    elif GFF_UNIGENE_GCLASS in column9keys:
                        # track belonging to an unigene
                        input[org]['gff-unigene'].append( gff )
                    elif GFF_GENE_GCLASS in column9keys:
                        # track belonging to current gene's annotation
                        input[org]['gff-gene'].append( gff )
                    elif gff[2].lower() in GFF_GENE_ELEMENT_FMETHOD_LIST and\
                    Set(GFF_TRANSCRIPT_GCLASS_LIST).intersection(column9keys):
                        # track belonging to current gene's annotation
                        input[org]['gff-gene'].append( gff )
                    else:
                        # other track..... Store to general gff's
                        input[org]['gffs'].append(gff)
                        pass

                # Check if we can/need to filter the unigene GFFs.
                # Required because my (Avdb) ABFGP/dbwarehouse implementation
                # In exceptional cases (~1:10.000), a unigene can be aligned
                # on multiple genomic loci. All these loci might be present
                # in the unigene file(s) in the GeneLocusDirectory. And this
                # will cause errors here. So, filter out these exceptions.
                if input[org]['gff-unigene'] and has_locusfile and\
                len(Set([gffTup[0] for gffTup in input[org]['gff-unigene']]))>1:
                    passedGff = []
                    coords = Set(range(0,length))
                    for gffTup in input[org]['gff-unigene']:
                        if gffTup[0] != locusGfft[0]:
                            continue
                        if len(coords.intersection(range(gffTup[3],gffTup[4])))==0:
                            continue
                        # if here -> append to passedGff
                        passedGff.append(gffTup)
                    # overwrite input[org]['gff-unigene']
                    input[org]['gff-unigene'] = passedGff

                # if strand is negative -> mirror the gff coordinates!
                if strand == '-':
                    input[org]['gff-gene']    = mirrorgfflist( input[org]['gff-gene'], length=length )
                    input[org]['gff-unigene'] = mirrorgfflist( input[org]['gff-unigene'], length=length )
                    input[org]['gffs']        = mirrorgfflist( input[org]['gffs'], length=length )

                # order the gff lists
                input[org]['gff-gene']    = order_gff_list(input[org]['gff-gene'])
                input[org]['gff-unigene'] = order_gff_list(input[org]['gff-unigene'])
                input[org]['gffs']        = order_gff_list(input[org]['gffs'])

    # return the incremented input dict
    return input

# end of function parseinputgff


def confirmcanonicalsplicesites(sequence,gfflist,exon_fmethod=None,verbose=False):
    """
    """
    onlycanonical = True
    warnings = []
    exons = filtergffs4fmethod(gfflist,exon_fmethod)
    for i in range(1,len(exons)):
        donor,accep = exons[i-1:i+1]
        dSeq = sequence[donor[4]:donor[4]+2].upper()
        aSeq = sequence[accep[3]-3:accep[3]-1].upper()
        if dSeq not in CANONICAL_DONOR_SITES:
            onlycanonical = False
            if dSeq not in NON_CANONICAL_DONOR_SITES:
                warnings.append( ImplausibleSpliceSiteWarning("Donor:"+dSeq) )
            else:
                warnings.append( NonCanonicalSpliceSiteWarning("Donor:"+dSeq) )
            if verbose: print warnings[-1]
        if aSeq not in CANONICAL_ACCEPTOR_SITES:
            warnings.append( ImplausibleSpliceSiteWarning("Acceptor:"+aSeq) )
            onlycanonical = False
            if verbose: print warnings[-1]

    # return the status of the observed splice sites
    return onlycanonical, warnings

# end of function confirmcanonicalsplicesites


def genestructurehasproperstartcodon(gfflist,dnasequence,fmethod,verbose=False):
    """
    """
    GENE_HAS_PROPER_START_CODON = False
    startcodon = []
    for gfftrack in gfflist:
        if gfftrack[2] == fmethod and len(startcodon)==0:
            if gfftrack[4] - gfftrack[3] +1 >= 3:
                startcodon = range(gfftrack[3]-1,gfftrack[3]-1+3)
                break
            else:
                startcodon = range(gfftrack[3]-1,gfftrack[4])
        elif gfftrack[2] == fmethod and len(startcodon)<3:
            offset = 0
            while len(startcodon)<3:
                pos = gfftrack[3]-1+offset
                startcodon.append(pos)
                offset+=1
        else:
            pass

    # now check the startcodon
    if "".join([dnasequence[pos] for pos in startcodon]).lower() == 'atg':
        GENE_HAS_PROPER_START_CODON = True

    if not GENE_HAS_PROPER_START_CODON and verbose and gfflist:
        print "Warning (uni)gene '%s' of sequence '%s' has no proper start codon (%s %s)" % (
                gfflist[0][8].split(";")[0],
                gfflist[0][0],
                "".join([dnasequence[pos] for pos in startcodon]).lower(),
                startcodon
                )

    # return the status
    return GENE_HAS_PROPER_START_CODON

# end of function genestructurehasproperstartcodon


def genestructurehasproperstopcodon(gfflist,dnasequence,fmethod,verbose=False):
    """
    """
    GENE_HAS_PROPER_STOP_CODON = False
    stopcodon = []
    for trackpos in range(len(gfflist)-1,-1,-1):
        gfftrack = gfflist[trackpos]
        if gfftrack[2] == fmethod and len(stopcodon)==0:
            if gfftrack[4] - gfftrack[3] +1 >= 3:
                stopcodon = range(gfftrack[4]-3,gfftrack[4])
                #if verbose:
                #    print "stopcodon (1)", stopcodon, fmethod,
                #    print "".join([dnasequence[pos] for pos in stopcodon])
                break
            else:
                stopcodon = range(gfftrack[3]-1,gfftrack[4])
                #if verbose: 
                #    print "stopcodon (2)", stopcodon, fmethod,
                #    print "".join([dnasequence[pos] for pos in stopcodon])
        elif gfftrack[2] == fmethod and len(stopcodon)<3:
            # the stop-codon is spliced !!!!
            offset = 0
            while len(stopcodon)<3:
                pos = gfftrack[4]-offset
                stopcodon.insert(0,pos)
                offset-=1
            #if verbose: 
            #    print "stopcodon (3)", stopcodon, fmethod,
            #    print "".join([dnasequence[pos] for pos in stopcodon])
        else:
            pass

    # now check the stop codon
    if "".join([dnasequence[pos] for pos in stopcodon]).lower() in ['tga','taa','tag']:
        GENE_HAS_PROPER_STOP_CODON = True

    if not GENE_HAS_PROPER_STOP_CODON and verbose and gfflist:
        for track in gfflist: print track
        print "Warning (uni)gene '%s' of sequence '%s' has no proper stop codon (%s %s)" % (
                gfflist[0][8].split(";")[0],
                gfflist[0][0],
                "".join([dnasequence[pos] for pos in stopcodon]).lower(),
                stopcodon
                )

    # return the status
    return GENE_HAS_PROPER_STOP_CODON

# end of function genestructurehasproperstopcodon


def create_gene_utrs(gene_gff_list):
    """
    Create UTR tracks for an (annotated) gene

    @type  gene_gff_list: list
    @param gene_gff_list: list with gff tuples of the (annotated) gene

    @rtype  utrs: list
    @return utrs: list with gff tuples of the UTR tracks

    @attention: requires global variables GFF_CDS_FMETHOD, GFF_EXON_FMETHOD
    @attention: requires global variables GFF_UTR5_FSOURCE, GFF_UTR5_FMETHOD
    @attention: requires global variables GFF_UTR3_FSOURCE, GFF_UTR3_FMETHOD
    """
    # make sets of unigene coordinates
    cdscoords  = gffs2coordset(gene_gff_list,fmethod=[GFF_CDS_FMETHOD])
    exoncoords = gffs2coordset(gene_gff_list,fmethod=[GFF_EXON_FMETHOD])

    # return list with UTR tracks
    utrs = []
    
    # create a list with 5'UTR coordinates
    if cdscoords and exoncoords:
        utr5p_coords = []
        utr3p_coords = []

        # find start codon in cdscoords: min()
        for coord in range( min(exoncoords), min(cdscoords) ):
            if coord in exoncoords:
                utr5p_coords.append(coord)

        # find stop codon in exoncoords: max()
        if max(cdscoords)+3 != max(exoncoords):
            # check if exon-end != cds-end+3 -> only a stop codon
            for coord in range( max(cdscoords)+1, max(exoncoords)+1 ):
                if coord in exoncoords:
                    utr3p_coords.append(coord)

        if utr5p_coords or utr3p_coords:
            # get list with coding exons for track backbone
            cexons = filtergffs4fmethod(gene_gff_list,fmethod=GFF_CDS_FMETHOD)

            # create 5'UTRs
            gfftrack    = list( cexons[0] )
            gfftrack[1] = GFF_UTR5_FSOURCE
            gfftrack[2] = GFF_UTR5_FMETHOD
            utrs.extend( coordset2gfftracks(utr5p_coords,gfftrack))

            # create 3'UTRs
            gfftrack    = list( cexons[0] )
            gfftrack[1] = GFF_UTR3_FSOURCE
            gfftrack[2] = GFF_UTR3_FMETHOD
            utrs.extend( coordset2gfftracks(utr3p_coords,gfftrack))

    # return the utr tracks
    return utrs

# end of function create_gene_utrs


def geneconfirmation(input,verbose=False):
    """
    Confirm given gene GFF structure with Orfs of the DNA sequence

    @type  input: dict
    @param input: input dictionary data structure

    @type  verbose: Boolean
    @param verbose: print discrepancies to STDOUT (True) or not (False, default)

    @rtype:  tuple 
    @return: tuple of ( input data structure, TRACKS_ARE_PROPERLY_MATCHED boolean )
    """
    input = readsequences(input)
    input = rungetorf(input)
    input = parseinputgff(input)

    GENE_TRACKS_ARE_PROPERLY_MATCHED = True
    for org in input.keys():

        # (A) find the genestructure_orfmodel of the annotated/known gene structure
        input[org]['orfid-genestructure'] = []
        _GENE_TRACKS_ARE_PROPERLY_MATCHED = False
        if input[org]['gff-gene']:
            # get only the CDS-type of tracks that define the coding sequence
            genecdstracks = filtergffs4fmethod( input[org]['gff-gene'], GFF_CDS_FMETHOD ) 

            # for comparison / later use get exon tracks too
            geneexontracks = filtergffs4fmethod( input[org]['gff-gene'], GFF_EXON_FMETHOD )

            if genecdstracks:
                # as expected...
                pass
            elif not genecdstracks and geneexontracks:
                # MANY GeneLoci miss CDS tracks but have exon tracks
                # In the BROAD annotation logic this represents partial genes:
                # based on similarity, but missing start and/or stop codon
                # In practice, these instances are likely the result of
                # sequence errors!
                # But, back to the problem we face here: what to do with these cases?
                # The choice that was made here is as follows:
                # - check the exon tracks as if they represent unigenes
                # - if OK -> convert the exon tracks to CDS tracks
                # - the GENE_TRACKS_ARE_PROPERLY_MATCHED will later in this
                #   function we recognized as False 

                # Append the Warning message
                message = "no tracks of type '%s'" % (GFF_CDS_FMETHOD)
                warn = IncompleteGeneStructureWarning(message)
                input[org]['warnings'].append( warn )

                # obtain exon track's DNA sequence
                command = """%s %s %s | grep -v ">" | tr -d "\\n" """ % ( PYTHON_PATH,
                                         EXECUTABLE_GFF2FASTA,
                                         input[org]['genomeseqfile'] )
                ci,co,ce = os.popen3(command)
                # make shure geneexontracks FREF matches input[org]['genomeseqfile'] FREF
                loci_dna_fref = open(input[org]['genomeseqfile']).readlines()[0].split(' ')[0].replace(">","") 
                ci.write( gffs2txt(geneexontracks).replace( geneexontracks[0][0],loci_dna_fref ) )
                ci.close()
                protseq = dna2protein( co.read().strip(), 0 )
                co.close()
                error = ce.read()
                ce.close()

                if error.find("unrecognized header:") > -1:
                    # rewrite GFF fref to genomic fref
                    print "XXXX",error
                    print command
                    print gffs2txt(geneexontracks)
                    print input[org]['genomefref'],input[org]['locusfref'],input[org]['proteinfref']
                    print open( input[org]['genomeseqfile']).readlines()[0]

                if protseq.count('*') == 0:
                    # The EXON type tracks have an ORF in frame 0
                    # This represents a partial gene structure.
                    # Copy the EXON tracks as CDS tracks and run
                    # the rest of this function.
                    for _track in geneexontracks:
                        track = list(_track)
                        track[2] = GFF_CDS_FMETHOD
                        track = tuple(track)
                        input[org]['gff-gene'].append( track )
                        genecdstracks.append( track ) 
                else:
                    message = "tracks of type '%s' have no ORF" % (GFF_EXON_FMETHOD)
                    warn = IncompleteGeneStructureWarning(message)
                    input[org]['warnings'].append( warn )
                    GENE_TRACKS_ARE_PROPERLY_MATCHED = False
                    continue

                ################################################################
                if verbose:
                    print org, "NO CDS TRACKS!!", input[org]['FREF']
                    print protseq[0:20]+"..."+protseq[-20:], len(protseq)
                    for track in genecdstracks: print track[0:7]
                ################################################################

            else:
                # neither types of tracks -> no gene annotation at all! 
                GENE_TRACKS_ARE_PROPERLY_MATCHED = False
                continue

            # correct the final exon for absence of STOP-codon
            # many given gene annotations do not include STOP-codon in final exon, some do
            # By definition:
            # A track of type CDS  should EXCLUDE the stop codon triplet itself
            # A track of type exon should INCLUDE the stop codon triplet itself
            # Not all annotations follow this principle, however...
            # And, in ABFGP, it was chosen to use ANNOTATED CDS tracks as input and
            # PREDICTED EXON tracks as output.
            # This was done to have the stop-codon included as a check-final-check for
            # the final exon (if the predicted CDS is indeed followed by a stop codon). 

            if not genestructurehasproperstopcodon(genecdstracks,
            input[org]['genomeseq'],GFF_CDS_FMETHOD,verbose=False):
                # Well, that is as expected -> we expect here CDS tracks
                # check if the NEXT triplet isa STOP-codon
                try:
                    triplet = input[org]['genomeseq'][ genecdstracks[-1][4] : genecdstracks[-1][4]+3 ].lower()
                except:
                    # here errors have occurred a lot! Some annotations mix up exon/CDS
                    # tracks (and in some cases I mixed them up during the conversion
                    # proces from annotation to GeneLoci). Print this exception ALWAYS!
                    print "Exception occurred:"
                    print org, "genecdstracks:", len(genecdstracks)
                    triplet = input[org]['genomeseq'][ genecdstracks[-1][4] : genecdstracks[-1][4]+3 ].lower()

                if triplet in ['tga','tag','taa']:
                    # update track by +3 length for STOP-codon
                    # this update is only valid within the scope of this function;
                    # input[org]['gff-gene'] remains unchanged.
                    track = list( genecdstracks[-1] )
                    track[4]+=3
                    genecdstracks[-1] = tuple(track)
                else:
                    # no CDS type of track followed by a stop-codon;
                    # ignore here, error message will follow later
                    pass

            elif genecdstracks and geneexontracks and\
            genecdstracks[-1][3] == geneexontracks[-1][3] and\
            genecdstracks[-1][4]+3 == geneexontracks[-1][4]:
                # do not correct anything; it's fine!
                pass

            else:
                # hmmm.... that is weird! We expect here CDS tracks, not EXON tracks
                # do a HARD correction on the provided annotation track!
                # TODO: this piece of code will fail for spliced stop codons:
                # Tgt.....agGA because it is spread over >1 track            
                for pos in range(0,len(input[org]['gff-gene'])):
                    track = input[org]['gff-gene'][pos]
                    if track == genecdstracks[-1]:
                        # HARD-replace in the input[org]['gff-gene'] list with gff tuples!
                        track = list(track)
                        track[4] = track[4]-3
                        input[org]['gff-gene'][pos] = tuple(track)
                        break
                else:
                    print "GeneStructureWarning: VERY UNUSUAL SPLITTED STOP CODON", org

            # convert exon/CDS tracks to UTR tracks
            utrs = create_gene_utrs(input[org]['gff-gene'])
            input[org]['gff-gene'].extend(utrs)

            # get the orfs of the genestructure
            orfids,_GENE_TRACKS_ARE_PROPERLY_MATCHED = genestructuregff2orfs(genecdstracks,input[org]['orfs'],verbose=False)
            if (not _GENE_TRACKS_ARE_PROPERLY_MATCHED or not orfids):
                print "_GENE_TRACKS_ARE_PROPERLY_MATCHED == False", org, orfids
                # redo just for logging purposes!
                orfids,_GENE_TRACKS_ARE_PROPERLY_MATCHED = genestructuregff2orfs(genecdstracks,input[org]['orfs'],verbose=True)

            # place the list with orfids in 'orfid-genestructure' dict key 
            input[org]['orfid-genestructure'] = orfids
            # check if the gene was properly matched.
            if not _GENE_TRACKS_ARE_PROPERLY_MATCHED or not input[org]['orfid-genestructure']:
                GENE_TRACKS_ARE_PROPERLY_MATCHED = False
                # append warning to warnings list
                warn = GeneStructureIsNotMappableOnOrfsWarning("not mappable on orfs")
                input[org]['warnings'].append( warn )
                # set known gene-tracks to empty list!
                input[org]['orfid-genestructure'] = []
                # check for a potential sequence error in the sequence
                potentialsequenceerror(genecdstracks,input[org]['orfs'],verbose=verbose)
            else:
                # all fine. Label all the **known** Orfs with IS_ANNOTATED_EXON_LABEL
                for orfid in input[org]['orfid-genestructure']:
                    orfObj = input[org]['orfs'].get_orf_by_id(orfid)
                    # add label to Orf object
                    setattr(orfObj,IS_ANNOTATED_EXON_LABEL,True)


        # vars for annotated start & stop codon check
        GENE_HAS_PROPER_START_CODON = False
        GENE_HAS_PROPER_STOP_CODON  = False

        # (B) check start codon of the gene structure
        if _GENE_TRACKS_ARE_PROPERLY_MATCHED and input[org]['gff-gene']:
            GENE_HAS_PROPER_START_CODON = genestructurehasproperstartcodon(
                    genecdstracks,
                    input[org]['genomeseq'],
                    GFF_CDS_FMETHOD,
                    verbose=verbose)
            # create warning message if a problem occured 
            if not GENE_HAS_PROPER_START_CODON:
                # create warning message
                try:
                    sta,end  = genecdstracks[0][3]-1, genecdstracks[0][3]+2
                    startdna = input[org]['genomeseq'][sta:end]
                    orf      = input[org]['orfs'].get_orf_by_id(
                                   input[org]['orfid-genestructure'][0] )
                    message  = "no proper start codon (%s [%s] %s) on orf %s" % (
                        startdna, dna2protein(startdna,0),range(sta,end), orf)
                except:
                    # coordinates are horribly wrong -> make simple warning message
                    message = "no proper start codon (Exception)"
                # append warning to warnings list
                warn = IncompleteGeneStructureWarning(message)
                input[org]['warnings'].append( warn )


        # (C) check stop codon of the gene structure
        if _GENE_TRACKS_ARE_PROPERLY_MATCHED and input[org]['gff-gene']:
            GENE_HAS_PROPER_STOP_CODON = genestructurehasproperstopcodon(
                    genecdstracks,
                    input[org]['genomeseq'],
                    GFF_CDS_FMETHOD,
                    verbose=verbose)
            if not GENE_HAS_PROPER_STOP_CODON:
                # create warning message
                try:
                    sta,end = genecdstracks[-1][4]-3, genecdstracks[-1][4]
                    stopdna = input[org]['genomeseq'][sta:end]
                    message = "no proper stop codon (%s [%s] %s) on orf %s" % (
                        stopdna, dna2protein(stopdna,0),range(sta,end),
                        input[org]['orfs'].get_orf_by_id(input[org]['orfid-genestructure'][-1])
                        )
                except:
                    # coordinates are horribly wrong -> make simple warning message
                    message = "no proper stop codon (Exception)"
                # append warning to warnings list
                warn = IncompleteGeneStructureWarning(message)
                input[org]['warnings'].append( warn )

        # (D) check if start_codon is present as annotated track
        if _GENE_TRACKS_ARE_PROPERLY_MATCHED and GENE_HAS_PROPER_START_CODON and\
        input[org]['gff-gene'] and not filtergffs4fmethod( input[org]['gff-gene'],
        GFF_GENESTART_FMETHOD ):
            # create start_codon track
            # TODO: write a generic function that does the job,
            # TODO: whick takes splices start codons into account
            startcodon    = list(genecdstracks[0])
            startcodon[2] = GFF_GENESTART_FMETHOD
            startcodon[4] = startcodon[3]+2
            startcodon[5] = "."
            startcodon[7] = "."
            startcodon[8] = startcodon[8].split("; ")[0]
            input[org]['gff-gene'].append( tuple(startcodon) )
            ################################################################
            if verbose:
                sta,end  = genecdstracks[0][3]-1, genecdstracks[0][3]+2
                startdna = input[org]['genomeseq'][sta:end]
                print "CREATED:", genecdstracks[0][0:7]
                print "CREATED:", tuple(startcodon)
                print "CREATED:", sta, end, startdna, dna2protein(startdna,0)
            ################################################################

        # (E) check if stop_codon is present as annotated track
        if _GENE_TRACKS_ARE_PROPERLY_MATCHED and GENE_HAS_PROPER_STOP_CODON and\
        input[org]['gff-gene'] and not filtergffs4fmethod( input[org]['gff-gene'],
        GFF_GENESTOP_FMETHOD ):
            # create stop_codon track
            stopcodon    = list(genecdstracks[-1])
            stopcodon[2] = GFF_GENESTOP_FMETHOD
            stopcodon[4] = stopcodon[4]
            stopcodon[3] = stopcodon[4]-2 
            stopcodon[5] = "."
            stopcodon[7] = "."
            stopcodon[8] = stopcodon[8].split("; ")[0]
            input[org]['gff-gene'].append( tuple(stopcodon) )
            ################################################################
            if verbose:
                sta,end = genecdstracks[-1][4]-3, genecdstracks[-1][4]
                stopdna = input[org]['genomeseq'][sta:end]
                print "CREATED:", genecdstracks[-1][0:7]
                print "CREATED:", tuple(stopcodon)
                print "CREATED:", sta, end, stopdna, dna2protein(stopdna,0)
            ################################################################



        # (F) label the first & final (annotated) Orfs in the gene structure
        if GENE_HAS_PROPER_START_CODON and input[org]['orfid-genestructure']:
            firstorfid = input[org]['orfid-genestructure'][0]
            firstorf = input[org]['orfs'].get_orf_by_id(firstorfid)
            # add label to Orf object
            setattr(firstorf,FIRST_ANNOTATED_EXON_LABEL,True)
        if GENE_HAS_PROPER_STOP_CODON and input[org]['orfid-genestructure']:
            finalorfid = input[org]['orfid-genestructure'][-1]
            finalorf = input[org]['orfs'].get_orf_by_id(finalorfid)
            # add label to Orf object
            setattr(finalorf,FINAL_ANNOTATED_EXON_LABEL,True)


        # check if START & STOP codon are properly matched
        if not GENE_HAS_PROPER_START_CODON or not GENE_HAS_PROPER_STOP_CODON:
            _GENE_TRACKS_ARE_PROPERLY_MATCHED = False
            GENE_TRACKS_ARE_PROPERLY_MATCHED  = False


        # print the negative outcome in verbose mode
        if not _GENE_TRACKS_ARE_PROPERLY_MATCHED and verbose:
            # set known gene-tracks to empty list!
            input[org]['orfid-genestructure'] = []
            print "# WARNING: Gene is not mappable on Orfs:", org

    # return the incremented input dict and error status (True or False)
    return input, GENE_TRACKS_ARE_PROPERLY_MATCHED

# end of function geneconfirmation


if __name__ == "__main__":
    # AbgpGeneLocus and AbgpGeneLocusDirectory imports
    from lib_optparse import (
            abgpoptparser,
            genelocustatsoptions,
            validate_genelocustatsoptions
            )
    from abgpgenelocusdirectory import AbgpGeneLocusDirectory
    from lib_sequencerepetitiveness import annotatedproteinsequencerepetitivenesscheck

    # construct the command line option parser
    parser = abgpoptparser()
    genelocustatsoptions(parser)
    parser.remove_option('--explain')
    # parse the command line & validate
    (OPTIONS, args) = parser.parse_args()
    validate_genelocustatsoptions(parser,OPTIONS)

    # change directory to the one specified by optparse option
    os.chdir(OPTIONS.DIRECTORY)
    locus = AbgpGeneLocusDirectory(OPTIONS.DIRECTORY)
    input = locus.toinputdict()
 
    # do geneconfirmation
    input,gene_status = geneconfirmation(input,verbose=True)

    # do a protein sequence repetitiveness check
    IS_REPETITIVE = annotatedproteinsequencerepetitivenesscheck(input)

    # do a check on tiny/small exons
    HAS_SMALL_OR_TINY_EXONS = annotatedgeneexonsizeevaluation(input)

    # do check on noncannocical splice sites
    for org in input.keys():
        status, warnings = confirmcanonicalsplicesites(input[org]['genomeseq'],
                input[org]['gff-gene'],exon_fmethod=GFF_CDS_FMETHOD,verbose=False)
        if not status:
            print "# WARNING: Non-canonical splice site(s) for organism/fref:", org, input[org]['FREF']
            for warn in warnings: print warn

    # go back to OPTIONS.CURRENTWORKINGDIR
    os.chdir(OPTIONS.CURRENTWORKINGDIR)

    if OPTIONS.verbose:
        print "# EOF geneconfirmation; if nothing printed, no abnormalities observed"


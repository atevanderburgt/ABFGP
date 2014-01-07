"""
AbgpGeneLocusDirectory class
"""

# Module metadata
__authors__ = "Ate van der Burgt"
__license__ = "MIT"

# Abgp Imports
from abgp_warnings import *
from gff import (
    parsegfffile, mirrorgfflist,
    gffs2txt, gff2db, dbcleanup,
    filtergffs4fmethods,
    exons2introns, filtergffs4fmethod,
    gffs2coordset, getseqbyfeaturetuple,
    apply_fref, gnamefromgfftuple,
    )
from abgp_geneconfirmation import (
    parsegfffile,
    mirrorgfflist,
    rungetorf,
    parseinputgff,
    geneconfirmation,
    )
from abgp_unigeneconfirmation import (
    correct_unigene_for_utrs,
    get_annotated_genes_start_codon,
    get_annotated_genes_stop_codon,
    )

from pythonlibs.uniqueness import get_random_string_tag, NONLETTERS_TRANS
from lib_tcode import obtaintcodedata 

from gene.start import (
    IC_TSS_PATTERN_OFFSET,
    TranslationalStartSite,
    score_tss,
    scan_pssm_tss,
    )
from gene.exon import (
    FirstExonOnOrf,
    FinalExonOnOrf,
    SingleExonOnOrf,
    ExonOnOrf,
    )
from gene.intron import (
    IntronConnectingOrfs,
    )
from gene.donor import (
    SpliceDonor,
    IC_DONOR_PATTERN_OFFSET,
    )
from gene.acceptor import (
    SpliceAcceptor,
    IC_ACCEPTOR_PATTERN_OFFSET,
    )
from gene.splicesite import (
    get_shared_nucleotides_at_splicesite,
    scan_pssm_splice_site,
    _score_splice_site,
    )
from gene.signalpeptide import SignalPSignalPeptide

from lib_fasta import (
    parseFasta,
    parseSingleFasta,
    parseSingleFastaHeaderFromFile,
    writeSingleFasta
    )
from lib_filevalidation import (
    fastafilesequencelength,
    gfffilesequencelength,
    )

# Python Imports
import os
from os import popen3 as osPopen3
from sets import Set
from copy import deepcopy

# Global variable Imports
from settings.abgp import ABGP_VERSION
from settings.dbwarehouse import ABGP_LOCUS_IDENTIFIER2ORGANISM_MAPPING
from settings.gff import (
    GFF_EXON_FMETHOD,
    GFF_CDS_FMETHOD,
    )
from settings.gff import (
    GFF_UGEXON_FMETHOD,
    GFF_UG3UTREXON_FMETHOD,
    GFF_UG5UTREXON_FMETHOD,
    GFF_TSS_FSOURCE,
    GFF_TSS_FMETHOD,
    GFF_SIGNALP_FSOURCE,
    GFF_SIGNALP_FMETHOD,
    GFF_TMHMM_FSOURCE,
    GFF_TMHMM_FMETHOD,
    )
from settings.executables import (
    EXECUTABLE_SIGNALP_VERSION,
    SIGNALP_MIN_D,
    SIGNALP_MIN_SPROB,
    EXECUTABLE_GFF2FASTA,
    EXECUTABLE_UNIGENEANNOTATION,
    EXECUTABLE_TRANSEQ,
    PYTHON_PATH,
    )
from settings.splicesites import (
    CANONICAL_DONOR_SITES,
    NON_CANONICAL_DONOR_SITES,
    CANONICAL_ACCEPTOR_SITES,
    )
from settings.translationalstartsites import (
    TSS_MIN_PSSM_SCORE
    )
from settings.ggb import (
    DB_NAME,DB_HOST,DB_CONNECTION,
    DB_USER,DB_PASS
    )


class AbgpGeneLocusDirectory:
    """ """
    def __init__(self,dirname):
        """ """
        if not IsAbgpGeneLocusDirectory(dirname):
            message = "directory '%s' is not a AbgpGeneLocusDirectory" % dirname
            raise NoAbgpGeneLocusDirectory, message
        # attributes for file pointers in locusdir
        self.locus_file   = None
        self.dna_file     = None
        self.gene_file    = None
        self.protein_file = None
        self.unigene_file = None

        # set absolute path of directory to self.dirname
        self.dirname = os.path.abspath(dirname)

        # attribute that can hold (miscellaneous) gff data
        self._gffdata = []

        # attributes that can hold (pre-)predicted gene components
        self._locus_objects_tss     = []
        self._locus_objects_tmhmm   = []
        self._locus_objects_signalp = []

        # Attribute that can contain an partial `input dict`
        # input dict is a dictionary data structure used in ABGP
        # self.input is the (sub)dictionary that is as
        # in the main ABGP algorithm located in input[FREF], where
        # FREF is the pointer to this organism/gene identifier
        self.input = {}

        # abstract filename(s) from locus directory
        self.abstract_locusfile()     # stored in self.locus_file
        self.abstract_dnafiles()      # stored in self.dna_file
        
        # get the header of the dna_file and use this as fref
        fref,descr = parseSingleFastaHeaderFromFile(os.path.join(self.dirname,self.dna_file))
        self.fref = fref

        # abstract gene's annotation
        self.abstract_genefiles()     # stored in self.gene_file
        self.abstract_proteinfiles()  # stored in self.protein_file
        self.abstract_unigenefiles()  # stored in self.unigene_file

        # check if locus gff file specifies the same length as the dna file is
        if fastafilesequencelength(self.dna_file) != gfffilesequencelength(self.locus_file):
            message = "dnafile '%s' (l=%s) and locusfile '%s' (l='%s') specify unidentical sequence lengths" % (
                    self.dna_file, fastafilesequencelength(self.dna_file),
                    self.locus_file, gfffilesequencelength(self.locus_file)
                    )
            raise UnequalGeneLocusLengths, message

        # get the header of the dna_file and use this as fref
        fref,descr = parseSingleFastaHeaderFromFile(os.path.join(self.dirname,self.dna_file))
        self.fref = fref

    # end of function __init__


    def __str__(self):
        """ """
        return "<%s %s %s>" % (self.__class__.__name__, self.fref, self.dirname)

    # end of function __str__

    def _files_unigene(self):
        """ """
        return [os.path.join(os.path.abspath(self.dirname),x) for x in os.listdir(self.dirname) if x.endswith('.gff') and x.find('.unigene') > 0]
    # end of function _files_unigene

    def _files_protein(self):
        """ """
        return [os.path.join(os.path.abspath(self.dirname),x) for x in os.listdir(self.dirname) if x.endswith('.fa') and x.find('.protein') > 0]
    # end of function _files_protein

    def _files_gene(self):
        """ """
        return [os.path.join(os.path.abspath(self.dirname),x) for x in os.listdir(self.dirname) if x.endswith('.gff') and x.find('.gene') > 0]
    # end of function _files_gene

    ############################################################################
    ### Functions for gene accessibility
    ############################################################################

    def abstract_genefiles(self):
        """ """
        # TMP TODO HACK: there can be more than 1 gene gff file in a locus directory!
        try:
            self.gene_file = self._files_gene()[0]
        except IndexError:
            self.gene_file = None 
    # end of function abstract_genefiles

    def _obtain_gene_gff(self,input={}):
        """
        @type  input: dict
        @param input: SINGLE SUBDICT of input dictionary data structure
        """
        # handle potentially applied input argument
        self._handle_input_subdict(input)
        if self.input and self.input.has_key('gff-gene') and self.input['gff-gene']:
            return self.input['gff-gene']
        elif self.gene_file:
            # read gff file to list of tuples
            locusgff   = parsegfffile(self.locus_file)
            # define offset and length from locusgff
            offset = locusgff[0][3] - 1
            length = locusgff[0][4] - locusgff[0][3] + 1
            genegff= parsegfffile(self.gene_file,offset=offset)
            if locusgff[0][6] == '-':
                genegff    = mirrorgfflist( genegff, length=length)
            return genegff
        else:
            return []

    # end of function _obtain_gene_gff


    ############################################################################
    ### Functions for protein accessibility
    ############################################################################

    def abstract_proteinfiles(self):
        """ """
        # TMP TODO HACK: there can be more than 1 protein file in a locus directory!
        try:
            self.protein_file = self._files_protein()[0]
        except IndexError:
            self.protein_file = None 
    # end of function abstract_proteinfiles


    def protein_fref(self):
        """ Get the header/identifier of the protein sequence """
        # if no (main) protein file is set -> return None
        if not self.protein_file: return None
        header, descr = parseSingleFastaHeaderFromFile(self.protein_file)
        return header

    # end of function protein_fref


    def get_protein_frefs(self):
        """ """
        return [ parseSingleFastaHeaderFromFile(fname)[0] for fname in self._files_protein() ]
    # end of function get_protein_frefs


    def set_protein_by_fref(self,fref):
        """ """
        for fname in self._files_protein():
            if fref == parseSingleFastaHeaderFromFile(fname)[0]:
                self.protein_file = fname
                return True
        else:
            # protein fref not found!
            return False

    # end of function set_protein_by_fref


    ############################################################################
    ### Functions for unigene accessibility
    ############################################################################

    def abstract_unigenefiles(self):
        """ """
        unigenefrefs = self.get_unigene_frefs()
        if not unigenefrefs:
            self.unigene_file = ""
        else:
            fref2annotation = []
            for fname in self._files_unigene():
                self.unigene_file = fname
                ident = self.are_gene_and_unigene_structure_identical()
                ugannotation = self.unigeneannotation()
                ugfref = self.unigene_fref()
                row = ( ident, int(ugannotation[1]), ugannotation[0], ugfref )
                fref2annotation.append(row)
            fref2annotation.sort()
            fref2annotation.reverse()
            # take longest unigene (first), potentially overruled by
            # an unigene which equals the annotated gene structure
            (ident,length,status,ugfref) = fref2annotation[0]
            self.set_unigene_by_fref(ugfref)

    # end of function abstract_unigenefiles


    def unigene_fref(self):
        """ Get the header/identifier of the unigene (-> TC00000, ... ) """
        # if no (main) unigene file is set -> return None
        if not self.unigene_file: return None
        return gnamefromgfftuple(parsegfffile(self.unigene_file)[0])

    # end of function unigene_fref


    def get_unigene_frefs(self):
        """ """
        return [ gnamefromgfftuple(parsegfffile(fname)[0]) for fname in self._files_unigene() ]
    # end of function get_unigene_frefs


    def set_unigene_by_fref(self,fref):
        """ """
        for fname in self._files_unigene():
            if fref == gnamefromgfftuple(parsegfffile(fname)[0]):
                self.unigene_file = fname
                return True
        else:
            # unigene fref not found!
            return False

    # end of function set_unigene_by_fref


    def _obtain_unigene_gff(self,input={}):
        """
        @type  input: dict
        @param input: SINGLE SUBDICT of input dictionary data structure
        """
        # handle potentially applied input argument
        self._handle_input_subdict(input)
        if self.input and self.input.has_key('gff-unigene') and self.input['gff-unigene']:
            return self.input['gff-unigene']
        elif self.unigene_file:
            # read gff file to list of tuples
            locusgff   = parsegfffile(self.locus_file)
            # define offset and length from locusgff
            offset = locusgff[0][3] - 1
            length = locusgff[0][4] - locusgff[0][3] + 1
            unigenegff = parsegfffile(self.unigene_file,offset=offset)
            if locusgff[0][6] == '-':
                unigenegff = mirrorgfflist( unigenegff, length=length)

            # correct unigene for utrs by the unigene annotation2DNA method
            unigenegff,typeofunigene = correct_unigene_for_utrs(
                unigenegff, dnaseqfname=self.dna_file
                )
            return unigenegff
        else:
            return []

    # end of function _obtain_unigene_gff


    ############################################################################
    ### Functions for locus & DNA file accessibility
    ############################################################################


    def abstract_locusfile(self):
        """ """
        self.locus_file = [os.path.join(os.path.abspath(self.dirname),x) for x in os.listdir(self.dirname) if x.endswith('.locus.gff') > 0][0]
    # end of function abstract_locusfile


    def abstract_dnafiles(self):
        """ """
        # TMP TODO HACK: there can be more than 1 dna file in a locus directory!
        self.dna_file   = [os.path.join(os.path.abspath(self.dirname),x) for x in os.listdir(self.dirname) if x.endswith('.fa') and x.find('.dna') > 0][0]
    # end of function abstract_dnafiles


    def dna_fref(self):
        """ Get the header/identifier of the DNA sequence file / e.g. the locus dir name """
        return os.path.basename( self.dirname )

    # end of function dna_fref


    def locus_fref(self):
        """ Get the header/identifier of the locus (-> chromosome, contig) """
        return open(self.locus_file).read().split('\t')[0]
    # end of function locus_fref


    def locus_fstrand(self):
        """ """
        return open(self.locus_file).read().split('\t')[6]
    # end of function locus_fstrand


    def locus_start(self):
        """ Get the start position of the gene locus """
        return int(open(self.locus_file).read().split('\t')[3])
    # end of function locus_start


    def locus_stop(self):
        """ Get the stop position of the gene locus """
        return int(open(self.locus_file).read().split('\t')[4])
    # end of function locus_stop



    def dnasequence(self):
        """
        Get *the* DNAsequence of this locus directory

        @rtype:  string
        @return: DNA sequence string

        @attention: TMP TODO HACK: there can be more than 1 DNA file !?
        """
        return parseSingleFasta(open(self.dna_file).readlines())[1]

    # end of function dnasequence


    def _splice_site_confirmation(self,introngfflist):
        """
        """
        warnings = []
        # get the DNA sequence belonging to this locus
        seq = self.dnasequence()
        for intron in [ getseqbyfeaturetuple(seq,gff) for gff in introngfflist ]:
            dsite, asite = intron[0:2].upper(), intron[-2:].upper()
            if dsite not in CANONICAL_DONOR_SITES:
                if dsite not in NON_CANONICAL_DONOR_SITES:
                    warnings.append( ImplausibleSpliceSiteWarning("Donor:"+dsite) )
                else:
                    warnings.append( NonCanonicalSpliceSiteWarning("Donor:"+dsite) )
            if asite not in CANONICAL_ACCEPTOR_SITES:
                warnings.append( ImplausibleSpliceSiteWarning("Acceptor:"+asite) )
        return warnings

    # end of function _splice_site_confirmation


    ############################################################################
    ### Functions for unigene DNA & protein sequence accessibility
    ############################################################################

    def get_unigene_dna_sequence(self,input={}):
        """
        """
        # handle potentially applied input argument
        self._handle_input_subdict(input)

        # if no (main) unigene file is set -> return None
        if not self.unigene_file: return None

        # read gff file to list of tuples
        locusgff   = parsegfffile(self.locus_file)
        # define offset and length from locusgff
        offset = locusgff[0][3] - 1
        length = locusgff[0][4] - locusgff[0][3] + 1
        unigenegff = parsegfffile(self.unigene_file,offset=offset)
        if locusgff[0][6] == '-':
            unigenegff = mirrorgfflist( unigenegff, length=length)

        # get unigene structure annotation; do not bother about UTRs
        unigeneexons = filtergffs4fmethod(unigenegff,GFF_UGEXON_FMETHOD)
        unigeneexons.sort()

        # apply protein_fref as fref for the unigene tracks
        unigeneexons = apply_fref(unigeneexons,self.dna_fref())

        # obtain unigene's full DNA sequence
        command = """ %s %s %s | grep -v ">" | tr -d "\n" """ % (
                PYTHON_PATH,
                EXECUTABLE_GFF2FASTA,
                self.dna_file,
                )
        ci,co,ce = osPopen3(command)  
        ci.write( gffs2txt( unigeneexons ) )
        ci.close()
        dnaseq = co.read().replace("\n","").strip()
        co.close()
        ce.close()

        # return the unigene's DNA sequence
        return dnaseq

    # end of function get_unigene_dna_sequence


    def get_unigene_cds_sequence(self,input={},correct_for_readingframe=True):
        """
        """
        # handle potentially applied input argument
        self._handle_input_subdict(input)

        # if no (main) unigene file is set -> return None
        if not self.unigene_file: return None

        # get tracks of unigene CODING exons
        ugcdsexons = filtergffs4fmethod( self._obtain_unigene_gff(), GFF_UGEXON_FMETHOD )
        ugcdsexons.sort()

        # apply protein_fref as fref for the unigene tracks
        ugcdsexons = apply_fref(ugcdsexons,self.dna_fref())

        # obtain unigene's full CDS sequence
        command = """ %s %s %s | grep -v ">" | tr -d "\n" """ % (
                PYTHON_PATH,
                EXECUTABLE_GFF2FASTA,
                self.dna_file,
                )
        ci,co,ce = osPopen3(command)            
        ci.write( gffs2txt( ugcdsexons ) )
        ci.close()
        cdsseq = co.read().replace("\n","").strip()
        co.close()
        ce.close()

        # return the unigene's CDS sequence
        return cdsseq

    # end of function get_unigene_cds_sequence


    def get_unigene_protein_sequence(self,input={}):
        """
        """
        # handle potentially applied input argument
        self._handle_input_subdict(input)

        # if no (main) unigene file is set -> return None
        if not self.unigene_file: return None

        # get unigene's CDS sequence
        cdsseq = self.get_unigene_cds_sequence()

        # get unigeneannotation tuple
        ugannotation = self.unigeneannotation()

        # define reading frame for translation
        if ugannotation[0] == "utr3p":  transeq_f_param = "F"
        elif ugannotation[0] == "fragm":transeq_f_param = "F"
        elif ugannotation[0] == "nonc": transeq_f_param = "1"
        else:                           transeq_f_param = "1"

        # obtain protein sequence from dna sequence
        command = """%s -filter -frame %s""" % (
                EXECUTABLE_TRANSEQ,
                transeq_f_param,
                )
        ci,co,ce = osPopen3(command)            
        ci.write( ">%s\n%s" % ("myheader", cdsseq) )
        ci.close()
        protseqs = parseFasta( co.readlines() )
        co.close()
        ce.close()
        if transeq_f_param == "1":
            protseq = protseqs.values()[0]
            frame = 0
        else:
            protseq = None
            frame = None
            for stop_cnt in [0,1]:
                for hdr,seq in protseqs.iteritems():
                    if seq.count("*") == stop_cnt:
                        if stop_cnt == 1 and not seq[-1] == '*':
                            # NO stop codon at the end -> incorrect frame!
                            continue
                        # if here, correct reading frame found
                        if frame == None:
                            protseq = seq
                            frame = int(hdr.split("_")[-1]) - 1
                        else:
                            # frame already assigned; multipel reading
                            # frames possible!! return None, None
                            # This can only happen for a unigene of
                            # type fragment (frag)
                            frame = None
                            protseq = None
                # if protein sequence assigned -> break out
                if protseq: break


        # return the protein sequence and its frame relative to the unigene
        return protseq, frame

    # end of function get_unigene_protein_sequence


    def has_unigene_implausible_splice_sites(self):
        """
        Is there a non-canonical splice site in the unigene alignment?
        """
        HAS_IMPLAUSIBLE_SPLICE_SITES = False
        warnings = []
        if self.unigene_file:
            # read gff file to list of tuples
            locusgff   = parsegfffile(self.locus_file)
            # define offset and length from locusgff
            offset = locusgff[0][3] - 1
            length = locusgff[0][4] - locusgff[0][3] + 1
            unigenegff = parsegfffile(self.unigene_file,offset=offset)
            if locusgff[0][6] == '-':
                unigenegff = mirrorgfflist( unigenegff, length=length)
            unigeneexons   = filtergffs4fmethod(unigenegff,GFF_UGEXON_FMETHOD)
            unigeneintrons = exons2introns(unigeneexons)
            warnings = self._splice_site_confirmation(unigeneintrons)
            if warnings: HAS_IMPLAUSIBLE_SPLICE_SITES=True
            return HAS_IMPLAUSIBLE_SPLICE_SITES, warnings
        else:
            return None, []

    # end of function has_unigene_implausible_splice_sites

 
    def gene_start_coordinate(self):
        """
        """
        try:
            return filtergffs4fmethod(
                    self._obtain_gene_gff(),
                    GFF_CDS_FMETHOD )[0][3]
        except IndexError:
            return False

    # end of function gene_start_coordinate


    def unigene_start_coordinate(self):
        """
        """
        try:
            return filtergffs4fmethod(
                    self._obtain_unigene_gff(),
                    GFF_UGEXON_FMETHOD )[0][3]
        except IndexError:
            return False

    # end of function unigene_start_coordinate


    def unigeneannotation(self):
        """
        """
        if self.unigene_file:
            # read gff file to list of tuples
            locusgff   = parsegfffile(self.locus_file)
            # define offset and length from locusgff
            offset = locusgff[0][3] - 1
            length = locusgff[0][4] - locusgff[0][3] + 1
            unigenegff = parsegfffile(self.unigene_file,offset=offset)
            if locusgff[0][6] == '-':
                unigenegff = mirrorgfflist( unigenegff, length=length)

            # get unigene structure annotation
            unigeneexons = filtergffs4fmethod(unigenegff,GFF_UGEXON_FMETHOD)
            unigeneexons.sort()

            # apply protein_fref as fref for the unigene tracks
            unigeneexons = apply_fref(unigeneexons,self.dna_fref())

            if min([gfftup[3] for gfftup in unigeneexons]) <= 0:
                # correct negative coordinates; can happen when
                # the unigene sticks out of the genelocus
                for pos in range(0,len(unigeneexons)):
                    if unigeneexons[pos][3] <= 0:
                        tmp = list(unigeneexons[pos])
                        tmp[3] = 1
                        unigeneexons[pos] = tuple(tmp)

            # run unigeneannotation command
            command = "%s %s %s" % (
                    PYTHON_PATH,
                    EXECUTABLE_UNIGENEANNOTATION,
                    self.dna_file )
            ci,co = os.popen2(command)
            ci.write( gffs2txt(unigeneexons) )
            ci.close()
            ugannotation = co.read().strip().split("\t")
            co.close()

            return ugannotation
        else:
            return None

    # end of function unigeneannotation


    def have_unigenes_alternative_splicing_evidence(self,verbose=False):
        """ """
        if not self._files_unigene(): return False
        if len(self._files_unigene()) == 1: return False

        has_alternative_splicing = False
        # more than 1 file. Read unigene coords into Sets
        current_unigene_file = deepcopy(self.unigene_file)

        thisunigenegff = self._obtain_unigene_gff()
        thisunigenecoordset = gffs2coordset(
                thisunigenegff, fmethod=[GFF_UGEXON_FMETHOD] )
        thisunigenearea = Set(range(min(thisunigenecoordset),max(thisunigenecoordset)+1))

        for fname in self._files_unigene():
            if fname == current_unigene_file: continue
            self.unigene_file = fname
            unigenegff = self._obtain_unigene_gff()
            unigenecoordset = gffs2coordset(
                unigenegff, fmethod=[GFF_UGEXON_FMETHOD] )
            unigenearea = Set(range(min(unigenecoordset),max(unigenecoordset)+1))

            if unigenecoordset.intersection(thisunigenearea).intersection(unigenearea).difference(thisunigenecoordset):
                has_alternative_splicing = True
                break
            if thisunigenecoordset.intersection(thisunigenearea).intersection(unigenearea).difference(unigenecoordset):
                has_alternative_splicing = True
                break

        # reset current original unigene
        self.unigene_file = current_unigene_file
        # return pointer for AS detection
        return has_alternative_splicing

    # end of function have_unigenes_alternative_splicing_evidence


    def are_gene_and_unigene_structure_identical(self,verbose=False):
        """
        """
        if self.unigene_file and self.gene_file:
            genegff    = self._obtain_gene_gff()
            unigenegff = self._obtain_unigene_gff()
            # make coord Set of gene coords
            gene_coordinate_set = gffs2coordset(
                genegff, fmethod=[GFF_CDS_FMETHOD] )

            # make coord Set of unigene coords
            unigene_coordinate_set = gffs2coordset(
                unigenegff, fmethod=[GFF_UGEXON_FMETHOD] )

            # calculate the difference
            difference = gene_coordinate_set.symmetric_difference(
                unigene_coordinate_set )

            if len(difference) == 0:
                # gene and unigene are identical!
                return True
            elif len(difference) == 3 and\
            ( max(difference) == max(gene_coordinate_set) or\
            max(difference) == max(unigene_coordinate_set) ):
                # gene and unigene are identical (stopcodon disambigouty)
                return True
            else:
                ################################################################
                if verbose:
                    print "dif:", len(difference), "gene:",
                    print len(gene_coordinate_set), "unigene:",
                    print len(unigene_coordinate_set),
                ################################################################
                if len(gene_coordinate_set)<=15:
                    ############################################################
                    if verbose: print "gene:", gene_coordinate_set
                    ############################################################
                elif len(unigene_coordinate_set)<=15:
                    ############################################################
                    if verbose: print "unigene:", unigene_coordinate_set
                    ############################################################
                else:
                    ############################################################
                    if verbose: print ""
                    ############################################################
                # return status False -> not identical
                return False
        else:
            return False

    # end of function are_gene_and_unigene_structure_identical 


    def _create_auto_key(self,identifier2organism={}):
        """
        Create a (logical) organism-based string-tag for this locus

        @rtype:  string
        @return: short string-tag that represents the organism of this locus
        """
        # when this is a locus in a dbwarehouse, abstract the genomedirname
        realdirname = os.path.realpath(self.dirname)
        if realdirname.find("/loci/") > 0:
            key = os.path.basename(realdirname[0:realdirname.find("/loci/")])
            if key: return key
        # if this point is reached, NOT a locus in dbwarehouse
        # check if we can map the gene's id to an organism ID
        if identifier2organism:
            for identifierpart,organism in identifier2organism.iteritems():
                if self.fref.find(identifierpart) == 0:
                    # succesfull mapping
                    return organism
            else:
                # mapping was not succesfull
                return self.fref
        else:
            return self.fref

    # end of function _create_auto_key


    def storetoggbdb(self,fref=None,getorf=True,gene=True,unigene=True,introns=True,splicesites=False,db_name=DB_NAME):
        """
        """
        # read gff file to list of tuples
        locusgff   = parsegfffile(self.locus_file)
        # define offset and length from locusgff
        offset = locusgff[0][3] - 1
        length = locusgff[0][4] - locusgff[0][3] + 1

        #list with known orf ids for visualization purposes
        known_orf_ids = []

        if unigene and self.unigene_file:
            #unigenegff = self._obtain_unigene_gff()
            #self._gffdata.extend( unigenegff )
            # abstract ALL unigenes from the locus directory i.s.o. only 1
            main_ugId = self.unigene_fref()
            for ugId in self.get_unigene_frefs():
                self.set_unigene_by_fref(ugId)
                self._gffdata.extend( self._obtain_unigene_gff() )
            # reset `main` unigene ID
            self.set_unigene_by_fref(main_ugId)


        if gene and self.gene_file:
            genegff = self._obtain_gene_gff()
            self._gffdata.extend( genegff )

        if introns:
            # we need abgp_geneconfirmation.geneconfirmation first!
            _retstruct = geneconfirmation( { self._create_auto_key(): self.input } )

            geneobjects = self.as_exons()
            known_introns = [ geneobjects[pos] for pos in range(1,len(geneobjects),2) ]
            known_exons   = [ geneobjects[pos] for pos in range(0,len(geneobjects),2) ]
            known_orf_ids = [ exon.orf.id for exon in known_exons ] 

            for intron in known_introns:
                intron.assign_bp_and_ppts()
                _gff = {'fref':fref,'fsource':'knownIntron','gclass':'KnownIntron'}
                knownintron_gfflines = intron.togff(gff=_gff,annotated_intron=True)
                # translate the actual intron track to 'IntronSpan'
                # Otherwise, it will pop up in the predicted IntronConnectingOrfs stack
                intron_gffline = list( knownintron_gfflines.pop(0) )
                intron_gffline[2] = 'IntronSpan'
                self._gffdata.append( tuple(intron_gffline) )
                self._gffdata.extend( knownintron_gfflines )

        if getorf:
            # perform getorf + tcode prediction(s)
            self.rungetorf()
            _retstruct = obtaintcodedata( { self._create_auto_key(): self.input } )

            # get gff lines for each ORF
            for orf in self.input['orfs'].orfs:
                if orf.id in known_orf_ids:
                    gfflist = list( orf.togff(fref=fref) )
                    gfflist[-1] = gfflist[-1] + "; knownOrf True"
                    self._gffdata.append( tuple(gfflist) )
                else:
                    self._gffdata.append( orf.togff(fref=fref) )

        # get TMHMM & SignalP locus data
        self._gffdata.extend( self.get_locus_tmhmm_gff() )
        self._gffdata.extend( self.get_locus_signalp_gff() )
        self._gffdata.extend( self.get_locus_proteomics_gff() )

        # scan the COMPLETE input sequence for high scoring TSS
        tsslist = scan_pssm_tss(self.input['genomeseq'],min_pssm_score=1.0)
        for tss in tsslist:
            _gff = {'fref':fref, 'fsource': 'tssPSSM4FUNGI', 'fmethod': 'tsspssm' }
            # format tss.pssm_score as a string-formatted float
            tss.pssm_score = "%2.1f" % tss.pssm_score
            self._gffdata.append( tss.togff( gff=_gff ) )

        if splicesites:
            # scan the COMPLETE input sequence for donor sites
            donorlist = scan_pssm_splice_site(self.input['genomeseq'],
                    splicetype="donor",min_pssm_score=1.0,
                    allow_non_canonical=True,
                    non_canonical_min_pssm_score=4.0)
            for dObj in donorlist:
                _gff = {'fref':fref, 'fsource': 'donorPSSMfungi', 'fmethod': 'donorpssm' }
                # format tss.pssm_score as a string-formatted float
                dObj.pssm_score = "%1.1f" % dObj.pssm_score
                self._gffdata.append( dObj.togff( gff=_gff ) )
            # scan the COMPLETE input sequence for acceptor sites
            acceptorlist = scan_pssm_splice_site(self.input['genomeseq'],
                    splicetype="acceptor",min_pssm_score=1.0,
                    allow_non_canonical=False)
            for aObj in acceptorlist:
                _gff = {'fref':fref, 'fsource': 'acceptorPSSMfungi', 'fmethod': 'acceptorpssm' }
                # format tss.pssm_score as a string-formatted float
                aObj.pssm_score = "%1.1f" % aObj.pssm_score
                self._gffdata.append( aObj.togff( gff=_gff ) )
                

        if self._gffdata:
            if not fref:
                fref = "locus_%s"  % os.path.basename(self.dirname)
            # update to this FREF
            for i in range(0,len(self._gffdata)):
                gffline = list(self._gffdata[i])
                gffline[0] = fref
                self._gffdata[i] = tuple(gffline)
            # create a Sequence track from the locus
            self._gffdata.insert( 0, ( fref, 'abfgp2.0', 'Sequence', 1, length, ".", "+", ".", "Sequence %s; GeneLocusDirectory %s" % (fref,fref) ) )
            # write (tmp) gff file
            tmp_gff_file = os.path.join(self.dirname,".tmp.storetoggbdb.gff")
            fh = open(tmp_gff_file, 'w')
            for line in self._gffdata:
                fh.write("\t".join([str(x) for x in line ])+"\n")
            fh.close()
            # write (tmp) fasta file -> adapt header to fref!
            tmp_fa_file = os.path.join(self.dirname,".tmp.storetoggbdb.fa")
            fh = open(tmp_fa_file, 'w')
            fh.write(">%s\n%s\n" % (fref,self.dnasequence()))
            fh.close()

            # cleanup the database 
            dbcleanup(fref=fref)
            # run storetoggbdb command
            gff2db( dbname=db_name,
                    dbhost=DB_HOST,
                    dbconnection=DB_CONNECTION,
                    dbuser=DB_USER,
                    dbpass=DB_PASS,
                    fnamegff=tmp_gff_file,
                    fnamefasta=tmp_fa_file)
            # and remove the tmp created gff file
            os.remove(tmp_gff_file)
            os.remove(tmp_fa_file)
            # return True because data is stored
            return True
        else:
            return False

    # end of function storetoggbdb


    def toinputdict(self,key=None):
        """ """
        if not key: key = self._create_auto_key()
        return { key: {
                'dirname':              self.dirname,
                'locusfile':            self.locus_file,
                'proteinseqfile':       self.protein_file,
                'genomeseqfile':        self.dna_file,
                'unigenefile':          self.unigene_file,
                'genefile':             self.gene_file,
                # backwards compatibilty to 'old' directory structure
                'gff':                  self.gene_file,
                'genomeseq':            '',
                'proteinseq':           '',
                'orfs':                 None,
                'tcode':                [],
                'blastdb':              {},
                'FREF':                 self.protein_fref(),
                'METHOD':               ABGP_VERSION,
                'genomefref':           self.locus_fref(), 
                'proteinfref':          self.protein_fref(),
                'locusfref':            self.locus_fref(),
                'fstrand':              self.locus_fstrand(),
                'gffs':                 [],
                'warnings':             [],
                'gldobj':               self,
                # keys for unigene informants
                'is_unigene':           False,
                'is_unigene_of':        None,
                'unigene-annotation':   None,
                } }
    # end of function toinputdict

    def _handle_input_subdict(self,input):
        """
        """
        if input:
            self.input = input
        elif not self.input:
            # create input data structure from this locusdir itself
            self.input = self.toinputdict().values()[0]
        else:
            pass
    # end of function _handle_input_subdict


    def rungetorf(self,input={}):
        """
        Run getorf and parse the output into Orf objects
    
        @type  input: dict
        @param input: SINGLE SUBDICT of input dictionary data structure
        """
        # handle potentially applied input argument
        self._handle_input_subdict(input)
        # check if orfs already predicted (in self.input data structure)
        if not self.input['orfs']:
            # predict ORFs on the dna sequence
            _retstruct = rungetorf( { self._create_auto_key(): self.input } )
        else:
            # ORFs already available in self.input['orfs']
            pass

    # end of function rungetorf


    def parseinputgff(self,input={}):
        """
        Parse the input gff file(s) into gene and the unigene tracks
    
        @type  input: dict
        @param input: SINGLE SUBDICT of input dictionary data structure
        """
        # handle potentially applied input argument
        self._handle_input_subdict(input)
        # check if gff is already parsed (in self.input data structure)
        if not self.input.has_key('gff-gene'):
            # parse gff file(s)
            _retstruct = parseinputgff({ self._create_auto_key(): self.input } )
        else:
            # gene gff etc data already parsed!
            pass

    # end of function parseinputgff


    def as_exons(self,input={}):
        """
        Get annotated gene structure as ExonOnOrf objects

        @type  input: dict
        @param input: SINGLE SUBDICT of input dictionary data structure

        @rtype:  list
        @return: list with ExonOnOrf objects
        """
        # handle potentially applied input argument
        self._handle_input_subdict(input)
        # parse data in the AbgpGeneLocusDir
        self.parseinputgff()
        self.rungetorf()
        # we need abgp_geneconfirmation.geneconfirmation first!
        geneconfirmation( { self._create_auto_key(): self.input } )

        # get only the CDS-type of tracks that define the coding sequence
        genecdstracks = filtergffs4fmethod( self._obtain_gene_gff(), GFF_CDS_FMETHOD ) 

        if len(genecdstracks) == 1:
            # deal with SingleExonOnOrf -> TSS + donor
            orf   = self.input['orfs'].get_orf_by_id(self.input['orfid-genestructure'][0])
            tss   = self._gene_cds_track_2_tss( genecdstracks[0], orf )
            return [ SingleExonOnOrf(tss,genecdstracks[-1][4],orf,gff={}) ]

        elif len(genecdstracks) == 0:
            # no tracks !?
            return []
        elif not self.input['orfid-genestructure']:
            # not mappable on Orfs / or no genestructure provided
            return []
        else:
            # list with exons,introns to return            
            exons = []
            introns = []
            exonsandintrons = []

            # deal with FirstExonOnOrf -> TSS + donor
            try:
                orf   = self.input['orfs'].get_orf_by_id(self.input['orfid-genestructure'][0])
            except:
                print self.input.keys(), self.input['proteinfref']
                orf   = self.input['orfs'].get_orf_by_id(self.input['orfid-genestructure'][0])

            tss   = self._gene_cds_track_2_tss( genecdstracks[0], orf )
            donor = self._gene_cds_track_2_donor( genecdstracks[0], orf )
            donor.phase = ( genecdstracks[0][4]-genecdstracks[0][3]-1 ) % 3
            exons.append( FirstExonOnOrf(tss,donor,orf,gff={}) )
            exonsandintrons.append( exons[-1] )

            # deal with internal ExonOnOrf(s): -> acceptor + donor
            for pos in range(1,len(genecdstracks)-1):
                orf   = self.input['orfs'].get_orf_by_id(self.input['orfid-genestructure'][pos])
                accep = self._gene_cds_track_2_acceptor( genecdstracks[pos], orf )
                accep.phase = exons[-1].donor.phase
                donor = self._gene_cds_track_2_donor( genecdstracks[pos], orf )
                donor.phase = ( genecdstracks[pos][4]-genecdstracks[pos][3]-1+accep.phase ) % 3
                exons.append( ExonOnOrf(accep,donor,orf,gff={}) )
                sharednts = get_shared_nucleotides_at_splicesite(
                    exons[-1].orf, exons[-2].orf,
                    exons[-1].acceptor, exons[-2].donor,
                    )
                intron = IntronConnectingOrfs(
                    exons[-2].donor, exons[-1].acceptor, sharednts,
                    exons[-2].orf, exons[-1].orf,
                    )
                introns.append(intron)
                exonsandintrons.append( introns[-1] )
                exonsandintrons.append( exons[-1] )

            # deal with FinalExonOnOrf -> acceptor + StopCodon
            orf   = self.input['orfs'].get_orf_by_id(self.input['orfid-genestructure'][-1])
            accep = self._gene_cds_track_2_acceptor( genecdstracks[-1], orf )
            accep.phase = exons[-1].donor.phase
            exons.append( FinalExonOnOrf(accep,genecdstracks[-1][4],orf,gff={}) )
            sharednts = get_shared_nucleotides_at_splicesite(
                exons[-1].orf,exons[-2].orf,
                exons[-1].acceptor,exons[-2].donor,
                )
            intron = IntronConnectingOrfs(
                exons[-2].donor, exons[-1].acceptor, sharednts,
                exons[-2].orf, exons[-1].orf,
                )
            introns.append(intron)
            exonsandintrons.append( introns[-1] )
            exonsandintrons.append( exons[-1] )

            # return list of exons&introns
            return exonsandintrons

    # end of function as_exons


    def _gene_cds_track_2_donor(self,track,orf,gff={}):
        """ """
        end   = track[4]
        offset= 1 # python/gff coordinate offset
        ccor  = -IC_DONOR_PATTERN_OFFSET[0]
        seqp  = orf.inputgenomicsequence[end+ccor:end+2+IC_DONOR_PATTERN_OFFSET[1]]
        score = _score_splice_site(seqp,splicetype='donor')
        donor = SpliceDonor(end+ccor,seqp,pssm_score=score,gff=gff)
        return donor
    # end of function _gene_cds_track_2_donor


    def _gene_cds_track_2_acceptor(self,track,orf,gff={}):
        """ """
        start = track[3]
        offset= 1 # python/gff coordinate offset
        ccor  = -IC_ACCEPTOR_PATTERN_OFFSET[0]-offset-2
        seqp  = orf.inputgenomicsequence[start+ccor:start-offset+IC_ACCEPTOR_PATTERN_OFFSET[1]]
        score = _score_splice_site(seqp,splicetype='acceptor')
        accep = SpliceAcceptor(start+ccor,seqp,pssm_score=score,gff=gff)
        return accep
    # end of function _gene_cds_track_2_acceptor


    def _gene_cds_track_2_tss(self,track,orf,gff={}):
        """ """
        start = track[3]
        offset= 1 # python/gff coordinate offset
        ccor  = -IC_TSS_PATTERN_OFFSET[0]-offset
        seqp  = orf.inputgenomicsequence[start+ccor:start-offset+3+IC_TSS_PATTERN_OFFSET[1]]
        score = score_tss(seqp)
        tss   = TranslationalStartSite(start+ccor,seqp,pssm_score=score,gff=gff)
        return tss
    # end of function _gene_cds_track_2_tss


    def printfiles(self):
        """ """
        print self
        print "LOCUS:", self.locus_file
        print "DNA  :", self.dna_file
        print "GENE :", self.gene_file
        print "PROT :", self.protein_file
        print "UNIGN:", self.unigene_file
    # end of function printfiles

    def filestats(self):
        """ """
        return "\n".join( [
            str(self),
            "LOCUS: %s" % self.locus_file,
            "DNA  : %s" % self.dna_file,
            "GENE : %s" % self.gene_file,
            "PROT : %s" % self.protein_file,
            "UNIGN: %s" % self.unigene_file,
            ] )
    # end of function filestats



    def get_locus_proteomics_gff(self):
        """ """

        """ """
        # TMP TODO HACK: there can be more than 1 dna file in a locus directory!
        self.dna_file   = [os.path.join(os.path.abspath(self.dirname),x) for x in os.listdir(self.dirname) if x.endswith('.fa') and x.find('.dna') > 0][0]

        proteomics_files = [os.path.join(os.path.abspath(self.dirname),x) for x in os.listdir(self.dirname) if x.endswith('.gff') and x.find('.proteomics.') > 0]
        # no gff for proteomics data
        if not proteomics_files: return []

        # read locus gff file to obtain offset/strand
        locusgff   = parsegfffile(self.locus_file)
        # define offset and length from locusgff
        offset = locusgff[0][3] - 1
        length = locusgff[0][4] - locusgff[0][3] + 1

        allproteomicsgffdata = []
        for proteomicsfile in proteomics_files:
            proteomicsgff = parsegfffile(proteomicsfile,offset=offset)
            if locusgff[0][6] == '-':
                proteomicsgff = mirrorgfflist( proteomicsgff, length=length)
            allproteomicsgffdata.extend( proteomicsgff )
        # return all parsed data
        return allproteomicsgffdata
    
    # end of function get_locus_proteomics_gff


    ############################################################################
    ### Functions for TSS objects accessibility
    ############################################################################

    def get_locus_tss_objects(self,input={},gff={}):
        """ """
        if not self._locus_objects_tss:
            # handle potentially applied input argument
            self._handle_input_subdict(input)

            # handle gff defaults for this object type
            gff_defaults = {'fref'   : self.fref,
                            'fsource': GFF_TSS_FSOURCE,
                            'fmethod': GFF_TSS_FMETHOD }
            gff_defaults.update(gff)

            for tssobj in scan_pssm_tss(self.dnasequence(),min_pssm_score=1.0):
                tssobj._gff.update(gff_defaults)
                self._locus_objects_tss.append( tssobj )

        # return list of objects
        return self._locus_objects_tss

    # end of function get_locus_tss_objects


    def get_locus_tss_gff(self,input={},gff={}):
        """ """
        return [ tss.togff( gff=gff ) for tss in self.get_locus_tss_objects(input=input) ]

    # end of function get_locus_tss_gff


    ############################################################################
    ### Functions for SignalP objects accessibility
    ############################################################################

    def get_locus_signalp_objects(self,input={}):
        """ """
        # handle potentially applied input argument
        self._handle_input_subdict(input)

        # return list of objects
        return self._locus_objects_signalp

    # end of function get_locus_signalp_objects


    def get_locus_signalp_gff(self,input={}):
        """ """
        # handle potentially applied input argument
        self._handle_input_subdict(input)

        try:
            # we need orf prediction to match SIGNALP on Orfs
            self.rungetorf()
            gfflines = []
            orfids2orf = {}
            theorf   = None
            theorfid = None
            fname = "%s.locus.signalp.out" % self.fref
            for line in open(os.path.join(self.dirname,fname)).readlines():
                parts = line.strip().split(" ")
                while '' in parts: parts.remove('')
                if SIGNALP_MIN_D and float(parts[12]) >= SIGNALP_MIN_D:
                    fscore  = float(parts[12])
                    splen   = int(parts[5])
                elif SIGNALP_MIN_SPROB and float(parts[-2]) >= SIGNALP_MIN_D:
                    fscore  = float(parts[-2])
                    splen   = int(parts[-4])
                else:
                    continue
                # fsource & method
                fsource = "SignalP3.0"
                fmethod = "predSignalPeptide"
                # obtain coordinates
                header_parts = parts[13].split("_")
                atg_nt_start = int(header_parts[-2])
                orf_nt_end   = int(header_parts[-1])

                # calculate start coords relative to locus
                signalp_start  = atg_nt_start
                signalp_end    = atg_nt_start + splen*3
                gffline = ( self.fref, fsource, fmethod,
                            signalp_start, signalp_end,
                            "%1.3f" % fscore, "+", ".", 
                            )
                gfflines.append( gffline )
            else:
                pass

            # return parsed file
            return gfflines

        except IOError:
            # signalp prediction not performed
            return []

    # end of function get_locus_signalp_gff


    def load_protein_signalp_to_orf(self,input={}):
        """ """
        # handle potentially applied input argument
        self._handle_input_subdict(input)

        try:
            fname = "%s.signalp.out" % self.fref
            line = open(os.path.join(self.dirname,fname)).readlines()[0]
            parts = line.strip().split(" ")
            while '' in parts: parts.remove('')
            # if no signalpeptide predicted -> return None
            if parts[-1] == 'N': return None
        except IOError:
            # signalp prediction not performed
            return False
        except IndexError:
            # file exists but is empty
            return False

        exonsandintrons = self.as_exons()
        exons   = [ exonsandintrons[pos] for pos in range(0,len(exonsandintrons),2) ]
        introns = [ exonsandintrons[pos] for pos in range(1,len(exonsandintrons),2) ]

        # no recognizable intron-exon gene structure for currently loaded gene
        if not exons: return False

        # calculate AA length & score of SignalPeptide
        if SIGNALP_MIN_D and float(parts[12]) >= SIGNALP_MIN_D:
            signalp_score  = float(parts[12])
            signalp_length = int(parts[5])*3
        elif SIGNALP_MIN_SPROB and float(parts[-2]) >= SIGNALP_MIN_D:
            signalp_score  = float(parts[-2])
            signalp_length = int(parts[-4])*3
        else:
            return False # SignalPeptide not of high enough score

        # calcuate coordinated and start creating object
        signalp_start  = exons[0].start
        signalp_end    = exons[0].start + signalp_length
        signalp_tss    = deepcopy(exons[0].acceptor) # acceptor == tss for first exon!
        signalp_tss._gff['fsource'] = 'tssPSSM4FUNGI'
        signalp_tss._gff['fmethod'] = 'tsspssmSignalPeptide'

        if signalp_length > exons[0].length:
            # SignalPeptide branched by an intron!
            sgp = SignalPSignalPeptide(signalp_start,signalp_end,
                    signalp_score,tss=signalp_tss,intron=introns[0])
            # load SignalPeptide to Orf object -> IT CANNOT EXIST
            # in the locus SignalP data because of the intron
            exons[0].orf._signalp_sites.append(sgp)
            exons[0].orf._has_signalp_sites_predicted = True

        else:
            # create SignalPSignalPeptide object
            sgp = SignalPSignalPeptide(signalp_start,signalp_end,signalp_score,
                    tss=signalp_tss)
            # load SignalPeptide to Orf object IF NOT EXISTS ALREADY!
            if (sgp.start,sgp.end) not in\
            [ (_sgp.start,_sgp.end) for _sgp in exons[0].orf._signalp_sites ]:
                exons[0].orf._signalp_sites.append(sgp)
                exons[0].orf._has_signalp_sites_predicted = True

        # return status == True for succesfull loading of SignalP object
        return True


    # end of function load_protein_signalp_to_orf


    def load_locus_signalp_to_orfs(self,input={}):
        """ """
        # handle potentially applied input argument
        self._handle_input_subdict(input)
        try:
            fname = "%s.locus.signalp.out" % self.fref
            lines = open(os.path.join(self.dirname,fname)).readlines()
        except IOError:
            # signalp prediction not performed
            return 0

        # we need predicted TSS sites
        tssobjects = self.get_locus_tss_objects()

        # we need orf prediction to match SIGNALP on Orfs
        self.rungetorf()
        signalpcnt = 0
        gfflines = []
        orfids2orf = {}
        theorf   = None
        theorfid = None
        for line in lines:
            parts = line.strip().split(" ")
            while '' in parts: parts.remove('')
            score  = float(parts[-2])
            header_parts = parts[13].split("_")
            atg_nt_start = int(header_parts[-2])
            orf_nt_end   = int(header_parts[-1])

            # get Orf object and do TSS prediction
            orfObj = self._get_orf_object_by_coords(None,orf_nt_end)
            if not orfObj:
                # this will result in a hard crash! Log message to debug it
                print "Warning:: no Orf object in %s with end position %s" % (self.fref,orf_nt_end)
                print line
                for orfObj in self.input['orfs'].orfs:
                    if atg_nt_start > orfObj.start and atg_nt_start < orfObj.end:
                        print orfObj
                orfObj = None

            orfObj.scan_orf_for_pssm_tss(min_pssm_score=TSS_MIN_PSSM_SCORE)

            # calculate PYTHON (not GGB!) start coords relative to locus
            signalp_start  = atg_nt_start -1
            signalp_end    = atg_nt_start -1 + int(parts[-4])*3
            # find TSS object belonging to the SignalPSignalPeptide
            tssObj = None
            for tsso in tssobjects:
                if tsso.pos == signalp_start:
                    tssObj = tsso
                    tssObj._gff['fmethod'] = "tsspssmSignalPeptide"
                    break
            else:
                # tss not found -> Methionine but NOT a likely tss
                tss_start = signalp_start - IC_TSS_PATTERN_OFFSET[0]
                tss_seqp  = "N"* ( sum(IC_TSS_PATTERN_OFFSET) + 3 )
                tssObj = TranslationalStartSite(tss_start,tss_seqp,pssm_score=0.0)
                tssObj._gff['fmethod'] = "tsspssmSignalPeptide"
                pass

            # create SignalPSignalPeptide object and store to Orf
            sgp = SignalPSignalPeptide(signalp_start,signalp_end,score,tss=tssObj)
            orfObj._signalp_sites.append(sgp)
            signalpcnt+=1
        else:
            pass

        # do not forget to label all Orf objects as _has_signalp_sites_predicted
        for orfobj in self.input['orfs'].orfs:
            orfobj._has_signalp_sites_predicted = True

        # return counter of how many SignalPeptides are loaded
        return signalpcnt

    # end of function load_locus_signalp_to_orfs


    ############################################################################
    ### Functions for TMHMM objects accessibility
    ############################################################################


    def get_locus_tmhmm_objects(self,input={}):
        """ """
        # handle potentially applied input argument
        self._handle_input_subdict(input)

        # return list of objects
        return self._locus_objects_tmhmm

    # end of function get_locus_tmhmm_objects


    def get_locus_tmhmm_gff(self,input={}):
        """ """
        # handle potentially applied input argument
        self._handle_input_subdict(input)

        try:
            # we need orf prediction to match TMHMMs on Orfs
            self.rungetorf()
            gfflines = []
            orfids2orf = {}
            theorf   = None
            theorfid = None
            fname = "%s.locus.tmhmm.out" % self.fref
            for line in open(os.path.join(self.dirname,fname)).readlines():
                while line.find("  ") >= 0:
                    line = line.replace("  "," ")
                while line.find(" \t") >= 0:
                    line = line.replace(" \t","\t")
                while line.find("\t ") >= 0:
                    line = line.replace("\t ","\t")
                line = line.replace(" ","\t")
                parts = line.strip().split("\t")
                fsource = parts[1]
                fmethod = parts[2]
                header_parts = parts[0].split("_")
                orf_nt_start = int(header_parts[-2])
                orf_nt_end   = int(header_parts[-1])
                # get Orf objects on which the TMHMM is predicted
                orf_key      = (orf_nt_start,orf_nt_end)
                if not orfids2orf.has_key(orf_key):
                    if theorf:
                        # TMHMMs on next ORFs -> make Orf tailing stop track
                        gfflines.append( (
                                self.fref, fsource, "OrgTailingStop",
                                theorf.end+1,theorf.end+3,
                                ".", "+", ".", "Orf %s" % theorfid
                                ) )
                    # find Orf by coords
                    for orfobj in self.input['orfs'].orfs:
                        if orfobj.start == orf_nt_start and orfobj.end == orf_nt_end:
                            theorf = orfobj
                            theorfid = orfobj.id
                            orfids2orf[orf_key] = theorf
                            # calulate `projected` leading stop gff lines
                            if theorf.start >= 4:
                                gfflines.append( (
                                        self.fref, fsource, "OrgLeadingStop",
                                        theorf.start-3,theorf.start-1,
                                        ".", "+", ".", "Orf %s" % theorfid
                                        ) )
                            # orf found -> break out
                            break

                    else:
                        theorf = None
                        theorfid = None
                else:
                    theorf = orfids2orf[orf_key]
                    theorfid = theorf.id

                # calculate start coords relative to locus
                tmhmm_start  = orf_nt_start + int(parts[3])*3
                tmhmm_end    = orf_nt_start + int(parts[4])*3
                gffline = ( self.fref, fsource, fmethod,
                            tmhmm_start, tmhmm_end,
                            ".", "+", ".", "Orf %s" % theorfid
                            )
                gfflines.append( gffline )
            else:
                # EOF lines -> store tailing stop of final TMHMM Orf
                if theorf:
                    # TMHMMs on next ORFs -> make Orf tailing stop track
                    gfflines.append( (
                            self.fref, fsource, "OrgTailingStop",
                            theorf.end+1,theorf.end+3,
                            ".", "+", ".", "Orf %s" % theorfid
                            ) )

            # return parsed file
            return gfflines

        except IOError:
            # tmhmm prediction not performed
            return []

    # end of function get_locus_tmhmm_gff


    def _get_orf_object_by_coords(self,start,end):
        """
        Get Orf object from input['orfs'] by nt coords

        @attention: only use when self.rungetorf() is performed!
        """
        for orfobj in self.input['orfs'].orfs:
            if orfobj.start == start and orfobj.end == end:
                return orfobj
            elif start == None and orfobj.end == end:
                # used in case of SignalP output -> no Orf start coord given!
                return orfobj
            else:
                pass
        else:
            return None

    # end of function _get_orf_object_by_coords


# end of class AbgpGeneLocusDirectory



########################################################################
#### Exceptions
########################################################################
class InproperlyAppliedArgument(Exception):
    """ """
    def __init__(self, message):
        self.message = message
    def __str__(self):
        return repr(self.message)

class NoAbgpGeneLocusDirectory(InproperlyAppliedArgument):
    pass


class UnequalGeneLocusLengths(InproperlyAppliedArgument):
    pass


########################################################################
#### Validators
########################################################################
def IsAbgpGeneLocusDirectory(dirname,verbose=False):
    """
    """
    validation_status = True
    if not os.path.isdir(dirname):
        if verbose: print "'%s' is not a recognized directory in: %s" % (dirname,os.getcwd())
        validation_status = False
    else:
        # check if a *.locus.gff AND a .dna*.fa file are present
        locus_files = [x for x in os.listdir(dirname) if x.endswith('.locus.gff')]
        dna_files   = [x for x in os.listdir(dirname) if x.endswith('.fa') and x.find('.dna')]
        if not locus_files:
            if verbose: print "no locus_file(s) '*.locus.gff' in directory '%s'" % dirname
            validation_status = False
        if len(locus_files) > 1:
            if verbose: print "more than 1 locus_file(s) '*.locus.gff' in directory '%s'" % dirname
            validation_status = False
        if not dna_files:
            if verbose: print "no dna_file(s) '*.dna*.fa' in directory '%s'" % dirname
            validation_status = False
    return validation_status

# end of validator function IsAbgpGeneLocusDirectory


def make_abgpgenelocusdirectory_from_fasta(hdr,sequence,outdir):
    """ Helper function which makes a (near-empty) AbgpGeneLocusDirectory from a provided hdr and sequence

    @type   hdr: string
    @param  hdr: fasta accession header

    @type   sequence: string
    @param  sequence: DNA sequence string of this informant

    @type   outdir: string
    @param  outdir: !existing! directory on the filesystem in which to make the AbgpGeneLocusDirectory
    
    @rtype:  AbgpGeneLocusDirectory
    @return: instantiated AbgpGeneLocusDirectory
    """
    safe_header   = hdr.translate(NONLETTERS_TRANS).replace("_","")
    unique_string = get_random_string_tag(10)
    if not safe_header:
        safe_header = unique_string
    locus_dirname = os.path.join(outdir,"%s_%s" % (safe_header,unique_string))
    if not os.path.exists(locus_dirname):
        os.mkdir(locus_dirname)
    fname_dna = os.path.join(locus_dirname,"%s.%s.dna.fa" % (safe_header,unique_string) )
    fname_gff = os.path.join(locus_dirname,"%s.%s.locus.gff" % (safe_header,unique_string) )
    writeSingleFasta(hdr,sequence,fname_dna)
    fh = open(fname_gff,'w')
    fh.write("%s\t%s\t%s\t1\t%s\t.\t+\t.\n" % (safe_header,ABGP_VERSION,"abgpGLD",len(sequence)))
    fh.close()
    gldObj = AbgpGeneLocusDirectory(locus_dirname)
    return gldObj

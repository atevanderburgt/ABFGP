#!/usr/bin/python
""" Quick test if all python imports performed from the various ABFGP packages,modules and scripts will work """
__authors__ = "Ate van der Burgt"
__license__ = "MIT"

# BioPython NCBIStandalone blast parser
import Bio
if Bio.__version__ != "1.60":
    print "Warning!! Biopython version != 1.60, but '%s'. Bio.Blast.NCBIStandalone will likely cause errors." % Bio.__version__
from Bio.Blast import NCBIStandalone

# Generic Python Imports
import gzip
import os
import os.path
import re
import sys
from sys import path as sysPath

#print sys.path

# function imports from os,re,sys,copy,time,random
from os import popen as osPopen
from os import popen2 as osPopen2
from os import popen3 as osPopen3
from os import popen2
from os import popen3
from os import remove as osRemove
from os import system
from os import system as osSystem
from os import getcwd as osGetcwd
from os import mkdir as osMkdir
from os import chdir as osChdir

from os.path import exists as osPathExists
from os.path import join as osPathJoin
from os.path import isfile as osPathIsfile
from os.path import dirname as osPathDirname
from os.path import abspath as osPathAbspath
from os.path import split as osPathSplit


from re import compile
from re import finditer
from re import IGNORECASE
from sys import exit as sysExit
from copy import deepcopy
from time import time
from random import randint
from random import sample
from math import log as mathlog
from math import sqrt

# Python imports that might be absent depending on the version of Python
from sets import Set # deprecated but backwards compatible in newer verions of python
from optparse import OptionParser
from optparse import OptionGroup
from optparse import OptionValueError
from StringIO import StringIO

# Most likely Numeric is not installed; this is solved by having it installed
# within the ABFGP code tree
# The path to Numeric is obtained either from the settings dir (absolute)
# or the path is reconstructed (relative)

# absolute path acquisition
from settings.abgp import MAIN_ABGP_PATH as BASEPATH
sys.path.append(osPathJoin(BASEPATH,"requiredmodules"))
from Numeric import zeros, where, greater, greater_equal

# relative path acquisition
sys.path.append(osPathJoin(osPathDirname(osPathAbspath(__file__)),"requiredmodules"))
from Numeric import zeros, where, greater, greater_equal

from abgp_etc import abgpsysexit
from abgp_etc import _blastdb_cleanup
from abgp_etc import _blastdb_cleanup, _file_cleanup
from abgp_exceptions import InproperlyAppliedArgument
from abgp_exceptions import NoCrossdataApplied, NoInputApplied
from abgp_geneconfirmation import *
from abgp_geneconfirmation import geneconfirmation
from abgpgenelocusdirectory import AbgpGeneLocusDirectory
from abgpgenelocusdirectory import IsAbgpGeneLocusDirectory
from abgpgenelocusdirectory import make_abgpgenelocusdirectory_from_fasta
from abgp_logging import logf
from abgp_logging import Logger
from abgp_unigeneconfirmation import geneandunigeneconfirmation
from abgp_warnings import *
from abgp_warnings import RepetitiveProteinSequenceWarning
from abgp_warnings import UniGeneStructureIsNotMappableOnOrfsWarning, GeneStructureIsNotMappableOnOrfsWarning

import dna2prot
from dna2prot import dna2protein
from dna2prot import dna2proteinbyframe

from gene.acceptor import SpliceAcceptor, ProjectedSpliceAcceptor
from gene.codingblock import CodingBlockStart, CodingBlockEnd
from gene.donor import SpliceDonor, ProjectedSpliceDonor
from gene.exon import ExonOnOrf
from gene.exon import FirstExonOnOrf
from gene.exon import FinalExonOnOrf
from gene.exon import SingleExonOnOrf
from gene.intron import IntronConnectingOrfs
from gene.orf import TcodeCodingOrf
from gene.pssm import pssmscore, parse_ic_file
from gene.sequenceerror import merge_orfs_with_sequenceerror
from gene.signalpeptide import SignalPSignalPeptide
from gene.splicesite import scan_pssm_splice_site
from gene.splicesite import _score_splice_site
from gene.start import scan_pssm_tss, TranslationalStartSite, IC_TSS_PATTERN_OFFSET
from gene.start import score_tss
from gene.start import TranslationalStartSite
from gene.stop import StopCodon
from gene.validators import IsProperStrand

from gff.exceptions import *
from gff import GffFeature
from gff import gffs2coordset
from gff import filtergffs4fmethod
from gff import gffs2txt
from gff import parsegfffile

import graphAbgp
import graphPlus
from graphAbgp.exceptions import OrganismNotPresentInGraph
from graphAbgp.graph_pacbpcollection import _delete_pacbp
from graphAbgp.graph_pacbpcollection import pacbporf2PCG
from graphAbgp.graph_codingblock import CodingBlockGraph
from graphAbgp import conversion
from graphAbgp import GeneTreeGraph
from graphAbgp import InwardsPointingCodingBlockGraph
from graphAbgp.ordering import order_list_by_attribute
from graphAbgp.subclass_sitealignment import sort_by_cumulative_score
from graphPlus.comparison import mutual_nodes
#from graphAbgp.exceptions import NoOverallMinimalSpanningRange


from lib_abgpgff import exons2introns
from lib_abgpgff import get_raw_abfgp_genestruture
from lib_blastp import formatdb, blastall_seq2db
from lib_cexpander import runcexpander
from lib_cexpander import cexpander_checkCBG4omsrbordergaps
from lib_cexpander import ZeroUniformlyAlignedPositions
from lib_clustalw import clustalw
from lib_clustalw import clustalw, strip_alignment_for_exterior_gaps
from lib_crossblast import iterativecrossblastp, createblastdbs
from lib_crossblast import _order_list_by_attribute
from lib_crossdatafunctions import *
from lib_crossdatafunctions import createcrossdata
from lib_crossdatafunctions import create_pacbpcollectiongraph_from_crossdata
from lib_fasta import parseFasta
from lib_fasta import parseSingleFasta
from lib_fasta import parseSingleFastaHeaderFromFile
from lib_fasta import IsSingleFastaDna
from lib_fasta import IsSingleFastaProtein
from lib_fasta import parseDecoratedFasta
from lib_fasta import writeMultiFasta
from lib_fasta import writeSingleFasta
from lib_filevalidation import fastafilesequencelength
from lib_filevalidation import gfffilesequencelength
from lib_genestructure import make_consensus_genestructure_from_compatible_pacb_graphs
from lib_hmm_pairwise import is_hmmpacbporf_conflicting_with_pacbporflist
from lib_intron import *
from lib_intron import _order_intron_list
from lib_orfset import OrfSet
from lib_orfset import GetorfOrfSet

from lib_pcg2blocks import *
from lib_pcg2blocks import PCG2inwpCBGS
from lib_proteinsimilaritymatrix import make_clustalw_alignment_match
from lib_proteinsimilaritymatrix import ProteinSimilarityMatrix
from lib_rip import RIPindex
from lib_sequenceperiodicity import coding_sequence_periodicity_test
from lib_sequenceperiodicity import importgeneticcode
from lib_sequenceperiodicity import test
from lib_sequencerepetitiveness import annotatedproteinsequencerepetitivenesscheck
from lib_sequencerepetitiveness import proteinsequencelowcomplexityscore
from lib_signalp import create_projected_signalp_for_informants
from lib_signalp import ProjectedSignalPSignalPeptide
from lib_stopwatch import StopWatch
from lib_synteny_pairwise import get_first_and_final_inwpcbg_pos
from lib_tcode import obtaintcodedata
from lib_tinyexononorf import scan_orf_for_tiny_exon
from lib_tmhmm import obtainlocustmhmmdata

from listofcodingblockgraphs import ListOfCodingBlockGraphs

import pacb
import pacb.conversion
import pacb.recombination
from pacb import PacbPCOORDS
from pacb import PacbP, PacbPORF
from pacb import swap_query_and_sbjct
from pacb.pacbp import PacbP
from pacb.connecting.functions import _update_kwargs
from pacb.connecting import project_splicesites_on_pacbporfs_with_lacking_intron_in_sbjct as PrjIntLsbj
from pacb.connecting import project_splicesites_on_pacbporfs_with_lacking_intron_in_query as PrjIntLqry
from pacb.connecting.mapping import _filter_aligned_introns_on_pssm_entropy_combination
from pacb.connecting.mapping import merge_pacbporfs_with_introns
from pacb.connecting.orfs import get_potential_tiny_exons_on_orf
from pacb.connecting.orfs import get_potention_first_exons_on_orf
from pacb.connecting.orfs import get_signalp_first_exon_on_orf
from pacb.connecting.projecting import merge_pacbporfs_by_intron_in_query
from pacb import conversion as pacbconversion
from pacb.conversion import pacbp2pacbporf
from pacb.conversion import pacbp_from_clustalw
from pacb.conversion import pacbporf2pacbp
from pacb.exceptions import CoordinateOutOfRange
from pacb.ordering import order_pacbp_list
from pacb.ordering import order_pacbporf_list
from pacb.overlap import correct_overlap_for_sbjct
from pacb.overlap import correct_overlap_for_query
from pacb.splitting import split_pacb_on_gaps

from parsers.cexpander import *
from parsers.clustalw import *
from parsers.fasta import *
from parsers.fasta import _parseFastaHeader
from parsers.dna import reversecomplement as _reversecomplement
from parsers.getorf import getorf
from parsers.hmm import hmmbuild_protein
from parsers.hmm import hmmsearch_protein
from parsers.proteinsimilaritymatrix import *
from parsers.proteinsimilaritymatrix import _generate_alignment_match
from parsers.tcode import *

from pythonlibs.optparsefunctions import *
from pythonlibs.ordering import order_list_by_attribute
from pythonlibs.sequencerepetitiveness import *
from pythonlibs.uniqueness import get_random_string_tag
from pythonlibs.stdoutmanaging import stdoutManagerClass

from settings import *
from settings.abgp import *
from settings.dbwarehouse import *
from settings.genestructure import *
from settings.genetreegraph import *
from settings.gff import *
from settings.ggb import *
from settings.gff.cbginterface import *

from settings.abgp import ABGP_VERSION
from settings.abgp import MAIN_ABGP_PATH as BASEPATH
from settings.abgp import USE_TCODE
from settings.blastp import blastoptions_iter2 as BLASTOPTIONS
from settings.dbwarehouse import ABGP_LOCUS_IDENTIFIER2ORGANISM_MAPPING
from settings.dbwarehouse import _get_organism_full_name
from settings.emboss import BLOSUM62_PATH
from settings.executables import EXECUTABLE_GETORF_VERSION
from settings.executables import EXECUTABLE_UNIGENEANNOTATION, PYTHON_PATH
from settings.executables import TCODE_MAX_NONCODING
from settings.executables import TCODE_MIN_CODING
from settings.genestructure import MAX_INTRON_NT_LENGTH
from settings.genestructure import MIN_INTRON_NT_LENGTH
from settings.genestructure import MIN_EXON_NT_LENGTH
from settings.genestructure import OPTIMAL_BRACNHPOINT_TO_ACCEPTOR_DISTANCE
from settings.genestructure import ORF_IS_UNIGENE_LABEL
from settings.gff.currentannotation import GFF_CDS_FMETHOD
from settings.gff import GFF_GENESTART_FMETHOD
from settings.gff import GFF_GENESTOP_FMETHOD
from settings.gff import GFF_CDS_FMETHOD
from settings.gff import GFF_UGEXON_FMETHOD
from settings.gff import GFF_UG3UTREXON_FMETHOD
from settings.gff import GFF_UG5UTREXON_FMETHOD
from settings.inframeintron import INFRAME_INTRON_MIN_AA_LENGTH
from settings.pacbp import LINEARIZATION_STARTWITHBEST
from settings.pacbp import LINEARIZATION_ACCEPTANCE_MARGIN
from settings.pacbp import LINEARIZATION_WEIGTHED
from settings.splicesites import IC_ACCEPTOR_DATA_FILE
from settings.splicesites import IC_ACCEPTOR_PATTERN_OFFSET
from settings.splicesites import IC_DONOR_DATA_FILE
from settings.splicesites import IC_DONOR_PATTERN_OFFSET
from settings.splicesites import IC_DONOR_NCGC_DATA_FILE
from settings.splicesites import KWARGS_TINYEXON_PAIRWISE
from settings.translationalstartsites import IC_TSS_DATA_FILE
from settings.translationalstartsites import IC_TSS_PATTERN_OFFSET
from settings.translationalstartsites import TSS_MIN_PSSM_SCORE
from settings.translationalstartsites import TSS_ALLOW_NON_CANONICAL
from settings.translationalstartsites import TSS_NON_CANONICAL_MIN_PSSM_SCORE

# if here, then none of the imports raised Exceptions
print "All imports succesfully performed!"

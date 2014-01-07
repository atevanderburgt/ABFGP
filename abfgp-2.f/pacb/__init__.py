"""
Pairwise-Aligned Coding Block package, used in ABFGP

This package contains 4 main classes

class PacbP
    Pairwise-Aligned Coding Block PROTEIN
    a blastp.hsp - like object of 2 aligned protein sequences
class PacbPDNA(PacbP)
    Pairwise-Aligned Coding Block PROTEIN & DNA
    a PacbP extended with its corresponding DNA sequeces
    call it a Protein-templated DNA-alignment
class PacbPORFS(PacbPDNA)
    Pairwise-Aligned Coding Block PROTEIN & ORFs
    a PacbPDNA objected with ORF objects of both query and sbjct
    a PacbP extended with its corresponding DNA sequeces,
    and the possibility to extend the alignment bejoined the
    currently aligned part (originating from a blastp-hsp)
    call it a Protein-templated DNA-alignment, with optional extension
    for Orf objects, see gene.orf classes
class AlignedProteinDnaPosition
    An aligned position: positions on sbjct & query, 2 AAs and 2x3 NTs
    each PacbPDNA and PacbPORF object is constructed from these

def swap_query_and_sbjct
    swaps query and sbjct in a PacbP, PacbPDNA or PacbPORF object

"""

# Python Imports
from copy import deepcopy

# Imports from pacb package
from alignedproteindnaposition import AlignedProteinDnaPosition
from alignedproteindnaposition import IsIdenticalAlignedProteinDnaPosition
from pacbp import PacbP
from pacbpdna import PacbPDNA
from pacbporf import PacbPORF
from pacbpcoords import PacbPCOORDS

# Import Functionality packages for PacbP(DNA/ORF) objects
import linearization
import ordering
import connecting
import splitting
import conversion
import merging
import overlap
import tcode
import validators           # validators accessable via pacb.validators
from validators import *    # validators directly accessable
from exceptions import * 

# ProteinSimilarityMatrix import
from lib_proteinsimilaritymatrix import ProteinSimilarityMatrix
DEFAULT_MATRIX = ProteinSimilarityMatrix(name="BLOSUM62")

##############################################################################
### some helper functions ####################################################
##############################################################################

def swap_query_and_sbjct(pacb):
    """
    This function can be used on PacbP and PacbPDNA classes
    """
    new_pacb = deepcopy(pacb)
    new_pacb._swap_query_and_sbjct()
    return new_pacb

# end of function swap_query_and_sbjct



def calculate_bitscoreratio(query,sbjct,matrix=DEFAULT_MATRIX):
    """
    Calculate the bitscore ratio between a given query and sbjct

    @type  query: string
    @param query: query AA string 

    @type  sbjct: string
    @param sbjct: sbjct AA string

    @type  matrix: ProteinSimilarityMatrix 
    @param matrix: ProteinSimilarityMatrix object (default BLOSUM62)

    @rtype:  float
    @return: bitscore ratio ( >>0.0 .. 1.0 )
    """
    bitQuery = float(matrix.scorealignment(query,query))
    bitSbjct = float(matrix.scorealignment(sbjct,sbjct))
    bitAlign = float(matrix.scorealignment(query,sbjct))
    try:
        return bitAlign / max([bitQuery,bitSbjct])
    except ZeroDivisionError:
        # query and sbjct bitscore == 0!
        return 1.0

# end of function calculate_bitscoreratio


def calculate_identityscore(alignment,verbose=False):
    """
    Calculate the identityscore% based on the alignment string

    @attention: based on an alignment string, NOT on a blast match string

    @type  alignment: string
    @param alignment: alignment string to calculate identityscore from

    @rtype:  float
    @return: identityscore ( >>0.0 .. 1.0 )
    """
    length      = float(len(alignment))
    if length == 0.0:
        if verbose:
            # zero-length alignments can occur!
            print "SeriousWarning, ZeroDivisionError in pacb.calculate_identityscore(alignment) ; zero-length alignment"
        # return identity of 0.0
        return 0.0
    else:
        identity    = float(alignment.count("*"))
        similarity  = length - float(alignment.count(" ")) - identity
        return ( identity + 0.5*similarity ) / length

# end of function calculate_identityscore


def calculate_identity(alignment):
    """
    Calculate the identity% based on the alignment string

    @attention: based on an alignment string, NOT on a blast match string

    @type  alignment: string
    @param alignment: alignment string to calculate identityscore from

    @rtype:  float
    @return: identity ( >>0.0 .. 1.0 )
    """
    length      = float(len(alignment))
    if length == 0.0:
        print "SeriousWarning, ZeroDivisionError in pacb.calculate_identity(alignment) ; zero-length alignment"
        return 0.0
    else:
        identity = float(alignment.count("*"))
        return identity / length

# end of function calculate_identity


"""
AlignedProteinDnaPosition class used in PacbPDNA and PacbPORF objects

class AlignedProteinDnaPosition
    An aligned position: positions on sbjct & query, 2 AAs and 2x3 NTs
    each PacbPDNA and PacbPORF object is constructed from these

function IsIdenticalAlignedProteinDnaPosition
    Are two AlignedProteinDnaPosition indentical?
 
"""

# Imports from pacb package
from exceptions import *

class AlignedProteinDnaPosition:
    """
    Object representing an aligned Amino Acid and its accompanying DNA sequence.
    So, the DNA alignment is guided by the Protein alignment!
    """
    def __init__(self,position,query_pos,sbjct_pos,alignment=("","",""),
        query_dna=("",0),sbjct_dna=("",0)):
        """
        Initialization function of AlignedProteinDnaPosition

        @type  position: zero or positive integer
        @param position: relative position in the protein alignment (zero based coordinates)

        @type  query_pos: positive integer
        @param query_pos: current amino acid position in the query sequence (zero based coordinates)

        @type  sbjct_pos: positive integer
        @param sbjct_pos: current amino acid position in the query sequence (zero based coordinates)

        @type  alignment: tuple of 3 elements (query,match,sbjct)
        @param alignment: tuple of 3 elements (
                            query   Single uppercase Amino Acid letter of query sequence
                            match   alignment match status, options are:
                                    identity    Single uppercase Amino Acid letter
                                    similarity  '+'
                                    gap         ' ' or '' (space) or empty string
                                    other       ' ' or '' (space) or empty string
                            sbjct   Single uppercase Amino Acid letter of sbjct sequence
                            )

        @type  query_dna: tuple of 2 elements (dna_triplet_sequence, abs_dna_pos)
        @param query_dna: tuple of 2 elements representing the query DNA sequence of currence query Amino Acid (
                            dna_triplet_sequence    DNA triplet sequence; ("",0) denotes a gap; not given denotes a gap
                            abs_dna_pos             absolute DNA position in the query sequence (zero based coordinates)
                            )

        @type  sbjct_dna: tuple of 2 elements (dna_triplet_sequence, abs_dna_pos)
        @param sbjct_dna: tuple of 2 elements representing the sbjct DNA sequence of current sbjct Amino Acid (
                            dna_triplet_sequence    DNA triplet sequence; ("",0) denotes a gap; not given denotes a gap
                            abs_dna_pos             absolute DNA position in the sbjct sequence (zero based coordinates)
                            )

        """
        query,match,sbjct= alignment
        self.position    = position
        self.query       = query
        self.match       = match
        self.sbjct       = sbjct
        self.query_pos   = query_pos
        self.sbjct_pos   = sbjct_pos
        self.isa_gap     = False
        self.isa         = ""       # identity, similarity, gap or None
        if self.match == '+':
            self.isa = 'similarity'
        elif self.query == self.sbjct:
            self.isa = 'identity' 
        else:
            self.isa = 'None'
        if query_dna[0] or query_dna[1]:
            self.query_dna_seq     = query_dna[0]
            if not self.query_dna_seq:
                self.query_dna_seq = "   "
                self.isa_gap = True
                self.isa     = 'gap'
            self.query_dna_start = query_dna[1]
            self.query_dna_end   = query_dna[1]+3
        else:
            self.query_dna_seq  = "   "
            self.isa_gap        = True
            self.isa     = 'gap'
            self.query_dna_start= 0
            self.query_dna_end  = 0
        if sbjct_dna[0] or sbjct_dna[1]:
            self.sbjct_dna_seq     = sbjct_dna[0]
            if not self.sbjct_dna_seq:
                self.sbjct_dna_seq = "   "
                self.isa_gap = True
                self.isa     = 'gap'
            self.sbjct_dna_start = sbjct_dna[1]
            self.sbjct_dna_end   = sbjct_dna[1]+3
        else:
            self.sbjct_dna_seq  = "   "
            self.isa_gap        = True
            self.isa     = 'gap'
            self.sbjct_dna_start= 0
            self.sbjct_dna_end  = 0

    # end of function __init__


    def __str__(self):
        """ """
        return "<%s %s P(%s,%s) D(%s,%s)>" % (
                self.__class__.__name__,
                self.position,
                self.query_pos,
                self.sbjct_pos,
                self.query_dna_start,
                self.sbjct_dna_start
                )

    # end of function __str__

    def _swap_query_and_sbjct(self):
        """
        USE WITH CARE!!!
        see PacbPDNA._swap_query_and_sbjct()
        """
        # swap actual query/sbjct characters and position
        self.query, self.sbjct = self.sbjct, self.query
        self.query_pos, self.sbjct_pos = self.sbjct_pos, self.query_pos
        # swap query/sbjct DNAsequence and coordinates
        (a,b,c) = (self.query_dna_seq, self.query_dna_start, self.query_dna_end)
        (d,e,f) = (self.sbjct_dna_seq, self.sbjct_dna_start, self.sbjct_dna_end)
        (self.query_dna_seq, self.query_dna_start, self.query_dna_end) = (d,e,f)
        (self.sbjct_dna_seq, self.sbjct_dna_start, self.sbjct_dna_end) = (a,b,c)

    # end of function _swap_query_and_sbjct


    def print_protein(self):
        """ """
        print "\n".join( self._return_print_protein() )

    def print_protein_and_dna(self):
        """ """
        print "\n".join( self._return_print_protein_and_dna() )

    def print_dna(self):
        """ """
        print "\n".join( self._return_print_dna() )

    def _return_print_protein(self):
        """ """
        return [ self.query, self.match, self.sbjct ]

    def _return_print_dna(self):
        """ """
        return [ self.query_dna_seq, self.sbjct_dna_seq ]

    def _return_print_protein_and_dna(self):
        """ """
        (_q, _m, _s)   = self._return_print_protein()
        (_qdna, _sdna) = self._return_print_dna()
        return [ _qdna, " %s " % _q, " %s " % _m, " %s " % _s, _sdna ]

# end of class AlignedProteinDnaPosition


def IsIdenticalAlignedProteinDnaPosition(objA,objB,quickanddirty=True):
    """
    Are 2 AlignedProteinDnaPosition identical? 

    @type  objA: AlignedProteinDnaPosition
    @param objA: AlignedProteinDnaPosition object

    @type  objB: AlignedProteinDnaPosition
    @param objB: AlignedProteinDnaPosition object

    @type  quickanddirty: Boolean
    @param quickanddirty: if True, check only AAs, not DNAs

    @rtype:  Boolean
    @return: Boolean
    """
    if objA == CoordinateOutOfRange or objB == CoordinateOutOfRange:
        return False
    elif objA.query_pos == objB.query_pos and objA.sbjct_pos == objB.sbjct_pos:
        # coordinates are the same; check if aa's are identical too
        if objA.query == objB.query and objA.sbjct == objB.sbjct:
            # objects are (presumably) identical.
            # If not param quickanddirty, check the exact DNA bases and coordinates too
            if quickanddirty:
                return True
            else:
                if objA.query_dna_start == objB.query_dna_start and objA.sbjct_dna_start == objB.sbjct_dna_start:
                    if objA.query_dna_seq == objB.query_dna_seq and objA.sbjct_dna_seq == objB.sbjct_dna_seq:
                        return True
                    else:
                        return False
                else:
                    return False
        else:
            return False
    else:
        return False

# end of function IsIdenticalAlignedProteinDnaPosition

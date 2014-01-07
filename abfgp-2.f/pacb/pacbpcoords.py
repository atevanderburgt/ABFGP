"""
"""

# Module metadata
__authors__ = "Ate van der Burgt"
__license__ = "MIT"


# Imports from pacb package
from pacbp import PacbP
from exceptions import * 
from validators import *

# Python Imports
from copy import deepcopy
from sets import Set


class PacbPCOORDS(PacbP):
    def __init__(self,input=("","",0,0),blastp_hsp=None,MATRIX=None,gff={}):
        """ """
        PacbP.__init__(self,input=input,blastp_hsp=blastp_hsp,MATRIX=MATRIX,gff=gff)

    # end of function __init__


    def __str__(self):
        """ """
        return """<%s Q:%s-%s S:%s-%s len:(%s,%s)>""" % (
            self.__class__.__name__,
            self.query_start, self.query_end,
            self.sbjct_start, self.sbjct_end,
            self.query_length, self.sbjct_length )

    # end of function __str__


    def _initialize(self,input):
        """
        Overwrites PacbP._initialize(input)
        """
        if self._hsp:
            # initialize from a blastp HSP
            self.query       = deepcopy(self._hsp.query)
            self.match       = deepcopy(self._hsp.match)
            self.sbjct       = deepcopy(self._hsp.sbjct)
            self.query_start = self._hsp.query_start-1
            self.sbjct_start = self._hsp.sbjct_start-1
            self.length      = len(self.query)
            if self.length == 0:
                raise ZeroSizedPacb
            self.source = 'blastp'
        else:
            # validate input data!
            self._validate_input(input)
            # initialize from scratch
            self.query       = input[0].upper()
            self.sbjct       = input[1].upper()
            self.query_start = input[2]
            self.sbjct_start = input[3]
            self.length      = max([ len(self.query),len(self.sbjct) ])

        # and set the now common attributes
        self.query_protein  = self.query.replace("-","")
        self.sbjct_protein  = self.sbjct.replace("-","")
        self.query_end      = self.query_start + self.length -\
                              self.query.count("-")
        self.query_length   = self.query_end - self.query_start
        self.sbjct_end      = self.sbjct_start + self.length -\
                              self.sbjct.count("-")
        self.sbjct_length   = self.sbjct_end - self.sbjct_start

    # end of function _initialize


    def _validate_input(self,input):
        """
        if anything False, an error is raised
        if all oke, nothing is raised or returned
        """
        if len(input) != 4:
            message = "`input` is not a 4-element tuple"
            raise InproperlyAppliedArgument, message
        if len(input[0]) == 0:
            # this MUST be raised to present ZeroDivision Exception
            message = "no input query or sbjct given"
            raise InproperlyAppliedArgument, message
        if type(input[2]) != type(int()) or input[2] < 0:
            raise CoordinateOutOfRange
        if type(input[3]) != type(int()) or input[3] < 0:
            raise CoordinateOutOfRange

    # end of function _validate_input


    def _make_alignment_match(self,q,s):
        pass

    def print_protein(self,_linesize=120):
        pass

    def set_matrix(self,matrix):
        pass

    def get_matrix_name(self):
        return None

    def alignment_has_gaps(self,gap_size=1):
        return None

    def query_has_gaps(self,gap_size=1,abs_coords=(0,0)):
        return None

    def sbjct_has_gaps(self,gap_size=1,abs_coords=(0,0)):
        return None

    def _score_alignment(self):
        pass

    def bitscore_slice_by_abs_protein_query(self,qstart,qstop):
        return 0

    def bitscore_slice_by_abs_protein_sbjct(self,sstart,sstop):
        return 0

    def identityscore_slice_by_abs_protein_query(self,qstart,qstop):
        return 0.0

    def identityscore_slice_by_abs_protein_sbjct(self,sstart,sstop):
        return 0.0

    def _match2alignment(self):
        pass

    def togff(self,gff={}):
        pass

    def strip_unmatched_ends(self):
        pass

    def _strip_one_left_position(self):
        pass

    def _strip_one_rigth_position(self):
        pass

    def slice(self, abs_start, abs_end, coords_on='query'):
        """ Not applicable on PacbPCOORD object """
        return None

    def returnslice(self, abs_start, abs_end, coords_on='query'):
        """ Not applicable on PacbPCOORD object """
        return None

    def potentially_contains_aligned_intron(self,**kwargs):
        return False

    def alignmentposition_by_query_pos(self,qposPY,forced_return=False):
        """ Not applicable on PacbPCOORD object """
        return None

    def alignmentposition_by_sbjct_pos(self,sposPY,forced_return=False):
        """ Not applicable on PacbPCOORD object """
        return None

    def get_aligned_protein_sequences(self):
        """ Not (really) applicable on PacbPCOORD object """
        return (self.query,self.match,self.sbjct)

    def get_unextended_aligned_protein_sequences(self):
        """ Not (really) applicable on PacbPCOORD object """
        return (self.query,self.match,self.sbjct)

# end of class PacbPCOORDS
"""
"""

# Module metadata
__authors__ = "Ate van der Burgt"
__license__ = "MIT"

# Imports from pacb package
import ordering
from exceptions import * 
from validators import *

# Python Imports
from copy import deepcopy
from sets import Set
from re import finditer

# Abgp Imports
from lib_proteinsimilaritymatrix import ProteinSimilarityMatrix

DEFAULT_MATRIX = ProteinSimilarityMatrix(name="BLOSUM62")


class PacbP:
    # Pairwise-Aligned Coding Block PROTEIN; something like an exon
    # alternative name: ProteinAlignment
    # this `flavour` only has an aligned Protein sequence
    def __init__(self,input=("","",0,0),blastp_hsp=None,MATRIX=DEFAULT_MATRIX,gff={}):
        """
        initialize an object of this class by either a:
        blastp_hsp      blastp HSP parsed by Bio.Blast.NCBIStandalone parser
        input           tuple (query,sbjct,query_start,sbjct_start)
                        query       must be an ALIGNED Protein sequence
                        sbjct       must be an ALIGNED Protein sequence
                        query_start in 0-based Protein coordinates
                        sbjct_start in 0-based Protein coordinates
        MATRIX          a ProteinSimilarityMatrix object
        """
        # fill input attribute data
        self._hsp           = blastp_hsp
        self._original_pacb = self      # TODO: is sucha construction okay??
        self.MATRIX         = MATRIX
        # set string attributes for the alignment
        self.query          = ""
        self.match          = ""
        self.sbjct          = ""
        self.alignment      = ""
        self.identity       = 0
        self.similarity     = 0
        self.identityscore  = 0
        self.bitscore       = 0
        self.bits           = self.bitscore
        # Source of where this PacbP(DNA/ORF) originates from
        # Feel free to introduce your own terminology; there is no
        # check upon the source of a PacbP(ORF) anywhere in the code.
        # The sources mentioned here are used in ABGP:
        #     unknown             unknown (default)
        #     blastp              alignment produced with blastp
        #     hmmsearch           alignment produced with a hmmprofile
        #     clustalw            alignment produced with (pairwise) clustalw
        #     clustalw-EXTENDED   alignment that is re-aligned with clustalw based on assumption of elongatable
        #     clustalw-OPTIMIZED  existing alignment that is OPTIMIZED with clustalw multiple alignment
        #     lsrPACBP            coordinate range that bridges the gap between 2 CBGs (used in ABGP)
        self.source         = 'unknown'
        # set integer attributes for the positions
        self.query_start    = 0
        self.query_end      = 0
        self.query_length   = 0
        self.sbjct_start    = 0
        self.sbjct_end      = 0
        self.sbjct_length   = 0
        self.query_protein  = ""
        self.sbjct_protein  = ""
        self.length         = 0
        # novel attributes defining if pacbp can be deleted/edited
        self._IS_DELETE_PROTECTED = False
        self._IS_EDIT_PROTECTED = False
        # set/update gff data
        self._gff = {
                'fref'      : None,
                'fsource'   : 'getorf-TBLASTX',
                'gclass'    : 'Similarity',
                'fphase'    : '.',
                'fstrand'   : '+',
                'column9data' : {},
            }
        self._gff.update(gff)
        # and finish initialization
        self._initialize(input)

    # end of function __init__

    def __str__(self):
        """
        """
        source = ' '
        if self.source != 'blastp':
            source = " src=%s " % self.source
        return """<pacbp%sQ:%s-%s S:%s-%s bits:%s iden:%1.3f len:(%s,%s)>""" % (
            source,
            self.query_start, self.query_end,
            self.sbjct_start, self.sbjct_end,
            self.bitscore, self.identityscore,
            self.query_end - self.query_start,
            self.sbjct_end - self.sbjct_start )

    # end of function __str__

    def __len__(self):
        """
        """
        return self.length

    # end of function __len__

    
    def barcode(self):
        """
        A unique tuple with integers describing this object
        """
        return ( self.query_start, self.query_end, self.sbjct_start, self.sbjct_end, self.bitscore )

    # end of function barcode


    def _initialize(self,input):
        """
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
            self.match       = self._make_alignment_match(
                                    self.query,
                                    self.sbjct )
            self.query_start = input[2]
            self.sbjct_start = input[3]
            self.length      = len(self.query)

        # and set the now common attributes
        self._score_alignment()
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
        if len(input[0]) != len(input[1]):
            message = "query and sbjct of input not equal in length"
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
        """
        Create a blastp alignment match string between
        a query and a sbjct sequence
        """
        match = []
        for pos in range(0,len(q)):
            score = self.MATRIX.scoreaapair( q[pos], s[pos] ) 
            if q[pos] == s[pos] and score > 0:
                # identical amino acids
                match.append(q[pos])
            elif q[pos] == s[pos] and score <= 0:
                # identical char, but no amino acid (stop-codon?)
                match.append(" ")
            else:
                if score > 0:
                    match.append("+")
                else:
                    match.append(" ")
        # and return
        return "".join(match)

    # end of function _make_alignment_match


    ########################################################################
    #### Functions for printing PacbP objects                           ####
    ########################################################################

    def print_protein(self,_linesize=120):
        """
        Print protein alignment of this PacbP object (query,match,sbjct)
        """
        for offset in range(0,len(self),_linesize):
            if offset > 0: print ""
            print self.query[offset:offset+_linesize]
            print self.match[offset:offset+_linesize]
            print self.sbjct[offset:offset+_linesize]

    # end of function print_protein


    def set_matrix(self,matrix):
        """ """
        self.MATRIX = matrix
        self._score_alignment()

    # end of function set_matrix


    def get_matrix_name(self):
        """ """
        return self.MATRIX.name

    # end of function get_matrix_name


    def construct_unique_key(self,node1,node2):
        """
        (Re)create the (semi!!)-unique key of this pacbp-type object

        @type  node1: tuple( gene_id, orf_id )
        @param node1: Node-identifier

        @type  node2: tuple( gene_id, orf_id )
        @param node2: Node-identifier

        @rtype:  tuple
        @return: Tuple of this form: (bits, length, orf-id-Q, orf-id-S)

        @attention: this key is semi-unique; in theory 2 pacbps with identical length and bitscore can be created!
        """
        # do *NOT* use self.length -> extended PacbPORF.length != PacbP.length
        return ( self.bits, self.get_unextended_length(), node1[1], node2[1] )

    # end of function construct_unique_key


    def alignment_has_gaps(self,gap_size=1):
        """
        Are there gaps (-) in the alignment (query or sbjct)?

        @type  gap_size: Integer (positive)
        @param gap_size: Length of minimal continious gap length

        @rtype:  Boolean
        @return: True or False
        """
        for strand in [self.query,self.sbjct]:
            if strand.find('-'*gap_size) >= 0:
                return True
        else:
            return False

    # end of function alignment_has_gaps


    def query_has_gaps(self,gap_size=1,abs_coords=(0,0)):
        """
        Are there gaps (-) in the aligned query sequence?

        @type  gap_size: Integer (positive)
        @param gap_size: Length of minimal continious gap length

		@type  abs_coords:  tuple of two positive integers
		@param abs_coords:  absolute AA start and stop position coordinates (on query)

        @rtype:  Boolean
        @return: True or False
        """
        if abs_coords != (0,0):
            sta = abs_coords[0] - self.query_start
            end = abs_coords[1] - self.query_start
            gap_offset = 0
            # compensate for gaps in the frontal part of the alignment
            while sta+gap_offset - self.query[0:sta+gap_offset].count('-') < sta:
                gap_offset+=1
            # and compensate for leading gaps
            while self.query[sta+gap_offset] == '-': gap_offset+=1
            # correct coordinates for LEFT gap_offset
            sta = sta+gap_offset
            end = end+gap_offset
            gap_offset = 0
            while end-sta+gap_offset - self.query[sta:end+gap_offset].count('-') < end-sta:
                gap_offset+=1
            # correct end coordinate for RIGHT gap_offset
            end = end+gap_offset
            # and now check if there is a gap
            if self.query[sta:end].find('-'*gap_size) >= 0:
                return True
            else:
                return False
        else:
            if self.query.find('-'*gap_size) >= 0:
                return True
            else:
                return False

    # end of function query_has_gaps


    def sbjct_has_gaps(self,gap_size=1,abs_coords=(0,0)):
        """
        Are there gaps (-) in the aligned sbjct sequence?

        @type  gap_size: Integer (positive)
        @param gap_size: Length of minimal continious gap length

		@type  abs_coords:  tuple of two positive integers
		@param abs_coords:  absolute AA start and stop position coordinates (on sbjct)

        @rtype:  Boolean
        @return: True or False
        """
        if abs_coords != (0,0):
            sta = abs_coords[0] - self.sbjct_start
            end = abs_coords[1] - self.sbjct_start
            gap_offset = 0
            # compensate for gaps in the frontal part of the alignment
            while sta+gap_offset - self.sbjct[0:sta+gap_offset].count('-') < sta:
                gap_offset+=1
            # and compensate for leading gaps
            while self.sbjct[sta+gap_offset] == '-': gap_offset+=1
            # correct coordinates for LEFT gap_offset
            sta = sta+gap_offset
            end = end+gap_offset
            gap_offset = 0
            while end-sta+gap_offset - self.sbjct[sta:end+gap_offset].count('-') < end-sta:
                gap_offset+=1
            # correct end coordinate for RIGHT gap_offset
            end = end+gap_offset
            # and now check if there is a gap
            if self.sbjct[sta:end].find('-'*gap_size) >= 0:
                return True
            else:
                return False
        else:
            if self.sbjct.find('-'*gap_size) >= 0:
                return True
            else:
                return False

    # end of function sbjct_has_gaps


    def gap_count(self,gap_size=1):
        """
        Count number of gaps of at least gap_size in query AND sbjct

        @type  gap_size: Integer (positive)
        @param gap_size: Length of minimal continious gap length

		@type  count:  Integer
		@param count:  Number of distinct gaps of at least gap_size in alignment
        """
        (q,m,s) = self.get_aligned_protein_sequences()
        cntQ = len([ match for match in finditer("-{%s,}" % int(gap_size),q) ])
        cntS = len([ match for match in finditer("-{%s,}" % int(gap_size),s) ])
        return cntQ+cntS

    # end of function gap_count


    def query_gap_count(self,gap_size=1):
        """
        Count number of gaps of at least gap_size in query

        @type  gap_size: Integer (positive)
        @param gap_size: Length of minimal continious gap length

		@type  count:  Integer
		@param count:  Number of distinct gaps of at least gap_size in query
        """
        (q,m,s) = self.get_aligned_protein_sequences()
        return len([ match for match in finditer("-{%s,}" % int(gap_size),q) ])

    # end of function query_gap_count


    def sbjct_gap_count(self,gap_size=1):
        """
        Count number of gaps of at least gap_size in query

        @type  gap_size: Integer (positive)
        @param gap_size: Length of minimal continious gap length

		@type  count:  Integer
		@param count:  Number of distinct gaps of at least gap_size in sbjct
        """
        (q,m,s) = self.get_aligned_protein_sequences()
        return len([ match for match in finditer("-{%s,}" % int(gap_size),s) ])

    # end of function sbjct_gap_count


    def gap_ratio_score(self,gap_size=1):
        """
        Score the relative presence of distinct gaps in alignment

        @type  gap_size: Integer (positive)
        @param gap_size: Length of minimal continious gap length

		@type  ratio:  Float
		@param ratio:  ratio of gap_count / max([query,sbjct]) length
        """
        # allow a single gap to be present in the alignment
        cnt = max([0, self.gap_count(gap_size=gap_size) - 1 ])
        query_len = len(self.alignment_protein_range_query())
        sbjct_len = len(self.alignment_protein_range_sbjct())
        # `punish` many gaps by taking 2th power of corrected gap count
        return float(cnt*cnt) /float(max([query_len,sbjct_len]))

    # end of function gap_ratio


    def query_gap_ratio_score(self,gap_size=1):
        """
        Score the relative presence of distinct gaps in aligned query sequence

        @type  gap_size: Integer (positive)
        @param gap_size: Length of minimal continious gap length

		@type  ratio:  Float
		@param ratio:  ratio of gap_count / query length
        """
        # allow a single gap to be present in the alignment
        cnt = max([0, self.query_gap_count(gap_size=gap_size) -1 ])
        # `punish` many gaps by taking 2th power of corrected gap count
        return float(cnt*cnt) / float(len(self.alignment_protein_range_query()))

    # end of function query_gap_ratio


    def sbjct_gap_ratio_score(self,gap_size=1):
        """
        Score the relative presence of distinct gaps in aligned sbjct sequence

        @type  gap_size: Integer (positive)
        @param gap_size: Length of minimal continious gap length

		@type  ratio:  Float
		@param ratio:  ratio of gap_count / sbjct length
        """
        # allow a single gap to be present in the alignment
        cnt = max([0, self.sbjct_gap_count(gap_size=gap_size) -1 ])
        # `punish` many gaps by taking 2th power of corrected gap count
        return float(cnt*cnt) / float(len(self.alignment_protein_range_sbjct()))

    # end of function sbjct_gap_ratio


    def _score_alignment(self):
        """
        Score several properties of the alignment.
        After a manual change of a Pacb object (splitting or extension),
        do not forget to re-call this function!
        This (re-calling) is taken care for in the functions:
            ``to do``
            ``to do``
        """
        # (re) set the clustalW-like alignment string
        self._match2alignment()
        if self.length == 0:
            self.identity       = 0
            self.similarity     = 0
            self.identityscore  = 0.0
            self.bitscore       = 0
            self.bits           = self.bitscore
        else:
            # do some counts on this string
            self.identity       = self.alignment.count("*")
            self.similarity     = len(self.alignment) - self.alignment.count(" ") -\
                              self.identity
            self.identityscore  = float(self.identity) / float(self.length) +\
                              (0.5*float(self.similarity)) / float(self.length) 
            self.bitscore       = self.MATRIX.scorealignment(self.query,self.sbjct)
            self.bits           = self.bitscore

    # end of function _score_alignment


    def bitscore_slice_by_abs_protein_query(self,qstart,qstop):
        """
        """
        rel_qstart = qstart - self.query_start 
        rel_qstop  = qstop - qstart + rel_qstart
        return self.MATRIX.scorealignment(
            self.query[rel_qstart:rel_qstop],
            self.sbjct[rel_qstart:rel_qstop]
            )

    # end of function bitscore_alignment_part_by_query


    def bitscore_slice_by_abs_protein_sbjct(self,sstart,sstop):
        """
        """
        rel_sstart = sstart - self.sbjct_start 
        rel_sstop  = sstop - sstart + rel_sstart
        return self.MATRIX.scorealignment(
            self.query[rel_sstart:rel_sstop],
            self.sbjct[rel_sstart:rel_sstop]
            )

    # end of function bitscore_alignment_part_by_query


    def identityscore_slice_by_abs_protein_query(self,qstart,qstop):
        """
        """
        rel_qstart = qstart - self.query_start 
        rel_qstop  = qstop - qstart + rel_qstart
        algmslice  = self.alignment[rel_qstart:rel_qstop]
        return calculate_identityscore(algmslice)

    # end of function identityscore_slice_by_abs_protein_query


    def identityscore_slice_by_abs_protein_sbjct(self,sstart,sstop):
        """
        """
        rel_sstart = sstart - self.sbjct_start 
        rel_sstop  = sstop - sstart + rel_sstart
        algmslice  = self.alignment[rel_sstart:rel_sstop]
        return calculate_identityscore(algmslice)

    # end of function identityscore_slice_by_abs_protein_sbjct


    def _match2alignment(self):
        """
        Translate blastp ``match`` string into clustalW like format
         - all AAs (identities)      into `*`
         - all +   (similarities)    into `+` (no translation)
         - all *   (aligned *,X,?)   into ` `
        """
        alignment = []
        for symbol in self.match:
            if symbol in ["*","X"]:
                alignment.append(" ")
            elif symbol in ["+"," "]:
                alignment.append(symbol)
            else:
                alignment.append("*")
        # and set to PACB object
        self.alignment = "".join(alignment)

    # end of function _match2alignment


    def _swap_query_and_sbjct(self):
        """
        USE WITH CARE!!!
        Swap query and sbjct sequences and positions.
        This enables easier comparison in some cases.
        It is done WITHIN THIS INITIALIZED object.
        That is what I mean with use with care.
        Recommended usage is only with the helper function
            ``swap_query_and_sbjct``
        Usage:
            ``swapped_pacb = swap_query_and_sbjct(pacb)``
        """
        # swap strings
        self.query, self.sbjct = self.sbjct, self.query
        self.query_protein, self.sbjct_protein =\
            self.sbjct_protein, self.query_protein

        # swap coordinates
        (a,b,c) = (self.query_start, self.query_end, self.query_length)
        (d,e,f) = (self.sbjct_start, self.sbjct_end, self.sbjct_length)
        (self.query_start, self.query_end, self.query_length) = (d,e,f)
        (self.sbjct_start, self.sbjct_end, self.sbjct_length) = (a,b,c)

        # swap the original blastp_hsp input object too!
        if self._hsp: self._swap_input_hsp_qands()

    # end of function _swap_query_and_sbjct


    def _swap_input_hsp_qands(self):
        """
        USE WITH CARE!!!
        See _swap_query_and_sbjct
        """
        # swap the original blastp_hsp input object
        self._hsp = deepcopy(self._hsp)
        self._hsp.query, self._hsp.sbjct = self._hsp.sbjct, self._hsp.query
        self._hsp.query_start, self._hsp.sbjct_start =\
            self._hsp.sbjct_start, self._hsp.query_start

    # end of function _swap_input_hsp_qands


    def relatively_positioned_towards(self,pacbp):
        """ """
        return ordering.relatively_positioned_towards(self,pacbp)

    # end of function relatively_positioned_towards


    def is_positioned_compatibly(self,pacbp):
        """ """
        pos = self.relatively_positioned_towards(pacbp)
        (q1,s1) = ordering.relpos2binaryrelpos(pos['Q1'],pos['S1'])
        (q2,s2) = ordering.relpos2binaryrelpos(pos['Q2'],pos['S2'])
        if q1 != s1:   return False
        elif q2 != s2: return False
        else:          return True

    # end of function is_positioned_compatibly


    def is_postioned_compatibly(self,pacbp):
        """ Type function name; keep intact untill fixed """
        return self.is_positioned_compatibly(pacbp)

    # end of function is_postioned_compatibly


    def distance_towards(self,pacbp):
        """
        Give the (AA) distance between two pacbps

        @attention: see pacb.ordering.distance_towards for documentation

        @rtype:  positive integer or zero
        @return: AA distance between both pacbps
        """
        return ordering.distance_towards(self,pacbp)

    # end of function distance_towards 


    def overlap(self,pacbp):
        """
        Give the overlap between this and another pacbp

        @attention: see pacb.ordering.overlap for documentation

        @rtype:  positive float 
        @return: highest overlap ratio of query and sbjct of these 2 pacbps
        """
        return ordering.overlap(self,pacbp)

    # end of function overlap


    def query_overlap(self,pacbp):
        """
        Give the overlap between the query of this and another pacbp

        @attention: see pacb.ordering.query_overlap for documentation

        @rtype:  positive float 
        @return: highest overlap ratio of overlap of query sequences
        """
        return ordering.query_overlap(self,pacbp)

    # end of function query_overlap


    def sbjct_overlap(self,pacbp):
        """
        Give the overlap between the sbjcts of this and another pacbp

        @attention: see pacb.ordering.sbjct_overlap for documentation

        @rtype:  positive float 
        @return: highest overlap ratio of overlap of sbjct sequences
        """
        return ordering.sbjct_overlap(self,pacbp)

    # end of function sbjct_overlap


    def togff(self,gff={}):
        """
        Return 8-element tuple of gff data.
        To be overwritten in the subclasses.
        """
        pass

    # end of function togff


    def check(self):
        """
        Function that checks if attribute values are (still) as expected.
        Object can be changed (coordinate-shuffeling, splitting, extending).
        After such events, a call to check() confirms the validity of data.
        """
        _tmp = []
        for attr in ['query','match','sbjct','alignment']:
            _tmp.append( len(getattr(self,attr)) )
        _tmp.append( self.length )
        if len(Set(_tmp)) != 1:
            # non-identical lengths of alignment strings
            return False
            raise "non-identical lengths of alignment strings"
        if self.query_end-self.query_start != self.query_len:
            # non-identical length of query string and coords
            return False
            raise "non-identical length of query string and coords"
        if self.sbjct_end-self.sbjct_start != self.sbjct_len:
            # non-identical length of sbjct string and coords
            return False
            raise "non-identical length of sbjct string and coords"
        for attr in dir(self):
            if attr[0] == "_": continue
            if type(getattr(self,attr)) == type(int()):
                if getattr(self,attr) < 0:
                    # all coordinates should be >= 0
                    return False
                    raise "coordindate < 0: '%s'" % attr
        # all tests succesfull!
        return True

    # end of function check


    def strip_unmatched_ends(self):
        """
        """
        if self.alignment:
            # okay, now strip the alignment from non-matched chars,
            # means leading and trailing spaces!
            while self.alignment and self.alignment[0] == " ":
                self._strip_one_left_position()
            while self.alignment and self.alignment[-1] == " ":
                self._strip_one_rigth_position()

    # end of function strip_unmatched_ends


    def strip_tailing_gaps(self):
        """
        """
        if self.alignment:
            # okay, now strip the alignment from leading and tailing gaps
            while self.alignment and (self.query[0] == "-" or self.sbjct[0] == "-"):
                self._strip_one_left_position()
            while self.alignment and (self.query[-1] == "-" or self.sbjct[-1] == "-"): 
                self._strip_one_rigth_position()

    # end of function strip_tailing_gaps 


    def strip_consistent_internal_gaps(self):
	"""
	Strip consistently aligned gaps present in both query and sbjct
	
	@attention: can originate when constructing pairwise from multiple alignments
	"""
	is_removed = False
	for pos in range(len(self.alignment)-1,-1,-1):
	    if not self.alignment:
		# Zerosized alignment; completely consumed
		break
		is_removed = True
	    if set([ self.query[pos],self.sbjct[pos] ]) == set("-"):
		# remove this position from the aligned sequences
		self.query = self.query[0:pos] + self.query[pos+1:]
		self.match = self.match[0:pos] + self.match[pos+1:]
		self.sbjct = self.sbjct[0:pos] + self.sbjct[pos+1:]
	        self.length -= 1
		is_removed = True
		# forwards compatibilty with PacbPDNA and PacbPORF objects
		# remove position from self._positions
		if hasattr(self,"_positions"): popped = self._positions.pop(pos)
	if is_removed:
	    # (re) score the aligment
	    self._score_alignment()
	    
    # end of function strip_consistent_internal_gaps


    def _strip_one_left_position(self):
        """
        """
        # no truncation possible -> PacbP does not exist anymore
        if self.length == 0: raise ZeroSizedPacb

        # correct coordinates and length
        if self.query[0] != "-":
            self.query_start+=1
            self.query_protein  = self.query_protein[1:]
        if self.sbjct[0] != "-":
            self.sbjct_start+=1
            self.sbjct_protein  = self.sbjct_protein[1:]

        # decrease length attribute
        self.length -= 1
        # get rid of FIRST char of all strings
        self.alignment      = self.alignment[1:]
        self.query          = self.query[1:]
        self.match          = self.match[1:]
        self.sbjct          = self.sbjct[1:]

        # and (re) score the aligment
        self._score_alignment()

    # end of function _strip_one_left_position


    def _strip_one_rigth_position(self):
        """
        """
        # no truncation possible -> PacbP does not exist anymore
        if self.length == 0: raise ZeroSizedPacb

        # correct coordinates and length
        if self.query[-1] != "-":
            self.query_end-=1
            self.query_protein  = self.query_protein[0:-1]
        if self.sbjct[-1] != "-":
            self.sbjct_end-=1
            self.sbjct_protein  = self.sbjct_protein[0:-1]

        # decrease length attribute
        self.length -= 1

        # get rid of LAST char of all strings
        self.alignment      = self.alignment[0:-1]
        self.query          = self.query[0:-1]
        self.match          = self.match[0:-1]
        self.sbjct          = self.sbjct[0:-1]

        # and (re) score the aligment
        self._score_alignment()

    # end of function _strip_one_rigth_position


    def slice(self, abs_start, abs_end, coords_on='query'):
        """
        get only a slice of current pacbp
        abs_start and abs_end are ABSOLUTE PROTEIN coordinates
        e.g. self.query_start = 123, self.query_end = 212
        if (abs_start, abs_end) == (154, 211), then
        31 left positions and 1 rigth position will be removed 
        
        if a set of conflicting coordinates is given (end <= start)
        a CoordinateOutOfRange exception is raised
        as soon as a 0-sized slice is about to be created, 
        a ZeroSizedPacb exception is raised
        """
        # check start and end coordinates
        if abs_end <= abs_start: raise CoordinateOutOfRange

        if coords_on == 'sbjct':
            while self.sbjct_start < abs_start:
                if self.length == 1: raise ZeroSizedPacb
                self._strip_one_left_position()
            while self.sbjct_end > abs_end:
                if self.length == 1: raise ZeroSizedPacb
                self._strip_one_rigth_position()
        else:
            while self.query_start < abs_start:
                if self.length == 1: raise ZeroSizedPacb
                self._strip_one_left_position()
            while self.query_end > abs_end:
                if self.length == 1: raise ZeroSizedPacb
                self._strip_one_rigth_position()

    # end of function slice


    def returnslice(self, abs_start, abs_end, coords_on='query'):
        """
        """
        new = deepcopy(self)
        new.slice(abs_start,abs_end,coords_on=coords_on)
        return new

    # end of function returnslice


    def pacbp_length_discrepancy(self):
        """
        Absolute length difference of the query and sbjct protein sequences

        @rtype:  positive integer
        @return: Absolute difference of alignment length (counted in amino acids)
        """
        ( qlen, slen ) = self._pacbp_lengths()
        return abs( slen - qlen )

    # end of pacbp_length_discrepancy


    def pacbp_relative_length_discrepancy(self):
        """
        Relative length difference of the query and sbjct protein sequences

        @rtype:  float
        @return: Relative difference of alignment length (smallest/longest)
        """
        ( qlen, slen ) = self._pacbp_lengths()
        return float( min([qlen, slen]) ) / float( max([qlen, slen]) )

    # end of pacbp_length_discrepancy


    def _pacbp_lengths(self):
        """
        Length of the originally aligned query and sbjct protein sequences

        @rtype:  tuple
        @return: ( protein_query_length, protein_sbjct_length )
        """
        if self.__class__.__name__ == 'PacbPORF':
            # PacbPORF or further inhertited object.
            # The alignment might be extended,so
            # so take original alignemnt coordinates
            spos = self._positions[self._original_alignment_pos_start]
            epos = self._positions[self._original_alignment_pos_end-1]
            # correct for +1 beacause coord is a python list slice coord!
            qlen = epos.query_pos+1 - spos.query_pos
            slen = epos.sbjct_pos+1 - spos.sbjct_pos
            return ( slen , qlen )
        else:
            # PacbP or PacbPDNA; no extension has taken place
            return ( self.query_length , self.sbjct_length )

    # end of function _pacbp_lengths


    def potentially_contains_aligned_intron(self,gap_size=5,length_discrepancy=10,gap_size_alone=10):
        """
        Does PacbP potentially contains an aligned inframe intron?

        @type  length_discrepancy: Integer (positive)
        @param length_discrepancy: Length of minimal length difference between query and sbjct

        @type  gap_size: Integer (positive)
        @param gap_size: Length of minimal continious gap length

        @type  gap_size_alone: Integer (positive)
        @param gap_size_alone: Length of minimal continious gap length alone to return a True

        @rtype:  Boolean
        @return: True or False

        @attention: both gap_size AND length_discrepancy must be furfilled to return a True
        """
        if self.pacbp_length_discrepancy() >= length_discrepancy and\
        self.alignment_has_gaps(gap_size=gap_size):
            return True
        elif self.alignment_has_gaps(gap_size=gap_size_alone):
            return True
        else:
            return False

    # end of function potentially_contains_aligned_intron


    def alignment_protein_range_query(self):
        """ """
        return range(self.query_start,self.query_end)

    # end of function alignment_protein_range_query


    def alignment_protein_range_sbjct(self):
        """ """
        return range(self.sbjct_start,self.sbjct_end)

    # end of function alignment_protein_range_sbjct


    def issubsetorsuperset(self,otherpacb):
        """ """
        return ordering.issubsetorsuperset(self,otherpacb)

    # end of function issubsetorsuperset


    def issubset(self,otherpacb):
        """ """
        return ordering.issubset(self,otherpacb)

    # end of function issubset


    def issuperset(self,otherpacb):
        """ """
        return ordering.issuperset(self,otherpacb)

    # end of function issuperset


    def alignmentposition_by_query_pos(self,qposPY,forced_return=False):
        """
        Returns the (integer) position in the alignment of the requested query protein position

        @attention: qposPY is ABSOLUTE position on query input protein sequence

        @type  qposPY: positive integer
        @param qposPY: protein position on query input protein sequence

        @type  forced_return: Boolean
        @param forced_return: ...

        @rtype:  integer
        @return: position in the alignemnt of this query protein position
        """
        if qposPY < self.query_start:
            if forced_return:   raise "TO DO... qposPY < self.query_start"
            else:               raise CoordinateOutOfRange
        elif qposPY >= self.query_end:
            if forced_return:   raise "TO DO... qposPY >= self.query_end"
            else:               raise CoordinateOutOfRange
        else:
            qpos = deepcopy(self.query_start)
            apos = 0
            for q in self.query:
                if qpos == qposPY: break
                # increase positional counters
                if q != '-':qpos+=1
                apos+=1                  
            # return the alignment position
            return apos

    # end if function alignmentposition_by_query_pos


    def alignmentposition_by_sbjct_pos(self,sposPY,forced_return=False):
        """
        Returns the (integer) position in the alignment of the requested sbjct protein position

        @attention: sposPY is ABSOLUTE position on sbjct input protein sequence

        @type  sposPY: positive integer
        @param sposPY: protein position on sbjct input protein sequence

        @type  forced_return: Boolean
        @param forced_return: ...

        @rtype:  integer
        @return: position in the alignemnt of this sbjct protein position
        """
        if sposPY < self.sbjct_start:
            if forced_return:   raise "TO DO... sposPY < self.sbjct_start"
            else:               raise CoordinateOutOfRange
        elif sposPY >= self.sbjct_end:
            if forced_return:   raise "TO DO... sposPY >= self.sbjct_end"
            else:               raise CoordinateOutOfRange
        else:
            spos = deepcopy(self.sbjct_start)
            apos = 0
            for s in self.sbjct:
                if spos == sposPY: break
                # increase positional counters
                if s != '-':spos+=1
                apos+=1                  
            # return the alignment position
            return apos

    # end if function alignmentposition_by_sbjct_pos

    def get_unextended_length(self):
        """ Forwards compatiblity with PacbPORF """
        return int(self.length)

    # end of function get_unextended_length



    def get_aligned_protein_sequences(self):
        """ Forwards compatiblity with PacbPORF """
        return (self.query,self.match,self.sbjct)

    # end of function get_aligned_protein_sequences


    def get_unextended_aligned_protein_sequences(self):
        """ Forwards compatiblity with PacbPORF """
        return (self.query,self.match,self.sbjct)

    # end of function get_unextended_aligned_protein_sequences


    def has_query_methionine(self):
        """ """
        _q,_m,_s = self.get_unextended_aligned_protein_sequences()
        if _q.count("M")+_q.count("m") > 0: return True
        else:                               return False

    # end of function has_query_methionine


    def has_sbjct_methionine(self):
        """ """
        _q,_m,_s = self.get_unextended_aligned_protein_sequences()
        if _s.count("M")+_s.count("m") > 0: return True
        else:                               return False

    # end of function has_sbjct_methionine


    def has_methionines(self):
        """ """
        _q,_m,_s = self.get_unextended_aligned_protein_sequences()
        if _s.count("M")+_s.count("m") == 0:    return False
        elif _q.count("M")+_q.count("m") == 0:  return False
        else:                                   return True

    # end of function has_methionines

# end of class PacbP


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


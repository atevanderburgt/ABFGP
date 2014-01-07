"""
"""

# Module metadata
__authors__ = "Ate van der Burgt"
__license__ = "MIT"


# Imports from pacb package
from pacbp import PacbP
from alignedproteindnaposition import *
from exceptions import * 
from validators import *

# Python Imports
from sets import Set

# Abgp Imports
from pythonlibs.uniqueness import get_random_string_tag
from lib_sequenceperiodicity import coding_sequence_periodicity_test
from lib_clustalw import clustalw
from lib_proteinsimilaritymatrix import ProteinSimilarityMatrix

DEFAULT_MATRIX = ProteinSimilarityMatrix(name="BLOSUM62")


class PacbPDNA(PacbP):
    # Pairwise-Aligned Coding Block PROTEIN and DNA; something like an exon
    # this `flavour` has an aligned Protein sequence and
    # its corresponding DNA sequence
    def __init__(self,input=("","",0,0),blastp_hsp=None,MATRIX={},
        pacbp=None,input_dna=("","",0,0)):
        """
        """
        # first, get to the point of an instantiated PacbP object
        if str(pacbp.__class__) in ['pacb.PacbP','pacb.pacbp.PacbP']:
            # loop over all attributes, and copy all true attributes
            # skip the `__`-starting attributes and the `instancemethods`
            for attr in dir(pacbp):
                if attr[0:2] == "__": continue
                if str(type(getattr(pacbp,attr))) != "<type 'instancemethod'>":
                    # copy this attribute from PacbP to PacbPDNA
                    setattr(self,attr,getattr(pacbp,attr))
            # and move the original hsp object
            self._hsp = pacbp._hsp
            # check if a non-default matrix was applied
            if MATRIX:
                self.MATRIX = MATRIX
        else:
            # call the __init__ function of the basal PacbP class
            if MATRIX:
                PacbP.__init__(self,input=input,blastp_hsp=blastp_hsp,MATRIX=MATRIX)
            else:
                PacbP.__init__(self,input=input,blastp_hsp=blastp_hsp)

        # set the PacbPDNA specific attributes
        self.query_dna          = ""
        self.query_dna_start    = 0
        self.query_dna_end      = 0
        self.query_dna_length   = 0
        self.sbjct_dna          = ""
        self.sbjct_dna_start    = 0
        self.sbjct_dna_end      = 0
        self.sbjct_dna_length   = 0
        self._positions         = []
        # and finish initialization
        self._initialize(input_dna)
        self._alignment2positions()
        # set start and end positions of the original alignment
        # this is needed to calculate the `entropy` of a position
        # in the alignment
        self._original_alignment_pos_start = 0
        self._original_alignment_pos_end   = len(self._positions)


    # end of function __init__


    def __str__(self):
        """
        """
        source = ' '
        if self.source != 'blastp':
            source = " src=%s " % self.source
        return """<pacbpdna%sQ:%s-%s S:%s-%s bits:%s iden:%1.3f len:(%s,%s)>""" % (
            source,
            self.query_start, self.query_end,
            self.sbjct_start, self.sbjct_end,
            self.bitscore, self.identityscore, 
            self.query_length, self.sbjct_length )

    # end of function __str__


    def full_report(self):
        """
        Usefull for debugging. A blueprint of the object
        """
        print "#$"*30
        print "#$"*30
        print self
        print self.print_protein_and_dna(rel_coords=(0,self.length))
        for attr in dir(self):
            item = getattr(self,attr)
            if str(item)[0:13] == "<bound method": continue
            printitem = item
            if len(str(item)) > 80 and str(item)[0]+str(item)[-1] != "<>":
                printitem = "len(%s)" % len(str(item))
            print attr, "\t", printitem
        print "#$"*30
        print "#$"*30

    # end of function full_report


    def _initialize(self,input_dna):
        """
        """
        # TODO TODO so some input integrity checks
        # self._input_validation(input_dna)

        # set DNA sequences
        self.query_dna          = input_dna[0]
        self.sbjct_dna          = input_dna[1]
        # coordinates of the DNA sequences
        self.query_dna_start    = input_dna[2]
        self.query_dna_length   = len(self.query_dna)
        self.query_dna_end      = self.query_dna_start + self.query_dna_length
        self.sbjct_dna_start    = input_dna[3]
        self.sbjct_dna_length   = len(self.sbjct_dna)
        self.sbjct_dna_end      = self.sbjct_dna_start + self.sbjct_dna_length

    # end of function _initialize


    def __len__(self):
        """
        """
        return self._get_original_alignment_pos_end().query_pos - self._get_original_alignment_pos_start().query_pos + 1 

    # end of function __len__


    def strip_unmatched_ends(self):
        """ OVERWRITES PacbP.strip_unmatched_ends() """
        # strip 5' side of the alignment
        spos = self._get_original_alignment_pos_start()
        epos = self._get_original_alignment_pos_end()
        while spos.position < epos.position:
            if spos.match == " ":
                self._original_alignment_pos_start+=1
                spos = self._get_original_alignment_pos_start()
            else:
                break

        # strip 3' side of the alignment
        spos = self._get_original_alignment_pos_start()
        epos = self._get_original_alignment_pos_end()
        while epos.position > spos.position:
            if epos.match == " ":
                self._original_alignment_pos_end-=1
                epos = self._get_original_alignment_pos_end()
            else:
                break

    # end of function strip_unmatched_ends


    def _get_original_alignment_pos_start(self):
        """
        """
        return self._positions[self._original_alignment_pos_start]

    # end of _get_original_alignment_pos_start


    def _get_original_alignment_pos_end(self):
        """
        """
        return self._positions[self._original_alignment_pos_end-1]

    # end of _get_original_alignment_pos_end


    def get_unextended_length(self):
        """ """
        return self._original_alignment_pos_end-self._original_alignment_pos_start

    # end of function get_unextended_length


    def get_at_ratio_averaged(self):
        """ """
        Qseq,Sseq = self.get_unextended_aligned_dna_sequences()
        atQ = Qseq.count("A")+Qseq.count("a")+Qseq.count("T")+Qseq.count("t")
        atS = Sseq.count("A")+Sseq.count("a")+Sseq.count("T")+Sseq.count("t")
        return float(atQ+atS)/(len(Qseq)+len(Sseq))

    # end of function get_at_ratio_averaged


    def get_gc_ratio_averaged(self):
        """ """
        Qseq,Sseq = self.get_unextended_aligned_dna_sequences()
        gcQ = Qseq.count("G")+Qseq.count("g")+Qseq.count("C")+Qseq.count("c")
        gcS = Sseq.count("G")+Sseq.count("g")+Sseq.count("C")+Sseq.count("c")
        return float(gcQ+gcS)/(len(Qseq)+len(Sseq))

    # end of function get_gc_ratio_averaged


    def get_at_ratio_query(self):
        """ """
        seq,Sseq = self.get_unextended_aligned_dna_sequences()
        at = seq.count("A")+seq.count("a")+seq.count("T")+seq.count("t")
        return float(at)/len(seq)

    # end of function get_at_ratio_query


    def get_at_ratio_sbjct(self):
        """ """
        Qseq,seq = self.get_unextended_aligned_dna_sequences()
        at = seq.count("A")+seq.count("a")+seq.count("T")+seq.count("t")
        return float(at)/len(seq)

    # end of function get_at_ratio_sbjct


    def get_gc_ratio_query(self):
        """ """
        seq,Sseq = self.get_unextended_aligned_dna_sequences()
        gc = seq.count("G")+seq.count("g")+seq.count("C")+seq.count("c")
        return float(gc)/len(seq)

    # end of function get_gc_ratio_query


    def get_gc_ratio_sbjct(self):
        """ """
        Qseq,seq = self.get_unextended_aligned_dna_sequences()
        gc = seq.count("G")+seq.count("g")+seq.count("C")+seq.count("c")
        return float(gc)/len(seq)

    # end of function get_gc_ratio_sbjct


    def alignment_dna_range_query(self):
        """ """
        return range(self.query_dna_start,self.query_dna_end)

    # end of function alignment_dna_range_query


    def alignment_dna_range_sbjct(self):
        """ """
        return range(self.sbjct_dna_start,self.sbjct_dna_end)

    # end of function alignment_dna_range_sbjct


    def _alignment2positions(self):
        """
        """
        query_pos     = self.query_start
        sbjct_pos     = self.sbjct_start
        query_dna_pos = self.query_dna_start
        sbjct_dna_pos = self.sbjct_dna_start

        # loop over each AA position in the alignment
        for pos in range(0,self.length):
            query_dna = ("",0)
            sbjct_dna = ("",0)
            # get AA-symbol from the alignment
            q = self.query[pos]
            s = self.sbjct[pos]
            m = self.match[pos]
            if q != "-":
                # query AA is not a gap, so get it's DNA sequence counterpart
                rel_query_dna_pos = query_dna_pos - self.query_dna_start
                query_dna = (self.query_dna[rel_query_dna_pos:rel_query_dna_pos+3],query_dna_pos)
            if s != "-":
                # sbjct AA is not a gap, so get it's DNA sequence counterpart
                rel_sbjct_dna_pos = sbjct_dna_pos - self.sbjct_dna_start
                sbjct_dna = (self.sbjct_dna[rel_sbjct_dna_pos:rel_sbjct_dna_pos+3],sbjct_dna_pos)
            # and append the position!
            self._positions.append( AlignedProteinDnaPosition(
                            pos, query_pos, sbjct_pos,
                            alignment=(q,m,s),
                            query_dna = query_dna,
                            sbjct_dna = sbjct_dna
                            ) )
            # check if we must increase the counters
            if q != "-":
                query_pos+=1
                query_dna_pos+=3
            if s != "-":
                sbjct_pos+=1
                sbjct_dna_pos+=3

    # end of function _alignment2positions


    def _strip_one_left_position(self):
        """
        """
        # no truncation possible -> PacbP does not exist anymore
        if self.length == 0: raise ZeroSizedPacb

        if self.query[0] != "-":
            self.query_dna          = self.query_dna[3:]
            self.query_dna_start    = self.query_dna_start+3
            self.query_dna_length   = self.query_dna_length-3

        if self.sbjct[0] != "-":
            self.sbjct_dna          = self.sbjct_dna[3:]
            self.sbjct_dna_start    = self.sbjct_dna_start+3
            self.sbjct_dna_length   = self.sbjct_dna_length-3

        # increase `position` counter for all existing positions
        for pos in self._positions: pos.position-=1

        # and pop leading aligned position
        self._positions.pop(0)

        # COMPENSATE self._original_alignment_pos_start and end
        if self._original_alignment_pos_start > 0: 
            self._original_alignment_pos_start-=1
            self._original_alignment_pos_end-=1
        elif self._original_alignment_pos_start==0:
            self._original_alignment_pos_end-=1
        else:
            raise "WHATEVER; coords for original alignment not oke"

        # correct all the protein coordinates and strings
        PacbP._strip_one_left_position(self)

    # end of function _strip_one_left_position


    def _strip_one_rigth_position(self):
        """
        """
        # no truncation possible -> PacbP does not exist anymore
        if self.length == 0: raise ZeroSizedPacb

        if self.query[-1] != "-":
            self.query_dna          = self.query_dna[0:-3]
            self.query_dna_end      = self.query_dna_end-3
            self.query_dna_length   = self.query_dna_length-3

        if self.sbjct[-1] != "-":
            self.sbjct_dna          = self.sbjct_dna[0:-3]
            self.sbjct_dna_end      = self.sbjct_dna_end-3
            self.sbjct_dna_length   = self.sbjct_dna_length-3

        # COMPENSATE self._original_alignment_pos_end
        if len(self._positions) == self._original_alignment_pos_end:
            self._original_alignment_pos_end-=1
        elif len(self._positions) > self._original_alignment_pos_end:
            pass
        else:
            raise "WHATEVER; coords for original alignment not oke"

        # and pop trailing aligned position
        self._positions.pop()

        # correct all the protein coordinates and strings
        PacbP._strip_one_rigth_position(self)

    # end of function _strip_one_rigth_position


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
        # first swap the data from the basal PacbP clas
        PacbP._swap_query_and_sbjct(self)

        # swap strings and orf objects
        self.query_dna, self.sbjct_dna = self.sbjct_dna, self.query_dna

        # swap coordinates
        (a,b,c) = (self.query_dna_start, self.query_dna_end, self.query_dna_length)
        (d,e,f) = (self.sbjct_dna_start, self.sbjct_dna_end, self.sbjct_dna_length)
        (self.query_dna_start, self.query_dna_end, self.query_dna_length) = (d,e,f)
        (self.sbjct_dna_start, self.sbjct_dna_end, self.sbjct_dna_length) = (a,b,c)

        # and swap the alignment on all positions
        for adp_pos in self._positions:
            adp_pos._swap_query_and_sbjct()

    # end of function _swap_query_and_sbjct


    def dnaposition_query(self,abs_dna_pos,forced_return=False):
        """
        Find the position of the AlignedProteinDnaPosition by absolute query DNA/nt position

        @attention: abs_dna_pos is ABSOLUTE position on query input genomic sequence

        @type  abs_dna_pos: positive integer
        @param abs_dna_pos: nt position on query input genomic sequence

        @type  forced_return: Boolean
        @param forced_return: True or False. If False and abs_dna_pos not in alignment,
                              a CoordinateOutOfRange exception is raised.
                              If True, abs_dna_pos is recalculated to as-if int the alignment

        @rtype:  tuple
        @return: (int_A,int_B) with int_A: position number, int_B: phase [012]
        """
        for apdpos in self._positions:
            ntpositions = range(apdpos.query_dna_start,apdpos.query_dna_end)
            if abs_dna_pos in ntpositions:
                #print "returningA", apdpos.position, abs_dna_pos, range(apdpos.query_dna_start,apdpos.query_dna_end),
                #print self.query_dna_start, self.query_dna_end, self._original_alignment_pos_start,
                #print self._original_alignment_pos_end, len(self._positions),
                #print self._original_alignment_pos_end - self._original_alignment_pos_start
                return ( apdpos.position, ntpositions.index(abs_dna_pos) )
        else:
            # Abs_dna_pos argument is not found in the for loop
            # In case of a forced_return, create a position as if
            # if exists. If no forced_return -> CoordinateOutOfRange
            if forced_return:
                # recalculate to a position as if in the alignment
                if abs_dna_pos < self.query_dna_start:
                    # abs_dna_pos is in front of the alignment
                    pos   = abs_dna_pos - self.query_dna_start
                    phase = pos % 3
                    pos   = pos / 3
                    #print "returningB"
                    return (pos,phase)
                else:
                    # abs_dna_pos is after the alignment
                    pos   = abs_dna_pos - (self.query_dna_end-1)
                    phase = (pos + 2) % 3
                    pos   = (pos / 3) + len(self._positions)-1
                    if phase != 2: pos+=1
                    #print "returningC"
                    return (pos,phase)
            else:
                #return ( -1, -1 )
                # requested position not in aligned pacbp
                #print "### ABOUT TO raise CoordinateOutOfRange"
                #print abs_dna_pos, "asked for, range is: ", self._positions[0].query_dna_start, "..", self._positions[-1].query_dna_end,
                #print "or", self.query_dna_start, "..", self.query_dna_end
                #print self
                #print self.print_protein_and_dna(rel_coords=(0,self.length))
                #print "query AA  range:", [ pos.query_pos for pos in self._positions ], len([ pos.query_pos for pos in self._positions ]), len(Set([ pos.query_pos for pos in self._positions ]))
                #print "query DNA range:", [ pos.query_dna_start for pos in self._positions ], len([ pos.query_dna_start for pos in self._positions ]), len(Set([ pos.query_dna_start for pos in self._positions ]))
                raise CoordinateOutOfRange

    # end of function dnaposition_query


    def dnaposition_sbjct(self,abs_dna_pos,forced_return=False):
        """
        Find the position of the AlignedProteinDnaPosition by absolute sbjct DNA/nt position

        @attention: abs_dna_pos is ABSOLUTE position on sbjct input genomic sequence

		@type  abs_dna_pos: positive integer
		@param abs_dna_pos: nt position on sbjct input genomic sequence

		@type  forced_return: Boolean
		@param forced_return: True or False. If False and abs_dna_pos not in alignment,
                              a CoordinateOutOfRange exception is raised.
                              If True, abs_dna_pos is recalculated to as-if int the alignment

		@rtype:  tuple
		@return: (int_A,int_B) with int_A: position number, int_B: phase [012]
        """
        for apdpos in self._positions:
            ntpositions = range(apdpos.sbjct_dna_start,apdpos.sbjct_dna_end)
            if abs_dna_pos in ntpositions:
                return ( apdpos.position, ntpositions.index(abs_dna_pos) )
        else:
            if forced_return:
                # recalculate to a position as if in the alignment
                if abs_dna_pos < self.sbjct_dna_start:
                    # abs_dna_pos is in front of the alignment
                    pos   = abs_dna_pos - self.sbjct_dna_start
                    phase = pos % 3
                    pos   = pos / 3
                    return (pos,phase)
                else:
                    # abs_dna_pos is after the alignment
                    pos   = abs_dna_pos - (self.sbjct_dna_end-1)
                    phase = (pos + 2) % 3
                    pos   = (pos / 3) + len(self._positions)-1
                    if phase != 2: pos+=1
                    return (pos,phase)
            else:
                #return ( -1, -1 )
                # requested position not in aligned pacbp
                #print abs_dna_pos, "asked for, range is: ", self._positions[0].sbjct_dna_start, "..", self._positions[-1].sbjct_dna_end, "or", self.sbjct_dna_start, "..", self.sbjct_dna_end
                raise CoordinateOutOfRange

    # end of function dnaposition_sbjct


    def get_distance_aligned_protein_positions(self,query=None,sbjct=None):
        """
        Get the distance between absolute query and sbjct amino acid positions.

		@type  query: positive integer
		@param query: absolute Amino-Acid position of query position

		@type  sbjct: positive integer
		@param sbjct: absolute Amino-Acid position of query position

		@rtype:  positive integer
		@return: distance in aligned amino acid positions
        """
        # some error checks for non-me users ;-)
        if query==None or sbjct==None:
            message = "specify both query=.. and sbjct=.. arguments"
            raise InproperlyAppliedArgument, message
        if not type(query) == type(int()):
            message = "specify `query` as an integer, not '%s'" % type(query)
            raise InproperlyAppliedArgument, message
        if not type(sbjct) == type(int()):
            message = "specify `sbjct` as an integer, not '%s'" % type(sbjct)
            raise InproperlyAppliedArgument, message

        # get relative query position
        if query < self.query_start:
            algQpos = query - self.query_start + self._positions[0].position
        elif query >= self.query_end:
            algQpos = query - self.query_end + self._positions[-1].position
        else:
            algQpos = self.alignmentposition_by_query_pos(query)

        # get relative sbjct position
        if sbjct < self.sbjct_start:
            algSpos = sbjct - self.sbjct_start + self._positions[0].position
        elif sbjct >= self.sbjct_end:
            algSpos = sbjct - self.sbjct_end + self._positions[-1].position
        else:
            algSpos = self.alignmentposition_by_sbjct_pos(sbjct)

        # and return (absolute) distance between relative positions
        return abs( algQpos - algSpos )

    # end of function get_distance_aligned_protein_positions


    def get_distance_aligned_nucleotide_positions(self,query=None,sbjct=None):
        """
        Get the distance between absolute query and sbjct nucleotide positions.

		@type  query: positive integer
		@param query: absolute nucleotide position of query position

		@type  sbjct: positive integer
		@param sbjct: absolute nucleotide position of query position

		@rtype:  positive integer
		@return: distance in aligned nucleotide positions
        """
        # some error checks for non-me users ;-)
        if query==None or sbjct==None:
            message = "specify both query=.. and sbjct=.. arguments"
            raise InproperlyAppliedArgument, message
        if not type(query) == type(int()):
            message = "specify `query` as an integer, not '%s'" % type(query)
            raise InproperlyAppliedArgument, message
        if not type(sbjct) == type(int()):
            message = "specify `sbjct` as an integer, not '%s'" % type(sbjct)
            raise InproperlyAppliedArgument, message

        # get relative query position
        if query < self.query_dna_start:
            algQpos = query - self.query_dna_start + self._positions[0].position*3
        elif query >= self.query_dna_end:
            algQpos = query - self.query_dna_end + self._positions[-1].position*3
        else:
            (pos,phase) = self.dnaposition_query(query)
            algQpos = pos*3 +phase

        # get relative sbjct position
        if sbjct < self.sbjct_dna_start:
            algSpos = sbjct - self.sbjct_dna_start + self._positions[0].position*3
        elif sbjct >= self.sbjct_dna_end:
            algSpos = sbjct - self.sbjct_dna_end + self._positions[-1].position*3
        else:
            (pos,phase) = self.dnaposition_sbjct(sbjct)
            algSpos = pos*3 +phase

        # and return (absolute) distance between relative positions
        return abs( algQpos - algSpos )

    # end of function get_distance_aligned_nucleotide_positions


    def alignmentposition_by_query_pos(self,qposPY,forced_return=False):
        """
        Returns the AlignedProteinDnaPosition of requested query protein position

        @attention: qposPY is ABSOLUTE position on query input protein sequence

		@type  qposPY: positive integer
		@param qposPY: protein position on query input protein sequence

		@rtype:  integer
        @return: relative position (AlignedProteinDnaPosition) of this query protein position
        """
        for offset in range(0,len(self._positions)):
            pos = self._positions[offset]
            if pos.query_pos == qposPY and pos.query != "-":
                return offset
        else:
            if forced_return:
                # recalculate to a position as if in the alignment
                if qposPY < self.query_start:
                    # qposPY is in front of the alignment
                    return qposPY - self.query_start
                else:
                    # qposPY is after the alignment
                    return qposPY - (self.query_end-1) + len(self._positions) - 1
            else:
                # requested position not in aligned pacbp
                #print "### ABOUT TO raise CoordinateOutOfRange"
                #print qposPY, "asked for, range is: ", self._positions[0].query_pos, "..", self._positions[-1].query_pos, "or", self.query_start, "..", self.query_end
                #print "### ABOUT TO raise CoordinateOutOfRange"
                #print self.full_report()
                #print "### ABOUT TO raise CoordinateOutOfRange"
                #print qposPY, "asked for, range is: ", self._positions[0].query_pos, "..", self._positions[-1].query_pos, "or", self.query_start, "..", self.query_end
                #print "### ABOUT TO raise CoordinateOutOfRange"

                #raise CoordinateOutOfRange
                return CoordinateOutOfRange

    # end of function alignmentposition_by_query_pos


    def alignmentposition_by_sbjct_pos(self,sposPY,forced_return=False):
        """
        Returns the AlignedProteinDnaPosition of requested sbjct protein position

        @attention: sposPY is ABSOLUTE position on sbjct input protein sequence

		@type  sposPY: positive integer
		@param sposPY: protein position on sbjct input protein sequence

		@rtype:  integer
        @return: relative position (AlignedProteinDnaPosition) of this sbjct protein position
        """
        for offset in range(0,len(self._positions)):
            pos = self._positions[offset]
            if pos.sbjct_pos == sposPY and pos.sbjct != "-":
                return offset
        else:
            if forced_return:
                # recalculate to a position as if in the alignment
                if sposPY < self.sbjct_start:
                    # sposPY is in front of the alignment
                    return sposPY - self.sbjct_start
                else:
                    # sposPY is after the alignment
                    return sposPY - (self.sbjct_end-1) + len(self._positions) - 1
            else:
                # requested position not in aligned pacbp
                ###print sposPY, "asked for, range is: ", self._positions[0].sbjct_pos, "..", self._positions[-1].sbjct_pos, "or", self.sbjct_start, "..", self.sbjct_end
                ###raise CoordinateOutOfRange
                return CoordinateOutOfRange

    # end of function alignmentposition_by_sbjct_pos


    def alignmentobject_by_queryaa(self,qposPY):
        """
        Returns the AlignedProteinDnaPosition object of requested ABSOLUTE query protein position

        @attention: qposPY is ABSOLUTE position on query input protein sequence

		@type  qposPY: positive integer
		@param qposPY: protein position on query input protein sequence

		@rtype:  AlignedProteinDnaPosition
        @return: AlignedProteinDnaPosition object
        """
        if self._positions[0].query_pos > qposPY:
            # requested position in front of (current) (extended) alignment 
            return CoordinateOutOfRange
        if self._positions[-1].query_pos < qposPY:
            # requested position after of (current) (extended) alignment 
            return CoordinateOutOfRange
        # loop over the AlignedProteinDnaPosition objects and return
        # it as soon as the coordinates match
        for offset in range(0,len(self._positions)):
            posobj = self._positions[offset]
            if posobj.query_pos == qposPY and posobj.query != "-":
                return posobj
        else:
            # position not in the (current) (extended) alignment
            return CoordinateOutOfRange

    # end of function alignmentobject_by_queryaa


    def alignmentobject_by_sbjctaa(self,sposPY):
        """
        Returns the AlignedProteinDnaPosition object of requested ABSOLUTE sbjct protein position

        @attention: sposPY is ABSOLUTE position on sbjct input protein sequence

		@type  sposPY: positive integer
		@param sposPY: protein position on sbjct input protein sequence

		@rtype:  AlignedProteinDnaPosition
        @return: AlignedProteinDnaPosition object
        """
        if self._positions[0].sbjct_pos > sposPY:
            # requested position in front of (current) (extended) alignment 
            return CoordinateOutOfRange
        if self._positions[-1].sbjct_pos < sposPY:
            # requested position after of (current) (extended) alignment 
            return CoordinateOutOfRange
        # loop over the AlignedProteinDnaPosition objects and return
        # it as soon as the coordinates match
        for offset in range(0,len(self._positions)):
            posobj = self._positions[offset]
            if posobj.sbjct_pos == sposPY and posobj.sbjct != "-":
                return posobj
        else:
            # position not in the (current) (extended) alignment
            return CoordinateOutOfRange

    # end of function alignmentobject_by_sbjctaa



    def alignmentpart_by_query(self,start,stop):
        """
        Returns a slice from the protein alignment, defined by query protein coordinates

        @attention: start end stop are ABSOLUTE positions on query input protein sequence

		@type  start: positive integer
		@param start: desired start position on query input protein sequence

		@type  stop:  positive integer
		@param stop:  desired stop position on query input protein sequence

		@rtype:  tuple
        @return: ( query, match, sbjct,(qs,qe,ss,se) )
        """
        positions = []
        for pos in self._positions:
            if pos.query_pos == stop:
                break
            if pos.query_pos >= start:
                positions.append( pos )
        # and make the return data!
        return self._part2returndata(positions)

    # end of function alignmentpart_by_query


    def alignmentpart_by_alignmentposition(self,start,stop):
        """
        coordinate system based on the self._positions !
        start and stop are python list-slice coordinates!

		@rtype:  tuple
        @return: ( query, match, sbjct,(qs,qe,ss,se) )
        """
        return self._part2returndata(self._positions[start:stop])

    # end of function alignmentpart_by_alignmentposition


    def alignmentpart_by_sbjct(self,start,stop):
        """
        Returns a slice from the protein alignment, defined by sbjct protein coordinates

        @attention: start end stop are ABSOLUTE positions on sbjct input protein sequence

		@type  start: positive integer
		@param start: desired start position on sbjct input protein sequence

		@type  stop:  positive integer
		@param stop:  desired stop position on sbjct input protein sequence

		@rtype:  tuple
        @return: ( query, match, sbjct,(qs,qe,ss,se) )
        """
        positions = []
        for pos in self._positions:
            if pos.sbjct_pos == stop:
                break
            if pos.sbjct_pos >= start:
                positions.append( pos )
        # and make the return data!
        return self._part2returndata(positions)

    # end of function alignmentpart_by_sbjct


    def _part2returndata(self,positions):
        """
        """
        q = "".join( [ pos.query for pos in positions ] )
        s = "".join( [ pos.sbjct for pos in positions ] )
        m = "".join( [ pos.match for pos in positions ] )
        try:
            qs = min([ pos.query_pos for pos in positions ] )
            qe = max([ pos.query_pos for pos in positions ] )
            ss = min([ pos.sbjct_pos for pos in positions ] )
            se = max([ pos.sbjct_pos for pos in positions ] )
        except:
            # can happen when it is the truly END of the HSP!
            (qs,qe,ss,se) = (None,None,None,None)
        return (q,m,s,(qs,qe,ss,se))

    # end of function _part2returndata


    ########################################################################
    #### Functions for printing PacbPDNA objects                        ####
    ########################################################################


    def print_protein(self,abs_coords=(0,0),rel_coords=(0,0),_linesize=40):
        """
        Print protein alignment of this PacbPDNA object (query,match,sbjct)

        @type  rel_coords:  tuple of two positive integers
        @param rel_coords:  start and stop position in self._positions to be printed

        @type  _linesize: positive integer
        @param _linesize: line length (in AA coordinates) to be printed

        @attention: argument abs_coords is NOT IMPLEMENTED YET!
        """
        aoa = [] # array of arrays of positions to print
        # if no length range given: print all
        if rel_coords == (0,0): rel_coords = (0,len(self._positions))
        for pos in range(rel_coords[0],rel_coords[1]):
            # check if this position exists!
            try:
                aoa.append(self._positions[pos]._return_print_protein())
            except:
                print "POSITION PRINTING ERROR !??!?!"
                pass
        # and print!
        print self._print_aoa(aoa,_linesize=_linesize)

    # end of function print_protein


    def print_dna(self,abs_coords=(0,0),rel_coords=(0,0),_linesize=40):
        """
        Print DNA alignment of this PacbPDNA object (queryDNA,sbjctDNA)

        @type  rel_coords:  tuple of two positive integers
        @param rel_coords:  start and stop position in self._positions to be printed

        @type  _linesize: positive integer
        @param _linesize: line length (in AA coordinates, so x3) to be printed

        @attention: argument abs_coords is NOT IMPLEMENTED YET!
        """
        aoa = [] # array of arrays of positions to print
        # if no length range given: print all
        if rel_coords == (0,0): rel_coords = (0,len(self._positions))
        for pos in range(rel_coords[0],rel_coords[1]):
            # check if this position exists!
            try:
                aoa.append(self._positions[pos]._return_print_dna())
            except:
                pass
        # and print!
        print self._print_aoa(aoa,_linesize=_linesize)

    # end of function print_dna


    def print_protein_and_dna(self,abs_coords=(0,0),rel_coords=(0,0),_linesize=40):
        """
        Print protein and DNA alignment of this PacbPDNA object (queryDNA,query,match,sbjct,sbjctDNA)

        @type  rel_coords:  tuple of two positive integers
        @param rel_coords:  start and stop position in self._positions to be printed

        @type  _linesize: positive integer
        @param _linesize: line length (in AA coordinates, so x3) to be printed

        @attention: argument abs_coords is NOT IMPLEMENTED YET!
        """
        aoa = [] # array of arrays of positions to print
        # if no length range given: print all
        if rel_coords == (0,0): rel_coords = (0,len(self._positions))
        for pos in range(rel_coords[0],rel_coords[1]):
            # check if this position exists!
            try:
                aoa.append(self._positions[pos]._return_print_protein_and_dna())
            except:
                pass
        # and print!
        print self._print_aoa(aoa,_linesize=_linesize)

    # end of function print_protein_and_dna


    def _print_aoa(self,aoa,_linesize=40):
        """
        """
        rows = []
        if aoa:
            totallength = len(aoa)
            batches = range(0,totallength,_linesize)
            for batch_start in batches:
                batch_end = min([batch_start+_linesize,totallength])
                subrows = [""]*len(aoa[0])
                for item in aoa[batch_start:batch_end]:
                    for i in range(0,len(item)):
                        subrows[i]+=item[i]
                rows.extend(subrows)
                rows.extend("")
        # and return
        return "\n".join(rows)

    # end of function _print_aoa


    ########################################################################
    #### End of functions for printing PacbPDNA objects                 ####
    ########################################################################

    def get_nt_identity(self):
        """ Get identity% of protein-guided DNA alignment """
        qseq,sseq = self.get_unextended_aligned_dna_sequences()
        cnt = 0
        for pos in range(0,len(qseq)):
            if qseq[pos].lower() == sseq[pos].lower():
                cnt+=1
        # return identity ratio
        return float(cnt)/len(qseq)

    # end of function get_nt_identity


    def get_ungapped_nt_identity(self):
        """ Get identity% of protein-guided DNA alignment, gaps ignored """
        qseq,sseq = self.get_unextended_aligned_dna_sequences()
        cnt    = 0
        length = 0
        for pos in range(0,len(qseq)):
            baseQ = qseq[pos].lower()
            baseS = sseq[pos].lower()
            # ignore gaps
            if baseQ == "-" or baseS == "-": continue
            length+=1
            if baseQ == baseS: cnt+=1

        # return identity ratio
        if length == 0:
            # can happen in freak cases; pacbpdna with wrong coordinates or zero length
            return 0.0
        else:
            return float(cnt)/float(length)

    # end of function get_ungapped_nt_identity


    def get_unguided_nt_identity(self):
        """ Get identity% of UNGUIDED DNA alignment """
        # if zerosized -> return 0.0
        if self.length == 0: return 0.0
        # get DNA sequences
        dnaQ,dnaS = self.get_aligned_dna_sequences()
        dnaQ,dnaS = dnaQ.replace("-",""), dnaS.replace("-","")
        # make (semi) unique headers
        uniqueid = get_random_string_tag()
        (qs,qe,ss,se) = self.barcode()[0:4]
        headerQ = "query%s%s%s" % (qs,qe,uniqueid)
        headerS = "sbjct%s%s%s" % (ss,se,uniqueid)
        # prepare & run clustalw
        seqs    = { headerQ: dnaQ, headerS: dnaS }
        out,alignment = clustalw( seqs=seqs )
        # get id% on aligned dna sequences
        cnt = 0
        for pos in range(0,len(out[headerQ])):
            if out[headerQ][pos] == out[headerS][pos]:
                cnt+=1
        # return relative ratio
        return float(cnt) / len(out[headerQ])

    # end of function get_unguided_nt_identity


    def get_aligned_protein_sequences(self):
        """
        get the aligned protein sequences
        return query, match and sbjct track
        """
        q,m,s = [],[],[]
        for posobj in self._positions:
            (_q,_m,_s) = posobj._return_print_protein()
            q.append(_q)
            m.append(_m)
            s.append(_s)
        return ( "".join(q), "".join(m), "".join(s) )

    # end of function get_aligned_protein_sequences


    def get_aligned_dna_sequences(self):
        """
        get the aligned DNA sequences
        return query and sbjct track
        """
        q,s = [],[]
        for posobj in self._positions:
            (_q,_s) = posobj._return_print_dna()
            q.append(_q)
            s.append(_s)
        return ( "".join(q).replace(" ","-"), "".join(s).replace(" ","-") )

    # end of function get_aligned_dna_sequences


    def get_unextended_aligned_protein_sequences(self):
        """ Forwards compatiblity with PacbPORF """
        return self.get_aligned_protein_sequences()

    # end of function get_unextended_aligned_protein_sequences


    def get_unextended_aligned_dna_sequences(self):
        """ Forwards compatiblity with PacbPORF """
        return self.get_aligned_dna_sequences()

    # end of function get_unextended_aligned_dna_sequences


    def _get_gff_group(self):
        """ """
        return "%s-%s" % (self.query_dna_start, self.query_dna_end)

    # end of function _get_gff_group


    def togff(self,gff={},with_alignment=True):
        """
        Return 8-element tuple of gff data.
        To be called from the subclasses, or
        to be overwritten in the subclasses!
        """
        # update self._gff with gff data
        self._gff.update(gff)
        if gff.has_key('column9data'):
            self._gff['column9data'].update( gff['column9data'] )
        if not self._gff.has_key('gname'):   self._gff['gname']   = "%s-%s" % (self.query_dna_start, self.query_dna_end)
        if not self._gff.has_key('fmethod'): self._gff['fmethod'] = self.__class__.__name__

        if with_alignment:
            ( query, match, sbjct,(qs,qe,ss,se) ) = self.nonextended_alignmentpart_by_alignmentposition()
            buffersize = 80
            alignment = ["### Aligned protein sequence coordinates: %s-%s (Q), %s-%s (S)<PRE>" % (qs,qe,ss,se) ]
            for pos in range(0,len(query),buffersize):
                alignment.append( query[pos:pos+buffersize] ) 
                alignment.append( match[pos:pos+buffersize] ) 
                alignment.append( sbjct[pos:pos+buffersize] ) 
                alignment.append( '' ) 
            alignment.append("</PRE>")
            self._gff['column9data']['Alignment'] = "<br>".join(alignment)

        column9string = ""
        # take care for the formatting of column9data
        if self._gff['column9data']:
            tmp = []
            for key,value in self._gff['column9data'].iteritems():
                if str(value).find(" ") > -1: value = "'%s'" % value
                tmp.append( "; %s %s" % (key,value) )
            column9string = "".join(tmp)

        # calculate start & end coordinates
        # method depends on if it is a PacbPDNA or PacbPORF
        if self.__class__.__name__ == 'PacbPORF':
            start = self._positions[self._original_alignment_pos_start].query_dna_start+1
            end   = self._positions[self._original_alignment_pos_end-1].query_dna_end
        else:
            start = self.query_dna_start+1,
            end   = self.query_dna_end,

        # and return the GFF tuple
        return (
            self._gff['fref'],
            self._gff['fsource'],
            self._gff['fmethod'],
            start,
            end,
            "%1.3f" % self.identityscore,
            self._gff['fstrand'],
            self._gff['fphase'],
            "%s %s%s" % ( self._gff['gclass'], self._gff['gname'], column9string )
        )

    # end of function togff


    ########################################################################
    #### Coordinate translation                                         ####
    ########################################################################

    def aapos_query2sbjct(self,pos):
        """ Map (absolute) query protein coordinate to sbjct coordinate """
        algPosCoord = self.alignmentposition_by_query_pos(pos)
        if algPosCoord == CoordinateOutOfRange:
            return CoordinateOutOfRange
        else:
            return self._positions[algPosCoord].sbjct_pos

    # end of function aapos_query2sbjct


    def dnapos_query2sbjct(self,pos):
        """ Map (absolute) query dna coordinate to sbjct coordinate """
        try:
            algPosCoord, phase = self.dnaposition_query(pos)
            # when accidentially hit on a gap -> find first 5p' non-gap
            while self._positions[algPosCoord].isa_gap:
                pos-=3
                algPosCoord, phase = self.dnaposition_query(pos)
            return self._positions[algPosCoord].sbjct_dna_start + phase
        except CoordinateOutOfRange:
            # recognized CoordinateOutOfRange exception
            return CoordinateOutOfRange
        except:
            raise CoordinateOutOfRange

    # end of function dnapos_query2sbjct


    def aapos_sbjct2query(self,pos):
        """ Map (absolute) sbjct protein coordinate to query coordinate """
        algPosCoord = self.alignmentposition_by_sbjct_pos(pos)
        if algPosCoord == CoordinateOutOfRange:
            return CoordinateOutOfRange
        else:
            return self._positions[algPosCoord].query_pos

    # end of function aapos_sbjct2query


    def dnapos_sbjct2query(self,pos):
        """ Map (absolute) sbjct dna coordinate to query coordinate """
        try:
            algPosCoord, phase = self.dnaposition_sbjct(pos)
            # when accidentially hit on a gap -> find first 5p' non-gap
            while self._positions[algPosCoord].isa_gap:
                pos-=3
                algPosCoord, phase = self.dnaposition_sbjct(pos)
            return self._positions[algPosCoord].query_dna_start + phase
        except CoordinateOutOfRange:
            # recognized CoordinateOutOfRange exception
            return CoordinateOutOfRange
        except:
            raise CoordinateOutOfRange

    # end of function dnapos_sbjct2query


    def is_coding(self):
        """ Does this aligment likely represent a coding (DNA) alignment? """
        Qseq,Sseq = self.get_unextended_aligned_dna_sequences()
        ####### ERRONEOUS! why??
        if len(Qseq) != len(Sseq):
            print "\n" , self
            print len(Qseq), "[UNEQUAL LENGTH]" ,len(Sseq)
            print Qseq[0:30], "...", Qseq[-30:]
            print Sseq[0:30], "...", Sseq[-30:]
        ####### ERRONEOUS! why??
        return coding_sequence_periodicity_test(Qseq,Sseq)

    # end of function is_coding

# end of class PacbPDNA

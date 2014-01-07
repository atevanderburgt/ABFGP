"""
"""

# Module metadata
__authors__ = "Ate van der Burgt"
__license__ = "MIT"

# Imports from pacb package
from pacbpdna import PacbPDNA
from pacbp import calculate_identityscore, calculate_identity
from alignedproteindnaposition import *
import conversion
import tcode
from exceptions import * 
from validators import *

# Python Imports
from copy import deepcopy
from sets import Set

# Abgp Imports
import dna2prot
from lib_proteinsimilaritymatrix import ProteinSimilarityMatrix
from pythonlibs.ordering import order_list_by_attribute as olba


# Stopless3n introns imports
from settings.splicesites import KWARGS_PACBPORF_ALIGNED_STOPLESS_3N_INTRONS
from connecting.functions import _update_kwargs
from connecting.orfs import find_stopless3n_introns_on_orf
from lib_stopless3n_intron import _filter_overlapping_stopless3n_introns


DEFAULT_MATRIX = ProteinSimilarityMatrix(name="BLOSUM62")


class PacbPORF(PacbPDNA):
    def __init__(self,pacbp,orfQ,orfS,MATRIX={}):
        """
        """
        # get nucleotide sequence of orfQuery
        orfQdnaseq = orfQ.getntseqslice(
                abs_aa_startPY = pacbp.query_start,
                abs_aa_endPY   = pacbp.query_end )
        # get nucleotide sequence of orfSbjct
        orfSdnaseq = orfS.getntseqslice(
                abs_aa_startPY = pacbp.sbjct_start,
                abs_aa_endPY   = pacbp.sbjct_end )

        if len(orfQdnaseq) % 3 != 0:
            # can happen when we are at the furthest 
            # exteriority of the sequence (Orfs NOT ended by STOP codons
            # solve by completion into triplets with n's
            correction = "NNN"[len(orfSdnaseq) % 3:]
            orfQdnaseq = orfQdnaseq + correction

        if len(orfSdnaseq) % 3 != 0:
            # can happen when we are at the furthest 
            # exteriority of the sequence (Orfs NOT ended by STOP codons
            # solve by completion into triplets with n's
            correction = "NNN"[len(orfSdnaseq) % 3:]
            orfSdnaseq = orfSdnaseq + correction

        # make input_dna tuple
        input_dna = (
            orfQdnaseq,
            orfSdnaseq,
            orfQ.aapos2dnapos(pacbp.query_start),
            orfS.aapos2dnapos(pacbp.sbjct_start)
            )
        # now instantiate the parental PacbPDNA object

        if MATRIX:
            PacbPDNA.__init__(self,pacbp=pacbp, input_dna=input_dna,MATRIX=MATRIX)
        else:
            PacbPDNA.__init__(self,pacbp=pacbp, input_dna=input_dna)

        # set PacbPORF specific attributes
        self.orfQ = orfQ
        self.orfS = orfS
        # set start and end positions of the original alignment
        # this is needed to calculate the `entropy` of a position
        # in the alignment
        self._original_alignment_pos_start = 0
        self._original_alignment_pos_end   = len(self._positions)

    # end of function __init__


    def barcode(self):
        """
        A unique tuple with integers describing this object
        """
        sta = self._get_original_alignment_pos_start()
        end = self._get_original_alignment_pos_end()
        return ( sta.query_pos, end.query_pos, sta.sbjct_pos, end.sbjct_pos, self.bitscore )

    # end of function barcode


    def _swap_query_and_sbjct(self):
        """
        USE WITH CARE!!!
        See PacbPDNA._swap_query_and_sbjct
        """
        # first swap the data from the basal PacbP clas
        PacbPDNA._swap_query_and_sbjct(self)

        # swap the orf objects
        self.orfQ, self.orfS = self.orfS, self.orfQ

    # end of function _swap_query_and_sbjct

    def __str__(self):
        """
        """
        source = ' '
        if self.source != 'blastp':
            source = " src=%s " % self.source
        if not self.is_extended():
            # pacbporf object is NOT extended
            return """<pacbporf%sQ:%s-%s (%s%s) S:%s-%s (%s%s) bits:%s iden:%1.3f len:(%s,%s)>""" % (
                source,
                self.query_start, self.query_end,
                self.orfQ.id, self.orfQ.tcode_symbolic(),
                self.sbjct_start, self.sbjct_end,
                self.orfS.id, self.orfS.tcode_symbolic(),
                self.bitscore, self.identityscore,
                self.query_end - self.query_start,
                self.sbjct_end - self.sbjct_start )
                ###self.query_length, self.sbjct_length )
        else:
            # pacbporf object is EXTENDED
            spos = self._get_original_alignment_pos_start()
            epos = self._get_original_alignment_pos_end()
            return """<pacbporf EXTENDED%sQ:(%s).%s..%s-%s..%s.(%s) (%s%s) S:(%s).%s..%s-%s..%s.(%s) (%s%s) bits:%s iden:%1.3f len:(%s,%s)>""" % (
                source,
                # coordinates of query: orf-start, extended-alignment-start, original-alignment-start
                # the same 3 in reversed order for the end coordinates
                self.orfQ.protein_startPY,
                self.query_start,
                spos.query_pos, epos.query_pos+1,   
                self.query_end,
                self.orfQ.protein_endPY,

                self.orfQ.id, self.orfQ.tcode_symbolic(),

                # coordinates of sbjct: orf-start, extended-alignment-start, original-alignment-start
                # the same 3 in reversed order for the end coordinates
                self.orfS.protein_startPY,
                self.sbjct_start,
                spos.sbjct_pos, epos.sbjct_pos+1,
                self.sbjct_end,
                self.orfS.protein_endPY,

                self.orfS.id, self.orfS.tcode_symbolic(),
                self.bitscore, self.identityscore,
                self.query_end - self.query_start,
                self.sbjct_end - self.sbjct_start )
                ###self.query_length, self.sbjct_length )

    # end of function __str__


    ########################################################################
    #### Functions for obtaining the actually aligned sequence          ####
    ########################################################################

    def get_unextended_aligned_protein_sequences(self):
        """
        get the unextended, acually aligned protein sequences
        return query, match and sbjct track
        """
        return self.get_aligned_protein_sequences_by_algpositions(
            self._original_alignment_pos_start,
            self._original_alignment_pos_end
            )
        #q,m,s = [],[],[]
        #for posobj in self._positions[self._original_alignment_pos_start:self._original_alignment_pos_end]:
        #    (_q,_m,_s) = posobj._return_print_protein()
        #    q.append(_q)
        #    m.append(_m)
        #    s.append(_s)
        #return ( "".join(q), "".join(m), "".join(s) )

    # end of function get_unextended_aligned_protein_sequences


    def get_unextended_aligned_dna_sequences(self):
        """
        get the unextended, acually aligned DNA sequences
        return query and sbjct track
        """
        return self.get_aligned_dna_sequences_by_algpositions(
            self._original_alignment_pos_start,
            self._original_alignment_pos_end
            )
        #q,s = [],[]
        #for posobj in self._positions[self._original_alignment_pos_start:self._original_alignment_pos_end]:
        #    (_q,_s) = posobj._return_print_dna()
        #    q.append(_q)
        #    s.append(_s)
        #return ( "".join(q).replace(" ","-"), "".join(s).replace(" ","-") )

    # end of function get_unextended_aligned_dna_sequences


    def get_aligned_protein_sequences_by_algpositions(self,sta,end):
        """
        @attention: sta & end must be positive integers, end > sta
        """
        q,m,s = [],[],[]
        for posobj in self._positions[sta:end]:
            (_q,_m,_s) = posobj._return_print_protein()
            q.append(_q)
            m.append(_m)
            s.append(_s)
        return ( "".join(q), "".join(m), "".join(s) )

    # end of function get_unextended_aligned_protein_sequences


    def get_aligned_dna_sequences_by_algpositions(self,sta,end):
        """
        @attention: sta & end must be positive integers, end > sta
        """
        q,s = [],[]
        for posobj in self._positions[sta:end]:
            (_q,_s) = posobj._return_print_dna()
            q.append(_q)
            s.append(_s)
        return ( "".join(q).replace(" ","-"), "".join(s).replace(" ","-") )

    # end of function get_aligned_dna_sequences_by_algpositions 

    ########################################################################
    #### Functions concerning TCODE data of Orf objects                 ####
    ########################################################################

    def tcode(self):
        """@attention: see pacb.tcode.tcode() for documentation """
        return tcode.tcode(self)
    
    # end of function tcode
  
 
    def tcode_query(self):
        """@attention: see pacb.tcode.tcode_query() for documentation """
        return tcode.tcode_query(self)
    
    # end of function tcode_query
    
    
    def tcode_sbjct(self):
        """@attention: see pacb.tcode.tcode_sbjct() for documentation """
        return tcode.tcode_sbjct(self)
    
    # end of function tcode_sbjct

    ############################################################################
    #### Function to extend PacbPORF matching tails (after PacbPORF trimming)
    ############################################################################

    def extend_matched_ends(self):
        """
        """
        # TODO: this function might be called with or without that
        # this PacbPORF is extended. But, for now, we only need this functionalty
        # on extended PacbPORF -> omit this extra coding effort
        if not self.is_extended(): return False

        spos = self._get_original_alignment_pos_start()
        epos = self._get_original_alignment_pos_end()

        IS_EXTENDED = False
        while spos.position > 0:
            prevpos = self._positions[spos.position-1]
            if "*" in [prevpos.query,prevpos.sbjct]: break
            if "-" in [prevpos.query,prevpos.sbjct]: break
            bitscore= self.MATRIX.scorealignment(prevpos.query,prevpos.sbjct)
            if bitscore <= 0: break
            IS_EXTENDED = True
            if prevpos.query.upper() == prevpos.sbjct.upper():
                prevpos.match = prevpos.query.upper()
                prevpos.isa = 'identity' 
            else:
                prevpos.match = '+'
                prevpos.isa = 'similarity' 

            # decrease self._original_alignment_pos_start position
            self._original_alignment_pos_start -=1

            # re-get current original_alignment_pos_start object
            spos = self._get_original_alignment_pos_start()

        while epos.position < len(self._positions)-1:
            nextpos = self._positions[epos.position+1]
            if "*" in [nextpos.query,nextpos.sbjct]: break
            if "-" in [nextpos.query,nextpos.sbjct]: break
            bitscore= self.MATRIX.scorealignment(nextpos.query,nextpos.sbjct)
            if bitscore <= 0: break
            IS_EXTENDED = True
            if nextpos.query.upper() == nextpos.sbjct.upper():
                nextpos.match = nextpos.query.upper()
                nextpos.isa = 'identity' 
            else:
                nextpos.match = '+'
                nextpos.isa = 'similarity' 

            # increase self._original_alignment_pos_end position
            self._original_alignment_pos_end +=1

            # re-get current original_alignment_pos_end object
            epos = self._get_original_alignment_pos_end()

        # prepare return status. If IS_EXTENDED, rescore alignment attributes
        if IS_EXTENDED:
            self._score_alignment()
            return True
        else:
            return False

    # end of function extend_matched_ends


    def is_extended(self):
        """
        Is PacbPORF object extended or not?
        """
        if self._original_alignment_pos_start == 0 and\
        self._original_alignment_pos_end == len(self._positions):
            return False
        else:
            return True

    # end of function is_extended


    def alignment_protein_range_query(self):
        """
        OVERRIDES PacbP.alignment_protein_range_query()
        """
        if not self.is_extended():
            return PacbPDNA.alignment_protein_range_query(self)
        else:
            spos = self._get_original_alignment_pos_start()
            epos = self._get_original_alignment_pos_end()
            # correct for +1 beacause coord is a python list slice coord!
            return range( spos.query_pos, epos.query_pos+1 )

    # end of function alignment_protein_range_query


    def alignment_protein_range_sbjct(self):
        """
        OVERRIDES PacbP.alignment_protein_range_sbjct()
        """
        if not self.is_extended():
            return PacbPDNA.alignment_protein_range_sbjct(self)
        else:
            spos = self._get_original_alignment_pos_start()
            epos = self._get_original_alignment_pos_end()
            # correct for +1 beacause coord is a python list slice coord!
            return range( spos.sbjct_pos, epos.sbjct_pos+1 )

    # end of function alignment_protein_range_query


    def alignment_dna_range_query(self):
        """
        OVERRIDES PacbPDNA.alignment_dna_range_query()
        """
        if not self.is_extended():
            return PacbPDNA.alignment_dna_range_query(self)
        else:
            spos = self._get_original_alignment_pos_start()
            epos = self._get_original_alignment_pos_end()
            # correct for +1 beacause coord is a python list slice coord!
            return range( spos.query_dna_start, epos.query_dna_end+1 )

    # end of function alignment_dna_range_query


    def alignment_dna_range_sbjct(self):
        """
        OVERRIDES PacbPDNA.alignment_dna_range_sbjct()
        """
        if not self.is_extended():
            return PacbPDNA.alignment_dna_range_sbjct(self)
        else:
            spos = self._get_original_alignment_pos_start()
            epos = self._get_original_alignment_pos_end()
            # correct for +1 beacause coord is a python list slice coord!
            return range( spos.sbjct_dna_start, epos.sbjct_dna_end+1 )

    # end of function alignment_dna_range_sbjct


    def alignment_entropy(self,position,method=None,window=10):
        """
        http://en.wikipedia.org/wiki/Faulhaber%27s_formula

        Entropy is calculated with the first order
        Faulhaber's formula ( F(n) = ( n^2 + n ) / 2
        over a range of size `window` in both side of
        the cutoff `position` on a specific side `method`

        See the example below; the alignment is converted
        from non-aligned (n) to aligned (a) positions
        into a binary vector.  A, B, C, D, E represent
        `positions` that are checked, in this case for the
        `method` right/donor.
        A window of size `window` (in the example 10) is taken
        on both sides on the `position`, empty positions (-)
        are treated as zeros.
        The sum of both `window`-vector is taked as `n` in the 
        first order Faulhaber's formula.
        

        nnnnnnnnnAnaaaaaaaaaaaaaaBaaaaaaCaaaaDnnnEnnnnn
        00000000000111111111111111111111111111000000000
        A000000000          
        A         0111111111
        A    0    -    45                                     = -45 => 0
        B               1111111111
        B                         1111111111
        B                   55    -   55                      = 0
        C                      1111111111
        C                                1111100000
        C                          55    -   15               = 40
        D                           1111111111
        D                                     000000000-
        D                               55    -   0           = 55
        E                               1111110000
        E                                         00000-----
        E                                   21    -   0       = 21


        For the given `window` of 10, this will result in
        the following pattern around the maximum (ratio's below it):

         0 10 19 27 34 40 45 49 52 54 55 45 36 28 21 15 10  6  3  1  0 
        .0 .2 .3 .5 .6 .7 .8 .9 .9 1. 1. .8 .7 .5 .4 .3 .2 .1 .1 .0 .0 

        Note the relative higher scores directly left from the
        entropy peak (55) compared to directly rigth of it.
        This correstponds nicely to the biological phenomenon that
        a few nt's downstream of a splice donor can be still alignable.
        On the contrary, the chance that the splice site itself and the
        amino acid directly upstream of it are not alignable anymore.
        """

        methods = ['left','acceptor','donor','right']
        if method not in methods:
            message = 'method must be in %s' % str(methods)
            raise InproperlyAppliedArgument, message
        if type(position) == type(int()):
            int_pos = position
        else:
            try:
                int_pos = position.position
            except:
                message = "'position' is not an integer or an AlignedProteinDnaPosition"
                raise InproperlyAppliedArgument, message
        # get the binary entropy array
        bea = [0]*len(self._positions)
        bea[self._original_alignment_pos_start:self._original_alignment_pos_end] = [1]*(self._original_alignment_pos_end-self._original_alignment_pos_start)
        # calculate maximum entropy score
        max_entropy = float(sum(range(1,window+1)))

        if method in ['donor','right']:
            l1,l2 = max([ (int_pos+1-window) , 0 ]), max([ (int_pos+1) , 0 ])
            r1,r2 = max([ (int_pos+1) , 0 ]), max([ (int_pos+1+window) , 0 ])
            score_left  = sum(range(1,sum(bea[l1:l2])+1))
            score_rigth = sum(range(1,sum(bea[r1:r2])+1))
            entropy = score_left - score_rigth
            entropy = float(entropy) / max_entropy
            ### CORRECTION AT 11/03/2009: allow for negative entropy!
            ###if entropy < 0.0: entropy = 0.0
            ### CORRECTION AT 11/03/2009: allow for negative entropy!
            return entropy
        if method in ['left','acceptor']:
            l1,l2 = max([ (int_pos-window) , 0 ]), max([ int_pos , 0 ])
            r1,r2 = max([ int_pos , 0 ]), max([ (int_pos+window) , 0 ])
            score_left  = sum(range(1,sum(bea[l1:l2])+1))
            score_rigth = sum(range(1,sum(bea[r1:r2])+1))
            entropy = score_rigth - score_left
            entropy = float(entropy) / max_entropy
            ### CORRECTION AT 11/03/2009: allow for negative entropy!
            ###if entropy < 0.0: entropy = 0.0
            ### CORRECTION AT 11/03/2009: allow for negative entropy!
            return entropy

    # end of function alignment_entropy

    #######################################################################

    def _append_one_left_position(self,_force_append=False):
        """
        Append one non-aligned position on the left (5') side
        of the alignment. Data (protein,dna) is taken from the Orf objects.

        @attention: USE WITH CARE!

        @rtype:  Boolean
        @return: True or False, depending on if a position is appended
        """

        # get to-be-appended protein and dna strings
        newQaa  = self.orfQ.getaa( abs_pos = self.query_start-1 )
        newMsym = " "
        newSaa  = self.orfS.getaa( abs_pos = self.sbjct_start-1 )
        newQdna = self.orfQ.getntseq( abs_aa_pos = self.query_start-1 )
        newSdna = self.orfS.getntseq( abs_aa_pos = self.sbjct_start-1 )

        # check if all strings are non-empty
        if (newQaa and newSaa and newQdna and newSdna) or _force_append:

            if _force_append or (len(newQdna)<3 or len(newSdna)<3):
                ## Q and/or S amino-acids are out-of-range of the
                ## original Orf. Get them from the nucleotidesequence
                #newQaa = dna2prot.dna2protein(newQdna)
                #newSaa = dna2prot.dna2protein(newSdna)
                # Q and/or S amino-acids are out-of-range of the
                # original Orf. Get them from the nucleotidesequence
                if len(newQdna) < 3:
                    # can happen when we are at the furthest 
                    # exteriority of the sequence -> complete triplet with n's
                    newQdna = newQdna + "nnn"
                    newQdna = newQdna[0:3]
                if len(newSdna) < 3:
                    # can happen when we are at the furthest 
                    # exteriority of the sequence -> complete triplet with n's
                    newSdna = newSdna + "nnn"
                    newSdna = newSdna[0:3]
                # translate in amino acid sequence
                newQaa = dna2prot.dna2protein(newQdna)
                newSaa = dna2prot.dna2protein(newSdna)


            # yep, we can add a position!
            newpos  = 0
            newAPDP = AlignedProteinDnaPosition(
                        newpos,self.query_start-1,self.sbjct_start-1,
                        alignment=(newQaa,newMsym,newSaa),
                        query_dna=(newQdna,self.query_dna_start-3),
                        sbjct_dna=(newSdna,self.sbjct_dna_start-3)
                        )

            # increase `position` counter for all existing positions
            for pos in self._positions: pos.position+=1
                
            # append new position by inserting on first position
            self._positions.insert(0,newAPDP)

            # increase the original alignment start position
            self._original_alignment_pos_start += 1
            self._original_alignment_pos_end += 1

            # and decrease all the absolute positions in the object
            self.query_start-=1
            self.sbjct_start-=1
            self.query_dna_start-=3
            self.sbjct_dna_start-=3

            # increase all lengths
            self.length+=1
            self.query_length+=1
            self.sbjct_length+=1
            self.query_dna_length+=3
            self.sbjct_dna_length+=3

            # append all sequence characters to the strings
            self.query      = newQaa  + self.query
            self.match      = newMsym + self.match
            self.sbjct      = newSaa  + self.sbjct
            self.alignment  = " " + self.alignment
            self.query_dna  = newQdna + self.query_dna
            self.sbjct_dna  = newSdna + self.sbjct_dna
            self.query_protein = newQaa + self.query_protein
            self.sbjct_protein = newSaa + self.sbjct_protein

            # and (re) score the aligment
            ### self._score_alignment()

            # return succesful status
            return True
        else:
            # no, something went wrong. End of Orf??
            return False

    # end of function _append_one_left_position


    def pop_alignment_object(self,aas=1):
        """
        Remove the final AlignedProteinDnaPosition object

        @type  aas: (small) positive integer
        @param aas: number of AAs to extend the alignment to the rigth/3p side

        @attention: USE WITH CARE!!!
        @attention: only possible on UNEXTENDED PacbPORF objects

        @rtype:  Boolean
        @return: True or False, depending on if position(s) is/are popped 
        """
        if self.is_extended():
            message = "pop_alignment_object() function only available on unextended PacbPORF objects"
            raise InproperlyAppliedArgument, message

        POPPED_CNT = 0 
        for iter in range(0,aas):
            # if an position is popped that is part of the actual
            # original alignment, reset the original alignment pointer
            if len(self._positions) == self._original_alignment_pos_end:
                self._original_alignment_pos_end -= 1
            # remove the AlignedProteinDnaPosition object
            try:
                self._positions.pop()
                POPPED_CNT+=1
            except IndexError:
                # nothing left to pop! Restore end coord
                self._original_alignment_pos_end += 1
                POPPED_CNT-=1
                break

            # get the new final AlignedProteinDnaPosition object
            endPos = self._positions[-1]

            # reset end coordinates
            self.query_end        = endPos.query_pos + 1
            self.sbjct_end        = endPos.sbjct_pos + 1

            # dna coordinate ends can be 0 in case a gap occurs on this exact popped position
            for i in range(-1,-len(self._positions),-1):
                if self._positions[i].query_dna_end:
                    self.query_dna_end    = self._positions[i].query_dna_end # end coord isa listslice coord!
                    break
            for i in range(-1,-len(self._positions),-1):
                if self._positions[i].sbjct_dna_end:
                    self.sbjct_dna_end    = self._positions[i].sbjct_dna_end # end coord isa listslice coord!
                    break

            # recreate all sequence strings by slicing
            # self.alignment is recreated by _score_alignment()
            self.query            = self.query[0:-1]
            self.match            = self.match[0:-1]
            self.sbjct            = self.sbjct[0:-1]
            self.query_dna        = self.query_dna[0:-3]
            self.sbjct_dna        = self.sbjct_dna[0:-3]
            self.query_protein    = self.query.replace("-","")
            self.sbjct_protein    = self.sbjct.replace("-","")
   
            # reset length attributes
            self._reset_length_attributes()
    
        # rescore the aligment
        self._score_alignment()

        # return status for if positions are popped 
        return POPPED_CNT > 0

    # end of function pop_alignment_object


    def append_alignment_object(self,algobjs):
        """
        Append AlignedProteinDnaPosition(s) to the alignment

        @attention: USE WITH CARE!
        @attention: only possible on UNEXTENDED PacbPORF objects

        @type  algobjs: AlignedProteinDnaPosition or list of AlignedProteinDnaPositions
        @param algobjs: AlignedProteinDnaPosition(s) to append to the alignment

        @rtype:  Boolean
        @return: True or False, depending on if position(s) is/are appended
        """
        if algobjs.__class__.__name__ == 'list':
            # assume a list of AlignedProteinDnaPosition objects
            pass
        elif algobjs.__class__.__name__ == 'AlignedProteinDnaPosition':
            algobjs = [ algobjs ]
        else:
            message = "argument isa '%s', not a (list of) AlignedProteinDnaPosition(s)" % algobjs.__class__.__name__
            raise InproperlyAppliedArgument, message

        if self.is_extended():
            message = "append_alignment_object() function only available on unextended PacbPORF objects"
            raise InproperlyAppliedArgument, message

        while algobjs:
            obj = algobjs.pop(0)
            self._positions.append(obj)
            # reset the original alignment pointer
            self._original_alignment_pos_end += 1

            # get the new final AlignedProteinDnaPosition object
            endPos  = self._positions[-1]
            prevPos = self._positions[-2]

            # reset the position attribute in this new AlignedProteinDnaPosition
            endPos.position = self._positions[-2].position + 1

            # reset all query coordinate attributes
            if endPos.query != "-":
                endPos.query_pos       = prevPos.query_pos + 1
                endPos.query_dna_start = prevPos.query_dna_end
                endPos.query_dna_end   = endPos.query_dna_start + 3 
            else: 
                endPos.query_pos       = prevPos.query_pos
                endPos.query_dna_start = 0
                endPos.query_dna_end   = 0

            # reset all sbjct coordinate attributes
            if endPos.sbjct != "-":
                endPos.sbjct_pos       = prevPos.sbjct_pos + 1
                endPos.sbjct_dna_start = prevPos.sbjct_dna_end
                endPos.sbjct_dna_end   = endPos.sbjct_dna_start + 3
            else:
                endPos.sbjct_pos       = prevPos.sbjct_pos 
                endPos.sbjct_dna_start = 0
                endPos.sbjct_dna_end   = 0

            # reset end coordinates
            self.query_end        = endPos.query_pos + 1
            self.sbjct_end        = endPos.sbjct_pos + 1
            self.query_dna_end    = endPos.query_dna_end  # end coord isa listslice coord! 
            self.sbjct_dna_end    = endPos.sbjct_dna_end  # end coord isa listslice coord! 

            # recreate all sequence strings by appending
            # self.alignment is recreated by _score_alignment()
            self.query            = self.query + endPos.query
            self.match            = self.match + endPos.match
            self.sbjct            = self.sbjct + endPos.sbjct
            self.query_dna        = self.query_dna + endPos.query_dna_seq
            self.sbjct_dna        = self.sbjct_dna + endPos.sbjct_dna_seq
            self.query_protein    = self.query.replace("-","")
            self.sbjct_protein    = self.sbjct.replace("-","")
    
            # reset length attributes
            self._reset_length_attributes()
    
        # rescore the aligment
        self._score_alignment()

        # return status for if positions are appended 
        return len(algobjs) > 0

    # end of function append_alignment_object


    def _reset_length_attributes(self):
        """
        Recalculate length attributes in object; use after an alignments update!
        """
        self.length           = len(self.query)
        self.query_length     = self.query_end - self.query_start
        self.sbjct_length     = self.sbjct_end - self.sbjct_start
        self.query_dna_length = len(self.query_dna)
        self.sbjct_dna_length = len(self.sbjct_dna)

    # end of function _reset_length_attributes
 


    def shift_alignment_object(self,aas=1):
        """
        Shift AlignedProteinDnaPosition(s) from the start of the alignment

        @attention: USE WITH CARE!
        @attention: only possible on UNEXTENDED PacbPORF objects

        @type  aas: (small) positive integer
        @param aas: number of AAs to shift from the left/5p side of the alignment

        @rtype:  Boolean
        @return: True or False, depending on if position(s) is/are shifted
        """
        if self.is_extended():
            message = "shift_alignment_object() function only available on unextended PacbPORF objects"
            raise InproperlyAppliedArgument, message

        shifted = 0
        for iter in range(0,aas):
            # remove the first AlignedProteinDnaPosition object
            try:
                removed = self._positions.pop(0)
                self._original_alignment_pos_end -= 1
                shifted+=1
            except IndexError:
                # nothing left to pop!
                break

            # decrease `position` counter for all existing positions
            for pos in self._positions: pos.position-=1

            # get the new frontal AlignedProteinDnaPosition object
            staPos = self._positions[0]

            # reset start coordinates
            self.query_start        = staPos.query_pos
            self.sbjct_start        = staPos.sbjct_pos
            self.query_dna_start    = staPos.query_dna_start
            self.sbjct_dna_start    = staPos.sbjct_dna_start

            # recreate all sequence strings by inserting
            # self.alignment is recreated by _score_alignment()
            self.query            = self.query[1:]
            self.match            = self.match[1:]
            self.sbjct            = self.sbjct[1:]
            self.query_dna        = self.query_dna[3:]
            self.sbjct_dna        = self.sbjct_dna[3:]
            self.query_protein    = self.query.replace("-","")
            self.sbjct_protein    = self.sbjct.replace("-","")

            # reset length attributes
            self._reset_length_attributes()

        # rescore the aligment
        self._score_alignment()

        # return status for if positions are shifted
        return shifted > 0

    # end of shift_alignment_object


    def unshift_alignment_object(self,algobjs):
        """
        Unshift AlignedProteinDnaPosition(s) to the alignment

        @attention: USE WITH CARE!
        @attention: only possible on UNEXTENDED PacbPORF objects

        @type  algobjs: AlignedProteinDnaPosition or list of AlignedProteinDnaPositions
        @param algobjs: AlignedProteinDnaPosition(s) to unshift to the alignment
        
        @rtype:  Boolean
        @return: True or False, depending on if position(s) is/are unshifted
        """
        if algobjs.__class__.__name__ == 'list':
            # assume a list of AlignedProteinDnaPosition objects
            pass
        elif algobjs.__class__.__name__ == 'AlignedProteinDnaPosition':
            algobjs = [ algobjs ]
        else:
            message = "argument isa '%s', not a (list of) AlignedProteinDnaPosition(s)" % algobjs.__class__.__name__
            raise InproperlyAppliedArgument, message

        if self.is_extended():
            message = "unshift_alignment_object() function only available on unextended PacbPORF objects"
            raise InproperlyAppliedArgument, message

        while algobjs:
            obj = algobjs.pop()

            # increase `position` counter for all existing positions
            for pos in self._positions: pos.position+=1
            # insert new object first in row
            self._positions.insert(0,obj)
            # set its position attribute to zero
            self._positions[0].position = 0
            # original alignment, reset the original alignment pointer
            # at the end of the alignment (it is increaded in length)
            self._original_alignment_pos_end += 1

            # get the new frontal AlignedProteinDnaPosition object
            staPos = self._positions[0]

            # reset start coordinates
            self.query_start        = staPos.query_pos
            self.sbjct_start        = staPos.sbjct_pos
            self.query_dna_start    = staPos.query_dna_start
            self.sbjct_dna_start    = staPos.sbjct_dna_start

            # recreate all sequence strings by inserting
            # self.alignment is recreated by _score_alignment()
            self.query            = staPos.query + self.query
            self.match            = staPos.match + self.match
            self.sbjct            = staPos.sbjct + self.sbjct
            self.query_dna        = staPos.query_dna_seq + self.query_dna
            self.sbjct_dna        = staPos.sbjct_dna_seq + self.sbjct_dna
            self.query_protein    = self.query.replace("-","")
            self.sbjct_protein    = self.sbjct.replace("-","")

            # reset length attributes
            self._reset_length_attributes()

        # rescore the aligment
        self._score_alignment()

        # return status for if positions are unshifted
        return len(algobjs) > 0

    # end of function unshift_alignment_object


    def _append_one_right_position(self,_force_append=False):
        """
        Append one non-aligned position on the right (3') side
        of the alignment. Data (protein,dna) is taken from the Orf objects.

        @attention: USE WITH CARE!

        @rtype:  Boolean
        @return: True or False, depending on if a position is appended
        """

        # get to-be-appended protein and dna strings
        newQaa  = self.orfQ.getaa( abs_pos = self.query_end )
        newMsym = " "
        newSaa  = self.orfS.getaa( abs_pos = self.sbjct_end )
        newQdna = self.orfQ.getntseq( abs_aa_pos = self.query_end )
        newSdna = self.orfS.getntseq( abs_aa_pos = self.sbjct_end )

        # check if all strings are non-empty
        if (newQaa and newSaa and newQdna and newSdna) or _force_append:

            if _force_append or (len(newQdna)<3 or len(newSdna)<3):
                ## Q and/or S amino-acids are out-of-range of the
                ## original Orf. Get them from the nucleotidesequence
                #newQaa = dna2prot.dna2protein(newQdna)
                #newSaa = dna2prot.dna2protein(newSdna)
                # Q and/or S amino-acids are out-of-range of the
                # original Orf. Get them from the nucleotidesequence
                if len(newQdna) < 3:
                    # can happen when we are at the furthest 
                    # exteriority of the sequence -> complete triplet with n's
                    newQdna = newQdna + "nnn"
                    newQdna = newQdna[0:3]
                if len(newSdna) < 3:
                    # can happen when we are at the furthest 
                    # exteriority of the sequence -> complete triplet with n's
                    newSdna = newSdna + "nnn"
                    newSdna = newSdna[0:3]
                # translate in amino acid sequence
                newQaa = dna2prot.dna2protein(newQdna)
                newSaa = dna2prot.dna2protein(newSdna)


            # yep, we can add a position!
            newpos  = len(self._positions)
            newAPDP = AlignedProteinDnaPosition(
                        newpos,self.query_end,self.sbjct_end,
                        alignment=(newQaa,newMsym,newSaa),
                        query_dna=(newQdna,self.query_dna_end),
                        sbjct_dna=(newSdna,self.sbjct_dna_end)
                        )

            # append new position by append in queu
            self._positions.append(newAPDP)

            # and increase all the absolute positions in the object
            self.query_end+=1
            self.sbjct_end+=1
            self.query_dna_end+=3
            self.sbjct_dna_end+=3

            # increase all lengths
            self.length+=1
            self.query_length+=1
            self.sbjct_length+=1
            self.query_dna_length+=3
            self.sbjct_dna_length+=3

            # append all sequence characters to the strings
            self.query      = self.query + newQaa
            self.match      = self.match + newMsym
            self.sbjct      = self.sbjct + newSaa 
            self.alignment  = self.alignment + " "
            self.query_dna  = self.query_dna + newQdna
            self.sbjct_dna  = self.sbjct_dna + newSdna
            self.query_protein = self.query_protein + newQaa
            self.sbjct_protein = self.sbjct_protein + newSaa

            # and (re) score the aligment
            ### self._score_alignment()

            # return succesful status
            return True
        else:
            # no, something went wrong. End of Orf??
            return False

    # end of function _append_one_rigth_position


    def _score_alignment(self):
        """
        Score several properties of the alignment.
        After a manual change of a Pacb object (splitting or extension),
        do not forget to re-call this function!

        OVERRIDES PacbP._score_alignment() because PacbPORF can be an extended

        """
        # (re) set the clustalW-like alignment string
        self._match2alignment()
        # take only the ORIGINAL alignment part to do the scoring upon
        origalignment = self.alignment[self._original_alignment_pos_start:self._original_alignment_pos_end]
        origlength    = len(origalignment)
        origquery     = self.query[self._original_alignment_pos_start:self._original_alignment_pos_end]
        origsbjct     = self.sbjct[self._original_alignment_pos_start:self._original_alignment_pos_end]
        # do some counts on this string
        self.identity       = origalignment.count("*")
        self.similarity     = origlength - origalignment.count(" ") -\
                              self.identity
        # calculate identity% and identityscore%, the latter counting similarities for 0.5
        self.identityscore  = calculate_identityscore(origalignment)
        self.identityperc   = calculate_identity(origalignment)
        # calculate bitscore/bits with MATRIX object
        self.bitscore       = self.MATRIX.scorealignment(origquery,origsbjct)
        self.bits           = self.bitscore

    # end of function _score_alignment


    def nonextended_alignmentpart_by_alignmentposition(self):
        """
        coordinate system based on the self._positions !
        start and stop are python list-slice coordinates!

		@rtype:  tuple
        @return: ( query, match, sbjct,(qs,qe,ss,se) )
        """
        start = self._original_alignment_pos_start
        stop  = self._original_alignment_pos_end
        return self._part2returndata(self._positions[start:stop])

    # end of function nonextended_alignmentpart_by_alignmentposition


    def extend_pacbporf_untill_stops(self):
        """
        """
        # extend on the left/5' side untill stop is reached
        status = True
        while status and self.query_start > self.orfQ.protein_startPY-1 and \
                self.sbjct_start > self.orfS.protein_startPY-1:
            status = self._append_one_left_position()
        # extend on the rigth/3' side untill stop is reached
        status = True
        #while status and "*" not in ( self.query[-1], self.sbjct[-1] ):
        while status and self.query_end <= self.orfQ.protein_endPY and \
                self.sbjct_end <= self.orfS.protein_endPY:
            status = self._append_one_right_position()

    # end of function extend_pacbporf_untill_stops


    def extend_pacbporf_after_stops(self,left_side=2,rigth_side=1):
        """
        left_side   2 positions because (longer) PSSM-acceptor must be fully covered
        right_side  1 position because (shorter) PSSM-donor must be fully covered
        """

        # check if extend_pacbporf_untill_stops is already performed
        check1 = self.query_start == self.orfQ.protein_startPY-1 or self.sbjct_start == self.orfS.protein_startPY-1
        check2 = self.query_end == self.orfQ.protein_endPY+1 or self.sbjct_end == self.orfS.protein_endPY+1

        if (check1,check2) != (True,True):
            # first extend UNTILL stops
            self.extend_pacbporf_untill_stops()

        # check if extend_pacbporf_after_stops is already performed
        check3 = self.query_start == self.orfQ.protein_startPY-1-left_side or self.sbjct_start == self.orfS.protein_startPY-1-left_side
        check4 = self.query_end == self.orfQ.protein_endPY+1+rigth_side or self.sbjct_end == self.orfS.protein_endPY+1+rigth_side

        while not check3:
            # and force the appending of one extra position on the left
            self._append_one_left_position(_force_append=True)
            # redo the check
            check3 = self.query_start == self.orfQ.protein_startPY-1-left_side or self.sbjct_start == self.orfS.protein_startPY-1-left_side

        while not check4:
            # and force the appending of one extra position on the rigth
            self._append_one_right_position(_force_append=True)
            # redo the check
            check4 = self.query_end == self.orfQ.protein_endPY+1+rigth_side or self.sbjct_end == self.orfS.protein_endPY+1+rigth_side

    # end of function extend_pacbporf_after_stops


    def unextend_pacbporf(self,left_side=True,rigth_side=True):
        """
        Unextend an extended PacbPORF object back to its original state
        """
        start  = self._original_alignment_pos_start
        end    = self._original_alignment_pos_end
        staPos = self._positions[start]
        endPos = self._positions[end-1]

        # reset start coordinates
        self.query_start      = staPos.query_pos
        self.sbjct_start      = staPos.sbjct_pos
        self.query_dna_start  = staPos.query_dna_start
        self.sbjct_dna_start  = staPos.sbjct_dna_start

        # reset end coordinates
        self.query_end        = endPos.query_pos + 1
        self.sbjct_end        = endPos.sbjct_pos + 1
        self.query_dna_end    = endPos.query_dna_end # dna coords are listslice coords! 
        self.sbjct_dna_end    = endPos.sbjct_dna_end # dna coords are listslice coords!

        # recreate all sequence strings by slicing
        self.query            = self.query[start:end]
        self.match            = self.match[start:end]
        self.sbjct            = self.sbjct[start:end]
        self.alignment        = self.alignment[start:end]
        self.query_dna        = self.query_dna[start*3:end*3]
        self.sbjct_dna        = self.sbjct_dna[start*3:end*3]
        self.query_protein    = self.query.replace("-","")
        self.sbjct_protein    = self.sbjct.replace("-","")

        # reset length attributes
        self.length           = len(self.query)
        self.query_length     = self.query_end - self.query_start
        self.sbjct_length     = self.sbjct_end - self.sbjct_start
        self.query_dna_length = len(self.query_dna)
        self.sbjct_dna_length = len(self.sbjct_dna)

        # reset _positions list
        self._positions = self._positions[start:end]
        self._original_alignment_pos_start = 0
        self._original_alignment_pos_end = len(self._positions)

        # reset position attributes in AlignedProteinDnaPosition list _positions
        for pos in range(0,len(self._positions)):
            self._positions[pos].position = pos

        # rescore the aligment
        self._score_alignment()

    # end of function unextend_pacbporf


    def detect_aligned_stopless3n_intron_in_query(self,verbose=False,**kwargs):
        """ """
        # edit **kwargs dictionary for some forced attributes
        _update_kwargs(kwargs,KWARGS_PACBPORF_ALIGNED_STOPLESS_3N_INTRONS)

        mind = self._get_original_alignment_pos_start().query_dna_start
        maxd = self._get_original_alignment_pos_end().query_dna_end
        mina = self._get_original_alignment_pos_start().query_dna_start
        maxa = self._get_original_alignment_pos_end().query_dna_end

        if maxd-mind < 75:
            # do not search for a stopless3n intron in a query of only 75nt
            # It is not likely that there is space for a stopless3n intron
            # in this alignment. It would represent two tiny exons separated
            # by an stopless3n intron (so in the same Orf!). Yeah right ;-)
            return []

        stopless3nsQ = find_stopless3n_introns_on_orf(self.orfQ,
                min_donor_pos = mind, max_donor_pos = maxd,
                min_acceptor_pos= mina, max_acceptor_pos= maxa,
                **kwargs)
        stopless3nsQ = _filter_overlapping_stopless3n_introns(stopless3nsQ)

        accepted_stopless3nsQ = []
        at_ratio_pacbporf =  self.get_at_ratio_query()
        for intron in stopless3nsQ:
            # re-check it's position towards the PacbPORF
            if intron.start <= mind: continue
            if intron.end   <= mind: continue
            if intron.start >= maxd: continue
            if intron.end   >= maxd: continue

            if kwargs['has_at_increase']:
                at_ratio_intron = intron.get_at_ratio()
                if at_ratio_pacbporf > at_ratio_intron:
                    if verbose: print "OMITTED AT% :", intron
                    continue


            # if here, get pacbp slice object of this putative intron
            pacbp = conversion.pacbporf2pacbp(self)

            try:

                slice = pacbp.returnslice(
                        self.orfQ.dnapos2aapos(intron.start),
                        self.orfQ.dnapos2aapos(intron.end),
                        coords_on='query'
                        )

            except:

                print self
                print pacbp
                print intron, intron.coords()
                print self.orfQ.dnapos2aapos(intron.start)
                print self.orfQ.dnapos2aapos(intron.end)

                slice = pacbp.returnslice(
                        self.orfQ.dnapos2aapos(intron.start),
                        self.orfQ.dnapos2aapos(intron.end),
                        coords_on='query'
                        )


            # calculate identityscore of the surrounding PacbP(ORF) parts
            # and ratios of identityscores.
            # ====================================== self.identityscore
            # ========================         ===== identity_score_surrounding
            #                         =========      slice.identityscore
            identity_score_surrounding = ( float(self.identity-slice.identity)+\
                    0.5*( self.similarity - slice.similarity ) ) /\
                    ( self.get_unextended_length() - slice.length )
            # calculate ids ratios; max([0.01,...]) to avoid ZeroDivisionError!
            ids_ratio_1 = identity_score_surrounding / self.identityscore
            ids_ratio_2 = self.identityscore / max([0.01,slice.identityscore])

            is_accepted = False
            if slice.bitscore <= 0:
                # append to accepted stopless3n intron list
                accepted_stopless3nsQ.append( intron )
                is_accepted = True
            elif ids_ratio_1 >= 1.10 and ids_ratio_2 >= 1.30:
                # append to accepted stopless3n intron list
                accepted_stopless3nsQ.append( intron )
                is_accepted = True
            else:
                pass

            ####################################################################
            if verbose and is_accepted:
                print "Q", self, self.get_at_ratio_query()
                print "Q", intron, intron.get_at_ratio()
                print slice
                sliceorf = conversion.pacbp2pacbporf(slice,self.orfQ,self.orfS)
                sliceorf._append_one_left_position()
                sliceorf._append_one_left_position()
                sliceorf._append_one_right_position()
                sliceorf._append_one_right_position()
                print "Q", sliceorf, "(%1.2f-%1.2f)" % (ids_ratio_1,ids_ratio_2)
                sliceorf.print_protein_and_dna()
            elif verbose:
                print "OMITTED bits:", intron, slice.bitscore,
                print "%1.2f-%1.2f-%1.2f" % (
                        self.identityscore,
                        identity_score_surrounding,
                        slice.identityscore ),
                print "(%1.2f-%1.2f)" % ( ids_ratio_1, ids_ratio_2 )
            else:
                pass
            ####################################################################

        # return the list with stopless3n introns
        return accepted_stopless3nsQ

    # end of function detect_aligned_stopless3n_intron_in_query


    def detect_aligned_stopless3n_intron_in_sbjct(self,verbose=False,**kwargs):
        """ """
        # edit **kwargs dictionary for some forced attributes
        _update_kwargs(kwargs,KWARGS_PACBPORF_ALIGNED_STOPLESS_3N_INTRONS)

        mind = self._get_original_alignment_pos_start().sbjct_dna_start
        maxd = self._get_original_alignment_pos_end().sbjct_dna_end
        mina = self._get_original_alignment_pos_start().sbjct_dna_start
        maxa = self._get_original_alignment_pos_end().sbjct_dna_end

        if maxd-mind < 75:
            # do not search for a stopless3n intron in a query of only 75nt
            # It is not likely that there is space for a stopless3n intron
            # in this alignment. It would represent two tiny exons separated
            # by an stopless3n intron (so in the same Orf!). Yeah right ;-)
            return []

        stopless3nsS = find_stopless3n_introns_on_orf(self.orfS,
                min_donor_pos = mind, max_donor_pos = maxd,
                min_acceptor_pos= mina, max_acceptor_pos= maxa,
                **kwargs)
        stopless3nsS = _filter_overlapping_stopless3n_introns(stopless3nsS)

        accepted_stopless3nsS = []
        at_ratio_pacbporf =  self.get_at_ratio_sbjct()
        for intron in stopless3nsS:
            # re-check it's position towards the PacbPORF
            if intron.start <= mind: continue
            if intron.end   <= mind: continue
            if intron.start >= maxd: continue
            if intron.end   >= maxd: continue

            if kwargs['has_at_increase']:
                at_ratio_intron = intron.get_at_ratio()
                if at_ratio_pacbporf > at_ratio_intron:
                    continue

            # if here, get pacbp slice object of this putative intron
            pacbp = conversion.pacbporf2pacbp(self)


            try:

                slice = pacbp.returnslice(
                        self.orfS.dnapos2aapos(intron.start),
                        self.orfS.dnapos2aapos(intron.end),
                        coords_on='sbjct'
                        )


            except:

                print self
                print pacbp
                pacbp.print_protein()
                print intron, intron.coords(), (mind,maxd),(mina,maxa)
                print self.orfS.dnapos2aapos(intron.start)
                print self.orfS.dnapos2aapos(intron.end)

                slice = pacbp.returnslice(
                        self.orfS.dnapos2aapos(intron.start),
                        self.orfS.dnapos2aapos(intron.end),
                        coords_on='sbjct'
                        )



            # calculate identityscore of the surrounding PacbP(ORF) parts
            # and ratios of identityscores.
            # ====================================== self.identityscore
            # ========================         ===== identity_score_surrounding
            #                         =========      slice.identityscore

            if  self.get_unextended_length() <= slice.length:
                # Slice is as long or even longer as the pacbporf itself !?
                # This is NOT what we are searching for. And, identical
                # lengths will cause ZeroDivisionError
                continue

            identity_score_surrounding = ( float(self.identity-slice.identity)+\
                    0.5*( self.similarity - slice.similarity ) ) /\
                    ( self.get_unextended_length() - slice.length )
            # calculate ids ratios; max([0.01,...]) to avoid ZeroDivisionError!
            ids_ratio_1 = identity_score_surrounding / self.identityscore
            ids_ratio_2 = self.identityscore / max([0.01,slice.identityscore])

            is_accepted = False
            if slice.bitscore <= 0:
                # append to accepted stopless3n intron list
                accepted_stopless3nsS.append( intron )
                is_accepted = True
            elif ids_ratio_1 >= 1.10 and ids_ratio_2 >= 1.30:
                # append to accepted stopless3n intron list
                accepted_stopless3nsS.append( intron )
                is_accepted = True
            else:
                pass

            ####################################################################
            if verbose and is_accepted:
                print "S", self, self.get_at_ratio_sbjct()
                print "S", intron, intron.get_at_ratio()
                print slice
                sliceorf = conversion.pacbp2pacbporf(slice,self.orfQ,self.orfS)
                sliceorf._append_one_left_position()
                sliceorf._append_one_left_position()
                sliceorf._append_one_right_position()
                sliceorf._append_one_right_position()
                print "S", sliceorf, "(%1.2f-%1.2f)" % (ids_ratio_1,ids_ratio_2)
                sliceorf.print_protein_and_dna()
            elif verbose:
                print "OMITTED bits:", intron, slice.bitscore,
                print "%1.2f-%1.2f-%1.2f" % (
                        self.identityscore,
                        identity_score_surrounding,
                        slice.identityscore ),
                print "(%1.2f-%1.2f)" % ( ids_ratio_1, ids_ratio_2 )
            else:
                pass
            ####################################################################

        # return the list with stopless3n introns
        return accepted_stopless3nsS

    # end of function detect_aligned_stopless3n_intron_in_sbjct


    def detect_conserved_stopless3n_intron(self,verbose=False,**kwargs):
        """ """
        # edit **kwargs dictionary for some forced attributes
        _update_kwargs(kwargs,KWARGS_PACBPORF_ALIGNED_STOPLESS_3N_INTRONS)

        # find stopless3n introns in query
        mindQ = self._get_original_alignment_pos_start().query_dna_start
        maxdQ = self._get_original_alignment_pos_end().query_dna_end
        minaQ = self._get_original_alignment_pos_start().query_dna_start
        maxaQ = self._get_original_alignment_pos_end().query_dna_end

        if maxdQ-mindQ < 75:
            # do not search for a stopless3n intron in a query of only 75nt
            # It is not likely that there is space for a stopless3n intron
            # in this alignment. It would represent two tiny exons separated
            # by an stopless3n intron (so in the same Orf!). Yeah right ;-)
            return []

        # do NOT perform _filter_overlapping_stopless3n_introns()
        stopless3nsQ = find_stopless3n_introns_on_orf(self.orfQ,
                min_donor_pos = mindQ, max_donor_pos = maxdQ,
                min_acceptor_pos= minaQ, max_acceptor_pos= maxaQ,
                **kwargs)

        # find stopless3n introns in sbjct
        mindS = self._get_original_alignment_pos_start().sbjct_dna_start
        maxdS = self._get_original_alignment_pos_end().sbjct_dna_end
        minaS = self._get_original_alignment_pos_start().sbjct_dna_start
        maxaS = self._get_original_alignment_pos_end().sbjct_dna_end

        if maxdS-mindS < 75:
            # do not search for a stopless3n intron in a query of only 75nt
            # It is not likely that there is space for a stopless3n intron
            # in this alignment. It would represent two tiny exons separated
            # by an stopless3n intron (so in the same Orf!). Yeah right ;-)
            return []

        # do NOT perform _filter_overlapping_stopless3n_introns()
        stopless3nsS = find_stopless3n_introns_on_orf(self.orfS,
                min_donor_pos = mindS, max_donor_pos = maxdS,
                min_acceptor_pos= minaS, max_acceptor_pos= maxaS,
                **kwargs)

        # test if any are found
        if not stopless3nsQ or not stopless3nsS:
            return []

        # return list for conserved_stopless3n_introns
        conserved_stopless3n_introns = []

        # max splice site offset
        offset = kwargs['max_nt_offset']
        at_ratio_pacbporf_query = self.get_at_ratio_query()
        at_ratio_pacbporf_sbjct = self.get_at_ratio_sbjct()
        for intronQ in stopless3nsQ:
            # re-check it's position towards the PacbPORF
            if intronQ.start <= mindQ: continue
            if intronQ.end   <= mindQ: continue
            if intronQ.start >= maxdQ: continue
            if intronQ.end   >= maxdQ: continue

            if kwargs['has_at_increase']:
                at_ratio_intron_query = intronQ.get_at_ratio()
                if at_ratio_pacbporf_query > at_ratio_intron_query:
                    continue

            # if here, get pacbp slice object of this putative intron
            pacbp = conversion.pacbporf2pacbp(self)
            slice = pacbp.returnslice(
                    self.orfQ.dnapos2aapos(intronQ.start),
                    self.orfQ.dnapos2aapos(intronQ.end),
                    coords_on='query'
                    )
            for intronS in stopless3nsS:
                # check phase compatibility (identity!)
                if intronS.phase != intronQ.phase: continue

                # re-check it's position towards the PacbPORF
                if intronS.start <= mindS: continue
                if intronS.end   <= mindS: continue
                if intronS.start >= maxdS: continue
                if intronS.end   >= maxdS: continue


                # do at ratio test if requested for
                if kwargs['has_at_increase']:
                    at_ratio_intron_sbjct = intronS.get_at_ratio()
                    if at_ratio_pacbporf_sbjct > at_ratio_intron_sbjct:
                        continue

                # get range of elegiable sbjct intron positions
                posS = self.orfS.aapos2dnapos(slice.sbjct_start)
                posE = self.orfS.aapos2dnapos(slice.sbjct_end)
                offset = kwargs['max_nt_offset']
                rangeS = range(posS-offset,posS+offset+1)
                rangeE = range(posE-offset,posE+offset+1)
                if intronS.start in rangeS and intronS.end in rangeE:
                    # calculate distance here and set to _distance attribute
                    distD = self.get_distance_aligned_nucleotide_positions(
                            intronQ.donor.pos,intronS.donor.pos)
                    distA = self.get_distance_aligned_nucleotide_positions(
                            intronQ.acceptor.pos,intronS.acceptor.pos)
                    intronQ._distance = distD+distA
                    intronS._distance = distD+distA

                    # append to list of conserved stopless3n introns
                    conserved_stopless3n_introns.append( ( intronQ, intronS ) )

                    ############################################################
                    if verbose:
                        print "QS", self, self.get_at_ratio_query()
                        print "Q.", intronQ, self.get_at_ratio_query(),
                        print intronQ.get_at_ratio()
                        print "S.", intronS, self.get_at_ratio_sbjct(),
                        print intronS.get_at_ratio()
                        print intronQ.coords(), intronQ.get_branchpoint_nt_distance(),
                        print intronS.coords(), intronS.get_branchpoint_nt_distance()
                        print slice
                        sliceorf = conversion.pacbp2pacbporf(slice,
                                self.orfQ,self.orfS)
                        sliceorf._append_one_left_position()
                        sliceorf._append_one_left_position()
                        sliceorf._append_one_right_position()
                        sliceorf._append_one_right_position()
                        print "QS", sliceorf
                        sliceorf.print_protein_and_dna()
                    ############################################################

        # return list with conserved_stopless3n_introns
        return conserved_stopless3n_introns

    # end of function detect_conserved_stopless3n_intron


    def has_query_orf_methionine(self):
        """ """
        return self.orfQ.has_methionine()

    # end of function has_query_orf_methionine


    def has_sbjct_orf_methionine(self):
        """ """
        return self.orfS.has_methionine()

    # end of function has_sbjct_orf_methionine


    def has_orfs_methionines(self):
        """ """
        if not self.orfQ.has_methionine():      return False
        elif not self.orfS.has_methionine():    return False
        else:                                   return True

    # end of function has_orfs_methionines


    def has_upstream_methionines(self):
        """ """
        eposobj = self._get_original_alignment_pos_end()
        posQ = eposobj.query_pos + 1
        posS = eposobj.sbjct_pos + 1
        if not self.orfQ.has_start_upstream_of(posQ):   return False
        elif not self.orfS.has_start_upstream_of(posS): return False
        else:                                           return True
       
    # end of function has_upstream_methionines


    def has_query_orf_tss(self):
        """ """
        return self.orfQ.has_translationalstartsite()

    # end of function has_query_orf_tss


    def has_sbjct_orf_tss(self):
        """ """
        return self.orfS.has_translationalstartsite()

    # end of function has_sbjct_orf_tss


    def has_orfs_tss(self):
        """ """
        if not self.orfQ.has_translationalstartsite():  return False
        elif not self.orfS.has_translationalstartsite():return False
        else:                                           return True

    # end of function has_orfs_tss


    def has_upstream_tss(self):
        """ """
        eposobj = self._get_original_alignment_pos_end()
        posQ = eposobj.query_pos + 1
        posS = eposobj.sbjct_pos + 1
        if not self.orfQ.has_tss_upstream_of(posQ):     return False
        elif not self.orfS.has_tss_upstream_of(posS):   return False
        else:                                           return True
       
    # end of function has_upstream_tss


    def get_sbjct_upstream_tss(self):
        """ """
        ordered_sites = []
        eposobj = self._get_original_alignment_pos_end()
        for site in self.orfS._tss_sites:
            if site.pos <= eposobj.sbjct_dna_end:
                ordered_sites.append( site )
        if ordered_sites:
            return olba(ordered_sites,order_by='pssm_score',reversed=True)[0]
        else:
            return None
       
    # end of function get_sbjct_upstream_tss


    def get_query_upstream_tss(self):
        """ """
        ordered_sites = []
        eposobj = self._get_original_alignment_pos_end()
        for site in self.orfQ._tss_sites:
            if site.pos <= eposobj.query_dna_end:
                ordered_sites.append( site )
        if ordered_sites:
            return olba(ordered_sites,order_by='pssm_score',reversed=True)[0]
        else:
            return None
       
    # end of function get_query_upstream_tss


    def get_projected_tailing_stop_aa_difference(self):
        """ """
        endpos = self._get_original_alignment_pos_end()
        qstop  = self.orfQ.protein_endPY - endpos.query_pos
        sstop  = self.orfS.protein_endPY - endpos.sbjct_pos
        return abs(qstop-sstop)

    # end of function get_projected_tailing_stop_aa_difference


    def get_projected_tailing_stop_nonaligned_aa_difference(self):
        """ """
        endpos = self._get_original_alignment_pos_end()
        qstop  = self.orfQ.protein_endPY - endpos.query_pos
        sstop  = self.orfS.protein_endPY - endpos.sbjct_pos
        return max([qstop,sstop])

    # end of function get_projected_tailing_stop_nonaligned_aa_difference


    def gff_projected_tailing_stop(self,gff={}):
        """
        Return gff tuple of ProjectedTailingStop of sbjct stop on query.
        """
        # update self._gff with gff data
        self._gff.update(gff)

        # project the stop-site of this sbjct organism
        endpos = self._get_original_alignment_pos_end()
        endpos_stop_dist = self.orfS.endPY - endpos.sbjct_dna_end
        projected_stop_start = endpos.query_dna_end + endpos_stop_dist
        sta,sto = projected_stop_start+1,projected_stop_start+3

        column9string = ""
        # take care for the formatting of column9data
        if self._gff['column9data']:
            tmp = []
            for key,value in self._gff['column9data'].iteritems():
                if str(value).find(" ") > -1: value = "'%s'" % value
                tmp.append( "; %s %s" % (key,value) )
            column9string = "".join(tmp)

        # and return the GFF tuple
        return ( self._gff['fref'], self._gff['fsource'], 'ProjectedTailingStop',
                 sta, sto, '.', self._gff['fstrand'], self._gff['fphase'],
                 "%s %s%s" % ( self._gff['gclass'], self._gff['gname'], column9string ) )

    # end of function gff_projected_tailing_stop


    def gff_projected_leading_stop(self,gff={}):
        """
        Return gff tuple of ProjectedLeadingStop of sbjct stop on query.
        """
        # update self._gff with gff data
        self._gff.update(gff)
        # project the leading stop-site of sbjct organism on query organism
        stapos = self._get_original_alignment_pos_start()
        stapos_stop_dist = stapos.sbjct_dna_start - self.orfS.startPY + 3
        projected_stop_start = stapos.query_dna_start - stapos_stop_dist
        sta,sto = projected_stop_start+1,projected_stop_start+3

        column9string = ""
        # take care for the formatting of column9data
        if self._gff['column9data']:
            tmp = []
            for key,value in self._gff['column9data'].iteritems():
                if str(value).find(" ") > -1: value = "'%s'" % value
                tmp.append( "; %s %s" % (key,value) )
            column9string = "".join(tmp)

        # and return the GFF tuple
        return ( self._gff['fref'], self._gff['fsource'], 'ProjectedLeadingStop',
                 sta, sto, '.', self._gff['fstrand'], self._gff['fphase'],
                 "%s %s%s" % ( self._gff['gclass'], self._gff['gname'], column9string ) )

    # end of function gff_projected_leading_stop

# end of class PacbPORF




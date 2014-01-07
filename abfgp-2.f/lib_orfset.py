"""
Parser functions for EMBOSS getorf and classes describing sets of Orf objects
"""

# Module metadata
__authors__ = "Ate van der Burgt"
__license__ = "MIT"

# Python Imports
from sets import Set

# Miscelaneous Imports
from lib_fasta import parseFasta
from parsers.getorf import getorf

# Gene package Imports
from gene.orf import TcodeCodingOrf
from gene.validators import IsProperStrand
from gene.gene_exceptions import (
    InproperlyAppliedArgument,
    MissingAttributeValue,
    )
    
# Global Variable Imports
from settings.executables import EXECUTABLE_GETORF_VERSION


class OrfSet:
    def __init__(self,sequence="",start=None):
        """
	Initialize OrfSet object
	
        @type  sequence: string
        @param sequence: DNA sequence string that served as input fot getorf

        @type  start: integer
        @param start: zero-positive integer, nt offset of DNA sequence slice
        """
        # sequence    genomic sequence that is fedded to getorf
        # accession   accession of this sequence
        # start       start coordinate (offset) of the sequence slice
        # stop        stop coordinate of the sequence slice (GFF coordinate)
	# end	    alias for self.stop (GFF coordinate)
        # strand      strand of the sequence slice
        # startPY     start coordinate in PYTHON coordinate system (==start-1)
        # stopPY      stop coordinate in PYTHON coordinate system (==stop)
	# endPY       alias for self.stopPY

	# main attributes of OrfSet
        self.orfs       = []
        self.sequence   = sequence

	# metadata attributes of OrfSet
        self.version    = None
        self.accession  = None

	# absolute coordinates relating to self.sequence
        self.strand     = None      # TODO: in use? not in __init__ or parse !?
        self.start      = start
        self.stop       = None      # TODO: in use? not in __init__ or parse !?
	self.end	= self.stop # TODO: in use? not in __init__ or parse !?
        self.startPY    = None      # TODO: in use? not in __init__ or parse !?
        self.stopPY     = None      # TODO: in use? not in __init__ or parse !?
	self.endPY	= self.stopPY

    # end of function __init__


    def get_orf_by_id(self,id):
        """
        Get Orf object from list with Orfs

        @type  id: * (integer)
        @param id: (integer) id of a Orf object

        @type  orf: Orf object (or None)
        @param orf: Orf object (or None)
        """
        for orf in self.orfs:
            if orf.id == id:
                return orf
        else:
            return None

    # end of function get_orf_by_id


    def tofasta(self):
        """
        Returns a text string of all orf's (protein) sequences

        @rtype:  text
        @return: fasta text string of all orfs in the getorf output
        """
        return "\n".join([ orf.tofasta() for orf in self.orfs ])

    # end of function tofasta


    def tofastadict(self):
        """
        Returns a dictionary all orf's (protein) sequences

        @rtype:  dictionary
        @return: dictionary with all orfs in the getorf output
        """
        mfa = {}
        for orf in self.orfs: 
            mfa[orf.fastaheader] = orf.protein_sequence
        return mfa

    # end of function tofastadict


    def tomaskedfasta(self,coords=[],orflist=[],header_prefix=""):
        """
        Returns a text string of all orf's (protein) sequences, masked by AA coord slices

        @type  coords: list
        @param coords: list of tuples in the shape op (mask_aa_start,mask_aa_end)

	@rtype:  text
	@return: fasta text string of all orfs in the getorf output
        """
        # create set of mask coordinates
        mask = Set()
        for (start,end) in coords: mask.update(range(start,end+1))
        seqs = []
        if orflist: orfstofasta = orflist
        else:       orfstofasta = self.orfs
        for orf in orfstofasta:
            orfcoords = Set( range( orf.protein_startPY, orf.protein_endPY ) )
            intersect = orfcoords.intersection(mask)
            header = "%(header)s_orf_%(orfid)s" % {
		'header': header_prefix, 'orfid': orf.id }
            if intersect:
                if len(intersect) == orf.protein_length:
                    # complete sequence is masked out!
                    continue
                letters = list( orf.protein_sequence )
                intersect = [ (coord - orf.protein_startPY) for coord in list(intersect) ]
                for coord in intersect: letters[coord] = 'X'
                seqstring = "".join(letters)
                if not seqstring or seqstring.count("X") == len(seqstring): continue
                seqs.append( ">%s\n%s" % (header,seqstring) )
            else:
                seqs.append( orf.tofasta(header=header) )
        return "\n".join(seqs)

    # end of function tomaskedfasta


    def get_elegiable_orfs(self,**kwargs):
        """
        Return a selection of elegiable orfs based on several parameters

        @attention: for function arguments, see get_elegiable_orfs() function

        @rtype:  list
        @return: list of elegiable Orfs
        """
        return get_elegiable_orfs(self.orfs,**kwargs)

    # end of function get_elegiable_orfs


    def filter(self,**kwargs):
        """
        Filter the Orf objects in this OrfSet (and remove failed Orf objects)
	
        @attention: for function arguments, see get_elegiable_orfs() function
        """
        self.orfs = get_elegiable_orfs(self.orfs,**kwargs)

    # end of function filter

# end of class OrfSet


class GetorfOrfSet(OrfSet):
    """
    """
    def __init__(self,sequence="",start=None,version=EXECUTABLE_GETORF_VERSION):
        """
	Initialize GetorfOrfSet object from EMBOSS getorf result
	
        @type  sequence: string
        @param sequence: DNA sequence string that served as input fot getorf

        @type  start: integer
        @param start: zero-positive integer, nt offset of DNA sequence slice

        @type  version: *
        @param version: EMBOSS Getorf version

        """
	OrfSet.__init__(self,sequence=sequence,start=start)
        self.version = version

    # end of function __init__


    def parse(self,filehandle,sequence="",accession="",start="",strand="+"):
        """
        Parse output of a EMBOSS getorf output file

        @attention: implemented for ORFs from STOP to STOP

        @type  filehandle: filehandle
        @param filehandle: filehandle to opened getorf output file

        @type  sequence: string
        @param sequence: DNA sequence string that served as input fot getorf

        @type  accession: string
        @param accession: accession to use for the DNA sequence and it orfs

        @type  start: integer
        @param start: zero-positive integer, nucletiotide offset of DNA sequence slice

        @type  strand: string
        @param strand: strand of the DNA sequence (default '+')
        """
        # input validation
        IsProperStrand(strand)
        if not self.sequence:
            raise MissingAttributeValue, "sequence not applied"

        self.sequence  = sequence
        self.accession = accession
        self.strand = strand

        if (start or start == 0) and type(start) == type(int()):
            self.start = start
	    self.end   = self.start+len(self.sequence)

        # parse the ouput file
        self._parsefilehandle(filehandle)

    # end of function parse


    def getorfs(self,strand="+"):
        """
        Run & Parse EMBOSS getorf from self.sequence 

        @attention: implemented for ORFs from STOP to STOP

        @type  strand: string
        @param strand: strand of the DNA sequence (default '+')
        """
        # input validation
        IsProperStrand(strand)
        self.strand = strand

        # parse the ouput file
        seqs = parseFasta( getorf( sequence=self.sequence ).split("\n") )
	self._dict_to_orfs( seqs )
	
    # end of function getorfs


    def _parsefilehandle(self,filehandle):
        """
	Parse getorf output file(handle) to Orf objects in OrfSet

        @type  filehandle: filehandle
        @param filehandle: filehandle to opened getorf output file
        """
        self._dict_to_orfs( parseFasta(filehandle.readlines() ) )

    # end of function _parsefilehandle


    def _dict_to_orfs(self,sequencedict,sequence=None):
	"""
	Convert a fasta sequence dictionary to Orf objects

        @type  sequencedict: dict
        @param sequencedict: fasta dictionary (values=sequence strings)

        @type  sequence: string (or None)
        @param sequence: sequence string on which Orfs are located
	"""
	if sequence: self.sequence = sequence
        for h,s in sequencedict.iteritems():
            force_id = len(self.orfs)
            self.orfs.append( TcodeCodingOrf(h,s,
		inputgenomicsequence=self.sequence, force_id=force_id) )

    # end of function _dict_to_orfs


    def add_novel_orf(self,start,end,aasequence):
        """
        Add a single Orf to this orfset object
        
        @type  start: int
        @param start: AA start coordinate of Orf

        @type  end: int
        @param end: AA end coordinate of Orf

        @type  aasequence: string
        @param aasequence: AA sequence string of this Orf

        @attention: use with care!
        @attention: apply only when a tiny Orf is regularly omitted
        """
        if self.orfs:
            new_orf_id = max([ orf.id for orf in self.orfs ])+1
            inputgenomicsequence = self.orfs[0].inputgenomicsequence
        else:
            new_orf_id = 0
            inputgenomicsequence = ""

        # append to OrfSet
        self.orfs.append( TcodeCodingOrf("",aasequence,
            inputgenomicsequence=inputgenomicsequence,
            force_id=new_orf_id, start = start, end = end ) )

        # return this created OrfObject
        return self.orfs[-1]
        
    # end of function add_novel_orf

# end of class GetorfOrfSet


def get_elegiable_orfs(orflist,min_orf_length=None,max_orf_length=None,
    min_orf_start=None,min_orf_end=None,max_orf_start=None,max_orf_end=None,
    acceptorfids=[],rejectorfids=[],has_starts=None,
    tcode_symbolic_in=[],tcode_is_coding=None,tcode_is_noncoding=None,
    is_internalexon=None,is_firstexon=None,
    is_finalexon=None,is_singleexon=None):
    """
    Return a selection of elegiable orfs based on several parameters

    @type  min_orf_length: positive integer (or None)
    @param min_orf_length: minimal nt length of Orf

    @type  max_orf_length: positive integer (or None)
    @param max_orf_length: maximal nt length of Orf

    @type  min_orf_start: positive integer (or None)
    @param min_orf_start: minimal start nt coordinate of Orf

    @type  max_orf_start: positive integer (or None)
    @param max_orf_start: maximal start nt coordinate of Orf

    @type  min_orf_end: positive integer (or None)
    @param min_orf_end: minimal end nt coordinate of Orf

    @type  min_orf_end: positive integer (or None)
    @param max_orf_end: maximal end nt coordinate of Orf

    @type  acceptorfids: list with integers
    @param acceptorfids: list of orf ids to accept

    @type  rejectorfids: list with integers
    @param rejectorfids: list of orf ids to reject 

    @type  has_starts: Boolean (or None)
    @param has_starts: True returns only Orfs with 'M' symbols, False all except these. 

    @type  tcode_symbolic_in: list
    @param tcode_symbolic_in: list of tcode symbolic statusses to return ['C','N','?','!','$']

    @type  tcode_is_coding: Boolean (or None)
    @param tcode_is_coding: True returns only Coding Orfs, False returns all except coding orfs

    @type  tcode_is_noncoding: Boolean (or None)
    @param tcode_is_noncoding: True returns only Non-Coding Orfs, False returns all except non-coding orfs

    @type  is_internalexon: NoneBoolean
    @param is_internalexon: only return Orfs than are potentionally internal exons

    @type  is_finalexon: NoneBoolean
    @param is_finalexon: only return Orfs than are potentionally final exons

    @type  is_firstexon: NoneBoolean
    @param is_firstexon: only return Orfs than are potentionally first exons

    @type  is_singleexon: NoneBoolean
    @param is_singleexon: only return Orfs than are potentionally single exons

    @rtype:  list
    @return: list of Orf objects
    """
    CHECK_FOR_TCODE_IS_APPLIED = True   # Tcode is applied on Orfs;
                                        # todo: write a function for this...
    retlist = []
    for orf in orflist:
        if min_orf_length:
            if orf.length < min_orf_length:     continue
            else:                               pass
        if max_orf_length:
            if orf.length > max_orf_length:     continue
            else:                               pass
        if min_orf_start:
            if orf.startPY >= min_orf_start:    pass
            else:                               continue
        if min_orf_end:
            if orf.endPY >= min_orf_end:        pass
            else:                               continue
        if max_orf_start:
            if orf.startPY <= max_orf_start:    pass
            else:                               continue
        if max_orf_end:
            if orf.endPY <= max_orf_end:        pass
            else:                               continue
        if acceptorfids:
            if orf.id in acceptorfids:          pass
            else:                               continue
        if rejectorfids:
            if orf.id in rejectorfids:          continue 
            else:                               pass
        if has_starts == False:
            if orf.has_start():                 continue
        elif has_starts == True:
            if not orf.has_start():             continue
        else:                                   pass

	if orf.__class__.__name__ in ['CodingOrf','TcodeCodingOrf']:

            # check for potential internal exon
            if is_internalexon == None:                 pass
            elif is_internalexon:
                if orf.can_orf_encode_internalexon():   pass
                else:                                   continue
            else:
                if orf.can_orf_encode_internalexon():   continue

            # check for potential first exon
            if is_firstexon == None:                    pass
            elif is_firstexon:
                if orf.can_orf_encode_firstexon():      pass
                else:                                   continue
            else:
                if orf.can_orf_encode_firstexon():      continue

            # check for potential final exon
            if is_finalexon == None:                    pass
            elif is_finalexon:
                if orf.can_orf_encode_finalexon():      pass
                else:                                   continue
            else:
                if orf.can_orf_encode_finalexon():      continue

            # check for potential single exon
            if is_singleexon == None:                   pass
            elif is_singleexon:
                if orf.can_orf_encode_singleexon():     pass
                else:                                   continue
            else:
                if orf.can_orf_encode_singleexon():     continue

        else:
            # hmm.. Orf object lack coding information
            pass


	if orf.__class__.__name__ in ['TcodeOrf','TcodeCodingOrf']:
            if tcode_is_coding == False:
                if orf.tcode_is_coding():
                    continue
            elif tcode_is_coding == True:
                if not orf.tcode_is_coding():
                    continue
            else:   pass
            if tcode_is_noncoding == False:
                if orf.tcode_is_noncoding():
                    continue
            elif tcode_is_noncoding == True:
                if not orf.tcode_is_noncoding():
                    continue
            else:   pass
            if tcode_symbolic_in:
                if orf.tcode_symbolic() not in tcode_symbolic_in:
                    continue 
            else:   pass
        else:
            # hmm.. Orf object lack tcode information
            pass

        # if this point is reached without hitting a `continue`
        # statement, this Orf is eligiable!
        retlist.append( orf )

    # return the subselection of orfs
    return retlist

# end of function get_elegiable_orfs



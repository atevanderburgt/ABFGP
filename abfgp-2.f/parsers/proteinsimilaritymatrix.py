"""
ProteinSimilarityMatrix parser function and classes
"""

# Module metadata
__authors__ = "Ate van der Burgt"
__license__ = "MIT"

##############################################################################
#
#       blast matrices information
#       http://flybase.net/static_pages/blast/matrix_info.html
#
#     Query length     Substitution matrix     Gap costs
#     ------------     -------------------     ---------
#     <35              PAM-30                  ( 9,1)
#     35-50            PAM-70                  (10,1)
#     50-85            BLOSUM-80               (10,1)
#     >85              BLOSUM-62               (11,1)
#
##############################################################################
#
# blosum62 = parse_emboss_matrix_file("/opt/ab/share/EMBOSS/data/EBLOSUM62")
# pam30 = parse_emboss_matrix_file("/opt/ab/share/EMBOSS/data/EPAM30")
#
##############################################################################

# Python Imports
from os.path import join as osPathJoin

# import paths to BLOSUM and PAM matrices
from executables import EMBOSS_DATA_DIRECTORY
BLOSUM45_PATH = osPathJoin(EMBOSS_DATA_DIRECTORY,"EBLOSUM45")
BLOSUM62_PATH = osPathJoin(EMBOSS_DATA_DIRECTORY,"EBLOSUM62")
BLOSUM80_PATH = osPathJoin(EMBOSS_DATA_DIRECTORY,"EBLOSUM80")
PAM30_PATH    = osPathJoin(EMBOSS_DATA_DIRECTORY,"EPAM30")
PAM70_PATH    = osPathJoin(EMBOSS_DATA_DIRECTORY,"EPAM70")


def parse_emboss_matrix_file(fname):
    """
    Parse an EMBOSS protein similarity matrix file
    Supported types are BLOSUM and PAM matrices
    """
    lines = []
    matrix = {}
    aas = []
    for line in open(fname).readlines():
        if line[0] == "#":
            continue
        elif line[0] == " ":
            line = line.strip()
            if not line: continue
            # the COLUMNS of amino acids
            aas = line.strip().replace("   "," ").replace("  "," ").split(" ")
        else:
            line = line.strip()
            if not line: continue
            # a ROW with scores
            line = line.replace("   "," ").replace("  "," ").strip()
            scores = line.split(" ")
            aa = scores.pop(0)
            scores = [ int(score) for score in scores ]
            for pos in range(0,len(aas)):
                matrix[aa+aas[pos]] = scores[pos]
    # ready, return
    return matrix

# end of function parse_emboss_matrix_file


def alignment2bitscore(_query,_match,_sbjct,matrix={},gap_opening=11,gap_extension=1):
    """
    Calculate the alignment (bit) score of
    two ALREADY ALIGNED protein sequences.
    REMARK: _match is deprecated
    """
    if not matrix:
        raise "no matrix dictionairy specified"
    if len(_query) != len(_sbjct):
        print 'Q', len(_query), "'%s'" % _query
        print 'S', len(_sbjct), "'%s'" % _sbjct
        raise "input sequence lengths NOT identical!"
    else:
        # define identical length
        length = len(_query)
    # take uppers
    _query = _query.upper()
    _sbjct = _sbjct.upper()
    # okay , enough error checking; do the alignment score!
    bitscore = 0
    for pos in range(0,length):
        q = _query[pos]
        s = _sbjct[pos]
        if q == "-" or s == "-":
            # a gap!
            if pos == 0 or (_query[pos-1] == "-" or _sbjct[pos-1] == "-"):
                # a gap opening
                bitscore -= gap_opening
            else:
                # a gap extension
                bitscore -= gap_extension
        else:
            bitscore += matrix[q.upper()+s.upper()]
    # and return
    return bitscore    

# end of function alignment2bitscore


def make_alignment_match(q,s,matrix={}):
    """
    Creates a BlastP alignment match string between a query and a sbjct sequence
    >0  Identity        '*'
    >0  Non-identity    '+'
    0   Non-identity    '.'
    <0  Other           ' '

    @attention: see _generate_alignment_match for documentation
    """

    return _generate_alignment_match(q,s,matrix=matrix,symbols=['*','+','.',' '])

    match = []
    q,s = q.upper(), s.upper()
    for pos in range(0,len(q)):
        if q[pos] == '-' or s[pos] == '-':
            # deal with gaps in the AA sequence
            match.append(" ")
            continue
        score = matrix[q[pos]+s[pos]]
        if q[pos] == s[pos] and score > 0:
            # identical amino acids
            match.append("*")
        elif q[pos] == s[pos] and score <= 0:
            # identical char, but no amino acid (stop-codon?)
            match.append(" ")
        elif score == 0:
            # indifferent aa substitution
            match.append(".")
        else:
            if score > 0:
                match.append("+")
            else:
                match.append(" ")
    # and return
    return "".join(match)

# end of function make_alignment_match


def make_blastp_alignment_match(q,s,matrix={}):
    """
    Creates a BlastP alignment match string between a query and a sbjct sequence

    @attention: see _generate_alignment_match for documentation
    @attention: alias for make_alignment_match 
    """
    return _generate_alignment_match(q,s,matrix=matrix,symbols=['*','+','.',' '])

# end of function make_blastp_alignment_match 


def make_clustalw_alignment_match(q,s,matrix={}):
    """
    Creates a ClustalW alignment match string between a query and a sbjct sequence
    >0  Identity        '*'
    >0  Non-identity    ':'
    0   Non-identity    '.'
    <0  Other           ' '

    """
    return _generate_alignment_match(q,s,matrix=matrix,symbols=['*',':','.',' '])


    match = []
    q,s = q.upper(), s.upper()
    for pos in range(0,len(q)):
        score = matrix[q[pos]+s[pos]]
        if q[pos] == s[pos] and score > 0:
            # identical amino acids
            match.append("*")
        elif q[pos] == s[pos] and score <= 0:
            # identical char, but no amino acid (stop-codon?)
            match.append(" ")
        elif score == 0:
            # indifferent aa substitution
            match.append(".")
        else:
            if score > 0:
                match.append(":")
            else:
                match.append(" ")
    # and return
    return "".join(match)

# end of function make_clustalw_alignment_match


def _generate_alignment_match(q,s,matrix={},symbols=['*','+','.',' ']):
    """
    Creates an alignment match string between a query and a sbjct protein sequence based on symbols

    @rtype:  q: string 
    @return: q: aligned(!) protein query

    @rtype:  s: string
    @return: s: aligned(!) protein subject

    @rtype:  matrix: dict 
    @return: matrix: similarity matrix (ProteinSimilarityMatrix.matrix)

    @type  symbols: list 
    @param symbols: 4 strings in a list: identity, similarity score >0, similarity score=0, other

    @attention: use symbols=['*','+','.',' '] for BlastP
    @attention: use symbols=['*',':','.',' '] for ClustalW 

    @rtype:  match: string
    @return: match: match string of the aligned protein sequence pair
    """
    match = []
    q,s = q.upper(), s.upper()
    for pos in range(0,len(q)):
        if q[pos] == '-' or s[pos] == '-':
            # deal with gaps in the AA sequence
            match.append( symbols[3] )
            continue
        # get aligned residue score from matrix
        score = matrix[q[pos]+s[pos]]
        if q[pos] == s[pos] and score > 0:
            # identical amino acids
            match.append( symbols[0] )
        elif q[pos] == s[pos] and score <= 0:
            # identical char, but no amino acid (stop-codon?)
            match.append( symbols[3] )
        elif score == 0:
            # indifferent aa substitution
            match.append( symbols[2] )
        else:
            if score > 0:
                match.append( symbols[1] )
            else:
                match.append( symbols[3] )
    # and return
    return "".join(match)

# end of function _generate_alignment_match


class ProteinSimilarityMatrix:
    """ """
    def __init__(self,name="BLOSUM62",fname="",gap_opening=11,gap_extension=1):
        """ """
        self.fname      = ""
        self._supported = { 'BLOSUM62'  : (BLOSUM62_PATH,11,1),
                            'BLOSUM80'  : (BLOSUM80_PATH,10,1),
                            'BLOSUM45'  : (BLOSUM45_PATH,11,1), # ?? are gap costs okay??
                            'PAM70'     : (PAM70_PATH,10,1),
                            'PAM30'     : (PAM30_PATH,9,1),
                            }
        if fname:
            self.fname          = fname
            self.name           = fname.split("/")[-1]
            self.gap_opening    = gap_opening
            self.gap_extension  = gap_extension
        elif name and name.upper() in self._supported.keys():
            self.name           = name.upper()
            self.fname          = self._supported[self.name][0]
            self.gap_opening    = self._supported[self.name][1]
            self.gap_extension  = self._supported[self.name][2]
        elif name and name.upper() not in self._supported:
            raise "Matrix '%s' not supported; not in %s" % (
                name, self._supported )
        else:
            raise "Specify a `name` or a filename `fname`"

        # and now parse the matrix file
        # if a non-existing filename is specified,
        # an error will be returned automatically
        self.matrix = self._parse_emboss_matrix_file()

    # end of function __init__

    def __str__(self):
        return "<ProtSimMatrix %s>" % self.name
    # end of function __str__


    def _parse_emboss_matrix_file(self):
        """ """
        return parse_emboss_matrix_file(self.fname)

    # end of function _parse_emboss_matrix_file


    def scorealignment(self,q,s):
        """
        score the alignment between a ALREADY ALIGNED (!)
        Query and Sbjct protein sequence of identical length
        """
        return alignment2bitscore(q, "", s, matrix=self.matrix,
            gap_opening=self.gap_opening,
            gap_extension=self.gap_extension )

    # end of function scorealignment


    def alignmentmatch(self,q,s):
        """ """
        return make_alignment_match(q,s,matrix=self.matrix)

    # end of function alignmentmatch


    def scoreaapair(self,aa1,aa2):
        """
        Get score for pair of amino acids
        This function does NOT add up to the
        total bitscore when iterated over the complete
        alignment, because gap_opening is not accounted for.
        """
        try:
            return self.matrix[aa1.upper()+aa2.upper()]
        except:
            # pair not present !? return lowest possible score
            return min( self.matrix.values() )

    # end of function scoreaapair


    def scansbjct(self,q,s,min_bitscore=None,min_bitscore_ratio=None):
        """
        Scan a sbjct sequence for a query sequence for an UNGAPPED ALIGNMENT

        @type  q: string
        @param q: (ungapped) protein sequence string

        @type  s: string
        @param s: (ungapped) protein sequence string

        @type  min_bitscore: number
        @param min_bitscore: minimal bitscore to report back (not recommended to use)

        @type  min_bitscore_ratio: float
        @param min_bitscore_ratio: float 0.0..1.0 of minimal ratio of bitscore vs perfect match to report back (recommended to use)
        """
        if not (min_bitscore_ratio or min_bitscore): raise Exception
        max_bitscore = float(self.scorealignment(q,q))
        qlen=len(q)
        results = []
        for pos in range(0,len(s)-qlen):
            sbjct = s[pos:pos+qlen]
            bitscore = self.scorealignment(q,sbjct)
            ratio = float(bitscore)/max_bitscore
            if min_bitscore_ratio:
                if ratio >= min_bitscore_ratio:
                    match = self.alignmentmatch(q,sbjct)
                    results.append( (ratio, pos, q, match, sbjct, bitscore) )
            else:
                if bitscore >= min_bitscore:
                    results.append( (ratio, pos, q, match, sbjct, bitscore) )
        # return the hits
        results.sort()
        results.reverse()
        return results

    # end of function scansbjct


# end of class ProteinSimilarityMatrix


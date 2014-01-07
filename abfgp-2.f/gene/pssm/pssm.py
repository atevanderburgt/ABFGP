# Module metadata
__authors__ = "Ate van der Burgt"
__license__ = "MIT"

# Imports
from math import log10

# Exceptions
class InproperlySizedPattern(Exception):
    """ applied pattern not of the proper length """
    pass

class NonDNALetter(Exception):
    """ Non-DNA letter occurring in sequence pattern in strict mode """
    pass

def log2(x):
    """ Base 2 logarithm. Example: log2(1024) = 10.0 """
    return log10(x) / log10(2)

# end of function log2

def pssmscore(pattern,IC,strict=False,ignore_unambiguity=False):
    """
    Get PSSM score for a sequence pattern given an IC data structure

    @type  pattern: string
    @param pattern: absolute filepath to an InformationContent file

    @type  IC: []
    @param IC: InformationContent list

    @type  strict: Boolean
    @param strict: if True, raise NonDNALetter if observed in pattern

    @rtype:  float
    @return: pssm score
    """
    pattern = pattern.upper()
    score = 0.0
    for i in range(0,len(pattern)):
        if ignore_unambiguity and round(max(IC[i].values()),1) == 2.0:
            # do not count unambigious positions
            continue
        base = pattern[i]
        if base in ['G','T','C','A']:
            score += IC[i][base]
        elif base in ['X','N']:
            if strict:
                raise NonDNALetter
            else:
                # return lowest possible score
                score += min(IC[i].values())
        elif base == 'U':
            score += IC[i]['T']
        else:
            raise NonDNALetter 
    # return summed PSSM score
    return score

# end of function pssmscore


def get_ic_max_score(ic,ignore_unambiguity=True):
    """ """
    scores = [ max(pos.values()) for pos in ic ]
    if ignore_unambiguity:
        for ii in range(len(scores)-1,-1,-1):
            if round(scores[ii],1) == 2.0:
                _popped = scores.pop(ii)
    return sum(scores)
# end of function get_ic_max_score


def get_ic_min_score(ic,ignore_unambiguity=True):
    """ """
    scores = [ min(pos.values()) for pos in ic ]
    if ignore_unambiguity:
        maxscores = [ max(pos.values()) for pos in ic ]
        for ii in range(len(maxscores)-1,-1,-1):
            if round(maxscores[ii],1) == 2.0:
                _popped = scores.pop(ii)
    return sum(scores)
# end of function get_ic_min_score



def print_ic(ic):
    """ Print human-readable representation of InformationContent list """
    print "\t",
    for base in list("GTAC"):
        print "%s\t\t\t" % base,
    print ""
    for i in range(0, len(ic) ):
        pssm = ic[i]
        print i, "\t",
        for base in list("GTAC"):
            print "%1.6f\t\t" % pssm[base],
        print ""

# end of function print_ic


def obtain_pssm_ic(seqs,start=None,end=None,report=False):
    """
    Generate a InformationContent list from a bunch of sequences

    @type  seqs: dict
    @param seqs: (fatsa) dictionary with sequences **of uniform length**

    @type  start: integer (or None)
    @param start: integer start offset to limit PSSM pattern

    @type  end: integer (or None)
    @param end: integer end offset to limit PSSM pattern

    @type  report: Boolean
    @param report: print status message to STDOUT if True

    @rtype:  ( [], integer )
    @return: InformationContent PSSM list, number of ignored sequences
    """
    _from, _to = 0, len(seqs.values()[0])
    if start: _from = start
    if end:   _to = end+1
    pattern_range = range(_from,_to)
    if report:
        print pattern_range
    # some variables we need
    data = []
    IC = []
    patlen = len(seqs.values()[0])
    for i in range(0,patlen):
        #data.append( {'A':0,'T':0,'C':0,'G':0 } )
        data.append( {'A':0,'T':0,'C':0,'G':0,'-':0 } )
    nns = 0 # sequences with 'N' bases
    for seq in seqs.values():
        if seq.upper().count("N") >= 1:
            nns+=1
            continue 
        for i in range(0,patlen):
            base = seq[i].upper()
            data[i][base]+=1
    for i in pattern_range:
        IC.append( {'A':0,'T':0,'C':0,'G':0 } )
        for base in list("ATGC"):
            freq = float(data[i][base])/float(len(seqs)-nns)
            if freq == 0.0:
                # about to encounter a Math Log error
                freq = 1.0 / float(len(seqs)-nns)
            #print base, "%2.0f" % (100.0*freq), "\t", 
            IC[-1][base] = 2.0-(-log2(freq))
    if report:
        print "%s sequences, %s ignored (N's), pattern length=%s" % (
            len(seqs),nns,len(pattern_range))
    return IC,nns

# end of function obtain_pssm_ic


def obtain_pssm_ic_from_cntdict(cntdict,report=False):
    """
    Generate a InformationContent list from a bunch of sequences

    @type  seqs: dict
    @param seqs: (fatsa) dictionary with sequences **of uniform length**

    @type  start: integer (or None)
    @param start: integer start offset to limit PSSM pattern

    @type  end: integer (or None)
    @param end: integer end offset to limit PSSM pattern

    @type  report: Boolean
    @param report: print status message to STDOUT if True

    @rtype:  ( [], integer )
    @return: InformationContent PSSM list, number of ignored sequences
    """
    patlen = len(cntdict.keys()[0])
    total = sum(cntdict.values())
    _from, _to = 0, patlen
    pattern_range = range(_from,_to)
    if report:
        print pattern_range
    # some variables we need
    data = []
    IC = []
    for i in range(0,patlen):
        data.append( {'A':0,'T':0,'C':0,'G':0 } )
    nns = 0 # sequences with 'N' bases
    for seq,cnt in cntdict.iteritems():
        if seq.upper().count("N") >= 1:
            nns+=cnt
            continue 
        for i in range(0,patlen):
            base = seq[i].upper()
            data[i][base]+=cnt
    for i in pattern_range:
        IC.append( {'A':0,'T':0,'C':0,'G':0 } )
        for base in list("ATGC"):
            freq = float(data[i][base])/float(total-nns)
            if freq == 0.0:
                # about to encounter a Math Log error
                freq = 1.0 / float(total-nns)
            #print base, "%2.0f" % (100.0*freq), "\t", 
            IC[-1][base] = 2.0-(-log2(freq))
    if report:
        print "%s sequences, %s ignored (N's), pattern length=%s" % (
            total,nns,len(pattern_range))
    return IC,nns

# end of function obtain_pssm_ic_from_cntdict



def parse_ic_data(x):
    """
    Parse InformationContent data into a PSSM IC list

    @type  x: string
    @param x: InformationContent data string (e.g. from print_ic() function)

    @attention: see README for InformationContent file format example

    @rtype:  list
    @return: InformationContent list
    """
    IC = []
    chars,data = x.strip().split("\n",1)
    chars = chars.strip().replace("\t"," ").split(" ")
    while '' in chars: chars.remove('')
    for line in data.split("\n"):
        scores = line.strip().replace("\t"," ").split(" ")
        while '' in scores: scores.remove('')
        empty_charset = {}
        for char in chars: empty_charset[char] = 0
        IC.append(empty_charset)
        for i in range(1,len(scores)):
            IC[-1][chars[i-1]] = float(scores[i])
    # and return
    return IC

# end of function parse_ic_data


def parse_ic_file(fname):
    """
    Parse InformationContent data file into a PSSM IC list

    @type  fname: string
    @param fname: absolute filepath to an InformationContent file

    @attention: see README for InformationContent file format example

    @rtype:  list
    @return: InformationContent list
    """
    content = []
    for line in open(fname).readlines():
        if line and line[0] != "#":
            content.append( line )
    return parse_ic_data("".join(content))

# end of function parse_ic_file


def ic2sfmpattern(ic,maxval=1000):
    """ Write InformationContent asa ScanForMatches quantitative pattern """
    sfmpattern = []
    for row in ic:
        sfmrow = []
        for base in ['A','C','G','T']:
            score = row[base]
            # calulate back to a frequency
            freq  = pow( 10, ( (score-2.0) * log10(2)) )
            # append as integer value to this sfmrow
            sfmrow.append( int( freq * maxval ) )
        sfmpattern.append( tuple(sfmrow) )
    return sfmpattern

# end of function ic2sfmpattern


class Pssm:
    """ """
    def __init__(self,ic=[],fname="",ignore_unambiguity=False,relativescore=False):
        """
        """
        if ic:
            self.ic     = ic
        elif fname:
            self.ic     = parse_ic_file(fname)
        else:
            # not properly initialized
            pass
        # initialize other arguments
        self.length     = len(self.ic)
        self.positions  = range(0,len(self.ic))
        self.ignore_unambiguity = ignore_unambiguity
        self.relativescore = relativescore 
        self._default_max = self.max(ignore_unambiguity=self.ignore_unambiguity,relativescore=False)
        self._default_min = self.min(ignore_unambiguity=self.ignore_unambiguity,relativescore=False)
        self._default_opt = self._default_max - self._default_min
    # end of function __init__

    def score(self,pattern,ignore_unambiguity=None,relativescore=None):
        """
        """
        if len(pattern) != self.length:
            raise InproperlySizedPattern
        if ignore_unambiguity == None: ignore_unambiguity = self.ignore_unambiguity
        if relativescore == None: relativescore = self.relativescore
        value = pssmscore(pattern,self.ic,ignore_unambiguity=ignore_unambiguity)
        if relativescore:
            return round( ( value - self._default_min ) / self._default_opt , 2 )
        else:
            return value

    # end of function score

    def ic2sfmpattern(self,maxval=1000):
        """ """
        return ic2sfmpattern(self.ic,maxval=maxval)

    # end of function ic2sfmpattern

    def max(self,ignore_unambiguity=None,relativescore=None):
        """ """
        if ignore_unambiguity == None: ignore_unambiguity = self.ignore_unambiguity
        if relativescore == None: relativescore = self.relativescore
        if relativescore:
            return 1.0
        else:
            return get_ic_max_score(self.ic,ignore_unambiguity=ignore_unambiguity)
    
    # end of function max
    
    def min(self,ignore_unambiguity=None,relativescore=None):
        """ """
        if ignore_unambiguity == None: ignore_unambiguity = self.ignore_unambiguity
        if relativescore == None: relativescore = self.relativescore
        if relativescore:
            return 0.0
        else:
            return get_ic_min_score(self.ic,ignore_unambiguity=ignore_unambiguity)
    
    # end of function min

    def pos_max(self,pos):
        """ Return highest score of current position in PSSM pattern """
        try:
            return max(self.ic[pos].values())
        except IndexError:
            raise InproperlySizedPattern
        except:
            raise "Unexpected Error in function PSSM.pos_highest()"

    # end of function pos_max

    def pos_min(self,pos):
        """ Return lowest score of current position in PSSM pattern """
        try:
            return min(self.ic[pos].values())
        except IndexError:
            raise InproperlySizedPattern
        except:
            raise "Unexpected Error in function PSSM.pos_lowest()"

    # end of function pos_min


    def get_strict_consensus(self):
        """ Return highest scoring dna sequence that matches to this PSSM """
        return "".join([ self.pos_max_base(pos) for pos in range(len(self.ic)) ])

    # end of function get_strict_consensus


    def abs_difference_with_other_pssm(self,other,ignore_unambiguity=None):
        """ Calculate the difference between 2 **IDENTICALLY SIZED** PSSMs """
        score = 0.0
        for pos in range(len(self.ic)):
            if (ignore_unambiguity or (ignore_unambiguity==None and self.ignore_unambiguity)) and\
            max(self.ic[pos].values()) == 2.0:
                continue
            for base in set(self.ic[pos].keys()).intersection(other.ic[pos].keys()):
                score+=abs(self.ic[pos][base]-other.ic[pos][base])
        return score

    def rel_difference_with_other_pssm(self,other,ignore_unambiguity=None):
        """ Calculate the difference between 2 **IDENTICALLY SIZED** PSSMs; averaged per columns taken into account """
        score = 0.0
        column_cnt = 0
        for pos in range(len(self.ic)):
            if (ignore_unambiguity or (ignore_unambiguity==None and self.ignore_unambiguity)) and\
            max(self.ic[pos].values()) == 2.0:
                continue
            column_cnt+=1
            for base in set(self.ic[pos].keys()).intersection(other.ic[pos].keys()):
                score+=abs(self.ic[pos][base]-other.ic[pos][base])
        return score/column_cnt


    def pos_max_base(self,pos):
        """ Return highest scoring base of current position in PSSM pattern """
        try:
            value = self.pos_max(pos)
            for k,v in self.ic[pos].iteritems():
                if v == value:
                    return k
        except IndexError:
            raise InproperlySizedPattern
        except:
            raise "Unexpected Error in function PSSM.pos_highest()"

    # end of function pos_max


    def pos_min_base(self,pos):
        """ Return lowest scoring base of current position in PSSM pattern """
        try:
            value = self.pos_min(pos)
            for k,v in self.ic[pos].iteritems():
                if v == value:
                    return k
        except IndexError:
            raise InproperlySizedPattern
        except:
            raise "Unexpected Error in function PSSM.pos_lowest()"

    # end of function pos_min

 
    def scan(self,sequence,min_score=None):
        """ """
        if min_score == None: min_score = self.min()
        pass

    # end of function scan


    def reversecomplement(self):
        """ """
        self.ic.reverse()
        for row in self.ic:
            row['A'], row['T'] = row['T'], row['A']
            row['G'], row['C'] = row['C'], row['G']

    # end of function reversecomplement

# end of class Pssm

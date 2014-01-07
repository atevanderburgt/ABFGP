"""
Code that can score the repetitiveness / low complexity of sequences
"""

# Module metadata
__authors__ = "Ate van der Burgt"
__license__ = "MIT"

# Python Imports
from sets import Set
from math import log as mathlog

# Exceptions
class NoAlphabetAppliedException(Exception):
    pass

def sequencerepetitiveness(seq,stringlengths=[1],alphabet=None):
    """
    """
    # initialize proper alphabet of aloued characters
    if type(alphabet) == type(Set()):
        alphabet = Set(alphabet)
    elif type(alphabet) == type(list()):
        alphabet = Set(alphabet)
    elif type(alphabet) == type(str()):
        alphabet = Set(list(alphabet.upper()))
    else:
        raise NoAlphabetAppliedException

    counts = {}
    seq = seq.upper()
    for n in stringlengths:
        counts[n] = {}
        for i in range(0,len(seq)-n):
            seqpart = seq[i:i+n]
            if Set(seqpart).difference(alphabet):
                continue
            if counts[n].has_key(seqpart):
                counts[n][seqpart]+=1
            else:
                counts[n][seqpart] =1
    # return counts data dictionary
    return counts

# end of function sequencerepetitiveness


def listproduct(ll):
    """ """
    if not ll: return 0
    val = ll[0]
    for item in ll[1:]: val = val*item
    return val

# end of function listproduct

def proteinsequencerepetitiveness(protseq,stringlengths=[1,2]):
    """
    """
    alphabet=Set(list("ARNDCEQGHILKMFPSTWYV"))
    counts = sequencerepetitiveness(protseq,stringlengths=stringlengths,alphabet=alphabet)
    for n in stringlengths:
        elements = sum(counts[n].values())
        items    = len(counts[n])
        maxitems = pow(len(alphabet),n)
        ratios   = [ float(v)/maxitems for v in counts[n].values() ]
        score    = listproduct(counts[n].values())
        expected = elements/float(min([elements,maxitems]))
        values   = counts[n].values()
        singleexp=elements/float(20)
        for i in range(len(values)-1,-1,-1):
            if values[i] < singleexp: values.pop(i)
        scoreval = listproduct(values)
        print n, score, scoreval, len(counts[n]), elements, maxitems, expected, len(values), items


# end of function proteinsequencerepetitiveness


def CalculateDiNucleotideComplexityScore(seqstring):
    """
    return fraction of two most occuring dinucleotide-combination
    """
    hpLength = len(seqstring)
    seqstring = seqstring.lower()
    dinucs = [ seqstring[i-1:i+1] for i in range (1,hpLength)]
    #dinucs_count = {}
    #for dinuc in Set(dinucs):
    #    dinucs_count[dinuc] = dinucs.count(dinuc)
    #return dinucs_count
    counts = []
    for dinuc in Set(dinucs):
        counts.append(dinucs.count(dinuc))
    counts.sort()
    return round( float(sum(counts[-2:]))/float(hpLength-1), 3 )
            
# end of CalculateDiNucleotideComplexityScore


def proteinsequencelowcomplexityscore(protseq):
    """
    """
    seqlen = len(protseq)
    protseq = protseq.upper()
    dinucs = [ protseq[i-1:i+1] for i in range (1,seqlen)]
    counts = [ 0 ]
    expected = max([seqlen/400,1])
    cutoff   = int(round(mathlog(seqlen,5)))
    for dinuc in Set(dinucs):
        occurrence = dinucs.count(dinuc) - expected
        if occurrence >= cutoff: 
            counts.append( occurrence - cutoff )
    return round( sum(counts) / float(seqlen), 3)

# end of function proteinsequencelowcomplexityscore


if __name__ == "__main__":
    # print some test data
    pass 

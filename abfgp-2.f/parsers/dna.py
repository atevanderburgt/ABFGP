""" functions on DNA sequence strings """

def reversecomplement(seq,seqtype='dna'):
    """
    """
    rev = []
    trans = {'a':'t','A':'T','t':'a','T':'A',
             'g':'c','G':'C','c':'g','C':'G',
             'u':'a','U':'A',
             'n':'n','N':'N',
             '-':'-','~':'~' }

    for base in list(seq):
        if trans.has_key(base):
            rev.append(trans[base])
        else:
            rev.append(base)
    # and reverse the sequence!
    rev.reverse()
    rev = "".join(rev)
    if seqtype.lower() == 'rna':
        rev = rev.replace("t","u").replace("T","U")
    # and return
    return rev

# end of function reversecomplement


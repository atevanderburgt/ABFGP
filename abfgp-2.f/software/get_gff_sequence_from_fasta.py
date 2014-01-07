"""
"""
import os, sys, gzip

##########################################################################

def parseFasta(content):
    """
    parses a array representing fasta (file)
    into a dict with headers and seqs
    content can be: sys.stdin.readlines()
    content can be: fh.readlines()
    """
    seqs = {}; name = ''; seqb = [];
    for line in content:
        line = line.strip()
        if not line: continue
        if line[0] == '>':
            #new fasta seq, store if there was a previous one
            if name and seqb: 
                seqs[name] = "".join(seqb)
            name = line[1:]; seqb = []
        else: seqb.append(line)
    if name and seqb: seqs[name] = "".join(seqb)
    return seqs

def _reversecomplement(seq,seqtype='dna'):
    """
    """
    rev = []
    trans = {'a':'t','A':'T','t':'a','T':'A',
             'g':'c','G':'C','c':'g','C':'G',
             'u':'a','U':'A',
             'n':'n','N':'N' }
    
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


##########################################################################


# read gff lines from STDIN
gfflines  = sys.stdin.readlines()
# obtain fasta filename from command line
fastafile = sys.argv[1]

_params = sys.argv[2:]
params = {}
for offset in range(0,len(_params),2):
    params[_params[offset]] = _params[offset+1]

# read fasta sequences
if fastafile[-3:] == ".gz":
    seqs = parseFasta(gzip.open(fastafile).readlines())
else:
    seqs = parseFasta(open(fastafile).readlines())

# all fasta headers
fastaheaders = seqs.keys()

def _freflookup(fref,headers,complexlookup=False):
    """ """
    # first simple string lookup
    if fref in headers: return fref
    # if no succeed: more complex lookup 
    seen = 0
    thefref = None
    for _fref in headers:
        if _fref.find(fref) >= 0 and (_fref.find(fref)+len(fref)==len(_fref) or _fref[_fref.find(fref)+len(fref)] == ' '):
            thefref = _fref
            seen+=1
            if seen>=2: break
    # check if only a single one seen
    if seen == 1: return thefref
    # if no succeed: even more complex lookup
    if not complexlookup and "_" in fref:
        return _freflookup(fref.replace("_"," "),headers,complexlookup=True)
    elif not complexlookup and " " in fref:
        return _freflookup(fref.replace(" ","_"),headers,complexlookup=True)
    else:
        pass
    # if here, then we definately cannot find the fref in the fasta headers 
    message = "unrecognized header: '%s' (format='%s',seen=%s)" % (fref, headers[0],seen)
    raise Exception, message

# end of function _freflookup 

for line in gfflines:
    if not line.strip(): continue
    gff = line.split("\t")
    # obtain header from gffline and do _freflookup
    fref = gff[0]
    thefref = _freflookup(fref,fastaheaders)
    # obtain coordinates and strand
    fstart,fstop,fstrand = int(gff[3]), int(gff[4]), gff[6]
    if params.has_key("-fstartextra"):
        fstart -= int(params["-fstartextra"])
    if params.has_key("-fstopextra"):
        fstop += int(params["-fstopextra"])
    seq = seqs[thefref][fstart-1:fstop]
    if fstrand == '-':
        seq = _reversecomplement(seq)
    # and print the header/sequence thingy
    print ">%s %s|%s|%s" % (fref,fstart,fstop,fstrand)
    print seq

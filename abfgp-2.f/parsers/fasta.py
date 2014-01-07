""" Functions for parsing,accessing and writing FASTA files """

# Module metadata
__authors__ = "Ate van der Burgt"
__license__ = "MIT"

# Python imports
import gzip
from subprocess import Popen, PIPE
from os.path import (
    join as osPathJoin,
    abspath as osPathAbspath,
    dirname as osPathDirname,
    )

# absolute paths to executables located within the same directory
EXECUTABLE_ISSINGLEFASTADNA     = osPathJoin(osPathDirname(osPathAbspath(__file__)),"issinglefastadna.sh")
EXECUTABLE_ISSINGLEFASTAPROTEIN = osPathJoin(osPathDirname(osPathAbspath(__file__)),"issinglefastaprotein.sh")

########################################################################
#### Exceptions
########################################################################

class NoSingleFastaFile(Exception):
    pass

########################################################################
#### Validation functions
########################################################################

def IsSingleFastaDna(fname,executable=EXECUTABLE_ISSINGLEFASTADNA):
    """
    Is this file a SINGLE DNA SEQUENCE fasta file?       

    @type  fname: string
    @param fname: (full) path to fasta file

    @type    executable: string
    @keyword executable: (full) path to issinglefastadna.sh executable (required)

    @rtype:  integer
    @return: 1 (True) or 0 (False)
    """
    return int(Popen("%s %s" % (executable,fname), shell=True, stdout=PIPE, close_fds=True).stdout.read().strip())
    #return int(ospopen("%s %s" % (executable,fname)).read())

# end of function IsSingleFastaDna


def IsSingleFastaProtein(fname,executable=EXECUTABLE_ISSINGLEFASTAPROTEIN):
    """
    Is this file a SINGLE PROTEIN SEQUENCE fasta file?
    
    @type  fname: string
    @param fname: (full) path to fasta file

    @type    executable: string
    @keyword executable: (full) path to issinglefastaprotein.sh executable (required)

    @rtype:  integer
    @return: 1 (True) or 0 (False)
    """
    return int(Popen("%s %s" % (executable,fname), shell=True, stdout=PIPE, close_fds=True).stdout.read().strip())
    #return int(ospopen("%s %s" % (executable,fname)).read())

# end of function IsSingleFastaProtein


########################################################################
#### Helper functions
########################################################################

def _validate_header(header):
    """
    validates and sanitizes a provided fasta accession header prior to writing to file

    @type  header: string
    @param header: fasta accession entry

    @rtype:  string
    @return: sanitized fasta accession entry
    
    @attention: removes potential present '>' characters
    @attention: removes potential present leading and tailing whitespace characters
    """
    return header.strip().strip(">").strip()

# end of function _validate_header
    
def _parseFastaHeader(header):
    """
    parse a fasta header line and split its on its non-spaced accession and description

    @type  header: string
    @param header: fasta accession entry
    
    @rtype:  ( string, None or string )
    @return: ( sanitized fasta accession entry, None or the fasta header description if available )
    """
    header = _validate_header(header)
    description = None
    for char in (" ","\t"):
        if header.find(char) >= 0:
            header, _description = header.split(char,1)
            if not description:
                description = _description
            else:
                # parsed description exists; concatenate 
                description = _description + char + description
    # return header & decription
    return header, description    

# end of function _parseFastaHeader
    
########################################################################
#### Main functions
########################################################################
    
def parseFasta(content,allow_decorated_headers=True):
    """
    Parse an iterable representing a fasta (file) into a dict with headers and seqs.
    
    @type  content: iterable
    @param content: iterable that represents a fasta file    

    @type    allow_decorated_headers: Boolean
    @keyword allow_decorated_headers: allow descriptions in fasta headers; removes description if False (required)   
    
    @rtype:  dict
    @return: fasta sequence dictionary
    
    @attention: content can be: sys.stdin.readlines()
    @attention: content can be: fh.readlines()
    @attention: The fasta header line(s) will be split on spaces and tabs
    """
    seqs = {}; header = ''; seqb = [];
    for line in content:
        line = line.strip()
        if not line: continue
        if line[0] == '>':
            #new fasta seq, store if there was a previous one
            if header and seqb: 
                seqs[header] = "".join(seqb)
                seqb = []
            # reinit/empty seqdb and parse new header
            seqb = []
            if allow_decorated_headers:
                header = _validate_header(line)
            else:
                header,_description = _parseFastaHeader(line)
        else: seqb.append(line)
    if header and seqb: seqs[header] = "".join(seqb)
    return seqs

# end of function parseFasta
    
def parseDecoratedFasta(content):
    """
    Parse an iterable representing a fasta (file) into a dict with headers and seqs, and additional fasta accession information

    @type  content: iterable
    @param content: iterable that represents a fasta file    

    @rtype:  ( dict, dict )
    @return: ( fasta sequence dictionary, dictionary with per fasta accession its header description )

    @attention: content can be: sys.stdin.readlines()
    @attention: content can be: fh.readlines()

    @attention: The fasta header line(s) will be split on spaces and tabs
    @attention: This fasta header meta-data will be returned as an additional dict
    """
    seqs = {}; header = ''; description = ''; headerinfo = {}; seqb = [];
    for line in content:
        line = line.strip()
        if not line: continue
        if line[0] == '>':
            #new fasta seq, store if there was a previous one
            if header and seqb:
                seqs[header] = "".join(seqb)
                seqb = []
                headerinfo[header] = description
            header,description = _parseFastaHeader(line)
        else: seqb.append(line)
    if header and seqb:
        seqs[header] = "".join(seqb)
        headerinfo[header] = description
    return seqs,headerinfo

# end of function parseDecoratedFasta

def parseSingleFasta(content):
    """ 
    Parse an iterable representing a fasta (file) with a single entry into its header, sequence and meta-description
    
    @type  content: iterable
    @param content: iterable that represents a fasta file    

    @rtype:  ( string, string, None or string )
    @return: ( fasta header, fasta sequence, None or the fasta header description if available )
    
    @attention: content can be: sys.stdin.readlines()
    @attention: content can be: fh.readlines()
    @attention: no sanity check will be performed: multi-fasta entries will be concatenated as-if sequence strings!
    @attention: The fasta header line(s) will be split on spaces and tabs
    @attention: This fasta header meta-data will be returned as an additional dict
    """
    try:
        header, description = _parseFastaHeader(content[0])
        sequence = "".join([ line.strip() for line in content[1:] ])
        return header, sequence, description
    except:
        raise NoSingleFastaFile

# end of function parseSingleFasta      
        
def parseSingleFastaHeaderFromFile(fname):
    """
    return the first occurring fasta header in a fasta file
    
    @type  fname: string
    @param fname: (full) path to existing fasta file
    
    @rtype:  ( string, string )
    @return: ( fasta header, None or the fasta header description if available )
    """
    command = """ head -n 1 %s | grep -m 1 '>' """ % fname
    p = Popen(command, shell=True, stdout=PIPE, stderr=PIPE, close_fds=True)
    error = p.stderr.read()
    if error:
        message = "ERROR: '%s'" % error
        raise IOError, message
    p.stderr.close()
    header = p.stdout.read()
    p.stdout.close()
    if not header:
        raise NoSingleFastaFile
    return _parseFastaHeader(header)
    
    ## outdated code using popen3
    #ci,co,ce = popen3("head -n 1 %s | grep -m 1 '>'" % fname)
    #ci.close()
    #error = ce.read()
    #if error:
    #    print "ERROR: '%s'" % error
    #    raise IOError
    #header = co.read()
    #if not header: raise NoSingleFastaFile
    #return _parseFastaHeader(header)

# end of function parseSingleFastaHeaderFromFile 


def writeMultiFasta(seqs,fname,linesize=None):
    """
    write a (multi) fasta file from a fasta dictionary

    @type  seqs: dictionary
    @param seqs: dictionary containing fasta accession entries
    
    @type  fname: string
    @param fname: (full) path to output fasta file
    
    @type    linesize: positive integer or None
    @keyword linesize: line width used for writing the sequence string to file.
    """
    fh = open(fname,'w')
    for h,s in seqs.iteritems():
        if linesize and type(linesize) == type(int()) and linesize > 0:
            fh.write(">%s\n" % (_validate_header(h)) )
            for offset in range(0,len(s),linesize):
                fh.write("%s\n" % s[offset:offset+linesize])
        else:
            fh.write(">%s\n%s\n" % (_validate_header(h),s) )
    fh.close() 

# end of function writeMultiFasta
    
# alias function name; no difference between single & multi fasta
def writeFasta(*args):
    writeMultiFasta(*args)

def writeSingleFasta(header,sequence,fname,linesize=None,description=""):
    """
    write a single fasta file, containing a meta description if desired

    @type  header: string
    @param header: fasta accession entry

    @type  sequence: string
    @param sequence: sequence string (dna|rna|protein)
    
    @type  fname: string
    @param fname: (full) path to output fasta file
    
    @type    description: string
    @keyword description: meta description on the fasta header line

    @type    linesize: positive integer or None
    @keyword linesize: line width used for writing the sequence string to file
    """
    fh = open(fname,'w')
    fh.write( ">%s" % _validate_header(header) )
    if description: fh.write(" %s" % description)
    if linesize and type(linesize) == type(int()) and linesize > 0:
        fh.write("\n%s\n" % "\n".join([ sequence[offset:offset+linesize] for offset in range(0,len(sequence),linesize) ]) )
    else:
        fh.write("\n%s\n" % sequence )
    fh.close()

# end of function writeSingleFasta


def getSeqFromFasta(header,fname):
    """
    obtain a sequence from a fasta file by its accession header
    
    @type  header: string
    @param header: fasta accession entry

    @type  fname: string
    @param fname: (full) path to existing fasta file
    
    @rtype:  string
    @return: fasta sequence string    
    """
    if fname[-3:] == ".gz":
        seqs = parseFasta(gzip.open(fname).readlines()) 
    else:
        seqs = parseFasta(open(fname).readlines())
    for h,seq in seqs.iteritems():
        if h == header or _validate_header(h) == header:
            return seq
    else:
        # EOF forloop reached; fasta header not found
        return False
        
# end of function getSeqFromFasta

def correctspacedfastaheaders(sequences):
    """
    corrects a fasta dictionary with descriptive fasta headers; provides a lookup dict for each adjusted accession

    @type  sequences: dictionary
    @param sequences: dictionary containing fasta accession entries, putatively with description in its headers

    @rtype:  ( dict, dict )
    @return: ( ADJUSTED fasta sequence dictionary, dictionary with lookup keys for the CHANGED ACCESSIONS ONLY )
    """
    hdr2hdrwithdescription = {}
    for hdr in sequences.keys():
        #### m = re.match("(?P<header>[^ \t]+)[ \t]{0,1}(?P<description>.*)",hdr)
        #### if m.groupdict()['description']:
        header, description = _parseFastaHeader(hdr)
        if description:
            hdr2hdrwithdescription[header] = hdr
    if hdr2hdrwithdescription:
        for newheader,fullheader in hdr2hdrwithdescription.iteritems():
            sequences[newheader] = sequences[fullheader]
            del( sequences[fullheader] )
    return sequences, hdr2hdrwithdescription

# end of function correctspacedfastaheaders


def __test(fname=None):
    """ unit test for the functions in this file """
    fastadict = {
        'entry01'                     : 'ATCGACGACTGACTAGCTACGTACGACTGATCC',
        'entry02 metadata entry02'    : 'ATCGACGACTGACTAGCTACGTACGACTGATCC',
        'entry03\tmetadata entry03'   : 'ATCGACGACTGACTAGCTACGTACGACTGATCC',
        ' entry04 incorrectly spaced' : 'ATCGACGACTGACTAGCTACGTACGACTGATCC',
        '> entry05 malformatted'      : 'ATCGACGACTGACTAGCTACGTACGACTGATCC',
        }
        
    def _layout_hdr(hdr):
        return (("'%s'" % hdr.replace("\t","        ")) + " "*40)[0:40]

    # get ordered fasta keys
    hdrs = fastadict.keys()
    hdrs.sort()
    hdrs.append( hdrs.pop(0) )
    hdrs.append( hdrs.pop(0) )
    # get iterable (pseudo) content
    content = []
    for hdr in hdrs: content.extend([(">"+hdr).replace(">>",">"),fastadict[hdr]])
    print "# testing: _validate_header(header)" 
    for hdr in hdrs:
        print "    %s -> '%s'" % (_layout_hdr(hdr),_validate_header(hdr))
    print "# testing: _parseFastaHeader(header)" 
    for hdr in hdrs:
        _hdr,etc = _parseFastaHeader(hdr)
        print "    %s -> '%s' '%s'" % (_layout_hdr(hdr),_hdr,etc)
    print "# testing: parseFasta(content,allow_decorated_headers=True)"
    print "    %s" % str(tuple(parseFasta(content,allow_decorated_headers=True).keys()))
    print "# testing: parseFasta(content,allow_decorated_headers=False)"
    print "    %s" % str(tuple(parseFasta(content,allow_decorated_headers=False).keys()))
    print "# testing: parseDecoratedFasta(content)"
    _seqs,_info = parseDecoratedFasta(content)
    print "    %s" % str(tuple(_seqs.keys()))
    print "    %s" % str(_info)
    print "# testing: correctspacedfastaheaders(fastadict)"
    print "    IN    : %s" % str(tuple(fastadict.keys()))
    fastadict, hdr2hdrwithdescription = correctspacedfastaheaders(fastadict)
    print "    OUT   : %s" % str(tuple(fastadict.keys()))
    print "    IN2OUT: %s" % str(hdr2hdrwithdescription)
    if fname:
        print "# testing: IsSingleFastaDna(fname)"
        print "    %s" % IsSingleFastaDna(fname)
        print "# testing: IsSingleFastaProtein(fname)"
        print "    %s" % IsSingleFastaProtein(fname)
        print "# testing: parseFasta(content,allow_decorated_headers=True)"
        _seqs = parseFasta(open(fname),allow_decorated_headers=True)
        print "    %s..." % str(tuple(_seqs.keys()[0:5]))
        print "# testing: parseFasta(content,allow_decorated_headers=False)"
        _seqs = parseFasta(open(fname),allow_decorated_headers=False)
        print "    %s..." % str(tuple(_seqs.keys()[0:5]))
        print "# testing: parseDecoratedFasta(content)"
        _seqs,_info = parseDecoratedFasta(open(fname))
        print "    %s..." % str(tuple(_seqs.keys()[0:5]))
        print "    %s..." % str(dict([(h,k) for h,k in _info.items()[0:5] ]))
    
# end of function __test

if __name__ == "__main__":
    import sys
    from os.path import isfile as osPathIsfile
    print """ Functions for parsing,accessing and writing FASTA files """
    print """ Unit testing: [ %s test ] [ %s test SINGLE_FASTA_FNAME ] """ % (sys.argv[0],sys.argv[0])
    if not 'test' in sys.argv: sys.exit()
    fname = None
    if len(sys.argv) == 3:
        fname = sys.argv[2]
        if not osPathIsfile(fname):
            fname = None
    # test all the functions
    __test(fname=fname)

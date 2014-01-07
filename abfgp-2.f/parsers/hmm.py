"""
HMM parser functions
Developed & tested for HMMER version 2.3.2
Currently only a small subset of protein HMMER functionality is covered.
"""

# Module metadata
__authors__ = "Ate van der Burgt"
__license__ = "MIT"

# Python Imports
import os

# Global variables / program paths
#EXECUTABLE_HMMPATH      = "/home/avdb/data/abfgp-dev/software/hmmer-2.3.2"
#EXECUTABLE_HMMSEARCH    = os.path.join(EXECUTABLE_HMMPATH,"hmmsearch")
#EXECUTABLE_HMMBUILD     = os.path.join(EXECUTABLE_HMMPATH,"hmmbuild")
from executables import (
    EXECUTABLE_HMMPATH,
    EXECUTABLE_HMMSEARCH,
    EXECUTABLE_HMMBUILD,
    )


def hmmbuild_protein(fname_fasta,fname_hmmbuild=""):
    """
    Run hmmbuild on a file with aligned protein sequences

    @type  fname_fasta: string
    @param fname_fasta: (absolute) file path of aligned protein sequence file
    
    @type  fname_hmmbuild: string
    @param fname_hmmbuild: specify desired name and path of the hmmbuild file

    @rtype:  string
    @return: (absolute) file path of the hmmbuild file
    """
    if not fname_hmmbuild:
        # get filename for the HMM build file
        if fname_fasta[-3:] == '.fa':
            fname_hmmbuild = fname_fasta[0:-3]+".hmmbuild"
        else:
            fname_hmmbuild = fname_fasta+".hmmbuild"

    # first remove hmmbuild file if it happens to be there
    try: os.remove(fname_hmmbuild)
    except OSError: pass 

    # run hmmbuild; -F for force-overwrite of output file, --amino for AA sequences
    command = "%s -F --amino %s %s" % (EXECUTABLE_HMMBUILD,fname_hmmbuild,fname_fasta)  
    ci,co,ce = os.popen3(command)
    ci.close()
    out = co.read()
    co.close()
    error = ce.read().strip()
    ce.close()
    if error:
        # messages are casted to STDER. Print them here to STDOUT
        print "HMMbuild ERROR:", command 
        print "HMMbuild ERROR:", error
        # set fname_hmmbuild to None
        fname_hmmbuild = None

    # return name of the created file
    return fname_hmmbuild

# end of function hmmbuild_protein


def hmmsearch_protein(fname_hmmbuild, fname_database,
    verbose = False,
    params= {'E':150,'A':5} ):
    """
    Run hmmsearch on a database with sequences with a hmmbuild pattern and
    return the best hits

    @type  fname_hmmbuild: string
    @param fname_hmmbuild: (absolute) file location of the hmmbuild file

    @type  fname_database: string
    @param fname_database: (absolute) file location of the (multifasta) database
                            to scan

    @type  verbose: Boolean
    @param verbose: print debugging messages to STDOUT (True) or not (False, default)

    @type  params: dict
    @param params: dictionary with alowed hmmsearch params (A,E,T,Z)

    @rtype:  list
    @return: list with hmmsearch hits
    """
    # each hmm hit spans 5 lines; add 5 more for some whitespace etc
    # BUT, when the alignment does not fit onto a single line...
    # depending on header length (?), ~47 characters are placed
    # Therefor, set spatiously large grep_offset
    grep_offset = 100000

    # hmmsearch command and some output pre-processing
    command = "%s %s %s %s | grep -A %s \"%s\" | grep -B %s \"%s\" | sed 1d" % (
            EXECUTABLE_HMMSEARCH,
            " ".join([ "-%s %s" % (k,v) for k,v in params.iteritems() ] ),
            fname_hmmbuild,
            fname_database,
            grep_offset,
            "^Alignments of top-scoring domains:$",
            grep_offset,
            "^Histogram of all scores:$"
            ) 


    # run the command
    ci,co = os.popen2(command)
    ci.close()
    output = "".join( co.readlines()[0:-1] ).strip()
    co.close()

    if output.find("no hits above thresholds") >= 0:
        # no hits found at all!!
        return []

    # parse the output
    results = []
    parts = output.split("\n\n")
    searchphrase = "[output cut off at A = %s" % params['A']
    if params.has_key('A') and parts[-1].find(searchphrase) >= 0:
        # output is trimmed at A = x best hits
        parts = parts[0:-1]


    check = [ len(part.split("\n")) for part in parts ]
    if 3 not in check:
        # just short, simple outputs with each hit sequence on a single line
        pass
    else:
        # Long hmmbuild causing the results to overrun the line lengths. 
        # Repair this by merging back on individual results
        parts = _convert_hmm_multiline_to_singleline(parts)

    for part in parts:
        # convert each hmmer protein hit to a tuple with data
        header,query_line,match_line,sbjct_line = part.split("\n") 
        ########################################################################
        if verbose:
            print header
            print query_line
            print match_line
            print sbjct_line
        ########################################################################

        # convert each hmmer protein hit to a tuple with data
        sbjct_header = header.split(':')[0]
        offsetA, offsetB = query_line.find('*->')+3, query_line.find('<-*')
        query = query_line[offsetA:offsetB]
        match = match_line[offsetA:offsetB]
        sbjct = sbjct_line[offsetA:offsetB]
        sbjct_start = int(sbjct_line[0:offsetA].strip().split(' ')[-1])
        sbjct_end   = int(sbjct_line[offsetB:].strip().split(' ')[0])
        query_start = 1
        query_end   = len(query) - query.count(".")
        score, expect = header.split(":")[-1].split(",")
        score = float( score.replace("score","").strip() )
        expect= float( expect.replace("E =","").strip() )
        # strip leading gaps from the sbjct alignment/sequence
        while sbjct[0] == '-':
            query = query[1:]
            match = match[1:]
            sbjct = sbjct[1:]
            query_start += 1
        # strip trailing gaps from the sbjct alignment/sequence
        while sbjct[-1] == '-':
            query = query[0:-1]
            match = match[0:-1]
            sbjct = sbjct[0:-1]
            query_end -= 1
        # append this hmmsearch protein hit to the list of results
        results.append( ( sbjct_header, sbjct_start, sbjct_end,
                          query_start, query_end,
                          query, match, sbjct, score, expect ) )

        #if sbjct_header == 'dotse_orf_98':
        #    print "#"*50
        #    print open(tmp_output_file).read() 
        #    print "#"*50
        #    import sys
        #    sys.exit()

    #### remove hmmsearch output file
    ###try: os.remove(tmp_output_file)
    ###except OSError: pass
       
    # return the results
    return results

# end of function hmmsearch_protein


def _convert_hmm_multiline_to_singleline(parts):
    """
    @type  parts: list
    @param parts: list with hmmer protein output blocks

    @attention: this function is only meaningfull as a sub in hmmsearch_protein

    @rtype  mergedparts: list
    @return mergedparts: list with hmmsearch hits joined in single-lined output
    """
    # list with line counts of parts; 4 == simple block, 3 means multilined
    # block that belongs to an above simple block (with 4 lines)
    check = [ len(part.split("\n")) for part in parts ]

    # Long hmmbuild causing the results to overrun the line lengths. 
    # Repair this by merging back on individual results
    mergelist = []
    mergedparts = []
    for i in range(0,len(check)):
        if check[i] == 4:
            mergelist.append([i])
        else:
            mergelist[-1].append(i)
            
    # mergelist now contains lists of part-id's
    # that belong to the same result block
    for pieces in mergelist:
        if len(pieces) == 1:
            # some outputs span a line, but this piece does not!
            mergedparts.append( parts[pieces[0]] )
            continue

        # start constructing single-lined output
        header,q,m,s = parts[pieces[0]].split("\n")
        while q[-1] == ' ': q = q[0:-1]
        offsetA = q.find('*->')
        offsetB = len(q) 
        reconstructedQ = q[0:offsetB]
        reconstructedM = m[0:offsetB]
        reconstructedS = s[0:offsetB]
        sbjctEndPos    = s[offsetA:].strip().split(' ')[-1]
        sbjctStartPos  = s[0:offsetA].strip().split(' ')[-1]
        for pos in pieces[1:]:
            q,m,s = parts[pos].split("\n")
            offsetB = len(q)
            reconstructedQ = reconstructedQ + q[offsetA:offsetB]
            reconstructedM = reconstructedM + m[offsetA:offsetB]
            reconstructedS = reconstructedS + s[offsetA:offsetB]
            if s[offsetB:].strip():
                # save sbjctEndPos in case last part is strangely truncated
                #endPos = s[offsetB:].strip()
                endPos = s[offsetA:].strip().split(' ')[-1]
                try:    sbjctEndPos = str(int(endPos))
                except: pass
            if sbjctStartPos == '-':
                # hmmmatch started with no sbjct coordinates (only ----)
                sbjctStartPos  = s[0:offsetA].strip().split(' ')[-1]

        # add sbjctEndPos to the Sbjct line
        reconstructedS = reconstructedS+"    "+sbjctEndPos
        if reconstructedS.strip()[-1] == '-':
            ### very exceptional case; break was through the <-* symbol, causing
            ### the existance of a nearly empty sbjct line:
            ###                    -*
            ###
            ###   ncu_orf_49     -    -   
            ### the absent coordinate is rescued in sbjctEndPos
            dashpos = reconstructedS.rfind('-')
            reconstructedS = reconstructedS[0:dashpos] + sbjctEndPos  
        if reconstructedS[0:offsetA].strip().split(' ')[-1] == '-':
            # very exceptional case; first sbjct coord wash a dash i.s.o an integer
            # that happens when there is not a single sbjct AA on the first line
            dashpos = reconstructedS[0:offsetA].rfind('-')
            intlen = len(sbjctStartPos)
            posA = reconstructedS[0:offsetA].rfind('-')+1-intlen 
            reconstructedS = reconstructedS[0:posA]+sbjctStartPos+' '+reconstructedS[offsetA:]

        # multi-lined part(s) reconstructed to single lined part
        mergedparts.append( "\n".join([header,reconstructedQ,reconstructedM,reconstructedS]) )

    # return mergedparts
    return mergedparts
 
# end of function _convert_hmm_multiline_to_singleline


class HmmsearchProteinHit:
    """ Class describing a hmmsearch obtained protein hit """ 
    def __init__(self):
        self.sbjct_header   = ""    # string
        self.sbjct_start    = None  # integer gte 0
        self.sbjct_end      = None  # integer gte 0
        self.query_start    = None  # integer gte 0
        self.query_end      = None  # integer gte 0
        self.query          = ""    # string
        self.match          = ""    # string
        self.sbjct          = ""    # string
        self.score          = 0     # integer
        self.expect         = 0.0   # float

    # end of function __init__


    def trim_unmatched_characters(self):
        """ Remove unmatched AA characters from both ends of the hit """
        while self.match and self.match[0] == ' ':
            self.query = self.query[1:]
            self.match = self.match[1:]
            self.sbjct = self.sbjct[1:]
            self.sbjct_start+=1
            self.query_start+=1
        while self.match and self.match[-1] == ' ':
            self.query = self.query[0:-1]
            self.match = self.match[0:-1]
            self.sbjct = self.sbjct[0:-1]
            self.sbjct_end-=1
            self.query_end-=1

    # end of function trim_unmatched_characters


    def trim_exterior_gaps(self):
        """ Remove leading and trailing gaps from the hit """
        while self.sbjct[0] == '-':
            self.query = self.query[1:]
            self.match = self.match[1:]
            self.sbjct = self.sbjct[1:]
            self.query_start += 1
        while self.sbjct[-1] == '-':
            self.query = self.query[0:-1]
            self.match = self.match[0:-1]
            self.sbjct = self.sbjct[0:-1]
            self.query_end -= 1

    # end of function trim_exterior_gaps


    def init_from_hmmsearch_tuple(self,hmmhit):
        """ Initialize objects values from a hmmsearch hit tuple """
        sbjct_header, sbjct_start, sbjct_end, query_start, query_end, query, match, sbjct, score, expect = hmmhit
        self.sbjct_header   = hmmhit[0]
        self.sbjct_start    = hmmhit[1]
        self.sbjct_end      = hmmhit[2]
        self.query_start    = hmmhit[3]
        self.query_end      = hmmhit[4]
        self.query          = hmmhit[5]
        self.match          = hmmhit[6]
        self.sbjct          = hmmhit[7]
        self.score          = hmmhit[8]
        self.expect         = hmmhit[9]

    # end of function init_from_hmmsearch_tuple


    def as_tuple(self):
        """ return the values in the hmmsearch protein hit as a tuple """
        return ( self.sbjct_header, self.sbjct_start, self.sbjct_end,
                 self.query_start, self.query_end,
                 self.query, self.match, self.sbjct,
                 self.score, self.expect )
    
    # end of function as_tuple
    
# end of class HmmsearchProteinHit

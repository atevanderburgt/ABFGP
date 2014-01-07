"""
CLUSTALW parser functions
Developed & tested for CLUSTALW version 1.83 
"""

# Module metadata
__authors__ = "Ate van der Burgt"
__license__ = "MIT"

# Imports
from os import popen2 as osPopen2, remove as osRemove
from sets import Set
from fasta import writeMultiFasta
from proteinsimilaritymatrix import ProteinSimilarityMatrix
from pythonlibs.uniqueness import get_random_string_tag

# Global variables
#EXECUTABLE_CLUSTALW = "/opt/ab/bin/clustalw"
from executables import EXECUTABLE_CLUSTALW

def clustalw(inputfile="",seqs={},remove_inputfile=True,params={}):
    """
    """
    if inputfile and seqs:
        raise "wrong usage!"
    elif inputfile and not seqs:
        # input is (hopefully) a filename
        pass
    elif not inputfile and seqs:
        # input is (hopefully) sequences
        # do a quick check if (sequence) strings are given
        ARE_ALL_STRINGS = True
        for header, seq in seqs.iteritems():
            if not seq:
                ARE_ALL_STRINGS = False
                break
        if not ARE_ALL_STRINGS:
            raise Exception, "no sequence string(s) specified: %s" % seqs
        # make a kind of semi-unique filename
        uniqueid = get_random_string_tag()
        inputfile = uniqueid+"_"+"_".join([ _nonstringheader2stringheader(hdr) for hdr in seqs.keys()[0:5] ])
        inputfile+=".mfa"
        writeMultiFasta(seqs,inputfile)
    else:
        # no input at all
        raise "no input specified"

    # okay, do the clustalw
    fname_in = inputfile
    # get hard-assigned parameters
    paramstring = " ".join([ "-%s=%s" % (k,v) for k,v in params.iteritems() ]) 
    ci,co = osPopen2("%s %s %s" % (EXECUTABLE_CLUSTALW,fname_in, paramstring))
    ci.close()
    clwout = co.read()
    co.close()
    # abstract output filenames from input filename
    if fname_in.find(".") == -1:
        fname_out  = fname_in+".aln"
        fname_tree = fname_in+".dnd"
    else:
        _base      = fname_in[0:fname_in.rfind(".")]
        fname_out  = _base+".aln"
        fname_tree = _base+".dnd"

    # parse alignment output file
    _seqs,_alignment = _parse_clustalw(fname_out)
    # and delete tmp. created files
    osRemove(fname_out)
    osRemove(fname_tree)
    if remove_inputfile: osRemove(fname_in)
    # check if the keys (headers) in _seqs correspont to those in seqs
    # differences can occur when non-string headers are used

    # and return
    return (_seqs,_alignment)

# end of function clustalw


def _nonstringheader2stringheader(hdr):
    """
    """
    if not hdr:
        raise "No fasta header specified!"
    elif type(hdr) == type(str()):
        return hdr.split(" ")[0]
    elif type(hdr) == type(tuple()):
        return str(hdr[0]) 
    elif type(hdr) == type([]):
        return str(hdr[0])
    else:
        return str(hdr).split(" ")[0]

# end of function _nonstringheader2stringheader


def _parse_clustalw(fname):
    """
    parse clustalw output file
    """
    seqs = {}
    alignmentstring = ""

    blocks = "".join(open(fname).readlines()[3:]).split("\n\n")
    first = blocks[0].split("\n")[0:-1]
    # get the headers and the offset (whitespace before sequence alignment) 
    for line in first:
        header,etc = line.split(" ",1)
        if header: seqs[header] = ""
    # get the offset (whitespace before sequence alignment)
    offset = len(" ".join(first[0].split(" ")[0:-1]))+1
    alignmwidth = len(first[0].split(" ")[-1])
    # and now parse all blocks!
    for block in blocks:
        lines = block.strip().split("\n")
        if not lines[-1]:
            lines.pop()
        # alignment of sequences
        alignmwidth = 0
        headerstosee = Set(seqs.keys())
        for i in range(0,len(seqs)):
        #for line in lines[0:-1]:
            line = lines[i]
            splitted = line.split(" ")
            header = splitted[0]
            headerstosee.remove(header)
            alignm = splitted[-1]
            alignmwidth = len(alignm)
            seqs[header]+=alignm
        if headerstosee:
            # hmmm spaces not present on alignment score string
            for header in headerstosee:
                seqs[header] += "-"*alignmwidth
            # and fix the alignment string
            alignmentstring += " "*alignmwidth
        elif len(lines) == len(seqs):
            # not a single aligned character on this line
            # append to alignmentstring
            alignmentstring += " "*alignmwidth 
        else:
            # nothing weird here
            # alignment score string
            part = lines[-1][offset:]
            if len(part) < alignmwidth:
                part += " "*(alignmwidth-len(part))
            # and append to alignmentstring
            alignmentstring += part
    return seqs,alignmentstring

# end of function _parse_clustalw


def strip_alignment_for_exterior_gaps(alignedseqs,alignment,coords):
    """
    @type  alignedseqs: dict
    @param alignedseqs: dict ASIS produced by clustalw.clustalw() function

    @type  alignment: string
    @param alignment: string ASIS produced by clustalw.clustalw() function

    @type  coords: dict
    @param coords: dict with identical keys as alignedseqs; values are lists (sta,end) coords 

    @attention: no error check on the coordinate ranges is performed; crap in == crap out!

    @rtype:  alignedseqs: dict
    @return: alignedseqs: corrected dict

    @rtype:  alignment: string 
    @return: alignment: corrected string

    @rtype:  coords: dict
    @return: coords: corrected dict
    """
    keys = alignedseqs.keys()
    # correct left side
    while True:
        if not alignedseqs[keys[0]]: break
        for key in keys:
            if alignedseqs[key][0] == '-':
                for _key in keys:
                    symbol = alignedseqs[_key][0]
                    if symbol != '-': coords[_key][0]+=1
                    alignedseqs[_key] = alignedseqs[_key][1:]
                # all seqs in this iteration corrected for!
                alignment = alignment[1:] 
                break # break the for loop
            else:
                pass
        else:
            break # end the While loop!
    # correct left side
    while True:
        if not alignedseqs[keys[0]]: break
        for key in keys:
            if alignedseqs[key][-1] == '-':
                for _key in keys:
                    symbol = alignedseqs[_key][-1]
                    if symbol != '-': coords[_key][-1]-=1
                    alignedseqs[_key] = alignedseqs[_key][0:-1]
                # all seqs in this iteration corrected for!
                alignment = alignment[0:-1]
                break # break the for loop
            else:
                pass
        else:
            break # end the While loop!


    # correct CONSISTENT internal gaps!
    (alignedseqs,alignment,coords) = strip_consistent_internal_gaps(
            alignedseqs,alignment,coords)

    # return corrected data structures
    return ( alignedseqs, alignment, coords )

# end of function strip_alignment_for_exterior_gaps


def strip_overall_nonaligned_residues(alignedseqs,alignment,coords):
    """
    Strip all residues which have no 'alignment' conservation

    @type  alignedseqs: dict
    @param alignedseqs: dict ASIS produced by clustalw.clustalw() function

    @type  alignment: string
    @param alignment: string ASIS produced by clustalw.clustalw() function

    @type  coords: dict
    @param coords: dict with identical keys as alignedseqs; values are lists (sta,end) coords 

    @attention: no error check on the coordinate ranges is performed; crap in == crap out!

    @rtype:  alignedseqs: dict
    @return: alignedseqs: corrected dict

    @rtype:  alignment: string 
    @return: alignment: corrected string

    @rtype:  coords: dict
    @return: coords: corrected dict
    """
    keys = alignedseqs.keys()
    # correct left side
    while True:
        if not alignedseqs[keys[0]] or not alignment: break
        if alignment[0] == ' ':
            alignment = alignment[1:] 
            for key in keys:
                symbol = alignedseqs[key][0]
                if symbol != '-': coords[key][0]+=1
                alignedseqs[key] = alignedseqs[key][1:]
        else:
            break # end the While loop!
    # correct rigth side
    while True:
        if not alignedseqs[keys[0]] or not alignment: break
        if alignment[-1] == ' ':
            alignment = alignment[0:-1] 
            for key in keys:
                symbol = alignedseqs[key][-1]
                if symbol != '-': coords[key][-1]-=1
                alignedseqs[key] = alignedseqs[key][0:-1]
        else:
            break # end the While loop!

    # return corrected data structures
    return ( alignedseqs, alignment, coords )

# end of function strip_overall_nonaligned_residues



def strip_poorly_supported_tails(alignedseqs,alignment,coords,ratio):
    """
    Strip all residues which have no 'alignment' conservation

    @type  alignedseqs: dict
    @param alignedseqs: dict ASIS produced by clustalw.clustalw() function

    @type  alignment: string
    @param alignment: string ASIS produced by clustalw.clustalw() function

    @type  coords: dict
    @param coords: dict with identical keys as alignedseqs; values are lists (sta,end) coords 

    @type  ratio: float 
    @param ratio: minimal presence of residues in a column to be maintained 

    @attention: no error check on the coordinate ranges is performed; crap in == crap out!

    @rtype:  alignedseqs: dict
    @return: alignedseqs: corrected dict

    @rtype:  alignment: string 
    @return: alignment: corrected string

    @rtype:  coords: dict
    @return: coords: corrected dict
    """
    keys = alignedseqs.keys()
    # correct left side
    while True:
        if not alignedseqs[keys[0]] or not alignment: break
        if alignment[0] == ' ':
            presence=0.0
            for key in keys:
                symbol = alignedseqs[key][0]
                if symbol != '-':
                    presence+=1.0
            if presence / float(len(keys)) < ratio:
                # strip this position
                alignment = alignment[1:]
                for key in keys:
                    symbol = alignedseqs[key][0]
                    if symbol != '-': coords[key][0]+=1
                    alignedseqs[key] = alignedseqs[key][1:]
            else:
                break # end the While loop!
        else:
            break # end the While loop!
    # correct rigth side
    while True:
        if not alignedseqs[keys[0]] or not alignment: break
        if alignment[-1] == ' ':
            presence=0.0
            for key in keys:
                symbol = alignedseqs[key][-1]
                if symbol != '-':
                    presence+=1.0
            if presence / float(len(keys)) < ratio:
                # strip this position
                alignment = alignment[0:-1]
                for key in keys:
                    symbol = alignedseqs[key][-1]
                    if symbol != '-': coords[key][-1]-=1
                    alignedseqs[key] = alignedseqs[key][0:-1]
            else:
                break # end the While loop!
        else:
            break # end the While loop!

    # return corrected data structures
    return ( alignedseqs, alignment, coords )

# end of function strip_poorly_supported_tails 


def strip_consistent_internal_gaps(alignedseqs,alignment,coords):
    """
    Strip consitently aligned gaps from (a subset of) ClustalW aligned sequences

    @attention: REQUIRED when pairwise/subsets of ClustalW aligned sequences are processed

    @type  alignedseqs: dict
    @param alignedseqs: dict ASIS produced by clustalw.clustalw() function

    @type  alignment: string
    @param alignment: string ASIS produced by clustalw.clustalw() function

    @type  coords: dict
    @param coords: dict with identical keys as alignedseqs; values are lists (sta,end) coords 

    @attention: no error check on the coordinate ranges is performed; crap in == crap out!

    @rtype:  alignedseqs: dict
    @return: alignedseqs: corrected dict

    @rtype:  alignment: string 
    @return: alignment: corrected string

    @rtype:  coords: dict
    @return: coords: corrected dict
    """
    for pos in range(len(alignment)-1,-1,-1):
        if Set([ seq[pos] for seq in alignedseqs.values() ]) == Set("-"):
            # remove this position from all the sequences
            for k in alignedseqs.keys():
                alignedseqs[k] = alignedseqs[k][0:pos] + alignedseqs[k][pos+1:]

    # return corrected data structures
    return ( alignedseqs, alignment, coords )

# end of function strip_consistent_internal_gaps


def find_substring_in_clustalw_alignment(substring,alignedseqs,key=None):
    """
    Find an (AA) substring in the ClustalW alignment, taking gaps into account
    
    @type  substring: string
    @param substring: amino acid substring to find in the ClustalW alignment

    @type  alignedseqs: dict
    @param alignedseqs: dict ASIS produced by clustalw.clustalw() function

    @type  key: string
    @param key: dictionary key in alignedseqs to search; None searches all (in random order)

    @rtype:  tuple
    @return: (substring_sta, substring_end) relative to the ClustalW alignment
    """
    substring_pos = None
    substring_sta = None
    substring_end = None
    substring_len = len(substring.replace("-",""))
    for k,algseq in alignedseqs.iteritems():
        if key and k != key: continue
        # find substring in algseq
        substring_pos = algseq.replace("-","").upper().find(substring.upper())
        if substring_pos < 0:
            # substring not found
            substring_pos = None
            continue
        # correct match position for gaps
        aa_cnt  = 0
        gap_cnt = 0
        for char in algseq:
            if char != "-":
                aa_cnt+=1
            else:
                gap_cnt+=1
            if aa_cnt == substring_pos:
                substring_sta = substring_pos + gap_cnt
            elif aa_cnt == substring_pos + substring_len:
                substring_end = substring_pos + gap_cnt + substring_len
                break
            else:
                pass
                
    # return substring_position
    return substring_sta,substring_end

# end of function find_substring_in_clustalw_alignment


def clustalw_alignment_score_string(alignedseqs,matrix=None):
    """
    Generates an alternative sequence similarity measurement string for ClustalW protein alignments

    @type  alignedseqs: dict
    @param alignedseqs: dict ASIS produced by clustalw.clustalw() function

    @type  matrix: ProteinSimilarityMatrix
    @param matrix: ProteinSimilarityMatrix object (None deaults to BLOSUM62)

    @rtype:  tuple
    @return: (scores_as_a_list, scores_in_a_string)
    """
    # input sanity check 
    seq_cnt    = len(alignedseqs)
    if not seq_cnt: return [], ""
    seq_length = len(alignedseqs.values()[0])
    if not seq_length: return [], ""
    # get ProteinSimilarityMatrix if no matrix was provided
    if not matrix: matrix = ProteinSimilarityMatrix()
    score_string_list = []
    for pos in range(0,seq_length):
        chars   = [ seq[pos] for seq in alignedseqs.values() ]
        charset = set(chars)
        if len(charset) == 1:
            score_string_list.append(1.0)
        else:
            scores = []
            # loop over all the unique chars, and see which char gives
            # the highest sequence conservation pattern
            for char in charset:
                # omit gap char as a starting point
                if char == "-":
                    # count as half score! Although a gap should be considered
                    # very deteriorous in the alignment, MANY gaps on the same
                    # position is signal for an insertion in a few cases!
                    scores.append(float(chars.count(char))/2)
                else:
                    scores.append(float(chars.count(char)))
                for otherchar in charset:
                    if char != otherchar:
                        if matrix.scoreaapair(char,otherchar) > 0:
                            scores[-1] = scores[-1] + float(chars.count(otherchar))/2
            score_string_list.append(max(scores)/seq_cnt)
    # convert score_string_list to string
    score_string = []
    for score in score_string_list:
        if score == 1.0:
            score_string.append("*")
        elif round(score*10) == 10.0:
            score_string.append("*")
        else:
            score_string.append(str(int(round(score*10))))
    # return score float list and concatenated string
    return score_string_list, "".join(score_string)

# end of function clustalw_alignment_score_string


def clustalw_gap_string(alignedseqs):
    """
    Generates a gap presence measurement string for ClustalW protein alignments

    @type  alignedseqs: dict
    @param alignedseqs: dict ASIS produced by clustalw.clustalw() function

    @rtype:  tuple
    @return: (scores_as_a_list, scores_in_a_string)
    """
    # input sanity check 
    seq_cnt    = len(alignedseqs)
    if not seq_cnt: return [], ""
    seq_length = len(alignedseqs.values()[0])
    if not seq_length: return [], ""
    gap_string_list = []
    for pos in range(0,seq_length):
        gaps = [ seq[pos] for seq in alignedseqs.values() ].count("-")
        gap_string_list.append( float(gaps) / seq_cnt ) 

    # convert gap_string_list to string
    gap_string = []
    for score in gap_string_list:
        if score == 1.0:
            gap_string.append("*")
        elif round(score*10) == 10.0:
            gap_string.append("*")
        else:
            gap_string.append(str(int(round(score*10))))
    # return score float list and concatenated string
    return gap_string_list, "".join(gap_string)

# end of function clustalw_gap_string


def fix_X_signs_clustalw(seqs,align):
    """
    """
    headers = seqs.keys()
    for pos in range(0,len(align)):
        symbols = list(Set([ seqs[h][pos] for h in headers ]))
        if symbols == ['X']:
            # yep, an aligned 'X' amino-acid!
            if align[pos] != " ":
                align = align[0:pos]+" "+align[pos+1:]
    return seqs,align                

# end of function fix_X_signs_clustalw


def _test():
    """
    test function for clustalw functionality
    """
    print "TESTING lib_clustalw.py"

    # do the multiple alignment
    # to do
    # parse output file
    seqs,align= _parse_clustalw(fname)
    # fix "X" signs in the alignment
    seqs,align= fix_X_signs_clustalw(seqs,align)
    
    step = 75
    for offset in range(step,2500,step):
        for h,s in seqs.iteritems():
            print (h+" "*20)[0:20], s[offset-step:offset]
            #print "foxg_-----_0", len(align), align[0:25], align[250:275], align[-25:]
        print ("score"+" "*20)[0:20], align[offset-step:offset].replace(" "," ")
        print " "

    for h,s in seqs.iteritems():
        print h, len(s)


if __name__ == "__main__":
    _test()

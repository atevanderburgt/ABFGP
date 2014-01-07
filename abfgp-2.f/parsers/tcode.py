"""
TCODE parser functions; Developed & tested for TCODE version ?? 
"""

# Module metadata
__authors__ = "Ate van der Burgt"
__license__ = "MIT"

# Parser Imports
from exceptions import InproperlyAppliedArgument

# Python Imports
from os import popen2, popen3
from re import finditer

# Global variables
from executables import EXECUTABLE_TCODE

def run_tcode(fname="",outfile="",sequence="",step=3,window=200,EXECUTABLE_TCODE=EXECUTABLE_TCODE):
    """
    """
    if not outfile: outfile = ".".join(fname.split(".")[0:-1])+".tcode"
    # okay, do the tcode run
    command = "%s %s -step %s -window %s -outfile %s" % (
            EXECUTABLE_TCODE,
            fname,
            step,
            window,
            outfile,
            )
    ci,co,ce = popen3(command)
    ci.close()
    out = co.read()
    co.close()
    error = ce.read()
    ce.close()

# end of function run_tcode


def tcode(sequence=None,fname=None,outputfile=None,step=3,window=200,EXECUTABLE_TCODE=EXECUTABLE_TCODE):
    """
    Run EMBOSS tcode and obtain results

    @type  sequence: string
    @param sequence: input DNA sequence

    @type  fname: string
    @param fname: full path to getorf input file (fasta)

    @type  outputfile: string
    @param outputfile: full path to tcode output file

    @type  EXECUTABLE_TCODE: string
    @param EXECUTABLE_TCODE: full path to tcode executable
    """
    # do some integrity checks
    if not sequence and not fname:
        raise InproperlyAppliedArgument, "specify `sequence` or `fname` variable, not neither"
    if sequence and fname:
        raise InproperlyAppliedArgument, "specify `sequence` or `fname` variable, not both"
    if not EXECUTABLE_TCODE:
        raise InproperlyAppliedArgument, "specify `EXECUTABLE_TCODE` variable"

    # create command line, execute with popen and parse
    command = "%s -step %s -window %s" % (
            EXECUTABLE_TCODE,
            step,
            window,
            )

    # handle input (sequence on STDIN or file
    if fname:      command = "cat %s | %s" % (fname, command)
    else:          command = "echo %s | %s" % (sequence, command)
    # handle output (STDOUT or file
    if outputfile: command = "%s -outseq %s" % (command, outputfile) 
    else:          command = "%s -filter" % (command)
    # run the command with popen3
    ci,co,ce = popen3(command)
    ci.close()
    output = co.read()
    co.close()
    error = ce.read()
    ce.close()
    # no error capturing done here!
    if not outputfile:
        return output
    else:
        return outputfile

# end of function tcode 


def parsetcodestdout(stdout,offset=0):
    """ """
    return _parsetcode(stdout=stdout,offset=offset)
# end of function parsetcodestdout


def parsetcodeoutputfile(fname,offset=0):
    """ """
    return _parsetcode(fname=fname,offset=offset)
# end of function parsetcodeoutputfile
 


def _parsetcode(stdout=None,fname=None,offset=0):
    """
    offset corrects for N-symbols at the start of the sequence...
    """
    # do some integrity checks
    if not stdout and not fname:
        raise InproperlyAppliedArgument, "specify `stdout` or `fname` variable, not neither"
    if stdout and fname:
        raise InproperlyAppliedArgument, "specify `stdout` or `fname` variable, not both"
    elif stdout:
        # split by lines, ignore '^$' lines,'^#' lines and when done ignore the first (header)
        out = [ line for line in stdout.split("\n") if line and line[0]!="#" ][1:]
    else:
        # split by lines, ignore '^$' lines,'^#' lines and when done ignore the first (header)
        out = [ line for line in open(fname).readlines() if line and line[0]!="#" ][1:]

    retlist = []
    # in newer versions of EMBOSS tcode, output style has changed!
    # <= 4.0.0, 4 elements on a line
    # >= 6.2.0, 5 elements on a line
    # in between 4.0 and 6.2, I did not check ;-)
    # check the first line which output style it is
    IS_NEW_OUTPUT_STYLE = False
    firstline = out[0].strip()
    while "  " in firstline: firstline = firstline.replace("  "," ")
    parts = firstline.split(" ")
    if parts[2] in ['+','-']: 
        IS_NEW_OUTPUT_STYLE = True

    for line in out:
        line = line.strip()
        while "  " in line: line = line.replace("  "," ")
        if IS_NEW_OUTPUT_STYLE:
            parts = line.split(" ",4)
            parts.pop(2)
        else:
            parts = line.split(" ",3)
        retlist.append( (int(parts[0])+offset, int(parts[1])+offset, float(parts[2]), parts[3] ) )
    return retlist

# end of function _parsetcode


def parse_tcode_outfile(fname,offset=0):
    """
    offset corrects for N-symbols at the start of the sequence...
    """
    command = "cat %s | grep -v '#' | sed '/^$/d' | sed 1d" % (fname) 
    ci,co = popen2(command)
    ci.close()
    out = co.readlines()
    co.close()
    retlist = []
    for line in out:
        line = line.strip()
        while "  " in line: line = line.replace("  "," ")
        parts = line.split(" ",3)
        retlist.append( (int(parts[0])+offset, int(parts[1])+offset, float(parts[2]), parts[3] ) )
    return retlist

# end of function parse_tcode_outfile


def correct_sequence_for_n_tracks(sequence,window=100):
    """
    Dirty hack function; replace N*(>window) by centered Adenosines
    """
    # set window to a minimum of 100nt
    window=int(min([100,window]))
    if sequence.upper().find("N"*window) == -1:
        return (True,sequence)
    else:
        seqlist = list(sequence)
        for match in finditer("N{%s,}" % window,sequence.upper()):
            ntracklength = match.end() - match.start()
            iters = 1 
            while True:
                piecelength = ntracklength / (iters+1) 
                if piecelength <= window:
                    for pos in range(match.start()+piecelength,match.start()+ntracklength,window):
                        seqlist[pos] = 'a'
                    # and break out of this while loop
                    break
                else:
                    # nope, track is longer...
                    iters+=1
        # return status False and joined seqlist as new sequence        
        return (False,"".join(seqlist))

# end of function correct_sequence_for_n_tracks


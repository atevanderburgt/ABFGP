"""
Parser functions for EMBOSS getorf
"""

# Module metadata
__authors__ = "Ate van der Burgt"
__license__ = "MIT"

# Python Imports
from os import popen3 as osPopen3

# Parser Imports
from exceptions import InproperlyAppliedArgument
from fasta import parseFasta

# Global variables
from executables import (
    EXECUTABLE_GETORF,
    EXECUTABLE_GETORF_VERSION,
    EXECUTABLE_GETORF_MINSIZE,
    )


def parseGetorfOutput(content):
    """
    Parse getorf output file(handle), returns dictionaty with orfs

    @type  content: open filehandle
    @param content: filehandle of a EMBOSS getorf output fasta file

    @rtype:  dictionary
    @return: dictionary with keys (original fastaheader, orf-start, orf-stop) and sequence as values

    @attention: content can be: sys.stdin.readlines()
    @attention: content can be: fh.readlines()
    """
    orfs = {}
    for header, sequence in parseFasta(content).iteritems():
	# below 3 lines retrieves start coordinates from getorf output
	(fref,whatever)= header.split("_",1)
	(start,stop)  = header.split(" ",1)[-1][1:-1].split(" - ")
	(start,stop)  = (int(start),int(stop))
	name = (fref,start,stop)
	orfs[name] = sequence
    # return fasta dictionary with openreadingframe info
    return orfs

# end of function parseGetorfOutput


def __parseGetorfOutput(content):
    """
    Parse getorf output file(handle), returns dictionaty with orfs

    @type  content: open filehandle
    @param content: filehandle of a EMBOSS getorf output fasta file

    @rtype:  dictionary
    @return: dictionary with keys (original fastaheader, orf-start, orf-stop) and sequence as values

    @attention: content can be: sys.stdin.readlines()
    @attention: content can be: fh.readlines()
    """
    seqs = {}; name = ''; seq = '';
    for line in content:
        line = line.strip()
        if not line: continue
        if line[0] == '>':
            #new fasta seq, store if there was a previous one
            if name and seq: 
                seqs[name] = seq
            name = line[1:]
            seq = ''
            # below 3 lines retrieves start coordinates from getorf output
            (org,whatever)= name.split("_",1)
            (start,stop)  = name.split(" ",1)[-1][1:-1].split(" - ")
            (start,stop)  = (int(start),int(stop))
            name = (org,start,stop)
        else: seq += line
    if name and seq: seqs[name] = seq
    return seqs

# end of function parseGetorfOutput


def run_getorf(inputfile="",outputfile="",
    executable_getorf=EXECUTABLE_GETORF,
    minsize=EXECUTABLE_GETORF_MINSIZE):
    """
    Run EMBOSS getorf and write results to file

    @type  inputfile: string
    @param inputfile: full path to getorf input file (fasta)

    @type  outputfile: string
    @param outputfile: full path to getorf output file (fasta)

    @type  EXECUTABLE_GETORF: string
    @param EXECUTABLE_GETORF: full path to getorf executable
    """
    # do some integrity checks
    if not inputfile:
        raise InproperlyAppliedArgument, "specify `inputfile` variable"
    if not outputfile:
        raise InproperlyAppliedArgument, "specify `outputfile` variable"
    if not executable_getorf:
        raise InproperlyAppliedArgument, "specify `executable_getorf` variable"
    # create command line, execute with popen and parse
    command = "%s -sequence %s -outseq %s -minsize %s -noreverse" % (
            executable_getorf,
            inputfile,
            outputfile,
            minsize
            ) 
    ci,co,ce = osPopen3(command)
    ci.close()
    output = co.read()
    co.close()
    error = ce.read()
    ce.close()
    if output:
        print "getorf OUTPUT:", output
    if error and error.strip() != "Finds and extracts open reading frames (ORFs)":
        print "getorf ERROR:", error

# end of function run_getorf


def getorf(sequence=None,fname=None,outputfile=None,
    executable_getorf=EXECUTABLE_GETORF,
    minsize=EXECUTABLE_GETORF_MINSIZE):
    """
    Run EMBOSS getorf and write results to file

    @type  sequence: string
    @param inputfile: full path to getorf input file (fasta)

    @type  outputfile: string
    @param outputfile: full path to getorf output file (fasta)

    @type  executable_getorf: string
    @param executable_getorf: full path to getorf executable
    """
    # do some integrity checks
    if not sequence and not fname:
	message = "specify `sequence` or `fname` variable, not neither"
        raise InproperlyAppliedArgument, message
    if sequence and fname:
	message = "specify `sequence` or `fname` variable, not both"
        raise InproperlyAppliedArgument, message
    if not executable_getorf:
	message = "specify `EXECUTABLE_GETORF` variable"
        raise InproperlyAppliedArgument, message

    # create command line, execute with popen and parse
    command = "%s -minsize %s -noreverse" % (
                executable_getorf,
                minsize
                )
    if fname:      command = "cat %s | %s" % (fname, command)
    else:          command = "echo %s | %s" % (sequence, command)
    if outputfile: command = "%s -outseq %s" % (command, outputfile) 
    else:          command = "%s -filter" % (command)
    ci,co,ce = osPopen3(command)
    ci.close()
    output = co.read()
    co.close()
    error = ce.read()
    ce.close()
    if error and error.strip() != "Finds and extracts open reading frames (ORFs)":
        print "getorf ERROR:", error
    if not outputfile:
        return output
    else:
        return outputfile

# end of function getorf

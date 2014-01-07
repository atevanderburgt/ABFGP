"""
Python-Cexpander parser & functions

A library for working with the Cexpander (by E. Severing)

The cexpander software package contains three algorithms:
 prep_launch.py
 cbalignp
 cexpander_dr

For the desired cexpander output, the three agorithms must be runned in serial.
The output of each algorithm is the input file for the next step.
An example of these serial command line calls:
 python prep_launch.py MYDATA.fasta MYDATA.pairwise MYDATA.report
 ./cbalignp_test_2/src/cbalignp -i MYDATA.pairwise -y > MYDATA.aligned
 ./cexpander_dr < cexpander_test_file.inp > MYDATA.aus

The settings file 'cexpander_test_file.inp' describes how the main
programm (./cexpander_dr) is operated, and on which input file it operates.
So, for a Python call from this parser script, the 'cexpander_test_file.inp'
must be created too.

The equivalent in Python-Cexpander functions would be
 run_make_all_vs_all_fasta("MYDATA.fasta","MYDATA.pairwise")
 run_cbalignp("MYDATA.pairwise","MYDATA.aligned")
 run_cexpander_dr("MYDATA.aligned","MYDATA.aus")
 # parse the output file into an CexpanderOutput object
 cxpdrOutput = parse_cexpander("input.cexpander")

 # all the above in shortcutted by this function:
 cxpdrOutput = runcexpander("MYDATA.fasta")

For complete documentation of the cexpander package algorithms and their
command line arguments, see the software itself.

In this Python-Cexander codebox, only the subset of the cexander package
functionalities are covered.
"""

# Module metadata
__authors__ = "Ate van der Burgt"
__license__ = "MIT"


# Other Imports from parsers module
from fasta import writeMultiFasta, parseFasta
from pythonlibs.uniqueness import get_random_string_tag

# Python Imports
from os.path import isfile as osPathIsfile
from os import (
    system as osSystem,
    popen3 as osPopen3,
    popen as  osPopen,
    )


# Global variables / program paths
from executables import (
    EXECUTABLE_CEXPANDER_PATH,
    EXECUTABLE_CEXPANDER_ALLVSALL,
    EXECUTABLE_CEXPANDER_CBALIGNP,
    EXECUTABLE_CEXPANDER_CEXPANDER,
    )


# cexpander exceptions
class ZeroUniformlyAlignedPositions(Exception):
    """ Not a single uniformly aligned Amino Acid position """
    pass

# end of class ZeroUniformlyAlignedPositions


def run_make_all_vs_all_fasta(infname,outfname,reportfname,verbose=False):
    """
    Make pairwise combinations of sequences in input file with prep_launch.py

    @type  infname: string
    @param infname: (absolute) path to protein multi fasta 

    @type  outfname: string
    @param outfname: (absolute) path to prep_launch.py output file

    @type  reportfname: string
    @param reportfname: (absolute) path to prep_launch.py report file

    @type  verbose: Boolean
    @param verbose: debugging messages to STDOUT (True) or not (False, default)

    @attention: requires global variable EXECUTABLE_cexpander_ALLVSALL
    @attention: this used to be called `make_all_vs_all_fasta.py`

    @rtype:  Boolean
    @return: True
    """
    command = "python %s %s %s %s" % (
        EXECUTABLE_cexpander_ALLVSALL,
        infname,
        outfname,
        reportfname,
        )
    if verbose: print command
    _retval = osSystem(command)
    return True

# end of function run_make_all_vs_all_fasta


def run_cbalignp(infname,outfname,commandline="",verbose=False):
    """
    Run the cbalignp program by specifying input & output filenames etc.

    @type  infname: string
    @param infname: (absolute) path to make_all_vs_all_fasta.py output file

    @type  outfname: string
    @param outfname: (absolute) path to cbalignp output file

    @type  commandline: string
    @param commandline: literal additional command line arguments for cbalignp

    @type  verbose: Boolean
    @param verbose: print debugging messages to STDOUT (True) or not (False, default)

    @attention: requires global variable EXECUTABLE_CEXPANDER_CBALIGNP
    @attention: see cbalignp for (additional) command line options
    @attention: only a subset of cbalignp commandline options are supported!

    @rtype:  Boolean
    @return: True
    """
    command = "%s -i %s %s > %s " % (
        EXECUTABLE_CEXPANDER_CBALIGNP,
        infname,
        commandline,
        outfname
        )
    if verbose: print command
    _retval = osSystem(command)
    return True

# end of function run_cbalignp


def run_cexpander_dr(settingsfile,outfname,verbose=False):
    """
    Run the cexpander_dr program by specifying input & output filenames etc.

    @type  infname: string
    @param infname: (absolute) path to cbalignp output file

    @type  outfname: string
    @param outfname: (absolute) path to cexpander_dr output file

    @type  commandline: string
    @param commandline: literal additional command line arguments for cexpander_dr

    @type  verbose: Boolean
    @param verbose: print debugging messages to STDOUT (True) or not (False, default)

    @attention: requires global variable EXECUTABLE_CEXPANDER_CEXPANDER
    @attention: see cexpander_dr for (additional) command line options
    @attention: only a subset of cexpander_dr commandline options are supported!

    @rtype:  Boolean
    @return: True
    """
    command = "%s < %s > %s " % (
        EXECUTABLE_CEXPANDER_CEXPANDER,
        settingsfile,
        outfname
        )
    if verbose: print command
    ci,co,ce = osPopen3(command)
    ci.close()
    co.close()
    ce.close()
    return True

# end of function run_cexpander_dr


def parse_cexpander(cexpanderdata,fname_fasta):
    """
    Parse the cexpander_dr output file into a CexpanderOutput class object

    @type  fname_cexpander: string
    @param fname_cexpander: (absolute) path to cexpander_dr output file

    @type  fname_fasta: string
    @param fname_fasta: (absolute) path to fasta input file

    @rtype:  CexpanderOutput object
    @return: CexpanderOutput object
    """

    # open file to txt string and initialize empty CexpanderOutput object
    cxpOut = CexpanderOutput()
    cxpOut.sequences = parseFasta(open(fname_fasta).readlines())
    data   = cexpanderdata.split("\n\n")[0:-1]

    # generate header list; omit first 5 lines (cexpander STDOUT messages)
    headers = data.pop(0)
    headers = [ line.split("\t")[2].replace(">","") for line in headers.split("\n")[5:] ]

    # loop over the `transfer blocks` in the file
    for pos in range(0,len(data)):
        cxpTrfblck = CexpanderTransferBlock()
        if not set([ line.split("\t")[1] for line in data[pos][1:].split("\n")[1:] ]).difference(['0','1']):
            # cexpander in binary mode
            # single-line string of zeros (0) and ones (1)
            #cxpTrfblck.binarystring = ''.join( [ line.split("\t")[1] for line in data[pos][1:].split("\n")[1:] ] )
            mode = "binary"
            cxpTrfblck.binarystring = "" 
            for line in data[pos][1:].split("\n")[1:]:
                for cell in line.strip().split("\t")[1:]:
                    cxpTrfblck.binarystring+=cell
        else:
            # cexpander in float mode
            mode = "float"
            cxpTrfblck.binarystring = []
            for line in data[pos][1:].split("\n")[1:]:
                for cell in line.strip().split("\t")[1:]:
                    cxpTrfblck.binarystring.append(float(cell))

        cxpTrfblck.header       = headers[pos]
        cxpTrfblck.sequence     = cxpOut.sequences[cxpTrfblck.header]
        cxpTrfblck.positions    = len(cxpTrfblck.sequence)
        cxpTrfblck.uniform      = cxpTrfblck.get_uniform_positions()
        cxpTrfblck.score        = len(cxpTrfblck.uniform)
        cxpTrfblck.ratio        = cxpTrfblck._binarystring2matchratio(cxpTrfblck.binarystring)
        cxpTrfblck.mode         = mode

        if cxpTrfblck.positions != len(cxpTrfblck.binarystring):
            print "CEXPANDER-PARSE-ERROR:", cxpTrfblck.positions, "!=", len(cxpTrfblck.binarystring)
            print "CEXPANDER-PARSE-ERROR:", cxpTrfblck.header, mode
            print "CEXPANDER-PARSE-ERROR:", cxpTrfblck.sequence
            print "CEXPANDER-PARSE-ERROR:", cxpTrfblck.binarystring
            print "#"*40
            print "".join(data[pos])
            print "#"*40

        # add CexpanderTransferBlock to CexpanderOutput object
        cxpOut.add_transferblock(cxpTrfblck)

    # return the created object
    return cxpOut

# end of function parse_cexpander


class CexpanderTransferBlock:
    """
    Class describing a TransferBlock in a cexpander_dr output file

    @attention: (only) used asa subclass in the CexpanderOutput class
    """
    def __init__(self):
        """ initialize an empty cexpander TransferBlock """
        self.binarystring   = ""
        self.sequence       = ""
        self.uniform        = ""
        self.header         = ""
        self.threshold      = 1.0
        self.positions      = 0
        self.score          = 0
        self.ratio          = 0.0
    # end of function __init__

    
    def get_uniform_positions(self,threshold=1.0):
        """ """
        uniform   = ""
        for pos in range(0,len(self.binarystring)):
            if float(self.binarystring[pos]) >= threshold:
                try:
                    uniform = uniform + self.sequence[pos]
                except:
                    print uniform
                    print self.sequence
                    print self.binarystring
                    print self.header, len(self.binarystring), len(self.sequence), self.positions, pos
                    print self.tmp 
                    uniform = uniform + self.sequence[pos]
                   

        return uniform

    # end of function get_uniform_positions


    def _binarystring2matchratio(self,binarystring,threshold=1.0):
        """
        Obtain the relative number of uniformly matched positions (1) in
        the binarystring (of a TransferBlock)

        @type  binarystring: string
        @param binarystring: TransferBlock string with only zeros (0) and ones (1)

        @rtype:  float
        @return: ratio 0.0 .. 1.0 ( no matched positions / only zeros ) ..
                 ( only matched positions / only ones )
        """
        if type(binarystring) == type(str()):
            if binarystring.count("1") == len(binarystring):
                return 1.0
            elif binarystring.count("0") == len(binarystring):
                return 0.0
            else:
                return float(binarystring.count("1")) / len(binarystring)
        else:
            # a `binarystring` with floats
            cnt_gte = 0.0
            for elem in binarystring:
                if elem >= threshold:
                    cnt_gte+=1.0
            return cnt_gte/len(binarystring)
            
    # end of function _binarystring2matchratio


# end of class CexpanderTransferBlock


class CexpanderOutput(CexpanderTransferBlock):
    """ """
    def __init__(self):
        """ initialize an empty CexpanderOutput objects """
        # inherit all attributes from the CexpanderTransferBlock
        # the FIRST CexpanderTransferBlock that is added to
        # the _transferblocks attribute is stored directly into
        # the CexpanderOutput object attributes too for easy access
        CexpanderTransferBlock.__init__(self)
        self.sequences      = {}
        self.uniforms       = {}
        self._transferblocks= []

    # end of function __init__


    def add_transferblock(self,block):
        """ """
        self._transferblocks.append(block)
        self.uniforms[block.header] = block.uniform
        if len(self._transferblocks) == 1:
            self.set_transferblock(block.header)
    # end of function add_transferblock


    def set_transferblock(self,header):
        """ """
        if header in self.sequences.keys():
            for block in self._transferblocks:
                if block.header == header:
                    self.binarystring   = block.binarystring
                    self.sequence       = block.sequence
                    self.uniform        = block.uniform
                    self.header         = block.header
                    self.threshold      = block.threshold
                    self.positions      = block.positions
                    self.score          = block.score
                    self.ratio          = block.ratio
            else:
                return False
        else:
            return False
    # end of function set_transferblock


    def get_transferblock(self,header):
        """ """
        if header in self.sequences.keys():
            for block in self._transferblocks:
                if block.header == header:
                    return block
        else:
            return None
    # end of function get_block

    
    def get_formatted_binarystring(self):
        """
        Format the (float) binarystring asa string
        """
        if type(self.binarystring) == type(str()):
            return self.binarystring
        else:
            items = []
            for elem in self.binarystring:
                value = int(round(elem*10))
                if value == 10: value = "*"
                items.append(str(value))
            return "".join(items)
    
    # end of function get_formatted_binarystring
    
    
    def uniformly_matched_ratio(self,header=None):
        """
        Obtain the relative number of uniformly matched positions of a TransferBlock

        @type  header: string (or None)
        @param header: fasta sequence header; None uses current TransferBlock

        @rtype:  float (or None)
        @return: ratio of uniformly matched positions vs. all positions;
                 when fasta header is not found or TransferBlock is empty
                 (output is not completely parsed), return value is None.
        """
        if header and header in self.sequences.keys():
            for trfblck in self._transferblocks:
                if trfblck.header == header:
                    if trfblck.binarystring:
                        return trfblck._binarystring2matchratio(
                            trfblck.binarystring)
                    else:
                        return None
        elif header:
            return None
        elif not hasattr(self,"binarystring"):
            return None
        elif not self.binarystring:
            return None
        else:
            return self._binarystring2matchratio(
                self.binarystring)

    # end of function uniformly_matched_ratio

# end of class CexpanderOutput


def runcexpander(fname_fasta,cbalignp_commandline=" -y",output='binary'):
    """
    Run the complete cascade of cexpander algorithms on an input multi fasta
    file and return the output as a CexpanderOutput object

    @type  fname_fasta: string
    @param fname_fasta: path to input multi fasta file

    @type  cbalignp_commandline: string
    @param cbalignp_commandline: (extra) command line for cbalignp

    @type  min_cols: integer
    @param min_cols: minimal number of uniformly matched positions (cols)
                     required to report transfer blocks for (>= 0)

    @type  projected_on: string
    @param projected_on: apply fasta seqeunce header which to use for projection;
                         apply ':::' to do projections on all input sequences

    @attention: requires global variable EXECUTABLE_cexpander_ALLVSALL
    @attention: requires global variable EXECUTABLE_CEXPANDER_CBALIGNP
    @attention: requires global variable EXECUTABLE_CEXPANDER_CEXPANDER
    @attention: see cexpander_dr for (additional) command line options
    @attention: only a subset of cexpander_dr commandline options are supported!
    
    @rtype:  CexpanderOutput object
    @return: CexpanderOutput object
    """
    if not fname_fasta: raise "NoProperFunctionArguments"
    if not osPathIsfile(fname_fasta): raise "FileDoesNotExist"

    # (0) create (~unique) filenames
    uniquetag = get_random_string_tag() 
    fname_allvsall  = ".".join([fname_fasta,uniquetag,"allvsall"])
    fname_report    = ".".join([fname_fasta,uniquetag,"report"])
    fname_aligned   = ".".join([fname_fasta,uniquetag,"aligned"])
    fname_settings  = ".".join([fname_fasta,uniquetag,"settings"])
    fname_cexpander = ".".join([fname_fasta,uniquetag,"cexpander"])

    # (1) create complete .fa -> cexpanderstring command
    command = """
        python %s %s %s %s;
        %s -i %s %s > %s;
        %s < %s;
        """ % (
        EXECUTABLE_CEXPANDER_ALLVSALL,
        fname_fasta,
        fname_allvsall,
        fname_report,
        EXECUTABLE_CEXPANDER_CBALIGNP,
        fname_allvsall,
        cbalignp_commandline,
        fname_aligned,
        EXECUTABLE_CEXPANDER_CEXPANDER,
        fname_settings,
        )


    # (2) create fname_settings file
    binorfloat = "$dumpcv"
    if output == "float": binorfloat = "$dumpcvc"
    fh = open(fname_settings,'w')
    content = "\n\n".join( [
        "$load\n%s\n%s" % (fname_report,fname_aligned),
        "$addquery\n-1",
        "$run",
        "$dumpentries",
        "$cv_linear",
        "%s" % ( binorfloat ), # BINARY == $dumpcv, FLOAT = $dumpcvc
        "$exit\n\n", 
        ] )
    fh.write(content)
    fh.close()


    # (3) run the command
    ci,co,ce = osPopen3(command)
    ci.close()
    # output of EXECUTABLE_CEXPANDER_ALLVSALL is cast to STDOUT as well!
    cexpanderdata = co.read()
    co.close()
    error = ce.read()
    ce.close()

    # (4) parse fname_cexpander to CexpanderOutput object
    cxpdr = parse_cexpander(cexpanderdata,fname_fasta)

    # (5) cleanup files
    osSystem("rm -f %s %s.%s.*" % ( fname_fasta, fname_fasta,uniquetag ) )

    # (6) return the output object
    return cxpdr

# end of function runcexpander


def _parse_cexpander_string(fname,verbose=False):
    """
    @attention: deprecated function, not used in any of the functions above
    """
    command = """ cat %s | grep "\$start_values" -A 1000 | grep "\$end_values" -m 1 -B 1000 | sed 's/^.*_values$//' | sed '/^$/d' | tr -d "\\t" | tr -d "\\n" """ % fname
    if verbose: print command
    co = osPopen(command)
    cexpanderstring = co.read().strip()
    co.close() 
    return cexpanderstring

# end of function _parse_cexpander_string



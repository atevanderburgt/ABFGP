import os

from settings.executables import (
    EXECUTABLE_FASTALENGTH,
    EXECUTABLE_GFFLENGTH,
    )

def fastafilesequencelength(fname,executable=EXECUTABLE_FASTALENGTH):
    """
    Get the sequence length of a SINGLE SEQUENCE fasta file

    @type  fname: string
    @param fname: (full) path to fasta file

    @type  executable: string
    @param executable: (full) path to fastalength.sh executable

    @rtype:  integer
    @return: sequence length
    """
    return int(os.popen("%s %s" % (executable,fname)).read())

# end of function fastafilesequencelength


def gfffilesequencelength(fname,executable=EXECUTABLE_GFFLENGTH):
    """
    Get the sequence track length of a SINGLE LINE gff file

    @type  fname: string
    @param fname: (full) path to gff file

    @type  executable: string
    @param executable: (full) path to gfflength.sh executable

    @rtype:  integer
    @return: sequence length
    """
    return int(os.popen("cat %s | bash %s" % (fname,executable)).read())

# end of function gfffilesequencelength


"""
"""

# Module metadata
__authors__ = "Ate van der Burgt"
__license__ = "MIT"

# import parent pacb class
import pacb
from validators import IsPacbPORF, IsPacbP
from exceptions import NoPacbPORF

def pacbporf2pacbp(pacbporf):
    """
    Convert PacbPORF object (backwards to) a PacbP object

    @type  pacbporf: pacb.PacbPORF instance
    @param pacbporf: pacb.PacbPORF instance 

    @rtype:  pacb.PacbP instance
    @return: pacb.PacbP instance 
    """
    if str(pacbporf.__class__) == 'pacb.PacbP':
        # object is already a PacbP, not a PacbPORF
        return pacbporf
    elif str(pacbporf.__class__) == 'pacb.pacbp.PacbP':
        # object is already a PacbP, not a PacbPORF
        return pacbporf
    else:
        # object isa PacbPDNA or PacbPDNA; convert to PacbP
        staPos = pacbporf._get_original_alignment_pos_start()
        endPos = pacbporf._get_original_alignment_pos_end()
        query  =  pacbporf.query[staPos.position:endPos.position+1]
        sbjct  =  pacbporf.sbjct[staPos.position:endPos.position+1]
        pacbpinput = (query,sbjct,staPos.query_pos,staPos.sbjct_pos)
        pacbpobj = pacb.PacbP(input=pacbpinput,MATRIX=pacbporf.MATRIX)
        # update required attributes
        pacbpobj.source = pacbporf.source
        pacbpobj._IS_DELETE_PROTECTED = pacbporf._IS_DELETE_PROTECTED
        pacbpobj._IS_EDIT_PROTECTED = pacbporf._IS_EDIT_PROTECTED
        # return the object
        return pacbpobj

# end of function pacbporf2pacbp


def pacbp2pacbporf(pacbp,orfQ,orfS):
    """
    Convert a PacbP object into a PacbPORF by applying its Orf objects

    @type  pacbp: pacb.PacbP instance 
    @param pacbp: pacb.PacbP instance 

    @type  orfQ: Orf instance 
    @param orfQ: Orf Query object 

    @type  orfS: Orf instance
    @param orfS: Orf Sbjct object 

    @rtype:  pacb.PacbPORF instance 
    @return: pacb.PacbPORF instance 

    @attention: new and only correct function name for earlier _return_into_pacbporf
    """
    try:
        IsPacbPORF(pacbp)
        return pacbp
    except NoPacbPORF:
        IsPacbP(pacbp)
        return pacb.PacbPORF(pacbp,orfQ,orfS)

# end of function pacbp2pacbporf


def _return_into_pacbporf(pacbp,orfQ,orfS):
    """
    Convert a PacbP object into a PacbPORF by applying its Orf objects

    @attention: DEPRECATED! use new function name pacbp2pacbporf
    """
    return pacbp2pacbporf(pacbp,orfQ,orfS)

# end of function _return_into_pacbporf


def pacbp2pacbpdna(pacbp,orfQ,orfS):
    """
    @attention: new alias name for make_PacbPDNA 
    """
    return make_PacbPDNA(pacbp,orfQ,orfS)

# end of function pacbp2pacbpdna


def make_PacbPDNA(pacbp,orfQ,orfS):
    """
    Convert a PacbP object into a PacbPDNA by applying its Orf objects

    @type  pacbp: pacb.PacbP instance
    @param pacbp: pacb.PacbP instance

    @type  orfQ: Orf instance
    @param orfQ: Orf Query object 

    @type  orfS: Orf instance
    @param orfS: Orf Sbjct object 

    @rtype:  pacb.PacbPDNA instance
    @return: pacb.PacbPDNA instance

    @attention: no extention on the aligned sequence range is performed!
    """
    # get nucleotide sequence of orfQuery
    orfQdnaseq = orfQ.getntseqslice(
            abs_aa_startPY = pacbp.query_start,
            abs_aa_endPY   = pacbp.query_end )
    # get nucleotide sequence of orfSbjct
    orfSdnaseq = orfS.getntseqslice(
            abs_aa_startPY = pacbp.sbjct_start,
            abs_aa_endPY   = pacbp.sbjct_end )
    # make hsp alignments of both extended hsp's
    pacbpdna = pacb.PacbPDNA(
            pacbp = pacbp,
            input_dna = (
                orfQdnaseq,
                orfSdnaseq,
                orfQ.aapos2dnapos(pacbp.query_start),
                orfS.aapos2dnapos(pacbp.sbjct_start)
                )
            )


    # update required attributes
    pacbpdna.source = pacbp.source
    pacbpdna._IS_DELETE_PROTECTED = pacbp._IS_DELETE_PROTECTED
    pacbpdna._IS_EDIT_PROTECTED = pacbp._IS_EDIT_PROTECTED
    # and return
    return pacbpdna

# end of function make_PacbPDNA


def pacbp_from_clustalw(alignment=(),coords=()):
    """
    Create a PacbP object from a ClustalW alignment

    @type  alignment: tuple (of 3 strings)
    @param alignment: (query,match,sbjct) in !ALIGNED! format (indentical string lengths with gaps)

    @type  coords: tuple (of 4 integers)
    @param coords: (qaastart,qaaend,saastart,saaend) integers, AA-coords!

    @rtype:  pacb.PacbP instance
    @return: pacb.PacbP instance
    """
    query,match,sbjct = alignment
    qaastart,qaaend,saastart,saaend = coords
    # trim leading and tailing gaps in the alignment
    while query and query[0] == '-':
        query = query[1:]
        match = match[1:]
        sbjct = sbjct[1:]
        saastart+=1
    while sbjct and sbjct[0] == '-':
        query = query[1:]
        match = match[1:]
        sbjct = sbjct[1:]
        qaastart+=1
    while query and query[-1] == '-':
        query = query[0:-1]
        match = match[0:-1]
        sbjct = sbjct[0:-1]
        saaend-=1
    while sbjct and sbjct[-1] == '-':
        query = query[0:-1]
        match = match[0:-1]
        sbjct = sbjct[0:-1]
        qaaend-=1
    # trim non-aligned characters as well!
    while match and match[0] == ' ':
        query = query[1:]
        match = match[1:]
        sbjct = sbjct[1:]
        qaastart+=1
        saastart+=1
    while match and match[-1] == ' ':
        query = query[0:-1]
        match = match[0:-1]
        sbjct = sbjct[0:-1]
        qaaend-=1
        saaend-=1

    if query and sbjct:
        # make a PacbPORF of this alignment
        pacbpinput = (query,sbjct,qaastart,saastart)
        pacbp      = pacb.PacbP(input=pacbpinput)
        pacbp.source = 'clustalw'
    else:
       # no proper clustalw alignable sequences -> no PacbP
       pacbp = None
    # return the created pacbp
    return pacbp

# end of function pacbp_from_clustalw


def exononorfs2pacbporf(exonQ,exonS,matrix=None):
    """
    Create a PacbPORF object from 2 ExonOnOrf objects

    @type  exonQ: ExonOnOrf object
    @param exonQ: ExonOnOrf object (query)

    @type  exonS: ExonOnOrf object
    @param exonS: ExonOnOrf object (sbjct)

    @attention: exonQ and exonS proteinsequence() MUST be identical in length!

    @rtype:  pacb.PacbPORF instance
    @return: pacb.PacbPORF instance
    """
    # prepare input data
    query = exonQ.proteinsequence()
    sbjct = exonS.proteinsequence()
    qaas  = exonQ.orf.dnapos2aapos(exonQ.start)
    saas  = exonS.orf.dnapos2aapos(exonS.start)
    # make pacbp and then pacbporf object
    if matrix:
        pacbpobj = pacb.PacbP(input=(query,sbjct,qaas,saas),MATRIX=matrix)
    else:
        pacbpobj = pacb.PacbP(input=(query,sbjct,qaas,saas))
    if pacbpobj.length == 0:
        return None
    pacbporfobj = pacbp2pacbporf(pacbpobj,exonQ.orf,exonS.orf)
    pacbporfobj.extend_pacbporf_after_stops()
    return pacbporfobj

# end of function exononorfs2pacbporf


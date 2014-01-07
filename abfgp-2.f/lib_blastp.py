""" blast parser to NCBI's blastall -p blastp used in ABFGP; version blastall2.2.8 is supported """

import sys
from Bio.Blast import NCBIStandalone
from StringIO import StringIO
from os import (
    system as OSsystem,
    popen3 as osPopen3,
    remove as osRemove,
    getcwd as OSgetcwd,
    )
from os.path import join as osPathJoin
from re import compile, finditer
from copy import deepcopy
from sets import Set

from pythonlibs.uniqueness import get_random_string_tag

from settings.executables import (
    EXECUTABLE_BLASTALL as BLASTALL_PATH,
    EXECUTABLE_FORMATDB as FORMATDB_PATH
    )

from settings.emboss import (
    BLOSUM62_PATH,
    BLOSUM80_PATH,
    BLOSUM45_PATH,
    PAM30_PATH,
    PAM70_PATH,  
    )


#
#       blast matrices information
#       http://flybase.net/static_pages/blast/matrix_info.html
#
#     Query length     Substitution matrix     Gap costs
#     ------------     -------------------     ---------
#     <35              PAM-30                  ( 9,1)
#     35-50            PAM-70                  (10,1)
#     50-85            BLOSUM-80               (10,1)
#     >85              BLOSUM-62               (11,1)

def formatdb(fastadata={},fname=""):
    """
    @seqs    dict in fasta style with input sequences
    @returns name of the blastdb(.*)
    """
    if not fastadata and fname:
        OSsystem("%s -i %s" % (FORMATDB_PATH,fname))
    elif fastadata and fname:
        pass
    else:
        raise "inproper input"
    return fname

# end of function formatdb


def blastall_seq2seq(fastadata=(),filenames=(),output="ncbiparsed",blastprogram="blastp",remove_files=True,extra_blastp_params={'F': 'F', 'e': '10'}):
    """
    choose proper input:
    fastadata   ( ( headerQUERY, seqQUERY ) , ( headerSBJCT, seqSBJCT ) )
     or
    filenames   ( filenameQUERY, filenameSBJCT )
    """
    input = None

    if blastprogram not in ['blastp','tblastn','tblastx','blastx']:
        raise "only blastp and tblastn are supported"
    elif blastprogram in ['tblastn','tblastx']:
        dna_or_prot = "F"
    else:
        dna_or_prot = "T"

    if fastadata and type(fastadata) == type(()) and len(fastadata) == 2 and not filenames:
        # input is fasta headers and sequence
        input = "fastadata"
        # write input filenames
        uniquetag = get_random_string_tag()
        fname_q = "_".join( [ uniquetag, str(fastadata[0][0]), 'Q.fa' ] )
        fname_s = "_".join( [ uniquetag, str(fastadata[1][0]), 'S.fa' ] )
        fh = open(fname_q,'w')
        fh.write(">%s\n%s" % (fastadata[0][0],fastadata[0][1]))
        fh.close()
        fh = open(fname_s,'w')
        fh.write(">%s\n%s" % (fastadata[1][0],fastadata[1][1]))
        fh.close()
    elif filenames and type(filenames) == type(()) and len(filenames) == 2 and not fastadata:
        # input is (supposed to be) filenames
        input = "filenames"
        # get filenames
        fname_q = filenames[0]
        fname_s = filenames[1]
    elif not filenames and not fastadata:
        raise "no input!"
    else:
        raise "inproper input!"

    # formatdb
    OSsystem("%s -i %s -p %s" % (FORMATDB_PATH,fname_s,dna_or_prot))
    # and blastall!
    extra_params = " ".join(["-%s %s" % (k,v) for k,v in extra_blastp_params.iteritems()])
    ci,co,ce = osPopen3("%s -p %s %s -i %s -d %s " % (BLASTALL_PATH,blastprogram,extra_params,fname_q,fname_s))
    ci.close()
    if output == "ncbiparsed":
        b_parser = NCBIStandalone.BlastParser()
        blastallout = b_parser.parse(co)
    else:
        blastallout = co.read()
    co.close()
    ce.close()
    if remove_files:
        OSsystem("rm %s.*" % fname_s)
        osRemove("%s" % fname_s)
        osRemove("%s" % fname_q)
    # and return!
    return blastallout
    
# end of function blastall_seq2seq


def blastall_seq2db(header,sequence,dbname="",blastprogram="blastp",output="ncbiparsed",extra_blastp_params={'F': 'F', 'e': '10'}):
    """
    """
    if blastprogram not in ['blastp','tblastn','blastn','blastx']:
        raise "only blastp and tblastn are supported"

    extra_params = " ".join(["-%s %s" % (k,v) for k,v in extra_blastp_params.iteritems()])
    # generate (semi ;-) unique filename
    uniquetag = get_random_string_tag()
    fname = "_".join( [ uniquetag, str(header).replace(" ","_"), sequence[0:10]+".fa" ] )
    fname = osPathJoin(OSgetcwd(),fname)
    fh = open(fname,'w')
    fh.write(">%s\n%s\n" % (header,sequence))
    fh.close()
    command = "%s -p %s %s -i %s -d %s " % (BLASTALL_PATH,blastprogram,extra_params,fname,dbname)
    try:
        ci,co,ce = osPopen3(command)
        ci.close()
        if output == "ncbiparsed":
            b_parser = NCBIStandalone.BlastParser()
            blastallout = b_parser.parse(co)
        else:
            blastallout = co.read()
        co.close()
        ce.close()
    except:
        # for some kind of - obvious or freak accident case -
        # Blast or parsing of the blast record failed
        # No debugging here; just cleanup and return False
        print "BLAST CRASHED::"
        print command
        blastallout = False

    # remove the created Query file
    osRemove(fname)
    # and return!
    return blastallout

# end of function blastall_seq2db



def blastall_file2db(fname,dbname="",blastprogram="blastp",output="ncbiparsed",extra_blastp_params={'F': 'F', 'e': '10'}):
    """
    """
    if blastprogram not in ['blastp','tblastn','blastn','tblastx']:
        raise "only blastp and tblastn are supported"

    extra_params = " ".join(["-%s %s" % (k,v) for k,v in extra_blastp_params.iteritems()])
    command = "%s -p %s %s -i %s -d %s " % (BLASTALL_PATH,blastprogram,extra_params,fname,dbname)
    try:
        ci,co,ce = osPopen3(command)
        ci.close()
        if output == "ncbiparsed":
            b_parser = NCBIStandalone.BlastParser()
            blastallout = b_parser.parse(co)
        else:
            blastallout = co.read()
        co.close()
        ce.close()
        # do NOT remove the input fname
    except:
        co.close()
        error = ce.read().strip()
        ce.close()
        print command
        print "ERROR: '%s'" % error
        raise "BLAST CRASHED...."
    # and return!
    return blastallout

# end of function blastall_file2db


def get_orfs_of_proteinseq(dbname,protein_fname):
    """
    # TODO: This step is not 100% accurate, it should be abstracted from the gene's gff and the reading frame.
    # TODO: Because we have gff in most cases, it can be solved in that fashion.
    # TODO: If no gff is available, this should be tackled by a dedicated algorithm as SCIPIO.
    """
    header, sequence = open(protein_fname).read()[1:].strip().split("\n",1)
    sequence = sequence.replace("\n","")
    params = {'F': 'F', 'g': 'F', 'M': 'BLOSUM80', 'e': '10'}
    brec = blastall_seq2db(header,sequence,dbname=dbname,output="ncbiparsed",extra_blastp_params=params)

    continuoussmatch = compile("[A-Z]{4,}")
    exons = []

    for alignment in brec.alignments:
        for hsp in alignment.hsps:
            if hsp.gaps != (None,None):
                # gapped alignment !? that is not okay!
                raise "GAPPED ALIGNMENT encountered in a hsp...."
            for mo in finditer(continuoussmatch,hsp.match):
                exon = deepcopy(hsp)
                offset = 0
                offset_end = len(hsp.match)
                if len(hsp.match) != len(mo.group()):
                    offset = hsp.match.find(mo.group())
                    exon.query_start+=offset
                    exon.sbjct_start+=offset
                    exon.query = exon.query[offset:]
                    exon.match = exon.match[offset:]
                    exon.sbjct = exon.sbjct[offset:]
                    if len(exon.query) != len(mo.group()):
                        offset_end = len(exon.query) - len(mo.group())
                        exon.query = exon.query[:-offset_end]
                        exon.match = exon.match[:-offset_end]
                        exon.sbjct = exon.sbjct[:-offset_end]
                if len(exon.match) != len(mo.group()):
                    raise "ERORORORORORORO"
                else:
                    exon.identities = ( len(exon.match), len(exon.match) )
                # and append to exons
                exon.sbjct_name = alignment.title.replace(">","")
                exons.append(exon)
    exons = _order_exons_by_length(exons)
    true_exons = []
    if exons:
        true_exons = [ exons[0] ]
        for exon in exons[1:]:
            extension = True
            exoncoords = Set(range(exon.query_start,(exon.query_start+len(exon.query)+1) ) )
            for true_exon in true_exons:
                true_exoncoords = Set(range(true_exon.query_start,(true_exon.query_start+len(true_exon.query)+1) ) )
                if exoncoords.issubset(true_exoncoords):
                    extension = False
                    break
            if extension:
                true_exons.append(exon)
    # and return the ORFids
    orfids = list(Set([ int(exon.sbjct_name.split("_")[-1]) for exon in true_exons ]) )
    return ( orfids, true_exons )


def _order_exons_by_length(exons):
    retexons = [ ( len(exon.query), exon ) for exon in exons ]
    retexons.sort()
    retexons.reverse()
    return [ exon for (ll, exon) in retexons ]

# end of function _order_exons_by_length

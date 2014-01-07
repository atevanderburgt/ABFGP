"""
Annotate an *unigene* on its reference *genome* sequence.

  *unigene* is represented by gff of its aligned exons (e.g. GenomeThreader) on STDIN
  *genome*  is represented by multi-fasta file as (first and only) command line argument

Annotation stati can be the following:

  full         5'UTR - M .. exon(s) .. - * - 3'UTR 
  cds5p        5'UTR - M .. exon(s) .. - *
  cds3p                M .. exon(s) .. - * - 3'UTR
  cds                  M .. exon(s) .. - *
  utr5p        5'UTR - M .. exon(s) .. - NO *        orf at least x AA
  utr3p      NO * / NO M .. exon(s) .. - * - 3'UTR   orf at least x AA
  fragm                   NO *
  nonc          no obvious open reading frame

  5'UTR        is defined as having a * symbol in it!

IMPORTANT NOTICE: gff of unigene exons must be ordered ASC by gff fstart coordinate!!

"""

# Module metadata
__authors__ = "Ate van der Burgt"
__license__ = "MIT"

# Python Imports
import sys, os
from re import compile, finditer

# Python library Imports
from lib_fasta import parseFasta

# Global variable import
from settings.executables import (
        PYTHON_PATH,           # path to python executable
        EXECUTABLE_TRANSEQ,    # path to transeq from EMBOSS
        EXECUTABLE_GFF2FASTA,  # path to get_gff_sequence_from_fasta.py tool
        )
MIN_ORF_AA_LENGTH = 50


# get python version (2.6.6 -> 2.6, 2.4.3 -> 2.4
# required for os.popen2 / subprocess.Popen import
python_version = float(".".join(sys.version.split(" ")[0].split(".",2)[0:2])) 

# obtain genome multi-fasta file name
fastafile = os.path.join( os.getcwd(), sys.argv[1] )
# check for other command line arguments
output_protein = False
transeq_frame = "F"
print_all_frames = False
if len(sys.argv) > 2:
    arguments = sys.argv[2:]
    if '--prot' in arguments:
        output_protein = True
    if '--forward' in arguments:
        transeq_frame = "F"
    elif '--reverse' in arguments:
        transeq_frame = "R"
    elif '--both' in arguments:
        transeq_frame = "6"
    else:
        transeq_frame = "F"
    if '--sixframe' in arguments:
        print_all_frames = True
        transeq_frame = "6"
        output_protein   = True
        MIN_ORF_AA_LENGTH = 30
    else:
        print_all_frames = False

# read unigene exon(s) from STDIN
unigenegffexons = sys.stdin.readlines()

# obtain strand of the alignment
fstrand = unigenegffexons[0].split("\t")[6]
# reverse unigenegffexon list in case of negative strand
if fstrand == '-': unigenegffexons.reverse()


def _obtain_unigene_gff_nt_start_coord(gffs,fstrand,readingframe,ugntlength,inputprotseqlength,aaoffset):
    """
    # TODO: CHECK IF COORD CORRECTION +/- readingframe IS CORRECT IN NEGATVIE STRAND!!!
    """
    if fstrand == '-':
        ntcoord = max([ int(gff[4]) for gff in gffs ]) - aaoffset*3 - readingframe
        return ntcoord
    else:
        ntcoord = min([ int(gff[3]) for gff in gffs ]) + aaoffset*3 + readingframe
        # correct for leading utr5p exons that are fully utr
        for i in range(1,len(gffs)):
            summedexonlength = sum([ (int(gff[4])-int(gff[3])+1) for gff in gffs[0:i] ])
            if summedexonlength < aaoffset*3:
                curintronlength = int(gffs[i][3])-1-int(gffs[i-1][4])
                ntcoord+= curintronlength
            else:
                break
        return ntcoord

# end of function _obtain_unigene_gff_nt_start_coord


def _obtain_unigene_gff_nt_end_coord(gffs,fstrand,readingframe,ugntlength,inputprotseqlength,aaoffset):
    """
    # TODO: CHECK IF COORD CORRECTION +/- readingframe IS CORRECT IN NEGATIVE STRAND!!!

    # WEEEEEIRD EXCEPTION; since shifting algorithm to calculation cluster (~feb 2010),
    # a newer version of EMBOSS.transeq is used. This causes very weird behavious in
    # reading frame 2 output. Suddenly, frame 2 protseq is sometimes 1 AA taller!
    # I suppose this was a bug in the older version of transeq, because this output
    # seems to be the correct output. To have backwards compatibility with ALL transeq
    # versions, check the length of the nt sequence with the protein sequence.
    # in case of a 'missing' AA -> correct for it!
    """
    if fstrand == '-':
        return min([ int(gff[3]) for gff in gffs ]) + aaoffset*3
    else:
        framecorrection = (ugntlength - readingframe) % 3
        if ugntlength%3 > readingframe: framecorrection -= 3 # (shift one AA position)

        # WEEEEEIRD EXCEPTION; found march 2010...
        if readingframe == 2 and ugntlength-(inputprotseqlength*3) == 1:
            framecorrection -= 3
        elif readingframe == 1 and ugntlength-(inputprotseqlength*3) == 0:
            framecorrection -= 3
        elif readingframe == 2 and ugntlength-(inputprotseqlength*3) == 0:
            framecorrection -= 3
        else:
            pass


        ntcoord = max([ int(gff[4]) for gff in gffs ]) - aaoffset*3 - framecorrection
        # correct for tailing utr3p exons that are fully utr
        for i in range(len(gffs)-1,-1,-1):
            summedexonlength = sum([ (int(gff[4])-int(gff[3])+1) for gff in gffs[i:] ])
            if summedexonlength < aaoffset*3:
                curintronlength = int(gffs[i][3])-1-int(gffs[i-1][4])
                ntcoord-= curintronlength
            else:
                break
        return ntcoord

# end of function _obtain_unigene_gff_nt_end_coord



# get fasta by gff, remove headers, concatenate into single DNA seq and make 3-frame translation
command = """%s %s %s | grep -v "^>" | tr -d " \n" | %s -filter -frame %s """ % (
        PYTHON_PATH,
        EXECUTABLE_GFF2FASTA,
        fastafile,
        EXECUTABLE_TRANSEQ,
        transeq_frame )

if python_version <= 2.5:
    # Python2.4 os.popen syntax
    ci,co = os.popen2(command)
    ci.write("\n".join(unigenegffexons))
    ci.close()
    frametrans = parseFasta( co.readlines() )
    co.close()
else:
    # Python2.6 subprocess.Popen syntax
    from subprocess import Popen, PIPE
    p = Popen(command,shell=True,stdin=PIPE, stdout=PIPE, close_fds=True)
    (child_stdin, child_stdout) = (p.stdin, p.stdout)
    child_stdin.write("\n".join(unigenegffexons))
    child_stdin.close()
    frametrans = parseFasta( child_stdout.readlines() )


# command to get DNA sequence of unigene exon(s)
command = """%s %s %s | grep -v "^>" """ % (
        PYTHON_PATH,
        EXECUTABLE_GFF2FASTA,
        fastafile,
        )

if python_version <= 2.5:
    # Python2.4 os.popen syntax
    ci,co = os.popen2(command)
    ci.write("\n".join(unigenegffexons))
    ci.close()
    _seq = co.read() 
    co.close()
else:
    # Python2.6 subprocess.Popen syntax
    p = Popen(command,shell=True,stdin=PIPE, stdout=PIPE, close_fds=True)
    (child_stdin, child_stdout) = (p.stdin, p.stdout)
    child_stdin.write("\n".join(unigenegffexons))
    child_stdin.close()
    _seq = child_stdout.read()


# redefine unigenegffexons; now as tuples of gff
for i in range(0,len(unigenegffexons)):
    parts = unigenegffexons[i].strip().split("\t")
    parts[3] = int(parts[3])
    parts[4] = int(parts[4])
    unigenegffexons[i] = tuple(parts)


# calculate total transcript length
transcriptlength = sum([ int(gff[4])-int(gff[3])+1 for gff in unigenegffexons ])

# now determine orfs in the 3frame translations and annotate protein
# case 1: status cds
for protseq in frametrans.values():
    protseqlength = len(protseq)
    if protseq[0] == 'M' and protseq[-1] == '*' and protseq.count('*') == 1:
        #print 'cds', protseq[0:20], "...", protseq[-20:]
        sta = _obtain_unigene_gff_nt_start_coord(unigenegffexons,fstrand,0,transcriptlength,protseqlength,0)
        end = _obtain_unigene_gff_nt_end_coord(unigenegffexons,fstrand,0,transcriptlength,protseqlength,0)
        if not output_protein:
            print "cds\t%s\t%s\t0\t0\t%s\t%s" % (len(protseq),len(protseq),sta,end)
        else:
            print "cds\t%s\t%s\t0\t0\t%s\t%s\t%s" % (len(protseq),len(protseq),sta,end,protseq)
        sys.exit()

# case 2: status fragm(end) 
for protseq in frametrans.values():
    if protseq.count('*') == 0:
        if not output_protein:
            print "fragm\t%s\t%s\t0\t0\tna\tna" % (len(protseq),len(protseq))
        else:
            print "fragm\t%s\t%s\t0\t0\tna\tna\t%s" % (len(protseq),len(protseq),protseq)
        sys.exit()


# no return? do pattern matching with regular expression
# UTR5p is defined by having an *, then a >=0 non-M AA's and then the Methionine/TSS
# UTR3p is defined by existance (any symbol after the * of the orf)
# the difference between cds and utr is that cds's strings must either start with ^M or end with *$ 
patterns = {
  'full'  : compile("(?P<utrus>\*[^M^\*]*)(?P<orf>M[A-Z]{%s,}\*)(?P<utrds>[A-Z\*]*)" % (MIN_ORF_AA_LENGTH-1) ),
  'cds5p' : compile("(?P<utrus>\*[^M^\*]*)(?P<orf>M[A-Z]{%s,}\*)$" % (MIN_ORF_AA_LENGTH-1) ),
  'cds3p' : compile("^(?P<orf>M[A-Z]{%s,}\*)(?P<utrds>[A-Z\*]*)" % (MIN_ORF_AA_LENGTH-1) ),
  'utr5p' : compile("(?P<utrus>\*[^M^\*]*)(?P<orf>M[A-Z]{%s,})$" % (MIN_ORF_AA_LENGTH-1) ),
  'utr3p' : compile("^(?P<orf>[^M^\*]{1}[A-Z]{%s,}\*)(?P<utrds>[A-Z\*]*)" % (MIN_ORF_AA_LENGTH-2) ),
}

length2status = {}
protseq2line = {}
for status in ['utr5p','utr3p','cds5p','cds3p','full']:
    pat = patterns[status]
    for hdr,protseq in frametrans.iteritems():
        # hdr is of format _1, _2 or _3; minus 1 for python 0-based coord system
        readingframe = int(hdr[-1])-1
        protseqlength = len(protseq)

        for mobj in finditer(pat,protseq):
            #print status, "\t", len(protseq), mobj.start(), mobj.end(), "\torflength:",
            #print len(mobj.group('orf'))  #mobj.group(), mobj.groups()
            # write to dict; when key exists, it is overwitten!
            # that is why the loop over the statusses is ordered, and not patterns.keys()

            sta, end = "na", "na"
            if status == 'full':
                utr5paasize = mobj.start()+len(mobj.group('utrus'))
                utr3paasize = len(protseq)-mobj.start()-len(mobj.group('utrus'))-len(mobj.group('orf'))
                #utr3paasize = len(protseq) - ( mobj.start()+len(mobj.group('utrus'))+len(mobj.group('orf')) ) 

                sta = _obtain_unigene_gff_nt_start_coord(unigenegffexons,fstrand,readingframe,transcriptlength,protseqlength,utr5paasize)
                end = _obtain_unigene_gff_nt_end_coord(  unigenegffexons,fstrand,readingframe,transcriptlength,protseqlength,utr3paasize)
                #print mobj.group('utrus')
                #print mobj.group('orf')
                #print mobj.group('utrds')
            elif status in ['cds5p','utr5p']:
                utr5paasize = mobj.start()+len(mobj.group('utrus'))
                utr3paasize = 0
                sta = _obtain_unigene_gff_nt_start_coord(unigenegffexons,fstrand,readingframe,transcriptlength,protseqlength,utr5paasize)
                if status == 'cds5p':
                    end = _obtain_unigene_gff_nt_end_coord(  unigenegffexons,fstrand,readingframe,transcriptlength,protseqlength,utr3paasize)
            else:
                utr5paasize = 0
                utr3paasize = len(protseq) - len(mobj.group('orf'))
                end = _obtain_unigene_gff_nt_end_coord(  unigenegffexons,fstrand,readingframe,transcriptlength,protseqlength,utr3paasize)
                if status == 'cds3p':
                    sta = _obtain_unigene_gff_nt_start_coord(unigenegffexons,fstrand,readingframe,transcriptlength,protseqlength,utr5paasize)

            resultline = "%s\t%s\t%s\t%s\t%s\t%s\t%s" % (status, len(mobj.group('orf')),
                       len(protseq), utr5paasize, utr3paasize, sta, end
                       )
            length2status[ len(mobj.group('orf')) ] = resultline
            protseq2line[resultline] = mobj.group('orf')

# add empty noncoding status too in case no comfirmation is found
length2status[0] = 'nonc\t0\t0\tna\tna\tna\tna'
protseq2line[0] = ""

# and print the status line of the longest orf evidence
if print_all_frames:
    keys = length2status.keys()
    keys.sort()
    keys.reverse()
    for key in keys:
        line = length2status[key]
        if line[0:4] != "nonc":
            protseq = protseq2line[line]
        else:
            protseq = ""
        print "%s\t%s" % (line,protseq)
else:
    if not output_protein:
        print length2status[max(length2status.keys())]
    else:
        line = length2status[max(length2status.keys())]
        if line[0:4] != "nonc":
            protseq = protseq2line[line]
        else:
            protseq = ""
        print "%s\t%s" % (line,protseq)


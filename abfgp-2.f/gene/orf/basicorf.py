"""
"""

# Module metadata
__authors__ = "Ate van der Burgt"
__license__ = "MIT"

# Python Imports
import re
from sys import path as sysPath

# make shure imports from the parental gene package will work,
# even if it is not specified (yet) in sys.path
from os.path import dirname as osPathDirname, abspath as osPathAbspath, split as osPathSplit
gene_package_dir = osPathSplit(osPathDirname(osPathAbspath(__file__)))[0]
if gene_package_dir not in sysPath: sysPath.append(gene_package_dir)

# Gene Imports
from gene_exceptions import InproperlyAppliedArgument

# compiled pattern match of coordinates and id of orfs in getorf ourput file
pat_getorf_coords = re.compile("\[(\d+) \- (\d+)\]")
pat_getorf_id =     re.compile("_(\d+) \[\d+ \- \d+\]")


class BasicOrf:
    def __init__(self,fastaheader,fastasequence,
	inputgenomicsequence="",force_id=None,start=None,end=None):
        """
        parse data from
        """
        self.fastaheader = fastaheader
        self.id = force_id
        if force_id or force_id == 0:
            self.id = force_id
        else:
            try:
                self.id  = int(re.search(pat_getorf_id,fastaheader).groups()[0])
            except:
                pass
        self.header = "_".join(fastaheader.split(" ")[0].split("_")[0:-1])
        self.inputgenomicsequence = inputgenomicsequence
        # get start & stop coordinates from arguments
        # start and stop are GFF/GGB DNA coordinates
        if start != None and end != None:
            self.start   = start
            self.end     = end
        else:
            # get start & stop coordinates
            # all coordinates are based on the genomic slice
            match = re.search(pat_getorf_coords,fastaheader)
            start, end = match.groups()
            start, end = int(start), int(end)
            self.start   = start
            self.end     = end

        # make python coordinates as well, compared to GGB coordinates
        # startPY and stopPY are python string DNA coordinates
        self.startPY = start-1
        self.endPY   = end
        self.frame   = self.startPY % 3
        self.length  = self.endPY - self.startPY
        # now data on the protein slice
        self.protein_sequence = fastasequence
        self.protein_length   = len(fastasequence)
        self.protein_startPY  = (self.startPY-self.frame) / 3
        self.protein_endPY    = (self.endPY-self.frame) / 3
        self.protein_start    = self.protein_startPY+1
        self.protein_end      = self.protein_endPY
        self.orf_type = ("stop","stop")

    # end of function __init__


    def __str__(self):
        """ """
        return """<Orf:%s %s %s-%s>""" % (
            self.id,
            self.header,
            self.startPY,
            self.endPY,
            )

    # end of function __str__


    def getaa(self,abs_pos=None):
        """
        Return a Amino acid from the Orf protein sequence

        @type abs_pos:  (positive) integer
        @param abs_pos: desired absolute Amino acid position

        @rtype:  string
        @return: requested Amino acid, stop-codon (*) or empty string
        """
        if abs_pos or abs_pos == 0:
            # recalculate abs_pos to relative pos
            rel_pos = abs_pos - self.protein_startPY
            if rel_pos < 0:
                if abs_pos+1 == self.protein_startPY:
                    return "*"
                else:
                    return ""
            if rel_pos == self.protein_length:
                return "*"
            try:
                return self.protein_sequence[rel_pos]
            except:
                # hmm...non-existing coordinate!
                return ""
        else:
            return ""

    # end of function getaa


    def getntseq(self,abs_aa_pos=None,rel_aa_pos=None):
        """ """
        if abs_aa_pos or abs_aa_pos == 0:
            ntstart = abs_aa_pos*3 + self.frame
            return self.inputgenomicsequence[ntstart:ntstart+3].lower()
        elif rel_aa_pos or rel_aa_pos == 0:
            ntstart = (self.protein_startPY+rel_aa_pos)*3 + self.frame
            return self.inputgenomicsequence[ntstart:ntstart+3].lower()
        else:
            return ""

    # end of function getntseq


    def getntseqslice(self,abs_aa_startPY=0,abs_aa_endPY=0):
        """ """
        return "".join([ self.getntseq(abs_aa_pos=pos) for pos in range(abs_aa_startPY,abs_aa_endPY)] )

    # end of function getntseqslice


    def dnapos2aapos(self,pos):
        """ """
        return ( pos - self.frame ) / 3

    # end of function dnapos2aapos


    def aapos2dnapos(self,pos):
        """
        ALIAS FOR proteinpos2dnapos
        TODO: delete on of these....
        """
        return pos*3 + self.frame

    # end of function aapos2dnapos


    def proteinpos2dnapos(self,pos):
        """
        ALIAS FOR aapos2dnapos
        TODO: delete on of these....
        """
        return pos*3 + self.frame

    # end of function proteinpos2dnapos


    def getaas(self,abs_pos_start=None,abs_pos_end=None):
        """
        """
        if abs_pos_start and abs_pos_end:
            # recalculate abs_pos to relative pos
            rel_pos_start = abs_pos_start - self.protein_startPY
            rel_pos_end   = abs_pos_end   - self.protein_startPY
            try:    return self.protein_sequence[rel_pos_start:rel_pos_end]
            except: return ""
        elif (abs_pos_start or abs_pos_start == 0) and not abs_pos_end:
            # no end pos
            rel_pos_start = abs_pos_start - self.protein_startPY
            try:    return self.protein_sequence[rel_pos_start:]
            except: return ""
        elif (abs_pos_end or abs_pos_end == 0) and not abs_pos_start:
            # no start pos
            rel_pos_end   = abs_pos_end   - self.protein_startPY
            try:    return self.protein_sequence[:rel_pos_end]
            except: return ""
        else:
            return ""

    # end of function getaas


    def nucleotidesequence(self,extra_start=0,extra_end=0):
        """
        get nucleotide sequence of this orf
        extra_start may be in [0,1,2,3] to recover the leading STOP codon
        extra_end   may be in [0,1,2,3,4]    to recover the trailing STOP codon + 1 nt
        the latter is for recovering splice sites in the STOP codon!
        """

        #if extra_start not in [0,1,2,3]:
        #    message = "extra_start (%s) not in [0,1,2,3]" %  extra_start
        #    raise InproperlyAppliedArgument, message
        #if extra_end not in [0,1,2,3,4]:
        #    message = "extra_end (%s) not in [0,1,2,3,4]" %  extra_end
        #    raise InproperlyAppliedArgument, message

        if not self.inputgenomicsequence:
            return ""
        else:
            # okay, return the genomic sequence slice
            start, end = ( self.startPY-extra_start, self.endPY+extra_end )
            # check for coordinates out of range of inputgenomicsequence
            # and correct those with `n`-nucleotides
            prefix, suffix = "", ""
            if start < 0 and end > len(self.inputgenomicsequence):
                prefix = "n"*abs(start)
                suffix = "n"*(end - len(self.inputgenomicsequence))
                start  = 0
            elif start < 0:
                prefix = "n"*abs(start)
                start  = 0
            elif end > len(self.inputgenomicsequence):
                suffix = "n"*(end - len(self.inputgenomicsequence))
            else:
                pass
            # and return the requested sequence
            return "".join([ prefix,self.inputgenomicsequence[start:end],suffix ])

    # end of function nucleotidesequence


    def printproteinanddna(self,_linesize=40):
        """ """
        ntseq = self.nucleotidesequence()
        pos=0
        aoa = []
        coordNT = "    %s" % self.startPY
        coordAA = "    %s" % self.protein_startPY
        coordNT = coordNT[-5:] + " "
        coordAA = coordAA[-5:] + " "
        aoa.append([coordAA,coordNT])
        aoapos=1
        for pos in range(0,self.protein_length):
            aoa.append( [ " %s " % self.protein_sequence[pos], ntseq[pos*3:(pos*3)+3] ] )
            aoapos += 1
            if (pos+1) % _linesize == 0:
                coordNT = "    %s" % ( self.startPY + (pos * 3 ) )
                coordAA = "    %s" % ( self.protein_startPY + pos )
                coordNT = " " + coordNT[-5:]
                coordAA = " " + coordAA[-5:]
                aoa.append([coordAA,coordNT])
                coordNT = "    %s" % ( self.startPY + (pos * 3 ) + 1)
                coordAA = "    %s" % ( self.protein_startPY + pos + 1)
                coordNT = coordNT[-5:] + " "
                coordAA = coordAA[-5:] + " "
                aoa.append([coordAA,coordNT])
                aoapos+=2
        else:
            if self.protein_length % _linesize != 0:
                coordNT = "    %s" % self.endPY
                coordAA = "    %s" % self.protein_endPY
                coordNT = " " + coordNT[-5:]
                coordAA = " " + coordAA[-5:]
                aoa.append([coordAA,coordNT])
        # and print!
        print self._print_aoa(aoa,_linesize=_linesize+2)

    # end of function printintron

    def _print_aoa(self,aoa,_linesize=50):
        """
        """
        rows = []
        if aoa:
            totallength = len(aoa)
            batches = range(0,totallength,_linesize)
            for batch_start in batches:
                batch_end = min([batch_start+_linesize,totallength])
                subrows = [""]*len(aoa[0])
                for item in aoa[batch_start:batch_end]:
                    for i in range(0,len(item)):
                        subrows[i]+=item[i]
                rows.extend(subrows)
                rows.extend("")
        # and return
        return "\n".join(rows)

    # end of function _print_aoa

    def has_dna_n_symbols(self):
        """ 
        Does Orf's DNA sequence contains N nucleotides (i.e. X-symbols in PROTEIN sequence)?
        
        @rtype:  Boolean
        @return: True or False
        """
        if self.nucleotidesequence().lower().find("n") >= 0:
            return True 
        else:
            return False

    # end of function has_dna_n_symbols 


    def has_protein_x_symbols(self):
        """
        Does Orf's PROTEIN sequence contains X symbols (i.e N-nucleotides in DNA sequence)?

        @rtype:  Boolean
        @return: True or False
        """
        if self.protein_sequence.find("X") >= -1:
            return True
        else:
            return False

    # end of function has_protein_x_symbols 


    def tofasta(self,header=""):
        """
        """
        if not header:
            header = self.fastaheader
        return ">%s\n%s" % (header,self.protein_sequence)

    # end of function tofasta


    def togff(self,fmethod='orf',fsource='GETORF',gclass='Orf',fstrand='+',fscore='.',fref='',one_group_per_frame=True):
        """
        Return gff tuple for current orf
        """
        if not fref:
            if self.fastaheader:
                fref = self.fastaheader.split(" ")[0].split("_")[0]
            else:
                raise "'fref' identifier MUST be specified"
        if one_group_per_frame:
            # identical gname for all orfs in the same stand
            # this enables nicer view in GGB
            gname = "frame%s" % self.frame
        else:
            gname = "%s-f%s-%s-(id%s)" % (fref,self.frame,self.start,self.id)
            gname = "%s-f%s-%s-(id%s)" % (fref,self.frame,self.start,self.id)

        # layout fscore
	fscore = self._get_gff_fscore()

        # return gff tuple
        return ( fref, fsource, fmethod, self.start, self.end,
		 fscore, fstrand, self.frame,
		 "%s %s; orfid %s" % (gclass,gname,self.id) )

    # end of function togff

    def _get_gff_fscore(self):
	return "."
    # end of function _get_gff_fscore

# end of class BasicOrf


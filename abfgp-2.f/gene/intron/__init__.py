"""
"""
# Module metadata
__authors__ = "Ate van der Burgt"
__license__ = "MIT"

# Python Imports
from os import popen as osPopen
from os import popen2 as osPopen2
from copy import deepcopy
from sys import path as sysPath

# make shure imports from the parental gene package will work,
# even if it is not specified (yet) in sys.path
from os.path import dirname as osPathDirname, abspath as osPathAbspath, split as osPathSplit
gene_package_dir = osPathSplit(osPathDirname(osPathAbspath(__file__)))[0]
if gene_package_dir not in sysPath: sysPath.append(gene_package_dir)

# Gene Imports
from gene_gff import BasicGFF, GffWithSequenceFunctionality
from splicesite import get_shared_nucleotides_at_splicesite
from donor import ProjectedSpliceDonor
from acceptor import ProjectedSpliceAcceptor
from polypyrimidinetract import PolyPirimidineTract
from branchpoint import BranchPoint 
from validators import *

# Other Imports
import dna2prot


# Import Global variables
from settings.executables import (
    EXECUTABE_SCAN_BRANCHPOINT,
    EXECUTABLE_SFM, # EXECUTABE_SCAN_BRANCHPOINT requires EXECUTABLE_SFM
    EXECUTABE_SCAN_POLYPYRIMIDINES,
    PYTHON_PATH,
    TCODE_MAX_NONCODING,
    TCODE_MIN_CODING,
    )
from settings.splicesites import OPTIMAL_BRACNHPOINT_TO_ACCEPTOR_DISTANCE

def parsegfftxt(gfftxt,offset=0):
    """
    Parse a gff text block into a list of tracks (as tuples)

    @type  gfftxt: string 
    @param gfftxt: gff text block

    @type  offset: integer
    @param offset: correct the positions in the gff tracks by this offset

    @rtype:  list
    @return: list with gff tracks as tuples (9 elements)
    """ 
    tracks = []
    for line in gfftxt.strip().split("\n"):
        if not line.strip(): continue
        tracks.append( _parse_line_2_track(line,offset=offset) )
    # return gff tracks list
    return tracks 

# end of function parsegfftxt


def _parse_line_2_track(line,offset=0,separator="\t"):
    """
    Parse a single gff file into a list of tracks (as tuples)
        
    @type  line: string 
    @param line: gff text line (of 8 columns!)
        
    @type  offset: integer
    @param offset: correct the positions in the gff tracks by this offset
        
    @type  separator: string
    @param separator: string which separates the gff line (default tab)
        
    @rtype:  tuple
    @return: single gff line as a tuple (9 elements)
    """ 
    tracks = []
    gff = line.strip().split(separator)
    gff[3] = int(gff[3])+offset
    gff[4] = int(gff[4])+offset
    if gff[5] != ".":
        if gff[5].find('.') > 0:
            gff[5] = float(gff[5])
        else:
            gff[5] = int(gff[5])
    return tuple(gff)

# end of function _parse_line_2_track


class Intron(BasicGFF):
    def __init__(self,donor,acceptor,shared_nts,gff={}):
        # input validation
        IsDonor(donor)
        IsAcceptor(acceptor)
        # initialization
        BasicGFF.__init__(self)
        self._gff.update(gff)
        self.donor      = donor
        self.acceptor   = acceptor
        # init by splice site objects
        self._init_by_splicesites()

        # forward compatibility with IntronConnectingOrfs
        self.nt_from_stop_donor = None
        self.nt_from_stop_acceptor = None

        # set shared nt/aa data
        self.shared_nts = shared_nts 
        self.shared_aa  = ""
        if self.shared_nts:
            self.shared_aa = dna2prot.dna2protein(self.shared_nts)

    # end of function __init__


    def _init_by_splicesites(self):
        """ """
        self.length     = self.acceptor.pos-self.donor.pos
        self.start      = self.donor.pos
        self.end        = self.acceptor.pos
        # check/set phase of the intron
        if self.donor.phase and self.acceptor.phase:
            if self.donor.phase == self.acceptor.phase:
                self.phase = self.donor.phase
            else:
                raise IncompatibleSpliceSitePhases
        elif self.donor.phase == 0 and self.acceptor.phase == 0:
            self.phase = 0
        else:
            # hmm... not recommended, but no phase specified on donor/acceptor
            self.phase = None
        # get aa positions of donor and acceptor site
        if self.phase != None:
            self.donor_aa_pos          = (self.donor.pos-self.donor.phase) / 3
            self.acceptor_aa_pos       = (self.acceptor.pos-self.acceptor.phase) / 3
        else:
            self.donor_aa_pos          = self.donor.pos / 3
            self.acceptor_aa_pos       = self.acceptor.pos / 3

    # end of function _init_by_splicesites


    def deepcopy(self):
        """ """
        return self.__class__(self.donor,self.acceptor,self.shared_nts,gff=self._gff)

    # end of function deepcopy

    def __str__(self):
        str_nt_from_stop_donor =    str(-self.nt_from_stop_donor)
        str_nt_from_stop_acceptor = str(-self.nt_from_stop_acceptor)
        if self.nt_from_stop_donor <= 0:
            str_nt_from_stop_donor = "+"+str_nt_from_stop_donor
        if self.nt_from_stop_acceptor <= 0:
            str_nt_from_stop_acceptor = "+"+str_nt_from_stop_acceptor
        return "<%s %s [%s] | GT .. [%s/%s] length=%s .. AG | [%s]  %s >" % (
            self.__class__.__name__,
            self.donor,
            str_nt_from_stop_donor,
            self.shared_nts,
            self.shared_aa,
            self.length, 
            str_nt_from_stop_acceptor,
            self.acceptor,
            )

    # end of function __str__

    def coords(self):
        """ """
        return (self.donor.pos,self.acceptor.pos)

    # end of function coords

    def is_3n_intron(self):
        """ Is the length of this intron devidable by 3 """
        return self.length % 3 == 0

    # end of function is_3n_intron

    def togff(self,gff={}):
        """ Overwrites BasicGFF.togff() """
        if not gff.has_key('gname') and not self._gff.has_key('gname'):
            gff['gname'] = "%s-%s" % (self.start, self.end)
        return BasicGFF.togff(self,gff=gff)

    # end of function togff

# end of class Intron


class IntronConnectingOrfs(Intron,GffWithSequenceFunctionality):
    def __init__(self,donor,acceptor,shared_nts,orfDonor,orfAcceptor,
        binary_entropy_donor=None,binary_entropy_acceptor=None,gff={}):
        """
        @type donor:    SpliceDonor object
        @param donor:   fully instantiated SpliceDonor object (on orfDonor)

        @type acceptor:  SpliceAcceptor object
        @param acceptor: fully instantiated SpliceAcceptor object (on orfAcceptor)

        @type shared_nts:   string
        @param shared_nts:  3-nt string or empty string

        @type orfDonor:     Orf object
        @param orfDonor:    Orf object that contains the SpliceDonor site

        @type orfAcceptor:  Orf object
        @param orfAcceptor: Orf object that contains the SpliceAcceptor site

        """
        # input validation
        IsOrf(orfDonor)
        IsOrf(orfAcceptor)

        # if shared_nts not given -> get them here
        # TODO: REMOVE shared_nts FUNCTION ARGUMENT; superfluously!
        if shared_nts == None:
            shared_nts = get_shared_nucleotides_at_splicesite(
                            orfAcceptor, orfDonor,
                            acceptor, donor )

        # initialise basal Intron class
        Intron.__init__(self,donor,acceptor,shared_nts,gff=gff)
        self.orfDonor              = orfDonor
        self.orfAcceptor           = orfAcceptor

        # init by splice site objects
        self._init_by_splicesites()

        self.binary_entropy_donor   = binary_entropy_donor
        self.binary_entropy_acceptor= binary_entropy_acceptor
        # optional attributes that can be obtained/filled
        # with active function calls
        self.branchpoint= None
        self.ppt5p      = None
        self.ppt3p      = None

    # end of function __init__


    def _init_by_splicesites(self):
        """ """
        Intron._init_by_splicesites(self)
        if hasattr(self,"orfDonor"):
            self.nt_from_stop_donor    = self.orfDonor.endPY - self.donor.pos
            self.nt_from_stop_acceptor = self.acceptor.pos - self.orfAcceptor.startPY

    # end of function _init_by_splicesites

    def deepcopy(self):
        """ """
        return self.__class__(self.donor,self.acceptor,self.shared_nts,
                self.orfDonor,self.orfAcceptor,
                binary_entropy_donor=self.binary_entropy_donor,
                binary_entropy_acceptor=self.binary_entropy_acceptor,
                gff=self._gff)

    # end of function deepcopy


    def is_stopless_3n_intron(self):
        """ Is this a stopless 3n intron? """
        if self.is_3n_intron() and self.orfDonor.id == self.orfAcceptor.id:
            return True
        else:
            return False

    # end of function is_stopless_3n_intron


    def togff(self,gff={},annotated_intron=False):
        """ """
        if not annotated_intron:
            return Intron.togff(self,gff=gff)
        else:
            gfflines = [ Intron.togff(self,gff=gff) ]
            gclass, gname = gfflines[0][-1].split("; ")[0].split(" ",1)
            if gname[0]+gname[-1] in ["''",'""']:
                gname = gname[1:-1]
            gff['gclass'] = gclass
            gff['gname'] = gname
            gfflines.append( self.donor.togff(gff=gff) )
            if self.ppt5p:
                gfflines.append( self.ppt5p.togff(gff=gff) )
            if self.branchpoint:
                gfflines.append( self.branchpoint.togff(gff=gff) )
            if self.ppt3p:
                gfflines.append( self.ppt3p.togff(gff=gff) )
            gfflines.append( self.acceptor.togff(gff=gff) )
            return gfflines

    # end of function togff


    def assign_bp_and_ppts(self):
        """ """
        bpdist = max(OPTIMAL_BRACNHPOINT_TO_ACCEPTOR_DISTANCE)
        bplength = 5
        bps = self.scan_branchpoint()
        ppts = self.scan_polypyrimidines()
        combis = []

        if ppts:
            dummy_ppt = list(ppts[0])
            dummy_ppt[3]=0
            dummy_ppt[4]=0
            dummy_ppt[5]=0
        else:
            dummy_ppt = ['intron', 'PPTregex', 'PPTregex', 0, 0, 0, '+', '.']

        # make dummy branchpoint (for absence of BP, only ppts)
        dummy_bp = ( None, None, None, self.acceptor.pos-bpdist-bplength, self.acceptor.pos-bpdist, -2 )
        # quit searhing after the first 10 branchpoints (should be located near the AG site!)
        bps = bps[0:10]
        bps.append( dummy_bp )
        for bp in bps:
            score = bp[5]
            # branchpoints are optimally ~15-20nt spaced from the acceptor site
            ####spacing = abs( bpdist - (self.acceptor.pos - bp[4]) )
            spacing = min([abs(offset - (self.acceptor.pos - bp[4])) for offset in OPTIMAL_BRACNHPOINT_TO_ACCEPTOR_DISTANCE])
            if bp[5] == -2:
                # this is the DUMMY branchpoint
                pptlen = 0
                ppt5p = None
                ppt3p = None
                for ppt in ppts:
                    # allow for slight overlap;
                    # [CT]T[AG]A[CT] easily overlaps with [CTG]+
                    if ppt[3] >= bp[4]-1:
                        ppt3p = ppt
                        score += 1
                        pptlen+=ppt[4]-ppt[3]+1
                        break
                for pos in range(len(ppts)-1,-1,-1):
                    ppt = ppts[pos]
                    # allow for slight overlap;
                    # [CT]T[AG]A[CT] easily overlaps with [CTG]+
                    if ppt[4]-2 <= bp[3]:
                        ppt5p = ppt
                        score += 1
                        pptlen+=ppt[4]-ppt[3]+1
                        break
                if not ppt5p and not ppt3p:
                    # no ppts available around putative dummy bp
                    continue
                else:
                    # correct score and set bp to None
                    bp = None
                    score += 1 # counts as the weakest possible branchpoint
                    # main ordering criterion is ppt-tract-length corrected for spacing of bp
                    combis.append( ( (pptlen-spacing+(score*10),score,pptlen,spacing), ppt5p, bp, ppt3p ) )
                    # break out of the for-loop
                    break

            # if here, a normal (non-dummy) branchpoint
            ppt5p = None
            ppt3p = None
            for pptA in ppts:
                # allow for slight overlap;
                # [CT]T[AG]A[CT] easily overlaps with [CTG]+
                if pptA[3] >= bp[4]-1:
                    ppt3p = pptA
                else:
                    ppt3p = None
                for pptB in ppts:
                    # allow for slight overlap;
                    # [CT]T[AG]A[CT] easily overlaps with [CTG]+
                    if pptB[4]-2 <= bp[3]:
                        ppt5p = pptB
                    else:
                        ppt5p = None

                    if not ppt3p and bp:
                        # search for short ppt3p overlapping with BP
                        offset = bp[3] - self.start - 1 + 4
                        minibpseq = self.dnasequence()[offset:offset+5]
                        if len(minibpseq) == 5 and\
                        minibpseq.upper().count("T") >= 3 and\
                        minibpseq.upper().count("A") == 0:
                            ppt3p = deepcopy(dummy_ppt)
                            ppt3p[3] = bp[4]
                            ppt3p[4] = bp[4]+4
                            ppt3p[5] = 5
                            ppt3p = tuple(ppt3p)

                    # calculate score of this combination
                    score  = bp[5]
                    pptlen = 0
                    if ppt3p:
                        pptlen+=ppt3p[4]-ppt3p[3]+1
                        score+=1
                    if ppt5p: 
                        pptlen+=ppt5p[4]-ppt5p[3]+1
                        score+=1

                    # check for presence of (poor) branchpoint AND both pptracts
                    # reward this by increasing the score (1) with +1
                    if bp and bp[5] == -1 and ppt5p and ppt3p:
                        score+=1

                    # append to found combinations;
                    ## main ordering criterion is score (presence/absence)
                    ## 2th ordering criterion is ppt-tract-length corrected for spacing of bp
                    ##combis.append( ( (score,pptlen-spacing,pptlen,spacing), ppt5p, bp, ppt3p ) )

                    # main ordering criterion is ppt-tract-length corrected for spacing of bp
                    combi = ( (pptlen-spacing+(score*10),score,pptlen,spacing), ppt5p, bp, ppt3p )
                    if combi not in combis: combis.append( combi )

        # if no combinations are found -> no bp/ppt data at all!
        if not combis: return False

        # order found combinations
        combis.sort()
        combis.reverse()

        scoretuple, ppt5p, bp, ppt3p = combis[0]
        if ppt5p:
            _gff = { 'fsource': ppt5p[1], 'fmethod': ppt5p[2] }
            self.ppt5p = PolyPirimidineTract(ppt5p[3]-1,ppt5p[4],gff=_gff)
        if bp:
            _gff = { 'fsource': bp[1], 'fmethod': bp[2], 'fscore': bp[5] }
            self.branchpoint = BranchPoint(bp[3]-1,bp[4],gff=_gff)
        if ppt3p:
            _gff = { 'fsource': ppt3p[1], 'fmethod': ppt3p[2] }
            self.ppt3p = PolyPirimidineTract(ppt3p[3]-1,ppt3p[4],gff=_gff)
        # return status message True
        return True
        
    # end of function assign_bp_and_ppts


    def get_branchpoint_nt_distance(self):
        """ Get nt distance between branchpoint and acceptor (yurAy .. aG  A->G) """
        if not self.branchpoint: self.assign_bp_and_ppts()
        if not self.branchpoint:
            return None
        else:
            return self.acceptor.pos - self.branchpoint.end + 2

    # end of function get_branchpoint_nt_distance


    def WAS_scan_branchpoint(self):
        """ """
        intronseq  = self.dnasequence()
        co = osPopen(""" echo ">intron\n%s\n" | %s %s """ % (
            intronseq,EXECUTABE_SCAN_BRANCHPOINT, EXECUTABLE_SFM ) )
        gfflines = parsegfftxt( co.read(), offset=self.start )
        co.close()
        orderable = [ ( int(line[3]), line ) for line in gfflines ]
        orderable.sort()
        orderable.reverse()
        return [ gffline for pos,gffline in orderable ]

    # end of function scan_branchpoint


    def scan_branchpoint(self):
        """ """
        intronseq  = self.dnasequence()
        command = """%s -e %s""" % (EXECUTABE_SCAN_BRANCHPOINT, EXECUTABLE_SFM)
        ci,co = osPopen2(command)
        ci.write(">intron\n%s\n" % intronseq)
        ci.close()
        gfflines = parsegfftxt( co.read(), offset=self.start )
        co.close()
        orderable = [ ( int(line[3]), line ) for line in gfflines ]
        orderable.sort()
        orderable.reverse()
        return [ gffline for pos,gffline in orderable ]

    # end of function scan_branchpoint    
    
    def scan_polypyrimidines(self):
        """ """
        intronseq  = self.dnasequence()
        co = osPopen(""" echo ">intron\n%s\n" | %s %s """ % (
            intronseq,PYTHON_PATH,EXECUTABE_SCAN_POLYPYRIMIDINES ) )
        gfflines = parsegfftxt( co.read(), offset=self.start )
        co.close()
        return gfflines

    # end of function scan_polypyrimidines


    def get_donor_nt_from_stop(self):
        """ """
        return self.orfDonor.endPY - self.donor.pos
    # end of function get_donor_nt_from_stop


    def get_acceptor_nt_from_stop(self):
        """ """
        return self.acceptor.pos - self.orfAcceptor.startPY
    # end of function get_acceptor_nt_from_stop


    def dnasequence(self):
        """ """
        return self.orfDonor.inputgenomicsequence[self.start:self.end]
    # end of function dnasequence


    def barcode(self):
        """
        A unique tuple with integers describing this object
        """
        return ( self.nt_from_stop_donor, self.orfDonor.id, self.start, self.phase, self.end, self.orfAcceptor.id, self.nt_from_stop_acceptor )

    # end of function barcode


    def tcode(self):
        """ """
        # doesn't matter if we obtain it from orfDonor or orfAcceptor...
        return self.orfDonor.average_tcode_score_of_range(
                    self.orfDonor._RAW_TCODE_DATA,
                    self.start,self.end
                    )

    # end of function tcode


    def is_tcode_coding(self):
        """ """
        return self.tcode() >= TCODE_MIN_CODING

    # end of function is_tcode_coding


    def is_tcode_noncoding(self):
        """ """
        return self.tcode() < TCODE_MAX_NONCODING

    # end of function is_code_noncoding


    def printintron(self,offset=20):
        """
        """
        aoa = []
        for pos in range(self.donor_aa_pos-offset,self.donor_aa_pos):
            item = []
            item.append( " %s " % self.orfDonor.getaa(pos) )
            item.append( self.orfDonor.getntseq(abs_aa_pos=pos) )
            item.append( self.orfDonor.getntseq(abs_aa_pos=pos) )
            item.append( "   " )
            aoa.append(item)
        if self.phase == 0:
            posA = self.acceptor_aa_pos-2
            posB = self.acceptor_aa_pos-1
            posC = self.acceptor_aa_pos
            layout = self.orfDonor.getntseq(self.donor_aa_pos)[0:2].upper() + self.orfDonor.getntseq(self.donor_aa_pos)[-1]  
            aoa.append( [ "-0-", "-0-", layout, "   " ] )
            aoa.append( [ "-0-", "-0-",  self.orfDonor.getntseq(self.donor_aa_pos+1), self.orfAcceptor.getntseq(posA) ] )
            layout = self.orfAcceptor.getntseq(self.acceptor_aa_pos)[0] + self.orfAcceptor.getntseq(posB)[1:3].upper() 
            aoa.append( [ "-0-", "-0-",  "   ", layout ] )
        elif self.phase == 1:
            posA = self.acceptor_aa_pos-1
            posB = self.acceptor_aa_pos
            posC = self.acceptor_aa_pos+1
            layout = self.orfDonor.getntseq(self.donor_aa_pos)[0] + self.orfDonor.getntseq(self.donor_aa_pos)[1:3].upper()
            aoa.append( [ "-1-", "-1-", layout, "   " ] ) 
            layout = self.orfAcceptor.getntseq(posA)[0:2] + self.orfAcceptor.getntseq(posA)[-1].upper() 
            aoa.append( [ " ? ", "%s%s" % ( 
                    self.orfDonor.getntseq(abs_aa_pos=self.donor_aa_pos)[0].upper(),
                    self.orfAcceptor.getntseq(posB)[1:3] 
                    ), 
                    self.orfDonor.getntseq(self.donor_aa_pos+1),
                    layout ] )
            layout =  self.orfAcceptor.getntseq(posB)[0].upper() + self.orfAcceptor.getntseq(posB)[1:3] 
            aoa.append( [ "-1-", "-1-", "   ", layout ] )

        elif self.phase == 2:
            posA = self.acceptor_aa_pos-1
            posB = self.acceptor_aa_pos
            posC = self.acceptor_aa_pos+1
            layout = self.orfDonor.getntseq(self.donor_aa_pos)[0:2] + self.orfDonor.getntseq(self.donor_aa_pos)[-1].upper() 
            aoa.append( [ "-2-", "-2-", layout, "   " ] )
            layout = self.orfDonor.getntseq(abs_aa_pos=self.donor_aa_pos+1)[0].upper() + self.orfDonor.getntseq(abs_aa_pos=self.donor_aa_pos+1)[1:3]  
            aoa.append( [ "-2-", "%s%s" % ( 
                    self.orfDonor.getntseq(abs_aa_pos=self.donor_aa_pos)[0:2].upper(),
                    self.orfAcceptor.getntseq(posB)[-1] 
                    ),
                    layout,
                    self.orfAcceptor.getntseq(posA) ] )
            layout = self.orfAcceptor.getntseq(posB)[0:2].upper() + self.orfAcceptor.getntseq(posB)[-1]
            aoa.append( [ "-2-", "-2-", "   ", layout ] )

        else:
            raise "UNEXPECTED PHASE!!"


        for pos in range(posC,self.acceptor_aa_pos+offset+1):
            item = []
            item.append( " %s " % self.orfAcceptor.getaa(pos) )
            item.append( self.orfAcceptor.getntseq(abs_aa_pos=pos) )
            item.append( "   " )
            item.append( self.orfAcceptor.getntseq(abs_aa_pos=pos) )
            aoa.append(item)


        # and print!
        print self._print_aoa(aoa,_linesize=50)

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


# end of class IntronConnectingOrfs


class IntronObtainedByAlignment(IntronConnectingOrfs):
    """ """
    def __init__(self,distance=None,query=None,sbjct=None,*args,**kwargs):
        """ """
        IntronConnectingOrfs.__init__(self,*args,**kwargs)
        # nt distance between aligned splice sites from two species
        # query & sbjct are attributes to easily check object's origine
        self.distance = distance
        self.query    = query
        self.sbjct    = sbjct
    # end of function __init__

# end of class IntronObtainedByAlignment


class IntronObtainedByProjection(IntronConnectingOrfs):
    """ """
    def __init__(self,distance=None,query=None,sbjct=None,*args,**kwargs):
        """ """
        IntronConnectingOrfs.__init__(self,*args,**kwargs)
        # nt distance between aligned splice sites from two species
        # query & sbjct are attributes to easily check object's origine
        self.distance = distance
        self.query    = query
        self.sbjct    = sbjct
    # end of function __init__

# end of class IntronObtainedByProjection


class ProjectedIntronConnectingOrfs:
    """
    """
    def __init__(self,orfProjection,projected_donor_pos,projected_acceptor_pos,gff={}):
        """
        """
        self.orfProjection          = orfProjection
        self.projected_donor_pos    = projected_donor_pos
        self.projected_acceptor_pos = projected_acceptor_pos
        self.distance               = self.projected_acceptor_pos -\
                                      self.projected_donor_pos
        self.projected_introns      = []


    def barcode(self):
        """
        """
        return ( self.orfProjection.id, self.projected_donor_pos, self.projected_acceptor_pos )

    # end of function barcode


    def add_projected_intron(self,intron_connecting_orfs):
        """
        """
        if intron_connecting_orfs.__class__.__name__ == 'IntronConnectingOrfs':
            self.projected_introns.append( intron_connecting_orfs )
        else:
            message = "argument is not a `IntronConnectingOrfs` class object (%s)" % (
                intron_connecting_orfs.__class__.__name__ )
            raise InproperlyAppliedArgument, message

    # end of function add_projected_intron

    def __str__(self):
        if self.projected_introns:
            phase = self.projected_introns[0].phase
        else:
            phase = None 

        return "<projectedintron %s (%1.2f %1.2f) [%s %sx phase %s] (%1.2f %1.2f) %s>" % (
            self.projected_donor_pos,
            sum([intron.binary_entropy_donor for intron in self.projected_introns]),
            sum([intron.donor.pssm_score for intron in self.projected_introns]),
            self.distance,
            len(self.projected_introns),
            phase,
            sum([intron.acceptor.pssm_score for intron in self.projected_introns]),
            sum([intron.binary_entropy_acceptor for intron in self.projected_introns]),
            self.projected_acceptor_pos,
            )

    # end of function __str__

    def get_projected_donor_site(self):
        """
        """
        example_donor   = self.projected_introns[0].donor
        cases           = len(self.projected_introns)
        total_entropy   = 0.0
        total_pssm_score= 0.0
        for intron in self.projected_introns:
            total_entropy   += intron.binary_entropy_donor
            total_pssm_score+= intron.donor.pssm_score
        # and return a ProjectedSpliceSite
        return ProjectedSpliceDonor(self.projected_donor_pos,"XX",
                example_donor,distance=self.distance,entropy=total_entropy,
                pssm_score=total_pssm_score,cases=cases)

    # end of function get_projected_donor_site


    def get_projected_acceptor_site(self):
        """
        """
        example_acceptor= self.projected_introns[0].acceptor
        cases           = len(self.projected_introns)
        total_entropy   = 0.0
        total_pssm_score= 0.0
        for intron in self.projected_introns:
            total_entropy   += intron.binary_entropy_acceptor
            total_pssm_score+= intron.acceptor.pssm_score
        # and return a ProjectedSpliceSite
        return ProjectedSpliceAcceptor(self.projected_acceptor_pos,"XX",
                example_acceptor,distance=self.distance,entropy=total_entropy,
                pssm_score=total_pssm_score,cases=cases)

    # end of function get_projected_acceptor_site


    def total_score(self):
        """ """
        return self.get_projected_donor_site().total_score() + self.get_projected_acceptor_site().total_score()

    # end of function total_score

    def togff(self):
        """ ProjectedIntronConnectingOrfs has no gff() function (yet)... """
        pass

    # end of function togff

# end of class ProjectedIntronConnectingOrfs


def _filter_intron_list(introns,filter_by='length',criterion=35,operator=None):
    """
    """
    filter_options = [
        None,           False,
        'length',       'total_pssm',
        'donor_pos',    'donor_pssm',
        'acceptor_pos', 'acceptor_pssm',
        ]
    if filter_by not in filter_options:
        # inproper filter_by argument
        message = "filter_by (%s) not in %s" % (
            filter_by, filter_options )
        raise InproperlyAppliedArgument, message

    operator_options = ['>','>=','==','!=','<','<=','%']

    if operator not in operator_options:
        # inproper operator argument
        message = "operator (%s) not in %s" % (
            operator, operator_options )
        raise InproperlyAppliedArgument, message

    filtered_introns = []
    if filter_by == 'length':
        for i in introns:
            if eval('i.length %s %s' % ( operator, criterion ) ):
                filtered_introns.append(i)
    elif filter_by == 'total_pssm':
        for i in introns:
            if eval('i.acceptor.pssm_score + i.donor.pssm_score %s %s' % ( operator, criterion ) ):
                filtered_introns.append(i)
    elif filter_by == 'donor_pos':
        for i in introns:
            if eval('i.donor.pos %s %s' % ( operator, criterion ) ):
                filtered_introns.append(i)
    elif filter_by == 'acceptor_pos':
        for i in introns:
            if eval('i.acceptor.pos %s %s' % ( operator, criterion ) ):
                filtered_introns.append(i)
    elif filter_by == 'donor_pssm':
        for i in introns:
            if eval('i.donor.pssm_score %s %s' % ( operator, criterion ) ):
                filtered_introns.append(i)
    elif filter_by == 'acceptor_pssm':
        for i in introns:
            if eval('i.acceptor.pssm_score %s %s' % ( operator, criterion ) ):
                filtered_introns.append(i)
    else:
        # no filtering requested
        filtered_introns = introns

    return filtered_introns

# end of function _filter_intron_list


def _order_intron_list(introns,order_by='length'):
    """
    """
    order_options = [
        None,           False,
        'length',       'total_pssm',
        'donor_pos',    'donor_pssm',
        'acceptor_pos', 'acceptor_pssm',
        ]
    if order_by not in order_options:
        # inproper order_by argument
        message = "order_by (%s) not in %s" % (
            order_by, order_options )
        raise InproperlyAppliedArgument, message

    # temporary list for ordering
    tmp = []

    if order_by == 'length':
        tmp = [ (i.length, -i.donor.pos, i) for i in introns ]
        tmp.sort()
        return [ i for (a,b,i) in tmp ]
    elif order_by == 'total_pssm':
        try:
            tmp = [ (-i.acceptor.pssm_score-i.donor.pssm_score, i.length, i) for i in introns ]
            tmp.sort()
            return [ i for (a,b,i) in tmp ]
        except:
            # hmm... pssm_scores not available in the splice site objects
            return _order_intron_list(introns,order_by='length')
    elif order_by == 'donor_pos':
        tmp = [ (-i.donor.pos, i.acceptor.pos, i) for i in introns ]
        tmp.sort()
        return [ i for (a,b,i) in tmp ]
    elif order_by == 'acceptor_pos':
        tmp = [ (i.acceptor.pos, -i.donor.pos, i) for i in introns ]
        tmp.sort()
        return [ i for (a,b,i) in tmp ]
    elif order_by == 'donor_pssm':
        try:
            tmp = [ (-i.donor.pssm_score, i.length, i) for i in introns ]
            tmp.sort()
            return [ i for (a,b,i) in tmp ]
        except:
            # hmm... pssm_scores not available in the splice site objects
            return _order_intron_list(introns,order_by='length')
    elif order_by == 'acceptor_pssm':
        try:
            tmp = [ (-i.acceptor.pssm_score, i.length, i) for i in introns ]
            tmp.sort()
            return [ i for (a,b,i) in tmp ]
        except:
            # hmm... pssm_scores not available in the splice site objects
            return _order_intron_list(introns,order_by='length')
    else:
        # no ordering requested
        return introns

# end of function _order_intron_list

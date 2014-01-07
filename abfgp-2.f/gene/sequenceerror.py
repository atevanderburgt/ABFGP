"""
"""
# Module metadata
__authors__ = "Ate van der Burgt"
__license__ = "MIT"

# Gene Imports
from gene.gff import BasicGFF
from gene.validators import IsSequenceErrorCoordinate, IsOrf

# Other Imports
import dna2prot

# Python Imports
from os import popen as osPopen

# Import Global variables


class SequenceErrorCoordinate:
    """ class describing a sequence error coordinate (donor / acceptor like) """
    def __init__(self,pos):
        self.pos = pos
        self.pssm_score = 0.0   # compatibilty with intron/splicesites
        self.phase = 0          # compatibilty with intron/splicesites

# end of class SequenceErrorCoordinate


class ExonStart:
    """ class describing the (unexpected) start point of an exon object """
    def __init__(self,pos):
        self.pos   = pos
        self.start = pos
        self.pssm_score = 0.0   # compatibilty with intron/splicesites
        self.phase = 0          # compatibilty with intron/splicesites

# end of class ExonStart


class ExonEnd:
    """ class describing the (unexpected) end point of an exon object """
    def __init__(self,pos):
        self.pos   = pos
        self.start = pos
        self.pssm_score = 0.0   # compatibilty with intron/splicesites
        self.phase = 0          # compatibilty with intron/splicesites

# end of class ExonEnd


class SequenceError(BasicGFF):
    def __init__(self,seqerrorstart,seqerrorend,gff={}):
        # input validation
        if type(seqerrorstart) == type(int()):
            seqerrorstart = SequenceErrorCoordinate(seqerrorstart)
        else:
            IsSequenceErrorCoordinate(seqerrorstart)
        if type(seqerrorend) == type(int()):
            seqerrorend = SequenceErrorCoordinate(seqerrorend)
        else:
            IsSequenceErrorCoordinate(seqerrorend)

        # initialization
        BasicGFF.__init__(self)
        self._gff.update(gff)
        self.donor      = seqerrorstart
        self.acceptor   = seqerrorend
        self.length     = self.acceptor.pos-self.donor.pos
        self.start      = self.donor.pos
        self.end        = self.acceptor.pos
        self.typeof     = None
        # Estimated coordinate of this sequence error
        # It is possible to pinpoint this SequenceError to a certain
        # position but (re)label the estimated coordinate somewhere is.
        # This option is needed in ABFGP in case of Intertions/Deletions
        self._estimated_nt_coord = self.donor.pos

        # set shared nt/aa data. This mimmics attributes from
        # the Intron class and maintains readingframe in a genestructure
        self.shared_nts = "nnn"
        self.shared_aa  = "X"

    # end of function __init__


    def __str__(self):
        """ """
        return "<%s %s-%s (%s,%s,~%s)>" % (
            self.__class__.__name__, self.donor.pos,
            self.acceptor.pos, self.length, self.typeof,
            self._estimated_nt_coord )

    # end of function __str__


    def coords(self):
        """ """
        return ( self.donor.pos, self.acceptor.pos )

    # end of function coords


    def deepcopy(self):
        """ """
        seqerror = self.__class__(self.donor,self.acceptor,gff=self._gff)
        seqerror.typeof = self.typeof
        return seqerror

    # end of function deepcopy


    def togff(self,gff={}):
        """ Overwrites BasicGFF.togff() """
        if not gff.has_key('gname') and not self._gff.has_key('gname'):
            gff['gname'] = "%s-%s" % (self.start, self.end)
        self._gff['column9data']['typeof'] = self.typeof
        retgff = {}
        retgff.update(self._gff)
        retgff.update(gff)
        return BasicGFF.togff(self,gff=retgff)

    # end of function togff


    # Functions for compatibility with Intron classes
    def is_stopless_3n_intron(self): return False
    def is_3n_intron(self): return False

# end of class SequenceError


class SequenceErrorConnectingOrfs(SequenceError):
    def __init__(self,seqerrorstart,seqerrorend,orfDonor,orfAcceptor,gff={}):
        """ """
        # input validation
        IsOrf(orfDonor)
        IsOrf(orfAcceptor)

        # initialise basal SequenceError class
        SequenceError.__init__(self,seqerrorstart,seqerrorend,gff=gff)
        self.orfDonor   = orfDonor
        self.orfAcceptor= orfAcceptor

    # end of function __init__


    def __str__(self):
        """ """
        return "<%s %s-%s (%s,%s,~%s,[%s-%s])>" % (
            self.__class__.__name__, self.donor.pos,
            self.acceptor.pos, self.length, self.typeof,
            self._estimated_nt_coord,
            self.orfDonor.id, self.orfAcceptor.id )

    # end of function __str__


    def deepcopy(self):
        """ """
        seqerror = self.__class__(self.donor,self.acceptor,
                self.orfDonor, self.orfAcceptor,
                gff=self._gff)
        seqerror.typeof = self.typeof
        return seqerror

    # end of function deepcopy


    def dnasequence(self):
        """ """
        return self.orfDonor.inputgenomicsequence[self.start:self.end]
    # end of function dnasequence


    def barcode(self):
        """ A unique tuple with integers describing this object """
        return (    self.get_donor_nt_from_stop(), self.orfDonor.id,
                    self.start, None, self.end,
                    self.orfAcceptor.id, self.get_acceptor_nt_from_stop() )

    # end of function barcode


    # Functions for compatibility with Intron classes
    def get_donor_nt_from_stop(self):
        return self.orfDonor.endPY - self.donor.pos
    def get_acceptor_nt_from_stop(self):
        return self.acceptor.pos - self.orfAcceptor.startPY

# end of class SequenceErrorConnectingOrfs


def merge_orfs_with_sequenceerror(orfD,orfA,allow_double_stopcodon=False):
    """ """
    if orfA.startPY - orfD.endPY == 3:
        # single base miscall resulting in a stopcodon
        se = SequenceErrorConnectingOrfs(orfD.endPY,orfA.startPY,orfD,orfA)
        se.typeof = "basecall"
        se._gff['column9data']['StopCodon'] = se.dnasequence()
        se._gff['column9data']['typeof'] = "basecall"
        # return this sequence error
        return se
    elif orfA.startPY - orfD.endPY == 6 and allow_double_stopcodon:
        # Orfs separated by double stop codon(s)
        se = SequenceErrorConnectingOrfs(orfD.endPY,orfA.startPY,orfD,orfA)
        se.typeof = "basecall"
        se._gff['column9data']['StopCodon'] = [ se.dnasequence()[0:3], se.dnasequence()[3:] ]
        se._gff['column9data']['typeof'] = "basecall"
        # return this sequence error
        return se
    elif (orfD.endPY > orfA.startPY and orfD.startPY < orfA.startPY) or\
    (orfD.endPY < orfA.endPY and orfD.startPY > orfA.startPY) or\
    (orfD.endPY > orfA.endPY and orfD.startPY < orfA.startPY):
        # overlapping/including Orfs
        # possible insertion or deletion
        frames2error = {
            (0,1):  "insertion",    (0,2):  "deletion",
            (1,0):  "deletion",     (1,2):  "insertion",
            (2,1):  "deletion",     (2,0):  "insertion",
        }
        if frames2error[(orfD.frame,orfA.frame)] == "deletion":
            se = SequenceErrorConnectingOrfs(orfD.endPY,orfD.endPY+2,orfD,orfA)
            se.typeof = "deletion"
            se._gff['column9data']['typeof'] = "deletion"
        else:
            se = SequenceErrorConnectingOrfs(orfD.endPY-1,orfD.endPY+3,orfD,orfA)
            se.typeof = "insertion"
            se._gff['column9data']['typeof'] = "insertion"
        # return this sequence error
        return se
    else:
        # no obvious (singe!) sequence error
        return False

# end of function merge_orfs_with_sequenceerror

"""

"""
# Module metadata
__authors__ = "Ate van der Burgt"
__license__ = "MIT"

# Imports
from gene_gff import BasicGFF, GffWithSequenceFunctionality
from start import StartCodon, _handle_startcodon
from stop import StopCodon, _handle_stopcodon

# Import Global variables
from settings.genestructure import (
    MAX_EXON_NT_LENGTH,
    MIN_EXON_NT_LENGTH
    )
from settings.executables import (
    TCODE_MAX_NONCODING,
    TCODE_MIN_CODING,
    )


class Exon(BasicGFF):
    """
    """
    def __init__(self,acceptor,donor,gff={}):
        BasicGFF.__init__(self)
        self.donor      = donor
        self.acceptor   = acceptor

        # init by splice site objects
        self._init_by_splicesites()
        # update gff attributes
        self._gff.update(gff)

    # end of function __init__


    def _init_by_splicesites(self):
        """ """
        self.length     = self.donor.pos-self.acceptor.pos
        self.start      = self.acceptor.pos
        self.end        = self.donor.pos
        try:
            self.pssm_score = self.donor.pssm_score + self.acceptor.pssm_score
        except:
            self.pssm_score = 0.0 # no pssm_score in gff object(s)

    # end of function _init_by_splicesites


    def __str__(self):
        """ """
        return "<%s %s .. length=%s .. %s >" % (
            self.__class__.__name__,
            self.acceptor,
            self.length, 
            self.donor,
            )

    # end of function __str__


    def deepcopy(self):
        """ """
        return self.__class__(self.acceptor,self.donor,gff=self._gff)

    # end of function deepcopy


# end of class Exon


class FirstExon(Exon):
    """ """
    def __init__(self,startpos,donor,gff={}):
        startcodon = _handle_startcodon(startpos)
        # initialise basal Exon class
        Exon.__init__(self,startcodon,donor,gff=gff)
    # end of function __init__

# end of class FirstExon


class FinalExon(Exon):
    """ """
    def __init__(self,acceptor,endpos,gff={}):
        stopcodon  = _handle_stopcodon(endpos)
        # initialise basal Exon class
        Exon.__init__(self,acceptor,stopcodon,gff=gff)
    # end of function __init__

# end of class FinalExon


class SingleExon(Exon):
    """ """
    def __init__(self,startpos,endpos,gff={}):
        startcodon = _handle_startcodon(startpos)
        stopcodon  = _handle_stopcodon(endpos)
        # initialise basal Exon class
        Exon.__init__(self,startcodon,stopcodon,gff=gff)

# end of class SingleExon


class ExonOnOrf(Exon,GffWithSequenceFunctionality):
    """ """
    def __init__(self,acceptor,donor,orf,gff={}):
        # initialise basal Exon class
        Exon.__init__(self,acceptor,donor,gff=gff)
        self.orf = orf
        # this is an intermediate exon
        self.IS_FIRST = False
        self.IS_FINAL = False
    # end of function __init__

    def __str__(self):
        return "<%s: %s..%s (%snt) [orf(%s): %s..%s (%snt)] [pssm=%2.2f-%2.2f]>" % (
            self.__class__.__name__,
            self.acceptor.pos,
            self.donor.pos,
            self.length, 
            self.orf.id,
            self.orf.start,
            self.orf.end,
            self.orf.length,
            self.acceptor.pssm_score,
            self.donor.pssm_score,
            )

    # end of function __str__

    def deepcopy(self):
        """ """
        return self.__class__(self.acceptor,self.donor,self.orf,gff=self._gff)

    # end of function deepcopy


    def coords(self):
        """ """
        return (self.acceptor.pos,self.donor.pos)

    # end of function coords

    def dnasequence(self):
        """ """
        return self.orf.inputgenomicsequence[self.start:self.end]
    # end of function dnasequence


    def proteinsequence(self):
        """ """
        startAA = self.protein_start() - self.orf.protein_startPY
        endAA   = self.protein_end() - self.orf.protein_startPY
        return self.orf.protein_sequence[startAA:endAA]

    # end of function proteinsequence


    def protein_start(self):
        """ """
        startAA = self.orf.dnapos2aapos(self.start)
        if self.acceptor.phase != 0: startAA+=1
        return startAA

    # end of function protein_start


    def protein_end(self):
        """ """
        return self.orf.dnapos2aapos(self.end)

    # end of function protein_end


    def protein_range(self):
        """ """
        return range(self.protein_start(),self.protein_end()+1)

    # end of function protein_range


    def tcode(self):
        """ """
        return self.orf.tcode_score(self.protein_start(),self.protein_end())

    # end of function tcode


    def is_tcode_coding(self):
        """ """
        return self.tcode() >= TCODE_MIN_CODING

    # end of function is_tcode_coding


    def is_tcode_noncoding(self):
        """ """
        return self.tcode() < TCODE_MAX_NONCODING

    # end of function is_code_noncoding

# end of class ExonOnOrf


class FirstExonOnOrf(ExonOnOrf):
    """ """
    def __init__(self,startpos,donor,orf,gff={}):
        startcodon = _handle_startcodon(startpos)
        # initialise basal ExonOnOrf class
        ExonOnOrf.__init__(self,startcodon,donor,orf,gff=gff)
        self.IS_FIRST = True 

    # end of function __init__

    def protein_start(self):
        """ """
        return self.orf.dnapos2aapos(self.start)
    # end of function protein_start

# end of class FirstExonOnOrf


class SingleExonOnOrf(ExonOnOrf):
    """ """
    def __init__(self,startpos,endpos,orf,gff={}):
        startcodon = _handle_startcodon(startpos)
        stopcodon  = _handle_stopcodon(endpos)
        # initialise basal ExonOnOrf class
        ExonOnOrf.__init__(self,startcodon,stopcodon,orf,gff=gff)
        self.IS_FIRST = True
        self.IS_FINAL = True

    # end of function __init__

    def protein_start(self):
        """ """
        return self.orf.dnapos2aapos(self.start)
    # end of function protein_start

# end of function SingleExonOnOrf


class FinalExonOnOrf(ExonOnOrf):
    """ """
    def __init__(self,acceptor,endpos,orf,gff={}):
        stopcodon  = _handle_stopcodon(endpos)
        # initialise basal ExonOnOrf class
        ExonOnOrf.__init__(self,acceptor,stopcodon,orf,gff=gff)
        self.IS_FINAL = True 

    # end of function __init__

# end of class FinalExonOnOrf

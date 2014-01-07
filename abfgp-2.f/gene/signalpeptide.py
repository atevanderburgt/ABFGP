"""
"""

# Module metadata
__authors__ = "Ate van der Burgt"
__license__ = "MIT"

# Gene Imports
from gene_gff import BasicGFF

# Import Global variables
from settings.executables import EXECUTABLE_SIGNALP_VERSION

class SignalPeptide(BasicGFF):
    """ basal SignalPeptide class """
    def __init__(self,start,end,score,tss=None,intron=None,gff={}):
        """
        @attention: when spliced by intron:
                    a) provide Intron object
                    b) provide end coordinate as if UNSPLICED (length*3nt)
        """
        BasicGFF.__init__(self)
        self.start      = start
        self.pos        = self.start
        self.end        = end
        self.length     = self.end - self.start
        self.pssm_score = score
        self.tss        = tss
        self.intron     = intron
        self._gff.update(gff)
        if not self._gff.has_key('fmethod'):
            self._gff['fmethod'] = 'predSignalPeptide'

    # end of function __init__


    def togff(self,gff={}):
        """ Overwrites BasicGFF.togff() functionality """
        if not self._gff.has_key('gname'):
            self._gff['gname'] = "%s_%s" % (self.start,self.end)
        if not self._gff.has_key('gclass'):
            self._gff['gclass'] = self.__class__.__name__
        # call basal BasicGFF.togff() function
        gfflines = [ BasicGFF.togff(self) ]
        # gff default for TSS/ATG
        tss_gff = { 'fref'  : self._gff['fref'],
                    'gname' : self._gff['gname'],
                    'gclass': self._gff['gclass'], }

        # if SignalPeptide is spliced -> make 2 tracks!
        if self.intron:
            # calculate lengths
            signalp_length    = self.end-self.start
            signalp_5p_length = self.intron.donor.pos - self.start
            signalp_3p_length = signalp_length - signalp_5p_length
            # create SignalP 5p and 3p tracks
            exon5p = list(gfflines[0])
            exon3p = list(gfflines[0])
            exon5p[4] = self.intron.donor.pos
            exon3p[3] = self.intron.acceptor.pos + 1
            exon3p[4] = self.intron.acceptor.pos + signalp_3p_length
            if self.tss: exon5p[3]+=3
            gfflines = [ tuple(exon5p), tuple(exon3p) ]
            if self.tss: gfflines.insert(0, self.tss.togff(gff=tss_gff) )
        elif self.tss:
            # correct SignalPeptide start coord to 1nt after the ATG codon
            track = list(gfflines[0])
            track[3]+=3
            gfflines = [ tuple(track) ]
            gfflines.insert(0, self.tss.togff(gff=tss_gff) )
        else:
            pass
        # return list with gff tuples
        return gfflines

    # end of function togff


    def __str__(self):
        """ print one-liner with SignalPeptide data """
        scoreSIGNALP = self.pssm_score
        try:    scoreSIGNALP = "%1.2f" % scoreSIGNALP
        except: pass
        try:    scoreTSS = self.tss.pssm_score
        except: scoreTSS = None
        try:    scoreTSS = "%1.2f" % scoreTSS
        except: pass

        return "<%s %s-%s (%s %s)>" % (
                self.__class__.__name__,
                self.start, self.end, scoreTSS,scoreSIGNALP, 
                )

    # end of function __str__


# end of class SignalPeptide



class SignalPSignalPeptide(SignalPeptide):
    """ SignalP predicted Signalpeptide """
    def __init__(self,*args,**kwargs):
        """
        @attention: see basal SignalPeptide for *args and **kwargs
        """
        SignalPeptide.__init__(self,*args,**kwargs)
        self._gff['fsource'] = EXECUTABLE_SIGNALP_VERSION

    # end of function __init__

# end of class SignalPSignalPeptide

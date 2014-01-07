"""

"""

# Module metadata
__authors__ = "Ate van der Burgt"
__license__ = "MIT"

# Imports
from re import finditer
import dna2prot
from gene.pssm import pssmscore, parse_ic_file

# Import Global variables
from settings.genestructure import MAX_INTRON_NT_LENGTH, MIN_INTRON_NT_LENGTH, MIN_EXON_NT_LENGTH
from settings.splicesites import IC_ACCEPTOR_DATA_FILE, IC_ACCEPTOR_PATTERN_OFFSET
from settings.splicesites import IC_DONOR_DATA_FILE, IC_DONOR_PATTERN_OFFSET, IC_DONOR_NCGC_DATA_FILE


class WrongPatternLengthApplied(Exception):
    pass

# parse IC PSSM files of cannonical sites
IC_ACCEPTOR             = parse_ic_file(IC_ACCEPTOR_DATA_FILE)
IC_DONOR                = parse_ic_file(IC_DONOR_DATA_FILE)


class SpliceSiteBase:
    """
    coordinate system in python coordinates
                GFF     Python
        start   x       x-1
        end     y       y

    """
    def __init__(self,start,phase=None,strand='+',pattern=None,
        pattern_offset=(0,0),pssm_score=None,gff={}):
        """
        Initialization function of Basal SpliceSite logic
        Recommended is to use only one of the inheriting classes

        @type  start: number
    	@param start: start coord of site (e.g GT) or pattern (e.g. tgtGTcgat)

        @type  phase: number
    	@param phase: [0,1,2] or None

        @type  phase: number
    	@param phase: [0,1,2] or None, default None

        @type  strand: string
    	@param strand: ['+','-'] or None, default '+'

        @type  pattern: string
    	@param pattern: splice site pattern sequence (e.g. tgtGTcgat)

        @type  pattern_offset: tuple
    	@param pattern_offset: tuple with 2 integers, specifying 5p and 3p
                               offset towards actual splice site.
                               e.g. (3,4) for tgtGTcgat

        @type  pssm_score: float
    	@param pssm_score: PSSM score of splice site pattern

        @type  gff: dictionary
    	@param gff: dictionary with gff data. Possible keys are:
                    'fref', 'fmethod', 'fstrand', 'fsource', 'gclass', 'gname'

        """
        self.start      = start
        self.strand     = strand
        self.phase      = phase
        self.end        = None  # to be filled in inheriting classes
        self.pos        = None  # to be filled in inheriting classes
        self.canonical_donor = "GT"
        self.canonical_acceptor = "AG"
        # pattern, pattern_offset and pssm_score are only truely
        # functional in inheriting classes SpliceDonor and SpliceAcceptor
        self.pattern    = pattern
        self.pssm_score = pssm_score
        self._offset_5p = pattern_offset[0]
        self._offset_3p = pattern_offset[1]
        # and some GFF data thingies
        self._gff       = {
                'fref'      : None,
                'fsource'   : 'undefined',
                'fstrand'   : '+',
                'fscore'    : self.pssm_score,
                'column9data' : {},
            }
        self._gff.update(gff)
        # and error-check the data
        self.check()

    # end of function __init__

    def check(self):
        """
        Basic error check of object's validity
        """
        if self.phase not in [None,0,1,2]:
            raise WrongPhaseApplied
        if self.pattern:
            if len(self.pattern) != self._offset_5p + self._offset_3p + 2:
                raise WrongPatternLengthApplied
        else:
            if ( self._offset_5p, self._offset_3p ) != (0,0):
                raise WrongPatternLengthApplied

    # end of function check

    def __str__(self):
        """ print one-liner with SpliceSite data """
        checkA = self.pssm_score != None
        checkB = self.pattern != None
        _pattern = ""
        _score   = ""
        if checkA:
            _pattern = ", %s%s%s" % (
                self.pattern[0:self._offset_5p].lower(),
                self.pattern[self._offset_5p:-self._offset_3p].upper(),
                self.pattern[-self._offset_3p:].lower(),
                )       
        if checkB:
            _score = ", score=%2.3f" % self.pssm_score
        if self.pssm_score != None:
            return "<%s %s-%s (%s)%s%s>" % (
                self.__class__.__name__,
                self.start, self.end, self.phase, 
                _score, _pattern
                )
        else:
            return "<%s %s-%s (%s)>" % (
                self.__class__.__name__,
                self.start, self.end, self.phase, 
                )

    # end of function __str__


    def is_canonical(self):
        """ is splicesite canonical or not """
        return True

    # end of function is_canonical


    def togff(self,gff={}):
        """
        Return 8-element tuple of gff data.
        To be called from the subclasses, or
        to be overwritten in the subclasses!
        """
        # update self._gff dictionary
        self._gff.update(gff)
        # set defaults if not set already
        if not self._gff.has_key('gclass'):  self._gff['gclass']  = self.__class__.__name__
        if not self._gff.has_key('gname'):   self._gff['gname']   = str(self.pos)
        if not self._gff.has_key('fmethod'): self._gff['fmethod'] = self.__class__.__name__
        if not self._gff.has_key('fstart'):  self._gff['fstart']  = self.start+1
        if not self._gff.has_key('fstop'):   self._gff['fstop']   = self.end
        column9string = ""
        # take care for the formatting of column9data
        if self._gff['column9data']:
            tmp = []
            for key,value in self._gff['column9data'].iteritems():
                if str(value).find(" ") > -1: value = "'%s'" % value
                tmp.append( "; %s %s" % (key,value) )
            column9string = "".join(tmp)

        return (
            self._gff['fref'],
            self._gff['fsource'],
            self._gff['fmethod'],
            self._gff['fstart'],
            self._gff['fstop'],
            self.pssm_score,
            self._gff['fstrand'],
            self.phase,
            "%s %s%s" % ( self._gff['gclass'], self._gff['gname'], column9string )
        )

    # end of function togff

# end of class SpliceSite


class SpliceDonorGT(SpliceSiteBase):
    def __init__(self,start,phase=None,strand=None,gff={}):
        # initialize object from parental SliceSiteBase class
        SpliceSiteBase.__init__(self,start,phase=phase,strand=strand,pattern=None,
            pattern_offset=(0,0),pssm_score=None,gff=gff)
        self.end        = self.start+2
        self.pos        = self.start
        self.donor      = "GT"

    # end of function __init__

# end of class SpliceDonorGT


class SpliceAcceptorAG(SpliceSiteBase):
    def __init__(self,start,phase=None,strand=None,gff={}):
        # initialize object from parental SliceSiteBase class
        SpliceSiteBase.__init__(self,start,phase=phase,strand=strand,pattern=None,
            pattern_offset=(0,0),pssm_score=None,gff=gff)
        self.end        = self.start+2
        self.pos        = self.end
        self.acceptor   = "AG"

    # end of function __init__

# end of class SpliceAcceptorAG


class SpliceDonor(SpliceSiteBase):
    def __init__(self,start,pattern,donor="GT",phase=None,strand=None,
        pattern_offset=IC_DONOR_PATTERN_OFFSET,pssm_score=None,gff={}):
        # initialize object from parental SliceSiteBase class
        SpliceSiteBase.__init__(self,start,phase=phase,strand=strand,gff=gff,
            pattern=pattern,pattern_offset=pattern_offset,pssm_score=pssm_score)
        self.end        = self.start + 2 + self._offset_5p + self._offset_3p
        self.pos        = self.start + self._offset_5p
        self.donor      = donor.upper()

    # end of function __init__

    def is_canonical(self):
        """ is splicesite canonical or not """
        if self.canonical_donor == self.donor:
            return True
        else:
            return False

    # end of function is_canonical

    def togff(self,gff={}):
        """
        Overrides SpliceSiteBase.togff()
        Takes care for proper fstart and fstop setting
        """
        self._gff['fstart']  = self.pos+1
        self._gff['fstop']   = self.pos+2
        return SpliceSiteBase.togff(self,gff=gff)

    # end of function togff

# end of class SpliceDonor


class SpliceAcceptor(SpliceSiteBase):
    def __init__(self,start,pattern,acceptor="AG",phase=None,strand=None,
        pattern_offset=IC_ACCEPTOR_PATTERN_OFFSET,pssm_score=None,gff={}):
        # initialize object from parental SliceSiteBase class
        SpliceSiteBase.__init__(self,start,phase=phase,strand=strand,gff=gff,
            pattern=pattern,pattern_offset=pattern_offset,pssm_score=pssm_score)
        self.end        = self.start + 2 + self._offset_5p + self._offset_3p
        self.pos        = self.end - self._offset_3p
        self.acceptor   = acceptor

    # end of function __init__

    def is_canonical(self):
        """ is splicesite canonical or not """
        if self.canonical_acceptor == self.acceptor:
            return True
        else:
            return False

    # end of function is_canonical

    def togff(self,gff={}):
        """
        Overrides SpliceSiteBase.togff()
        Takes care for proper fstart and fstop setting
        """
        self._gff['fstart']  = self.pos-1
        self._gff['fstop']   = self.pos
        return SpliceSiteBase.togff(self,gff=gff)

    # end of function togff


# end of class SpliceAcceptor


def _projected_site_total_score(self):
    """
    Function in use as total_score() in 
    ProjectedSpliceDonor
    ProjectedSpliceAcceptor
    The total score of a projected site is a
    balanced average of PSSM, entropy, distance and
    number of cases this projection is supported by
    """
    case_correction = max([ 0.75*float(self.cases) , 1.0 ])
    pssm_part    = self.total_pssm_score / case_correction
    entropy_part = self.total_entropy / case_correction

    weight_pssm_part    = 1.0  # a single pssm score is in range -3.0..9.0
    weight_entropy_part = 1.0  # a single entropy score has range 0.0..1.0
                               # ATTENTION!! this lower range is corrected
                               # for in the ABFGP code at a later stage
                               # so, leave here as 1.0

    return ( ( pssm_part * weight_pssm_part ) +\
             ( entropy_part * weight_entropy_part ) ) /\
             float( abs(self.distance)/3 + 1 )

# end of function _projected_site_total_score


# end of function _projected_site_total_score


class ProjectedSpliceDonor(SpliceDonor):
    def __init__(self,start,pattern,splicedonor,distance=None,
        entropy=None,pssm_score=None,cases=1):
        # initialize object from parental SpliceDonor class
        SpliceDonor.__init__(self,start,phase=splicedonor.phase,
            strand=splicedonor.strand, gff=splicedonor._gff,
            pattern=pattern,pattern_offset=(0,0) )
        self.end        = self.start + 2 + self._offset_5p + self._offset_3p
        self.pos        = self.start + self._offset_5p
        self.donor      = pattern.upper()
        self.distance   = distance
        self.cases      = cases
        # needed for compatibilty with SpliceSiteBase object
        self.pssm_score = pssm_score
        self.entropy    = entropy
        # special thingies for Projected flavour
        self.total_entropy      = entropy
        self.total_pssm_score   = pssm_score
        self.average_entropy    = entropy / float(cases)
        self.average_pssm_score = pssm_score / float(cases)

    # end of function __init__
    
    def is_canonical(self):
        """ is splicesite canonical or not """
        return False

    # end of function is_canonical

    def total_score(self):
        """ total score of this pojected splice site"""
        return _projected_site_total_score(self)

    # end of function total_score

    def __str__(self):
        """
        Overrides SpliceDonor.__str__()
        """
        retstr = SpliceDonor.__str__(self)
        retstr = retstr.replace(", score=", " [%sx], b.e.t.w. %1.2f, score=" % (self.cases,self.total_entropy) )
        return retstr

    # end of function __str__

    def togff(self,gff={}):
        """
        Overrides SpliceSiteBase.togff()
        Takes care for proper fstart and fstop setting
        """
        self._gff['column9data'] = {
                'Entropy'       : self.total_entropy,
                'Alignedsites'  : self.cases,
                'Distance'      : self.distance,
            }
        return SpliceDonor.togff(self,gff=gff)

    # end of function togff

    
# end of class ProjectedSpliceDonor


class ProjectedSpliceAcceptor(SpliceAcceptor):
    def __init__(self,start,pattern,spliceacceptor,distance=None,
        entropy=None,pssm_score=None,cases=1):
        # initialize object from parental SpliceAcceptor class
        SpliceAcceptor.__init__(self,start,phase=spliceacceptor.phase,
            strand=spliceacceptor.strand, gff=spliceacceptor._gff,
            pattern=pattern,pattern_offset=(0,0) )
        self.end        = self.start + 2 + self._offset_5p + self._offset_3p
        self.pos        = self.start + self._offset_5p
        self.acceptor   = pattern.upper()
        self.distance   = distance
        self.cases      = cases
        # needed for compatibilty with SpliceSiteBase object
        self.pssm_score = pssm_score
        self.entropy    = entropy
        # special thingies for Projected flavour
        self.total_entropy      = entropy
        self.total_pssm_score   = pssm_score
        self.average_entropy    = entropy / float(cases)
        self.average_pssm_score = pssm_score / float(cases)

    # end of function __init__
    
    def is_canonical(self):
        """ is splicesite canonical or not """
        return False

    # end of function is_canonical

    def total_score(self):
        """ total score of this pojected splice site"""
        return _projected_site_total_score(self)

    # end of function total_score

    def __str__(self):
        """
        Overrides SpliceAcceptor.__str__()
        """
        retstr = SpliceAcceptor.__str__(self)
        retstr = retstr.replace(", score=", " [%sx], b.e.t.w. %1.2f, score=" % (self.cases,self.total_entropy) )
        return retstr

    # end of function __str__

    def togff(self,gff={}):
        """
        Overrides SpliceSiteBase.gff()
        Takes care for proper fstart and fstop setting
        """
        self._gff['column9data'] = {
                'Entropy'       : self.total_entropy,
                'Alignedsites'  : self.cases,
                'Distance'      : self.distance,
            }
        return SpliceAcceptor.togff(self,gff=gff)

    # end of function togff

# end of class ProjectedSpliceAcceptor


class Intron:
    def __init__(self,donor,acceptor,shared_nts,gff={}):
        self.donor      = donor
        self.acceptor   = acceptor
        self.length     = acceptor.pos-donor.pos
        self.start      = donor.pos
        self.end        = acceptor.pos
        self.phase      = None
        self.nt_from_stop_donor = None
        self.nt_from_stop_acceptor = None

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
            pass
        # set/update gff data
        self._gff       = {
                'fref'      : None,
                'fscource'  : 'undefined',
                'fstrand'   : '+',
                'fscore'    : ".",
                'column9data' : {},
            }
        self._gff.update(gff)
        # and set some shared nt/aa data
        self.shared_nts = shared_nts 
        self.shared_aa  = ""
        if self.shared_nts:
            self.shared_aa = dna2prot.dna2protein(self.shared_nts)

    # end of function __init__

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


    def togff(self,gff={}):
        """
        Return 8-element tuple of gff data.
        To be called from the subclasses, or
        to be overwritten in the subclasses!
        """
        # update self._gff dictionary
        self._gff.update(gff)
        # set defaults if not set already
        if not self._gff.has_key('gclass'):  self._gff['gclass']  = self.__class__.__name__
        if not self._gff.has_key('gname'):   self._gff['gname']   = "%s-%s" % (self.start, self.end)
        if not self._gff.has_key('fmethod'): self._gff['fmethod'] = self.__class__.__name__
        column9string = ""
        # take care for the formatting of column9data
        if self._gff['column9data']:
            tmp = []
            for key,value in self._gff['column9data'].iteritems():
                if str(value).find(" ") > -1: value = "'%s'" % value
                tmp.append( "; %s %s" % (key,value) )
            column9string = "".join(tmp)

        return (
            self._gff['fref'],
            self._gff['fsource'],
            self._gff['fmethod'],
            self.start+1,
            self.end,
            self._gff['fscore'],
            self._gff['fstrand'],
            self.phase,
            "%s %s%s" % ( _gff['gclass'], self._gff['gname'], column9string )
        )

    # end of function togff

# end of class Intron


class IntronConnectingOrfs(Intron):
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
        # initialise basal Intron class
        Intron.__init__(self,donor,acceptor,shared_nts,gff=gff)
        self.orfDonor              = orfDonor
        self.orfAcceptor           = orfAcceptor
        self.nt_from_stop_donor    = self.orfDonor.endPY - self.donor.pos
        self.nt_from_stop_acceptor = self.acceptor.pos - self.orfAcceptor.startPY
        self.donor_aa_pos          = (donor.pos-donor.phase) / 3
        self.acceptor_aa_pos       = (acceptor.pos-acceptor.phase) / 3
        self.binary_entropy_donor   = binary_entropy_donor
        self.binary_entropy_acceptor= binary_entropy_acceptor

    # end of function __init__


    def barcode(self):
        """
        A unique tuple with integers describing this object
        """
        return ( self.nt_from_stop_donor, self.orfDonor.id, self.start, self.phase, self.end, self.orfAcceptor.id, self.nt_from_stop_acceptor )

    # end of function barcode


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


class BasicGFF:

    def __str__(self):
        """ """
        return "<%s %s>" % ( self.__class__.__name__,self.pos)

    # end of function __str__


    def _get_gff_defaults(self):
        """ """
        return {
                'fref'      : None,
                'fsource'   : 'undefined',
                'fstrand'   : '+',
                'fscore'    : '.',
                'column9data' : {},
            }

    # end of function _get_gff_defaults

    def togff(self,gff={}):
        """
        Return 8-element tuple of gff data.
        To be called from the subclasses, or
        to be overwritten in the subclasses!
        """
        # update self._gff dictionary
        self._gff.update(gff)
        # set defaults if not set already
        if not self._gff.has_key('gclass'):  self._gff['gclass']  = self.__class__.__name__
        if not self._gff.has_key('gname'):   self._gff['gname']   = str(self.pos)
        if not self._gff.has_key('fmethod'): self._gff['fmethod'] = self.__class__.__name__
        if not self._gff.has_key('fstart'):  self._gff['fstart']  = self.start+1
        if not self._gff.has_key('fstop'):   self._gff['fstop']   = self.end
        column9string = ""
        # take care for the formatting of column9data
        if self._gff['column9data']:
            tmp = []
            for key,value in self._gff['column9data'].iteritems():
                if str(value).find(" ") > -1: value = "'%s'" % value
                tmp.append( "; %s %s" % (key,value) )
            column9string = "".join(tmp)

        return (
            self._gff['fref'],
            self._gff['fsource'],
            self._gff['fmethod'],
            self._gff['fstart'],
            self._gff['fstop'],
            self.pssm_score,
            self._gff['fstrand'],
            self.phase,
            "%s %s%s" % ( self._gff['gclass'], self._gff['gname'], column9string )
        )

    # end of function togff    


# end of class BasicGFF


class StartCodon(BasicGFF):
    def __init__(self,pos,gff={}):
        """ """
        self.pos        = pos
        self.start      = self.pos
        self.end        = self.start+3
        self.pssm_score = 1.0  # default, dummy value
        self._gff       = self._get_gff_defaults()
        self.phase      = 0
    # end of function __init__

# end of class StartCodon


class StopCodon(BasicGFF):
    def __init__(self,pos,gff={}):
        """ """
        self.pos        = pos
        self.start      = self.pos
        self.end        = self.start+3
        self.pssm_score = 1.0  # default, dummy value
        self._gff       = self._get_gff_defaults()
        self.phase      = 0
    # end of function __init__

# end of class StopCodon


class CodingBlockStart(BasicGFF):
    def __init__(self,pos,gff={}):
        """ """
        self.pos        = pos
        self.start      = self.pos
        self.end        = self.start+2
        self.pssm_score = 1.0  # default, dummy value
        self.phase      = None # no phase for this feature!
        self._gff       = self._get_gff_defaults()

    # end of function __init__

    def __str__(self):
        """ """
        return "<%s pos=%s>" % (self.__class__.__name__,self.pos)
    # end of function __str__

# end of class CodingBlockStart


class CodingBlockEnd(BasicGFF):
    def __init__(self,pos,gff={}):
        """ """
        self.pos        = pos
        self.start      = self.pos-2
        self.end        = self.pos
        self.pssm_score = 1.0  # default, dummy value
        self.phase      = None # no phase for this feature!
        self._gff       = self._get_gff_defaults()

    # end of function __init__

    def __str__(self):
        """ """
        return "<%s pos=%s>" % (self.__class__.__name__,self.pos)
    # end of function __str__

# end of class CodingBlockEnd


class Exon:
    def __init__(self,acceptor,donor,gff={}):
        self.donor      = donor
        self.acceptor   = acceptor
        self.length     = donor.pos-acceptor.pos
        self.start      = acceptor.pos
        self.end        = donor.pos
        self._gff       = gff
        self.pssm_score = self.donor.pssm_score + self.acceptor.pssm_score

    # end of function __init__

    def __str__(self):
        return "<%s %s .. length=%s .. %s >" % (
            self.__class__.__name__,
            self.acceptor,
            self.length, 
            self.donor,
            )

    def gff(self):
        """ return 8-element tuple of gff data """
        pass

# end of class Exon


class FirstExon(Exon):
    def __init__(self,startpos,donor,gff={}):
        # initialise basal Exon class
        startcodon = StartCodon(startpos)
        Exon.__init__(self,startcodon,donor,gff=gff)
    # end of function __init__

# end of class FirstExon

class FinalExon(Exon):
    def __init__(self,acceptor,endpos,gff={}):
        # initialise basal Exon class
        stopcodon = StopCodon(endpos)
        Exon.__init__(self,acceptor,stopcodon,gff=gff)
    # end of function __init__

# end of class FinalExon

class ExonOnOrf(Exon):
    def __init__(self,acceptor,donor,orf,gff={}):
        # initialise basal Exon class
        Exon.__init__(self,acceptor,donor,gff=gff)
        self.orf = orf
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

# end of class ExonOnOrf

class FirstExonOnOrf(ExonOnOrf):
    def __init__(self,startpos_or_tsspssm_site,donor,orf,gff={}):
        # do some integrity check
        if startpos_or_tsspssm_site.__class__.__name__ in ['TranslationalStartSiteATG','TranslationalStartSite']:
            startcodon = startpos_or_tsspssm_site
        elif type(startpos_or_tsspssm_site) == type(int()):
            startcodon = StartCodon(startpos_or_tsspssm_site)
        else:
            raise InproperlyAppliedArgument, "``startpos_or_tsspssm_site`` not an expected object or integer"
        # initialise basal Exon class
        ExonOnOrf.__init__(self,startcodon,donor,orf,gff=gff)
    # end of function __init__

    def protein_start(self):
        """ """
        return self.orf.dnapos2aapos(self.start)
    # end of function protein_start

# end of class FirstExonOnOrf


class SingleExonOnOrf(ExonOnOrf):
    def __init__(self,startpos_or_tsspssm_site,endpos,orf,gff={}):
        # do some integrity check
        if startpos_or_tsspssm_site.__class__.__name__ in ['TranslationalStartSiteATG','TranslationalStartSite']:
            startcodon = startpos_or_tsspssm_site
        elif type(startpos_or_tsspssm_site) == type(int()):
            startcodon = StartCodon(startpos_or_tsspssm_site)
        else:
            raise InproperlyAppliedArgument, "``startpos_or_tsspssm_site`` not an expected object or integer"
        # initialise basal Exon class
        stopcodon = StopCodon(endpos)
        ExonOnOrf.__init__(self,startcodon,stopcodon,orf,gff=gff)
    # end of function __init__

    def protein_start(self):
        """ """
        return self.orf.dnapos2aapos(self.start)
    # end of function protein_start

# end of function SingleExonOnOrf

class FinalExonOnOrf(ExonOnOrf):
    def __init__(self,acceptor,endpos,orf,gff={}):
        # initialise basal Exon class
        stopcodon = StopCodon(endpos)
        ExonOnOrf.__init__(self,acceptor,stopcodon,orf,gff=gff)
    # end of function __init__

# end of class FinalExonOnOrf


def find_GT_donor_sites(seq):
    """
    Find canonical GT donor sites on input sequence

    @type seq:  string
    @param seq: DNA sequence

    @rtype:  list
    @return: list with GT =start= coordinates (0-based Python string coords)
    """
    sites = []
    for item in finditer( "GT", seq.upper().replace("U","T") ):
        sites.append( item.start() )
    sites.reverse()
    return sites

# end of function find_GT_donor_sites


def find_AG_acceptor_sites(seq):
    """
    Find canonical AG acceptor sites on input sequence

    @type seq:  string
    @param seq: DNA sequence

    @rtype:  list
    @return: list with AG =end= coordinates (0-based Python string coords)
    """
    sites = []
    for item in finditer( "AG", seq.upper() ):
        sites.append( item.start()+2 )
    return sites

# end of function find_AG_acceptor_sites


def _score_splice_site(seq,splicetype="donor"):
    """
    Return PSSM splice site score

    @type seq:  string
    @param seq: DNA sequence of EXACT length of the PSSM

    @type splicetype:   string
    @param splicetype:  'donor' or 'acceptor'

    @rtype:     float
    @return:    PSSM score for this splice site
    """
    if splicetype=='acceptor':
        return pssmscore(seq,IC_ACCEPTOR)
    else:
        return pssmscore(seq,IC_DONOR)

# end of _score_splice_site


def scan_pssm_splice_site(seq,splicetype="donor",override_pattern_offset=(),min_pssm_score=None,
    allow_non_canonical=False,non_canonical_min_pssm_score=0.0):
    """
    Find splice sites by a PSSM on input sequence

    @type seq:  string
    @param seq: DNA sequence of EXACT length of the PSSM

    @type splicetype:   string
    @param splicetype:  'donor' or 'acceptor'

    @type min_pssm_score:   float
    @param min_pssm_score:

    @type allow_non_canonical:  boolean
    @param allow_non_canonical: True of False

    @type non_canonical_min_pssm_score:   float
    @param non_canonical_min_pssm_score:

    @type override_pattern_offset:  tuple
    @param override_pattern_offset: tuple with 2 integers; use cautiously!!

    @rtype:     list
    @return:    list with SpliceDonors or SpliceAcceptors
    """
    if splicetype=='acceptor':
        PSSM_MATRIX     = IC_ACCEPTOR
        pattern_offset  = IC_ACCEPTOR_PATTERN_OFFSET 
        canonical       = "AG"
    else:
        PSSM_MATRIX     = IC_DONOR
        pattern_offset  = IC_DONOR_PATTERN_OFFSET 
        canonical       = "GT"

    if allow_non_canonical:
        # obtain PSSM_IC for non-canonical (GC) donors
        IC_NCGC_DONOR = parse_ic_file(IC_DONOR_NCGC_DATA_FILE)

    # hmm... somebody knows what he or she is doing ;-)
    if override_pattern_offset:
        pattern_offset = override_pattern_offset

    pssmlength = len(PSSM_MATRIX)
    sites = []
    for offset in range(0, len(seq) - pssmlength + 1 ):
        # get sequence slice of pattern and actual splice site
        seqpart = seq[offset:offset+pssmlength].upper()
        splicesite = seqpart[pattern_offset[0]:-pattern_offset[1]]

        # continue if non-canonical sites if not requested for
        if not allow_non_canonical and splicesite != canonical:
            continue
        elif splicesite == canonical:
            # score this splicesite
            score = _score_splice_site(seqpart,splicetype=splicetype)
            # check if site must be stored
            if min_pssm_score or min_pssm_score == 0.0:
                if score < min_pssm_score:
                    continue
        elif splicesite != canonical and splicetype == 'donor':
            # score non-canonical donor site
            score = pssmscore(seqpart,IC_NCGC_DONOR) 
            # check if site must be stored
            if non_canonical_min_pssm_score or non_canonical_min_pssm_score == 0.0:
                if score < non_canonical_min_pssm_score:
                    continue
        else:
            continue

        if splicetype=='acceptor':
            a = SpliceAcceptor(offset,seqpart,acceptor=splicesite,pssm_score=score)
            sites.append(a)
        else:
            d = SpliceDonor(offset,seqpart,donor=splicesite,pssm_score=score)
            sites.append(d)

    # return sites for Donor
    if splicetype == 'donor':
        sites.reverse()

    # and return
    return sites

# end of function scan_pssm_splice_site


def scan_orf_for_pssm_splice_sites(orf,splicetype="donor",min_pssm_score=None,
    allow_non_canonical=False,non_canonical_min_pssm_score=0.0):
    """
    """
    if splicetype=='acceptor':
        pattern_offset  = IC_ACCEPTOR_PATTERN_OFFSET 
        offset_5p       = pattern_offset[0]+4 
        offset_3p       = 0
    else:
        pattern_offset  = IC_DONOR_PATTERN_OFFSET 
        offset_5p       = 0
        offset_3p       = pattern_offset[1]+4 

    # get (elongated) sequence of this orf
    seqslice = orf.nucleotidesequence(extra_start=offset_5p,extra_end=offset_3p)

    # scan for occurrences
    sites = scan_pssm_splice_site(seqslice,splicetype=splicetype,
        min_pssm_score=min_pssm_score,allow_non_canonical=allow_non_canonical,
        non_canonical_min_pssm_score=non_canonical_min_pssm_score)

    # correct site positions to absolute coords
    # and set correct phase of the splice site
    for site in sites:
        site.start = site.start + orf.startPY - offset_5p
        site.end   = site.end   + orf.startPY - offset_5p
        site.pos   = site.pos   + orf.startPY - offset_5p
        site.phase = (site.pos - orf.startPY) % 3

    # and return the splice sites on this orf
    return sites

# end of function scan_orf_for_pssm_splice_sites


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


def get_shared_nucleotides_at_splicesite(orfA,orfD,acceptor,donor):
    """
    """
    # do splice site phase compatibility check
    if acceptor.phase != donor.phase: raise IncompatibleSpliceSitePhases

    if donor.phase == 0:
        shared_nts = ""
    elif donor.phase == 1:
        shared_nts = "%s%s" % (
            orfD.inputgenomicsequence[donor.pos-1:donor.pos],
            orfA.inputgenomicsequence[acceptor.pos:acceptor.pos+2]
            )
    elif donor.phase == 2:
        shared_nts = "%s%s" % (
            orfD.inputgenomicsequence[donor.pos-2:donor.pos],
            orfA.inputgenomicsequence[acceptor.pos:acceptor.pos+1]
            )
    else:
        # what else !? Unexpected error!
        raise UnexpectedSpliceSitePhase

    # return shared_nucleotides
    return shared_nts

# end of function get_shared_nucleotides_at_splicesite

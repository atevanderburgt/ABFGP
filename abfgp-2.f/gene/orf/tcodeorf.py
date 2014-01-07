"""
Class describing an OpenReadingFrame (stop->stop) on a DNA sequence with
EMBOSS TCODE data attached to it.
"""

# Module metadata
__authors__ = "Ate van der Burgt"
__license__ = "MIT"

# Orf Imports
from basicorf import BasicOrf

# Python Imports

# Global variable Imports
from settings.executables import (
    EXECUTABLE_TCODE_WINDOW,
    EXECUTABLE_TCODE_STEP,
    TCODE_MAX_NONCODING,
    TCODE_MIN_CODING,
    )

class TcodeFunctions:
    """ """
    def tcode_status(self):
        """
        """
        if not self._tcode_stats:
            status = []
            for statid in (1,0,-1):
                status.append( self._tcode_status_list.count(statid) )
            statsum = float(sum(status))
            if statsum == 0.0:
                self._tcode_stats = ( None, None, None, 0 )
            else:
                rettup = [ float(stat)/statsum for stat in status ]
                rettup.append(int(statsum))
                self._tcode_stats = tuple( rettup )
        # and return the stats
        return self._tcode_stats

    # end of function tcode_status

    def tcode_coding_ratio(self):
        return self.tcode_status()[0]
    def tcode_noopinion_ratio(self):
        return self.tcode_status()[1]
    def tcode_noncoding_ratio(self):
        return self.tcode_status()[2]

    def tcode_is_coding(self,CODING_NONCODING_RATIO=10):
        """ """
        if self.tcode_status()[0] == None:
            return None
        elif self.tcode_status()[0] and self.tcode_status()[0] >\
            self.tcode_status()[2]*CODING_NONCODING_RATIO:
            return True
        else:
            return False

    # end of function tcode_is_coding


    def tcode_is_noncoding(self):
        """ """
        if self.tcode_status()[2] == None:
            return None
        elif self.tcode_status()[2] and not self.tcode_status()[0] and\
            self.tcode_status()[2] > self.tcode_status()[1]*2:
            return True
        else:
            return False

    # end of function tcode_is_noncoding


    def tcode_symbolic(self):
        """
        C for Coding
        N for Non-coding
        ? for No opinion
        ! for absent status / not set
        $ for unexpected -> error!
        """
        if self.tcode_is_coding():
            return "C"
        elif self.tcode_is_noncoding():
            return "N"
        elif self.tcode_is_coding() == None:
            return "!"
        elif self.tcode_is_coding() == False and self.tcode_is_noncoding() == False:
            return "?"
        else:
            # what else!? this should be non-existing
            return "$"

    # end of function tcode_symbolic


    def average_tcode_score_of_range(self,raw_tcode_data,startCoord,endCoord):
        """
        """
        scores = []
        for (start,stop,score,status) in raw_tcode_data:
            coord = ((stop-start-1)/2) + start-1
            if coord < startCoord:
                continue   # in front of orf
            elif coord >= endCoord:
                break     # after orf; break!
            scores.append(score)
        # Remove None points (TCODE ~= window-based -> no score at tails of the sequence!
        while None in scores: scores.remove(None)
        # return the average of the summed scores
        if not scores: scores = [0.0]
        return sum(scores) / float(len(scores))

    # end of function average_tcode_score_of_range


    def find_lowest_scoring_tcode_stretch(self,raw_tcode_data,startCoord,endCoord,window_size=EXECUTABLE_TCODE_WINDOW):
        """
        Find the area with the lowest EMBOSS-TCODE scores, probably a non-coding area

        @type  startCoord: Integer
        @param startCoord: Absolute nt coordinate to take as offset start

        @type  endCoord: Integer
        @param endCoord: Absolute nt coordinate to take as offset end
        """

        selectedrange = []
        for (start,stop,score,status) in raw_tcode_data:
            center = start + ((stop-start)/2)
            if center  < startCoord:
                continue   # in front of range
            elif center > endCoord:
                break     # behind range; break!
            selectedrange.append(score)

        tcode_score_batch = window_size / EXECUTABLE_TCODE_STEP
        tcodescores = []

        for offset_sta in range(0,len(selectedrange)-tcode_score_batch):
            offset_end = offset_sta + tcode_score_batch
            if None in selectedrange[offset_sta:offset_end]: continue
            average = sum(selectedrange[offset_sta:offset_end]) / tcode_score_batch
            tcodescores.append(average)

        if not tcodescores:
            return None
        else:
            return min(tcodescores)

    # end of function find_lowest_scoring_tcode_stretch


# end of class TcodeFunctions


class TcodeOrf(BasicOrf,TcodeFunctions):
    """
    """
    def __init__(self,*args,**kwargs):
        """
        """
        BasicOrf.__init__(self,*args,**kwargs)
        # list for tcode scores on certain positions
        self._tcode_score_list  = [None]*self.length
        # list for tcode scores on certain positions
        # status options are
        # None  no score for this position)
        # -1    Non-coding
        # 0     No opinion
        # 1     Coding
        self._tcode_status_list = [None]*self.length
        self._tcode2status = { 'Non-coding': -1, 'No opinion': 0, 'Coding': 1 }
        self._tcode_stats = ()
        self._RAW_TCODE_DATA = []

    # end of function __init__


    def _get_gff_fscore(self):
        """
        """
        tcs = self.tcode_symbolic()
        if tcs == 'C':   fscore = 1
        elif tcs == 'N': fscore = -1
        elif tcs == '?': fscore = 0
        else:            fscore = 0 # this cannot happen!
        return fscore
    # end of function _get_gff_fscore

    def __str__(self):
        """ """
        return """<Orf:%s %s %s %s-%s>""" % (
            self.id,
            self.tcode_symbolic(),
            self.header,
            self.startPY,
            self.endPY,
            )

    # end of function __str__


    ####################################################################
    #### Functions for EMBOSS-TCODE data accessibility in Orfs      ####
    ####################################################################

    def add_tcode_data(self,tcodeoutput):
        """
        """
        # store complete sequence TCODE data to _RAW_TCODE_DATA
        self._RAW_TCODE_DATA = tcodeoutput
        for (start,stop,score,status) in tcodeoutput:
            coord = ((stop-start-1)/2) + start-1
            if coord < self.startPY:
                continue   # in front of orf
            elif coord >= self.endPY:
                break     # after orf; break!
            else:
                # yep; append tcode data to orf
                orf_rel_coord = coord-self.startPY
                self._tcode_score_list[orf_rel_coord]=score
                self._tcode_status_list[orf_rel_coord]=self._tcode2status[status]

    # end of function add_tcode_data


    def tcode_entropy_around_pos(self,pos,aa_step_size=33,
        range_of_steps=range(-5,3),window_left=None,window_right=None,
        default_window_size=201):
        """
        Get ratios of averaged EMBOSS-TCODE score of a window left and right at positions around a requested position

        @type  pos: Integer (positive)
        @param pos: Absolute AA coordinate to take as center

        @type  aa_step_size: Integer (positive)
        @param aa_step_size: Size of step (in AA) to next position to take into account

        @type  range_of_steps: list 
        @param range_of_steps: range of steps to calculate; default range(-5,3); 0 == window around `pos` 

        @type  window_left: Integer (positive)
        @param window_left: nt length to take as a left/5p window

        @type  window_rigth: Integer (positive)
        @param window_rigth: nt length to take as a rigth/3p window

        @type  default_window_size: Integer (positive)
        @param default_window_size: nt length of default window size when left/rigth window is omitted

        @rtype:  tuple of (positive) floats
        @return: tuple of ratios of averaged tcode window scores around a requested position
        """
        scores = []
        for position in [ pos+(step*aa_step_size) for step in range_of_steps ]:
            (scoreL,scoreR) = self.tcode_entropy_of_pos(position,
                    window_left=window_left,
                    window_right=window_right,
                    default_window_size=default_window_size
                    )
            scores.append( (scoreL,scoreR) )
        return scores

    # end of function tcode_entropy_around_pos


    def tcode_entropy_of_pos(self,pos,window_left=None,window_right=None,default_window_size=201):
        """
        Get average EMBOSS-TCODE score of a window left and right of a requested position

        @type  pos: Integer (positive)
        @param pos: Absolute AA coordinate to take as center

        @type  window_left: Integer (positive)
        @param window_left: nt length to take as a left/5p window

        @type  window_rigth: Integer (positive)
        @param window_rigth: nt length to take as a rigth/3p window

        @type  default_window_size: Integer (positive)
        @param default_window_size: nt length of default window size when left/rigth window is omitted

        @rtype:  tuple of two (positive) floats
        @return: tuple of average tcode score of a window left (5p) and right (3p) of the requested position
        """

        # Notice when pos is (way) out of leage of this orf
        # If so, tails are never the boundaries of the orf
        # but just the windows
        POS_IS_OUTSIDE_ORF = False
        if pos < self.protein_startPY: POS_IS_OUTSIDE_ORF = True
        if pos >= self.protein_endPY:  POS_IS_OUTSIDE_ORF = True

        if POS_IS_OUTSIDE_ORF:
            if not window_right: window_right = default_window_size
            if not window_left:  window_left = default_window_size


        # calculate RIGHT window
        if not POS_IS_OUTSIDE_ORF and not window_right:
            tcodeR = self.tcode_score(start=pos)
        elif not POS_IS_OUTSIDE_ORF and window_right and (window_right/3)+pos < self.protein_endPY:
            tcodeR = self.tcode_score(start=pos,end=(window_right/3)+pos)
        else:
            startDNA = pos * 3
            endDNA   = ( (window_right/3)+pos ) * 3
            tcodeR = self.average_tcode_score_of_range(
                    self._RAW_TCODE_DATA,
                    startDNA,endDNA
                    )

        # calculate LEFT window
        if not POS_IS_OUTSIDE_ORF and not window_left:
            tcodeL = self.tcode_score(end=pos)
        elif not POS_IS_OUTSIDE_ORF and window_left and pos-(window_left/3) >= self.protein_startPY:
            tcodeL = self.tcode_score(start=pos-(window_left/3),end=pos)
        else:
            endDNA   = pos * 3
            startDNA = ( pos - (window_left/3) ) * 3
            tcodeL = self.average_tcode_score_of_range(
                    self._RAW_TCODE_DATA,
                    startDNA,endDNA
                    )
        # return average tcode scores on both sides of this position
        return ( tcodeL, tcodeR )

    # end of function tcode_entropy_of_pos


    def tcode_score(self,start=None,end=None):
        """
        Get average tcode score of (part of) this Orf

        @type  start: Integer (or None)
        @param start: Absolute AA coordinate to take as offset start

        @type  end: Integer (or None)
        @param end: Absolute AA coordinate to take as offset end

        @rtype:  float (in range ~0.49 .. ~1.26 )
        @return: average tcode score for this sequence part
        """
        startDNA = self.startPY
        endDNA   = self.endPY
        if start:   startDNA = self.aapos2dnapos(start)
        if end:     endDNA   = self.aapos2dnapos(end)
        # check coordinates; are they not boud to get out-of-range?
        # in some weird (or bogus?) cases (hmmsearches?), PacbPORF can
        # be a few nt longer than the Orf it self (e.g. case mgg0014)
        # temporarily, correct the coordinates here.
        if endDNA <= (self.startPY + len(self._tcode_score_list)):
            pass
        elif endDNA == (self.startPY + len(self._tcode_score_list) +1):
            # exact EOF Orf position / donor IN stopcodon!
            endDNA = self.startPY + len(self._tcode_score_list)
        else:
            ###if endDNA > (self.startPY + len(self._tcode_score_list)):
            print "Warning: end coordinate changed in Orf.tcode_score (%s->%s) %s" % (
                    endDNA, self.startPY + len(self._tcode_score_list), self)
            endDNA = self.startPY + len(self._tcode_score_list) 
        # get tcode position score list
        points = [ self._tcode_score_list[coord-self.startPY] for coord in range(startDNA,endDNA) ]
        # Remove None points (TCODE ~= window-based -> no score at tails of the sequence!
        while None in points: points.remove(None)
        # Error check for empty lists
        if len(points) == 0: points = [0.0]
        # return the average score
        return sum(points) / float(len(points))

    # end of function tcode_score


    def aa_pos_tcode_score(self,pos):
        """
        Average tcode score of a single amino acid of this Orf

        @type  pos: Integer (or None)
        @param pos: Absolute AA coordinate
 
        @rtype:  float (in range ~0.49 .. ~1.26 )
        @return: average tcode score for this sequence part
        """
        return  self.tcode_score(start=pos,end=pos+1)

    # end of function aa_pos_tcode_score


    def has_lowscroring_tcode_stretch(self,start=None,end=None,aa_window_size=15,min_tcode_score=0.740):
        """
        Is there an area with low EMBOSS-TCODE scores, probably non-coding, in this orf?

        @type  start: Integer (or None)
        @param start: Absolute AA coordinate to take as offset start

        @type  end: Integer (or None)
        @param end: Absolute AA coordinate to take as offset end

        """
        if start: start = self.proteinpos2dnapos(start) - self.startPY
        else:     start = 0
        if end:   end   = self.proteinpos2dnapos(end) - self.startPY
        else:     end   = self.endPY - self.startPY

        scores = []
        if start+(aa_window_size*3) > end:
            points  = [ self._tcode_score_list[coord] for coord in range(start,end) ]
            avtcode = sum(points) / float(len(points))
            scores.append( avtcode )
        else:
            for rel_offset_start in range(start,end-(aa_window_size*3),3):
                rel_offset_end = rel_offset_start + (aa_window_size*3)
                points  = [ self._tcode_score_list[coord] for coord in range(rel_offset_start,rel_offset_end) ]
                avtcode = sum(points) / float(len(points))
                scores.append( avtcode )

        if min(scores) < min_tcode_score:
            return True, min(scores)
        else:
            return False, min(scores)

    # end of function has_lowscroring_tcode_stretch

# end of class TcodeOrf

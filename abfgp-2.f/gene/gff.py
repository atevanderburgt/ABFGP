"""
Basic GFF functionality for gene entity objects in the python gene package
"""
# Module metadata
__authors__ = "Ate van der Burgt"
__license__ = "MIT"

class BasicGFF:
    """ """
    def __init__(self):
        self._gff = self._get_gff_defaults()
        self.phase = "."
    # end of function __init__


    def __str__(self):
        """ """
        try:
            return "<%s %s>" % ( self.__class__.__name__,self.pos)
        except:
            return "<%s %s-%s>" % ( self.__class__.__name__,self.start,self.end)


    # end of function __str__


    def dna_range(self):
        """ """
        return range(self.start,self.end+1)

    # end of function dna_range


    def _get_gff_defaults(self):
        """ """
        # return dict({}) to prevent call-by-reference
        return dict({
                'fref'      : None,
                'fsource'   : 'undefined',
                'fstrand'   : '+',
                'fscore'    : '.',
                'column9data' : {},
            })

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
        if not self._gff.has_key('gname'):   self._gff['gname']   = str(self.start)
        if not self._gff.has_key('fstart'):  self._gff['fstart']  = self.start+1
        if not self._gff.has_key('fstop'):   self._gff['fstop']   = self.end
        if not self._gff.has_key('fmethod'): self._gff['fmethod'] = self.__class__.__name__
        if hasattr(self,"pssm_score") and getattr(self,"pssm_score") != None:
            self._gff['fscore'] = self.pssm_score
            if type(self._gff['fscore']) == type(float()):
                self._gff['fscore'] = round(self._gff['fscore'],2)

        column9string = ""
        # take care for the formatting of column9data
        if self._gff.has_key('column9data') and self._gff['column9data']:
            tmp = []
            for key,value in self._gff['column9data'].iteritems():
                if type(value) == type([]):
                    valuelist = value
                else: 
                    valuelist = [ value ]
                for v in valuelist: 
                    if str(v).find(" ") > -1:  v = "'%s'" % v
                    tmp.append( "; %s %s" % (key,v) )
            column9string = "".join(tmp)

        return (
            self._gff['fref'],
            self._gff['fsource'],
            self._gff['fmethod'],
            self._gff['fstart'],
            self._gff['fstop'],
            # TODO!! pssm_score and _gff['fscore'] must be converged
            # TODO!! if float value: round(score,2)
            self._gff['fscore'],
            self._gff['fstrand'],
            self.phase,
            "%s %s%s" % ( self._gff['gclass'], self._gff['gname'], column9string )
        )

    # end of function togff    

# end of class BasicGFF



class GffWithSequenceFunctionality:
    """ @attention: only inherit when the object has a dnasequence() function """
    
    def get_at_ratio(self):
        """ """
        seq = self.dnasequence()
        at = seq.count("A")+seq.count("a")+seq.count("T")+seq.count("t")
        return float(at)/len(seq)

     # end of function get_at_ratio


    def get_gc_ratio(self):
        """ """
        seq = self.dnasequence()
        gc = seq.count("G")+seq.count("g")+seq.count("C")+seq.count("c")
        return float(gc)/len(seq)

     # end of function get_gc_ratio

# end of class GffWithSequenceFunctionality

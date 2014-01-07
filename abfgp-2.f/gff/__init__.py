"""
Currect basal GFF functionality (GFF files and GFF features) used from
http://sourceforge.net/projects/pythongff

Other - much more elaborate - implementations of GFF in Python
http://www.biopython.org/DIST/docs/api/Bio.GFF-pysrc.html

IMPORTANT:

Coordinate system is in python coordinates
                GFF     Python
        start   x       x-1
        end     y       y

"""

# Imports
from pythongff.gff import Feature as PythonGffFeature
from pythongff.gff import GFF
from gff_functions import *
from os import system

# Global variables (for load_gff.pl)
PERL_PATH     = 'perl'
LOAD_GFF_PATH = '/home/avdb/code/abfgp-dev/gff/load_gff.pl'    


# Exceptions
class WrongPatternLengthApplied(Exception):
    """ """
    pass

class WrongPhaseApplied(Exception):
    """ """
    pass

class InproperlyAppliedArgument(Exception):
    """ """
    def __init__(self, message):
        self.message = message
    def __str__(self):
        return repr(self.message)


class GffFeature(PythonGffFeature):
    """
    """
    def __init__(self,start,end,strand="+",score=".",frame=".",
        fref=None,fsource=None,fmethod=None,attribute_string="",gff={}):
        """
        @type  start: number
    	@param start: start coord of Feature in PYTHON COORDINATE (GFF-1)

        @type  end: number
    	@param end: end coord of Feature

        @type  strand: string
    	@param strand: ['+','-','.'] or None, default '+'

        @type  score: float
    	@param score: score of Feature, default "."

        @type  frame: number
    	@param frame: [0,1,2,'.'] or None, default "."

        @type  phase: USE FRAME IN STEAD OF PHASE
    	@param phase: USE FRAME IN STEAD OF PHASE

        @type  fref: string
    	@param fref: default None ('None'), to be set lateron with seqname

        @type  fsource: string
    	@param fsource: default None ('None'), mostly inherited

        @type  fmethod: string
    	@param fmethod: default None ('None'), mostly inherited

        @type  gff: dictionary
    	@param gff: dictionary with (default) gff data. Possible keys are:
                    'fref', 'fmethod', 'fstrand', 'fsource', 'gclass', 'gname'

        @attention: TODO gff dictionary application is not coded yet! 

        """
        if fsource==None: fsource = self.__class__.__name__
        if fmethod==None: fmethod = self.__class__.__name__

        gfflist = [fref,fsource,fmethod,start,end,score,strand,frame]
        PythonGffFeature.__init__( self,"\t".join([str(elem) for elem in gfflist]) )

    # end of function __init__


    def __str__(self):
        """ Human-readable xml-like representation of this Feature.
            This function is likely to be overwritten in inheriting Features.
        """
        scorestring = ""
        if self.score() != ".":
            scorestring = " (score=%2.3f)" % self.score()
        return "<%s %s-%s%s>" % (
            self.__class__.__name__,
            self.start(), self.end(),
            scorestring 
            )

    # end of function __str))


    def togff(self):
        """ Original __str__ function of Pyhtongff.gff.feature.Feature """
        if not self.gclass():  self.gclass(self.__class__.__name__)
        if not self.gname():   self.gname(str(self.start()+1))
        return PythonGffFeature.__str__(self)

    # end of function togff


    ## access functions ##
    def fref(self,name = ""):
        return self.seqname(name)

    def fsource(self,name = ""):
        return self.source(name)

    def fmethod(self,name = ""):
        return self.feature(name)

    def gclass(self,name = ""):
        return self.attribute(name)

    def gname(self,name = ""):
        return self.comment(name)


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
        if not self._gff.has_key('gname'):   self._gff['gname']   = str(self.pos+1)
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

# end of class Feature


def dbcleanup(dbname=None,dbhost=None,dbconnection=None,dbuser=None,dbpass=None,fref=None):
    """
    """
    if not fref: return False
    if not dbname:
        DB_NAME = 'abfgp'
    else:
        DB_NAME = dbname
    DB_HOST       = '##########'
    DB_USER       = '##########'
    DB_PASS       = '##########'
    sql = "DELETE FROM fdata WHERE FREF = '%s'" % (fref)
    import MySQLdb
    connection = MySQLdb.connect(host=DB_HOST, user=DB_USER, db=DB_NAME,passwd=DB_PASS)
    cursor = connection.cursor(MySQLdb.cursors.DictCursor)
    _ret = cursor.execute(sql)
    return True

# end of function dbcleanup


def gff2db(dbname=None,dbhost=None,dbconnection=None,dbuser=None,dbpass=None,
    fnamefasta=None,fnamegff=None,
    pathperl=PERL_PATH,pathloadgff=LOAD_GFF_PATH):
    """
    Load a gff file (and its fasta) into a GGB database
    """
    # load to GGB database
    # example commandline
    # perl /home/avdb/code/abfgp-2.0/gff/load_gff.pl --dsn dbi:mysql:<DB_NAME>:<DB_HOST> --user <DB_USER> --pass <DB_PASS> $fname --fasta $fasta_fname;
    # DB_CONNECTION = 'dbi:mysql'
    # DB_NAME       = '##########'
    # DB_HOST       = '##########'
    # DB_USER       = '##########'
    # DB_PASS       = '##########'
    # PERL_PATH     = 'perl'
    # LOAD_GFF_PATH = '/home/avdb/data/abfgp-dev/gff/load_gff.pl'
    # make the command line to load_gff.pl
    command = "%s %s --dsn %s:%s:%s --user %s --pass %s %s" % (
        pathperl,pathloadgff,
        dbconnection,dbname,dbhost,dbuser,dbpass,
        fnamegff
        )
    # check if fasta filename is applied
    if fnamefasta:
        command = "%s --fasta %s" % (command,fnamefasta)
    # os.system -> command!
    try:
        system(command)
        return True
    except:
        return False

# end of function gff2db

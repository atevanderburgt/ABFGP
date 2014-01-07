"""
(General) functions for logging used in Aligment Based Gene Prediction
"""

# Module metadata
__authors__ = "Ate van der Burgt"
__license__ = "MIT"

VERBOSE = False
QUIET   = False
SILENT  = False

# N.B. various function have been removed after the development of ABFGP was done:
# logfGSG(GSG,**kwargs)
# printdetailedGSGreport(GSG)

def logf(*args):
    """
    Log results of intermediate steps to STDOUT
    Change this function when more/less information is required/needed
    Major milestone messages start with a '##' symbol
    Milestones messages start with a '#' symbol
    """
    if SILENT: pass
    else:
        x = " ".join([str(item) for item in args])
        if VERBOSE:
            print x
        else:
            if not x:
                pass
            elif x[0:2] == "##":
                print x
            elif x[0] == "#" and x[1] != "#":
                if not QUIET:
                    print x
            else:
                pass
# end of function logf


class Logger:
    def __init__(self,verbose=False,quiet=False,silent=False):
        self.VERBOSE = verbose
        self.SILENT = silent
        self.QUIET = quiet

    def logf(self,*args):
        """
        Log results of intermediate steps to STDOUT
        Change this function when more/less information is required/needed
        Major milestone messages start with a '##' symbol
        Milestones messages start with a '#' symbol
        """
        if self.SILENT: pass
        else:
            x = " ".join([str(item) for item in args])
            if self.VERBOSE:
                print x
            else:
                if not x:
                    pass
                elif x[0:2] == "##":
                    print x
                elif x[0] == "#" and x[1] != "#":
                    if not self.QUIET:
                        print x
                else:
                    pass

# end of class Logger

def debug(crossdata,input):
    """
    Debugging function feeded with `crossdata` and `input`.
    After each step, it can be checked what is in what-where-how.
    If nothing needed, return True immediately!
    """
    return True
    # below here, type functionality for debugging
    for k,v in crossdata[('fgsg', 'mgg')].iteritems():
        if type(v) != type({}): continue
        for (a,b,c,d),pacbp in v.iteritems():
            if c == 85 and d == 1:
                print k, (a,b,c,d), pacbp

# end of function debug



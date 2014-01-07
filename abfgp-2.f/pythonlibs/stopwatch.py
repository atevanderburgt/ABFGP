"""
StopWatch class that enables timing of range of processes (steps and total time)
"""

# Imports (Python)
from time import time

class AlgorithmStep:
    def __init__(self,name="",substep=list("abcdefghijklmnopqrstuvwxyz")):
        """ """
        self.name     = name
        self._substep = substep
        self._pos     = 0
        if not self.name:
            self.name = self.__class__.__name+".%s"
    # end of function __init__

    def __str__(self):
        """ """
        try:
            stepchar = self._substep[self._pos]
        except:
            stepchar = ""
        try:
            return self.name % stepchar
        except TypeError:
            return str(self.name)+str(stepchar)
    # end of function __str__ 

    def next(self):
        """ """
        self._pos+=1
    # end of function next

# end of class AlgorithmStep 


class StopWatch:
    def __init__(self,name="",stepname=None,RETURN_ON_PRESS=True,PRINT_ON_PRESS=False):
        """ """
        self._laptime = []
        self._start   = None
        self.name     = str(name)
        self.stepname = stepname
        if not name:
            self.name = self.__class__.__name__
        self.PRINT_ON_PRESS  = PRINT_ON_PRESS
        self.RETURN_ON_PRESS = RETURN_ON_PRESS
    # end of function __init__

    def __len__(self):
        return len(self._laptime)
    # end of function __len__

    def __str__(self,*args):
        """ """
        # format self.stepname; is it a string or a AlgorithmStep instance?
        try:
            stepname = str(self.stepname)
            stepname = " %s " % stepname
            self.stepname.next()
            # stepname isa AlgorithmStep!
        except:
            if self.stepname == None:
                stepname = " "
            else:
                stepname = " %s " % str(self.stepname)
        # make the provided message
        message = " ".join([str(arg) for arg in args])
        if message: message = " "+message
        return "<%s %7s %7s ( % 3s ) %s>%s" % (
                self.name,
                "%3.3f" % self.gettime_lap(),
                "%5.1f" % self.gettime_total(),
                len(self)-1,
                stepname,
                message,
                )
    # end of function __str__

    def start(self,*args):
        """ """
        self._start   = time()
        self._laptime = [ 0.0 ]
        if self.PRINT_ON_PRESS:
            print self.__str__(*args)
        if self.RETURN_ON_PRESS:
            return self.__str__(*args)
    # end of function start

    def stop(self,*args):
        """ """
        return self._add_time(*args) 
    # end of function stop

    def lap(self,*args):
        """ """
        return self._add_time(*args)
    # end of funtion lap

    def _add_time(self,*args):
        """ """
        if self._start == None:
            # user forgot to press on start
            self.start(*args)
        self._laptime.append( time()-self._start )
        if self.PRINT_ON_PRESS:
            print self.__str__(*args)
        if self.RETURN_ON_PRESS:
            return self.__str__(*args)
        else:
            pass # return nothing! 
    # end of function _add_time

    def _get_lap_time(self,lapnumber):
        """ """
        return self._laptime[lapnumber] - self._laptime[lapnumber-1]
    # end of function _get_lap_time

    def gettime_lap(self,*args):
        """ """
        if not args: return self._get_lap_time(len(self)-1)
        else:        return self._get_lap_time(args[0])
    # end of function gettime_lap

    def gettime_total(self,*args):
        """ """
        if not args: return self._laptime[-1]
        else:        return self._laptime[args[0]]
    # end of function gettime_total

# end of class Stopwatch


class SilentStopWatch(StopWatch):
    """ Default silent version of the StopWatch object """
    def __init__(self,*args,**kwargs):
        StopWatch.__init__(self,*args,**kwargs)
        self.RETURN_ON_PRESS=False
        self.name="StopWatch"

    # end of funcion __init__

# end of class SilentStopWatch


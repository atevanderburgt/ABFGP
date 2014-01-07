import sys
from StringIO import StringIO

class stdoutManagerClass:
    def __init__(self):
        self._orig_stdout = None
        self._buffered_stdout = None
    def buffer(self):
        """ Catch stdout into a StringIO object """
        self._orig_stdout = sys.stdout
        self._buffered_stdout = StringIO()
        sys.stdout = self._buffered_stdout
    def reset(self):
        """ Reset stdout to its original state """
        sys.stdout = self._orig_stdout
        self._buffered_stdout = None
    def clean(self):
        """ Clean all data currently gathered in the buffered stdout in StringIO() """
        if self._buffered_stdout:
            # overwrite asa novel empty StringIO() object
            self._buffered_stdout = StringIO()
    def log(self):
        """ Print (not clean) all data currently gathered in the buffered stdout in StringIO() """
        if self._buffered_stdout:
            sys.stdout = self._orig_stdout 
            print self._buffered_stdout.getvalue()
            sys.stdout = self._buffered_stdout
    def getlogtxt(self):
        """ Get (not clean) all data currently gathered in the buffered stdout in StringIO() """
        txt = ""
        if self._buffered_stdout:
            sys.stdout = self._orig_stdout
            txt = self._buffered_stdout.getvalue()
            sys.stdout = self._buffered_stdout
        return txt
    def logclean(self):
        """ Print and successively clean all data currently gathered in the buffered stdout in StringIO() """
        self.log()
        self.clean()
    def getlogtxtclean(self):
        """ Get and successively clean all data currently gathered in the buffered stdout in StringIO() """
        txt = self.getlogtxt()
        self.clean()
        return txt
# end of class stdoutManagerClass


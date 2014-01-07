"""
"""

# Module metadata
__authors__ = "Ate van der Burgt"
__license__ = "MIT"

# Orf Imports
from basicorf import BasicOrf
from codingorf import CodingOrf
from tcodeorf import TcodeOrf

# Python Imports

# Global variable Imports

class TcodeCodingOrf(TcodeOrf,CodingOrf,BasicOrf):

    def __init__(self,*args,**kwargs):
        """
        """
        BasicOrf.__init__(self,*args,**kwargs)
        TcodeOrf.__init__(self,*args,**kwargs)
        CodingOrf.__init__(self,*args,**kwargs)

    # end of function __init__

# end of class TcodeCodingOrf

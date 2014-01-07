"""
Exceptions used in parser scripts & functions
"""

# Module metadata
__authors__ = "Ate van der Burgt"
__license__ = "MIT"

# Exceptions

class InproperlyAppliedArgument(Exception):
    """ """
    def __init__(self, *args):
        self.message = " ".join([str(item) for item in args])
    def __str__(self):
        return repr(self.message)


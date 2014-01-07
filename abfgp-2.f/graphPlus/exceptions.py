"""
Exceptions for graphPlus class
"""

# Module metadata
__authors__ = "Ate van der Burgt"
__license__ = "MIT"

class InproperlyAppliedArgument(Exception):
    """ """
    def __init__(self, message):
        self.message = message
    def __str__(self):
        return repr(self.message)


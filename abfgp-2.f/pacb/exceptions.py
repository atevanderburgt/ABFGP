"""
Exceptions for python-pacb.
"""

# Module metadata
__authors__ = "Ate van der Burgt"
__license__ = "MIT"

# Exceptions

class CoordinateOutOfRange(Exception):
    """ """
    pass

class ZeroSizedPacb(Exception):
    """ """
    pass

class InproperlyAppliedArgument(Exception):
    """ """
    def __init__(self, *args):
        self.message = " ".join([str(item) for item in args])
    def __str__(self):
        return repr(self.message)


class UnidenticalOrfs(InproperlyAppliedArgument):
    """ unidentical Orf(s) in PacbPORF merging/comparison """
    pass

class NoPacbPORF(InproperlyAppliedArgument):
    """ attribute is not a PacbPORF object """
    pass

class NoPacbPDNA(InproperlyAppliedArgument):
    """ attribute is not a PacbP object """
    pass

class NoPacbP(InproperlyAppliedArgument):
    """ attribute is not a PacbP object """
    pass

class NoOrf(InproperlyAppliedArgument):
    """ attribute is not an Orf object """
    pass

"""
Exceptions for the python gene package
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

class MissingAttributeValue(InproperlyAppliedArgument):
    """ Function's attribute or variable which is required, is not given """
    pass

class WrongPatternLengthApplied(Exception):
    """ incorrect pattern length applied for a PSSM element """
    pass

class WrongStrandApplied(Exception):
    """ incorrect strand applied (not + or -) """
    pass


class WrongPhaseApplied(Exception):
    """ incorrect or no phase applied """
    pass

class UnexpectedSpliceSitePhase(Exception):
    """ incorrect or no splice site phase applied """
    pass

class IncompatibleSpliceSitePhases(Exception):
    """ Splice site phases do not match """
    pass

class NoOrf(InproperlyAppliedArgument):
    """ attribute is not an Orf object """
    pass

class NoDonor(InproperlyAppliedArgument):
    """ attribute is not a Donor object """
    pass

class NoAcceptor(InproperlyAppliedArgument):
    """ attribute is not an Acceptor object """
    pass

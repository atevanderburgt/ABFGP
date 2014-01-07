"""
Exceptions in use if the gff class GffFeature and Pythongff classes
Currect basal GFF functionality (GFF files and GFF features) used from
http://sourceforge.net/projects/pythongff
"""

# Exceptions

class WrongPatternLengthApplied(Exception):
    """ """
    pass


class WrongLengthApplied(Exception):
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


class WrongStrandApplied(InproperlyAppliedArgument):
    """ """
    def __init__(self):
        self.message = "Strand must be any of '+','-','.'"


class IncompatibleSpliceSitePhases(Exception):
    """ """
    pass


class UnexpectedSpliceSitePhase(Exception):
    """ """
    pass


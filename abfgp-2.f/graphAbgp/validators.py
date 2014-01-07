"""
Validators for graphAbgp classes in Aligment Based Gene Prediction
"""
from exceptions import InproperlyAppliedArgument

def IsProperPhaseValidator(phase):
    """ """
    if phase not in [0,1,2]:
        message = "phase is not in [0,1,2], but '%s'" % phase
        raise InproperlyAppliedArgument, message

# end of function IsProperPhaseValidator

def IsCodingBlockEndValidator(cbgend):
    """ """
    if cbgend.__class__.__name__ != 'CodingBlockEnd':
        message = "first argument `cbgend` is not a CodingBlockEnd, but a `%s`" % cbgend.__class__.__name__
        raise InproperlyAppliedArgument, message

# end of function IsCodingBlockEndValidator

def IsCodingBlockStartValidator(cbgstart):
    """ """
    if cbgstart.__class__.__name__ != 'CodingBlockStart':
        message = "first argument `cbgstart` is not a CodingBlockStart, but a `%s`" % cbgstart.__class__.__name__
        raise InproperlyAppliedArgument, message

# end of function IsCodingBlockStartValidator


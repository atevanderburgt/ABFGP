"""
Validators for the python gene package
"""
from gene_exceptions import *

def IsOrf(obj):
    """ Is applied attribute an Orf object? """
    orfs = ['Orfs','BasicOrf','TcodeOrf','CodingOrf','TcodeCodingOrf']
    if obj.__class__.__name__ not in orfs:   
        raise NoOrf, obj

# end of function IsOrf


def IsDonor(obj):
    """ Is applied attribute a Donor object? """
    if obj.__class__.__name__.find("SpliceDonor") != 0:
        raise NoDonor, obj

# end of function IsDonor


def IsProperStrand(strand):
    """ Is strand a recognized and proper strand (+ or -) """
    if strand not in ['+','-']:
        raise WrongStrandApplied

# end of function IsProperStrand


def IsProperIntronPhase(phase):
    """ Is phase a proper intron phase? """
    if phase not in [None,0,1,2]:
        raise WrongPhaseApplied, "phase (%s) not in [None,0,1,2]" % phase

# end of function IsProperIntronPhase


def IsProperPhase(phase):
    """ Is phase (gff column 8) properly applied? """
    if phase not in ['.',None,0,1,2]:
        raise WrongPhaseApplied, "phase (%s) not in ['.',None,0,1,2]" % phase

# end of function IsProperPhase


def IsAcceptor(obj):
    """ Is applied attribute an Acceptor object? """
    if obj.__class__.__name__.find("SpliceAcceptor") != 0:
        raise NoAcceptor, obj

# end of function IsAcceptor


def IsSequenceErrorCoordinate(obj):
    """ Is applied attribute a SequenceErrorCoordinate object? """
    if obj.__class__.__name__ != "SequenceErrorCoordinate":
        message = "No SequenceErrorCoordinate but '%s'"  % obj.__class__.__name__
        raise InproperlyAppliedArgument, message

# end of function IsSequenceErrorCoordinate

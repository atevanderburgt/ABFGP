"""
Validators for python-pacb.
"""
from exceptions import *

def IsIdenticalOrfs(orfA,orfB):
    """ Are Orfs identical? """
    if orfA.id != orfB.id:
        raise UnidenticalOrfs, "Orf objects not identical"
    elif not orfA.id and not orfB.id and str(orfA) != str(orfB):
        # no id given; compare strings
        raise UnidenticalOrfs, "Orf objects not identical"
    else:
        pass

# end of function IsIdenticalOrfs


def IsPacbPORF(pacbpobj):
    """ Is applied attribute a PacbPORF object? """
    if pacbpobj.__class__.__name__ != "PacbPORF":
        raise NoPacbPORF, pacbpobj

# end of function IsPacbPORF


def IsPacbP(pacbpobj):
    """ Is applied attribute a PacbP object? """
    if pacbpobj.__class__.__name__ != "PacbP":
        raise NoPacbP, pacbpobj

# end of function IsPacbP


def IsPacbPDNA(pacbpobj):
    """ Is applied attribute a PacbPDNA (or PacbPORF) object? """
    if pacbpobj.__class__.__name__ not in ["PacbPDNA","PacbPORF"]:
        raise NoPacbPDNA, pacbpobj

# end of function IsPacbPDNA


def IsOrf(obj):
    """ Is applied attribute an Orf object? """
    orfs = ['Orfs','BasicOrf','TcodeOrf','CodingOrf','TcodeCodingOrf']
    if obj.__class__.__name__ not in orfs: 
        raise NoOrf, obj

# end of function IsOrf

"""
(General) Exceptions used in Aligment Based Gene Prediction
"""

# Module metadata
__authors__ = "Ate van der Burgt"
__license__ = "MIT"


class InproperlyAppliedArgument(Exception):
    """ """
    def __init__(self, *args):
        self.message = " ".join([str(item) for item in args])
    def __str__(self):
        return repr(self.message)

class MissingArgumentException(InproperlyAppliedArgument):
    pass

class NoCrossdataApplied(MissingArgumentException):
    pass

class NoInputApplied(MissingArgumentException):
    pass

class NoCodingBlockGraphApplied(MissingArgumentException):
    """ No attribute applied or attribute is not a CBG """
    pass

class NoGenestructureOfCodingBlockGraphsApplied(MissingArgumentException):
    """ No attribute applied or attribute is not a GSG """
    pass

class NoGeneTreeGraphApplied(MissingArgumentException):
    """ No attribute applied or attribute is not a GTG """
    pass

class OrganismNotPresentInGraph(InproperlyAppliedArgument):
    """ Organism Identifier is absent in this graph """
    pass

class NodeNotPresentInGraph(Exception):
    """ Node Identifier is absent in this graph """
    pass

class EdgeNotPresentInGraph(Exception):
    """ Requested Edge is absent in this graph """
    pass

class NoOverallMinimalSpanningRange(InproperlyAppliedArgument):
    """ CodingBlockGraph has no Overall Minimal Spanning Range (OMSR) """
    pass

class UnexptectedAttributeDataFormat(InproperlyAppliedArgument):
    """ Function's attribute or variable is not of correct data type or format """
    pass

class MissingAttributeValue(Exception):
    """ Function's attribute or variable which is required, is not given """ 
    def __init__(self, *args):
        self.message = " ".join([str(item) for item in args])
    def __str__(self):
        return repr(self.message)

class IncompatibleSpliceSitePhases(Exception):
    pass

class InValidStringRepresentationOfNode(Exception):
    """ """
    def __init__(self, *args):
        self.message = " ".join([str(item) for item in args]) 
    def __str__(self):
        return repr(self.message)


class InproperGeneStructureOfCodingBlockGraphs(Exception):
    """ Inproperly ordered CBGs and/or lsrCBGs in a GeneStructureOfCodingBlockGraphs """
    def __init__(self, *args):
        self.message = " ".join([str(item) for item in args])
    def __str__(self):
        return repr(self.message)


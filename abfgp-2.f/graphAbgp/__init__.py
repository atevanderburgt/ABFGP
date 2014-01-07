# Copyright (c) 2008 Ate van der Burgt <ate.vanderburgt@wur.nl>
#
# Permission is hereby granted, free of charge, to any person
# obtaining a copy of this software and associated documentation
# files (the "Software"), to deal in the Software without
# restriction, including without limitation the rights to use,
# copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the
# Software is furnished to do so, subject to the following
# conditions:

# The above copyright notice and this permission notice shall be
# included in all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
# OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
# HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
# WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
# OTHER DEALINGS IN THE SOFTWARE.


"""
graph classes for Aligment Based Gene Predictions
All classes inherit from graphPlus, a Plus-version of the python-graph package of Pedro Matiello

class BasalSiteAlignmentFunctions
class BasalPSSMObjectGraphFunctions
class OrganismGraph(graphPlus.graphPlus)
class CodingBlockGraph(OrganismGraph)
class AlignedStopCodonGraph(OrganismGraph)
class AlignedStartCodonGraph(OrganismGraph,BasalSiteAlignmentFunctions)
class AlignedPSSMTSSGraph(AlignedStartCodonGraph,BasalSiteAlignmentFunctions,BasalPSSMObjectGraph)
class AlignedSpliceSiteGraph(OrganismGraph,BasalSiteAlignmentFunctions,BasalPSSMObjectGraph)
class AlignedDonorSiteGraph(AlignedSpliceSiteGraph)
class AlignedAcceptorSiteGraph(AlignedSpliceSiteGraph)


"""


# Module metadata
__authors__ = "Ate van der Burgt"
__license__ = "MIT"


# Imports
from graph_organism import OrganismGraph
from graph_genetree import GeneTreeGraph
from graph_pacbpcollection import PacbpCollectionGraph
from graph_codingblock import CodingBlockGraph
from graph_inwardspointingcodingblock import InwardsPointingCodingBlockGraph 
from graph_lowsimilaritycodingblock import LowSimilarityRegionCodingBlockGraph
from graph_genestructure import GenestructureOfCodingBlockGraphs
from graph_pssmcollections import *
from graph_alignedsites import *
from graph_exoncollection import ExonCollectionGraph

#from conversion import *

###from sitealignment import BasalSiteAlignmentFunctions
###from pssmobjects import BasalPSSMObjectGraphFunctions
###from tcodedata import GraphTcodeDataAccesFunctions
###from codingblockcollectionharvesting import *
###from codingblockprinting import *
###import codingblocksplitting

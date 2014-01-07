"""
Functions for identification of complete K(s) graphs in a larger network for python-graph.
"""

# Module metadata
__authors__ = "Ate van der Burgt"
__license__ = "MIT"

# Imports
from exceptions import InproperlyAppliedArgument


def makecompletegraph(gra,wt=0.0):
    """
    Make the graph complete, i.e. set all non-existing edges to wt(=0.0)

    @type  gra: graph
    @param gra: Graph object

    @type  wt: * (presumably/preferably float)
    @param wt: weight to assign to new edges

    @attention: only 100% correct when no arrows are present
    """
    if gra.connectivitysaturation() < 1.0:
        combis = []
        for n1 in gra.get_nodes():
            for n2 in gra.get_nodes():
                if n1 == n2: continue
                if (n1,n2) in combis or (n2,n1) in combis: continue
                # set this edge if it does not exist
                if not gra.has_edge(n1,n2):
                    gra.add_edge(n1,n2,wt=wt)

# end of function makecompletegraph


def is_complete(gra):
    """
    Is this graph complete / fully connected ?

    @type  gra: graph
    @param gra: Graph object

    @rtype:  Boolean
    @return: True or False
    """
    if gra.connectivitysaturation() == 1.0:
        return True
    else:
        return False

# end of function is_complete


def remove_incomplete_components(gra,s=None):
    """
    Remove all nodes and edges that are not part of a Complete K(s) graph

    @type  gra: graph
    @param gra: Graph object

    @type  s: integer
    @param s: integer value (gte 3)

    @attention: function is overwritten in OrganismGraph class (a coloured graph)
    """
    # some input sanity checks
    if type(s) != type(int()):
        raise InproperlyAppliedArgument, "variable 's' not integer but '%s'" % type(s)
    if s <= 2:
        raise InproperlyAppliedArgument, "variable 's' not gte 3 but %s" % s

    # delete all the nodes with a to low degree
    gra.minimal_node_degree(s-1)

# end of function remove_incomplete_components


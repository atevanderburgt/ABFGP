"""
Subclass with functions for CodingBlockGraphs (etc...)
Related to Maximal Spanning Range (MAXSR)
"""

# Module metadata
__authors__ = "Ate van der Burgt"
__license__ = "MIT"

# Python Imports
from sets import Set

# graphAbgp import
from graphAbgp.exceptions import OrganismNotPresentInGraph

# MINSR imports
from MINSR import (
    _size_dict,
    _min_dict,
    _max_dict,
    _calculate_cbg_distance,
    )

class MaximalSpanningRange:
    """
    Subclass with MaximalSpanningRange functions for CodingBlockGraphs
    
    @attention: not functional on its own, required to be inherited from
    @attention: requires MINSR._get_ranges_by_nodes() function
    """

    def maximal_spanning_range(self,organism=None,node=None):
        """
        TODO description of maximal_spanning_range

        @attention: only functional if pacbp's are present in the graph
                    (perform the function `harvest_pacbps_from_crossdata`)

        @attention: only meaningfull for (nearly) Fully Connected Graphs

        @type  organism: *
        @param organism: organism identifier (or None)

        @type  node: *
        @param node: Node identifier (or None)

		@rtype:  dictionary (or Set if organism is specified)
		@return: dictionary with nodes (keys) and spanning range Sets
                 (values), or only the spanning range Set if an organism
                 or node identifier was specified
        """
        if organism and organism not in self.organism_set():
            raise OrganismNotPresentInGraph
        if node and node not in self.get_nodes():
            raise NodeNotPresentInGraph
        # try to get from the objects' cached attributes
        if self._maxsr:
            if organism:    return self._maxsr[self.node_by_organism(organism)]
            elif node:      return self._maxsr[node]
            else:           return self._maxsr

        # create Maximal Spanning Range dictionary
        ranges = {}
        for (node1,node2) in self.pairwisecrosscombinations_node():
            # get ranges for these nodes
            thisrangeQ, thisrangeS = self._get_ranges_by_nodes(node1,node2)
            # set to ranges
            if not ranges.has_key(node1):
                ranges[node1] = thisrangeQ
            else:
                ranges[node1].union_update(thisrangeQ)
            if not ranges.has_key(node2):
                ranges[node2] = thisrangeS
            else:
                ranges[node2].union_update(thisrangeS)

        # return all ranges or the range of a specific organism
        if organism:
            for node, orgrange in ranges.iteritems():
                if self._organism_from_node(node) == organism:
                    return orgrange
        else:
            return ranges

    # end of function maximal_spanning_range


    def has_maximal_spanning_range(self):
        """
        Has this (nearly Fully Connected) Graph a minimal spanning range?

        @attention: only functional if pacbp's are present in the graph
                    (perform the function `harvest_pacbps_from_crossdata`)

        @attention: only meaningfull for (nearly) Fully Connected Graphs
        @attention: in practice, always True

		@rtype:  Boolean
		@return: True or False
        """
        if self.maximal_spanning_range_sizes().values():
            if min(self.maximal_spanning_range_sizes().values()) == 0:
                return False
            else:
                return True
        else:
            return False

    # end of function has_maximal_spanning_range


    def has_maxsr(self):
        """
        @attention: alias function name of has_maximal_spanning_range
        """
        return self.has_maximal_spanning_range()

    # end of function has_msr


    def maximal_spanning_range_sizes(self):
        """
        Get the sizes of the maximal spanning range of all pacbp's for each node

        @attention: only functional if pacbp's are present in the graph
                    (perform the function `harvest_pacbps_from_crossdata`)

        @attention: only meaningfull for (nearly) Fully Connected Graphs

		@rtype:  dictionary
		@return: dictionary with nodes (keys) and size of spanning range Sets (values)
        """
        return _size_dict(self.maximal_spanning_range())

    # end of function maximal_spanning_range_sizes


    def maxsr_starts(self):
        """
        Return a dictionary of MAXSR start coodinates (keys==nodes)

        @rtype:  dictionary 
        @return: dict[Node] = MAXSR start coords (integer) as values
        """
        return _min_dict(self.maximal_spanning_range())

    # end of function maxsr_starts


    def maxsr_ends(self):
        """
        Return a dictionary of MAXSR end coodinates (keys==nodes)

        @rtype:  dictionary 
        @return: dict[Node] = MAXSR end coords (integer) as values
        """
        return _max_dict(self.maximal_spanning_range())

    # end of function maxsr_ends


    def maxsr_distance_between_codingblocks(self,other,**kwargs):
        """
        Distance in AA between the MAXSR of two CodingBlockGraphs

        @attention: see MINSR._calculate_cbg_distance() for **kwargs
        """
        coords1 = self.maximal_spanning_range()
        coords2 = other.maximal_spanning_range()
        return _calculate_cbg_distance(self,other,coords1,coords2,**kwargs)

    # end of function maxsr_distance_between_codingblocks

# end of class MaximalSpanningRange
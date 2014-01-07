"""
Subclass with functions for CodingBlockGraphs (etc...)
Related to Overall Minimal Spanning Range (OMSR)
"""

# Module metadata
__authors__ = "Ate van der Burgt"
__license__ = "MIT"

# Python Imports
from sets import Set
from copy import deepcopy

# graphAbgp import
from graphAbgp.exceptions import OrganismNotPresentInGraph

# MINSR imports
from MINSR import (
    _size_dict,
    _min_dict,
    _max_dict,
    _calculate_cbg_distance,
    )

class OverallMinimalSpanningRange:
    """
    Subclass with OverallMinimalSpanningRange functions for CodingBlockGraphs
    
    @attention: not functional on its own, required to be inherited from
    @attention: requires MINSR._get_ranges_by_nodes() function
    """

    def overall_minimal_spanning_range(self,organism=None,node=None):
        """
        TODO description of overall_minimal_spanning_range

        @attention: only functional if pacbpORFs (not Pacbps) are present in the graph
                    (perform the function `harvest_pacbps_from_crossdata`)

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
        if self._omsr:
            if organism:    return self._omsr[self.node_by_organism(organism)]
            elif node:      return self._omsr[node]
            else:           return self._omsr

        # FIRST iteration; get normal minimal spanning range
        initial_ranges = self.minimal_spanning_range()
        final_ranges = deepcopy(initial_ranges)
        # SECOND iteration: get minimal spanning range by overall comparison
        for (node1,node2) in self.pairwisecrosscombinations_node():
            # check if these nodes are indeed present as an edge
            if not self.has_edge(node1,node2): continue
            # get ranges for these nodes
            thisrangeQ, thisrangeS = self._get_ranges_by_nodes(
                    node1,node2,
                    _limit_range_node1=initial_ranges[node1],
                    _limit_range_node2=initial_ranges[node2],
                    )
            # set to ranges
            final_ranges[node1].intersection_update(thisrangeQ)
            final_ranges[node2].intersection_update(thisrangeS)

        # return ranges or the range of a specific organism or node
        if organism:
            return final_ranges[self.node_by_organism(organism)]
        elif node:
            return final_ranges[node]
        else:
            return final_ranges

    # end of function overall_minimal_spanning_range


    def has_overall_minimal_spanning_range(self):
        """
        Has this (nearly Fully Connected) Graph an overall minimal spanning range?

        @attention: only functional if pacbp's are present in the graph
                    (perform the function `harvest_pacbps_from_crossdata`)

        @attention: only meaningfull for (nearly) Fully Connected Graphs

        @rtype:  Boolean
        @return: True or False
        """
        # TODO TODO temporarily set in try/except clause; this must be fixed!!
        try:
            if not self.are_all_edges_covered_by_pacbps():
                return False
            elif self.overall_minimal_spanning_range_sizes().values():
                # 8/4/09 smallest OMSR size 3 (means 4 AAs overall aligned)
                if min(self.overall_minimal_spanning_range_sizes().values()) <\
                self.MINIMAL_OVERAL_SPANNING_RANGE_SIZE:
                    return False
                else:
                    return True
            else:
                return False
        except:
            return False

    # end of function has_overall_minimal_spanning_range


    def has_omsr(self):
        """
        @attention: alias function name of has_overall_minimal_spanning_range
        """
        return self.has_overall_minimal_spanning_range()

    # end of function has_omsr


    def overall_minimal_spanning_range_sizes(self):
        """
        Get the sizes of the overall_minimal spanning range of all pacbp's for each node

        @attention: only functional if pacbp's are present in the graph
                    (perform the function `harvest_pacbps_from_crossdata`)

        @attention: only meaningfull for (nearly) Fully Connected Graphs

		@rtype:  dictionary
		@return: dictionary with nodes (keys) and size of spanning range Sets (values)

        """
        return _size_dict(self.overall_minimal_spanning_range())

    # end of function overall_minimal_spanning_range_sizes


    def omsr_starts(self):
        """
        Return a dictionary of OMSR start coodinates (keys==nodes)

        @rtype:  dictionary 
        @return: dict[Node] = OMSR start coords (integer) as values
        """
        return _min_dict(self.overall_minimal_spanning_range())

    # end of function omsr_starts


    def omsr_ends(self):
        """
        Return a dictionary of OMSR end coodinates (keys==nodes)

        @rtype:  dictionary 
        @return: dict[Node] = OMSR end coords (integer) as values
        """
        return _max_dict(self.overall_minimal_spanning_range())

    # end of function omsr_ends


    def omsr_distance_between_codingblocks(self,other,**kwargs):
        """
        Distance in AA between the OMSR of two CodingBlockGraphs

        @attention: alias for distance_between_codingblocks() function
        """
        coords1 = self.overall_minimal_spanning_range()
        coords2 = other.overall_minimal_spanning_range()
        return _calculate_cbg_distance(self,other,coords1,coords2,**kwargs)

    # end of function omsr_distance_between_codingblocks


    def omsr_identityscore(self):
        """
        Get the CBG's identityscore in the OMSR region
   
        @attention: a different measure as total_weight()
        """
        omsr = self.overall_minimal_spanning_range()
        wts = []
        for (key,nodeQ,nodeS), pacbporf in self.pacbps.iteritems():
            sta,end = min(omsr[nodeQ]), max(omsr[nodeQ])
            wt = pacbporf.identityscore_slice_by_abs_protein_query(sta,end)
            wts.append(wt)
        return float(sum(wts))/self.edge_count()
    
    # end of function omsr_identityscore


    def omsr_length(self):
        """
        Get the (smallest) OMSR length for this CBG
    
        @attention: no OMSR will return 0
    
        @rtype:  integer
        @return: smallest OMSR length of this CBG in AA coords
        """
        if self.has_overall_minimal_spanning_range(): 
            return min(self.overall_minimal_spanning_range_sizes().values())
        else:
            return 0

    # end of function omsr_length

    
    def omsrlength(self): 
        """ 
        @attention: alias for omsr_length()
        """
        return self.omsr_length()
    
    # end of function omsrlength


    def omsr2orfstart(self,organism=None,node=None):
        """
        """
        omsrends = self.omsr_ends()
        omsr2start = {}
        for _node in omsrends.keys():
            org = self.organism_by_node(_node)
            if node and _node != node: continue
            if organism and org != organism: continue
            orf = self.get_orfs_of_graph(organism=org)[0]
            if organism or node:
                omsr2start = Set(range(orf.protein_startPY,omsrends[_node]+1))
            else:
                omsr2start[_node] = Set(range(orf.protein_startPY,omsrends[_node]+1))

        # return omsr2start dict with Sets / single Set
        return omsr2start

    # end of function omsr2orfstart


    def omsr2orfend(self,organism=None,node=None):
        """
        """
        omsrstarts = self.omsr_starts()
        omsr2end = {}
        for _node in omsrstarts.keys():
            org = self.organism_by_node(_node)
            if node and _node != node: continue
            if organism and org != organism: continue
            orf = self.get_orfs_of_graph(organism=org)[0]
            if organism or node:
                omsr2end = Set(range(omsrstarts[_node],orf.protein_endPY))
            else:
                omsr2end[_node] = Set(range(omsrstarts[_node],orf.protein_endPY))

        # return omsr2end dict with Sets / single Set
        return omsr2end

    # end of function omsr2orfend

# end of class OverallMinimalSpanningRange


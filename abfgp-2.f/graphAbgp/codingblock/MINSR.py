"""
Subclass with functions for CodingBlockGraphs (etc...)
Related to Minimal Spanning Range ( MINSR / MSR )
"""

# Python Imports
from sets import Set

# graphAbgp import
from graphAbgp.exceptions import OrganismNotPresentInGraph

# Module metadata
__authors__ = "Ate van der Burgt"
__license__ = "MIT"

class MinimalSpanningRange:
    """
    Subclass with MinimalSpanningRange functions for CodingBlockGraphs
    
    @attention: not functional on its own, required to be inherited from
    """

    def _get_ranges_by_nodes(self,node1,node2,_limit_range_node1=[],_limit_range_node2=[]):
        """
        Function in use in minimal_spanning_range() and maximal_spanning_range()

        @type  node1: *
        @param node1: One Node

        @type  node2: *
        @param node2: One Node

        @type  _limit_range_node1: list
        @param _limit_range_node1: list of integers, used as allowed interval for node1

        @type  _limit_range_node2 list
        @param _limit_range_node2: list of integers, used as allowed interval for node2

		@rtype:  ( Set, Set )
		@return: tuple of 2 sets, range for Query and Sbjct

        @attention: node1 MUST be the query node, node2 MUST be the sbjct node!
        @attention: THIS FUNCTION IS ONLY 100% FUNCTIONAL FOR PacbPORF AND PacbPDNA objects!

        """
        # set known ranges by pre-limit known ranges
        thisrangeQ = Set()
        thisrangeS = Set()

        # IMPORTANT REMINDER!!! Most of the classes that use this function,
        # are (semi-)complete Graphs, with a single PacbPobj per edge.
        # In that case the iteration over ``get_pacbps_by_nodes`` will obviously
        # result in a single pacbpobj to be returned.
        # However, this function is called as well from more ``collection`` types
        # of graphs: (semi-)complete Graphs with >1 PacbPobj per edge. In the older
        # version of this function, this caused a - for long - undetected flaw in
        # the function's output. All PacbPobject of a single edge must here
        # be treated as a `single` PacbPobject. So, its range is the SUM of all
        # the individual ranges, and not the cross-section!

        for pacbpobj in self.get_pacbps_by_nodes(node1=node1,node2=node2):
            if pacbpobj.__class__.__name__ in ['PacbPORF','PacbPDNA']:
                # PacbPORF or further inhertited object.
                # The alignment might be extended,so
                # so take original alignment coordinates
                spos = pacbpobj._positions[pacbpobj._original_alignment_pos_start]
                epos = pacbpobj._positions[pacbpobj._original_alignment_pos_end-1]
                startQ, endQ = spos.query_pos, epos.query_pos + 1
                startS, endS = spos.sbjct_pos, epos.sbjct_pos + 1

            else:
                # a Pacbp or PacbpDNA object, take start & end
                startQ, endQ = pacbpobj.query_start, pacbpobj.query_end
                startS, endS = pacbpobj.sbjct_start, pacbpobj.sbjct_end

            # if _limit_range_node1 and DNA-type PacbP,
            # then limit by allowed interval _limit_range_node1 
            if pacbpobj.__class__.__name__ in ['PacbPORF','PacbPDNA']:
                if _limit_range_node1:
                    if startQ < min(_limit_range_node1):
                        spos = pacbpobj._positions[ pacbpobj.alignmentposition_by_query_pos( min(_limit_range_node1) ) ]
                        startQ, startS = spos.query_pos, spos.sbjct_pos
                    if endQ-1 > max(_limit_range_node1):
                        #print "max(_limit_range_node1)", max(_limit_range_node1)
                        #print "obj.a_b_q_p()", pacbpobj.alignmentposition_by_query_pos( max(_limit_range_node1) )
                        #print node1,node2, _limit_range_node1, _limit_range_node2
                        #print self.get_nodes()
                        #print startQ, endQ
                        #print startS, endS

                        epos = pacbpobj._positions[ pacbpobj.alignmentposition_by_query_pos( max(_limit_range_node1) ) ]
                        endQ, endS = epos.query_pos+1, epos.sbjct_pos+1

                if _limit_range_node2:
                    if startS < min(_limit_range_node2):
                        spos = pacbpobj._positions[ pacbpobj.alignmentposition_by_sbjct_pos( min(_limit_range_node2) ) ]
                        startQ, startS = spos.query_pos, spos.sbjct_pos
                    if endS-1 > max(_limit_range_node2):
                        epos = pacbpobj._positions[ pacbpobj.alignmentposition_by_sbjct_pos( max(_limit_range_node2) ) ]
                        endQ, endS = epos.query_pos+1, epos.sbjct_pos+1

            else:
                # A Pacbp or PacbpDNA object.
                # Limiting ranges by preset ranges (_limit_range_node1 and node2)
                # is less easy; check the sequence for occurring gaps etc.
                # For now, leave it here as such. This means that the function
                # overall_minimal_spanning_range() has a different result for
                # PacbP objects than for PacbPORF, PacbPDNA objects
                pass


            # and update the ranges; this can only enlarge the final range
            thisrangeQ.update( range(startQ, endQ) )
            thisrangeS.update( range(startS, endS) )

        # and return the ranges
        return thisrangeQ, thisrangeS

    # end of function _get_ranges_by_nodes


    def minimal_spanning_range(self,organism=None,node=None):
        """
        TODO description of minimal_spanning_range

        @attention: only functional if pacbp's are present in the graph
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
        if self._msr:
            if organism:    return self._msr[self.node_by_organism(organism)]
            elif node:      return self._msr[node]
            else:           return self._msr

        initial_ranges = {}
        # first iteration: get minimal spanning range by mutual comparison
        for (node1,node2) in self.pairwisecrosscombinations_node():
            # check if these nodes are indeed present as an edge
            if not self.has_edge(node1,node2): continue
            # get ranges for these nodes
            thisrangeQ, thisrangeS = self._get_ranges_by_nodes(node1,node2)
            # set to ranges
            if not initial_ranges.has_key(node1):
                initial_ranges[node1] = thisrangeQ
            else:
                initial_ranges[node1].intersection_update(thisrangeQ)
            if not initial_ranges.has_key(node2):
                initial_ranges[node2] = thisrangeS
            else:
                initial_ranges[node2].intersection_update(thisrangeS)

        # return ranges or the range of a specific organism or node
        if organism:
            return initial_ranges[self.node_by_organism(organism)]
        elif node:
            return initial_ranges[node]
        else:
            return initial_ranges

    # end of function minimal_spanning_range


    def has_minimal_spanning_range(self):
        """
        Has this (nearly Fully Connected) Graph a minimal spanning range?

        @attention: only functional if pacbp's are present in the graph
                    (perform the function `harvest_pacbps_from_crossdata`)

        @attention: only meaningfull for (nearly) Fully Connected Graphs

		@rtype:  Boolean
		@return: True or False
        """
        if self.minimal_spanning_range_sizes().values():
            if min(self.minimal_spanning_range_sizes().values()) == 0:
                return False
            else:
                return True
        else:
            return False

    # end of function has_minimal_spanning_range


    def has_msr(self):
        """
        @attention: alias function name of has_minimal_spanning_range
        """
        return self.has_minimal_spanning_range()

    # end of function has_msr


    def minimal_spanning_range_sizes(self):
        """
        Get the sizes of the minimal spanning range of all pacbp's for each node

        @attention: only functional if pacbp's are present in the graph
                    (perform the function `harvest_pacbps_from_crossdata`)

        @attention: only meaningfull for (nearly) Fully Connected Graphs

		@rtype:  dictionary
		@return: dictionary with nodes (keys) and size of spanning range Sets (values)

        """
        return _size_dict(self.minimal_spanning_range())

    # end of function minimal_spanning_range_sizes


    def minsr_starts(self):
        """
        Return a dictionary of MSR / MINSR start coodinates (keys==nodes)

        @rtype:  dictionary 
        @return: dict[Node] = MINSR start coords (integer) as values
        """
        return _min_dict(self.minimal_spanning_range())

    # end of function minsr_starts


    def minsr_ends(self):
        """
        Return a dictionary of MSR / MINSR end coodinates (keys==nodes)

        @rtype:  dictionary 
        @return: dict[Node] = MINSR end coords (integer) as values
        """
        return _max_dict(self.minimal_spanning_range())

    # end of function minsr_ends


    def minsr_distance_between_codingblocks(self,other,**kwargs):
        """
        Distance in AA between the MSR/MINSR of two CodingBlockGraphs

        @attention: see MINSR._calculate_cbg_distance() for **kwargs
        """
        coords1 = self.minimal_spanning_range()
        coords2 = other.minimal_spanning_range()
        return _calculate_cbg_distance(self,other,coords1,coords2,**kwargs)

    # end of function minsr_distance_between_codingblocks


    def msr_distance_between_codingblocks(self,other,**kwargs):
        """
        Distance in AA between the MSR/MINSR of two CodingBlockGraphs

        @attention: alias of minsr_distance_between_codingblocks
        """
        return self.minsr_distance_between_codingblocks(other,**kwargs)

    # end of function msr_distance_between_codingblocks

# end of class MinimalSpanningRange

################################################################################
# Helper functions for coordinate recalculations etc.
################################################################################

def _size_dict(rangedict):
    """ """
    return dict([(node,len(vlist)) for node, vlist in rangedict.iteritems() ])

# end of function _size_dict


def _min_dict(rangedict):
    """ """
    return dict([(node,min(vlist)) for node, vlist in rangedict.iteritems() ])

# end of function _min_dict


def _max_dict(rangedict):
    """ """
    return dict([(node,max(vlist)) for node, vlist in rangedict.iteritems() ])

# end of function _max_dict



def _calculate_cbg_distance(cbg1,cbg2,coords1,coords2,organism=None,node=None):
    """ "" 

    Distance in AA between two CodingBlockGraphs

    @type  cbg1: CodingBlockGraph (or related) object
    @param cbg1: CodingBlockGraph (or related) object

    @type  cbg2: CodingBlockGraph (or related) object
    @param cbg2: CodingBlockGraph (or related) object

    @type  coords1: dict
    @param coords1: MSR/MINSR/MAXSR/OMSR coordinate range dictionary

    @type  coords2: dict
    @param coords2: MSR/MINSR/MAXSR/OMSR coordinate range dictionary

    @type  organism: *
    @param organism: Organism identifier (or None)

    @type  node: *
    @param node: Nodeidentifier (or None)


    @rtype:  dictionary (or integer if Organism or Node is specified)
    @return: Organisms (keys) and AA-distance between CBGs (values),
             or only a distance if Organism or Node identifier was specified

    @attention: cbg2 is supposed to be 3p/rigth of cbg1
    """

    distances = {}
    for org in cbg1.organism_set().intersection(cbg2.organism_set()):
        if organism and org != organism: continue
        nodeA = cbg1.node_by_organism(org)
        nodeB = cbg2.node_by_organism(org)
        if node and node != nodeA: continue
        if node and node != nodeB: continue
        distA = min( coords1[nodeA] ) - max( coords2[nodeB] )
        distB = min( coords2[nodeB] ) - max( coords1[nodeA] )
        distances[org] = max( [ distA, distB ] ) - 1

    # return distance only of a specific organism is requested for
    if organism:
        return distances[organism]
    else:
        return distances

# end of function _calculate_cbg_distance

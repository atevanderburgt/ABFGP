"""
Pacbp accessibility functions joined in a class that is inherited in
other graph classes for Aligment Based Gene Predictions
"""

# Pacb imports
import pacb

# graphAbgp imports
import ordering
from exceptions import *

# Abfgp imports
from lib_stopwatch import StopWatch

# Python Imports

class PacbpAccesibilityFunctions:
    """
    Class with PacbP(ORF) accesibility functions

    @attention: class cannot be used stand-alone!
    @attention: used in CodingBlockGraph class by inheritance
    @attention: used in PacbpCollectionGraph class by inheritance
    """
    def get_pacbps_by_nodes(self,node1=None,node2=None,order_by='bits'):
        """
        Get the pacbp(s) from the CodingBlockGraph by node, optional two nodes

        @type  node1: *
        @param node1: node identifier

        @type  node2: *
        @param node2: node identifier (or None)

        @type  order_by: string
        @param order_by: 'length' (DESC), 'bits' (DESC), 'start' (ASC);
			 default 'bits'

        @rtype:  list
        @return: list of pacbps objects

        @attention: when only a single node is requested, pacbp(s) are swapped
		    to make the requested node the query node
        """
        if node1 not in self.get_nodes():
            message = "node1 `%s` not in graph: %s" % ( node1,self.get_nodes())
            raise InproperlyAppliedArgument, message
        if node2 and node2 not in self.get_nodes():
            message = "node2 `%s` not in graph: %s" % ( node2,self.get_nodes())
            raise InproperlyAppliedArgument, message
        if order_by not in ['bits','length']:
            order_by ='bits'

        # if no pacbps are stored into the object yet, return []
        if not self.pacbps: return []

        thepacbps = []        
        for (key,nodeA,nodeB),pacbporf in self.pacbps.iteritems():
            if nodeA == node1 or nodeB == node1:
                if not node2:
                    if nodeA == node1:
                        thepacbps.append( pacbporf )
                    else:
                        # swap query and sbjct!
                        thispacbporf = pacb.swap_query_and_sbjct(pacbporf)
                        thepacbps.append( thispacbporf )
                else:
                    if nodeA == node2 or nodeB == node2:
                        thepacbps.append( pacbporf )
                    else:
                        pass

        # order the requested pacbps
        if order_by == 'bits':
            thepacbps = ordering.order_list_by_attribute(
		thepacbps,"bits",reversed=True)
        else:
            thepacbps = ordering.order_list_by_attribute(
		thepacbps,"length",reversed=True)

        # return the requested pacbps
        return thepacbps

    # end of function get_pacbps_by_nodes


    def get_pacbp_by_organisms(self,orgA,orgB):
        """
        @type  orgA: * (string)
        @param orgA: Organism identifier

        @type  orgB: * (string)
        @param orgB: Organism identifier

        @rtype:  PacbP(ORF) object 
        @return: PacbP(ORF) object or None

        @attention: requires presence of a single Orf node per
		    Organism node in the graph ( K(s) graph, e.g. CBG )!
        """
        # check if requested organisms are present in this graph
        if orgA not in self.organism_set() or orgB not in self.organism_set():
            raise OrganismNotPresentInGraph

        # if no pacbps are stored into the object yet, return None 
        if not self.pacbps: return None 

        thepacbp = None
        for (key,(org1,orf1),(org2,orf2)),pacbporf in self.pacbps.iteritems():
            if orgA == org1 and orgB == org2:
                thepacbp = pacbporf
                break
            if orgB == org1 and orgA == org2:
                thepacbp = pacbporf
                break
        # return the requested pacbp
        return thepacbp

    # end of function get_pacbp_by_organisms


    def get_pacbps_by_organisms(self,orgA,orgB,order_by='bits'):
        """
        @type  orgA: * (string)
        @param orgA: Organism identifier

        @type  orgB: * (string)
        @param orgB: Organism identifier

        @type  order_by: string
        @param order_by: 'length' (DESC), 'bits' (DESC), 'start' (ASC);
                         default 'bits'

        @rtype:  list
        @return: list of pacbps objects
        """
        # check if requested organisms are present in this graph
        if orgA not in self.organism_set():
            raise OrganismNotPresentInGraph, orgA
        if orgB not in self.organism_set(): 
            raise OrganismNotPresentInGraph, orgB

        # if no pacbps are stored into the object yet, return empty list 
        if not self.pacbps: return [] 

        # if no pacbps are stored into the object yet, return []
        if not self.pacbps: return []

        thepacbps = []
        for (key,(org1,orf1),(org2,orf2)),pacbporf in self.pacbps.iteritems():
            if orgA == org1 and orgB == org2:
                thepacbps.append( pacbporf )
            elif orgB == org1 and orgA == org2:
                thepacbps.append( pacbporf )
            else:
                continue

        # order the requested pacbps
        if order_by == 'bits':
            thepacbps = ordering.order_list_by_attribute(
                thepacbps,"bits",reversed=True)
        else:
            thepacbps = ordering.order_list_by_attribute(
                thepacbps,"length",reversed=True)

        # return the requested pacbps
        return thepacbps

    # end of function get_pacbps_by_organisms


    def get_pacbps_by_organism(self,organism,order_by=None):
        """
        Get the pacbp(s) from the CodingBlockGraph of a single organism

        @type  organism: * (string)
        @param organism: Organism identifier

        @type  order_by: string
        @param order_by: 'length' (DESC), 'bits' (DESC), 'node' or None
			 (on node); default None

	@rtype:  list
	@return: list of pacbps objects

        @attention: pacbps are swapped such that `organism` is always the query!
        @attention: pacbps are ordered by their sbjct nodes
        """
        # check if requested organism is present in this graph
        if organism not in self.organism_set():
            raise OrganismNotPresentInGraph

        # if no pacbps are stored into the object yet, return []
        if not self.pacbps: return []

        # reset order_by if falsely assigned
        if order_by not in [None,'bits','length']:
            order_by = None

        thepacbps = []
        for (key,(org1,orf1),(org2,orf2)),pacbporf in self.pacbps.iteritems():
            if organism ==  org1:
                thepacbps.append( ( (org2,orf2), pacbporf ) )
            elif organism ==  org2:
                # swap query and sbjct!
                thispacbporf = pacb.swap_query_and_sbjct(pacbporf)
                thepacbps.append( ( (org1,orf1), thispacbporf ) )
            else:
                pass
        # sort the requested pacbps on Node
        thepacbps.sort()
        thepacbps = [ pacbporf for node,pacbporf in thepacbps ]

        # order the requested pacbps if requested for
        if order_by == 'bits':
            thepacbps = ordering.order_list_by_attribute(
		thepacbps,"bits",reversed=True)
        if order_by == 'length':
            thepacbps = ordering.order_list_by_attribute(
		thepacbps,"length",reversed=True)

        # return the requested pacbps
        return thepacbps

    # end of function get_pacbps_by_organism


    def remove_pacbp(self,pacbp,nodeQ,nodeS):
        """
        Remove a Pacbp(DNA/ORF) object from this graph

        @type  pacbp: PacbP, PacbPDNA or PacbPORF object
        @param pacbp: object instance

        @type  nodeQ: *
        @param nodeQ: PacbP Query Node identifier

        @type  nodeS: *
        @param nodeS: PacbP Sbjct Node identifier

        @attention: USE WITH CARE; only use when you know what you are intending
        """
        searchkey = pacbp.construct_unique_key(nodeQ,nodeS)
        keys = self.pacbps.keys()
        is_deleted = False
        for (pacbptuplekey,n1,n2) in keys:
            if (n1,n2) == (nodeQ,nodeS):
		thispacbp = self.pacbps[(pacbptuplekey,n1,n2)]
                thiskey   = thispacbp.construct_unique_key(nodeQ,nodeS)
                if thiskey == searchkey:
                    del( self.pacbps[(pacbptuplekey,n1,n2)] )
                    is_deleted = True
                    break
        # return weather or not pacbp is deleted
        return is_deleted

    # end of function remove_pacbp


    def extend_pacbporfs(self,input):
        """
        Convert PacbP -> PacbPORF objcects and extend them

        @type  input: dict
        @param input: input <data structure> that contains lists of Orfs
        """
        try:
            # pacbps are already (unextended) PacPORFs
            for (k,n1,n2), pacbporf in self.pacbps.iteritems():
                pacbporf.extend_pacbporf_after_stops()
                self.pacbps[(k,n1,n2)] = pacbporf
        except:
            # no, pacbps are still PacbPs
            self.pacbps2pacbporfs(input)
            for (k,n1,n2), pacbporf in self.pacbps.iteritems():
                pacbporf.extend_pacbporf_after_stops()
                self.pacbps[(k,n1,n2)] = pacbporf 

    # end of function extend_pacbporfs


    def pacbps2pacbporfs(self,input):
        """
        Convert PacbP -> PacbPORF objects, but DO NOT extend them

        @type  input: dict
        @param input: input <data structure> that contains lists of Orfs
        """
        for (k,n1,n2), pacbp in self.pacbps.iteritems():
            orgQ = self.organism_by_node(n1)
            orfQ = input[orgQ]['orfs'].get_orf_by_id( n1[1] )
            orgS = self.organism_by_node(n2)
            orfS = input[orgS]['orfs'].get_orf_by_id( n2[1] )
            pacbporf = pacb.conversion.pacbp2pacbporf(pacbp,orfQ,orfS)
            self.pacbps[(k,n1,n2)] = pacbporf 

    # end of function pacbps2pacbporfs


    def pacbporfs2pacbps(self):
        """
        Convert PacbPORF -> PacbP objects, a BACKWARDS conversion
	(e.g. for splitting == quicker!)
        """
        for (k,n1,n2), pacbporf in self.pacbps.iteritems():
            pacbp = pacb.conversion.pacbporf2pacbp(pacbporf)
            self.pacbps[(k,n1,n2)] = pacbp

    # end of function pacbporfs2pacbps


    def _recrute_pacbporfs_from_parental_cbg(self,parentcbg,
        create_cache=True,
        ignore_nonexisting_edges=False,verbose=False):
        """
	Harvest PacbPORFs from (parental) CodingBlockGraph

	@attention: alternative for harvest_pacbps_from_crossdata()
        @attention: required in _place_cbg_in_partialgsg() function
        @attention: use create_cache=False with care!

	@type  parentcbg: CodingBlockGraph
	@param parentcbg: CodingBlockGraph that has to delived PacbPORFs

	@type  create_cache: Boolean
	@param create_cache: run the create_cache() function on the CBG (self)

	@type  ignore_nonexisting_edges: Boolean
	@param ignore_nonexisting_edges: when False, do not create edges in the
					 CBG (self) that are absent (but present
					 in the parentcbg)

        @type  verbose: Boolean
        @param verbose: print debugging information to STDOUT when True
        """
        replacements = {}
        substituted = 0

        ####################################################################
        if verbose:
            stw = StopWatch("recruteParentalPacbps")
            print stw.start()
            print "target:", self
            print "source:", parentcbg
        #################################################################### 

        for (node1,node2) in self.pairwisecrosscombinations_node():
            # if this edge is not present in the parent, ignore it
            if not parentcbg.has_edge(node1,node2): continue
            # get PacbPORF of the parent
            origpacbporf = parentcbg.get_pacbps_by_nodes(
		    node1=node1,node2=node2)[0]
            curpacbporf = None
            replace_pacbporf = False
            if not self.has_edge(node1,node2):
                if ignore_nonexisting_edges:
                    # if ignore_nonexisting_edges -> do not recrute this pacbp
                    continue 
                else:
                    # replace this Pacbporf if it exists and
		    # simultaniously create novel edge
                    replace_pacbporf = True
            elif self.has_edge(node1,node2) and not\
	    self.get_pacbps_by_nodes(node1=node1,node2=node2):
                replace_pacbporf = True
            else:
                curpacbporf  = self.get_pacbps_by_nodes(
			node1=node1,node2=node2)[0]
                if pacb.comparison.IsIdenticalPacbPORF(
		    origpacbporf,curpacbporf):
                    # Pacbporfs are already identical; not relevant to copy
                    continue
                if origpacbporf.issuperset(curpacbporf):
                    # store to replacements dict
                    replacements[(node1,node2)] = curpacbporf
                    # remove from the CBG -> replacement in progress
                    self.remove_pacbp(curpacbporf,node1,node2)
                    replace_pacbporf = True

            # check if replace_pacbporf is set to True
            if replace_pacbporf:
                ################################################################
                if verbose:
                    print stw.lap(), "REPLACING PacbPORF Source->Target:"
                    print "T:", curpacbporf, "(current)"
                    print "S:", origpacbporf
                    origpacbporf.print_protein(_linesize=100) 
                ################################################################
                newkey = origpacbporf.construct_unique_key(node1,node2)
                self.set_edge_weight( node1, node2, wt=origpacbporf.bitscore )
                self.pacbps[(newkey,node1,node2)] = origpacbporf
                substituted+=1

        # check if substitutions have been taken place
        if create_cache and substituted:
            #####################################################################
            if verbose:
                print stw.lap(), "CREATE_CACHE & substituted PacbPORFS:",
                print substituted, "edges:", len(self.weights)/2,
                print "pacbps:", len(self.pacbps)
                ####for k,pacbporf in self.pacbps.iteritems():
                ####    print k,"\n",pacbporf
            #####################################################################
            self.clear_cache()
            # check if there is an OMSR upon recreation; in very
            # exceptional cases, OMSR can get lost in this step
            if self.has_overall_minimal_spanning_range():
                self.create_cache()
                self.update_edge_weights_by_minimal_spanning_range()
            else:
                #############################################################
                if verbose:
                    print stw.lap(), "OMSR got lost!",
                    print "replacements:",  len(replacements)
                    for (n1,n2), curpacbporf in replacements.iteritems():
                        print "REP:", curpacbporf, n1, n2
                #############################################################
                # OMSR got lost! Restore replacements dict and as such
                # restore the original PacbPs one by one (in random order)
                # and quit as soon as an OMSR is restored
                for (node1,node2),curpacbporf in replacements.iteritems():
                    newkey = curpacbporf.construct_unique_key(node1,node2)
                    tobereplpacbporf = self.get_pacbps_by_nodes(node1=node1,node2=node2)[0]
                    # remove from the CBG
                    self.remove_pacbp(tobereplpacbporf,node1,node2)
                    # and place back the original one
                    self.set_edge_weight( node1, node2, wt=curpacbporf.bitscore )
                    self.pacbps[(newkey,node1,node2)] = curpacbporf
                    substituted-=1
                    if self.has_overall_minimal_spanning_range():
                        self.create_cache()
                        self.update_edge_weights_by_minimal_spanning_range()
                        #########################################################
                        if verbose:
                            print stw.lap(), "OMSR restored, substitutions:",
                            print substituted 
                            print "T:", self
                        ##########################################################
                        # break out of the for loop of PacbP replacement
                        break

        # return number of replaced/added pacbporfs
        return substituted

    # end of function _recrute_pacbporfs_from_parental_cbg

# end of class PacbpAccesFunctions

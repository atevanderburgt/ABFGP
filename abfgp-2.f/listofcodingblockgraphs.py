"""
Class that stores an (unordered) list of CodingBlockGraphs on which several
operations can be performed on. Used in Alignment Based Gene Prediction.
"""

# Module metadata
__authors__ = "Ate van der Burgt"
__license__ = "MIT"

# Function Imports
from graphAbgp.ordering import (
    order_list_by_attribute,
    order_graphlist_by_total_weight,
    order_graphlist_by_total_weight_and_identity,
    )

# Exception Imports
from abgp_exceptions import NoCrossdataApplied, NoInputApplied

class ListOfCodingBlockGraphs:
    """ """
    def __init__(self,listofcbgs,input={},crossdata={}):
        """
        @type  listofcbgs: list
        @param listofcbgs: list of CodingBlockGraphs (and potential other types)

        @type  input: dict
        @param input: input <dict data structure>

        @type  crossdata: dict
        @param crossdata: crossdata <dict data structure>
        """
        self.codingblockgraphs = []
        self.codingblockgraphs.extend( listofcbgs )
        self._input     = input
        self._crossdata = crossdata

    # end of function __init__

    ############################################################################
    ### Build-in functions to enshure behavious of LOCBGS asa list of CBGs  ####
    ############################################################################

    def __len__(self):
        """ Return the length of list self.codingblockgraphs """
        return len(self.codingblockgraphs)
    # end of function __len__


    def __str__(self):
        """ Nicely formatted onleliner of this ListOfCodingBlockGraphs """
        return "<%s (%s) N,E,P: %s,%s,%s i,c: %s,%s >" % (
            self.__class__.__name__, len(self),
            sum([cbg.node_count() for cbg in self]),
            sum([cbg.edge_count() for cbg in self]),
            sum([len(cbg.pacbps) for cbg in self]),
            str(len(self._input) > 0)[0],
            str(len(self._crossdata) > 0)[0], 
            )
    # end of function __str__

    def __iter__(self):
        """ Return an iterable of self.codingblockgraphs """
        return iter(self.codingblockgraphs)
    # end of function __iter__

    def __getitem__(self,x):
        """ Return an xth element of self.codingblockgraphs """
        return self.codingblockgraphs[x]
    # end of function __getitem__

    #def __setitem__(self,x,item):
    #    """ Add xth element of self.codingblockgraphs """
    #    return self.codingblockgraphs.insert(x,item)
    ## end of function __setitem__

    #def __delitem__(self,x):
    #    """ Remove xth element of self.codingblockgraphs """
    #    return self.codingblockgraphs.pop(x)
    ## end of function __delitem__

    ############################################################################
    ### functions that do not change the number of CBGs in the list         ####
    ############################################################################

    def harvest_pacbps_from_crossdata(self,crossdata={}):
        """
        @attention: see CodingBlockGraph.harvest_pacbps_from_crossdata()
        """
        if not crossdata: crossdata = self._crossdata
        if not crossdata: raise NoCrossdataApplied
        for cbg in self: cbg.harvest_pacbps_from_crossdata(crossdata)

    # end of function harvest_pacbps_from_crossdata


    def harvest_pacbps_from_pacbpcollection(self,pcg):
        """
        Harvest PacbP(ORF) objects from a PacbpCollectionGraph

        @type  pcg: PacbpCollectionGraph
        @param pcg: PacbpCollectionGraph
        
        @attention: see CodingBlockGraph.harvest_pacbps_from_pacbpcollection()
        """
        for cbg in self: cbg.harvest_pacbps_from_pacbpcollection(pcg)

    # end of function harvest_pacbps_from_pacbpcollection


    def update_cbgs_by_removing_alternatives_from_pacbps_dict(self):
        """
        @attention: requires presence of PacbP(ORF) objects
        @attention: see CodingBlockGraph.update_cbgs_by_removing_alternatives_
                    from_pacbps_dict()
        """
        for cbg in self: cbg.remove_alternatives_from_pacbps_dict()

    # end of function update_cbgs_by_removing_alternatives_from_pacbps_dict


    def keep_only_best_alternative_from_pacbps_dict(self):
        """
        @attention: requires presence of PacbPORF objects
        @attention: see CodingBlockGraph.keep_only_best_alternative_from_
                    pacbps_dict()
        """
        for cbg in self: cbg.keep_only_best_alternative_from_pacbps_dict()

    # end of function


    def extend_pacbporfs(self,input={}):
        """
        @attention: requires presence of PacbPORF objects
        @attention: see CodingBlockGraph.extend_pacbporfs()
        """
        if not input: input = self._input
        if not input: raise NoInputApplied
        for cbg in self: cbg.extend_pacbporfs(input)

    # end of function extend_pacbporfs


    def make_pacbps_for_missing_edges(self):
        """
        @attention: requires presence of PacbP(ORF) objects
        @attention: see CodingBlockGraph.make_pacbps_for_missing_edges()
        """
        for cbg in self: cbg.make_pacbps_for_missing_edges()

    # end of function make_pacbps_for_missing_edges


    def update_edge_weights_by_minimal_spanning_range(self):
        """
        @attention: requires presence of PacbPORF objects
        @attention: see CodingBlockGraph.update_edge_weights_by_
                    minimal_spanning_range()
        """
        for cbg in self: cbg.update_edge_weights_by_minimal_spanning_range()

    # end of function update_edge_weights_by_minimal_spanning_range


    def try_create_omsr_by_edge_deletion(self):
        """
        @attention: see CodingBlockGraph.try_create_omsr_by_edge_deletion()
        """
        for cbg in self: cbg.try_create_omsr_by_edge_deletion()

    # end of function try_create_omsr_by_edge_deletion


    def recover_rejected_pacbps(self,crossdata={}):
        """
        @attention: see CodingBlockGraph.recover_rejected_pacbps()
        """
        if not crossdata: crossdata = self._crossdata
        if not crossdata: raise NoCrossdataApplied
        for cbg in self: cbg.recover_rejected_pacbps(crossdata)

    # end of function recover_rejected_pacbps 


    def order_list_by_attribute(self,order_by='',reversed=False):
        """
        @attention: see graphAbgp.ordering.order_list_by_attribute
        """
        self.codingblockgraphs = order_list_by_attribute(
            self.codingblockgraphs,order_by=order_by,reversed=reversed)

    # end of function order_list_by_attribute


    def order_graphlist_by_total_weight(self):
        """
        @attention: see graphAbgp.ordering.order_graphlist_by_total_weight
        """
        self.codingblockgraphs = order_graphlist_by_total_weight(
            self.codingblockgraphs)

    # end of function order_graphlist_by_total_weight 


    def order_graphlist_by_total_weight_and_identity(self):
        """
        @attention: see graphAbgp.ordering.order_graphlist_by_total_weight
        """
        self.codingblockgraphs = order_graphlist_by_total_weight_and_identity(
            self.codingblockgraphs)

    # end of function order_graphlist_by_total_weight_and_identity

    ############################################################################
    ### functions that can change the number of CBGs in the list            ####
    ############################################################################

    def remove_all_but_cbgs(self):
        """
        Remove all non-CodingBlockGraph objects from the list

        @rtype:  list
        @return: list with removed non-CodingBlockGraph objects
        """
        deleted = []
        for i in range(len(self)-1,-1,-1):
            try:
                if self[i].__class__.__name__ != 'CodingBlockGraph':
                    deleted.append( self.codingblockgraphs.pop(i) )
                elif self[i].node_count() <= 2:
                    # a CBG must have at least 3 nodes...
                    deleted.append( self.codingblockgraphs.pop(i) )
                else:
                    pass
            except:
                # element in list without __class__.__name__ attribute
                deleted.append( self.codingblockgraphs.pop(i) )
        # return list of deleted elements
        return deleted

    # end of function remove_all_but_cbgs


    def remove_all_but_complete_cbgs(self):
        """
        Remove all incomplete CodingBlockGraphs

        @attention: function removes as well all non-CodingBlockGraph objects

        @rtype:  list
        @return: list with removed incomplete CodingBlockGraph objects
        """
        # verify that all objects in list are CBGs
        deleted = self.remove_all_but_cbgs()
        for i in range(len(self)-1,-1,-1):
            if self[i].connectivitysaturation() < 1.0:
                deleted.append( self.codingblockgraphs.pop(i) )
        # return list of deleted elements
        return deleted

    # end of function remove_all_but_complete_cbgs

        
    def remove_cbgs_without_msr(self):
        """
        Remove all CodingBlockGraphs with no MSR

        @attention: function removes as well all non-CodingBlockGraph objects

        @rtype:  list
        @return: list with removed CodingBlockGraph objects without MSR
        """
        # verify that all objects in list are CBGs
        deleted = self.remove_all_but_cbgs()
        for i in range(len(self)-1,-1,-1):
            if not self[i].has_minimal_spanning_range():
                deleted.append( self.codingblockgraphs.pop(i) )
        # return list of deleted elements
        return deleted

    # end of function remove_cbgs_without_msr


    def remove_cbgs_without_omsr(self):
        """
        Remove all CodingBlockGraphs with no OMSR

        @attention: function removes as well all non-CodingBlockGraph objects

        @rtype:  list
        @return: list with removed CodingBlockGraph objects without OMSR
        """
        # verify that all objects in list are CBGs
        deleted = self.remove_all_but_cbgs()
        for i in range(len(self)-1,-1,-1):
            if not self[i].has_overall_minimal_spanning_range():
                deleted.append( self.codingblockgraphs.pop(i) )
        # return list of deleted elements
        return deleted

    # end of function remove_cbgs_without_omsr


    def remove_cbgs_with_lt_nodes(self,x):
        """
        Remove all CodingBlockGraphs with less than x nodes

        @type  x: integer
        @param x: minimal number of nodes in the CBG in order to be retained

        @attention: function removes as well all non-CodingBlockGraph objects

        @rtype:  list
        @return: list with removed CodingBlockGraph objects with to few nodes
        """
        # border values for minimal node count in order to stay a CBG
        x = max([x,3])
        # verify that all objects in list are CBGs
        deleted = self.remove_all_but_cbgs()
        for i in range(len(self)-1,-1,-1):
            if self[i].node_count() < x:
                deleted.append( self.codingblockgraphs.pop(i) )
        # return list of deleted elements
        return deleted

    # end of function remove_cbgs_with_lt_nodes


    def remove_cbgs_with_lt_edges(self,x):
        """
        Remove all CodingBlockGraphs with less than x edges

        @type  x: integer
        @param x: minimal number of edges in the CBG in order to be retained

        @attention: function removes as well all non-CodingBlockGraph objects

        @rtype:  list
        @return: list with removed CodingBlockGraph objects with to few edges
        """
        # border values for minimal edge count in order to stay a CBG
        x = max([x,2])
        # verify that all objects in list are CBGs
        deleted = self.remove_all_but_cbgs()
        for i in range(len(self)-1,-1,-1):
            if self[i].edge_count() < x:
                deleted.append( self.codingblockgraphs.pop(i) )
        # return list of deleted elements
        return deleted

    # end of function remove_cbgs_with_lt_edges


    def remove_incompatible_cbgs(self,
        minimal_node_count=None,
        minimal_edge_count=None,
        filter_for_msr=True,
        filter_for_omsr=False):
        """
        Remove CodingBlockGraphs according to different criteria

        @type  minimal_node_count: integer (or None)
        @param minimal_node_count: minimal number of nodes in the CBG
                                   in order to be retained

        @type  minimal_edge_count: integer (or None)
        @param minimal_edge_count: minimal number of edges in the CBG
                                   in order to be retained

        @type  filter_for_msr: Boolean
        @param filter_for_msr: when True, only CBGs with MSR will pass

        @type  filter_for_omsr: Boolean
        @param filter_for_omsr: when True, only CBGs with OMSR will pass

        @attention: function removes as well all non-CodingBlockGraph objects

        @rtype:  list
        @return: list with removed CodingBlockGraph objects
        """
        # border values for minimal node & edge count in order to stay a CBG
        if minimal_node_count:
            minimal_node_count = max([minimal_node_count,3])
        if minimal_edge_count:
            minimal_edge_count = max([minimal_edge_count,2])
        
        # verify that all objects in list are CBGs
        deleted = self.remove_all_but_cbgs()
        if minimal_node_count:
            deleted.extend( self.remove_cbgs_with_lt_nodes(minimal_node_count) )
        else:
            if self._input:
                deleted.extend(
                    self.remove_cbgs_with_lt_nodes(len(self._input))
                    )
        if minimal_edge_count:
            deleted.extend( self.remove_cbgs_with_lt_edges(minimal_edge_count) )
        else:
            if self._crossdata:
                deleted.extend(
                    self.remove_cbgs_with_lt_edges(len(self._crossdata))
                    )


        if filter_for_msr:
            deleted.extend( self.remove_cbgs_without_msr() )
        if filter_for_omsr:
            deleted.extend( self.remove_cbgs_without_omsr() )
        # return list of deleted elements
        return deleted

    # end of function remove_incompatible_cbgs


    def remove_reversecomplement_cbgs(self,overlap_ratio = None):
        """
        Remove all CBGs that are ReversecomplementCodingBlockGraph

        @attention: function removes as well all non-CodingBlockGraph objects

        @rtype:  ( list, list )
        @return: [ deleted CBGs ], [ ReversecomplementCodingBlockGraphs ]
        """
        # verify that all objects in list are CBGs
        deleted = self.remove_all_but_cbgs()
        reversedcbgs = []
        for i in range(len(self)-1,-1,-1):
            cbg = self[i]
            IS_DELETED = False
            for earlier_obtained_revcbg in reversedcbgs:
                if cbg.is_reversecomplement(revcbg=earlier_obtained_revcbg):
                    deleted.append( self.codingblockgraphs.pop(i) )
                    IS_DELETED = True
                    break
            # continue if cbg is deleted
            if IS_DELETED: continue
            # obtain reverse complement CBG
            revcbg = cbg.get_reversecomplement()
            if cbg.is_reversecomplement(revcbg=revcbg):
                reversedcbgs.append( revcbg )
                deleted.append( self.codingblockgraphs.pop(i) )
                
        # return list with (unique) reversedcbgs & list of deleted CBGs
        return deleted, reversedcbgs

    # end of function remove_reversecomplement_cbgs


    def split_codingblock_on_alternatives_in_pacbps_dict(self,
        filter_for_msr=True,
        filter_for_omsr=False):
        """
        """
        for i in range(len(self)-1,-1,-1):
            splits = self[i].split_codingblock_on_alternatives_in_pacbps_dict(
                    filter_for_msr=filter_for_msr,
                    filter_for_omsr=filter_for_omsr
                    )

            if len(splits) > 1:
                # update current CBG, append new CBG(s)
                self.codingblockgraphs[i] = splits[0]
                for newcbg in splits[1:]: self.codingblockgraphs.append(newcbg)
            elif len(splits) == 1:
                # update current CBG
                self.codingblockgraphs[i] = splits[0]
            else:
                # no CBGs in splits !? Should be impossible, but just for shure
                self.codingblockgraphs.pop(i)

    # end of function split_codingblock_on_alternatives_in_pacbps_dict


    def try_create_omsr_by_edge_deletion(self):
        """
        TODO: move to CodingBlockGraph
        """
        if self.has_overall_minimal_spanning_range():
            return None

        # clear_cache() but do not create_cache() !
        self.clear_cache()

        for nodeQ,nodeS in self.pairwisecrosscombinations_node():
            # (temporarily) remove this edge
            wt = self.get_edge_weight(nodeQ,nodeS)
            sg.del_edge(nodeQ,nodeS)

            # check if this CBG has a minimal spanning range now
            if self.has_overall_minimal_spanning_range():
                # Yes, now it is a CBG with a OMSR. Remove the PacbP(s)
                # that belonged to the removed edge
                self._update_pacbps()
                # break the loop and exit the function with a True
                return True
            else:
                # nope, still not an overall minimal spannign range
                # set back the edge and try another edge....
                sg.add_edge(nodeQ,nodeS,wt=wt)
                continue
            
        # if this point is reached, edge deletion did not yield an OMSR
        return False

    # end of function try_create_omsr_by_edge_deletion

# end of class ListOfCodingBlockGraphs

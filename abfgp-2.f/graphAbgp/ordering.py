"""
Ordering lists of abgp-type of graphs in Aligment Based Gene Prediction
"""

def order_graphlist_by_total_weight(gralist):
    """ 
    """ 
    # sort the subgraphs based on their total_weight
    tmp = []
    for sg in gralist:
        tmp.append( ( sg.total_weight(), sg ) )
    tmp.sort()
    tmp.reverse()
    return [ sg for (tw,sg) in tmp ]
        
# end of function order_graphlist_by_total_weight


def order_graphlist_by_total_weight_and_identity(cbgralist):
    """
    """
    # sort the subgraphs based on their total_weight & identity score
    tmp = []
    for cbg in cbgralist:
        tmp.append( ( cbg.total_weight()*cbg.omsr_identityscore(), cbg ) )
    tmp.sort()
    tmp.reverse()
    return [ cbg for (tw,cbg) in tmp ]

# end of function order_graphlist_by_total_weight_and_identity 


def order_cbgraphlist(cbgralist):
    """
    """
    # sort the subgraphs based on their total_weight
    tmp = []
    for cbg in cbgralist:
        ## CHANGED 21/11/2009 cbg.total_weight() * cbg.omsr_identityscore() !!!
        tmp.append( ( cbg.total_weight()*cbg.omsr_identityscore(), cbg ) )
    tmp.sort()
    tmp.reverse()
    return [ cbg for (tw,cbg) in tmp ]

# end of function order_cbgraphlist 



def order_graphlist_by_identity(gralist):
    """
    """
    # sort the CBGs in the list on identity 
    tmp = []
    for sg in gralist:
        tmp.append( ( sg.omsr_identityscore(), sg ) )
    tmp.sort()
    tmp.reverse()
    return [ sg for (tw,sg) in tmp ]

# end of function order_graphlist_by_identity



def reorder_cbgs_on_node_occurrence(cbgs,prev=None,next=None):
    """
    """
    ordering = []
    for i in range(0,len(cbgs)):
        cbg = cbgs[i]
        cnt = 0
        if prev: cnt+=len(prev.node_set().intersection(cbg.get_nodes()))
        if next: cnt+=len(next.node_set().intersection(cbg.get_nodes()))
        ordering.append( ( cnt, -i, cbg ) )
    # now sort and return
    ordering.sort()
    ordering.reverse()
    return [ cbg for (cnt,i,cbg) in ordering ]

# end of function reorder_cbgs_on_node_occurrence 


def order_list_by_attribute(inputlist,order_by='',reversed=False):
    """
    Helper function to order objects in a list by any given attribute
    
    @type  inputlist: list
    @param inputlist: list of items to be ordered on a attribute
    
    @type  order_by: string
    @param order_by: attribute of items in inputlist to be sorted upon
        
    @type  reversed: Boolean 
    @param reversed: reverse the ordered list (True) or not (False)
            
    @rtype:  list
    @return: ordered list of items, ordered on order_by attribute
    """
    if not inputlist: return inputlist
    if hasattr(inputlist[0],str(order_by)):
        orderedlist = [ ( getattr(o,str(order_by)), o ) for o in inputlist ]
        orderedlist.sort()
        if reversed: orderedlist.reverse()
        return [ o for (attr,o) in orderedlist ]
    else:
        return inputlist

# end of function order_list_by_attribute


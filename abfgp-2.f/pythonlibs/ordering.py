""" Generic functions for ordering lists (of objects) """

# Module metadata
__authors__ = "Ate van der Burgt"
__license__ = "MIT"

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


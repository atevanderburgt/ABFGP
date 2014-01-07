"""
Generic functions for recombining lists
"""
from copy import deepcopy

def pairwise(iterable):
    """
    Make a half cross through all (ordered) pairwise combinations in an iterable

    @type  iterable: *
    @param iterable: python class instance that is iterable (list,tuple,...)

    @rtype:  list
    @return: list with all unique, non-self pairwise ordered combinations
    """
    allcombis = []
    for i in range(0,len(iterable)):
        for j in range(len(iterable)-1,-1,-1):
            if i==j: break
            combi = [iterable[i],iterable[j]]
            combi.sort() 
            allcombis.append(tuple(combi))
    allcombis.sort()
    return allcombis 

# end of function pairwise
    

def cross(iterableofiterables):
    """
    Make a full cross through all iterables in a list of iterables
    
    @type  iterableofiterables: *
    @param iterableofiterables: iterable python class instance with iterables as elements (list,tuple,...)

    @rtype:  list
    @return: list with all unique combinations through all iterables
    """
    allcombis = []
    prev_len = 0
    for items in iterableofiterables:
        if not allcombis:
            for i in range(0,len(items)):
                allcombis.append( [ items[i] ] )
            prev_len = len(allcombis)
            continue
        for duplicate in range(1,len(items)):
            duplpart = deepcopy( allcombis[0:prev_len] )
            allcombis.extend(duplpart)
        for i in range(0,len(items)):
            for combi in allcombis[i*prev_len:(i+1)*prev_len]:
                combi.append(items[i])
        # and update the length of the list for the next iteration
        prev_len = len(allcombis)
    return allcombis

# end of function cross


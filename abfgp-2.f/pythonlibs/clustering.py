""" Generic and simple functions for clustering various types of data """

# Module metadata
__authors__ = "Ate van der Burgt"
__license__ = "MIT"

def cluster_coordinates(coordinatelist,offset):
    """  Cluster a list of (unique) coordinates based on offset """
    if not coordinatelist: return {}
    clustered_coords = []
    coordinatelist = list(set(coordinatelist))
    coordinatelist.sort()
    clustered_coords = [ [ coordinatelist[0] ] ]
    for pos in coordinatelist[1:]:
        if pos - clustered_coords[-1][-1] <= offset:
            clustered_coords[-1].append(pos)
        else:
            clustered_coords.append([pos])
    # return the clustered coordinates
    return clustered_coords
    
# end of function cluster_coordinates


def split_in_continious_ranges(coordinatelist):
    """ Split a list of coordinates into continious ranges of (start,end) coordinates """
    return [ (locus[0],locus[-1]+1) for locus in cluster_coordinates(coordinatelist,1) ]

# end of function split_in_continious_ranges

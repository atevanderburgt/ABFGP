""" Generic statistics & mathematics related functions  """

# Module metadata
__authors__ = "Ate van der Burgt"
__license__ = "MIT"

# Python Imports
from math import sqrt

def meanstdv(x):
    """ Calculate mean and standard deviation of values in an iterable """
    n, mean, std = len(x), 0, 0
    for a in x:
        mean = mean + a
    mean = mean / float(n)
    for a in x:
        std = std + (a - mean)**2
    std = sqrt(std / float(n-1))
    return mean, std
# end of function meanstdv

def product(values):
    """ Calculate the product of values in an iterable """
    retval = values[0]
    for v in values[1:]: retval = retval*v
    return retval
# end of function product


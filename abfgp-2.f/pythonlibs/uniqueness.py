"""
Functions to obtain (semi) unique strings/integers/tags 
"""
# Module metadata
__authors__ = "Ate van der Burgt"
__license__ = "MIT"

# Python Imports
from random import sample

# Global variables
ALPHABET = "abcdefghijklmnopqrstuvwxyz"

# translate string to get rid of non-alphanumerical characters in strings
# use like: "abcde1234!@%^".translate(NONLETTERS_TRANS)
# see: http://stackoverflow.com/questions/2171095/python-efficient-method-to-remove-all-non-letters-and-replace-them-with-unders
NONLETTERS_TRANS=''.join(chr(c) if chr(c).isupper() or chr(c).islower() else '_' for c in range(256))

def get_random_string_tag(length=8):
    """ """
    return "".join(sample(ALPHABET,length))

# end of function get_random_string_tag


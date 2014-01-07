"""
Functions for obtaining TMHMM data for Alignment Based Gene Prediction
"""

# Module metadata
__authors__ = "Ate van der Burgt"
__license__ = "MIT"

# Imports

# Global Variable Imports
USE_TMHMM = True

def obtainlocustmhmmdata(input):
    """ """
    for org in input.keys():
        if USE_TMHMM:
            # store the locus TMHMM predictions to the Orf objects
            input[org]['gldobj'].load_locus_tmhmm_to_orfs()

    # return the incremented input dict
    return input

# end of function obtainlocustmhmmdata
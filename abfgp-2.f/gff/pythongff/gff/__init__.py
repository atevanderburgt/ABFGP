"""
    BiRC GFF library   
"""
__version__ = "$Revision: 0.1 $"[11:-2]

# Imports
import sys, string, operator, math, os
from feature import Feature, FeatureInputError
from gff import GFF, GFF_IOError
#__all__ = ["Feature","GFF"]
__all__ = ["Feature","FeatureInputError","GFF","GFF_IOError"]


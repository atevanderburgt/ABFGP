"""
Scanning functions for Orf objects with PSSMs
"""

# Module metadata
__authors__ = "Ate van der Burgt"
__license__ = "MIT"

# Python Imports
from sys import path as sysPath

# make shure imports from the parental gene package will work,
# even if it is not specified (yet) in sys.path
from os.path import dirname as osPathDirname, abspath as osPathAbspath, split as osPathSplit
gene_package_dir = osPathSplit(osPathDirname(osPathAbspath(__file__)))[0]
if gene_package_dir not in sysPath: sysPath.append(gene_package_dir)

# Gene Imports
from gene_exceptions import InproperlyAppliedArgument
from splicesite import scan_pssm_splice_site
from start import scan_pssm_tss

# Import Global variables
from settings.splicesites import (
    IC_DONOR_PATTERN_OFFSET,
    IC_ACCEPTOR_PATTERN_OFFSET,
    )
from settings.translationalstartsites import (
    IC_TSS_PATTERN_OFFSET
    )

def scan_orf_for_pssm_splice_sites(orf,splicetype="donor",min_pssm_score=None,
    allow_non_canonical=False,non_canonical_min_pssm_score=0.0):
    """
    """
    if splicetype=='acceptor':
        pattern_offset  = IC_ACCEPTOR_PATTERN_OFFSET 
        offset_5p       = pattern_offset[0]+4 
        offset_3p       = 0
    elif splicetype=='donor':
        pattern_offset  = IC_DONOR_PATTERN_OFFSET 
        offset_5p       = 0
        offset_3p       = pattern_offset[1]+4 
    else:
        message = "'splicetype' argument '%s', not acceptor/donor" % splicetype
        raise InproperlyAppliedArgument, message

    # get (elongated) sequence of this orf
    seqslice = orf.nucleotidesequence(extra_start=offset_5p,extra_end=offset_3p)

    # scan for occurrences
    sites = scan_pssm_splice_site(seqslice,splicetype=splicetype,
        min_pssm_score=min_pssm_score,allow_non_canonical=allow_non_canonical,
        non_canonical_min_pssm_score=non_canonical_min_pssm_score)

    # correct site positions to absolute coords
    # and set correct phase of the splice site
    for site in sites:
        site.start = site.start + orf.startPY - offset_5p
        site.end   = site.end   + orf.startPY - offset_5p
        site.pos   = site.pos   + orf.startPY - offset_5p
        site.phase = (site.pos - orf.startPY) % 3

    # and return the splice sites on this orf
    return sites

# end of function scan_orf_for_pssm_splice_sites


def scan_orf_for_pssm_tss(orf,override_pattern_offset=(),min_pssm_score=None,
    allow_non_canonical=False,non_canonical_min_pssm_score=0.0):
    """
    """
    pattern_offset  = IC_TSS_PATTERN_OFFSET
    canonical       = "ATG"

    # hmm... somebody knows what he or she is doing ;-)
    if override_pattern_offset:
        pattern_offset = override_pattern_offset

    offset_5p       = pattern_offset[0]+2 
    offset_3p       = pattern_offset[1]+2

    # get (elongated) sequence of this orf
    seqslice = orf.nucleotidesequence(extra_start=offset_5p,extra_end=offset_3p)

    # scan for occurrences of TranslationalStartSites
    # REMEMBER: all 3 frames of this ORF are scanned
    # filter out wrong phases/frames afterwards
    sites = scan_pssm_tss(seqslice,override_pattern_offset=override_pattern_offset,
        min_pssm_score=min_pssm_score,allow_non_canonical=allow_non_canonical,
        non_canonical_min_pssm_score=non_canonical_min_pssm_score)

    correctsites = []
    # only return sites in the correct reading frame
    for site in sites:
        # check if the start of the site matches with current frame
        # TSS (ATG) are allways located in frame 0
        # site.start == position of ATG (0-based)
        #               minus offset_5p 
        if (site.start - offset_5p) % 3 == 0:
            # yep, this site is in the correct reading frame
            # and thus falls on this orf
            site.start = site.start + orf.startPY - offset_5p
            site.end   = site.end   + orf.startPY - offset_5p
            site.pos   = site.pos   + orf.startPY - offset_5p
            site.phase = orf.frame
            correctsites.append(site)

    # and return the TranslationalStartSites on this orf
    return correctsites

# end of function scan_orf_for_pssm_tss


"""
Functions concerning GFF files & lists of GFF tuples
"""

# Module metadata
__authors__ = "Ate van der Burgt"
__license__ = "MIT"

# Function Imports
from lib_fasta import _reversecomplement

# Pyhton Imports
from sets import Set
from copy import deepcopy
    

def __parsegfffile(gfffile,offset=0):
    """
    Parse a single gff file into a list of tracks (as tuples)

    @type  gfffile: string 
    @param gfffile: (absolute) gff filename

    @type  offset: integer
    @param offset: correct the positions in the gff tracks by this offset

    @rtype:  list
    @return: list with gff tracks as tuples (9 elements)
    """ 
    tracks = []
    for line in open(gfffile).readlines():
        if not line.strip(): continue
        tracks.append( _parse_line_2_track(line) )
    # return gff tracks list
    return tracks 

# end of function parsegfffile


def parsegfftxt(gfftxt,offset=0):
    """
    Parse a gff text block into a list of tracks (as tuples)

    @type  gfftxt: string 
    @param gfftxt: gff text block

    @type  offset: integer
    @param offset: correct the positions in the gff tracks by this offset

    @rtype:  list
    @return: list with gff tracks as tuples (9 elements)
    """ 
    tracks = []
    for line in gfftxt.strip().split("\n"):
        if not line.strip(): continue
        tracks.append( _parse_line_2_track(line) )
    # return gff tracks list
    return tracks 

# end of function parsegfftxt


def _parse_line_2_track(line,offset=0,separator="\t"):
    """
    Parse a single gff file into a list of tracks (as tuples)

    @type  line: string 
    @param line: gff text line (of 8 columns!)

    @type  offset: integer
    @param offset: correct the positions in the gff tracks by this offset

    @type  separator: string
    @param separator: string which separates the gff line (default tab)

    @rtype:  tuple
    @return: single gff line as a tuple (9 elements)
    """ 
    _gff = line.strip().split(separator)
    _gff[3] = int(_gff[3])-offset
    _gff[4] = int(_gff[4])-offset
    return tuple(_gff)

# end of function _parse_line_2_track


def parsegfffile(gfffile,offset=0):
    """
    Parse a single gff file into a list of tracks (as tuples)

    @type  gfffile: string 
    @param gfffile: (absolute) gff filename

    @type  offset: integer
    @param offset: correct the positions in the gff tracks by this offset

    @rtype:  list
    @return: list with gff tracks as tuples (9 elements)
    """ 
    tracks = []
    for line in open(gfffile).readlines():
        gff = line.strip().split("\t")
        gff[3] = int(gff[3])-offset
        gff[4] = int(gff[4])-offset
        tracks.append( tuple(gff) )
    # return gff tracks list
    return tracks 

# end of function parsegfffile


def mirrorgfflist(gfflist,length=0):
    """
    """
    strandtrans={ '+':'-', '-':'+', '.': '-' }
    for i in range(0,len(gfflist)):
        gffl = list(gfflist[i])
        gffl[4], gffl[3] = length-gffl[3]+1, length-gffl[4]+1
        gffl[6] = strandtrans[ gffl[6] ]
        gfflist[i] = tuple(gffl)
    # return the (ordered) mirrored gfflist
    gfflist.sort()
    return gfflist

# end of function mirrorgfflist


def order_gff_list(gfflist,reversed=False):
    """
    Order a list of gff tuples by their start coordinates

    @type  gfflist: list
    @param gfflist: list with gff tuples

    @type  reversed: Boolean
    @param reversed: reversed order when True (highest coordinate first)

    @rtype:  list
    @return: input gfflist ordered on start (4th element in list)
    """
    tmp = []
    for gff in gfflist: tmp.append( ( gff[0], int(gff[3]), gff ) )
    tmp.sort()
    if reversed: tmp.reverse()
    # return ordered list
    return [ gff for (fref, start, gff) in tmp ] 

# end of function order_gff_list


def filtergffs4fmethods(gffs,fmethods=[]):
    """
    Filter list of gff tuples for certain fmethods

    @type  gffs: list
    @param gffs: list with gff tuples

    @type  fmethods: list
    @param fmethods: list of strings of fmethods (column 3) to take into account

    @rtype:  list 
    @return: list with filtered gff tuples
    """
    retlist = []
    for gff in gffs:
        if gff[2] in fmethods:
            retlist.append(gff)
    return retlist

# end of function filtergffs4fmethods


def filtergffs4fmethod(gffs,fmethod=""):
    """
    Filter list of gff tuples for a specific fmethod

    @type  gffs: list
    @param gffs: list with gff tuples

    @type  fmethod: string
    @param fmethod: fmethod (column 3) to take into account
    
    @rtype:  list
    @return: list with filtered gff tuples
    """
    retlist = []
    for gff in gffs:
        if gff[2] == fmethod:
            retlist.append(gff)
    return retlist

# end of function filtergffs4fmethod


def filtergffs4fref(gffs,fref=""):
    """
    Filter list of gff tuples for certain fref

    @type  gffs: list
    @param gffs: list with gff tuples

    @type  fref: string 
    @param fref: fref (column 1) to take into account

    @rtype:  list 
    @return: list with filtered gff tuples
    """
    retlist = []
    for gff in gffs:
        if gff[0] == fref:
            retlist.append(gff)
    return retlist

# end of function filtergffs4fref



def gffs2coordset(gffs,fmethod=[],gclass=[]):
    """
    Return a coordinate set of a (subset of) list of gff tuples

    @type  gffs: list
    @param gffs: list with gff tuples

    @type  fmethod: list
    @param fmethod: list of strings of fmethods (column 3) to take into account

    @type  gclass: list
    @param gclass: list of strings of gclass's (column 9 first item) to take into account

    @rtype:  sets.Set
    @return: Set with coordinates
    """
    retset = Set()
    for gff in gffs:
        if gff[2] in fmethod:
            retset.update( range( int(gff[3]), int(gff[4])+1 ) )
            continue
        if gff[8].split(" ")[0] in gclass:
            retset.update( range( int(gff[3]), int(gff[4])+1 ) )
            continue
    # and return the Set of coordinates
    return retset

# end of function gffs2coordset


def coordset2gfftracks(coords,gfftrackbackbone):
    """
    Convert a list of (semi-continiuous) coordinates to GFF tracks

    @type  gfftrackbackbone: tuple
    @param gfftrackbackbone: backbone gff tuple in which coordinates are placed

    @type  coords: list
    @param coords: list with nt integer coordinates
    
    @rtype  gfftracks: list
    @return gfftracks: list with created gff tuples
    """
    # check if coords list is not empty
    if not coords: return []
    # return list for tracks
    gfftracks = []
    # data preparation (ordering & setting of first track first coord)
    coords.sort()
    track_coords  = [ [ coords[0] ] ]
    for coord in coords[1:]:
        if coord == max(track_coords[-1])+1:
            track_coords[-1].append(coord)
        else:
            track_coords.append( [ coord ] )
            
    for track in track_coords:
        # make a deepcopy of the gfftrackbackbone
        newgff = list( deepcopy( gfftrackbackbone ) )
        # update the coordinates
        newgff[3] = min(track)
        newgff[4] = max(track)
        gfftracks.append( tuple(newgff) )

    # return list with tracks
    return gfftracks

# end of function coordset2gfftracks


def exons2introns(exons,fmethod='intron',gclass="Intron"):
    """
    Return a coordinate set of a (subset of) list of gff tuples

    @type  exons: list
    @param exons: list with exons ordered on fstart

    @type  fmethod: string
    @param fmethod: gff fmethod string for the introns

    @type  gclass: string
    @param gclass: gff gclass string for the introns

    @rtype:  list
    @return: list with gff intron tuples

    @attention: Exon list is supposed to be ordered!!
    """
    introns = []
    for pos in range(1,len(exons)):
        exonA  = exons[pos-1]
        exonB  = exons[pos]
        fstart = int(exonA[4])+1 # end of previous exon
        fstop  = int(exonB[3])-1 # start of next exon
        gname = "Intron_%s_%s" % (fstart,fstop)
        intron = ( exonA[0], exonA[1], fmethod, fstart, fstop,
                    ".", exonA[6], ".", "%s %s" % (gclass,gname) )
        introns.append(intron)
    return introns

# end of function exons2introns


def gffs2txt(gffs,separator='\t'):
    """   """
    #data = []
    #for line in gffs:
    #     data.append( separator.join([str(x) for x in line ]) )
    #return "\n".join(data)
    return "\n".join([separator.join([str(x) for x in line ]) for line in gffs])

# end of function gffs2txt


def getseqbyfeaturetuple(sequence,gfftuple):
    """ """
    seqpart = sequence[int(gfftuple[3])-1:int(gfftuple[4])]
    if gfftuple[5] == '-':
        return _reversecomplement(seqpart)
    else:
        return seqpart

# end of function getseqbyfeaturetuple


def apply_fref(gffs,fref,fmethods=[]):
    """
    Apply given fref to all gff tracks (that have fiven fmethod(s))

    @type  gffs: list
    @param gffs: list with gff tuples

    @type  fref: string
    @param fref: fref string

    @type  fmethods: list
    @param fmethods: list of strings of fmethods (column 3) to take into account

    @rtype:  list 
    @return: list with rewritten gff tuples
    """
    retlist = []
    for gff in gffs:
        if not fmethods or gff[2] in fmethods:
            gff = list(gff)
            gff[0] = fref
            retlist.append(tuple(gff))
    return retlist

# end of function apply_fref


def merge_overlapping_features(gffs,identical_fmethod=True):
    """ """
    retlist = []
    if identical_fmethod:
        splittedgffs = []
        for fmethod in set([ gff[2] for gff in gffs ]):
            splittedgffs.append( order_gff_list( filtergffs4fmethod(gffs,fmethod=fmethod) ) )
    else:
        splittedgffs = [ order_gff_list( gffs ) ]
    for splitted in splittedgffs:
        mergedlist = [ splitted[0] ]
        for gff in splitted[1:]:
            if gff[0] != mergedlist[-1][0]:
                mergedlist.append(gff)
            elif gff[3] >= mergedlist[-1][3] and gff[3] <= mergedlist[-1][4] and gff[4] > mergedlist[-1][4]:
                merged = list(mergedlist[-1])
                merged[4] = gff[4]
                try:
                    merged[5] = max([merged[5],gff[5]])
                except:
                    merged[5] = "."
                while len(merged) < 8: merged.append(".")
                if len(merged) == 8: merged.append("")
                if merged[8].find("merged_feature \"True\";") == -1:
                    if merged[8] == "":
                        merged[8] = "merged_feature \"True\";"
                    else:
                        merged[8]+= "; merged_feature \"True\";"
                        merged[8].replace(";;",";")
                mergedlist[-1] = tuple(merged)
            else:
                mergedlist.append(gff)
        retlist.extend(mergedlist) 
    # return complete merged list
    return retlist

# end of function merge_overlapping_features

def gnamefromgfftuple(gfftuple):
    """
    @type  gfftuple: tuple 
    @param gfftuple: 9 element gff tuple 

    @rtype:  string 
    @return: gname (or empty string)
    """
    if len(gfftuple) != 9: return ''
    
    gname = gfftuple[8].split("; ")[0].split(" ",1)[1]
    if gname[0]+gname[1] == "''":
        gname = gname[1:-1]
    if gname[0]+gname[1] == '""':
        gname = gname[1:-1]
    return gname

# end of function gnamefromgfftuple

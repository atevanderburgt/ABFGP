"""
Splitting algorithms for python-pacb.

split_pacb_on_coordinates
_has_splitter
split_pacb_on_gaps
split_pacb_on_lowsimilarities
strip_unmatched_ends

"""

# Module metadata
__authors__ = "Ate van der Burgt"
__license__ = "MIT"

# Imports
import pacb
from re import search
from copy import deepcopy

# Global variables
MINIMUM_PACB_ALIGNMENT_AA_LENGTH = 4
SPLITTERS = [
        (  8, 0, 0 ),
        ( 10, 1, 1 ),
        ( 12, 2, 0 ),
        ( 13, 3, 2 ),
        ( 12, 1, 2 ),
        ( 12, 0, 4 ),
    ]



def split_pacb_on_coordinates(thepacbp,coords,returnside='left'):
    """
    thepacbp    pacbp object to be splitted at splitter coordinates
    coords      2-sized tuple with split coordinates
                first  int: end of left splitted pacb part
                second int: start of rigth splitted pacb part
    returnside  left    left part of splitted pacb, belonging to first splitter int coord
                rigth   rigth part of splitted pacb, belonging to second splitter int coord
    """
    (splitter_start,splitter_end) = coords
    if splitter_start and returnside == 'left' and \
    splitter_start >= MINIMUM_PACB_ALIGNMENT_AA_LENGTH:

        # check for occurrence of the splitter in the actually aligned region
        if thepacbp.__class__.__name__ == 'PacbPORF':
           if splitter_start <= thepacbp._original_alignment_pos_start:
               # a split BEFORE this pacbporf -> left side == None
               return None
           elif splitter_start >= thepacbp._original_alignment_pos_end:
               # a split AFTER this pacbporf -> left side == thepacbp
               return thepacbp
           else:
               # a split IN this pacbporf -> do actual splitting!
               pass

        # make new left side of the alignment
        pacbL           = deepcopy(thepacbp)
        pacbL.query     = pacbL.query[0:splitter_start]
        pacbL.match     = pacbL.match[0:splitter_start]
        pacbL.sbjct     = pacbL.sbjct[0:splitter_start]
        pacbL.alignment = pacbL.alignment[0:splitter_start]
        pacbL.length    = splitter_start
        # set the new end positions, take "-" (gaps) into account!
        pacbL.query_end = pacbL.query_start+pacbL.length-pacbL.query.count("-")
        pacbL.sbjct_end = pacbL.sbjct_start+pacbL.length-pacbL.sbjct.count("-")

        # split PacbPDNA and PacbPORF dna-sequences too!
        if pacbL.__class__.__name__ in ['PacbPDNA','PacbPORF']:
            # take a slice of the _positions list
            pacbL._positions = pacbL._positions[0:splitter_start]

            # query/sbjct_dna_end must be changed
            pacbL.query_dna_end     = pacbL._positions[-1].query_dna_end
            pacbL.sbjct_dna_end     = pacbL._positions[-1].sbjct_dna_end
            pacbL.query_dna_length  = pacbL.query_dna_end - pacbL.query_dna_start
            pacbL.sbjct_dna_length  = pacbL.sbjct_dna_end - pacbL.sbjct_dna_start
            pacbL.query_dna         = pacbL.query_dna[0:pacbL.query_dna_length]
            pacbL.sbjct_dna         = pacbL.sbjct_dna[0:pacbL.sbjct_dna_length]

        if pacbL.__class__.__name__ == 'PacbPORF':
            if len(pacbL._positions) < pacbL._original_alignment_pos_end:
                pacbL._original_alignment_pos_end = len(pacbL._positions)

        # and score the alignment!
        pacbL.strip_unmatched_ends()

        if pacbL.length >= MINIMUM_PACB_ALIGNMENT_AA_LENGTH:
            pacbL._score_alignment()
            # and return
            return pacbL
        else:
            return None

    elif splitter_end and returnside == 'rigth' and \
    (thepacbp.length-splitter_end) >= MINIMUM_PACB_ALIGNMENT_AA_LENGTH:

        # check for occurrence of the splitter in the actually aligned region
        if thepacbp.__class__.__name__ == 'PacbPORF':
           if splitter_start <= thepacbp._original_alignment_pos_start: 
               # a split BEFORE this pacbporf -> rigth side == thepacbp
               return thepacbp 
           elif splitter_start >= thepacbp._original_alignment_pos_end:
               # a split AFTER this pacbporf -> rigth side == None 
               return None 
           else:
               # a split IN this pacbporf -> do actual splitting!
               pass

        # make new rigth side of the alignment
        pacbR           = deepcopy(thepacbp)
        pacbR.query     = pacbR.query[splitter_end:]
        pacbR.match     = pacbR.match[splitter_end:]
        pacbR.sbjct     = pacbR.sbjct[splitter_end:]
        pacbR.alignment = pacbR.alignment[splitter_end:]
        pacbR.length    = len(pacbR.query)
        # set the new start positions, take "-" (gaps) into account!
        pacbR.query_start = pacbR.query_end-pacbR.length+pacbR.query.count("-")
        pacbR.sbjct_start = pacbR.sbjct_end-pacbR.length+pacbR.sbjct.count("-")
        # update the query_protein & sbjct_protein strins
        pacbR.query_protein = pacbR.query.replace("-","")
        pacbR.sbjct_protein = pacbR.sbjct.replace("-","")


        # split PacbPDNA and PacbPORF dna-sequences too!
        if pacbR.__class__.__name__ in ['PacbPDNA','PacbPORF']:
            # take a slice of the _positions list
            pacbR._positions  = pacbR._positions[splitter_end:]

            # renumber `position` counter for all existing positions
            for i in range(0,len(pacbR._positions)): pacbR._positions[i].position = i

            # query/sbjct_dna_start must be changed
            pacbR.query_dna_start   = pacbR._positions[0].query_dna_start
            pacbR.sbjct_dna_start   = pacbR._positions[0].sbjct_dna_start
            pacbR.query_dna_length  = pacbR.query_dna_end - pacbR.query_dna_start
            pacbR.sbjct_dna_length  = pacbR.sbjct_dna_end - pacbR.sbjct_dna_start
            pacbR.query_dna         = pacbR.query_dna[-pacbR.query_dna_length:]
            pacbR.sbjct_dna         = pacbR.sbjct_dna[-pacbR.sbjct_dna_length:]


        if pacbR.__class__.__name__ == 'PacbPORF':
            if splitter_end < thepacbp._original_alignment_pos_start:
                pacbR._original_alignment_pos_start -= splitter_end
                pacbR._original_alignment_pos_end -= splitter_end
            elif splitter_end < thepacbp._original_alignment_pos_end:
                pacbR._original_alignment_pos_start = 0
                pacbR._original_alignment_pos_end -= splitter_end
            else:
                print spitter_end, thepacbp._original_alignment_pos_start, thepacbp._original_alignment_pos_end
                raise "WHATEVER; splitting of pacbR failed because unexpected coordinates!"

        # and score the alignment!
        pacbR.strip_unmatched_ends()

        if pacbR.length >= MINIMUM_PACB_ALIGNMENT_AA_LENGTH:
            pacbR._score_alignment()
            # and return
            return pacbR
        else:
            return None


    elif splitter_start and returnside == 'left' and \
    splitter_start < MINIMUM_PACB_ALIGNMENT_AA_LENGTH:
        # to small left pacb to return!
        return None

    elif splitter_end and returnside == 'rigth' and \
    (thepacbp.length-splitter_end) < MINIMUM_PACB_ALIGNMENT_AA_LENGTH:
        # to small rigth pacb to return!
        return None

    elif returnside not in ['left','rigth']:
        message = "argument `returnside` must be either 'left' or 'rigth'"
        raise pacb.InproperlyAppliedArgument, message

    else:
        # if we reach this point, nothing to return!
        # means an unexpected/undesired splitting action
        return None

# end of function split_pacb_on_coordinates


def _has_splitter(pacbp,splitter,may_contain_methione=True):
    """
    splitter is composed of:
    (windowsize,identity,similarity)
    """
    splitter_start = -1
    splitter_end   = -1
    (windowsize,identity,similarity) = splitter
    for i in range(windowsize,len(pacbp.alignment)):
        window = pacbp.alignment[i-windowsize:i]
        if window.count("*") <= identity and window.count("+") <= similarity:
            if not may_contain_methione and pacbp.match[i-windowsize:i].find('M') > -1:
                # no aligned Methionies aloued in teh to-be-splitted region!
                # this is required in the potenial first pacbp
                continue
            # splitted position found
            (splitter_start, splitter_end) = ( i-windowsize, i )
            # correct left side for potentially last matched character
            windowlist = list(window)
            if splitter_start > 0:
                for char in windowlist:
                    if char != " ": splitter_start += 1
                    else:           break
            # correct right side for potentially last matched character
            windowlist.reverse()
            for char in windowlist:
                if char != " ": splitter_end -= 1
                else:           break
            # Splitter is fully prepared now.
            # Break out of for loop
            break
    # Return the splitter outcome
    return (splitter_start,splitter_end)

# end of function _has_splitter


def split_pacb_on_lowsimilarities(pacbp,splitters=SPLITTERS,may_contain_methione=True):
    """
    """
    # return list
    unsplittable_pacbs = []
    putatively_splittable_pacbs = [ pacbp ]
    # and a checker if something was splitted after all
    pacb_is_splitted = False
    # now do the split!
    cnt = 1
    while len(putatively_splittable_pacbs) >= 1:
        # list of new proposals for the next iteration
        new_putatively_splittable_pacbs = []
        for _this_pacb in putatively_splittable_pacbs:
            # checker for if there has been a split
            is_splitted = False
            for splitter in splitters:
                if is_splitted == True:
                    #already splitted, continue!
                    continue
                has_splitter = _has_splitter(_this_pacb,splitter,may_contain_methione=may_contain_methione)
                if has_splitter != (-1,-1):
                    # set the splitting-checkers to True
                    pacb_is_splitted = True
                    is_splitted = True
                    # make new left side of pacb
                    pacbL = split_pacb_on_coordinates(_this_pacb,has_splitter,returnside='left')
                    if pacbL and pacbL.length >= MINIMUM_PACB_ALIGNMENT_AA_LENGTH:
                        new_putatively_splittable_pacbs.append( pacbL )
                    # make new rigth side of pacb
                    pacbR = split_pacb_on_coordinates(_this_pacb,has_splitter,returnside='rigth')
                    if pacbR and pacbR.length >= MINIMUM_PACB_ALIGNMENT_AA_LENGTH:
                        new_putatively_splittable_pacbs.append( pacbR )
                    # irregardedly of what happened, jump out after the split!
                    break
                else:
                    # no split proposed for current split; no problem
                    pass
            else:
                # if we reach this point,
                # no split has been performed for this pacb
                # append this pacb to the non-splittable list
                unsplittable_pacbs.append( _this_pacb )

        # End of this round.
        # Set new_putatively_splittable_pacbs to putatively_splittable_pacbs
        # for the next iteration
        putatively_splittable_pacbs = deepcopy(new_putatively_splittable_pacbs)
        new_putatively_splittable_pacbs = []

    # and return all that can not be further splitted
    return (unsplittable_pacbs, pacb_is_splitted)

# end of function split_pacb_on_lowsimilarities


def pacb_slice_abs_coords(pacb,start=None,end=None,expressed_on="query"):
    """
    Create a slice from a PacbP (or inheriting) object on absolute coordinates

    @attention: parameters start and end must fall within PacbP's range

    @type pacbp:    PacbP (or inheriting) object
    @param pacbp:   PacbP (or inheriting) object

    @type start:    number
    @param start:   absolute AA start position of the slice

    @type end:  number
    @param end: absolute AA end position of the slice

    @type expressed_on:     string
    @param expressed_on:    'query' or 'sbjct'

    @rtype:  PacbP (or inheriting) object
    @return: PacbP (or inheriting) object
    """
    # do some parameter sanity checks
    if expressed_on == 'query':
        if start == None: start = pacb.query_start
        if end   == None: end   = pacb.query_end
        for coord in [start,end]:
            if coord not in range(pacb.query_start,pacb.query_end+1):
                raise CoordinateOutOfRange
        measure_coords = {'start': pacb.query_start, 'end': pacb.query_end+1}
        measure_alignm = list(pacb.query)
    elif expressed_on == 'sbjct':
        if start == None: start = pacb.query_start
        if end   == None: end   = pacb.query_end
        for coord in [start,end]:
            if coord not in range(pacb.sbjct_start,pacb.sbjct_end+1):
                raise CoordinateOutOfRange
        measure_coords = {'start': pacb.sbjct_start, 'end': pacb.sbjct_end+1}
        measure_alignm = list(pacb.sbjct)
    else:
        message = "`expressed_on` not in ['query','sbjct']: '%s'" % expressed_on
        raise InproperlyAppliedArgument(message)

    # now translate absolute coords to relative coords
    rel_coord_start, rel_coord_end = 0, 0
    while start >= measure_coords['start']:
        aa = measure_alignm.pop(0)
        rel_coord_start += 1
        if aa != "-":
            measure_coords['start']+=1
            if measure_coords['start'] > start:
                # this is the final step upon breaking the while clause
                rel_coord_start -= 1
    while end <= measure_coords['end']:
        aa = measure_alignm.pop()
        rel_coord_end -= 1
        if aa != "-":
            measure_coords['end']-=1
            if measure_coords['end'] < end:
                # this is the final step upon breaking the while clause
                rel_coord_end += 1

    # and make the slice by the relative coordinates
    return pacb_slice_rel_coords(pacb,start=rel_coord_start,end=-rel_coord_end)

# end of function pacb_slice_abs_coords


def pacb_slice_rel_coords(pacb,start=0,end=0):
    """
	Create a slice from a PacbP (or inheriting) object on relative coordinates

	@attention: parameters start and end must fall within PacbP's range

	@type  pacbp: PacbP (or inheriting) object
	@param pacbp: PacbP (or inheriting) object

	@type  start: number
	@param start: relative AA start position of the slice: 0 .. len(PacbP)

	@type  end: number
	@param end: relative AA end position of the slice: -(len(PacbP)-1) .. 0

	@rtype:  PacbP (or inheriting) object
	@return: PacbP (or inheriting) object
    """
    # do some parameter sanity checks
    if start not in range(0,pacb.length):    raise CoordinateOutOfRange
    if end not in range(-(pacb.length-1),1): raise CoordinateOutOfRange

    if start == 0 and end == 0:
        return pacb
    elif start == 0:
        coords = (pacb.length+end,pacb.length+end)
        pacbL = split_pacb_on_coordinates(pacb,coord,returnside='left')
        return pacbL
    elif end == 0:
        coords = (start,start)
        pacbR = split_pacb_on_coordinates(pacb,coords,returnside='rigth')
        return pacbR
    else:
        pacbL = split_pacb_on_coordinates(pacb,(end,end),returnside='left')
        pacbR = split_pacb_on_coordinates(pacbL,(start,start),returnside='rigth')
        return pacbR

# end of function pacb_slice_rel_coords


def split_pacb_on_gaps(thepacbp,gapsize=2):
    """
    NEW function for splitting pacbps on gaps
    OLD function had errors.....
    """
    # find all positions that has to be splitted
    if not thepacbp.alignment_has_gaps(gap_size=gapsize):
        return ( [ thepacbp ], False )
    else:
        splitted_pacbs = []
        gap_split_coords = []

        # list of to-be-splitted and unspittable ones
        tobesplitted = [ thepacbp ]
        returnunsplittable = []

        while tobesplitted:
            furthersplittable  = []
            for part in tobesplitted:
                splitcoords = ()
                for strand in [ part.query, part.sbjct ]:
                    match = search("-{%s,}" % gapsize,strand)
                    if match:
                        splitcoords = ( match.start(), match.end() )
                        break
                # check if there is a split possible
                if splitcoords:
                    # create left and right pacbp
                    pacbL = split_pacb_on_coordinates(part,splitcoords,returnside='left')
                    pacbR = split_pacb_on_coordinates(part,splitcoords,returnside='rigth')

                    ####print "pacbL", (pacbL._original_alignment_pos_start, pacbL._original_alignment_pos_end), min([ pos.position for pos in pacbL._positions ]), max([ pos.position for pos in pacbL._positions ]), len(pacbL._positions)
                    ####print "pacbL", [ pos.position for pos in pacbL._positions ]
                    ####print ""
                    ####print "pacbR", (pacbR._original_alignment_pos_start, pacbR._original_alignment_pos_end), min([ pos.position for pos in pacbR._positions ]), max([ pos.position for pos in pacbR._positions ]), len(pacbR._positions)
                    ####print "pacbR", [ pos.position for pos in pacbR._positions ]
                    ####print ""

                    # evaluate the LEFT part
                    if pacbL and pacbL.length >= MINIMUM_PACB_ALIGNMENT_AA_LENGTH:
                        if pacbL.alignment_has_gaps(gap_size=gapsize) and pacbL != part:
                            # proceed for the next iteration of splitting
                            furthersplittable.append( pacbL )
                        else:
                            # nope, now this part is unsplittable!
                            returnunsplittable.append( pacbL )
                    # evaluate the RIGHT part
                    if pacbR and pacbR.length >= MINIMUM_PACB_ALIGNMENT_AA_LENGTH:
                        if pacbR.alignment_has_gaps(gap_size=gapsize) and pacbR != part:
                            # proceed for the next iteration of splitting
                            furthersplittable.append( pacbR )
                        else:
                            # nope, now this part is unsplittable!
                            returnunsplittable.append( pacbR )
                else:
                    # this part is not splittable any further
                    returnunsplittable.append( part )
            # prepare for the next iteration
            tobesplitted = furthersplittable

        # return list with splitted pacb's
        if returnunsplittable:  splitted = True
        else:                   splitted = False
        return ( returnunsplittable, splitted )

# end of function split_pacb_on_gaps


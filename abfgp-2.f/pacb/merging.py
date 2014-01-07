"""
Merging algorithms for python-pacb.
Only possible for PacbPORF class
"""

# Module metadata
__authors__ = "Ate van der Burgt"
__license__ = "MIT"

# Pacb Imports
from alignedproteindnaposition import IsIdenticalAlignedProteinDnaPosition
from comparison import IsPacbPORFSubstringOf
from exceptions import * 

# Python Imports
from copy import deepcopy
from lib_stopwatch import StopWatch

def merge_pacbporfs(old,new,side,verbose=False):
    """
    Amalgate 2 (near) identical pacbporfs with a slighly different alignment on the right/3p side

    @type  old: PacbPORF
    @param old: PacbPORF instance

    @type  new: PacbPORF
    @param new: PacbPORF instance to take the rigth side from to amalgate onto the other PacbPORF

    @attention: assumes 'new' PacbPORF to be different/longer on the specified side

    @rtype:  PacbPORF
    @return: merged/amalgamated PacbPORF object

    @rtype:  Boolean
    @return: status weather or not the PacbPORFs are merged 
    """
    if side not in ['left','rigth']:
        message = "'side' not left or rigth but '%s'" % side
        raise InproperlyAppliedArgument, message

    stw = StopWatch()
    stw.start()

    # make temporarily unextended PacbPORFs
    old_was_extended, new_was_extended = False, False
    if old.is_extended():
        old.unextend_pacbporf()
        old_was_extended = True
    if new.is_extended():
        new.unextend_pacbporf()
        new_was_extended = True

    # call the merging helper function for the proper side
    if side == 'rigth':
        # do not get log messages from the helper function; it works fine...
        merged, status = _merge_pacbporfs_on_rigth_side(old,new,verbose=False)
        ###if verbose:
        ###    if status: print merged, "MERGED RIGTH"
        ###    else: print "FALSE MERGING rigth STATUS"
    else:
        # do not get log messages from the helper function; it works fine...
        merged, status = _merge_pacbporfs_on_left_side(old,new,verbose=False)
        ###if verbose:
        ###    if status: print merged, "MERGED LEFT"
        ###    else: print "FALSE MERGING left STATUS"

    # reset old and new PacbPORF object to its extention status
    if old_was_extended: merged.extend_pacbporf_after_stops()
    if old_was_extended: old.extend_pacbporf_after_stops()
    if new_was_extended: new.extend_pacbporf_after_stops()
    # verbose logging message when merging was succesfull
    if verbose and status:
         print old, "STARTING SITUATION", side
         #old.print_protein(_linesize=150)
         print merged, "RE-EXTENDED", side
         #merged.print_protein(_linesize=150)

    # return the (freshly merged) pacbporf and the status if it was succesfull
    return merged, status

# end of function merge_pacbporfs


def _merge_pacbporfs_on_rigth_side(old,new,verbose=False):
    """
    Helper function for merge_pacbporfs; rigth side of `new` is placed on `old`

    @rtype:  PacbPORF
    @return: merged/amalgamated PacbPORF object

    @rtype:  Boolean 
    @return: status weather merging was succesfull
    """

    # get the original alignment end objects of the current PacbPORF
    epos = old._get_original_alignment_pos_end()

    # get the AlignedProteinDnaPosition of the old PacbPORF end positions
    # in the new/alternative PacbPORF object
    eposInNewAlignment = new.alignmentobject_by_queryaa(epos.query_pos)

    # find positions from where they are identical
    pointers = (None,None)
    for oldpointer in range(-1,-len(old._positions)-1,-1):
        if old._positions[oldpointer].isa == 'identity':
            oldAlgObj = old._positions[oldpointer]
            for newpointer in range(-1,-len(new._positions)-1,-1):
                if new._positions[newpointer].isa == 'identity':
                    newAlgObj = new._positions[newpointer]
                    if IsIdenticalAlignedProteinDnaPosition(oldAlgObj,newAlgObj):
                        pointers = (oldpointer,newpointer)
                        break
                    elif newAlgObj.query_pos < oldAlgObj.query_pos:
                        break
                    else:
                        continue
                else:
                    continue
            if pointers != (None,None):
                break
        else:
            continue
    # if 100% identical pointers are found based on identical AAs,
    # try to move this merging point to the rigth (similar or dissimilar,
    # but the same AAs)
    if pointers != (None,None):
        while True:
            oldpointer, newpointer = pointers
            try:
                oldAlgObj = old._positions[oldpointer+1]
                newAlgObj = new._positions[newpointer+1]
                if IsIdenticalAlignedProteinDnaPosition(oldAlgObj,newAlgObj):
                    pointers = ( oldpointer+1, newpointer+1 )
                else:
                    break
            except:
                break


    # if pointers are found -> merge. If not, return the unmergable PacbPORF object
    if pointers == (None,None):
        # no merging is possible
        return old, False

    # if this point is reached -> start merging pacbporf objects
    if verbose:
        print "MERGING POINTERS rigth:", pointers
        print old, "OLD UNEXTENDED"
        print new, "NEW UNEXTENDED"

    # first, make a deepcopy of 'old' because it is about to be changed drastically! 
    merged = deepcopy(old)
    oldpointer, newpointer = pointers

    # start merging by removing the non-identical part on the rigth side of 'old'
    merged.pop_alignment_object(aas=abs(oldpointer)-1)
    if verbose: print merged, "POPPED"

    if IsPacbPORFSubstringOf(merged,new): print "substring PACBPORFS!!!!! before merging"

    # now append the rigth part of 'new' onto the current alignment
    # take list-slice of new._positions; INcrease the negative variable
    # newpointer by one to point to the next AlignedProteinDnaPosition object
    merged.append_alignment_object(  new._positions[newpointer+1:] )
    if verbose: print merged, "APPENDED"

    # return the freshly merged pacbporf
    return merged, True

# end of function _merge_pacbporfs_on_rigth_side


def _merge_pacbporfs_on_left_side(old,new,verbose=False):
    """
    Helper function for merge_pacbporfs; left side of `new` is placed on `old`

    @rtype:  PacbPORF
    @return: merged/amalgamated PacbPORF object

    @rtype:  Boolean
    @return: status weather merging was succesfull
    """

    # get the original alignment start object of the current PacbPORF
    spos = old._get_original_alignment_pos_start()

    # get the AlignedProteinDnaPosition of the old PacbPORF start positions
    # in the new/alternative PacbPORF object
    sposInNewAlignment = new.alignmentobject_by_queryaa(spos.query_pos)

    # find positions from where they are identical
    pointers = (None,None)
    for oldpointer in range(0,len(old._positions)):
        if old._positions[oldpointer].isa == 'identity':
            oldAlgObj = old._positions[oldpointer]
            for newpointer in range(0,len(new._positions)):
                if new._positions[newpointer].isa == 'identity':
                    newAlgObj = new._positions[newpointer]
                    if IsIdenticalAlignedProteinDnaPosition(oldAlgObj,newAlgObj):
                        pointers = (oldpointer,newpointer)
                        break
                    elif newAlgObj.query_pos > oldAlgObj.query_pos:
                        break
                    else:
                        continue
                else:
                    continue
            if pointers != (None,None):
                break
        else:
            continue

    # if 100% identical pointers are found based on identical AAs,
    # try to move this merging point to the left (similar or dissimilar,
    # but the same AAs)
    if pointers != (None,None):
        while True:
            oldpointer, newpointer = pointers
            try:
                oldAlgObj = old._positions[oldpointer-1]
                newAlgObj = new._positions[newpointer-1]
                if IsIdenticalAlignedProteinDnaPosition(oldAlgObj,newAlgObj):
                    pointers = ( oldpointer-1, newpointer-1 )
                else:
                    break
            except:
                break

    # if pointers are found -> merge. If not, return the unmergable PacbPORF object
    if pointers == (None,None):
        # no merging is possible
        return old, False

    # if this point is reached -> start merging pacbporf objects
    if verbose:
        print "MERGING POINTERS left:", pointers
        print old, "OLD UNEXTENDED"
        print new, "NEW UNEXTENDED"

    # first, make a deepcopy of 'old' because it is about to be changed drastically!
    merged = deepcopy(old)
    oldpointer, newpointer = pointers

    # start merging by removing the non-identical part on the left side of 'old'
    status = merged.shift_alignment_object(aas=oldpointer-1)
    if verbose: print merged, "SHIFTED", status

    if IsPacbPORFSubstringOf(merged,new): print "substring PACBPORFS!!!!! before merging"

    # now append the rigth part of 'new' onto the current alignment
    # take list-slice of new._positions; INcrease the negative variable
    # newpointer by one to point to the next AlignedProteinDnaPosition object
    status = merged.unshift_alignment_object(  new._positions[0:newpointer] )
    if verbose: print merged, "UNSHIFTED", status

    # return the freshly merged pacbporf
    return merged, True

# end of function _merge_pacbporfs_on_left_side

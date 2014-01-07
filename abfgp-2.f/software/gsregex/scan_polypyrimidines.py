#!/usr/bin/env python

import sys
from re import finditer,compile,IGNORECASE
# according to Kupfer e.a. 2004; Introns and splicing elements in five diverse fungi
REPAT_PYRIMIDINES = compile("[^A]{6,}",IGNORECASE)
MIN_TRACT_LENGTH  = 6
MIN_TRACT_T_COUNT = 3
SCAN_PYRIMIDINES_FMETHOD = 'PPTregex'
SCAN_PYRIMIDINES_FSOURCE = 'PPTregex'


def scan_polypyrimidine(seq,min_tract_length=MIN_TRACT_LENGTH,min_tract_t_count=MIN_TRACT_T_COUNT):
    """ """
    coords = []
    for m in finditer(REPAT_PYRIMIDINES,seq):
        if m.end()-m.start() >= min_tract_length and (m.group().count("T")+m.group().count("t")) >= min_tract_t_count:
            coords.append((m.start(),m.end()))
    return coords
# end of function scan_polypyrimidine



#REPAT_PYRIMIDINES = compile("[^A]*T[^A]*T[^A]*T[^A]*")
#MIN_TRACK_LENGTH  = 6
#def scan_polypyrimidine(seq,min_track_length=MIN_TRACK_LENGTH):
#    """ """
#    coords = []
#    for m in finditer(REPAT_PYRIMIDINES,seq):
#        if m.end()-m.start() >= min_track_length:
#            coords.append((m.start(),m.end()))
#    return coords
## end of function scan_polypyrimidine

def results_to_gff(header,coords):
    """ """
    for start,end in coords:
        print "%s\t%s\t%s\t%s\t%s\t%s\t+\t." % (
            header,
            SCAN_PYRIMIDINES_FSOURCE,
            SCAN_PYRIMIDINES_FMETHOD,
            start+1,
            end,
            end-start,
            )

# end of function results_to_gff

def main(min_tract_length=MIN_TRACT_LENGTH):
    """ """
    # check if MIN_TRACT_LENGTH is applied as command line argument
    try:    min_tract_length = int(sys.argv[1])
    except: pass
    # init header and seqparts variables
    header,seqparts = None,[]
    line = sys.stdin.readline()
    while line: 
        if line and line[0] == '>':
            if header and seqparts:
                results_to_gff( header, scan_polypyrimidine(
                    "".join(seqparts).upper().replace("U","T"),
                    min_tract_length=min_tract_length
                    ) )
                header,seqparts = None,[]
            # reset header
            header = line.strip().split(' ')[0].replace(">","")
        else:
            seqparts.append(line.strip())
        line = sys.stdin.readline()
    else:
        if header and seqparts:
            results_to_gff( header, scan_polypyrimidine(
                "".join(seqparts).upper().replace("U","T"),
                min_tract_length=min_tract_length
                ) )
# end of function main

if __name__ == "__main__":
    main()


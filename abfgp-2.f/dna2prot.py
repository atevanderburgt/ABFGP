from sets import Set

def dna2protein(seq,error_check=False):
    """
    """
    # sequence string in proper format
    seq = seq.upper().replace("U","T")

    # define the standard_genetic_code dictionary
    standard_genetic_code =     {"TTT" : "F" ,
                                 "TCT" : "S" ,
                                 "TAT" : "Y" ,
                                 "TGT" : "C" ,
                                 "TTC" : "F" ,
                                 "TCC" : "S" ,
                                 "TAC" : "Y" ,
                                 "TGC" : "C" ,
                                 "TTA" : "L" ,
                                 "TCA" : "S" ,
                                 "TAA" : "*" ,
                                 "TGA" : "*" ,
                                 "TTG" : "L" ,
                                 "TCG" : "S" ,
                                 "TAG" : "*" ,
                                 "TGG" : "W" ,
                                 "CTT" : "L" ,
                                 "CCT" : "P" ,
                                 "CAT" : "H" ,
                                 "CGT" : "R" ,
                                 "CTC" : "L" ,
                                 "CCC" : "P" ,
                                 "CAC" : "H" ,
                                 "CGC" : "R" ,
                                 "CTA" : "L" ,
                                 "CCA" : "P" ,
                                 "CAA" : "Q" ,
                                 "CGA" : "R" ,
                                 "CTG" : "L" ,
                                 "CCG" : "P" ,
                                 "CAG" : "Q" ,
                                 "CGG" : "R" ,
                                 "ATT" : "I" ,
                                 "ACT" : "T" ,
                                 "AAT" : "N" ,
                                 "AGT" : "S" ,
                                 "ATC" : "I" ,
                                 "ACC" : "T" ,
                                 "AAC" : "N" ,
                                 "AGC" : "S" ,
                                 "ATA" : "I" ,
                                 "ACA" : "T" ,
                                 "AAA" : "K" ,
                                 "AGA" : "R" ,
                                 "ATG" : "M" ,
                                 "ACG" : "T" ,
                                 "AAG" : "K" ,
                                 "AGG" : "R" ,
                                 "GTT" : "V" ,
                                 "GCT" : "A" ,
                                 "GAT" : "D" ,
                                 "GGT" : "G" ,
                                 "GTC" : "V" ,
                                 "GCC" : "A" ,
                                 "GAC" : "D" ,
                                 "GGC" : "G" ,
                                 "GTA" : "V" ,
                                 "GCA" : "A" ,
                                 "GAA" : "E" ,
                                 "GGA" : "G" ,
                                 "GTG" : "V" ,
                                 "GCG" : "A" ,
                                 "GAG" : "E" ,
                                 "GGG" : "G" ,
                                 }

    # and translate!
    prot = []
    for pos in range(0,len(seq),3):
        try:
            prot.append( standard_genetic_code[ seq[pos:pos+3] ] )
        except:
            # triplet not in genetic code! append 'X'
            prot.append( 'X' )

    # and return!
    return "".join(prot)

# end of function dna2protein


def dna2proteinbyframe(_sequence,_frame):
    """
    Translate a DNA sequence into AAs for a given frame

    @type  _sequence: string
    @param _sequence: DNA seqeunce string

    @type  _frame: integer
    @param _frame: desired reading frame (0,1,2)

    @rtype:  string
    @return: Protein sequence
    """
    _untill = (len(_sequence)-_frame) % 3
    if not _untill:  _untill = len(_sequence)
    else:            _untill = -_untill
    # and do the actual translation
    return dna2protein(_sequence[_frame:_untill])

# end of function dna2proteinbyframe


"""
Functions needed to determine the sequence periodicity

for length in [30,31,32]:
  x = "x"*length
  for frame in [0,1,2]:
    length, frame, range(frame,len(x)-((len(x)-frame) % 3),3)

from lib_sequenceperiodicity import test
test()

"""

from sets import Set



def test():
    """
    """
    seqA = """
        GTATGGTTCTTGGGTAAGACCGTTGGCCTCGTTGAGCAACAGTATAATGTCATCAAAGGACAAATGCCAT
        CGGTGCCCATGAGACTGCTTACGGGGCAGCTCAACATCGATGCATGGTCGCCA---GGGGTTTGGCCCCA
        AATTCTCGAGGGCACTCGCATCATCGTATCTACATTCTCCATCTTGCGTGATGCTCTTGACCACGCCTTT
        GTCACAATGGACATGCTCGCTTTAATAGTAGTTGACGAA"""
    seqB = """
        ATCTGGTTCCTGACGCCAACGGTCGCCCTGGCTCGCCAACAACACCGAGTTTTGCAATCACAAATACCTT
        CCGTCAAGGCTATAATGCTCTGTGGCCAGGATGGGGTGGACAGCTGGTCTGAGCAGGCTGTCTGGGATGC
        CGTCTTGCTCAACGTTCGTATCGTGGTGTCGACATACCAGATCTTGTTCGATGCAAATGCCCACAGCTTC
        GTCCGCTTGGACTCCTTGAGCTTGATTGTCATCGATGAA"""
    seqA = seqA.replace(" ","").replace("\n","")
    seqB = seqB.replace(" ","").replace("\n","")

    buffer_seqA = []
    buffer_cons = []
    buffer_seqB = []
    buffer_aaA  = []
    buffer_aaB  = []
    buffer_mtch = []

    geneticcode = importgeneticcode()
    seqlen = len(seqA)
    frame = 0
    for offset in range(frame,seqlen-((seqlen-frame) % 3),3):
        tripletA = seqA[offset:offset+3]
        tripletB = seqB[offset:offset+3]
        try:    aaA = geneticcode[tripletA]
        except: aaA = "+"
        try:    aaB = geneticcode[tripletB]
        except: aaB = "-"
        tripletA = list(tripletA)
        tripletB = list(tripletB)
        buffer_aaA.append( " " )
        buffer_aaB.append( " " )
        buffer_mtch.append( " " )
        is_conserved = is_conserved_aa(aaA,aaB)
        if is_conserved:
            if aaA == aaB:  buffer_mtch.append(aaA)
            else:           buffer_mtch.append("+")
        else:
            buffer_mtch.append(" ")
        buffer_aaA.append( aaA.replace("+"," ") )
        buffer_aaB.append( aaB.replace("-"," ") )
        buffer_aaA.append( " " )
        buffer_aaB.append( " " )
        buffer_mtch.append( " " )
        for pos in [0,1,2]:
            buffer_seqA.append( tripletA[pos] )
            buffer_seqB.append( tripletB[pos] )
            if tripletA[pos] == tripletB[pos]:
                buffer_cons.append( str(pos+1) )
            else:
                buffer_cons.append( " " )
                
    printsize = 99
    for offset in range(0,seqlen,printsize):
        for line in (buffer_seqA,buffer_cons,buffer_seqB,buffer_aaA,buffer_mtch,buffer_aaB):
            print "".join(line[offset:offset+printsize])
        print ""

    print "# codon_nucleotide_conservation(seqA,seqB)"
    cnc = codon_nucleotide_conservation(seqA,seqB)
    for frame in [0,1,2]:
        print "FRAME:", frame
        for item in cnc[frame]:
            print "(%1.4f, %1.4f, %1.4f, %s)" % item

    print "# coding_sequence_periodicity_test(seqA,seqB)"
    print coding_sequence_periodicity_test(seqA,seqB)


def coding_sequence_periodicity_test(seqA,seqB):
    """
    """
    cnc = codon_nucleotide_conservation(seqA,seqB,allframes=True)
    is_lowest_conservation_of_3_positions = False
    is_highest_ratio_of_3_frames = False

    if cnc[0][1][2] == min(cnc[0][1][0:3]):
        is_lowest_conservation_of_3_positions = True
    divide_by = [ cnc[frame][1][2] for frame in [0,1,2] ]
    if 0.0 in divide_by:
        # hmm, ZeroDivisionError is about to happen
        try:
            while True:
                pos = divide_by.index(0.0)
                # set to a very small offset
                divide_by[pos] = 0.05
        except ValueError:
            pass
    # and calulate the ratios
    ratios = [ (cnc[frame][1][0]+cnc[frame][1][1])/(2.0*divide_by[frame]) for frame in [0,1,2] ]
    if ratios[0] == max(ratios):
        is_highest_ratio_of_3_frames = True
    return ( is_lowest_conservation_of_3_positions, is_highest_ratio_of_3_frames )



def codon_nucleotide_conservation(seqA,seqB,verbose=True,allframes=True):
    """
    seqA    input dna sequence A
    seqB    input dna sequence B
    Requirements
    seqA    contains [A,T,G,C,-]
    seqB    contains [A,T,G,C,-]
    len(seqA) == len(seqB)
    """

    geneticcode = importgeneticcode()

    seqA = seqA.upper().strip().replace(" ","").replace("\n","")
    seqB = seqB.upper().strip().replace(" ","").replace("\n","")
    seqlen = len(seqA)
    framescores = {}
    doframes = [0,1,2]
    if not allframes: doframes = [0]
    for frame in doframes:
        pos_conservation_counts = [ 0, 0, 0 ]
        pos_conserv_counts    = [ 0, 0, 0, 0 ]
        pos_nonconserv_counts = [ 0, 0, 0, 0 ]
        for offset in range(frame,seqlen-((seqlen-frame) % 3),3):
            tripletA = seqA[offset:offset+3]
            tripletB = seqB[offset:offset+3]
            try:    aaA = geneticcode[tripletA]
            except: aaA = "+"
            try:    aaB = geneticcode[tripletB]
            except: aaB = "-"
            tripletA = list(tripletA)
            tripletB = list(tripletB)
            is_conserved = is_conserved_aa(aaA,aaB)
            for pos in [0,1,2]:
                if tripletA[pos] == tripletB[pos]:
                    pos_conservation_counts[pos]+=1
                    # and now check if this is a conserved
                    # position or not
                    if is_conserved:
                        pos_conserv_counts[pos] += 1
                    else:
                        pos_nonconserv_counts[pos] += 1
            # and increase the main counter
            if is_conserved:
                pos_conserv_counts[-1] += 1
            else:
                pos_nonconserv_counts[-1] += 1

        # and make the score tuple
        if (seqlen-frame)/3:
            # sequence of ample length; short seqlenghts <=5 give obvious problems (ZeroDivisionError)
            pos_conservation_score = [ float(cnt)/float((seqlen-frame)/3) for cnt in pos_conservation_counts ]
            pos_conservation_score.append( (seqlen-frame)/3 )
            pos_conservation_score = tuple( pos_conservation_score )
        else:
            # very short sequence <= 5 nt -> no counts in this frame!
            pos_conservation_score = ( 0.0, 0.0, 0.0, 0 )
        # and make the conserved/nonconserved tuples
        for i in range(0,3):
            try:
                pos_conserv_counts[i]    = float(pos_conserv_counts[i])/float(pos_conserv_counts[-1])
            except:
                pos_conserv_counts[i]    = 0.0
            try:
                pos_nonconserv_counts[i] = float(pos_nonconserv_counts[i])/float(pos_nonconserv_counts[-1])
            except:
                pos_nonconserv_counts[i]    = 0.0

        # and make tuples
        pos_conserv_counts    = tuple(pos_conserv_counts)
        pos_nonconserv_counts = tuple(pos_nonconserv_counts)
        framescores[frame] = ( pos_conservation_score, pos_conserv_counts, pos_nonconserv_counts )

    # and return all framescores
    if not allframes:
        return framescores[0]
    else:
        return framescores


def is_conserved_aa(aaA,aaB):
    """ """
    # check if combi is identical
    if aaA == aaB:
        return True
    else:
        # import conserved ones
        conservedsubstitutions = importconservedsubstitutions()
        # check if this combi is a conserved one
        for alias, substitutionset in conservedsubstitutions.iteritems():
            if Set([aaA,aaB]).issubset(substitutionset):
                return True
        else:
            return False


def importconservedsubstitutions():
    """ """
    conservedsubstitutions = {
        'Small (small+ hydrophobic (incl.aromatic -Y))':
            Set(list('AVFPMILW')),
        'Acidic':
            Set(list('DE')),
        'Basic - H':
            Set(list('RK')),
        'Hydroxyl + sulfhydryl + amine + G':
            Set(list('STYHCNGQ')), 
    }
    return conservedsubstitutions


def get_triplets_by_aa(aa,gencode={}):
    """
    """
    if not gencode: gencode = importgeneticcode()
    aa = aa.upper()
    triplets = []
    if aa not in gencode.values():
        return triplets
    for triplet,aatrans in gencode.iteritems():
        if aatrans == aa:
            triplets.append(triplet)
    return triplets

# end of function get_triplets_by_aa

def tripletsdegeneracy(triplets):
    """
    """
    f0,f1,f2 = [],[],[]
    for triplet in triplets:
        if triplet[0] not in f0: f0.append(triplet[0])
        if triplet[1] not in f1: f1.append(triplet[1])
        if triplet[2] not in f2: f2.append(triplet[2])
    f0.sort()
    f1.sort()
    f2.sort()
    return ( tuple(f0), tuple(f1), tuple(f2) )
    
# end of function tripletsdegeneracy

    
def gencodeperiodicity(dnaseq):
    """
    """
    geneticcode = importgeneticcode()
    dnaseq = dnaseq.upper()
    seqlen = len(dnaseq)
    frame = 0
    at = float( dnaseq.count('A') + dnaseq.count('T') ) / seqlen
    gc = float( dnaseq.count('G') + dnaseq.count('C') ) / seqlen
    # PeriodicityMaskedSequence
    pms = list(dnaseq)
    gcp = []
    for base in ['A','T','G','C']:
        for offset in range(frame,seqlen-((seqlen-frame) % 3),3):
            triplet  = dnaseq[offset:offset+3]
            try:
                aa   = geneticcode[triplet]
            except:
                # N or other abbarant nucleotides
                continue
            triplets = get_triplets_by_aa(aa,gencode=geneticcode)
            tripdeg  = tripletsdegeneracy(triplets)
            for pos in [0,1,2]:
                if triplet[pos] == base:
                    if len(tripdeg[pos]) == 1:
                        # obligate BASE -> leave upper
                        pass
                    else:
                        # non-obligate BASE -> lower()
                        pms[offset+pos] = pms[offset+pos].lower()
                else:
                    if len(tripdeg[pos]) == 1:
                        # obligate non-BASE -> leave upper
                        pass
                    elif base not in tripdeg[pos]:
                        # no BASE possibility -> obligate non-BASE -> leave upper
                        pass
                    else:
                        # non-obligate non-BASE -> lower()
                        pms[offset+pos] = pms[offset+pos].lower()
        #printsize = 99
        #for offset in range(0,seqlen,printsize):
        #    print base,"\t", "".join(pms[offset:offset+printsize])
        stats         = pmscounts("".join(pms),base)
        ratio         = float(stats[0])/(stats[0]+stats[3])
        ratioGenCode1 = float(stats[0]+stats[2])/(stats[0]+stats[3]-stats[4])
        ratioGenCode2 = ratioGenCode1 * ( float(stats[2]) / stats[0] )
        ratioGenCode3 = ( float(stats[2])/float(stats[0]) ) / ( float(stats[5]) / float(stats[3]) )
        #print stats, ratio, ratioGenCode1, ratioGenCode2, ratioGenCode3
        gcp.append( ratioGenCode3 )
    return gcp

# end of function gencodeperiodicity


def pmscounts(pms,base):
    """
    """
    alBaseCnt    = pms.count(base.upper()) + pms.count(base.lower())
    obBaseCnt    = pms.count(base.upper())   
    noBaseCnt    = pms.count(base.lower())
    alNonBaseCnt = len(pms) - alBaseCnt
    obNonBaseCnt = 0   
    noNonBaseCnt = 0
    for b in ['A','T','G','C']:
        if b != base:
            obNonBaseCnt += pms.count(b)
            noNonBaseCnt += pms.count(b.lower())
    # return the counters
    return ( alBaseCnt, obBaseCnt, noBaseCnt,
             alNonBaseCnt, obNonBaseCnt, noNonBaseCnt )
    
# end of function pmscounts


def importgeneticcode():
    """ """
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
    return standard_genetic_code

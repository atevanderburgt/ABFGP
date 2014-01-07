from lib_sequenceperiodicity import importgeneticcode
from sets import Set

def importgeneticcodedegeneracy(geneticcode=None):
    """
    """
    if not geneticcode: geneticcode = importgeneticcode()
    degeneracy = {}
    for triplet, aa in geneticcode.iteritems():
        if degeneracy.has_key(aa):
            degeneracy[aa].append(triplet)
        else:
            degeneracy[aa] = [ triplet ]
    return degeneracy

# end of function importgeneticcodedegeneracy

def importtripletbiasaward(degeneracy=None):
    """
    """
    if not degeneracy: degeneracy = importgeneticcodedegeneracy()
    tripletbiasaward = {}
    for triplets in degeneracy.values():
        for triplet in triplets:
            tripletbiasaward[triplet] = { 0: {}, 1: {}, 2: {} }
            for frame in [0,1,2]:
                bases = [ trip[frame] for trip in triplets ]
                if len(Set(bases)) == 1:
                    # no choice in base for this frame for this AA!
                    pass
                else:
                    for base in Set(bases):
                        tripletbiasaward[triplet][frame][base] = 1
    # and return this dict
    return tripletbiasaward

# end of function importtripletbiasaward

GENETICCODE   = importgeneticcode()
DEGENERACY    = importgeneticcodedegeneracy(GENETICCODE)
TRIPBIASAWARD = importtripletbiasaward(DEGENERACY)

def basebias(seq,GENETICCODE=GENETICCODE,DEGENERACY=DEGENERACY,TRIPBIASAWARD=TRIPBIASAWARD):
    """
    """
    biasperframe={
        0: { 'A': 0, 'T': 0, 'G': 0, 'C': 0 },
        1: { 'A': 0, 'T': 0, 'G': 0, 'C': 0 },
        2: { 'A': 0, 'T': 0, 'G': 0, 'C': 0 },
     }

    cnt = 0
    for triplet in [ seq[offset:offset+3].upper() for offset in range(0,len(seq),3) ]:
        try:
            aa = GENETICCODE[triplet]
            cnt+=1
        except KeyError:
            continue
        except:
            raise "UnExpected Non-KeyError while lookup in dictionary 'GENETICCODE'"
        # store the degeneracy of this triplet
        for alternative in DEGENERACY[aa]:
            if alternative == triplet:
                for frame in TRIPBIASAWARD[triplet]:
                    for base in TRIPBIASAWARD[triplet][frame].keys():
                        biasperframe[frame][base] += TRIPBIASAWARD[triplet][frame][base]
            else:
                for frame in [0,1,2]:
                    if alternative[frame] != triplet[frame]:
                        biasperframe[frame][alternative[frame]] -= 1
    # done! return the biasperframe dict
    for frame in [0,1,2]:
        for base in ['A','T','G','C']:
            biasperframe[frame][base] = float(biasperframe[frame][base]) / cnt
    return biasperframe

# end of function basebias


def basebiascheck(seq,GENETICCODE=GENETICCODE,DEGENERACY=DEGENERACY,TRIPBIASAWARD=TRIPBIASAWARD):
    """
    """
    cnt = 0
    biasperbase = {'A':[0,0],'T':[0,0],'G':[0,0],'C':[0,0]}
    for triplet in [ seq[offset:offset+3].upper() for offset in range(0,len(seq),3) ]:
        try:
            aa = GENETICCODE[triplet]
            cnt+=1
        except KeyError:
            continue
        except:
            raise "UnExpected Non-KeyError while lookup in dictionary 'GENETICCODE'"
        # store the degeneracy of this triplet

        for frame in [0,1,2]:
            base = triplet[frame].upper()
            for variablebase in TRIPBIASAWARD[triplet][frame].keys():
                if variablebase == base:
                    biasperbase[base][0]+=1
                else:
                    biasperbase[base][1]-=1
    return biasperbase

# end of function basebiascheck

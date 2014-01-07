import sys
from pssm import parse_ic_data, parse_ic_file

if len(sys.argv) == 2:
    IC = parse_ic_file(sys.argv[1])
else:
    ic_data = []
    for line in sys.stdin.readlines():
        if line[0] == '#': continue
        ic_data.append(line.strip())
    IC = parse_ic_data("\n".join(ic_data))


buffer = []    
num_seqs = 1000
for col in range(len(IC)):
    buffer.append([])
    vdict = IC[col]
    for base,value in vdict.iteritems():
        freq = pow(2,value-2.0)
        cnt = int(round(freq*num_seqs))
        buffer[-1].extend( [ base ]*cnt )
    while len(buffer[-1]) < num_seqs:
        buffer[-1].append("n")

        
seqs = {}
cnt=1
while buffer[0]:
    seq = "".join( [ buffer[col].pop() for col in range(len(IC)) ] )
    seqs[cnt] = seq
    cnt+=1

print "\n".join(seqs.values())

import sys
from pssm import obtain_pssm_ic, print_ic

seqs = {}
cnt=0
for line in sys.stdin.readlines():
    line = line.strip()
    if not line: continue
    if line[0] == ">": continue
    seqs[cnt] = line
    cnt+=1

IC, nns = obtain_pssm_ic(seqs)
print_ic(IC)

import random
import sys
import subprocess as sp

# Supply the total number of sequences for choosing a random one by index
n_seqs = int(sp.getoutput("grep -c '>' refs/HA_A_allsp_nfu_nospaces.fa"))

seen = set()
# Build commands for the case without simulated quasispecies noise
for j in [0.01, 0.02, 0.03]:
    for i in range(500):
        roll = random.randint(0, n_seqs-1)
        while roll in seen:
            roll = random.randint(0, n_seqs-1)
        seen.add(roll)
        print("python3 code/simulate_and_classify.py refs/HA_A_allsp_nfu_nospaces.fa %s %s -nq > results/%s_%s.out" % (roll, j, roll, j))

seen = set()
# Build commands for the case with simulated quasispecies noise
for j in [0.03]:
    for i in range(500):
        roll = random.randint(0, n_seqs-1)
        while roll in seen:
            roll = random.randint(0, n_seqs-1)
        seen.add(roll)
        print("python3 code/simulate_and_classify.py refs/HA_A_allsp_nfu_nospaces.fa %s %s -q > results/%s_%s_q.out" % (roll, j, roll, j))

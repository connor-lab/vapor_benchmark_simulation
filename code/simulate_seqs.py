import random

def mutate(s, p):
    mut = ""
    for k in s:
        alt = random.choice([b for b in "ATCG" if b != k])
        roll = random.uniform(0,1)
        if roll < p:
            mut += alt
        else:
            mut += k
    return mut

def generate_quasispecies(s, counts, p):
    # Generates quasispecies variants in given counts
    # Based on sequence s
    # With mutation probability p per base
    # Specified by counts
    # The original being the first count
    # Initialize with the first count, that of the major variant
    print(counts)
    population = [s for i in range(counts[0])]
    for c in counts[1:]:
        # For every other count c, mutate a new sequence
        mut = mutate(s,p)
        for i in range(c):
            population.append(mut)
    return population

def readize(seqs,L=200,p=0.0005,C=5):
    # given s, generate reads of length L
    # with error rate p
    # and coverage C
    ret = []
    for s in seqs:
        for j in range(C):
            for i in range(len(s)-L+1):
                seq = list(s[i:i+L])
                for k in range(len(seq)):
                    if random.uniform(0,1) < p:
                        alt = random.choice([b for b in "ATCG" if b != s[k]])
                        seq[k] = alt
                ret.append("".join(seq))
    return ret

def readize_random(seqs,L=200,p=0.0005,N=1000):
    # given s, generate reads of length L
    # with error rate p
    # and coverage C
    ret = []

    for i in range(N):
        rs = random.choice(seqs)
        ri = random.randint(0, len(rs)-L)
        tmps = ""
        for j in range(L):
            roll = random.uniform(0,1)
            b = rs[ri+j]
            if roll < p:
                tmps += random.choice([c for c in "ATCG" if b != c])
            else:
                tmps += b
        ret.append(tmps)
    return ret

def write_reads(reads, tag, fname):
    f = open(fname, "w")
    for ri,read in enumerate(reads):
        one = "@" + tag + "_" + str(ri) + "\n"
        two = read + "\n"
        three = "+\n"
        four = "I"*len(read) + "\n"
        f.write(one + two+three+four)
    f.close()
        
def write_reads_fasta(reads, tag, fname):
    f = open(fname, "w")
    for ri,read in enumerate(reads):
        one = ">" + tag + "_" + str(ri) + "\n"
        two = read + "\n"
        f.write(one + two)
    f.close()

if __name__ == "__main__":
	seqs = ["AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT", "AAAAAAAAAAAAAAAAAAAAAAAA", "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG"]
	print(readize_random(seqs,L=10,p=0.1,N=10))

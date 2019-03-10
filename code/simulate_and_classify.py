""" Simulates read sets from reference, then classifies with MASH, BLAST (see classification scheme), and VAPOR """
""" Firstly, takes a reference, mutates (substitutions) it to a given divergence, generates reads, adds substitutions into reads """
""" Also takes:"""
""" \t the Levenshtein distance of the original chosen sequence to its mutated counterpart """
""" \t the Levenshtein distance of each retrieved sequence to the mutated sequence """

import subprocess as sp
import sys
from Bio import SeqIO
import random
from Bio import pairwise2
from simulate_seqs import *

def levenshtein(s1, s2):
    # This will return the negative Levenshtein distance, take -1 of it
    aln = pairwise2.align.globalms(s1, s2, 0,-1,-1,-1, one_alignment_only=True)
    return aln

# Firstly load in database sequences to choose from
db = sys.argv[1]
seqs = [seq for seq in SeqIO.parse(db, "fasta")]
seqsd = {">"+seq.description : str(seq.seq) for seq in seqs}

# Take the ith sequence, and the per-base mutation probability as arguments
i = int(sys.argv[2])
mutp = float(sys.argv[3])
qflag = sys.argv[4]

# Retrieve the sequence
ref_record = seqs[i]
ref = str(seqs[i].seq)

# Set the coverage
CMP=200

# Mutate the sequence
mut = mutate(ref, mutp)
if qflag == "-q":
    # Reduce the Coverage
    CMP /= 100.
    mutref = ""
    for si, s in enumerate(generate_quasispecies(mut, [100,5,1,1,1], 0.01):)
        mutref += ">tmpseq"+str(si)+"\n"+s+"\n"
else:
    mutref = ">tmpseq\n"+mut

tmpfa = "tmp/tmp"+str(i)+"_"+qflag+".fa"
tmpfq = "tmp/tmp"+str(i)+"_"+qflag+".fq"

# Write the fasta
with open(tmpfa, "w") as f:
    f.write(mutref)

# Generate reads
sp.call("""java -jar code/afg/ArtificialFastqGenerator.jar -O tmp/%s_afg_out_0 -R %s -F1 refs/template_R1_001.fastq -F2 refs/template_R2_001.fastq -RL 151 -S '>tmpseq' -URQS true -SE true -CMP %d""" % (str(i), tmpfa, CMP), shell=True)

# Mutate the reads additionally to account for experimental error and experimental noise
sp.call("python3 code/read_mutator.py tmp/%s_afg_out_0.1.fastq > %s" % (str(i), tmpfq), shell=True)
sp.call("python3 code/read_mutator.py tmp/%s_afg_out_0.2.fastq > %s" % (str(i), tmpfq), shell=True)

n_reads = int(sp.getoutput("wc -l %s" % tmpfq).split()[0])/4.0

# Now classify the reads with MASH
sp.call("mash sketch -r %s" %tmpfq, shell=True)
msh_out = [l.split("\t") for l in sp.check_output("mash dist %s %s" % (db+".msh", tmpfq), shell=True).decode("utf-8").split("\n") if l != '']
msh_out = sorted(msh_out, key=lambda x:float(x[2]))
min_msh_dist = float(min(msh_out, key = lambda x:float(x[2]))[2])
msh_choices = [msh_out[i][0] for i in range(len(msh_out)) if float(msh_out[i][2]) == min_msh_dist]

# Now classify the reads with BLAST
read_tmpfa = tmpfq.split(".")[0]+".reads.fa"
sp.call("python3 code/fastq2fa.py %s > %s" % (tmpfq, read_tmpfa), shell=True)
blast_choices = [l.strip() for l in sp.check_output("python3 code/blast_classify.py %s %s" % (read_tmpfa, db), shell=True).decode("utf-8").split("\n") if l != '']

# Now classify the reads with VAPOR
sys.stderr.write("vapor.py -q -fa %s -fq %s" % (db, tmpfq))
vapor_raw_out = sp.getoutput("vapor.py -q -fa %s -fq %s" % (db, tmpfq))
sys.stderr.write("??"+vapor_raw_out + "/\n")
vapor_out = [l for l in vapor_raw_out.split("\n") if l != ''][0]
sys.stderr.write(vapor_raw_out + "\n")
vapor_choices = [vapor_out.split("\t")[5]]

sys.stderr.write(str(n_reads)+","+str(len(vapor_choices))+str(len(blast_choices))+"\n")

# Calculate Levenshtein distances
mashlevs = []
mashpids = []
vaporlevs = []
vaporpids = []
blastlevs = []
blastpids = []
for mc in msh_choices:
    alnm = levenshtein(mut, seqsd[">"+mc])[0]
    levm = -1*alnm[2]
    pidm = 1-levm/len(alnm[0])
    mashlevs.append(str(levm))
    mashpids.append(str(pidm))
for vc in vapor_choices:
    alnv = levenshtein(mut, seqsd[vc])[0]
    levv = -1*alnv[2]
    pidv = 1-levv/len(alnv[0])
    vaporlevs.append(str(levv))
    vaporpids.append(str(pidv))
for bc in blast_choices:
    alnb = levenshtein(mut, seqsd[">"+bc])[0]
    levb = -1*alnb[2]
    pidb = 1-levb/len(alnb[0])
    blastlevs.append(str(levb))
    blastpids.append(str(pidb))

# Get the true Levenshtein distance of the major variant (if quasi, if not, the only sequence)
alnt = levenshtein(mut, str(seqs[i].seq))[0]
levt = -1*alnt[2]
pidt = 1-levt/len(alnt[0])
true_score = str(levt)
true_pid = str(pidt)

msh_choices_str = ",".join(msh_choices)
mashlevsstr = ",".join([str(m) for m in mashlevs])
mashpidsstr = ",".join([str(m) for m in mashpids])
vapor_choices_str = ",".join(vapor_choices)
vaporlevsstr = ",".join([str(m) for m in vaporlevs])
vaporpidsstr = ",".join([str(m) for m in vaporpids])
blast_choices_str = ",".join(blast_choices)
blastlevsstr = ",".join([str(m) for m in blastlevs])
blastpidsstr = ",".join([str(m) for m in blastpids])
covstr = str(151*n_reads/float(len(mut)))

# Output results
sys.stderr.write(str(i) + " "+ ref_record.description + " "+ str(n_reads) + " " + str(len(mut)) + " " + covstr + " " + true_score + " " + true_pid + " " + msh_choices_str + " " + mashlevsstr + " " + mashpidsstr + " " + vapor_choices_str + " " + vaporlevsstr + " " + vaporpidsstr + " " + blast_choices_str + " " + blastlevsstr + " " + blastpidsstr + "\n")
print(str(i) + " "+ ref_record.description + " "+ str(n_reads) + " " + str(len(mut)) + " " + covstr + " " + true_score + " " + true_pid + " " + msh_choices_str + " " + mashlevsstr + " " + mashpidsstr + " " + vapor_choices_str + " " + vaporlevsstr + " " + vaporpidsstr + " " + blast_choices_str + " " + blastlevsstr + " " + blastpidsstr)

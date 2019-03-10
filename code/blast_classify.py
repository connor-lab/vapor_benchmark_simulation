import sys
import subprocess as sp
import itertools as it

""" Classifies a set of short reads, in .fa format, by taking the best hit for each read """
""" Best is defined as that with the lowest e-value, followed by the highest bit score """

readres = {}
db, query = sys.argv[1:]
# Perform BLAST
blstr = "blastn -query %s -db %s -outfmt '6 qseqid sseqid evalue score'" % (db, query)
sys.stderr.write(blstr+"\n")
blastres = [l.strip() for l in sp.check_output(blstr, shell=True).decode("utf-8").split("\n") if l !=  '']
# Parse the results on a per-read basis
for l in blastres:
	spl = l.split("\t")
	read, ref, evalue, score = spl
	evalue = float(evalue)
	score = float(score)
	if read in readres:
		readres[read].append([ref, evalue, score])
	else:
		readres[read] = [[ref,evalue,score]]
sys.stderr.write("Got results\n")

# Now get the best ref for each read, store counts in a dictionary
maxrefs = {}
for read, res in readres.items():
        # For each read, get the result with both 1. the minimum evalue
        # And in any ties, the maximum bit score
        # Achieved by sorting the tuple; here negative bit score is taken since 
        # The minimum is taken
	maxres = min(res, key = lambda x: (x[1], -x[2]))[1:]
	for q in res:
		ref = q[0]
		if q[1:] == maxres:
			if ref in maxrefs:
				maxrefs[ref] += 1
			else:
				maxrefs[ref] = 1

# Finally, take the reference(s) that are selected the greatest number of times
maxofmax = max(maxrefs.items(), key = lambda x:x[1])[1]
for i,q in maxrefs.items():
	if q == maxofmax:
		print(i)

			

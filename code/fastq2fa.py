import sys
from Bio import SeqIO

for r in SeqIO.parse(sys.argv[1], "fastq"):
	print(">"+r.description)
	print(str(r.seq))

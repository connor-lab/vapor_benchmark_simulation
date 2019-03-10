import random
import sys

def mutate_read(r, p):
	read = list(r)
	for i in range(len(read)):
		if random.uniform(0,1) < p:
			alt = random.choice([b for b in "ATCG" if read[i] != b])
			read[i] = alt
	return "".join(read)
	
c = 0
with open(sys.argv[1]) as f:
	for l  in f:
		if c == 1:
			print(mutate_read(l.strip(), 0.01))
		else:
			print(l.strip())
		c += 1
		if c > 3:
			c = 0
		
		

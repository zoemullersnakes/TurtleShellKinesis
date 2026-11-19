# python3 get_chrom_lengths.py genome.fasta

import sys

def getfasta(fastahandle):
	fastadata = {}
	name = None
	sequence = ''
	
	for line in fastahandle:
		
		if not line.startswith(">"):		# this is continuing a sequence
			line = line.rstrip()
			sequence += line
		
		elif line.startswith(">"):			# start of a new sequence
			if name is None:				# the first sequence in the file!
				line = line.rstrip()
				name = line[1:len(line)]
			else:							# report the previous name and sequence, then refresh the name and sequence
				fastadata[name] = sequence
				line = line.rstrip()
				name = line[1:len(line)]
				sequence = ''				# restart the sequence
				
	# report the last name and sequence
	fastadata[name] = sequence
	return(fastadata)

genome = getfasta(open(sys.argv[1]))

genome_size = 0

genome_size_file = sys.argv[1] + ".sizes"
SIZES = open(genome_size_file, "w")

for name in genome.keys():
	SIZES.write(name + "\t" + str(len(genome[name])) + "\n")



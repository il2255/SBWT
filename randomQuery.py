#! /usr/bin/python

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from random import randint, random, choice
import os, pprint, re, sys, getopt, const

def cleanFiles():
	print 'Clean all query generation related files'
	prefixToDelete = [const.QUERY_ORI_PREFIX, const.QUERY_MUT_PREFIX]
	path = "./"
	for f in os.listdir(path):
		fileName = os.path.splitext(f)[0]
		for prefix in prefixToDelete:
			if fileName.startswith(prefix):
				os.remove(f)

def mutateSequence(original_sequence, queryStartPos, mutationRate_1, mutationRate_2, mutationRate_3):
	# https://www.biostars.org/p/12417/
	count = 0
	sequence = list(original_sequence)
	for i, s in enumerate(sequence):
		curr_pos = i + queryStartPos
		val = random()
		if curr_pos % 3 == 0:
			MutationRate = mutationRate_1
		elif curr_pos % 3 == 1:
			MutationRate = mutationRate_2
		else:
			MutationRate = mutationRate_3
		if val < MutationRate:
			count += 1
			sequence[i] = choice([x for x in "actg" if x != s.lower()])

	# Forcefully generate mutation if none is found.  This alters the statistics thus turn it off.
	# But in the future this might be helpful
	# if count == 0:
	# 	# Randomly select at least one mutation
	# 	mutationPoint = randint(0,len(sequence)-1)
	# 	sequence[mutationPoint] = choice([x for x in "actg" if x != sequence[mutationPoint].lower()])
	mutatedSequence = "".join(sequence)
	return mutatedSequence	

def main():
	try:
		opts, args = getopt.getopt(sys.argv[1:], "1:2:3:n:l:i:ch", ["help","clean"])
	except getopt.GetoptError as err:
		# print help information and exit:
		print str(err) # will print something like "option -a not recognized"
		usage()
		sys.exit(2)

	# Initialize variable
	mut_rate_1 = const.MUTATION_RATE_1
	mut_rate_2 = const.MUTATION_RATE_2
	mut_rate_3 = const.MUTATION_RATE_3
	num_sequence = const.RANDOM_SEQ_NUMBER
	len_sequence = const.RANDOM_SEQ_LENGTH
	len_codon = const.RANDOM_SEQ_CODON_LENGTH
	inputFileName = None

	# Retrieve corresponding command line values
	for o, a in opts:
		if o == "-1":
			if len(a) > 0: 
				mut_rate_1 = float(a)
			else: 
				assert False, "No first codon mutation rate found"
		elif o in "-2":
			if len(a) > 0: 
				mut_rate_2 = float(a)
			else: 
				assert False, "No second codon mutation rate found"
		elif o in "-3":
			if len(a) > 0: 
				mut_rate_3 = float(a)
			else: 
				assert False, "No third codon mutation rate found"
		elif o in "-n":
			if len(a) > 0: 
				num_sequence  = int(a)
			else: 
				assert False, "No number of random sequence found"     
		elif o in "-l":
			if len(a) > 0: 
				len_codon  = int(a)
				len_sequence = len_codon * 3
			else: 
				assert False, "No codon length found"  
		elif o == "-i":
			if len(a) > 0: 
				inputFileName = a
			else: 
				assert False, "no input file name found"
		elif o in ["--clean", "-c"]:
			cleanFiles()
			sys.exit(2)    
		elif o in ["--help", "-h"]:
			usage()
			sys.exit(2)    
		else:
			assert False, "unhandled option"

	# Check if input file is found
	if inputFileName is None:
		print 'No Reference file sequence is found'
		usage()
		sys.exit(2)
     
	# Open the reference file
	input_f = open(inputFileName, 'rb')
	refFileName = os.path.splitext(inputFileName)[0]
	ref_sequences = list(SeqIO.parse(input_f,'fasta'))
	input_f.close()

	# Filter out the sequence which has less than the desired length
	sequences = filter(lambda k: len(k.seq) >= len_sequence, ref_sequences)
	if len(sequences) == 0:
		print 'References do not have any sequence greater than random sequence length'
		sys.exit(2)
	
	# Open files to write randomly selected segment
	outputName_Ori = const.QUERY_ORI_PREFIX + refFileName + const.FASTA_SUFFIX
	output_f_original = open(outputName_Ori, 'wb')	

	# Open files to write mutated randomly selected segment
	outputName_Mut = const.QUERY_MUT_PREFIX + refFileName + const.FASTA_SUFFIX
	output_f_mutated = open(outputName_Mut, 'wb')	

	# Create random sequences
	for i in xrange(num_sequence):
		# Pick a random sequence
		singleSequence = choice(sequences)
		refSeqName = singleSequence.id
		refSeq = singleSequence.seq
		refSeqLength = len(refSeq)

		# Get the random starting point.
		randEnd = refSeqLength - len_sequence
		start = randint(0,randEnd)
		end = start + len_sequence

		# Randomly selected sequence
		random_seq = str(refSeq[start:end])

		# Mutate the sequence according to the probablity
		mutatedSeq = mutateSequence(random_seq, start, mut_rate_1, mut_rate_2, mut_rate_3)

		# Generate reverse complement randomly as well
		val = random()
		if val <= const.REVERSE_COMP_RATE:
			random_seq = str(Seq(random_seq).reverse_complement())
			mutatedSeq = str(Seq(mutatedSeq).reverse_complement())
			start = start * -1

		# Generate the ID string for the fasta format
		id_string = str(refSeqName) + const.ID_DELIMIETER + str(i) + const.ID_DELIMIETER + str(start) + const.ID_DELIMIETER + str(0)
		# Print to the output original file. 
		original_seq = SeqRecord(Seq(random_seq),id= id_string,name='',description='START:'+str(start))
		SeqIO.write(original_seq, output_f_original, 'fasta')
		# Print to the output of mutated file
		mutated_seq = SeqRecord(Seq(mutatedSeq),id=id_string,name='',description='START:'+str(start))
		SeqIO.write(mutated_seq,output_f_mutated, 'fasta')		
		if i % 1000 == 0:
			print str(i) + ' sequences are generated'	
	# Close the file and exit
	print 'Random sequence generation done'
	output_f_original.close()
	output_f_mutated.close()

# Usage information. 
def usage():
	print "Usage:"
	print " -i <input file name for reference sequence>"
	print "\t The input file name in fasta format which contains the reference sequence"
	print " -1 <mutation rate for first codon position>"
	print "\t The values between 0 and 1. Default is set to 0.05>"
	print " -2 <mutation rate for second codon position>"
	print "\t The values between 0 and 1. Default is set to 0.05>"
	print " -3 <mutation rate for third codon position>"
	print "\t The values between 0 and 1. Default is set to 0.05>"
	print " -n <number of random sequence to generate.>"
	print "\t The number of random sequence to generate. Default is set to 10000"
	print " -l <number of codons for each random sequence.>"
	print "\t This is number of codons not nucleotides. Default is set to 20 and its nucleotides length will be 60"
	print "--clean | -c"
	print "\t Remove all random query generation related files"
	print "--help | -h"
	print "\t Help text which is what you see :)"

if __name__ == '__main__':
	main()
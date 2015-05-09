#! /usr/bin/python

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from random import randint, random, choice
import os, pprint, re, sys, getopt, const

def cleanFiles():
	print 'Clean all cleaned reference file'
	prefixToDelete = [const.REF_PREFIX_CLEAN]
	path = "./"
	for f in os.listdir(path):
		fileName = os.path.splitext(f)[0]
		for prefix in prefixToDelete:
			if fileName.startswith(prefix):
				os.remove(f)

def main():
	try:
		opts, args = getopt.getopt(sys.argv[1:], "i:ch", ["help", "clean"])
	except getopt.GetoptError as err:
		# print help information and exit:
		print str(err) # will print something like "option -a not recognized"
		usage()
		sys.exit(2)

	# Initialize variable
	inputFileName = None

	# Retrieve corresponding command line values
	for o, a in opts:
		if o == "-i":
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

	# Open the file using the SeqIO
	input_f = open(inputFileName, 'rb')
	refFileName = os.path.splitext(inputFileName)[0]
	outputName = const.REF_PREFIX_CLEAN + refFileName + const.FASTA_SUFFIX
	output_f = open(outputName, 'wb')
	for singleSequence in SeqIO.parse(input_f,'fasta'):
		# Retrieve the reference sequence.
		refSeqName = singleSequence.id
		refSeq = singleSequence.seq
		refSeqLength = len(refSeq)
		# Make a mutable sequence and remove repeats
		cleaned_Ref = re.sub('[^ACGTacgt]', '', str(refSeq))
		# Print the second only codon information from reference sequence	
		refSeq_clean_record = SeqRecord(Seq(cleaned_Ref),id=const.REF_PREFIX_CLEAN+refSeqName,name=refSeqName,description='This is  Cleaned Reference Seq')
		SeqIO.write(refSeq_clean_record,output_f, 'fasta')
	input_f.close()
	output_f.close()

# Usage information.
def usage():
	print "Usage:"
	print " -i <input file name for reference sequence>"
	print "\t The input file name in fasta format which contains the reference sequence to be cleaned"
	print "\t This will remove mutable sequence and remove repeat"
	print "\t Assumes lower character as mutable sequence such as a,c,g,t"
	print "\t Assumes N as the repeat sequence"
	print "--clean | -c"
	print "\t Remove all cleaned reference files"
	print "--help | -h"
	print "\t Help text which is what you see :)"

if __name__ == '__main__':
	main()
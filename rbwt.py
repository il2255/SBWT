#! /usr/bin/python
import getopt, sys, os, subprocess, const, pysam

def cleanFiles(excludeSAM = False):
	suffixToDelete = ['.sai', '.amb', '.ann', '.bwt', '.pac', '.sa', '.sam']
	if excludeSAM:
		del suffixToDelete[-1]
	path = "./"
	for f in os.listdir(path):
		suffix = os.path.splitext(f)
		if suffix[1] in suffixToDelete:
			os.remove(f)

def main():
	try:
		opts, args = getopt.getopt(sys.argv[1:], "r:q:hc", ["help", "clean"])
	except getopt.GetoptError as err:
		# print help information and exit:
		print str(err) # will print something like "option -a not recognized"
		usage()
		sys.exit(2)

	referenceFile = None
	queryFile = None

	# Retrieve corresponding command line values
	for o, a in opts:
		if o == "-r":
			if len(a) > 0: 
				referenceFile = a
			else: 
				assert False, "no query text found"
		elif o in "-q":
			if len(a) > 0: 
				queryFile = a
			else: 
				assert False, "no query text found"
		elif o in ["--clean", "-c"]:
			cleanFiles(False)
			sys.exit(2)
		elif o in ["--help", "-h"]:
			usage()
			sys.exit(2)        
		else:
			assert False, "unhandled option"

	# Initialize classes	
	if referenceFile == None or queryFile == None:
		usage()
		sys.exit(2)

	# Build the index from the reference file. 
	print '1. BUILD THE INDEX FOR THE REFERENCE SEQUENCE'
	index_command = const.COMMAND_BUILD_INDEX + referenceFile
	subprocess.call(index_command, shell=True)

	# Run a single end alignment usnig aln
	print '\n\n2. RUN THE BWA ON THE QUERY TO FIND THE ALIGNMENT'
	queryFileName = os.path.splitext(queryFile)[0]
	aln_outputFileName = const.ALN_OUTPUT_PREFIX + queryFileName + const.ALN_OUTPUT_SUFFIX
	query_command = const.COMMAND_FIND_ALN + aln_outputFileName + ' ' + referenceFile + ' ' + queryFile
	subprocess.call(query_command, shell=True)

	# Convert the aln to the sam using single end.
	print '\n\n3. CONVERT THE BWT RESULT TO THE SAM FORMAT'
	sam_outputFileName = const.SAM_OUTPUT_PREFIX + queryFileName + const.SAM_OUTPUT_SUFFIX
	sam_command = const.COMMAND_SAM_FORMAT + sam_outputFileName + ' ' + referenceFile + ' ' +  aln_outputFileName + " " + queryFile
	subprocess.call(sam_command, shell=True)
	
	print '\n\n4. END OF THE PROGRAM'
	

# Usage information. Copied from the reference version.
def usage():
	print "Usage:"
	print " -r <reference file>"
	print "\t Reference file sequence in fasta format"
	print " -q <query file>"
	print "\t Query file sequence in fasta format"
	print "--clean | -c"
	print "\t Remove all BWT realted files"
	print "--help | -h"
	print "\t Help text which is what you see :)"

if __name__ == '__main__':
	main()

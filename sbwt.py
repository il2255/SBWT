#! /usr/bin/python
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from bitarray import bitarray
from bitarray import bitdiff
import os, pprint, re, sys, getopt, const, subprocess, pysam, checkResult

DEBUG = 0

def cleanFiles(excludeSAM = False):
	suffixToDelete = ['.sai', '.amb', '.ann', '.bwt', '.pac', '.sa', '.sam']
	if excludeSAM:
		del suffixToDelete[-1]
	path = "./"
	for f in os.listdir(path):
		suffix = os.path.splitext(f)
		if suffix[1] in suffixToDelete:
			os.remove(f)
	prefixToDelete = [const.REF_PREFIX_SEC, const.QUERY_SEC_PREFIX]
	path = "./"
	for f in os.listdir(path):
		fileName = os.path.splitext(f)[0]
		for prefix in prefixToDelete:
			if fileName.startswith(prefix):
				os.remove(f)

def referencePos(codonPosition, reverse, start_pos):
	# Get the first codon position query sequence
	if codonPosition == 3 and reverse:
		reference_start = (abs(start_pos) + 1) 
	elif codonPosition == 1 and not reverse:
		reference_start = abs(start_pos) + 1
	else:
		reference_start = abs(start_pos)

	if reverse:
		reference_start = reference_start * -1
	return reference_start	

def checkFirstCodonPos(refSequence, refSeqPosition, querySequence, queryLength, mask_array, threshold, reverse):
	if DEBUG: print "Ref Pos   :", refSeqPosition
	#if DEBUG: print "Reference :", refSequence
	if DEBUG: print "Query1    :", querySequence

	if reverse:
		# Get the reverse complementary
		querySequence = str(Seq(querySequence).reverse_complement())
		refSeqPosition = abs(refSeqPosition)
	# Convert to bitarray
	querySequenceBitarray = bitarray()
	querySequenceBitarray.encode(const.BIT_DICT, querySequence)
	
	bitStart = abs(refSeqPosition * 4)
	bitEnd = bitStart + queryLength * 4
	
	refSeqPosition = refSeqPosition
	refBitArray = refSequence[0][refSeqPosition*4:(refSeqPosition+queryLength)*4]
	if DEBUG: print "Query2    :", querySequence
	if DEBUG: print "Ref Corres:", refSequence[1][refSeqPosition:refSeqPosition+queryLength]
	if DEBUG: print "Ref Bit   :", refBitArray
	if DEBUG: print "Query Bit :", querySequenceBitarray
	try:
		diffPercent = float(queryLength - bitdiff(refBitArray, querySequenceBitarray)/2)/queryLength
		if DEBUG: print "Diff Bit  :", diffPercent, bitdiff(refBitArray, querySequenceBitarray)/2
		if diffPercent >= threshold:
			return True
		else:
			return False
	except:
		return False

def main():
	try:
		opts, args = getopt.getopt(sys.argv[1:], "r:q:1t:ch", ["help","clean"])
	except getopt.GetoptError as err:
		# print help information and exit:
		print str(err) # will print something like "option -a not recognized"
		usage()
		sys.exit(2)

	# Initialize variable
	orgRefFileName = None
	orgQueryFileName = None
	checkFirstCodon = const.CHECK_FIRST_CODON
	firstCodonThreshold = const.FIRST_CODON_THRESHOLD

	# Retrieve corresponding command line values
	for o, a in opts:
		if o == "-1":
			checkFirstCodon = True
		elif o in "-t":
			if len(a) > 0: 
				firstCodonThreshold = float(a)
			else: 
				assert False, "No first codon accuracy threshold found"
		elif o == "-r":
			if len(a) > 0: 
				orgRefFileName = a
			else: 
				assert False, "no input file name found"
		elif o == "-q":
			if len(a) > 0: 
				orgQueryFileName = a
			else: 
				assert False, "no query file name found"
		elif o in ["--clean", "-c"]:
			cleanFiles(False)
			sys.exit(2)    
		elif o in ["--help", "-h"]:
			usage()
			sys.exit(2)    
		else:
			assert False, "unhandled option"

	# Check if input file is found
	if orgRefFileName is None:
		print 'No reference sequence file is found'
		usage()
		sys.exit(2)

	if orgQueryFileName is None:
		print 'No query sequence file is found'
		usage()
		sys.exit(2)

	# Open the reference file and print 2nd Position only Reference or BitArray for 1st Codon Position.
	print '\n\n1. 2ND CODON POSITION OF THE REFERENCE SEQUENCE'
	input_f = open(orgRefFileName, 'rb')
	refFileName = os.path.splitext(orgRefFileName)[0]
	outputName = const.REF_PREFIX_SEC + refFileName + const.FASTA_SUFFIX
	referenceFile = outputName
	output_f = open(outputName, 'wb')

	# If first codon position check is on, create a dictionry of bitarry of first codon position.
	if checkFirstCodon:
		firstCodonDict = dict()
	for singleSequence in SeqIO.parse(input_f, 'fasta'):
		refSeqName = singleSequence.id
		refSeq = singleSequence.seq
		refSeqLength = len(refSeq)
		# Get second codon position only
		refSeq_second = (refSeq[1:])[::3]		
		refSeq_second_record = SeqRecord(Seq(str(refSeq_second)),id=refSeqName,name=refSeqName,description='This is Reference Seq with 2nd Pos')
		SeqIO.write(refSeq_second_record,output_f, 'fasta')
		# Get first codon position only
		if checkFirstCodon:
			refSeq_first = str(refSeq[::3])
			refSeq_first_bitarray = bitarray()
			refSeq_first_bitarray.encode(const.BIT_DICT, refSeq_first)
			firstCodonDict[str(refSeqName)] = [refSeq_first_bitarray, refSeq_first]
	input_f.close()
	output_f.close()
	
	print '\n\n2. 2ND CODON POSITION OF THE QUERY SEQUENCE'
	# Open the query file and generate the query in 3 directions
	input_f = open(orgQueryFileName, 'rb')
	queryFileName = os.path.splitext(orgQueryFileName)[0]
	# Open file for the output for query with codon based
	outputName_2nd = const.QUERY_SEC_PREFIX + queryFileName + const.FASTA_SUFFIX
	queryFile = outputName_2nd
	output_f_query = open(outputName_2nd, 'wb')
	firstCodonDictQuery = dict()
	for singleSequence in SeqIO.parse(input_f, 'fasta'):
		refSeqName = singleSequence.id
		refSeq = singleSequence.seq
		refDesc = singleSequence.description
		refSeqLength = len(refSeq)
		# Generate three shift queries
		first_seq = str(refSeq[::3])
		second_seq = str((refSeq[1:])[::3])
		third_seq = str((refSeq[2:])[::3])		
		# Generate three Seq Record
		query_first = SeqRecord(Seq(first_seq),id=refSeqName[:-1] + '1',name='',description=refDesc)
		query_second = SeqRecord(Seq(second_seq),id=refSeqName[:-1] + '2',name='',description=refDesc)
		query_third = SeqRecord(Seq(third_seq),id=refSeqName[:-1] + '3',name='',description=refDesc)
		# Write to the query file
		SeqIO.write(query_first,output_f_query, 'fasta')
		SeqIO.write(query_second,output_f_query, 'fasta')
		SeqIO.write(query_third,output_f_query, 'fasta')
		# Generate First Codon Check dictionary
		if checkFirstCodon:
			firstCodonDictQuery[refSeqName[:-1]] = [first_seq, second_seq, third_seq]
	output_f_query.close()
	input_f.close()	

	# Get the seq length if it is checking first codon
	if checkFirstCodon:
		query_length = len(first_seq)
		second_pos_mask = bitarray('01'*query_length)

	# Now run the regular bwt with the generated files. 
	# Build the index from the reference file. 
	print '\n\n3. BUILD THE INDEX FOR THE REFERENCE SEQUENCE'
	index_command = const.COMMAND_BUILD_INDEX + referenceFile
	subprocess.call(index_command, shell=True)

	# Run a single end alignment usnig aln
	print '\n\n4. RUN THE BWA ON THE QUERY TO FIND THE ALIGNMENT'
	queryFileName = os.path.splitext(queryFile)[0]
	aln_outputFileName = const.ALN_OUTPUT_PREFIX + queryFileName + const.ALN_OUTPUT_SUFFIX
	query_command = const.COMMAND_FIND_ALN + aln_outputFileName + ' ' + referenceFile + ' ' + queryFile
	subprocess.call(query_command, shell=True)

	# Convert the aln to the sam using single end.
	print '\n\n5. CONVERT THE BWT RESULT TO THE SAM FORMAT'
	sam_outputFileName = const.SAM_OUTPUT_PREFIX + queryFileName + const.SAM_OUTPUT_SUFFIX
	sam_command = const.COMMAND_SAM_FORMAT + sam_outputFileName + ' ' + referenceFile + ' ' +  aln_outputFileName + " " + queryFile
	subprocess.call(sam_command, shell=True)


	print '\n\n6. CHECK FIRST CODON POSITION'
	truePositive = 0
	falsePositive = 0
	truePositiveMultiple = 0
	falsePositiveMultiple = 0
	notMapped = 0
	# Check the first codon position
	if checkFirstCodon:
		samfile = pysam.AlignmentFile(sam_outputFileName, 'r')
		iterator = samfile.fetch()

		# Initialize the variable for three positions.
		count = 0
		aln_pos = []
		ref_pos = None

		for x in iterator:
			count = count + 1			
			seqID = x.query_name
			ID_words = seqID.split(const.ID_DELIMIETER)
			refSeqName = seqID[:-1]
			codon_position = int(seqID[-1])
			reference_start = (abs(x.reference_start))
			ref_pos = int(ID_words[2])
			offset = int(ID_words[3])
			ref_pos_ans = tuple([ID_words[0],ref_pos])

			# Get the first codon position query sequence
			if codon_position == 3:
				first_codon_seq_reverse = firstCodonDictQuery[refSeqName][0]
				first_codon_seq_forward = firstCodonDictQuery[refSeqName][1]
			elif codon_position == 2:
				first_codon_seq_reverse = firstCodonDictQuery[refSeqName][2]
				first_codon_seq_forward = firstCodonDictQuery[refSeqName][0]					
			else:
				first_codon_seq_reverse = firstCodonDictQuery[refSeqName][1]
				first_codon_seq_forward = firstCodonDictQuery[refSeqName][2]

			reference_start = referencePos(codon_position, x.is_reverse, reference_start)

			# Check for the first codon position match
			if x.is_reverse:
				first_codon_seq = first_codon_seq_reverse
			else:
				first_codon_seq = first_codon_seq_forward

			if not x.is_unmapped:					
				if checkFirstCodonPos(firstCodonDict[samfile.getrname(x.reference_id)], reference_start, first_codon_seq, query_length, second_pos_mask, firstCodonThreshold, x.is_reverse):
					# Get the alignment position.
					temp_aln_pos = checkResult.getAlignedPosition(x.reference_start + 1, 1, offset)				
					if x.is_reverse:
						temp_aln_pos = (temp_aln_pos + offset - 4 + offset) * -1
					aln_pos.append(tuple([samfile.getrname(x.reference_id),temp_aln_pos]))
					
			# Check the multiple position check
			if x.has_tag('XA'):
				XA_TAG = x.get_tag('XA')
				tags = XA_TAG.split(const.TAG_DELIMIETER)
				for singleTag in tags:
					tag_words = singleTag.split(',')
					if len(tag_words) > 1:
						sam_seq_id = tag_words[0]
						sam_pos_val = int(tag_words[1])
						# Check if reverse
						reverseAln = (sam_pos_val < 0)
						reference_start = referencePos(codon_position, reverseAln, abs(sam_pos_val) - 1)
						# If reverse get differnet first codon position
						if reverseAln:
							first_codon_seq = first_codon_seq_reverse
						else:
							first_codon_seq = first_codon_seq_forward							
						# If first codon position still matches quite high append to the possible solution
						if checkFirstCodonPos(firstCodonDict[sam_seq_id], reference_start, first_codon_seq, query_length, second_pos_mask, firstCodonThreshold, reverseAln):
							temp_aln_pos = checkResult.getAlignedPosition(sam_pos_val, 1, offset) 
							aln_pos.append(tuple([sam_seq_id,temp_aln_pos]))
								
			if count%3 == 0:
				# Check if the reference is in the one of the answer.
				if len(aln_pos) == 0:
					notMapped = notMapped + 1
				elif ref_pos_ans in aln_pos:
					truePositive = truePositive + 1
					if len(aln_pos) > 1:
						truePositiveMultiple = truePositiveMultiple + 1
				else:
					falsePositive = falsePositive + 1
					if len(aln_pos) > 1:
						falsePositiveMultiple = falsePositiveMultiple + 1
				aln_pos = []
				ref_pos = None

		total = notMapped + truePositive + falsePositive
		print 'True Positive    					 	:',truePositive
		print 'True Positive with multiple alignment  				:',truePositiveMultiple
		print 'False Positive   						:',falsePositive
		print 'False Positive with multiple alignment 				:',falsePositiveMultiple
		print 'Not Mapped       						:', notMapped
		print 'Total            						:', total		
	# End of the program
	print '\n\nEND OF THE PROGRAM'
	return
	
# Usage information. 
def usage():
	print "Usage:"
	print " -r <input file name for reference sequence>"
	print "\t The input file name in fasta format which contains the reference sequence"
	print " -q <input file name for query sequence>"
	print "\t The input file name in fasta format which contains the query sequences"
	print " -1"
	print "\t Set this flag to turn on the first codon position match. It is off by default"
	print " -t <threshold for the 1st codon position match>"
	print "\t The values between 0 and 1. Default is set to 1 to match all>"
	print "--clean | -c"
	print "\t Remove all random query generation related files"
	print "--help | -h"
	print "\t Help text which is what you see :)"

if __name__ == '__main__':
	main()
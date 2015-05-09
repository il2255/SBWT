import pysam
import os, pprint, re, sys, getopt, const


SAM_TYPE = 0

def cleanFiles():
	print 'Clean all query generation related files'
	prefixToDelete = ['.sam']
	path = "./"
	for f in os.listdir(path):
		fileName = os.path.splitext(f)[0]
		for prefix in prefixToDelete:
			if fileName.startswith(prefix):
				os.remove(f)

def getAlignedPositions(XA_TAGS, sam_type, offset):
	als_positions = []
	tags = XA_TAGS.split(const.TAG_DELIMIETER)
	for singleTag in tags:
		tag_words = singleTag.split(',')
		if len(tag_words) > 1:
			sam_seq_id = tag_words[0]
			sam_pos_val = int(tag_words[1])
			aln_pos = getAlignedPosition(sam_pos_val, sam_type, offset) 
			als_positions.append(tuple([sam_seq_id,aln_pos]))
	return als_positions
	
def getAlignedPosition(start_position, sam_type, offset):
	if start_position < 0:
		offset = 4 - offset
		start_position = abs(start_position) -1
		if sam_type == 1:
			ref_pos = start_position * 3 + 1 - offset + 1
		else:
			ref_pos = start_position
		ref_pos = ref_pos * -1		
	else:
		start_position = abs(start_position) -1
		if sam_type == 1:
			ref_pos = start_position * 3 + 1 - offset + 1
		else:
			ref_pos = start_position
	return ref_pos

def main():
	try:
		opts, args = getopt.getopt(sys.argv[1:], "i:t:hc", ["help", "clean"])
	except getopt.GetoptError as err:
		# print help information and exit:
		print str(err) # will print something like "option -a not recognized"
		usage()
		sys.exit(2)

	sam_type = None
	sam_file = None

	# Retrieve corresponding command line values
	for o, a in opts:
		if o == "-i":
			if len(a) > 0: 
				sam_file = a
			else: 
				assert False, "no query text found"
		elif o in "-t":
			if len(a) > 0: 
				sam_type = int(a)
			else: 
				assert False, "no query text found"
		elif o in ["--clean", "-c"]:
			cleanFiles()
			sys.exit(2)
		elif o in ["--help", "-h"]:
			usage()
			sys.exit(2)        
		else:
			assert False, "unhandled option"

	# Initialize classes	
	if sam_type == None or sam_file == None:
		print 'Missing SAM Type or SAM file'
		usage()
		sys.exit(2)

	samfile = pysam.AlignmentFile(sam_file, 'r')
	iterator = samfile.fetch()

	truePositive = 0
	falsePositive = 0
	truePositiveMultiple = 0
	falsePositiveMultiple = 0
	notMapped = 0
	if sam_type == 1:
		withinLoop = 3
	else:
		withinLoop = 1
	count = 0
	aln_pos = []
	ref_pos = None
	# Collect answers for a single sequence for regular but three for second codon position
	for x in iterator:		
		count = count + 1
		seqID = x.query_name
		ID_words = seqID.split(const.ID_DELIMIETER)
		offset = int(ID_words[3])
		# Get the reference position (Original position)
		ref_pos = int(ID_words[2])
		# Get the reference solution
		ref_pos_ans = tuple([ID_words[0],ref_pos])		

		# Skip if it is unmapped
		if not x.is_unmapped:
			# Get the alignment position.
			temp_aln_pos = getAlignedPosition(x.reference_start + 1, sam_type, offset)				
			if x.is_reverse:
				if sam_type == 1:
					temp_aln_pos = (temp_aln_pos + offset - 4 + offset) * -1
				else:
					temp_aln_pos = temp_aln_pos * -1
			aln_pos.append(tuple([samfile.getrname(x.reference_id),temp_aln_pos]))

		# Check the multiple position check
		if x.has_tag('XA'):
			XA_TAG = x.get_tag('XA')
			ID_words = x.query_name.split(const.ID_DELIMIETER)
			offset = int(ID_words[3])			
			aln_pos.extend(getAlignedPositions(XA_TAG, sam_type, offset))

		if count%withinLoop == 0:
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
			seqID = None
			aln_pos = []
			ref_pos = None

	total = notMapped + truePositive + falsePositive
	print 'True Positive    					 	:',truePositive
	print 'True Positive with multiple alignment  				:',truePositiveMultiple
	print 'False Positive   						:',falsePositive
	print 'False Positive with multiple alignment 				:',falsePositiveMultiple
	print 'Not Mapped       						:', notMapped
	print 'Total            						:', total	
	return

# Usage information. Copied from the reference version.
def usage():
	print "Usage:"
	print " -i <input sam file to be analyzed>"
	print "\t SAM file to look at"
	print " -t <1>|<2>"
	print "\t 1 for the second codon SAM file with second codon BWT"
	print "\t 2 for the regular SAM file with regular BWT"
	print "--clean | -c"
	print "\t Remove all SAM realted files"
	print "--help | -h"
	print "\t Help text which is what you see :)"

if __name__ == '__main__':
	main()

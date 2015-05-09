from bitarray import bitarray

FASTA_SUFFIX = '.fa'

REF_PREFIX_CLEAN = 'Cleaned_'
REF_PREFIX_SEC = 'SBWT_Reference_2nd_'
QUERY_SEC_PREFIX = 'SBWT_Query_2nd_'
QUERY_ORI_PREFIX = 'Random_Original_'
QUERY_MUT_PREFIX = 'Random_Mutated_'

RANDOM_SEQ_CODON_LENGTH = 20
RANDOM_SEQ_LENGTH = RANDOM_SEQ_CODON_LENGTH * 3
RANDOM_SEQ_NUMBER = 10000

ID_DELIMIETER = "|"
TAG_DELIMIETER = ';'

MUTATION_RATE_1 = float(0.05)
MUTATION_RATE_2 = float(0.05)
MUTATION_RATE_3 = float(0.05)
REVERSE_COMP_RATE = float(0.30)

CHECK_FIRST_CODON = False
FIRST_CODON_THRESHOLD = 1

ALN_OUTPUT_PREFIX = 'ALN_'
ALN_OUTPUT_SUFFIX = '.sai'
SAM_OUTPUT_PREFIX = 'SAM_'
SAM_OUTPUT_SUFFIX = '.sam'
COMMAND_BUILD_INDEX = './bwa index '
COMMAND_FIND_ALN = './bwa aln -n 0 -o 0 -l 10000000 -d 0 -i 0 -k 0 -f '
COMMAND_SAM_FORMAT = './bwa samse -n 1000000 -f '

BIT_DICT = {'A':bitarray('0001'), 'C':bitarray('0010'), 'G':bitarray('0100'), 'T':bitarray('1000'), 'a':bitarray('0001'), 'c':bitarray('0010'), 'g':bitarray('0100'), 't':bitarray('1000')}

# Second Codon Position Burrows Wheeler Transform Alignment
This is a program which implements the second codon position borrows wheeler transform alignment. 
It heavily relies on other package and existing BWA (Borrows Wheeler Aligner).
The algorithm basically uses the second codon position only to perform the Burrows Wheeler Transform. Also, it performs the first codon position similarity check when requested.

##Packages or other sources to install
* pysam 
  * A Sam file formatting package from python
  * https://pypi.python.org/pypi/pysam
* biopython 
  * Biopython package for bioinformatics in python
  * http://biopython.org/
* biarray 
  * Bitarray python package.
  * https://pypi.python.org/pypi/bitarray
* bwa
  * Burrows Wheeler Alignment program in C 
  * Python program needs the bwa binary to be placed in same directory to work properly. 
  * Included bwa was compiled in MAC OS-X, try to download the source and compile for corresponding machine if the included binary source does not work.
  * http://bio-bwa.sourceforge.net/

##Data Used
* I used the human DNA chromosome 1 data from UCSC for testing after cleaning it to remove N which is unknown section. (http://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/chr1.fa.gz)
* I used the refMrna from UCSC genome. (http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/refMrna.fa.gz)

##Program Instruction
### cleanSeq.py
This is a program to clean out the sequence. USCS data contains N for the unknown nucleotide which is handled as random character in the **bwa**. To remove all of the unknown variable, this program is created to clean out the sequence. It removes all other characters other than [acgt]. It creates an output file with prefix of (Cleaned_) 

* Output file form:
  * **Cleaned_(input_file_name).fa**
```
Usage:
 -i <input file name for reference sequence>
	 The input file name in fasta format which contains the reference sequence to be cleaned
	 This will remove mutable sequence and remove repeat
	 Assumes lower character as mutable sequence such as a,c,g,t
	 Assumes N as the repeat sequence
--clean | -c
	 Remove all cleaned reference files
--help | -h
	 Help text which is what you see :)
```

### randomQuery.py
This is a program to generate random sequences to use for the alignment testing given the reference sequence. It takes various parameters such as number of sequences, number of codons for each sequences, each codon position mutation rate. It creates two output files, one with the original randomly selected sequences with prefix of (Random_Original_), the other with mutated randomly selected sequences with prefix of (Random_Mutated_)
Each sequences ID contains its original position as well as its original sequence ID thus it can be retrieved easily.
[Sequence ID] | [Sequence Number] | [Original Starting Position] |[Offset]


* Output file form:
  * **Random_Original_(input_file_name).fa**
  * **Random_Mutated_(input_file_name).fa**

```
Usage:
 -i <input file name for reference sequence>
	 The input file name in fasta format which contains the reference sequence
 -1 <mutation rate for first codon position>
	 The values between 0 and 1. Default is set to 0.05>
 -2 <mutation rate for second codon position>
	 The values between 0 and 1. Default is set to 0.05>
 -3 <mutation rate for third codon position>
	 The values between 0 and 1. Default is set to 0.05>
 -n <number of random sequence to generate.>
	 The number of random sequence to generate. Default is set to 10000
 -l <number of codons for each random sequence.>
	 This is number of codons not nucleotides. Default is set to 20 and its nucleotides length will be 60
--clean | -c
	 Remove all random query generation related files
--help | -h
	 Help text which is what you see :)
```

### rbwt.py
This is a program to process regular **bwa** with specifically predefined parameters. To test with naive algorithm I defined following parameter for the **bwa** program.

'./bwa aln -n 0 -o 0 -l 10000000 -d 0 -i 0 -k 0 -f '

To change **bwa** parameter, use the constant parameter COMMAND_FIND_ALN defined in **const.py**. These predefined parameter will disabled seed, gap, indel, mismatches and use **bwa aln** which is commonly used for the short segments to enforce only burrows wheeler transform without any local alignment algorithm.

The program performs an indexing procedure, alignment procedure and converting the result to SAM file procedure. It creates multiple output files related to the **bwa** program but the important output is the result output in the SAM format which will be used for the analysis.

* Output file form:
  * **SAM_(query file).sam**

```
Usage:
 -r <reference file>
	 Reference file sequence in fasta format
 -q <query file>
	 Query file sequence in fasta format
--clean | -c
	 Remove all BWT realted files
--help | -h
	 Help text which is what you see :)
```

### sbwt.py
This is a program to process second codon position **bwa**. It uses the same **bwa** program with same predefiend parameter. 

'./bwa aln -n 0 -o 0 -l 10000000 -d 0 -i 0 -k 0 -f '

To change **bwa** parameter, use the constant parameter COMMAND_FIND_ALN defined in **const.py**. These predefined parameter will disabled seed, gap, indel, mismatches and use **bwa aln** which is commonly used for the short segments to enforce only burrows wheeler transform without any local alignment algorithm.

This is similar to the regular bwt but performs two preprocessing steps for the second codon position **bwa**

First, it performs to retrieve only second codon position from the reference sequence. 
Secondly, it performs a three shift on the reference seqeunce to take care of all of the possible frames from the query.  Since the **bwa** creates an index in forward and reverse complementary, just three different frames in query seqeunces covers all six frames. The program assumes that the sequence starts with first codon position.

The program performs an indexing procedure, alignment procedure and converting the result to SAM file procedure like regular **bwa** algorithm does after two preprocessing steps.

The program performs an indexing procedure, alignment procedure and converting the result to SAM file procedure. It creates multiple output files related to the **bwa** program but the important output is the result output in the SAM format which will be used for the analysis.

Furthermore, there is a parameter to perform additional 1st codon position similarity check with similarity threshold. If this option is used, its statistical result is also calculated from the program which will printed in the result without modifying any SAM result file.

* Output file form:
  * **SBWT_Query_2nd_(query file).sam**

```
Usage:
 -r <input file name for reference sequence>
	 The input file name in fasta format which contains the reference sequence
 -q <input file name for query sequence>
	 The input file name in fasta format which contains the query sequences
 -1
	 Set this flag to turn on the first codon position match. It is off by default
 -t <threshold for the 1st codon position match>
	 The values between 0 and 1. Default is set to 1 to match all>
--clean | -c
	 Remove all random query generation related files
--help | -h
	 Help text which is what you see :)
```

### checkResult.py
This is a program to check if the returned SAM result is correctly aligned or not. It uses the sequence ID which contains the original position to check if the returned alignment is correctly found or not. It considers the alignemnt as true positive if any of the found positions contain the original position. It considers the alignment as false positive if any of the found positions do not contain the original position. If the alignment was not found, it is counted as Not Found and did not included in True or False Positive.
The program also returns the number of True and False positive with multiple alignment.
The output is printed on the stdout and no output file is generated.

The sequence ID section is quite different depending on the algorithm that is used. 
Thus specify the type of the alignment (regular **bwa** or second codon position **bwa**) carefullly for this program to get correct result.

```
Usage:
 -i <input sam file to be analyzed>
	 SAM file to look at
 -t <1>|<2>
	 1 for the second codon SAM file with second codon BWT
	 2 for the regular SAM file with regular BWT
--clean | -c
	 Remove all SAM realted files
--help | -h
	 Help text which is what you see :)
```

### const.py
Contains constant parameter used throughout the program such as prefixes for the output file and default values for the parameters and programs.

## Example Runs
I will show the example run with the given the refMrnaSMALL.fa

1. Clean the reference sequence
  * >python cleanSeq.py -i refMrnaSMALL.fa
2. Use the cleaned reference sequence to generate random query sequences
  * >python randomQuery.py -i Cleaned_refMrnaSMALL.
3. Use the mutated random sequence to run the regular **bwa** algorithm
  * >python rbwt.py -r Cleaned_refMrnaSMALL.fa -q Random_Mutated_Cleaned_refMrnaSMALL.fa
4. Use the mutated random sequence to run the second codon position **bwa** algorithm
  * >python sbwt.py -r Cleaned_refMrnaSMALL.fa -q Random_Mutated_Cleaned_refMrnaSMALL.fa
5. Check the result.
  * >python checkResult.py -i SAM_Random_Mutated_Cleaned_refMrnaSMALL.sam -t 2
  * >python checkResult.py -i SAM_SBWT_Query_2nd_Random_Mutated_Cleaned_refMrnaSMALL.sam -t 1
6. Run the second codon position **bwa** with first codon similarity to compare as well.
  * >python sbwt.py -r Cleaned_refMrnaSMALL.fa -q Random_Mutated_Cleaned_refMrnaSMALL.fa -1
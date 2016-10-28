#cp-blast README File

This is a blast wrapper for identifying novel proteins by BLAST with a seed query and comparison to a database of existing intellecual property.

#Objective

In many cases, existing enzymes will have sub-optimal characteristics for the desired final pathway. Mining existing databases for alternatives is an inexpensive way to generate new options that can be tested in vivo. However, standard command line blast makes it difficult to determine which hits are closely related to existing, potentially patented homologs.

The direct_blast script is designed to search using BLAST and then screen the resulting hits against a secondary database to allow for quick and dirty identification of possible IP problems. Results are outputted as a TSV file.

The shift_blast script is designed to take advantage of bacterial operons, and to identify potential pathway components by BLAST against a genomic DNA database with a seed sequence and then identification of candidate ORFs nearby. Again, it can screen the resulting hits against a secondary database to allow for quick and dirty IP analysis. Results are outputted as a TSV file.

The pairwise_identity script is designed to do repeated pairwise alignments between all sequences in a file, outputting a TSV file which contains the percent identity for each unique combination. 

#Requirements

-	Python 3.4, 3.5

-	BioPython 1.67

	https://github.com/biopython/biopython
	
-	NCBI BLAST+ 2.2.31-4

	https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/

-	MUSCLE 3.8.31

	http://drive5.com/muscle/
	
	Edgar, R.C. (2004). MUSCLE: a multiple sequence alignment method with reduced time and space complexity.
    BMC bioinformatics 5, 113.
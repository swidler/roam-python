# roam-python
RoAM project translated into Python

Required modules:

numpy, math, scipy, copy, pysam, Bio, gzip, datetime, glob, re, pickle, sys, itertools, pybedtools

Input files and other variables should be specified in config.py.

Input files:

bam file with ancient genome and corresponding bai file

    examples:
		Ust Ishim (single stranded) can be found at https://www.ebi.ac.uk/ena/browser/view/PRJEB6622
		SF12 (double stranded) can be found at https://www.ebi.ac.uk/ena/browser/view/PRJEB21940
		If there is no bai file, it can be created using samtools: samtools index <file>.bam  <file>.bam.bai
		
		
genome assembly sequence file

    hg19.fa.gz can be downloaded from http://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/
		
modern sample

    this is currently loaded from a text file
    
To run RoAM, start by editing the variables in the config.py file. These include input and output filenames, sample name, 
	a list of chromosomes, their respective lengths, various flags, and the the parts of the script to run. When done, 
	run the run_roam.py script.
  
Warnings

    When running run_roam.py, there will be warnings about invalid values and divide by zeros. They can be safely ignored.

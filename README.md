# roam-python
RoAM project translated into Python

Required modules:

numpy, math, scipy, copy, pysam, Bio, gzip, datetime, glob, re, pickle,itertools, pybedtools

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
  


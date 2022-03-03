# roam-python
RoAM project translated into Python

Required modules:

numpy, math, scipy, copy, pysam, Bio, gzip, datetime, glob, re, pickle, sys, itertools, pybedtools, matplotlib

Input files and other variables should be specified in config.py.

Input files:

bam file with ancient genome and corresponding bai file

    examples:
		Ust Ishim (single stranded) can be found at https://www.ebi.ac.uk/ena/browser/view/PRJEB6622
		SF12 (double stranded) can be found at https://www.ebi.ac.uk/ena/browser/view/PRJEB21940
		If there is no bai file, it can be created using samtools: samtools index <file>.bam  <file>.bam.bai
		(In later versions of samtools, don't specify output file name.)

To index all bam files in a directory, use the indexing script:

		./create_bai.sh <directory>/*.bam
		
genome assembly sequence file

    hg19.fa.gz can be downloaded from http://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/
		
modern sample

    This can be loaded from a text file:
    Download bone5.zip from http://carmelab.huji.ac.il/data.html and unzip into script directory.
    Alternatively, use a Bismark result file in the BedGraph or Cov format and specify in the config file.

CpG coordinates
	
	cpg_coords.P is a pickled object containing the coordinates of all CpGs in the genome.
	Download the CpG coordinate file from http://carmelab.huji.ac.il/data.html and unzip into 
	objects subdirectory of script directory.
    
Running the scripts

    To run RoAM, start by editing the variables in the config.py file. These include input and output filenames, 
    sample name, a list of chromosomes, their respective lengths, various flags, and the the parts of the script to run. 
    The stages are "bam", "diagnose", "filter", "drate", and "meth".
    The first stage, bam, is the conversion of bam file(s) to Amsample object. It is a prerequisite to the other stages.
    It can be run by itself or with the other stages, but need not be run more than once.
    The diagnose stage computes basic statistics on each input chromosome, and recommends what thresholds to use when 
    excluding PCR duplicates and true mutations. It also generates various plots for sanity checks. These are stored as 
    .png files in the directory specified in the config.py file.
    The next step is filter, which removes information from CpG sites that did not pass various quality control tests.
    Next, drate, estimates the deamination rate.
    The last step, meth, computes methylation from c_to_t data, based on some function of the C->T ratio (no_t/no_ct).
    When done adjusting the config file, run the run_roam.py script.
    
    The DMR process has a similar flow. First edit the variables in the config_DMR.py file. These include directory and
    filenames, samples and group names, parameters for grouping the DMRs, and the parts of the script to run.
    The stages are "create_files", "DMR", "permute", "permutstat", and "plotmethylation".
    The first stage, create_files, takes text files from the RoAM process and converts them into pickled files that the
    rest of the process uses. It is a prerequisite to the other stages. It can be run by itself or with the other 
    stages, but need not be run more than once for a give group of samples.
    The DMR step compares the methylation of the samples and finds differentially methylated regions. It creates 
    both a text file and a pickled file of the DMR object, for use by the subsequent stages.
    The permute step scrambles up the DMR groups to validate the DMRs found and the permutstat step calculates the 
    statistical significance of the results.
    The plotmethylation step creates a plot of a specific DMR to compare its methylation between samples and groups.
    When done adjusting the config file, run the run_DMRs.py script.
  
Warnings

    When running run_roam.py, there will be warnings about invalid values and divide by zeros. They can be safely ignored.

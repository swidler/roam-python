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
    The first 4 variables, input (.bam) file, name, abbreviation, and library (single or double stranded) are required. 
    The rest have defaults loaded, but be sure to put all necessary files in the proper directory (default is the current directory).
    
    The script can also be run from the command line, using flags for the required parameters, as follows:
    run_roam.py -f "path to bam file" -n "sample name" -a "sample abbreviation" -l "library--single or double"
    The rest of the parameters can be specified as well, to override defaults (quoted strings, except where noted):
    -le chrom lengths, specified with no quotes or commas, eg -le 12345 23456 (must correspond with chromosome list, below)
	-s sample species
	-c list of chromosomes to use, a list specified with no quotes or commas, eg -c chr1 chr2 chr3
	-t set True to trim ends during processing
	-cf set True for dir with one file per chromosome
	-m mapping quality for read 
	-q mapping quality for position 
	-st stages of process to be run, a list specified with no quotes or commas, eg -st bam diagnose
	-o directory for saved (pickled) object files 
	-fd directory for multiple input files 
	-od output directory 
	-ld directory for log files 
	-pd directory for images 
	-g path for the assembly file, eg hg19.fa.gz 
	-tf path for the object saved in text format 
	-mo text file for modern genome 
	-b .cov or .bedGraph file for modern genome 
	-gc CpG file 
	-mn modern sample name 
	-ma modern sample abbreviation 
	-ms modern sample species 
	-mr modern sample reference genome 
	-mm modern sample sequencing method
	-bed flag for bed file (use this flag when bed file output is not desired) 
    
    The stages of the process (which can be specified in the config file) are "bam", "diagnose", "filter", "drate", and "meth".
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

	When running from the command line, parameters can be specified as follows:
	-s sample names, a list specified with no quotes or commas
	-g group names--should correspond with samples, a list specified with no quotes or commas
	-o directory for pickled objects
	-d directory for data files from RoAM process
	-du directory for output txt files and pics
	-gc CpG file
	-ge sorted text file with genes
	-cg CGI file
	-di pickled DMR file for use in fdr, permutations, and plots
	-dpi pickled DMR permutation file for use in permutstat
	-t template to match any extra text in sample filename
	-st stages of process to be run, a list specified with no quotes or commas
	-de minimum methylation difference between the two groups
	-mc DMRs whose number of CpGs is less than min_CpGs are filtered out
	-mq DMRs with Qt < min_qt are filtered out
	-mf minimum number of ancient samples for which we require data
	-w window size for smoothing
	-l low coverage factor
	-p number of permutations to run
	-dmi index of DMR
	-dmc chromosome of DMR
	-b .cov or .bedGraph file for modern genome
	-mo text file for modern genome
	-r reference genome for use in histogram matching
	-re flag for logging info--use when logging not desired

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

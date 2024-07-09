RoAM is a flexible tool with two algorithms allowing (1) reconstruction of ancient methylation and
(2) detection of differentially methylated regions (DMRs) distinguishing groups. Detailed
information can be found in THE PAPER

The following modules are required for RoAM and should be downloaded before using it:

numpy, math, scipy, copy, pysam, Bio, gzip, datetime, glob, re, pickle, sys, itertools, pybedtools, matplotlib, screeninfo

RoAM can be operated either from the command line or by specifying the inputs in config.ini as detailed in this file.

# Reconstructing ancient methylation algorithm

Input files:

bam file with ancient genome and corresponding bai file

    examples:
		Ust Ishim (single stranded) can be found at 
		https://www.ebi.ac.uk/ena/browser/view/PRJEB6622
		SF12 (double stranded) can be found at 
		https://www.ebi.ac.uk/ena/browser/view/PRJEB21940
		If there is no bai file, it can be created using samtools: samtools index <file>.bam  
		<file>.bam.bai (In later versions of samtools, don't specify output file name.)

To index all bam files in a directory, use the indexing script:

		./create_bai.sh <directory>/*.bam
		
genome assembly sequence file (used only when creating a new CpG coordinates file)

    hg19.fa.gz can be downloaded from http://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/
		
CpG coordinates
	
	If you are using the Hg19 genome, a pickled object containing the coordinates of all CpGs 
	in the 	format fit for the code already exists. Download it from 
	http://carmelab.huji.ac.il/data.html and unzip into the objects subdirectory of script 
	directory. To create a CpG file that corresponds to another reference genome, make sure 
	to download the assembly in fa.gz format. When running RoAM, use the --create_cpg flag 
	and include the path of the assembly file in the config.ini or in the input parameters. 
	The cpg file will be created and stored in the object directory, using the format
	<species_with_underscores>_cpg_coords.P. In order to create the file without running the 
	rest of RoAM, add the --no_roam flag to the input parameters.
    
Modern reference sample

	This is a present-day methylation map generated from the same tissue as the ancient 
	remains. When reconstructing a bone sample with Hg19 reference genome, the reference can 
	be	loaded from a text file:
	Download bone5.zip from http://carmelab.huji.ac.il/data.html and unzip into script 
	directory. Alternatively, use a Bismark result file in the .cov format and specify in the 
	config file or in the input parameters. The modern reference is not strictly required, 
	but it is highly recommended, since the more accurate methods of deamination estimation 
	and methylation reconstruction depend on it.
	
Running the scripts

    To run the process, start by editing the variables in the config.ini file (or run using 
    command line parameters, as described below). These include input and output filenames, 
    sample name, a list of chromosomes, their respective lengths, various flags, and the the 
    parts of the script to run. 
    
    The first 3 variables, input (.bam) file, name, and library (single or double stranded) 
    are required. If instead of one bamfile, the sample has a directory with all the files, 
    the filename param can be skipped. In this case, change the file_per_chrom flag to true 
    if the directory contains (exactly) one file for each chromosome. The filenames should 
    be in the format <your label>_chr<chrom>.bam. 
    Please note: if there is more than one file for a given chromosome, only one will be read. 
    Alternatively, the directory can have multiple bamfiles that do not correspond to 
    specific chromosomes, in which case this flag should not be changed from its default value 
    (false).
    
    The rest of the parameters have defaults loaded, but be sure to put all necessary files in 
    the proper directory (default is the current directory). 
    
    When done adjusting the config file, run the run_roam.py script.
    
    The script can also be run from the command line, using flags for the required parameters, 
    as follows:
    
    	run_roam.py -f "path to bam file" -n "sample name" -l "library--single or double"
    
    IMPORTANT: If the sample is aligned to any genome but hg19, it is important to modify the chromosome 
    lengths, either in the config or by using the -le parameter.
    The rest of the parameters can be specified as well, to override defaults (strings, except 
    where noted):
    -co path of config file (different config files can be used for different runs)
    -le chromosome lengths, specified with no quotes or commas, eg -le 12345 23456 (must 
    correspond with chromosome list, below) (default: chromosome lengths of hg19).
	-s sample species (default: human)
	-c list of chromosomes to use, a list specified with no quotes or commas, 
	eg -c chr1 chr2 chr3 	(default: chromosomes 1 to 22, X)
	-t set True to trim ends during processing (default: false)
	-cf set True for dir with exactly one file per chromosome (default: false)
	-m minimum mapping quality for read (default: 20)
	-q minimum mapping quality for position (default: 20)
	-st stages of process to be run, a list specified with no quotes or commas, 
	eg -st bam diagnose Further details later in this document (default: all five stages)
	-o directory for saved (pickled) object files (include final / for all directories) 
	(default: current directory)
	-fd directory for multiple input files (default: current directory)
	-od output directory (default: current directory)
	-ld directory for log files (default: current directory)
	-pd directory for images (default: current directory)
	-g path for the assembly file, eg hg19.fa.gz (default: hg19.fa.gz in current directory)
	-tf path for the object saved in text format 
	-mo text file for modern genome 
	-b .cov file for modern genome 
	-gc CpG file 
	-mn modern reference sample name (default:empty)
	-ms modern reference sample species (default:empty)
	-mr modern sample reference genome (default:empty)
	-mm modern reference sample sequencing method (default:empty)
	-bed flag for bed file (use this flag when bed file output is not desired) (default: false)
	-cpg flag for creating CpG file (default: false)
	-cr reference genome assembly for CpG file
	-no flag for not running the rest of RoAM after creating CpG file (default: false)
	-u determines which threshold to use--from diagnose (if true), or user-entered (default: true)
	-ct user-entered threshold used to identify sites with a true C->T mutation (used only if 
	use_max_TsPerCoverage is false) (default: 0.25)
	-mg merge two consecutive coordinates of every CpG position (default: true)
	-ga threshold used to identify sites with a true C->T mutation. Only for library='single' (default: 0.25)
	-me method used to remove true C->T mutations (default: c_to_t for library='double', otherwise both) 
	-dm method of deamination rate calculation (can be reference [default and 
	highly recommended] or global)
	-mc minimum coverage of sites for deamination rate calculation (default: 1)
	-mb minimum beta value to consider for the estimation in the reference method of deamination rate calculation (default: 1)
	-rm method of reconstruction (can be histogram, linear, or logistic; default is histogram) 
	-gm global methylation value for use with global method of deamination rate calculation
	-lcf low coverage factor for methylation reconstruction (default: 0.05)
	-sl slope for linear/logistic methods of methylation reconstruction (default: 1/d_rate)
	-in intercept for linear method of methylation reconstruction (default: 0)
	-w window size for reconstruction of methylation--'auto' or list of 1 val or 1 for each 
	chromosome (default: auto)
	-wm window size calculation method--'prob' or 'relerror' (default: prob)
	-min minimum methylation level to detect (default: 0.2)
	-p param used in winsize calculation for win_method = prob (default: 0.01)
	-k reciprocal of param used in winsize calculation for win_method = relerror 
	(default: 1/2.5)
	-max maximum window size (default: 31)
	
Outputs
	
	If the process runs to the final stage, two main outputs are given: a bed file 
	(<sample>_meth.bed) with the methylation values per position, and a pickled file that 
	can be 	used for the DMR detection process. The process has five stages: "bam", "diagnose", 
	"filter", "drate", and "meth". These stages can be 	specified in the config file and from 
	the command line (input -st), allowing the user to choose whether to run the full algorithm, 
	or just several stages. When the script ends (after the last requested stage), it outputs 
	a text file for the last stage completed. This file, in the format <sample>_<stage>.txt, 
	can be found in the specified output directory and can be used as input for later stages. 
    
The stages
    
    The first stage, bam, is the conversion of bam file(s) to Amsample object. It is a 
    prerequisite to the other stages. It can be run by itself or with the other stages, but 
    need not be run more than once.
    
    The diagnose stage computes basic statistics on each input chromosome, and recommends 
    the thresholds to use to remove PCR duplicates and true mutations. The information is 
    recorded in a file in the specified log directory, in the format <sample>_diagnostics.txt. 
    It also generates various plots for sanity checks. These are stored as .png files in the 
    directory specified in the config.ini file.
    
    The next step is filter, which removes information from CpG sites that did not pass 
    various quality control tests, and merges the counts on the plus and minus strands in 
    each CpG. The details of what was removed are stored in a file in the specified log 
    directory, in the format <sample>_filter.txt.
    
    Next, drate, estimates the deamination rate.
    
    The last step, meth, reconstructs methylation from c_to_t data, based on some function of 
    the C->T ratio (no_t/no_ct).
    
    When done adjusting the config file, run the run_roam.py script.
    
Warnings

    When running run_roam.py, there will be warnings about invalid values and divide by zeros. 
    They can be safely ignored.
    
# Detect DMRs

	This process receives samples divided into two groups and searches for DMRs between them.

Preparing modern samples

	If one of the groups is made of modern samples, then they need to be formatted properly 
	for the code before running the DMR process. The files are created from .cov files 
	(Bismark output). Run the create_files_for_DMRs.py script, using the following input 
	parameters (quoted strings, 	unless otherwise noted):
	
	-ms sample names, a list specified with no quotes or commas (mandatory)
	-gc CpG file (mandatory)
	-o directory for saved output (pickled) object files (include final /) (optional)
	-d directory for data files from RoAM process (include final /) (optional)
	-mt template to match any extra text in modern sample filename (optional)
	-msp modern sample species (optional)
	-mr modern sample reference genome (optional)
	-mm modern sample sequencing method (optional)
		
	If any of the last 3 parameters differ for different modern samples, they will need to be 
	run separately.

Running the DMR detection process

    The DMR process has a similar flow to the reconstruction process. First edit the 
    variables in the config_DMR.ini file. These include directory and filenames, samples 
    and group names, parameters for grouping the DMRs, and the parts of the script to run.
    
    	run_DMRs.py -s sample1 sample2 -g group1 group2 -gc “CpG file path”
    
Inputs

	Three inputs are required for the process to work:
	a. Samples (-s): a list of all the samples that the process should compare. They should 
	be	specified with no quotes or commas. The algorithm will use the pickled files with the 
	same name (or with any extra text specified in the templ parameter) from the object 
	directory for the corresponding samples.
	b. Groups (-g): a list with identical size to “samples” with allocation of the samples 
	to the	groups that the process detect DMRs between. Groups should include ancient samples 
	or	modern samples, but not both. 
	For example: -s ust_ishim SF12 Altai_Neanderthal Denisovan -g modern_human modern_human archaic_human archaic_human.
	c. The CpG file (-gc), identical to the file in methylation reconstruction.
	
	In addition to these three, modern reference methylation map: this will be used to 
	generate simulations. Just like the methylation reconstruction process this can be 
	provided as a .cov file from Bismark (-b) or as a text file (-mo) in the same format as 
	“bone5” described earlier (bone5 can be found in http://carmelab.huji.ac.il/data.html). 
	If no -b or -mo input is provided, the	process will use the reference defined in the 
	config file of the methylation reconstruction	process, and its path can be specified 
	with the parameter -rco. If no such config file exist, the process will not be able 
	to run.
    
    Other optional parameters can be specified as follows:
	-co path of config file (different config files can be used for different runs)
	-rco path of RoAM config file
	-b .cov file for modern genome
	-mo text file for modern genome
	-ms modern sample names (list like samples)
	-mn modern reference sample name 
	-msp modern reference sample species 
	-mr modern sample reference genome 
	-mm modern reference sample sequencing method
	-c list of chromosome names
	-o directory for saved output (pickled) object files (include final / for all directories)
	-du directory for output txt files and pics
	-gc CpG file
	-ge sorted text file with genes and their coordinates for annotation
	-cg BED file with CGI coords for annotation
	-c1 first custom bed file for annotation
	-c2 second custom bed file for annotation
	-pd promoter definition around TSS: a list of 2 values [before, after], where before 
	is the 	number of nucleotides into the intergenic region, and after is the number of 
	nucleotides (default: [5000, 1000])
	-di pickled DMR file for use in fdr, permutations, and plots
	-t template to match any extra text in sample filename
	-st stages of process to be run, a list specified with no quotes or commas
	-de minimum methylation difference between the two groups to be considered a DMR (default 0.5)
	-mc (min_CpGs) DMRs whose number of CpGs is less than this value are filtered out 
	(default: 10)
	-mq (min_qt) DMRs with Qmax < this value are filtered out (default: 0)
	-mb (min_bases) DMRs shorter than min_bases are filtered out (default: 100)
	-mad (max_adj_dist) max distance between adjacent CpGs within the same DMR. If the 
	distance between consecutive CpG positions is larger than max_adj_dist, the algorithm sets 
	Qt to 0 (default: 1000)
	-mf (min_fin) list with 2 elements (1 for each group), each representing the minimum 
	number of samples for which we require data at a given position. If in a position there 
	are not enough samples with data, a NaN is	substituted in this position. It can also be a 
	fraction between 0 and 1 (0 <= min_fin < 1), in which case it is understood as the minimum 
	fraction of the total number of samples in the group (default: 1 for each group).
	-w window size for smoothing (default: calculated automatically)
	-l low coverage factor, positions with coverage below this value will be filtered out for 
	the corresponding sample (default: calculated automatically)
	-sp number of permutations to run for fdr
	-th FDR threshold (default: 0.05)
	-r flag for using histogram matching in the pooled methylation function (default: true)
	-re flag for logging info--use when logging not desired
	-an flag for running annotation--use when annotation not desired
    
    The stages are DMR and fdr. (Extra functions that can be used here are permute, 
    permutstat, and plotmethylation. For more details see “extra functions”).
    
    The DMR step compares the methylation of the samples and finds differentially methylated 
    regions. It creates 2 text files in the specified dump directory (DMR_gen_<timestamp>.txt 
    lists input parameters and results, while DMRs_<timestamp> has the data in tab separated 
    format), and a pickled file of the DMR object in the specified object directory 
    (DMR_obj_<timestamp>), for use by the subsequent stages. It also outputs a log file in 
    the current dir (DMR_log.txt) with info about the parameters used and DMRs found.
    
    The fdr (False Detection Rate) step filters the DMRs to the desired FDR level. It also 
    notes the thresholds used at the end of the DMR_log file created in the DMR step. It 
    outputs 2 text files in the dump directory (filtered_DMR_gen_<timestamp>.txt and
    filtered_DMRs_<timestamp>, 
    as above), and a log file in the current dir with info about the thresholds and the ratio 
    obtained for each combination.
    
    When done adjusting the config file, run the run_DMRs.py script.    
    
Outputs
    
	The algorithm provides a table with a list of all the DMRs detected and details about 
	them. It also	outputs a pickle file that can be used later for the extra functions.

Extra functions

	We provide extra functions for users who wish deeper analysis of the data and/or to plot 
	mean methylation of the DMRs. To operate them, change the parameter “stages” to the 
	corresponding	term.
	
	The permute function scrambles up the DMR groups and detects “permuted DMRs” to validate
	the DMRs found. The input for the function is the pickled output of the DMR detection 
	process. It can be set via the config file or input parameters. The number for 
	permutations can also be set using the num_permutations parameter in the config file or 
	the input param:
	
	-p number of permutations to run
	
	For each permutation, it creates 2 text files in the dump directory as above
	(DMR_gen_<timestamp>.<permutation #>.txt and DMRs_<timestamp>.<permutation #>.txt), and a 
	log file in the current dir (p<permutation #>_groupDMRs.txt) with info on parameters used 
	and DMRs found. This step also outputs a pickled file of the DMR object in the specified 
	object directory (DMP_obj_<timestamp>) for use in the permutstat function.
    
    Permute is followed by the permutstat function, which calculates the statistical 
    significance of the results of the permute step. 
    
    -dpi pickled DMR permutation file for use in permutstat
	
    It outputs a text file in the dump 
    directory (pstat_<timestamp>.txt) with the resulting statistics and p-values.
    
    The plotmethylation function creates a scatter plot of the mean methylation in each 
    sample for a specific DMR chosen for this function (meth_plot<chromosome>_<DMR index>.png, 
    in the dump directory) to compare its methylation between samples and groups. The inputs 
    for the function are the pickled file of the detected DMRs, the chromosome of the DMR that 
    the user wishes to plot, and its (0 based) index number within the chromosome, which can be 
    found in the DMR table output. In this file DMRs are ordered by chromosome. If, for example, 
    the selected DMR is the third DMR found in chromosome 5, then the value for this parameter 
    should be “2”.
	
	-dmi index of DMR
	-dmc chromosome of DMR
	

  

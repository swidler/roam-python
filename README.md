# roam-python
RoAM project translated into Python

Required modules:

numpy, math, scipy, copy, pysam, Bio, gzip, datetime, glob, re, pickle, sys, itertools, pybedtools, matplotlib

Input files and other variables should be specified in config.ini.

Input files:

bam file with ancient genome and corresponding bai file

    examples:
		Ust Ishim (single stranded) can be found at https://www.ebi.ac.uk/ena/browser/view/PRJEB6622
		SF12 (double stranded) can be found at https://www.ebi.ac.uk/ena/browser/view/PRJEB21940
		If there is no bai file, it can be created using samtools: samtools index <file>.bam  <file>.bam.bai
		(In later versions of samtools, don't specify output file name.)

To index all bam files in a directory, use the indexing script:

		./create_bai.sh <directory>/*.bam
		
genome assembly sequence file (used only when creating a new CpG coordinates file)

    hg19.fa.gz can be downloaded from http://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/
		
modern sample

    This can be loaded from a text file:
    Download bone5.zip from http://carmelab.huji.ac.il/data.html and unzip into script directory.
    Alternatively, use a Bismark result file in the Cov format and specify in the config file.
    The modern reference is not strictly required, but it is highly recommended, since the more accurate
    methods of deamination estimation and methylation reconstruction depend on it.

CpG coordinates
	
	cpg_coords.P is a pickled object containing the coordinates of all CpGs in the human genome (Hg19).
	Download the CpG coordinate file from http://carmelab.huji.ac.il/data.html and unzip into 
	objects subdirectory of script directory. To using a custom CpG file, make sure to download the 
	assembly in fa.gz format. When running RoAM, use the --create_cpg flag and include the path of the 
	assembly file in the config.ini or in the input parameters. The cpg file will be created and stored in 
	the object directory, using the format <species_with_underscores>_cpg_coords.P. In order to create
	the file without running the rest of RoAM, add the --no_roam flag to the input parameters.
    
Running the scripts

    To run RoAM, start by editing the variables in the config.ini file. These include input and output filenames, 
    sample name, a list of chromosomes, their respective lengths, various flags, and the the parts of the script to run. 
    The first 4 variables, input (.bam) file, name, abbreviation, and library (single or double stranded) are required. 
    The rest have defaults loaded, but be sure to put all necessary files in the proper directory 
    (default is the current directory).
    
    The script can also be run from the command line, using flags for the required parameters, as follows:
    run_roam.py -f "path to bam file" -n "sample name" -a "sample abbreviation" -l "library--single or double"
    The rest of the parameters can be specified as well, to override defaults (strings, except where noted):
    -co path of config file (different config files can be used for different runs)
    -le chrom lengths, specified with no quotes or commas, eg -le 12345 23456 (must correspond with chromosome list, below)
	-s sample species
	-c list of chromosomes to use, a list specified with no quotes or commas, eg -c chr1 chr2 chr3
	-t set True to trim ends during processing
	-cf set True for dir with one file per chromosome
	-m mapping quality for read 
	-q mapping quality for position
	-rm method of reconstruction (can be histogram, linear, or logistic; default is histogram) 
	-st stages of process to be run, a list specified with no quotes or commas, eg -st bam diagnose
	-o directory for saved (pickled) object files (include final / for all directories)
	-fd directory for multiple input files 
	-od output directory 
	-ld directory for log files 
	-pd directory for images 
	-g path for the assembly file, eg hg19.fa.gz 
	-tf path for the object saved in text format 
	-mo text file for modern genome 
	-b .cov file for modern genome 
	-gc CpG file 
	-mn modern reference sample name 
	-ma modern reference sample abbreviation 
	-ms modern reference sample species 
	-mr modern sample reference genome 
	-mm modern reference sample sequencing method
	-bed flag for bed file (use this flag when bed file output is not desired) 
	-cpg flag for creating CpG file
	-cr reference genome assembly for CpG file
	-no flag for not running the rest of RoAM after creating CpG file
	-dm method of deamination rate calculation (can be reference [highly recommended] or global)
	-mc minimum coverage of sites for deamination rate calculation
	-mb minimum beta value for reference method of deamination rate calculation
	-gm global methylation value for use with global method of deamination rate calculation
	-lcf low coverage factor for methylation reconstruction
	-sl slope for linear/logistic methods of methylation reconstruction
	-in intercept for linear method of methylation reconstruction
	-w window size for reconstruction of methylation--'auto' or list of 1 val or 1 for each chrom
	-wm window size calculation method--'prob' or 'relerror'
	-min minimum methylation level to detect
	-p param used in winsize calculation for win_method = prob
	-k reciprocal of param used in winsize calculation for win_method = relerror
	-max maximum window size
    
    The stages of the process (which can be specified in the config file) are "bam", "diagnose", "filter", "drate", and "meth". 
    When the script ends (after the last requested stage), it outputs a text file for the last stage completed. 
    This file, in the format <sample>_<stage>.txt, can be found in the specified output directory and can 
    be used as input for later stages. Additionally, after the last stage (reconstruct methylation), a bed
    file (<sample>_meth.bed) is also produced, unless the user has specifically prevented this, using the
    parameter in the config file or from the command line.
    
    The first stage, bam, is the conversion of bam file(s) to Amsample object. It is a prerequisite to the other stages.
    It can be run by itself or with the other stages, but need not be run more than once.
    
    The diagnose stage computes basic statistics on each input chromosome, and recommends what thresholds to use when 
    excluding PCR duplicates and true mutations. The information is recorded in a file in the specified log directory, 
    in the format <sample>_diagnostics.txt. It also generates various plots for sanity checks. These are stored as 
    .png files in the directory specified in the config.ini file.
    
    The next step is filter, which removes information from CpG sites that did not pass various quality control tests.
    The details of what was removed are stored in a file in the specified log directory, in the format <sample>_filter.txt.
    
    Next, drate, estimates the deamination rate.
    
    The last step, meth, computes methylation from c_to_t data, based on some function of the C->T ratio (no_t/no_ct).
    
    When done adjusting the config file, run the run_roam.py script.

	Before running the DMR process, the output files from the RoAM process must be converted into pickled files. 
	If using modern samples as well, these are created at this stage from .cov files (Bismark output).
	Run the create_files_for_DMRs.py script, using the following input parameters (quoted strings, unless otherwise noted):
	
	-s (ancient) sample names, a list specified with no quotes or commas
	-ms modern sample names, a list specified with no quotes or commas
	-o directory for saved (pickled) object files (include final /)
	-d directory for data files from RoAM process (include final /)
	-gc CpG file (for modern samples)
	-ttemplate to match any extra text in sample filename
	-mt template to match any extra text in modern sample filename
	-ma modern sample abbreviation, a list specified with no quotes or commas (same length as modern samples)
	-msp modern sample species
	-mr modern sample reference genome
	-mm modern sample sequencing method
		
	If any of the last 3 parameters differ for different modern samples, they will need to be run separately.

    The DMR process has a similar flow. First edit the variables in the config_DMR.ini file. These include
    directory and filenames, samples and group names, parameters for grouping the DMRs, and the parts of 
    the script to run.
    
    When running from the command line, parameters can be specified as follows:
	-co path of config file (different config files can be used for different runs)
	-rco path of RoAM config file
	-s sample names, a list specified with no quotes or commas
	-ms modern reference sample names (list like samples)
	-mn modern reference sample name 
	-ma modern reference sample abbreviation 
	-msp modern reference sample species 
	-mr modern sample reference genome 
	-mm modern reference sample sequencing method
	-g group names--should correspond with samples, a list specified with no quotes or commas
	-o directory for saved (pickled) object files (include final / for all directories)
	-du directory for output txt files and pics
	-gc CpG file
	-ge sorted text file with genes
	-cg CGI file
	-c1 first custom bed file
	-c2 second custom bed file
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
	-sp number of permutations to run for fdr
	-p number of permutations to run
	-dmi index of DMR
	-dmc chromosome of DMR
	-b .cov file for modern genome
	-mo text file for modern genome
	-r reference genome for use in histogram matching
	-re flag for logging info--use when logging not desired
	-an flag for running annotation--use when annotation not desired
    
    The stages are "DMR", "fdr", "permute", "permutstat", and "plotmethylation".
    
    The DMR step compares the methylation of the samples and finds differentially methylated regions. It creates 
    2 text files in the specified dump directory (DMR_gen_<timestamp>.txt lists input parameters and results,
    while DMRs_<timestamp> has the data in tab separated format), and a pickled file of the DMR object in the specified 
    object directory (DMR_obj_<timestamp>), for use by the subsequent stages. It also outputs a log file in the current dir
    (DMR_log.txt) with info about the parameters used and DMRs found.
    
    The fdr (False Detection Rate) step filters the DMRs to the desired FDR level. It also notes the thresholds used at the end 
    of the DMR_log file created in the DMR step. It outputs 2 text files in the dump directory (filtered_DMR_gen_<timestamp>.txt
    and filtered_DMRs_<timestamp>, as above), and a log file in the current dir with info about the thresholds and the ratio 
    obtained for each combination.
    
    The permute step scrambles up the DMR groups to validate the DMRs found. For each permutation, it creates 2 text files 
    in the dump directory as above (DMR_gen_<timestamp>.<permutation #>.txt and DMRs_<timestamp>.<permutation #>.txt), and a 
    log file in the current dir (p<permutation #>_groupDMRs.txt) with info on parameters used and DMRs found. This step also 
    outputs a pickled file of the DMR object in the specified object directory (DMP_obj_<timestamp>) for use
    in the permutstat stage.
    
    The permutstat step calculates the statistical significance of the results of the permute step. It outputs a text file 
    in the dump directory (pstat_<timestamp>.txt) with the resulting statistics and p-values.
    
    The plotmethylation step creates a scatter plot of a specific DMR (meth_plot_<chromosome>_<DMR index>.png, in the 
    dump directory) to compare its methylation between samples and groups.
    
    When done adjusting the config file, run the run_DMRs.py script.
  
Warnings

    When running run_roam.py, there will be warnings about invalid values and divide by zeros. They can be safely ignored.

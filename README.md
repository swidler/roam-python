RoAM is a flexible tool including two pipelines allowing (1) the reconstruction of ancient DNA methylation, and
(2) detection of differentially methylated regions (DMRs) distinguishing two groups of samples. Detailed
information can be found in BIORXIV_LINK.

The following modules are required for RoAM and should be downloaded before using it:
numpy, math, scipy, copy, pysam, Bio, gzip, datetime, glob, re, pickle, sys, itertools, pybedtools, matplotlib, screeninfo

RoAM can be operated either from the command line or by specifying input parameters in config.ini.

# Pipeline 1: Reconstructing ancient DNA methylation

Input files:

1. BAM file with the reads of an ancient genome, and the corresponding BAI file.
   	Examples:
   		Ust Ishim (single stranded) can be found at https://www.ebi.ac.uk/ena/browser/view/PRJEB6622
		SF12 (double stranded) can be found at https://www.ebi.ac.uk/ena/browser/view/PRJEB21940
   	If there is no BAI file, it can be created using samtools:
   		samtools index <file>.bam <file>.bam.bai
   	(In later versions of samtools, don't specify output file name.)

   	To index all bam files in a directory, use the indexing script:
		./create_bai.sh <directory>/*.bam
		
3. Genome assembly sequence FASTA file (used only when creating a new CpG coordinates file).
   	hg19.fa.gz can be downloaded from http://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/
		
4. CpG coordinates
	If you are using the Hg19 human genome version, a pickled object containing all CpG
	coordinates in the right format already exists. Download it from
	http://carmelab.huji.ac.il/data.html and unzip into the objects subdirectory of the script
	directory.
	To create a CpG file that corresponds to another reference genome version, make sure 
	to download the assembly file in fa.gz format. When running RoAM, use the --create_cpg flag 
	and include the path of the assembly file in the config.ini or in the command line. 
	The cpg file will be created and stored in the objects directory, using the format
	<species_with_underscores>_cpg_coords.P. In order to create the file without running the 
	rest of RoAM, add the --no_roam flag to the input parameters.
    
5. Reference DNA methylation map
	Typically, this would be a present-day methylation map generated from the same tissue
	as the ancient sample.
	When reconstructing methylation from a bone sample using the Hg19 human genome version,
	a reference DNA methylation map (called bone1.zip) can be downloaded from
	http://carmelab.huji.ac.il/data.html and unzip into the script directory.
	Alternatively, use a Bismark output file in COV format which you should specify in the 
	config file or in the command line. The modern reference is not strictly required, 
	but it is highly recommended, as it allows for more accurate estimation of deamination
	rate.
	
Running the scripts:

	To run the process, start by editing the variables in the config.ini file (or edit 
	the command line parameters). These include input and output filenames, sample name,
 	a list of chromosomes, their respective lengths, various flags, and the the parts
  	of the script to run. 
    
	The first 3 variables, filename (input .bam file), name (identifier of the sample),
 	and library (single or double stranded) are required. If instead of one BAM file, the
  	directory contains (exactly) one file per chromosome, the filename param can be
   	skipped. In this case, change the file_per_chrom flag to TRUE. The filenames should
    	be in the format <your label>_chr<chrom>.bam. 
	Please note: if there is more than one file for a given chromosome, only one will be read. 
	Alternatively, the directory can have multiple bamfiles that do not correspond to
 	specific chromosomes, in which case this flag should not be changed from its default value
  	(FALSE). ** It is unclear what would happen in this last scenario. Which of the BAM files would be read? **
    
	The rest of the parameters have defaults loaded, but be sure to put all necessary files in 
	the proper directory (default is the current directory).
    
	When done adjusting the config file, run the run_roam.py script. The script can also be run
 	from the command line, using flags for the required parameters,	as follows:
    
    	run_roam.py -f "path to bam file" -n "sample name" -l "library--single or double"
    
	IMPORTANT: If the sample is aligned to any genome version but Hg19, it is important to
 	modify the chromosome lengths, either in the config file or by using the -le parameter.
	The rest of the parameters can be specified as well, to override defaults (strings, except 
	where noted):
 
 	-co path of config file (different config files can be used for different runs).
    	-le chromosome lengths, specified with no quotes or commas, eg -le 12345 23456 (must 
		correspond to chromosome list below) (default: chromosome lengths of hg19).
	-s taxonomical classification of the ancient sample (default: human).
	-c list of chromosomes to use, a list specified with no quotes or commas, 
		e.g., -c chr1 chr2 chr3 (default: chromosomes 1 to 22, X).
	-t set to TRUE to trim read ends during processing (default: FALSE).
	-cf set TRUE for directory with exactly one file per chromosome (default: FALSE).
	-m minimum mapping quality for a read (default: 20).
	-q minimum mapping quality for a position (default: 20).
	-st stages of process to be run, a list specified with no quotes or commas, 
		e.g., -st bam diagnose.
  		Further details later in this document (default: all five stages)
	-o directory for saved (pickled) object files (include final / for all directories) **What does this mean (the part in parenthesis)?**
		(default: current directory)
	-fd directory for multiple input BAM files (default: current directory).
	-od output directory (default: current directory).
	-ld directory for log files (default: current directory).
	-pd directory for images (default: current directory).
	-g path for the reference genome assembly file, e.g., hg19.fa.gz
 		(default: hg19.fa.gz in current directory).
	-tf path for the object saved in text format. ** default **
	-mo file for Reference DNA methylation map (TXT format).
	-b file for Reference DNA methylation map (COV format). 
	-gc CpG coordinates file.
	-mn sample name of reference DNA methylation map (default:empty).
	-ms taxonomical classification of reference DNA methylation map (default:empty).
	-mr reference genome of reference DNA methylation map (default:empty).
	-mm sequencing method of reference DNA methylation map (default:empty).
	-bed flag for bed file (use this flag when bed file output is not desired)
 		(default: FALSE). **unclear**
	-cpg flag for creating CpG coordinates file (default: FALSE). **above it says that the flag is called --create_cpg**
	-cr genome assembly verstion for CpG coordinates file. ** when would this be
 		different than the regular genome assembly? **
	-no flag for not running the rest of RoAM after creating CpG file (default: FALSE). ** above this is indicated as --no_roam **
	-u determines whether filtering thresholds would be taken from diagnose (TRUE),
 		or would be user-defined (FALSE) (default: TRUE).
	-ct user-defined threshold to identify sites with a true C->T mutation (used only
 		if use_diagnose_filter is FALSE) (default: 0.25). ** what is use_diagnose_filter? **
	-mg merge information from the two strands for each CpG position (default: TRUE).
	-ga threshold used to identify sites with a true C->T mutation. Only relevant when
 		library='single' (default: 0.25).
	-me method used to remove true C->T mutations (default: c_to_t for library='double',
 		otherwise both).
	-dm method of estimating deamination rate (can be reference [default and 
		highly recommended] or global).
	-mc minimum coverage of sites used for estimation of deamination rate
 		(default: 1).
	-mb minimum beta value of sites used to estimate deamination rate in the
 		reference method (default: 1).
	-rm method of DNA methylation reconstruction (can be histogram, linear, or
 		logistic; default is histogram).
	-gm global methylation value. Used when method of estimating deamination rate
 		is global.
	-lcf low coverage factor for methylation reconstruction (default: 0.05).
	-sl slope for linear/logistic methods of methylation reconstruction
 		(default: 1/deamination_rate).
	-in intercept for linear method of methylation reconstruction (default: 0).
	-w window size used for DNA methylation reconstruction. Can be
 		--'auto' or a list of one value for each chromosome (default: auto).
	-wm method for computing window size. Can be --'prob' or 'relerror'
 		(default: prob).
	-min parameter for computing window size. The minimum methylation level
 		to detect (default: 0.2).
	-p parameter for computing window size when the method is prob (default: 0.01).
	-k parameter for computing window size when the method is relerror 
		(default: 1/2.5)
	-max maximum window size (default: 31).
	
Outputs:
	
	If the pipeline runs to the final stage, two main outputs are produced: a BED file 
	(<sample>_meth.bed) with the methylation values per position, and a pickled file that 
	can be used for the DMR detection process.
 	The pipeline has five stages: "bam", "diagnose", "filter", "drate", and "meth".
  	These stages can be specified in the config file and from the command line (input -st),
   	allowing the user to choose whether to run the full pipeline, or just several stages.
    	When the script ends (after the last requested stage), it outputs a text file for the
     	last stage completed. This file, in the format <sample>_<stage>.txt, 
	can be found in the specified output directory and can be used as input for later stages. 
    
The stages:
    
    	The first stage, "bam", performs the conversion of BAM file(s) to amSample object(s). It is a 
	prerequisite to the other stages. It can be run by itself or with the other stages, and would
	typically not be be run more than once per sample.
    
	The "diagnose" stage evaluates basic statistics of each input chromosome, and computes
    	recommended values of parameters that are later used (in the "filter" stage) to remove PCR
    	duplicates and positions that are suspected as true mutations. The information is recorded
    	in a text file in the specified log directory, in the format <sample>_diagnostics.txt. It
    	also generates various plots for sanity checks. These are stored as PNG files in the directory
    	specified in the config.ini file.
    
    	The "filter" stage excludes from the analysis CpG sites are suspected to not have gone through
    	deamination, but rather to represent true C->T mutations or to have biased counts due
    	to PCR duplication errors. In this stage, counts from the two strands of each CpG are merged.
    	The details of what sites have been removed are stored in a file in the specified log 
    	directory, in the format <sample>_filter.txt.
    
    	The "drate" stage estimates the deamination rate.
    
    	The "meth" stage carries out the actual reconstruction of the premortem DNA methylation
    	from the C->T ratio, which is simply the quotient #Ts/(#Ts + #Cs).
    
    	When done adjusting the config file, run the run_roam.py script.
    
Warnings:

    	When running run_roam.py, there will be warnings about invalid values and divide by zeros. 
    	They can be safely ignored.
    
# Pipeline 2: Detect DMRs separating two groups of samples

This pipeline receives a list of DNA methylation samples, divided into two groups, and searches
for DMRs between them. Each group can contain ancient samples or modern samples, but not a mixture
of both. 

Preparing modern samples:

	If one of the groups is made of modern samples, these samples should first be
 	converted from COV (Bismark output) format to a format that is appropriate for RoAM.
  	To this end, run the create_files_for_DMRs.py script, using the following input 
	parameters (quoted strings, unless otherwise noted):
	
	-ms sample names, a list specified with no quotes or commas (mandatory).
	-gc CpG coordinates file (mandatory).
	-o directory for saved output (pickled) object files (include final /) (optional) ** what is the meaning of the parenthesis. What is the default? **
	-d directory for data files generated by RoAM (include final /) (optional) ** what is the meaning of the parenthesis. Default? **
	-mt template to match any extra text in modern sample filename (optional). ** I didn't understand **
	-msp modern samples taxonomic classification (optional).
	-mr modern samples reference genome version (optional).
	-mm modern samples sequencing method (optional).
		
	The last three parameters should all be identical for all all modern samples, otherwise the pipeline
 	cannot be run.

Running the DMR detection pipeline:

    	The DMR detection pipeline has a similar flow to the DNA methylation reconstruction pipeline.
     	First edit the parameters in the config_DMR.ini file, or designate them in the command line.
      	These include directory and filenames, samples and group names, parameters for grouping the
	DMRs, and the parts of the script to run.

 	When done adjusting the config file, run the run_DMRs.py script. The script can also be run
 	from the command line, using flags for the required parameters,	as follows:
    
    	run_DMRs.py -s sample1 sample2 -g group1 group2 -gc “CpG file path”
     	(and additional parameter below).
    
Inputs:

	Three input variables are required:
	1. 	Samples (-s): a list of all the samples that are used in the analysis. They
 		should be specified with no quotes or commas. The algorithm will use the pickled
   		files with the same name (or with any extra text specified in the template parameter)
     		from the object directory.
	2.	Groups (-g): a list with identical size to “samples”, with allocation of the samples 
		to the two groups that are compared. For example:
  			-s ust_ishim SF12 Altai_Neanderthal Denisovan -g modern_human modern_human
     									archaic_human archaic_human
    		Each group may contain either ancient samples or modern samples, but not a
      		mixture of both.
	3.	The CpG coordinates file (-gc), identical to this input in the DNA methylation
 		reconstruction pipeline.
	
	Apart from these three mandatory paramters, a modern reference methylation map can be
 	provided, similar to the DNA methylation reconstruction pipeline. This reference is used
  	for the simulations. Lust like in the DNA methylation reconstruction pipeline, this file
   	can be provided as a COV Bismark output (-b) or as a text file (-mo), as described earlier
    	(the bone1 reference modern human bone DNA methylation can be found in
     	http://carmelab.huji.ac.il/data.html). If no reference is provided, RoAM will use the
      	reference defined in the config file of the DNA methylation reconstruction pipeline, and
       	its path can be specified with the parameter -rco. If no such config file exists, simulations
	would not be able to run.
    
    	Other optional parameters can be specified as follows:
     
	-co path of config file (different config files can be used for different runs).
	-rco path of the config file of the DNA methylation reconstruction pipeline.
	-b COV file of modern reference methylation map.
	-mo text file of modern reference methylation map.
	-ms modern reference methylation sample names (list like samples) ** why is this a list?**
	-mn modern reference methylation sample name ** what is the difference between this option and -ms? **
	-msp modern reference methylation sample species.
	-mr modern reference methylation sample genome version.
	-mm modern reference methylation sample sequencing method.
	-c list of chromosome names.
	-o directory for saved output (pickled) object files (include final / for all directories). ** What is in the parenthesis? **
	-du directory for output txt files and pictures.
	-gc CpG cooridantes file.
	-ge sorted text file with genes and their coordinates (for DMR annotation). ** what format? BED? **
	-cg BED file with CpG island coordinates (for DMR annotation). ** what format? BED? **
	-c1 first custom BED file (for DMR annotation).
	-c2 second custom BED file (for DMR annotation).
	-pd promoter definition: a list of two values [before, after], defining the promotoer as
 		running from "before" bases upstream to the TSS, up to "after" bases downstream to it 
		(default: [5000, 1000]).
	-di pickled DMR file for use in fdr, permutations, and plots.
	-t template to match any extra text in sample filename. ** what is this? **
	-st stages of the pipelibe to be run, a list specified with no quotes or commas.
	-de minimum methylation difference between the two groups to be considered as a DMR (default 0.5).
	-mc (min_CpGs) DMRs with less CpGs than this value are filtered out (default: 10).
	-mq (min_qt) DMRs with Qmax lower than this value are filtered out (default: 0).
	-mb (min_bases) DMRs with less than this number of bases are filtered out (default: 100).
	-mad (max_adj_dist) maximum distance between adjacent CpGs within the same DMR. If the 
		distance between consecutive CpG positions is larger than max_adj_dist, the
  		algorithm forces them to belong to different DMRs (default: 1000).
	-mf (min_fin) a list with two elements (one for each group), representing the minimum 
		number of samples for which we require data at a given position. If there are not
  		enough samples with data in a certain position, a NaN is substituted in this position.
    		This parameter can also be a fraction between 0 and 1 (0 <= min_fin < 1), in which
      		case it is understood as the minimum fraction of the total number of samples in the
		group (default: 1 for each group).
	-w window size for smoothing (default: calculated automatically).
	-l low coverage factor, positions with coverage below this value will be filtered out for 
		the corresponding sample (default: calculated automatically).
	-sp number of permutations to run for fdr. ** permutations or simulations? **
	-th FDR threshold (default: 0.05).
	-r flag for using histogram matching in the pooled methylation function (default: TRUE).
	-re flag for logging info, use when logging not desired. ** set FALSE when logging is not required? **
	-an flag for running annotation, use when annotation not desired. ** same as above **
    
The stages:
	This pipeline is comprised of "DMR" and "fdr". (Extra functions that can be used here are permute, 
    	permutstat, and plotmethylation. For more details see “extra functions”).
    
    	The "DMR" stage uses the methylation values in all samples to detect differentially
     	methylated regions. It creates two text files in the specified dump directory
      	(DMR_gen_<timestamp>.txt lists the input parameters and the results, while
       	DMRs_<timestamp>.txt has the data in tab separated format), and a pickled file of
	the DMRs object in the specified object directory (DMR_obj_<timestamp>), for use by
 	subsequent stages and functions. It also outputs a log file in the current dir
  	(DMR_log.txt) with info about the parameters used and DMRs found.
    
    	The "fdr" (false detection rate) stage filters the DMRs to the match a desired FDR level.
     	It also records the thresholds used at the end of the DMR_log file created in the DMR stage.
      	This stage outputs two text files in the dump directory (filtered_DMR_gen_<timestamp>.txt and
    	filtered_DMRs_<timestamp>.txt, as above), and a log file in the current dir with info about 
    	the thresholds and the ratio of the number of observed to the mean number of simulated, DMRs,
     	obtained for each combination of parameters.
    
    	When done adjusting the config file, run the run_DMRs.py script.    
    
Outputs:
    
	The algorithm provides a table with a list of all the DMRs detected and details about 
	them. It also outputs a pickle file that can be used later for the extra functions.

Extra functions:

	We provide several functions along with RoAM, that can be used by users who wish to
 	carry out deeper analysis or to generate plots. To use these functions, change the
  	parameter “stages” to the corresponding term.
	
	1.	The permute function scrambles up the sample group labels and repeat the DMR
 		detection process to obtain “permuted DMRs” that are detected under the null
   		hypothesis of no methylation differences between the groups. The input for the
     		function is the pickled object of the DMR detection process. The number of
       		permutations can be set using the num_permutations parameter in the config
	 	file or in the command line:
		-p number of permutations to run
	
		For each permutation, two text files are created in the dump directory as above
		(DMR_gen_<timestamp>.<permutation #>.txt and DMRs_<timestamp>.<permutation #>.txt),
  		and a log file is created in the current dir (p<permutation #>_groupDMRs.txt) with
    		information on parameters used and DMRs found. This step also outputs a pickled file
      		of the DMR object in the specified object directory (DMP_obj_<timestamp>) for use by
		the permutstat function.
    
    	2.	The permutstat function should follow "permute". It carries out the actual permutation
     		test and outputs the statistical significance of the results of the permute step. It
       		outputs a text file in the dump directory (pstat_<timestamp>.txt) with the resulting
	 	statistics and p-values.
    			-dpi pickled DMR permutation file for use in permutstat
	    
    	3.	The plotmethylation function creates a scatter plot of the mean methylation in each 
   		sample for a specific DMR. It generates a PNG file in the dump directory
     		(meth_plot<chromosome>_<DMR index>.png).
       		The inputs for this function are the pickled file of the detected DMRs, as well as the
	 	unique identifier of the DMR, comprised of the number of chromosome along with the
   		serial number of this DMR among all DMRs in this chromosome (this information can be 
    		found in the DMR table output). If, for example, we wish to plot methylation in the
      		third DMR along chromosome 5, we would use 5 as chromose number and 3 as the DMR index.
			-dmi index of DMR
			-dmc chromosome of DMR
	

  

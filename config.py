# user-supplied vars
name = "sample_name"
# eg: name = "Vindija"
abbrev = "sample_abbrev"
# eg: abbrev = "Vin"
# change these chrom lengths if necessary (these are the lengths for hg19). The list should correspond with
# the list of chromosomes, below.
chr_lengths = [249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,141213431,135534747,135006516,133851895,115169878,107349540,102531392,90354753,81195210,78077248,59128983,63025520,48129895,51304566,155270560,59373566,16571]
#chr_lengths = [51304566,155270560]
library = "strands"  # 'single' or 'double'
# eg: library = "single"
species = "species_name"
# eg: species = "Homo sapiens"
chroms = ["chr"+str(x) for x in range(1,23)] + ["chrX", "chrY"]  # list of chromosomes
# eg: chroms = ["chr22", "chrX"]
trim_ends = False  # if this has already been done (change to True to trim ends during processing)
file_per_chrom = False  # change to True for dir with one file per chrom
# chrom names in filename should match input names in chroms, above

stages = ["bam", "diagnose", "filter", "drate", "meth"]  # stages of process. Change to reflect the stage(s) you want to run.
#  bam is the first stage--conversion of bam file(s) to Amsample object. It is a prerequisite to the other stages.
#  It can be run by itself or with the other stages, but need not be run more than once.
#  When running in stages, the last part will be dumped to a text file in the specified directory, below.
#  Text files will be identified by object name and stage. When running later stages without the first part,
#  specify the text file name in the text_infile var, below.

# input files
filedir = "mult_file_dir"  # the dir for multiple files, or None

filename = "path_to_bam_file"  # the path for the bam file (if files in filedir, above, filename = None)

genome_file = "path_to_seq_file"  # the path for the assembly file, eg hg19.fa.gz
# eg: genome_file = "hg19.fa.gz"
text_infile = "path_to_object_textfile"  # the path for the object saved in text format
# eg: text_infile = "ams_vin.txt"
modern_infile = "bone5.txt"  # sample text file for modern genome
bismark_infile = "Bismark file"  # this file is a .cov or .bedGraph (or None if using sample text file)
# eg: bismark_infile = "sample5.bedGraph"
gc_object = "objects/cpg_coords.P"  # CpG file
mod_name = "Bone 5"  # sample (individual) name for modern genome
mod_abbrev = "Bone5"  # sample abbreviation
mod_spec = "Homo sapiens"  # sample species
mod_ref = "hg19"  # sample reference genome
mod_method = "WGBS"  # sample sequencing method

# output
outdir = "dir_for_textfiles"
# eg: outdir = "data/python_dumps/"
logdir = "dir for logfiles"
# eg: logdir = "data/logs/"
picdir = "dir for output figures"
# eg: picdir = "data/figures/"

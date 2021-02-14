# user-supplied vars
name = "sample_name"  # eg "Ust_Ishim"
abbrev = "sample_abbrev"  # eg "Ust"
# change these chrom lengths if nec
chr_lengths = [249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,141213431,135534747,135006516,133851895,115169878,107349540,102531392,90354753,81195210,78077248,59128983,63025520,48129895,51304566,155270560,59373566]
library = "strands"  # single or double
species = "species_name"  # eg Homo sapiens
chroms = list(range(23))  # list of chrom indices
trim_ends = False  # if this has already been done (change to True to trim ends during processing)
file_per_chrom = False  # change to True for dir with one file per chrom

# input files
filedir = "mult_file_dir"  # the dir for multiple files, or None
filename = "path_to_bam_file"  # the path for the bam file
genome_file = "path_to_seq_file"  # the path for the assembly file, eg hg19.fa.gz
text_infile = "path_to_object_textfile"  # the path for the object saved in text format
pickle_infile = "path_to_pickled_object"  # the path for the pickled object

# output

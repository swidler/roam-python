# dirs
object_dir = "directory for pickled objects"
# eg object_dir = "objects/"
dump_dir = "directory for output txt files and pics"
# eg dump_dir = "data/python_dumps/DMR/"

# files
gc_object = object_dir + "cpg_coords.P"  # CpG file
gene_file = "../../UCSC_Genes_sorted.txt"  # text file  sorted on the cmdline: sort -k1,1n -k2,2n -k3,3n fn>fn_sorted
cgi_file = "../../UCSC_CGIs.txt"  # CGI file
DMR_obj_infile = "pickled DMR file for use in permutations and plots"
# eg DMR_obj_infile = object_dir + "DMR_obj_27-02-2022_10.33"
DMP_obj_infile = "pickled DMR permutation file for use in permutstat"
# eg DMP_obj_infile = object_dir + "DMP_obj27-02-2022_12.21"

# filename vars
templ = "template to match any extra text in sample file"
# eg templ = "filt_"

scripts = ["DMR", "permute", "permutstat", "plotmethylation"] # scripts in the process. Change to reflect what you want to run.
# DMR runs the groupDMR process and the annotation. It can be run independently or in conjunction with permutation scripts.
# Permutation scripts (permute runs permutations, permutstat calculates permutation statistics) can be run independently as 
# long as they have a DMR file to work with. The plotmethylation script also requires a DMR file.


# groupDMRs vars
# change sample and group names to reflect your samples and their file names
samples = ["1116","1496","2520","4873","4875","4877","4878","4914","5077","5233","5235","5236"]
group_names = ["Farmers","Farmers","Farmers","HGs","HGs","HGs","HGs","HGs","Farmers","HGs","HGs","HGs"]
delta = 0.15
min_CpGs = 40

# permute vars
num_permutations = 10

# plotmethylation vars
# index number of desired DMR (indexing is zero-based!) and chromosome name
DMR_idx = 3
DMR_chrom = "chr11"

    
# dirs
object_dir = ""  # directory for pickled objects--defaults to ""
# eg object_dir = "objects/"
data_dir = ""  # directory for data files from RoAM process--defaults to ""
# eg data_dir = "data/python_dumps/"
dump_dir = ""  # directory for output txt files and pics--defaults to ""
# eg dump_dir = data_dir + "DMR/"


# files
gc_object = object_dir + "cpg_coords.P"  # CpG file
gene_file = "UCSC_Genes_sorted.txt"  # text file  sorted on the cmdline: sort -k1,1n -k2,2n -k3,3n fn>fn_sorted
cgi_file = "UCSC_CGIs.txt"  # CGI file
DMR_obj_infile = "pickled DMR file for use in permutations and plots"
# eg DMR_obj_infile = object_dir + "DMR_obj_27-02-2022_10.33"
DMP_obj_infile = "pickled DMR permutation file for use in permutstat"
# eg DMP_obj_infile = object_dir + "DMP_obj27-02-2022_12.21"

# filename vars
#templ = "template to match any extra text in sample file"
# eg 
templ = "_meth"

stages = ["create_files", "DMR", "fdr", "permute", "permutstat", "plotmethylation"] # stages in the process. Change to reflect what you want to run.
# DMR runs the groupDMR process and the annotation. It can be run independently or in conjunction with permutation stages.
# Permutation stages (permute runs permutations, permutstat calculates permutation statistics) can be run independently as 
# long as they have a DMR file to work with. The plotmethylation method also requires a DMR file.
ref = True  # for histogram matching. When ref is True, histogram matching is done (ref used is bone5 by default, and changeable in modern file params)

# groupDMRs vars
# change sample and group names to reflect your samples and their file names
samples = ["1116","1496","2520","4873","4875","4877","4878","4914","5077","5233","5235","5236"]
group_names = ["Farmers","Farmers","Farmers","HGs","HGs","HGs","HGs","HGs","Farmers","HGs","HGs","HGs"]
delta = 0.5
min_CpGs = 40
min_fin = [1,1]
min_Qt = 0
win_size = "meth"  # for pre-methylation-stage samples (eg DMR), change to "auto"
lcf = "meth"  # for pre-methylation-stage samples (eg DMR), change to value (eg 0.05)

# permute vars
num_permutations = 10
sim_permutations = 5
# plotmethylation vars
# index number of desired DMR (indexing is zero-based!) and chromosome name
DMR_idx = 0
DMR_chrom = "chr1"

    
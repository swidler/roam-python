# dirs
object_dir = "objects/"
dump_dir = "data/python_dumps/"

# files
gc_object = object_dir + "cpg_coords.P"  # CpG file
gene_file = "../../UCSC_Genes_sorted.txt"  # text file  sorted on the cmdline: sort -k1,1n -k2,2n -k3,3n fn>fn_sorted
cgi_file = "../../UCSC_CGIs.txt"
DMR_obj_infile = object_dir + "DMR_obj_27-02-2022_10.33"
DMP_obj_infile = object_dir + "DMP_obj27-02-2022_12.21"

# filename vars
templ = "filt_"
num = 1

scripts = ["DMR", "permute", "permutstat", "plotmethylation"]
scripts = ["DMR", "permutstat"]

# groupDMRs vars
samples = ["1116","1496","2520","4873","4875","4877","4878","4914","5077","5233","5235","5236"]
group_names = ["Farmers","Farmers","Farmers","HGs","HGs","HGs","HGs","HGs","Farmers","HGs","HGs","HGs"]
delta = 0.15
min_CpGs = 40

# permute vars
num_permutations = 3

# plotmethylation vars
DMR_idx = 3
DMR_chrom = "chr11"


    
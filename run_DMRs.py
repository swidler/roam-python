#!/usr/bin/python3

import DMRs as d
import tools as t

scripts = ["DMR", "permute", "permutstat"]
scripts = ["permutstat"]
if "DMR" in scripts or "permute" in scripts:
    infile = "objects/cpg_coords"
    gc = t.load_object(infile)
    chr_names = gc.chr_names 
    chr_names = [x for x in chr_names if "X" not in x]  # remove chrx from list
    samples = ["1116","1496","2520","4873","4875","4877","4878","4914","5077","5233","5235","5236"]
    group_names = ["Farmers","Farmers","Farmers","HGs","HGs","HGs","HGs","HGs","Farmers","HGs","HGs","HGs"]
    delta = 0.15
    min_CpGs = 40
    samplist = []
    templ = "filt_"
    #for the next step, all input files must be pickled
    indir = "objects/"
    for sample in samples:
        infile = indir + templ + sample
        print(f"loading sample I{sample}")
        input_obj = t.load_object(infile)
        samplist.append(input_obj)
#print("loading modern")
#mod = "ref_bone5_take2"
##samples.append(mod)
#group_names.append("Farmers")
#infile = indir + mod
#input_obj = t.load_object(infile)
#input_obj.merge(False)
#samplist.append(input_obj)
    print(samples)
if "DMR" in scripts:
    dms = d.DMRs()
else:
    infile = "objects/DMR_obj"
    dms = t.load_object(infile)
if "DMR" in scripts:
    (qt_up, qt_down) = dms.groupDMRs(samples=samplist, sample_groups=group_names, coord=gc, chroms=chr_names, min_finite=[0.8,0.8], min_CpGs=min_CpGs, delta=delta)
    gene_file = "../../UCSC_Genes_sorted.txt"  # text file  sorted on the cmdline: sort -k1,1n -k2,2n -k3,3n fn>fn_sorted
    cgi_file = "../../UCSC_CGIs.txt"
    dms.annotate(gene_file, cgi_file)
    dms.dump_DMR()
    t.save_object(f"{indir}/DMR_obj", dms)
if "permute" in scripts:
    num_permutations = 10
    dmp = dms.permute(num_permutations, samplist, gc)
    t.save_object(f"{indir}/DMP_obj", dmp)
if "permutstat" in scripts:
    if "permute" not in scripts:
        infile = "objects/DMP_obj"
        dmp = t.load_object(infile)
    pstat = dms.permutstat(dmp)
    dms.dump_pstat(pstat)
    
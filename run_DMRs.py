#!/usr/bin/python3

import DMRs as d
import tools as t
import datetime

scripts = ["DMR", "permute", "permutstat"]  # currently hardcoded vars should be in config file
scripts = ["DMR"]
if "DMR" in scripts or "permute" in scripts:
    infile = "objects/cpg_coords.P"
    gc = t.load_object(infile)
    chr_names = gc.chr_names 
    chr_names = [x for x in chr_names if "X" not in x]  # remove chrx from list
    #chr_names = chr_names[6:11]
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
    print(samples)
if "DMR" in scripts:
    dms = d.DMRs()
else:
    infile = "objects/DMR_obj" #use fn from config
    dms = t.load_object(infile)
if "DMR" in scripts:
    import cProfile
    #cProfile.run("(qt_up, qt_down) = dms.groupDMRs(samples=samplist, sample_groups=group_names, coord=gc, chroms=chr_names, min_finite=[0.8,0.8], min_CpGs=min_CpGs, delta=delta)", "data/logs/DMR_profile")
    (qt_up, qt_down) = dms.groupDMRs(samples=samplist, sample_groups=group_names, coord=gc, chroms=chr_names, min_finite=[0.8,0.8], min_CpGs=min_CpGs, delta=delta)
    gene_file = "../../UCSC_Genes_sorted.txt"  # text file  sorted on the cmdline: sort -k1,1n -k2,2n -k3,3n fn>fn_sorted
    cgi_file = "../../UCSC_CGIs.txt"
    dms.annotate(gene_file, cgi_file)
    #cProfile.run("dms.annotate(gene_file, cgi_file)", "data/logs/annotate_profile")
    dms.dump_DMR()
    time = datetime.datetime.now()
    time = time.strftime("%d-%m-%Y_%H.%M")
    t.save_object(f"{indir}/DMR_obj_{time}", dms)  # save more specifically?
if "permute" in scripts:
    num_permutations = 10  # use var from config
    dmp = dms.permute(num_permutations, samplist, gc)
    i = 1
    for d in dmp:  # bad use of d. change and test
        d.dump_DMR(num=i)
        i += 1
    t.save_object(f"{indir}/DMP_obj", dmp)
if "permutstat" in scripts:
    if "permute" not in scripts:
        infile = "objects/DMP_obj"
        dmp = t.load_object(infile)
    pstat = dms.permutstat(dmp)
    dms.dump_pstat(pstat)
if "plotmethylation" in scripts:
    chrom = "chr11"  # use vars from config
    DMR_idx = 3
    dms.plotmethylation(chrom, DMR_idx)  # add file name to save pic
    
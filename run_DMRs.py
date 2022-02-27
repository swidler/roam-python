#!/usr/bin/python3

import DMRs as d
import tools as t
import datetime
from config_DMR import *

time = datetime.datetime.now()
time = time.strftime("%d-%m-%Y_%H.%M")
if "DMR" in scripts or "permute" in scripts:
    gc = t.load_object(gc_object)
    chr_names = gc.chr_names 
    chr_names = [x for x in chr_names if "X" not in x]  # remove chrx from list
    samplist = []
    #for the next step, all input files must be pickled
    for sample in samples:
        infile = object_dir + templ + sample
        print(f"loading sample I{sample}")
        input_obj = t.load_object(infile)
        samplist.append(input_obj)
    print(samples)
if "DMR" in scripts:
    dms = d.DMRs()
else:
    dms = t.load_object(DMR_obj_infile)
if "DMR" in scripts:
    import cProfile
    #cProfile.run("(qt_up, qt_down) = dms.groupDMRs(samples=samplist, sample_groups=group_names, coord=gc, chroms=chr_names, min_finite=[0.8,0.8], min_CpGs=min_CpGs, delta=delta)", "data/logs/DMR_profile")
    (qt_up, qt_down) = dms.groupDMRs(samples=samplist, sample_groups=group_names, coord=gc, chroms=chr_names, min_finite=[0.8,0.8], min_CpGs=min_CpGs, delta=delta)
    dms.annotate(gene_file, cgi_file)
    #cProfile.run("dms.annotate(gene_file, cgi_file)", "data/logs/annotate_profile")
    fname = dump_dir + f"DMRs_{time}.txt"
    dms.dump_DMR(fname)
    t.save_object(f"{object_dir}DMR_obj_{time}", dms)  
if "permute" in scripts:
    dmp = dms.permute(num_permutations, samplist, gc)
    i = 1
    for x in dmp:  
        fname = dump_dir + f"DMRs_{time}.{i}.txt"
        x.dump_DMR(fname)
        i += 1
    t.save_object(f"{object_dir}DMP_obj{time}", dmp) 
if "permutstat" in scripts:
    if "permute" not in scripts:
        dmp = t.load_object(DMP_obj_infile)
    pstat = dms.permutstat(dmp)
    fname = dump_dir + f"pstat_{time}.txt"
    dms.dump_pstat(pstat, fname)
if "plotmethylation" in scripts:
    fname = dump_dir + f"meth_plot_{DMR_chrom}_{DMR_idx}"
    dms.plotmethylation(DMR_chrom, DMR_idx, fname)  
    
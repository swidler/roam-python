#!/usr/bin/python3

import DMRs as d
import amsample as a
import tools as t
import datetime
from config_DMR import *

time = datetime.datetime.now()
time = time.strftime("%d-%m-%Y_%H.%M")
if "create_files" in stages:
    for sample in samples:
        ams = a.Amsample()
        filename = data_dir + sample + templ + ".txt"
        ams.parse_infile(filename)
        outfile = object_dir + sample + templ
        t.save_object(outfile, ams)
if "DMR" in stages or "permute" in stages:
    gc = t.load_object(gc_object)
    chr_names = gc.chr_names 
    chr_names = [x for x in chr_names if "X" not in x]  # remove chrx from list
    samplist = []
    #for the next step, all input files must be pickled
    for sample in samples:
        infile = object_dir + sample + templ
        print(f"loading sample {sample}")
        input_obj = t.load_object(infile)
        samplist.append(input_obj)
    print(samples)
if "DMR" in stages:
    dms = d.DMRs()
else:
    dms = t.load_object(DMR_obj_infile)
if "DMR" in stages:
    import cProfile
    #cProfile.run("(qt_up, qt_down) = dms.groupDMRs(samples=samplist, sample_groups=group_names, coord=gc, chroms=chr_names, min_finite=[0.8,0.8], min_CpGs=min_CpGs, delta=delta)", "data/logs/DMR_profile")
    (qt_up, qt_down) = dms.groupDMRs(samples=samplist, sample_groups=group_names, coord=gc, chroms=chr_names, min_finite=[0.8,0.8], min_CpGs=min_CpGs, delta=delta)
    dms.annotate(gene_file, cgi_file)
    #cProfile.run("dms.annotate(gene_file, cgi_file)", "data/logs/annotate_profile")
    fname = dump_dir + f"DMRs_{time}.txt"
    dms.dump_DMR(fname)
    t.save_object(f"{object_dir}DMR_obj_{time}", dms)  
if "permute" in stages:
    dmp = dms.permute(num_permutations, samplist, gc)
    i = 1
    for x in dmp:  
        fname = dump_dir + f"DMRs_{time}.{i}.txt"
        x.dump_DMR(fname)
        i += 1
    t.save_object(f"{object_dir}DMP_obj{time}", dmp) 
if "permutstat" in stages:
    if "permute" not in stages:
        dmp = t.load_object(DMP_obj_infile)
    pstat = dms.permutstat(dmp)
    fname = dump_dir + f"pstat_{time}.txt"
    dms.dump_pstat(pstat, fname)
if "plotmethylation" in stages:
    fname = dump_dir + f"meth_plot_{DMR_chrom}_{DMR_idx}"
    dms.plotmethylation(DMR_chrom, DMR_idx, fname)  
    
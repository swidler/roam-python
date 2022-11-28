#!/usr/bin/python3

import DMRs as d
import amsample as a
import mmsample as m
import tools as t
import datetime
from config_DMR import *  # import as cfg for disambiguation?

time = datetime.datetime.now()
time = time.strftime("%d-%m-%Y_%H.%M")
if "create_files" in stages:
    for sample in samples:
        ams = a.Amsample()
        filename = data_dir + sample + templ + ".txt"
        ams.parse_infile(filename)
        outfile = object_dir + sample + templ
        t.save_object(outfile, ams)
if "DMR" in stages or "permute" in stages or "plot" in stages:
    gc = t.load_object(gc_object)
    chr_names = gc.chr_names 
    chr_names = [x for x in chr_names if "X" not in x]  # remove chrx from list
    samplist = []
    #for the next step, all input files must be pickled
    for sample in samples:
        infile = object_dir + templ + sample
        #infile = object_dir + sample + templ
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
    min_finite = min_fin[:]
    #cProfile.run("(qt_up, qt_down) = dms.groupDMRs(samples=samplist, sample_groups=group_names, coord=gc, chroms=chr_names, min_finite=[0.8,0.8], min_CpGs=min_CpGs, delta=delta)", "data/logs/DMR_profile")
    (qt_up, qt_down) = dms.groupDMRs(win_size=win_size, lcf=lcf, samples=samplist, sample_groups=group_names, coord=gc, chroms=chr_names, min_finite=min_finite, min_CpGs=min_CpGs, delta=delta)
    #dms.annotate(gene_file, cgi_file)  # move annotation to after fdr? or call here only when dmr, but not fdr in stages.
    #cProfile.run("dms.annotate(gene_file, cgi_file)", "data/logs/annotate_profile")
    fname = dump_dir + f"DMRs_{time}.txt"
    #fname = dump_dir + "python_DMRs.txt"
    #with open(fname, "w") as fid:
    #    for chrom in range(len(dms.chromosomes)):
    #        for dmr in range(len(dms.cDMRs[chrom].gen_start)):
    #            fid.write(f"{chrom+1}\t{dms.cDMRs[chrom].gen_start[dmr]}\t{dms.cDMRs[chrom].gen_end[dmr]}\n")

    #dms.dump_DMR(fname)  # dump doesn't currently work w/out annotate
    t.save_object(f"{object_dir}DMR_obj_{time}", dms) 
if "fdr" in stages:
    #create Mmsample object
    mms = m.Mmsample()
    if cfg.bismark_infile:
        mms.create_mms_from_bismark_file()
    else:
        mms.create_mms_from_text_file()
    samplist = []  # if dmr in stages, samplist already loaded
    for sample in samples:
        #use filtered files
        infile = object_dir + templ + sample
        #infile = object_dir + sample + templ
        print(f"loading sample {sample}")
        input_obj = t.load_object(infile)
        samplist.append(input_obj)
    sim_obj_list = []
    dmr_obj_list = []
    for perm in range(sim_permutations):
        sim_obj_list.append({})
        dmr_obj_list.append(d.DMRs())
    for sample in samplist:
        for perm in range(sim_permutations):
            samp_name = sample.name
            sample.simulate(mms)
            sample.estimate_drate(ref=mms)
            sample.reconstruct_methylation(ref=mms)
            #dump object to text file
            sim_obj_list[perm][samp_name] = sample
            #sample.dump(f"simeth_{perm}")  
    gc = t.load_object(gc_object)
    chr_names = gc.chr_names 
    chr_names = [x for x in chr_names if "X" not in x]  # remove chrx from list
    for perm in range(sim_permutations):
        min_finite = min_fin[:]
        dmr_obj_list[perm].groupDMRs(win_size=win_size, lcf=lcf, samples=list(sim_obj_list[perm].values()), sample_groups=group_names, coord=gc, chroms=chr_names, min_finite=min_finite, min_CpGs=min_CpGs, delta=delta)
    filtered_DMR = dms.run_fdr(dmr_obj_list)
    print("done")
    
        
    
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
if "plot" in stages:
    dms.plot(DMR_chrom, DMR_idx, gc, samplist, gene_file, cgi_file, widenby=5000)
    
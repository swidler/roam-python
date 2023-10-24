#!/usr/bin/python3

import DMRs as d
import amsample as a
import mmsample as m
import tools as t
import datetime
import copy
import config_DMR as cfg  
import config as rcfg
import argparse
import sys
import os

# params from config can be specified on the command line
argParser = argparse.ArgumentParser()
argParser.add_argument("-s", "--samples", nargs="+", help="sample names")
argParser.add_argument("-ms", "--mod_samples", nargs="+", help="modern sample names")
argParser.add_argument("-g", "--groups", nargs="+", help="group names--should correspond with samples")
argParser.add_argument("-o", "--object_dir", help="directory for pickled objects")
argParser.add_argument("-d", "--data_dir", help="directory for data files from RoAM process")
argParser.add_argument("-du", "--dump_dir", help="directory for output txt files and pics")
argParser.add_argument("-gc", "--gc_file", help="CpG file")
argParser.add_argument("-ge", "--gene_file", help="sorted text file with genes")
argParser.add_argument("-cg", "--cgi_file", help="CGI file")
argParser.add_argument("-di", "--dmr_infile", help="pickled DMR file for use in fdr, permutations, and plots")
argParser.add_argument("-dpi", "--dmp_infile", help="pickled DMR permutation file for use in permutstat")
argParser.add_argument("-t", "--templ", help="template to match any extra text in sample filename")
argParser.add_argument("-st", "--stages", nargs="+", help="stages of process to be run")
argParser.add_argument("-de", "--delta", type=float, help="minimum methylation difference between the two groups")
argParser.add_argument("-mc", "--min_cpgs", type=int, help="DMRs whose number of CpGs is less than min_CpGs are filtered out")
argParser.add_argument("-mq", "--min_qt", type=int, help="DMRs with Qt < min_qt are filtered out")
argParser.add_argument("-mf", "--min_finite", nargs="+", help="minimum number of ancient samples for which we require data")
argParser.add_argument("-w", "--win_size", help="window size for smoothing")
argParser.add_argument("-l", "--lcf", help="low coverage factor")
argParser.add_argument("-p", "--permutations", type=int, help="number of permutations to run")
argParser.add_argument("-dmi", "--dmr_idx", type=int, help="index of DMR")
argParser.add_argument("-dmc", "--dmr_chrom", help="chromosome of DMR")
argParser.add_argument("-b", "--bismark", help=".cov or .bedGraph file for modern genome")
argParser.add_argument("-mo", "--modern", help="text file for modern genome")
argParser.add_argument("-r", "--ref", help="reference genome for use in histogram matching")
argParser.add_argument("-re", "--noreport", help="flag for logging info", action="store_true")

args = argParser.parse_args()
keys = [x for x in vars(args).keys() if vars(args)[x] != None]
vals = [vars(args)[x] for x in keys]
parameters = dict(zip(keys, vals))
samples = parameters["samples"] if "samples" in parameters else cfg.samples
mod_samples = parameters["mod_samples"] if "mod_samples" in parameters else cfg.mod_samples
templ = parameters["templ"] if "templ" in parameters else cfg.templ
object_dir = parameters["object_dir"] if "object_dir" in parameters else cfg.object_dir
stages = parameters["stages"] if "stages" in parameters else cfg.stages
gc_object = parameters["gc_file"] if "gc_file" in parameters else cfg.gc_object
min_fin = parameters["min_finite"] if "min_finite" in parameters else cfg.min_fin
win_size = parameters["win_size"] if "win_size" in parameters else cfg.win_size
lcf = parameters["lcf"] if "lcf" in parameters else cfg.lcf
group_names = parameters["groups"] if "groups" in parameters else cfg.group_names
min_CpGs = parameters["min_cpgs"] if "min_cpgs" in parameters else cfg.min_CpGs
min_Qt = parameters["min_qt"] if "min_qt" in parameters else cfg.min_Qt
delta = parameters["delta"] if "delta" in parameters else cfg.delta
bismark_infile = parameters["bismark"] if "bismark" in parameters else rcfg.bismark_infile
modern = parameters["modern"] if "modern" in parameters else rcfg.modern_infile
gene_file = parameters["gene_file"] if "gene_file" in parameters else cfg.gene_file
cgi_file = parameters["cgi_file"] if "cgi_file" in parameters else cfg.cgi_file
dump_dir = parameters["dump_dir"] if "dump_dir" in parameters else cfg.dump_dir
report = False if parameters["noreport"] else cfg.report
# add param for permutations in permute?

time = datetime.datetime.now()
time = time.strftime("%d-%m-%Y_%H.%M")
if "create_ancient_files" in stages:
    data_dir = parameters["data_dir"] if "data_dir" in parameters else cfg.data_dir
    for sample in samples:
        ams = a.Amsample()
        filename = data_dir + sample + templ + ".txt"
        ams.parse_infile(filename)
        outfile = object_dir + sample + templ
        t.save_object(outfile, ams)
if "create_modern_files" in stages:
    data_dir = parameters["data_dir"] if "data_dir" in parameters else cfg.data_dir
    for sample in mod_samples:
        mms = m.Mmsample()
        bisfile = ""
        filename = data_dir + sample + ".bedGraph"
        if os.path.isfile(filename):
            bisfile = filename
        else:
            filename = data_dir + sample + ".cov"
            if os.path.isfile(filename):
                bisfile = filename
        if bisfile:
            mms.bismark_to_mm(bisfile, gc_object, sample, rcfg.mod_abbrev, rcfg.mod_spec, "", rcfg.mod_method)
            outfile = object_dir + sample + templ
            t.save_object(outfile, mms)
        else:
            raise Exception(f"No file in {data_dir} matches {sample}")
if "DMR" in stages or "permute" in stages or "plot" in stages:
    gc = t.load_object(gc_object)
    chr_names = gc.chr_names  # assumes user wants all chroms (or all but x)
    chr_names = [x for x in chr_names if "X" not in x]  # remove chrx from list
    samplist = []
    #for the next step, all input files must be pickled
    for samps in (samples, mod_samples):
        for sample in samps:
            infile = object_dir + sample + templ
            print(f"loading sample {sample}")
            input_obj = t.load_object(infile)
            samplist.append(input_obj)
    print(samples)  # print mod_samples too
if "DMR" in stages:
    dms = d.DMRs()
elif "fdr" in stages or "permute" in stages or "permutstat" in stages or "plotmethylation" in stages or "plot" in stages:
    DMR_obj_infile = parameters["dmr_infile"] if "dmr_infile" in parameters else cfg.DMR_obj_infile
    dms = t.load_object(DMR_obj_infile)
if "DMR" in stages or "fdr" in stages:
    mms = m.Mmsample()
    if bismark_infile:
        mms.create_mms_from_bismark_file(bismark_infile, gc_object, rcfg.mod_name, rcfg.mod_abbrev, rcfg.mod_spec, rcfg.mod_ref, rcfg.mod_method)
    else:
        mms.create_mms_from_text_file(modern)
if "DMR" in stages:
    import cProfile
    ref = parameters["ref"] if "ref" in parameters else cfg.ref
    if ref:
        ref = mms
    min_finite = min_fin[:]
    #cProfile.run("(qt_up, qt_down) = dms.groupDMRs(samples=samplist, sample_groups=group_names, coord=gc, chroms=chr_names, min_finite=[0.8,0.8], min_CpGs=min_CpGs, delta=delta)", "data/logs/DMR_profile")
    (qt_up, qt_down) = dms.groupDMRs(win_size=win_size, lcf=lcf, samples=samplist, sample_groups=group_names, coord=gc, chroms=chr_names, min_finite=min_finite, min_CpGs=min_CpGs, delta=delta, ref=ref)
    #cProfile.run("dms.annotate(gene_file, cgi_file)", "data/logs/annotate_profile")
    
    t.save_object(f"{object_dir}DMR_obj_{time}", dms) 
    if report:
        dms.annotate(gene_file, cgi_file)  
        fname = dump_dir + f"DMRs_{time}.txt"
        dms.dump_DMR(fname)  # dump doesn't currently work w/out annotate
if "fdr" in stages:
    #create Mmsample object
    sim_permutations = parameters["permutations"] if "permutations" in parameters else cfg.sim_permutations
    
    samplist = []  # if dmr in stages, samplist already loaded
    for sample in samples:  
    #for samps in (samples, mod_samples):
        #for sample in samps:
        #use filtered files
        #infile = object_dir + templ + sample
        infile = object_dir + sample + templ
        print(f"loading sample {sample}")
        input_obj = t.load_object(infile)
        samplist.append(input_obj)
    mod_samplist = []
    for sample in mod_samples:
        infile = object_dir + sample + templ
        print(f"loading sample {sample}")
        input_obj = t.load_object(infile)
        mod_samplist.append(input_obj)
    dmr_obj_list = []
    gc = t.load_object(gc_object)
    chr_names = gc.chr_names  # assumes user wants all chroms (or all but x)
    chr_names = [x for x in chr_names if "X" not in x]  # remove chrx from list
    for perm in range(sim_permutations):
        sim_obj_list = {}
        dmr_obj_list.append(d.DMRs())
        for sample in samplist:
            print(f"sample {sample.name}, perm {perm}\n")
            samp_name = sample.name
            samp_copy = copy.deepcopy(sample)
            samp_copy.simulate(mms)
            samp_copy.estimate_drate(ref=mms)
            samp_copy.reconstruct_methylation(ref=mms)
            #dump object to text file
            sim_obj_list[samp_name] = samp_copy
            del samp_copy  # will this fix memory errors?
            #sample.dump(f"simeth_{perm}")  
        min_finite = min_fin[:]
        dmr_obj_list[perm].groupDMRs(win_size=win_size, lcf=lcf, samples=list(sim_obj_list.values(), mod_samplist), sample_groups=group_names, coord=gc, chroms=chr_names, min_finite=min_finite, min_CpGs=min_CpGs, delta=delta)
    adjusted_DMR = dms.adjust_params(dmr_obj_list)
    #if not ("DMR" in stages and report):
    adjusted_DMR.annotate(gene_file, cgi_file)  
    fname = dump_dir + f"filtered_DMRs_{time}.txt"
    adjusted_DMR.dump_DMR(fname)  # dump doesn't currently work w/out annotate
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
        DMP_obj_infile = parameters["dmp_infile"] if "dmp_infile" in parameters else cfg.DMP_obj_infile
        dmp = t.load_object(DMP_obj_infile)
    pstat = dms.permutstat(dmp)
    fname = dump_dir + f"pstat_{time}.txt"
    dms.dump_pstat(pstat, fname)
    
DMR_chrom = parameters["dmr_chrom"] if "dmr_chrom" in parameters else cfg.DMR_chrom
DMR_idx = parameters["dmr_idx"] if "dmr_idx" in parameters else cfg.DMR_idx
    
if "plotmethylation" in stages:
    fname = dump_dir + f"meth_plot_{DMR_chrom}_{DMR_idx}"
    dms.plotmethylation(DMR_chrom, DMR_idx, fname)  
if "plot" in stages:
    dms.plot(DMR_chrom, DMR_idx, gc, samplist, gene_file, cgi_file, widenby=5000)
    
#!/usr/bin/python3

import DMRs as d
import mmsample as m
import tools as t
import datetime
import copy
import argparse
import sys
import os
import configparser as cp

# params from config can be specified on the command line
argParser = argparse.ArgumentParser()
argParser.add_argument("-co", "--config", help="path of config file")
argParser.add_argument("-rco", "--rconfig", help="path of RoAM config file")
argParser.add_argument("-s", "--samples", nargs="+", help="sample names")
argParser.add_argument("-ms", "--mod_samples", nargs="+", help="modern sample names")
argParser.add_argument("-c", "--chroms", nargs="+", help="list of chromosomes")
argParser.add_argument("-g", "--groups", nargs="+", help="group names--should correspond with samples")
argParser.add_argument("-o", "--object_dir", help="directory for saved (pickled) object files (include final /)")
argParser.add_argument("-du", "--dump_dir", help="directory for output txt files and pics")
argParser.add_argument("-lo", "--log_dir", help="directory for logging txt files")
argParser.add_argument("-gc", "--gc_file", help="CpG file")
argParser.add_argument("-ge", "--gene_file", help="sorted text file with genes")
argParser.add_argument("-cg", "--cgi_file", help="CGI file")
argParser.add_argument("-c1", "--cust_file1", help="first custom bed file")
argParser.add_argument("-c2", "--cust_file2", help="second custom bed file")
argParser.add_argument("-pd", "--prom_def", nargs="+", help="promoter definition around TSS: a list of 2 values [before, after], where before is the number of nucleotides into the intergenic region, and after is the number of nucleotides")
argParser.add_argument("-di", "--dmr_infile", help="pickled DMR file for use in fdr, permutations, and plots")
argParser.add_argument("-dpi", "--dmp_infile", help="pickled DMR permutation file for use in permutstat")
argParser.add_argument("-t", "--templ", help="template to match any extra text in sample filename")
argParser.add_argument("-st", "--stages", nargs="+", help="stages of process to be run")
argParser.add_argument("-de", "--delta", type=float, help="minimum methylation difference between the two groups")
argParser.add_argument("-mc", "--min_cpgs", type=int, help="DMRs whose number of CpGs is less than min_CpGs are filtered out")
argParser.add_argument("-mq", "--min_qt", type=int, help="DMRs with Qt < min_qt are filtered out")
argParser.add_argument("-mb", "--min_bases", type=int, help="DMRs shorter than min_bases are filtered out")
argParser.add_argument("-mad", "--max_adj_dist", type=int, help="max distance between adjacent CpGs within the same DMR")
argParser.add_argument("-mf", "--min_finite", nargs="+", help="minimum number of ancient samples for which we require data")
argParser.add_argument("-w", "--win_size", help="window size for smoothing")
argParser.add_argument("-l", "--lcf", help="low coverage factor")
argParser.add_argument("-sp", "--sim_permutations", type=int, help="number of permutations to run for fdr")
argParser.add_argument("-th", "--thresh", type=float, help="FDR threshold")
argParser.add_argument("-p", "--permutations", type=int, help="number of permutations to run")
argParser.add_argument("-dmi", "--dmr_idx", type=int, help="index of DMR")  # index or name??
argParser.add_argument("-dmc", "--dmr_chrom", help="chromosome of DMR")
argParser.add_argument("-b", "--bismark", help=".cov file for modern genome")
argParser.add_argument("-mo", "--modern", help="text file for modern genome")
argParser.add_argument("-mn", "--mname", help="modern reference sample name")
argParser.add_argument("-msp", "--mspecies", help="modern reference sample species")
argParser.add_argument("-mr", "--mref", help="modern sample reference genome")
argParser.add_argument("-mm", "--mmethod", help="modern reference sample sequencing method")
argParser.add_argument("-r", "--noref", help="flag for using reference genome for histogram matching", action="store_true")
argParser.add_argument("-re", "--noreport", help="flag for logging info", action="store_true")
argParser.add_argument("-an", "--noannot", help="flag for running annotation", action="store_true")

args = argParser.parse_args()
keys = [x for x in vars(args).keys() if vars(args)[x] != None]
vals = [vars(args)[x] for x in keys]
parameters = dict(zip(keys, vals))
confile = parameters["config"] if "config" in parameters else "config_DMR.ini"
config = cp.ConfigParser(interpolation=cp.ExtendedInterpolation())
config.read(confile)
rconfile = parameters["rconfig"] if "rconfig" in parameters else "config.ini"
rconfig = cp.ConfigParser(interpolation=cp.ExtendedInterpolation())
rconfig.read(rconfile)

samples = parameters["samples"[:]] if "samples" in parameters else config["basic"]["samples"].split(",")
mod_samples = parameters["mod_samples"[:]] if "mod_samples" in parameters else config["basic"]["mod_samples"].split(",")
chroms = parameters["chroms"[:]] if "chroms" in parameters else config["basic"]["chromosomes"].split(",")
templ = parameters["templ"] if "templ" in parameters else config["files"]["templ"]
object_dir = parameters["object_dir"] if "object_dir" in parameters else config["paths"]["object_dir"]
stages = parameters["stages"[:]] if "stages" in parameters else config["basic"]["stages"].split(",")
gc_object = parameters["gc_file"] if "gc_file" in parameters else config["files"]["gc_object"]
min_fin = parameters["min_finite"] if "min_finite" in parameters else list(map(int, config["basic"]["min_fin"].split(",")))
win_size = parameters["win_size"] if "win_size" in parameters else config["basic"]["win_size"]
lcf = parameters["lcf"] if "lcf" in parameters else config["basic"]["lcf"]
lcf = lcf if lcf == "meth" else float(lcf)  # if lcf isn't "meth" convert to float
group_names = parameters["groups"[:]] if "groups" in parameters else config["basic"]["group_names"].split(",")
min_CpGs = parameters["min_cpgs"] if "min_cpgs" in parameters else config["basic"].getint("min_CpGs")
min_Qt = parameters["min_qt"] if "min_qt" in parameters else config["basic"].getint("min_Qt")
min_bases = parameters["min_bases"] if "min_bases" in parameters else config["basic"].getint("min_bases")
max_adj_dist = parameters["max_adj_dist"] if "max_adj_dist" in parameters else config["basic"].getint("max_adj_dist")
delta = parameters["delta"] if "delta" in parameters else float(config["basic"]["delta"])
bismark_infile = parameters["bismark"] if "bismark" in parameters else rconfig["files"]["bismark_infile"]
modern = parameters["modern"] if "modern" in parameters else rconfig["files"]["modern_infile"]
mod_name = parameters["mname"] if "mname" in parameters else rconfig["modern"]["mod_name"]
mod_species = parameters["mspecies"] if "mspecies" in parameters else rconfig["modern"]["mod_spec"]
mod_ref = parameters["mref"] if "mref" in parameters else rconfig["modern"]["mod_ref"]
mod_method = parameters["mmethod"] if "mmethod" in parameters else rconfig["modern"]["mod_method"]
gene_file = parameters["gene_file"] if "gene_file" in parameters else config["files"]["gene_file"]
cgi_file = parameters["cgi_file"] if "cgi_file" in parameters else config["files"]["cgi_file"]
cust_file1 = parameters["cust_file1"] if "cust_file1" in parameters else config["files"]["cust_file1"]
cust_file2 = parameters["cust_file2"] if "cust_file2" in parameters else config["files"]["cust_file2"]
prom_def = parameters["prom_def"[:]] if "prom_def" in parameters else config["basic"]["prom_def"].split(",")
dump_dir = parameters["dump_dir"] if "dump_dir" in parameters else config["paths"]["dump_dir"]
log_dir = parameters["log_dir"] if "log_dir" in parameters else config["paths"]["log_dir"]
report = False if parameters["noreport"] else config["options"].getboolean("report")
annot = False if parameters["noannot"] else config["options"].getboolean("annot")

if not os.path.exists(dump_dir):
        os.makedirs(dump_dir)
if not os.path.exists(log_dir):
        os.makedirs(log_dir)
    

time = datetime.datetime.now()
time = time.strftime("%d-%m-%Y_%H.%M.%S")
if not samples or samples[0] == "":
    print("Samples is a required parameter")
    sys.exit(1)
if not group_names or group_names[0] == "":
    print("Group names is a required parameter")
    sys.exit(1)
if "DMR" in stages or "permute" in stages or "plot" in stages:
    gc = t.load_object(gc_object)
    chr_names = chroms if chroms and chroms[0] != "" else gc.chr_names  # assumes user wants all chroms (or all but x)
    chr_names = [x for x in chr_names if "X" not in x]  # remove chrx from list
    samplist = []
    #for the next step, all input files must be pickled
    for samps in (samples, mod_samples):
        for sample in samps:
            if sample:
                infile = object_dir + sample + templ
                print(f"loading sample {sample}")
                input_obj = t.load_object(infile)
                samplist.append(input_obj)
    print(samples)  # print mod_samples too
if "DMR" in stages:
    dms = d.DMRs()
elif "fdr" in stages or "permute" in stages or "permutstat" in stages or "plotmethylation" in stages or "plot" in stages:
    DMR_obj_infile = parameters["dmr_infile"] if "dmr_infile" in parameters else config["files"]["DMR_obj_infile"]
    dms = t.load_object(DMR_obj_infile)
if "DMR" in stages or "fdr" in stages:
    mms = m.Mmsample()
    if bismark_infile:
        mms.create_mms_from_bismark_file(bismark_infile, gc_object, mod_name, mod_species, mod_ref, mod_method)
    else:
        mms.create_mms_from_text_file(modern)
if "DMR" in stages:
    import cProfile
    ref = False if parameters["noref"] else config["basic"].getboolean("ref")
    if ref:
        ref = mms
    min_finite = min_fin[:]
    logfile = log_dir + f"DMR_log_{time}.txt"
    (qt_up, qt_down) = dms.groupDMRs(win_size=win_size, lcf=lcf, samples=samplist, sample_groups=group_names, coord=gc, chroms=chr_names, min_finite=min_finite, min_CpGs=min_CpGs, delta=delta, ref=ref, max_adj_dist=max_adj_dist, min_bases=min_bases, min_Qt=min_Qt, fname=logfile)
    
    t.save_object(f"{object_dir}DMR_obj_{time}", dms) 
    if report:
        if annot:
            dms.annotate(gene_file, cgi_file, cust_bed1=cust_file1, cust_bed2=cust_file2, prom_def=prom_def)  
        fname = dump_dir + f"DMRs_{time}.txt"
        dms.dump_DMR(fname)
if "fdr" in stages:
    win_size_orig = win_size  # save win_size from original groupDMRs process
    #create Mmsample object
    
    samplist = []  # if dmr in stages, samplist already loaded
    for sample in samples:  
        infile = object_dir + sample + templ
        print(f"loading sample {sample}")
        input_obj = t.load_object(infile)
        samplist.append(input_obj)
    mod_samplist = []
    for sample in mod_samples:
        if sample:
            infile = object_dir + sample + templ
            print(f"loading sample {sample}")
            input_obj = t.load_object(infile)
            mod_samplist.append(input_obj)
    dmr_obj_list = []
    gc = t.load_object(gc_object)
    chr_names = chroms if chroms and chroms[0] != "" else gc.chr_names  # assumes user wants all chroms (or all but x)
    chr_names = [x for x in chr_names if "X" not in x]  # remove chrx from list
    sim_permutations = parameters["sim_permutations"] if "sim_permutations" in parameters else config["permute"].getint("sim_permutations")
    thresh = parameters["thresh"] if "thresh" in parameters else float(config["basic"]["thresh"])
    logfile = log_dir + f"fdr_DMR_log_{time}.txt"
    for perm in range(sim_permutations):
        sim_obj_list = {}
        dmr_obj_list.append(d.DMRs())
        for sample in samplist:
            print(f"sample {sample.name}, perm {perm}\n")
            samp_name = sample.name
            samp_copy = copy.deepcopy(sample)
            samp_copy.simulate(mms)
            method = sample.d_rate["method"]
            min_cov = int(sample.d_rate["min_coverage"])
            d_params = {}
            if method == "reference":
                d_params["ref"] = mms
                d_params["min_beta"] = sample.d_rate["min_beta"]
            else:
                d_params["global_meth"] = sample.d_rate["global_meth"]
            samp_copy.estimate_drate(method=method, min_cov=min_cov, **d_params)
            method = sample.methylation["algorithm"]
            win_size = sample.methylation["win_size"]
            lcf = sample.methylation["lcf"]
            if isinstance(lcf, list):
                lcf = lcf[0]
            m_params = {}
            if method == "linear" or method == "lin" or method == "logistic" or method == "log":
                m_params["slope"] = sample.methylation["slope"]
                if method == "linear" or method == "lin":
                    m_params["intercept"] = sample.methylation["intercept"]
            samp_copy.reconstruct_methylation(ref=mms, function=method, win_size=win_size, lcf=lcf, **m_params)
            #dump object to text file
            sim_obj_list[samp_name] = samp_copy
            del samp_copy  # will this fix memory errors?
            #sample.dump(f"simeth_{perm}")  
        for sample in mod_samplist:
            print(f"sample {sample.name}, perm {perm}\n")
            samp_name = sample.name
            samp_copy = copy.deepcopy(sample)
            samp_copy.simulate_modern(mms)
            #dump object to text file
            sim_obj_list[samp_name] = samp_copy
            del samp_copy  # will this fix memory errors?       
        ref = False if parameters["noref"] else config["basic"].getboolean("ref")
        if ref:
            ref = mms
        min_finite = min_fin[:]
        samps = list(sim_obj_list.values())
        dmr_obj_list[perm].groupDMRs(win_size=win_size_orig, lcf=lcf, samples=samps, sample_groups=group_names, coord=gc, chroms=chr_names, min_finite=min_finite, min_CpGs=min_CpGs, delta=delta, ref=ref, max_adj_dist=max_adj_dist, min_bases=min_bases, min_Qt=min_Qt, fname=logfile)
    statfile = log_dir + f"fdr_stats_{time}.txt"
    adjusted_DMR = dms.adjust_params(dmr_obj_list, thresh=thresh, fname=logfile, statfile=statfile)
    if annot:
        adjusted_DMR.annotate(gene_file, cgi_file, cust_bed1=cust_file1, cust_bed2=cust_file2, prom_def=prom_def)  
    fname = dump_dir + f"filtered_DMRs_{time}.txt"
    adjusted_DMR.dump_DMR(fname)
    print("done")
    

    
    
        
    
if "permute" in stages:
    num_permutations = parameters["permutations"] if "permutations" in parameters else config["permute"].getint("num_permutations")
    dmp = dms.permute(num_permutations, samplist, gc)
    i = 1
    for x in dmp:  
        fname = dump_dir + f"DMRs_{time}.{i}.txt"
        x.dump_DMR(fname)
        i += 1
    t.save_object(f"{object_dir}DMP_obj_{time}", dmp) 
if "permutstat" in stages:
    if "permute" not in stages:
        DMP_obj_infile = parameters["dmp_infile"] if "dmp_infile" in parameters else config["files"]["DMP_obj_infile"]
        dmp = t.load_object(DMP_obj_infile)
    pstat = dms.permutstat(dmp)
    fname = dump_dir + f"pstat_{time}.txt"
    dms.dump_pstat(pstat, fname)
    
DMR_chrom = parameters["dmr_chrom"] if "dmr_chrom" in parameters else config["plotmethylation"]["DMR_chrom"]
DMR_idx = parameters["dmr_idx"] if "dmr_idx" in parameters else config["plotmethylation"].getint("DMR_idx")
    
if "plotmethylation" in stages:
    fname = dump_dir + f"meth_plot_{DMR_chrom}_{DMR_idx}"
    dms.plotmethylation(DMR_chrom, DMR_idx, fname)  
if "plot" in stages:  # deal with custom files here?
    dms.plot(DMR_chrom, DMR_idx, gc, samplist, gene_file, cgi_file, widenby=5000)
    
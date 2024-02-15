#!/usr/bin/python3

import amsample as a
import mmsample as m
import tools as t
#from config import *
#import config as cfg
import glob
import argparse
import sys
import configparser as cp
import gcoordinates as gcoord
import numpy as np



# params from config can be specified on the command line
argParser = argparse.ArgumentParser()
argParser.add_argument("-co", "--config", help="path of config file")
argParser.add_argument("-f", "--filename", help="path of bam input file")
argParser.add_argument("-l", "--library", help="single or double stranded")
argParser.add_argument("-n", "--name", help="sample name")
argParser.add_argument("-a", "--abbrev", help="sample abbreviation")
argParser.add_argument("-le", "--lengths", nargs="+", help="chrom lengths--should correspond with the list of chromosomes")
argParser.add_argument("-s", "--species", help="sample species")
argParser.add_argument("-c", "--chroms", nargs="+", help="list of chromosomes to use")
argParser.add_argument("-t", "--trim", help="set True to trim ends during processing")
argParser.add_argument("-cf", "--chrom_file", help="set True for dir with one file per chromosome")
argParser.add_argument("-m", "--mapq", help="mapping quality for read")
argParser.add_argument("-q", "--qual", help="mapping quality for position")
argParser.add_argument("-st", "--stages", nargs="+", help="stages of process to be run (bam diagnose filter drate meth)")
argParser.add_argument("-o", "--objectdir", help="directory for saved (pickled) object files (include final /)")
argParser.add_argument("-fd", "--filedir", help="directory for multiple input files (include final /)")
argParser.add_argument("-od", "--outdir", help="output directory (include final /)")
argParser.add_argument("-ld", "--logdir", help="directory for log files (include final /)")
argParser.add_argument("-pd", "--picdir", help="directory for images (include final /)")
argParser.add_argument("-g", "--genome", help="path for the assembly file, eg hg19.fa.gz")
argParser.add_argument("-tf", "--text_in", help="path for the object saved in text format")
argParser.add_argument("-mo", "--modern", help="text file for modern genome")
argParser.add_argument("-b", "--bismark", help=".cov or .bedGraph file for modern genome")
argParser.add_argument("-gc", "--gc_file", help="CpG file")
argParser.add_argument("-mn", "--mname", help="modern sample name")
argParser.add_argument("-ma", "--mabbrev", help="modern sample abbreviation")
argParser.add_argument("-ms", "--mspecies", help="modern sample species")
argParser.add_argument("-mr", "--mref", help="modern sample reference genome")
argParser.add_argument("-mm", "--mmethod", help="modern sample sequencing method")
argParser.add_argument("-bed", "--nobed", help="flag for bed file", action="store_true")
argParser.add_argument("-cpg", "--create_cpg", help="flag for creating cpg file", action="store_true")
argParser.add_argument("-cr", "--cpg_ref", help="reference genome assembly for CpG file")
argParser.add_argument("-no", "--no_roam", help="flag for not running the rest of RoAM", action="store_true")
argParser.add_argument("-dm", "--dmethod", help="method of deamination rate calculation (can be reference [highly recommended] or global)")
argParser.add_argument("-mc", "--min_cov", help="minimum coverage of sites for deamination rate calculation")
argParser.add_argument("-mb", "--min_beta", help="minimum beta value for reference method of deamination rate calculation")
argParser.add_argument("-gm", "--global_meth", help="global methylation value for use with global method of deamination rate calculation")
argParser.add_argument("-rm", "--rmethod", help="method of reconstruction (can be histogram, lin, or log)")
argParser.add_argument("-lcf", "--lcf", help="low coverage factor for methylation reconstruction")
argParser.add_argument("-sl", "--slope", nargs="+", help="slope for linear/logistic methods of methylation reconstruction")
argParser.add_argument("-in", "--intercept", nargs="+", help="intercept for linear method of methylation reconstruction")
argParser.add_argument("-w", "--win_size", help="window size for reconstruction of methylation--'auto' or list of 1 val or 1 for each chrom")
argParser.add_argument("-wm", "--win_method", help="window size calculation method--'prob' or 'relerror'")
argParser.add_argument("-min", "--min_meth", help="minimum methylation level to detect")
argParser.add_argument("-p", "--p0", help="param used in winsize calculation for win_method = prob")
argParser.add_argument("-k", "--k", help="reciprocal of param used in winsize calculation for win_method = relerror")
argParser.add_argument("-max", "--max_width", help="maximum window size")

args = argParser.parse_args()
keys = [x for x in vars(args).keys() if vars(args)[x] != None]
vals = [vars(args)[x] for x in keys]
parameters = dict(zip(keys, vals))
confile = parameters["config"] if "config" in parameters else "config.ini"
config = cp.ConfigParser(interpolation=cp.ExtendedInterpolation())
config.read(confile)

filedir = parameters["filedir"] if "filedir" in parameters else config["paths"]["filedir"]
file_per_chrom = parameters["chrom_file"] if "chrom_file" in parameters else config["basic"].getboolean("file_per_chrom")
bed = False if parameters["nobed"] else config["files"].getboolean("bed")
species = parameters["species"] if "species" in parameters else config["basic"]["species"]
chroms = parameters["chroms"] if "chroms" in parameters else config["basic"]["chroms"].split(",")
genome_file = parameters["genome"] if "genome" in parameters else config["files"]["genome_file"]
object_dir = parameters["objectdir"] if "objectdir" in parameters else config["paths"]["object_dir"]
lengths = parameters["lengths"] if "lengths" in parameters else config["basic"]["chr_lengths"].split(",")
lengths = [int(x) for x in lengths]
if parameters["create_cpg"]:
    outfile = object_dir + species.replace(" ","_") + "_cpg_coords.P"
    ref = parameters["cpg_ref"] if "cpg_ref" in parameters else config["basic"]["cpg_ref"]
    if not ref or ref == "":
        ref = genome_file.split("/")[-1].split(".")[0]
    gc_obj = gcoord.Gcoordinates(chr_names=chroms, species=species, name="CpG coordinates", reference=ref)
    gc_obj.create_cpg_file(genome_file, outfile, lengths)
    if parameters["no_roam"]:
        print("Finished creating CpG file. Exiting.")
        sys.exit(0)
    gc = outfile
else:
    gc = parameters["gc_file"] if "gc_file" in parameters else config["files"]["gc_object"]
               
#def roam_pipeline(filename=cfg.filename, name=cfg.name, abbrev=cfg.abbrev, library=cfg.library):
def roam_pipeline(**params):
    name=params["name"] if "name" in params else config["required"]["name"]
    abbrev=params["abbrev"] if "abbrev" in params else config["required"]["abbrev"]
    if name == "name" or not name:
        print("Name is a required parameter")
        sys.exit(1)
    if abbrev == "abbrev" or not abbrev:
        print("Abbreviation is a required parameter")
        sys.exit(1)
    #create Amsample object
    ams = a.Amsample(name=name, abbrev=abbrev)
    # eg: ams = a.Amsample(name="Ust_Ishim", abbrev="Ust")
    stages = params["stages"[:]] if "stages" in params else config["basic"]["stages"].split(",") 
    stage = stages[0]
    if stage == "bam":
        filename = params["filename"] if "filename" in params else config["required"]["filename"]
        if not file_per_chrom:
            if filename == "filename_path" or not filename:
                print("Filename is a required parameter")
                sys.exit(1)
        library = params["library"] if "library" in params else config["required"]["library"]
        if library == "strands" or not library:
            print("library is a required parameter")
            sys.exit(1)
        trim = params["trim"] if "trim" in params else config["basic"].getboolean("trim_ends")
        mapq = int(params["mapq"]) if "mapq" in params else int(config["basic"]["mapq_thresh"])
        qual = int(params["qual"]) if "qual" in params else int(config["basic"]["qual_thresh"])
        
#populate object from bam file
        ams.bam_to_am(filename=filename, library=library, chr_lengths=lengths, species=species, trim_ends=trim, chroms=chroms, filedir=filedir, file_per_chrom=file_per_chrom, mapq_thresh=mapq, qual_thresh=qual, gc_object=gc)
       # eg: ams.bam_to_am(filename="../../ust_ishim.bam", library="single", chr_lengths=ust_chr_lengths, genome_seq="../../hg19.fa.gz", species="Homo sapiens")
        stages = stages[1:]  # remove bam stage from list 
    else:
        #get object info from text file
        text_in = params["text_in"] if "text_in" in params else config["files"]["text_infile"]
        
        ams.parse_infile(text_in)                                   
        # eg: ams.parse_infile("data/python_dumps/ust_ishim_bam.txt")                                       
    mm_flag = 0
    for stage in stages:
        picdir = params["picdir"] if "picdir" in params else config["paths"]["picdir"]
        logdir = params["logdir"] if "logdir" in params else config["paths"]["logdir"]
        if stage == "diagnose":
            ams.diagnose(picdir=picdir, logdir=logdir)
        elif stage == "filter":
            ams.filter(logdir=logdir)
        elif stage == "drate" or stage == "meth":
            if not mm_flag:
                #create Mmsample object
                bismark = params["bismark"] if "bismark" in params else config["files"]["bismark_infile"]
                modern = params["modern"] if "modern" in params else config["files"]["modern_infile"]
                mod_name = params["mname"] if "mname" in params else config["modern"]["mod_name"]
                mod_abbrev = params["mabbrev"] if "mabbrev" in params else config["modern"]["mod_abbrev"]
                mod_species = params["mspecies"] if "mspecies" in params else config["modern"]["mod_spec"]
                mod_ref = params["mref"] if "mref" in params else config["modern"]["mod_ref"]
                mod_method = params["mmethod"] if "mmethod" in params else config["modern"]["mod_method"]
                drate_method = params["dmethod"] if "dmethod" in params else config["drate"]["deamination_method"]
                recon_method = params["rmethod"] if "rmethod" in params else config["meth"]["reconstruction_method"]
                if recon_method == "histogram" or drate_method == "reference":
                    mms = m.Mmsample()
                    if bismark:
                        mms.create_mms_from_bismark_file(bismark, gc, mod_name, mod_abbrev, mod_species, mod_ref, mod_method)
                    elif modern:
                        mms.create_mms_from_text_file(modern)
                    else:
                        print("No modern reference file entered. Using default for hg19. Please make sure bone5.txt is in the current directory or specify a reference file.")
                        mms.create_mms_from_text_file("bone5.txt")
                    mm_flag += 1
                elif drate_method == "global" and stage == "drate":
                    print("We highly recommend using the reference method for deamination rate calculation.")
                else:
                    mms = ""
            if stage == "drate":
                min_cov = params["min_cov"] if "min_cov" in params else config["drate"]["min_cov"]
                min_beta = params["min_beta"] if "" in params else config["drate"]["min_beta"]
                global_meth = params["global_meth"] if "global_meth" in params else config["drate"]["global_meth"]
                drate_params = {}
                drate_params["min_cov"] = int(min_cov)
                if drate_method == "reference":
                    drate_params["ref"] = mms
                    drate_params["min_beta"] = float(min_beta)
                else:
                    drate_params["global_meth"] = float(global_meth)
                drate_params["method"] = drate_method
                ams.estimate_drate(**drate_params)
            elif stage == "meth":
                lcf = params["lcf"] if "lcf" in params else float(config["meth"]["lcf"])
                slope = params["slope"] if "slope" in params else config["meth"]["slope"]
                intercept = params["intercept"] if "intercept" in params else config["meth"]["intercept"].split(",")
                win_size = params["win_size"] if "win_size" in params else config["meth"]["win_size"]
                win_method = params["win_method"] if "win_method" in params else config["meth"]["win_method"]
                min_meth = params["min_meth"] if "min_meth" in params else config["meth"]["min_meth"]
                p0 = params["p0"] if "p0" in params else config["meth"]["p0"]
                k_recip = np.reciprocal(float(params["k"])) if "k" in params else np.reciprocal(float(config["meth"]["k"]))
                max_width = params["max_width"] if "max_width" in params else config["meth"]["max_width"]
                win_params = {}
                if win_method: win_params["method"] = win_method
                if min_meth: win_params["min_meth"] = float(min_meth)
                if p0: win_params["p0"] = float(p0)
                if k_recip: win_params["k_recip"] = float(k_recip)
                if max_width: win_params["max_width"] = int(max_width)
                
                ams.reconstruct_methylation(ref=mms, function=recon_method, win_size=win_size, lcf=lcf, slope=slope, intercept=intercept, winsize_alg=win_params)
                
    
    #dump object to text file
    outdir = params["outdir"] if "outdir" in params else config["paths"]["outdir"]
    ams.dump(stage, dir=outdir, bed=bed, gc_object=gc)
    
if filedir and not file_per_chrom:
    filenames = glob.glob(filedir+"/*.bam")
    for fn in filenames:
        name = abbrev = fn.split("/")[-1].split(".")[0]
        #roam_pipeline(filename=fn, name=name, abbrev=abbrev, library=library)
        parameters["filename"] = fn
        parameters["name"] = name
        parameters["abbrev"] = abbrev
        roam_pipeline(**parameters)
else:
    #roam_pipeline(filename=fn, name=name, abbrev=abbrev, library=library)
    roam_pipeline(**parameters)

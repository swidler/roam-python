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
argParser.add_argument("-st", "--stages", nargs="+", help="stages of process to be run")
argParser.add_argument("-o", "--objectdir", help="directory for saved (pickled) object files")
argParser.add_argument("-fd", "--filedir", help="directory for multiple input files")
argParser.add_argument("-od", "--outdir", help="output directory")
argParser.add_argument("-ld", "--logdir", help="directory for log files")
argParser.add_argument("-pd", "--picdir", help="directory for images")
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
    gc = ""
    # eg: ams = a.Amsample(name="Ust_Ishim", abbrev="Ust")
    stages = params["stages"[:]] if "stages" in params else config["basic"]["stages"].split(",") 
    stage = stages[0]
    if stage == "bam":
        filename = params["filename"] if "filename" in params else config["required"]["filename"]
        if filename == "filename_path" or not filename:
            print("Filename is a required parameter")
            sys.exit(1)
        library = params["library"] if "library" in params else config["required"]["library"]
        if library == "strands" or not library:
            print("library is a required parameter")
            sys.exit(1)
        lengths = params["lengths"] if "lengths" in params else config["basic"]["chr_lengths"].split(",")
        lengths = [int(x) for x in lengths]
        species = params["species"] if "species" in params else config["basic"]["species"]
        trim = params["trim"] if "trim" in params else config["basic"].getboolean("trim_ends")
        chroms = params["chroms"] if "chroms" in params else config["basic"]["chroms"].split(",")
        mapq = int(params["mapq"]) if "mapq" in params else int(config["basic"]["mapq_thresh"])
        qual = int(params["qual"]) if "qual" in params else int(config["basic"]["qual_thresh"])
        gc = params["gc_file"] if "gc_file" in params else config["files"]["gc_object"]
        
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
                mms = m.Mmsample()
                bismark = params["bismark"] if "bismark" in params else config["files"]["bismark_infile"]
                gc = params["gc_file"] if "gc_file" in params else config["files"]["gc_object"]
                modern = params["modern"] if "modern" in params else config["files"]["modern_infile"]
                mod_name = params["mname"] if "mname" in params else config["modern"]["mod_name"]
                mod_abbrev = params["mabbrev"] if "mabbrev" in params else config["modern"]["mod_abbrev"]
                mod_species = params["mspecies"] if "mspecies" in params else config["modern"]["mod_spec"]
                mod_ref = params["mref"] if "mref" in params else config["modern"]["mod_ref"]
                mod_method = params["mmethod"] if "mmethod" in params else config["modern"]["mod_method"]
                
                if bismark:
                    mms.create_mms_from_bismark_file(bismark, gc, mod_name, mod_abbrev, mod_species, mod_ref, mod_method)
                else:
                    mms.create_mms_from_text_file(modern)
                mm_flag += 1
            if stage == "drate":
                ams.estimate_drate(ref=mms)
            elif stage == "meth":
                ams.reconstruct_methylation(ref=mms)
                gc = params["gc_file"] if "gc_file" in params else config["files"]["gc_object"]
    
    
    #dump object to text file
    outdir = params["outdir"] if "outdir" in params else config["paths"]["outdir"]
    ams.dump(stage, dir=outdir, bed=bed, gc_object=gc)
    
if filedir and not file_per_chrom:
    filenames = glob.glob(filedir+"/*.bam")
    for fn in filenames:
        name = abbrev = fn.split("/")[-1].split(".")[0]
        #roam_pipeline(filename=fn, name=name, abbrev=abbrev, library=library)
        roam_pipeline(**parameterss)
else:
    #roam_pipeline(filename=fn, name=name, abbrev=abbrev, library=library)
    roam_pipeline(**parameters)

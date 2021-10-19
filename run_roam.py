#!/usr/bin/python3

import amsample as a
import mmsample as m
import tools as t
#from config import *
import config as cfg
import glob

def roam_pipeline(filename=cfg.filename, name=cfg.name, abbrev=cfg.abbrev):
    #create Amsample object
    ams = a.Amsample(name=name, abbrev=abbrev)
    # eg: ams = a.Amsample(name="Ust_Ishim", abbrev="Ust")
       
    stage = cfg.stages[0]
    if stage == "bam":
       #populate object from bam file
       ams.bam_to_am(filename=filename, library=cfg.library, chr_lengths=cfg.chr_lengths, genome_seq=cfg.genome_file, species=cfg.species, trim_ends=cfg.trim_ends, chroms=cfg.chroms, filedir=cfg.filedir, file_per_chrom=cfg.file_per_chrom)
       # eg: ams.bam_to_am(filename="../../ust_ishim.bam", library="single", chr_lengths=ust_chr_lengths, genome_seq="../../hg19.fa.gz", species="Homo sapiens")
       stages = cfg.stages[1:]  # remove bam stage from list
    else:
        #get object info from text file
        ams.parse_infile(cfg.text_infile)                                   
        # eg: ams.parse_infile("data/python_dumps/ust_ishim_bam.txt")                                       
    mm_flag = 0
    for stage in cfg.stages:
        if stage == "diagnose":
            ams.diagnose()
        elif stage == "filter":
            ams.filter()
        elif stage == "drate" or stage == "meth":
            if not mm_flag:
                #create Mmsample object
                mms = m.Mmsample()
                mms.create_mms_from_text_file()
                mm_flag += 1
            if stage == "drate":
                ams.estimate_drate(ref=mms)
            elif stage == "meth":
                ams.reconstruct_methylation(ref=mms)
    
    
    #dump object to text file
    ams.dump(stage)
    
if cfg.filedir and not cfg.file_per_chrom:
    filenames = glob.glob(cfg.filedir+"/*.bam")
    for fn in filenames:
        name = abbrev = fn.split("/")[-1].split(".")[0]
        roam_pipeline(filename=fn, name=name, abbrev=abbrev)
else:
    roam_pipeline(filename=cfg.filename)

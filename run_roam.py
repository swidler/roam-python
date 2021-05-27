#!/usr/bin/python3

import amsample as a
import mmsample as m
import tools as t
from config import *

#create Amsample object
ams = a.Amsample(name=name, abbrev=abbrev)
# eg: ams = a.Amsample(name="Ust_Ishim", abbrev="Ust")
   
stage = stages[0]
if stage == "bam":
   #populate object from bam file
   ams.bam_to_am(filename=filename, library=library, chr_lengths=chr_lengths, genome_seq=genome_file, species=species, trim_ends=trim_ends, chroms=chroms, filedir=filedir, file_per_chrom=file_per_chrom)
   # eg: ams.bam_to_am(filename="../../ust_ishim.bam", library="single", chr_lengths=ust_chr_lengths, genome_seq="../../hg19.fa.gz", species="Homo sapiens")
   stages = stages[1:]  # remove bam stage from list
else:
    #get object info from text file
    ams.parse_infile(text_infile)                                   
    # eg: ams.parse_infile("data/python_dumps/ust_ishim_bam.txt")                                       
mm_flag = 0
for stage in stages:
    if stage == "diagnose":
        ams.diagnose()
    elif stage == "filter":
        ams.filter()
    elif stage == "drate" or stage == "meth":
        if not mm_flag:
            #create Mmsample object
            mms = m.Mmsample()
            mms.create_mms_from_file()
            mm_flag += 1
        if stage == "drate":
            ams.estimate_drate(ref=mms)
        elif stage == "meth":
            ams.reconstruct_methylation(ref=mms)


#dump object to text file
ams.dump(stage)

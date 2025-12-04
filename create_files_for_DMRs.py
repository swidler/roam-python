#!/usr/bin/python3

import amsample as a
import mmsample as m
import tools as t
import datetime
import argparse
import os

argParser = argparse.ArgumentParser()
#argParser.add_argument("-s", "--samples", nargs="+", help="sample names")
argParser.add_argument("-ms", "--mod_samples", nargs="+", help="modern sample names")
argParser.add_argument("-o", "--object_dir", help="directory for saved (pickled) object files (include final /)")
argParser.add_argument("-d", "--data_dir", help="directory for data files from RoAM process (include final /)")
argParser.add_argument("-gc", "--gc_file", help="CpG file")
#argParser.add_argument("-t", "--templ", help="template to match any extra text in sample filename")
argParser.add_argument("-mt", "--mtempl", help="template to match any extra text in modern sample filename")
argParser.add_argument("-msp", "--mspecies", help="modern sample species")
argParser.add_argument("-mr", "--mref", help="modern sample reference genome")
argParser.add_argument("-mm", "--mmethod", help="modern sample sequencing method")
argParser.add_argument("-sm", "--smoothed", help="is sample already smmothed (use methylation and not counts). Default False")

args = argParser.parse_args()
keys = [x for x in vars(args).keys() if vars(args)[x] != None]
vals = [vars(args)[x] for x in keys]
parameters = dict(zip(keys, vals))

#samples = parameters["samples"[:]] if "samples" in parameters else []
mod_samples = parameters["mod_samples"[:]] if "mod_samples" in parameters else []
#templ = parameters["templ"] if "templ" in parameters else ""
mtempl = parameters["mtempl"] if "mtempl" in parameters else ""
object_dir = parameters["object_dir"] if "object_dir" in parameters else ""
gc_object = parameters["gc_file"] if "gc_file" in parameters else ""
mod_species = parameters["mspecies"] if "mspecies" in parameters else ""
mod_ref = parameters["mref"] if "mref" in parameters else ""
mod_method = parameters["mmethod"] if "mmethod" in parameters else ""
alsm = parameters["smoothed"] if "smoothed" in parameters else False


"""if samples:
    data_dir = parameters["data_dir"] if "data_dir" in parameters else ""
    for sample in samples:
        ams = a.Amsample()
        filename = data_dir + sample + templ + ".txt"
        ams.parse_infile(filename)
        outfile = object_dir + sample + templ
        t.save_object(outfile, ams)"""
#if mod_samples:
data_dir = parameters["data_dir"] if "data_dir" in parameters else ""
i = 0
for sample in mod_samples:
    print(f"Processing {sample}")
    if alsm:
        print("Using given methylation values (already smoothed)") 
    else:
        print("Using counts to calculate methylation")
    mms = m.Mmsample()
    bisfile = data_dir + sample + ".cov"
    if os.path.isfile(bisfile):
        mms.bismark_to_mm(bisfile, gc_object, sample, mod_species, mod_ref, mod_method, alsm)
        mms.scale()
        outfile = object_dir + sample + mtempl
        t.save_object(outfile, mms)
    else:
        raise Exception(f"No file in {data_dir} matches {sample}")
    i+=1

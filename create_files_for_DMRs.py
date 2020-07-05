#!/usr/bin/python3

import tools as t
import amsample as a

samples = ["1116","1496","2520","4873","4875","4877","4878","4914","5077","5233","5235","5236"] 
indir = "data/matlab_dumps/filtered/"
outdir = "objects/"
templ = "filt_"
for sample in samples:
    ams = a.Amsample()
    filename = indir + sample + ".txt"
    ams.parse_infile(filename)
    outfile = outdir + templ + sample
    t.save_object(outfile, ams)


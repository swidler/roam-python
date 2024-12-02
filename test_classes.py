#!/usr/bin/python3

import mmsample as m
import amsample as a
import gcoordinates as g
import tools as t
import gintervals as gint
import cProfile

#mms = m.Mmsample()
#mms.parse_infile("/mnt/x/ref_bone5_take2.txt")
#outfile = "objects/ref_bone5_take2"
#t.save_object(outfile, mms)
#gc = g.Gcoordinates(chr_names=["chr2", "chr3", "chr4", "chr8"], coords=[[234,453, 563], [222, 223, 224], [121, 122, 123, 124, 125, 126, 455, 456, 457, 458], [123, 234, 235]])
#gc = g.Gcoordinates()
#gc.parse_infile("/mnt/x/cpg_coords.txt")
#infile = "objects/bone_5"
#mms = t.load_object(infile)
ams = a.Amsample()
#ams.parse_infile("../../I1116.txt")
ams.parse_infile("data/matlab_dumps/I1116_meth.txt")
outfile = "objects/I1116"
#t.save_object(outfile, ams)
#infile = "objects/I1116"
#ams = t.load_object(infile)
#infile = "objects/cpg_coords.P"
#gc = t.load_object(infile)
#region = "chr5:89887654-987645387"
#meth = mms.region_methylation(region, gc)
#meth = ams.region_methylation(region, gc)
#print(f"meth is {meth}")
#chroms = ["chr1", "chr2", "chr3"]
#gin = gint.Gintervals(chr_names=chroms)
mms = m.Mmsample()
#mms.create_mms_from_text_file()
#mms.create_mms_from_bismark_file()
#mms.dump("bed_test")

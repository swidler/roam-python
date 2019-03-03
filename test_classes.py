#!/usr/bin/python3

import mmsample as m
import gcoordinates as g

mms = m.Mmsample()
mms.parse_infile("/mnt/x/bone_5_short.txt")
#gc = g.Gcoordinates(chr_names=["chr2", "chr3", "chr4", "chr8"], coords=[[234,453, 563], [222, 223, 224], [121, 122, 123, 124, 125, 126, 455, 456, 457, 458], [123, 234, 235]])
gc = g.Gcoordinates()
gc.parse_infile("/mnt/x/cpg_coords.txt")
region = "chr4:10100-10330"
meth = mms.region_methylation(region, gc)
print(f"meth is {meth}")

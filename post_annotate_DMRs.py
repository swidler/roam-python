#!/usr/bin/python3

import amsample as a
import mmsample as m
import tools as t
import datetime
import argparse
import os
import numpy as np
import itertools
import cDMRs as c
import pybedtools as pbt
import gintervals as gint
import gcoordinates as gcoord
import copy
import math
import random
import re


argParser = argparse.ArgumentParser()
argParser.add_argument("-di", "--dmr_infile", help="DMR object (pre-filtered)")
argParser.add_argument("-ge", "--gene_file", help="sorted text file with genes")
argParser.add_argument("-cg", "--cgi_file", help="CGI file")
argParser.add_argument("-pd", "--prom_def", type=str, help="promoter definition around TSS: a list of 2 values [before, after], where before is the number of nucleotides into the intergenic region, and after is the number of nucleotides")
#argParser.add_argument("-pd", "--prom_def", nargs="+", help="promoter definition around TSS: a list of 2 values [before, after], where before is the number of nucleotides into the intergenic region, and after is the number of nucleotides")
argParser.add_argument("-tqt", "--thresh_Qt", type=int, help="thresh_Qt from FDR run")
argParser.add_argument("-tcg", "--thresh_CpG", type=int, help="thresh_CpG from FDR run")
argParser.add_argument("-ot", "--outfile", help="Output file for FDR-filtered and annotated DMRs table file")

args = argParser.parse_args()
keys = [x for x in vars(args).keys() if vars(args)[x] != None]
vals = [vars(args)[x] for x in keys]
parameters = dict(zip(keys, vals))

DMR_obj_infile = parameters["dmr_infile"] if "dmr_infile" in parameters else config["files"]["DMR_obj_infile"]
thresh_Qt = parameters["thresh_Qt"] if "thresh_Qt" in parameters else 0
thresh_CpG = parameters["thresh_CpG"] if "thresh_CpG" in parameters else 0
cgi_file = parameters["cgi_file"] if "cgi_file" in parameters else config["files"]["cgi_file"]
gene_file = parameters["gene_file"] if "gene_file" in parameters else config["files"]["gene_file"]
fname = parameters["outfile"] if "outfile" in parameters else "annotated." + os.path.basename(DMR_obj_infile) + ".txt"
#prom_def = parameters["prom_def"[:]] if "prom_def" in parameters else config["basic"]["prom_def"].split(",")
prom_def = [int(item) for item in parameters["prom_def"].split(',')] if "prom_def" in parameters else config["basic"]["prom_def"].split(",")


def post_annotate_DMRs(self, thresh_Qt, thresh_CpG):
    adjusted_dm = copy.deepcopy(self)
    cdm = [c.cDMR() for i in range(self.no_chromosomes)]
    for chrom in range(self.no_chromosomes):
        chromosome = self.chromosomes[chrom]
        idx = sorted(list(set(list(np.where(np.array(self.cDMRs[chrom].no_CpGs) >= thresh_CpG)[0])).intersection(list(np.where(np.array(self.cDMRs[chrom].max_Qt) >= thresh_Qt)[0]))))
        if idx:
            cdm[chrom].chromosome = chromosome
            cdm[chrom].CpG_start = np.array(self.cDMRs[chrom].CpG_start)[idx]
            cdm[chrom].CpG_end = np.array(self.cDMRs[chrom].CpG_end)[idx]
            cdm[chrom].gen_start = np.array(self.cDMRs[chrom].gen_start)[idx]
            cdm[chrom].gen_end = np.array(self.cDMRs[chrom].gen_end)[idx]
            cdm[chrom].no_bases = np.array(self.cDMRs[chrom].no_bases)[idx]
            cdm[chrom].no_CpGs = np.array(self.cDMRs[chrom].no_CpGs)[idx]
            cdm[chrom].max_Qt = np.array(self.cDMRs[chrom].max_Qt)[idx]
            cdm[chrom].methylation = np.array([np.array(self.cDMRs[chrom].methylation[x])[idx] for x in range(len(self.cDMRs[chrom].methylation))])
            cdm[chrom].no_DMRs = len(idx)
            cdm[chrom].grp_methylation_statistic = np.array(self.cDMRs[chrom].grp_methylation_statistic)[idx]
    adjusted_dm.cDMRs = cdm
    return(adjusted_dm)

dmrs = t.load_object(DMR_obj_infile)

dms = post_annotate_DMRs(dmrs, thresh_Qt, thresh_CpG)

dms.annotate(gene_file, cgi_file, prom_def=prom_def) 

dms.dump_DMR(fname)


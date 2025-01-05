#!/usr/bin/python3

import tools as t
import numpy as np
from chroms import Chrom
#from config import *
import sys
import re

class Mmsample(Chrom):
    """Modern methylation sample class

    This class inherits from Chroms superclass and has attributes name, reference,
    method, no_chrs, metadata, chr_names, coord_per_position, and methylation. metadata, chr_names, and methylation are
    lists. no_chrs is determined based on the length of chr_names. An mmsample object is 
    created (with empty defaults): mms = Mmsample(). The attributes can then be populated.
    
    Methods:
    """
    
    def __init__(self, name="unknown", species="unknown", reference="", method="", metadata=[], chr_names=[], coord_per_position="", methylation=[], coverage=[]):
        self.name = name
        self.species = species
        self.reference = reference
        self.method = method
        self.metadata = metadata
        self.chr_names = chr_names
        self.coord_per_position = coord_per_position
        self.methylation = methylation
        self.coverage = coverage
        self.no_chrs = len(chr_names)   
        
    def __repr__(self): #defines print of object
        return "name: %s\nspecies: %s\nreference: %s\nmethod: %s\nmetadata: %s\nchr_names: %s\ncoord_per_position: %s\nmethylation: %s\nno_chrs: %s" % (self.name, self.species, self.reference, self.method, self.metadata, self.chr_names, self.coord_per_position, self.methylation, self.no_chrs)
    
    def parse_infile(self, infile):
        """Parses mmsample object from a text file
        
        The text file is in the format created by the Matlab dumps script associated with the mmsample class.
        """
        i = 0
        meth_flag = 0
        cov_flag = 0
        with open(infile, "rt") as mmfile:
            for line in mmfile:
                line = line.rstrip("\n") #remove line feeds
                fields = line.split(": ") 
                if i < 7: #first 7 lines are headers
                    if fields[0] == "Name":
                        self.name = fields[1]
                    elif fields[0] == "Species":
                        self.species = fields[1]
                    elif fields[0] == "Reference":
                        self.reference = fields[1]
                    elif fields[0] == "Method":
                        self.method = fields[1]
                    elif fields[0] == "Coordinates per position":
                        coord = fields[1].split(" ")
                        coord = [int(x) for x in coord]
                        self.coord_per_position = coord[0]
                    elif fields[0] == "Chromosomes":
                        self.chr_names = re.sub("[\[\]\']", "", fields[1]).split(", ")
                else:
                    if fields[0] == "Methylation:":
                        meth_flag = 1
                        continue
                    if fields[0] == "Coverage:":
                        meth_flag = 0
                        cov_flag = 1
                        continue
                    if meth_flag:
                        chrom = fields[0]
                        meth = fields[1].split(" ")
                        meth = [float(x) for x in meth] #convert all numbers to float (NaN is a float)
                        #self.chr_names.append(chrom)
                        self.methylation.append(meth)
                    elif cov_flag:
                        chrom = fields[0]
                        cov = fields[1].split(" ")
                        cov = [float(x) for x in cov] #convert all numbers to float (NaN is a float)
                        #self.chr_names.append(chrom)
                        self.coverage.append(cov)
                i += 1
        self.no_chrs = len(self.chr_names) #reassign chrom num based on new info
    
    def bismark_to_mm(self, bisfile, gc_object, mod_name, mod_spec, mod_ref, mod_method):
        """Converts Bismark result file (.cov) into mmsample object
        
        Input: empty Mmsample object, Bismark file name, gcoordinates object (reference) filename, modern sample specifics
        Output: populated Mmsample object
        """
        gc = t.load_object(gc_object)
        if bisfile.endswith(".cov"):
            file_type = "cov"
        else:
            print("File for modern sample must be .cov")
            sys.exit(1)
        self.name = mod_name
        self.species = mod_spec
        self.reference = mod_ref
        self.method = mod_method
        chr_lengths = [len(x) for x in gc.coords]
        chr_names = {}
        coord_map = []
        for chrom in range(len(gc.coords)):  # create map of coords keyed by coord, with index as value
            d = dict(enumerate(gc.coords[chrom]))
            dr = {v: k for k, v in d.items()}
            coord_map.append(dr)
        coord_per_pos = "1"
        self.coord_per_position = coord_per_pos
        meth = []
        cov = []
        for chrom in chr_lengths:
            meth.append(np.empty(chrom))  # build list of nan arrays corresponding
            meth[-1].fill(np.nan)         # to chrom lengths for methylation and for coverage
            cov.append(np.empty(chrom))
            cov[-1].fill(np.nan)
        with open(bisfile, "rt") as mmfile:
            save_line = None
            unmatched_counter = np.zeros(gc.no_chrs)
            matched_counter = np.zeros(gc.no_chrs)
            chromosome_not_in_data=0
            while True:
                if save_line:
                    line1 = save_line
                    line2 = mmfile.readline()
                else:
                    line1 = mmfile.readline()
                    line2 = mmfile.readline()
                if not line1:
                    break
                if "track" in line1:
                    line1 = line2
                    line2 = mmfile.readline()
                elif "track" in line2:
                    line2 = mmfile.readline()
                line1 = line1.rstrip("\n")
                line2 = line2.rstrip("\n")
                fields1 = line1.split("\t")
                fields2 = line2.split("\t")
                chr_name1 = fields1[0]
                chr_name2 = fields2[0]
                chr_ind = gc.index([chr_name1])[0]
                if np.isnan(chr_ind):
                    chromosome_not_in_data+=1
                    save_line = line2
                    continue
                chr_names[chr_name1] = chr_ind
                start1 = int(fields1[1])
                m1 = int(fields1[4])
                um1 = int(fields1[5])
                if len(fields2) <=1:
                    start2 = None
                    m2 = None
                    um2 = None
                else:    
                    start2 = int(fields2[1])
                    m2 = int(fields2[4])
                    um2 = int(fields2[5])
                if start1 + 1 != start2 or not line2 or chr_name1 != chr_name2:
                    if line2:
                        save_line = line2
                else:
                    save_line = None
                if start1 in coord_map[chr_ind]:
                    i = coord_map[chr_ind][start1]
                    matched_counter[chr_ind] += 1
                    if start1 + 1 == start2:
                        tot_cov = m1 + m2 + um1 + um2
                        cov[chr_ind][i] = tot_cov
                        meth[chr_ind][i] = (m1+m2)/tot_cov
                    else:
                        tot_cov = m1 + um1
                        cov[chr_ind][i] = tot_cov
                        meth[chr_ind][i] = (m1)/tot_cov
                elif (start1 - 1) in coord_map[chr_ind]:
                    i = coord_map[chr_ind][start1 - 1]
                    matched_counter[chr_ind] += 1
                    tot_cov = m1 + um1
                    cov[chr_ind][i] = tot_cov
                    meth[chr_ind][i] = (m1)/tot_cov
                else:
                    unmatched_counter[chr_ind] += 1
                if len(fields2) <=1:
                        break
                     
        print(f"found {unmatched_counter} unmatched positions and {matched_counter} matched positions") 
        print(f"proportion matched: {matched_counter/chr_lengths}")   
        self.chr_names = gc.chr_names  # list chrom names as they appear in coord file
        self.methylation = meth
        self.coverage = cov  

    def scale(self):
        """Converts all values to be between 0 and 1
        """
        for i in range(0,len(self.methylation)):
            biggest = max([x for x in self.methylation[i] if ~np.isnan(x)]) #get biggest number in list (ignore string vals)
            if biggest > 1:
                self.methylation[i] = [x/100 if ~np.isnan(x) else x for x in self.methylation[i]] #divide numeric vals by 100

    def merge(self, report=True):
        """Merges pairs of consecutive CpG positions by averaging their values
        """
        for ind in range(0,len(self.chr_names)):
            if self.coord_per_position == "1": #if coord_per_position has 1, meth vals are merged
                if report:
                    print(f"{self.chr_names[ind]} is already merged")
                continue
            if not self.coverage:
                self.coverage = [[] for x in range(self.no_chrs)] 
            self.methylation[ind] = t.nanmerge(self.methylation[ind], "average")
            self.coverage[ind] = t.nanmerge(self.coverage[ind], "sum")
        self.coord_per_position = "1"

    def region_methylation(self, region, gc, standardize=True): 
        """Computes methylation in a region as a simple average of the values in all CpGs in the region
        
        Input: genomic region, gcoordinates object of CpG positions
        Output: methylation value in the region
        """
        if standardize:
            region = t.standardize_region(region) #take region input and change to dictionary
        chrom = region["chrom"] #get chrom from region dict
        chr_ind = gc.index([chrom])[0] #find index of region chrom in gc object
        
        cpg_start = np.where(gc.coords[chr_ind] >= region["start"])[0][0] #get index of first
        cpg_end = np.where(gc.coords[chr_ind] <= region["end"])[0][-1]    #and last coords in region
        meth = self.get_methylation(chrom) #get methylation for chrom
        meth = np.nanmean(meth[1][cpg_start:cpg_end+1]) #compute average methylation
        return meth

    def smooth(self, chrom, winsize, name=None):
        """Provides smoothed methylation using a sliding window.
        
        Input: chromosome, window size.
        Output: smoothed vector, weights (as defined in algorithm.docx)
        """
        if chrom == None:  # to account for partial chrom lists, send name rather than index
            chrom = self.index([name])[0]  # get chrom index
        meth = self.get_methylation(chrom)[1] #get methylation for chrom
        MIN_VAR = 0.01**2 #1% error
        (smooth_vec, zero_for_nan) = t.nansmooth(meth, winsize, "same")
        tpl = np.ones(winsize)
        meth_sq = np.square(meth)
        m2 = t.nanconv(meth_sq, tpl, "same")
        smooth_sq = np.square(smooth_vec)
        variance = (m2 - zero_for_nan*smooth_sq) / (zero_for_nan-1) #raises divide by zero runtime warnings
        variance[variance<MIN_VAR] = MIN_VAR #raises invalid value (for nan) runtime warnings 
        weights = zero_for_nan/variance
        return (smooth_vec, weights)
    
    def create_mms_from_text_file(self, mod_infile):
        self.parse_infile(mod_infile)
        
    def create_mms_from_bismark_file(self, bis_infile, gc, mod_name, mod_species, mod_ref, mod_method):
        self.bismark_to_mm(bis_infile, gc, mod_name, mod_species, mod_ref, mod_method)
        
    def to_m(self, chroms=None):
        """Transforms methylation beta-values to M-values.
        
        Input: list of chromosome names (or indices?). If none sent, uses all chromosomes from object.
        Output: list containing M values for each chromosome
        """
        if chroms == None:
            chroms = self.chr_names
        chr_ind = self.index(chroms)
        # make sure object is merged and scaled
        if self.coord_per_position == "2":
            self = self.merge()
        self.scale()
        # loop on chromosomes
        M = []
        for chrom in chr_ind:
            vec = self.methylation[chrom]
            M.append(np.log(vec/(1-vec)))  # might cause divide by 0 warning
        return M
    
    def dump(self, desc):
        """Dumps Mmsample object to text file.
        
        Input: desc    description for use in file name.
        Output: text file in format <object_name>_<desc>.txt (directory currently hard-coded).
        """
        mname = self.name
        mname = re.sub("\s", "_", mname)
        fname = outdir + mname + "_" + desc + ".txt"
        with open(fname, "w") as fid:
            fid.write(f"Name: {mname}\nSpecies: {self.species}\nReference: {self.reference}\nMethod: {self.method}\n")
            fid.write(f"Coordinates per position: {self.coord_per_position}\n")
            fid.write(f"Chromosomes: {self.chr_names}\n")
            fid.write("Methylation:\n")
            for chrom in range(len(self.methylation)):
                meth = self.methylation[chrom]
                fid.write(f"\t{self.chr_names[chrom]}: {' '.join(map(str, meth))}\n")
            fid.write("Coverage:\n")
            for chrom in range(len(self.coverage)):
                cov = self.coverage[chrom]
                fid.write(f"\t{self.chr_names[chrom]}: {' '.join(map(str, cov))}\n")
            
        


if __name__ == "__main__":
    mms = Mmsample()
    print(mms)
    #mms2 = Mmsample(name="First Attempt", abbrev="one", coord_per_position="2")
    #print(mms2)
    mms.parse_infile("data/python_dumps/Bone_5_cov_test.txt")
    print(mms)
    chr = ["chr4", "chr1", "chr7"]
    ind = mms.index(chr)
    print(f"{chr} is at index {ind}")
    mms.scale()
    print(mms)
    chrom = "chr4"
    meth = mms.get_methylation(chrom)
    print(meth)
    chrom = 2
    meth = mms.get_methylation(chrom)
    print(meth)
    meth = mms.get_methylation()
    print(meth)
    mms.merge()
    print(mms)
    mms.merge()

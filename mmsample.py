#!/usr/bin/python3

import tools as t
import numpy as np
from chroms import Chrom
from config import *

class Mmsample(Chrom):
    """Modern methylation sample class

    This class inherits from Chroms superclass and has attributes name, abbrev, reference,
    method, no_chrs, metadata, chr_names, coord_per_position, and methylation. The last 4 are
    lists. no_chrs is determined based on the length of coord_per_position. An mmsample object is 
    created (with empty defaults): mms = Mmsample(). The attributes can then be populated.
    
    Methods:
    """
    
    def __init__(self, name="unknown", abbrev="unk", species="unknown", reference="", method="", metadata=[], chr_names=[], coord_per_position=[], methylation=[]):
        self.name = name
        self.abbrev = abbrev
        self.species = species
        self.reference = reference
        self.method = method
        self.metadata = metadata
        self.chr_names = chr_names
        self.coord_per_position = coord_per_position
        self.methylation = methylation
        self.no_chrs = len(coord_per_position) 

    def __repr__(self): #defines print of object
        return "name: %s\nabbrev: %s\nspecies: %s\nreference: %s\nmethod: %s\nmetadata: %s\nchr_names: %s\ncoord_per_position: %s\nmethylation: %s\nno_chrs: %s" % (self.name, self.abbrev, self.species, self.reference, self.method, self.metadata, self.chr_names, self.coord_per_position, self.methylation, self.no_chrs)
    
    def parse_infile(self, infile):
        """Parses mmsample object from a text file
        
        The text file is in the format created by the Matlab dumps script associated with the mmsample class.
        """
        i = 0
        with open(infile, "rt") as mmfile:
            for line in mmfile:
                line = line.rstrip("\n") #remove line feeds
                fields = line.split(": ") 
                if i < 6: #first 6 lines are headers
                    if fields[0] == "Name":
                        self.name = fields[1]
                    elif fields[0] == "Abbreviation":
                        self.abbrev = fields[1]
                    elif fields[0] == "Species":
                        self.species = fields[1]
                    elif fields[0] == "Reference":
                        self.reference = fields[1]
                    elif fields[0] == "Method":
                        self.method = fields[1]
                    elif fields[0] == "Coordinates per position":
                        self.coord_per_position = fields[1].split(" ")
                        self.coord_per_position = [int(x) for x in self.coord_per_position]
                elif i > 6: #7th line is blank, rest are chroms
                    chrom = fields[0]
                    meth = fields[1].split(" ")
                    meth = [float(x) for x in meth] #convert all numbers to float (NaN is a float)
                    self.chr_names.append(chrom)
                    self.methylation.append(meth)
                i += 1
        self.no_chrs = len(self.coord_per_position) #reassign chrom num based on new info

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
            if self.coord_per_position[ind] == 1: #if coord_per_position has 1, meth vals are merged
                if report:
                    print(f"{self.chr_names[ind]} is already merged")
                continue
            self.methylation[ind] = t.nanmerge(self.methylation[ind], "average")
            self.coord_per_position[ind] = 1

    def region_methylation(self, region, gc): 
        """Computes methylation in a region as a simple average of the values in all CpGs in the region
        
        Input: genomic region, gcoordinates object of CpG positions
        Output: methylation value in the region
        """
        region = t.standardize_region(region) #take region input and change to dictionary
        chrom = region["chrom"] #get chrom from region dict
        chr_ind = gc.indexofchr([chrom])[0] #find index of region chrom in gc object
        
        cpg_start = np.where(gc.coords[chr_ind] >= region["start"])[0][0] #get index of first
        cpg_end = np.where(gc.coords[chr_ind] <= region["end"])[0][-1]    #and last coords in region
        meth = self.get_methylation(chrom) #get methylation for chrom
        meth = np.nanmean(meth[1][cpg_start:cpg_end+1]) #compute average methylation
        return meth

    def smooth(self, chrom, winsize):
        """Provides smoothed methylation using a sliding window.
        
        Input: chromosome, window size.
        Output: smothed vector, weights (as defined in algorith.docx)
        """
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
    
    def create_mms_from_file(self):
        self.parse_infile(modern_infile)


if __name__ == "__main__":
    mms = Mmsample()
    print(mms)
    #mms2 = Mmsample(name="First Attempt", abbrev="one", coord_per_position=[2,2,2,2,2])
    #print(mms2)
    mms.parse_infile("/mnt/x/bone_5_short.txt")
    print(mms)
    chr = ["chr4", "chr1", "chr7"]
    ind = mms.indexofchr(chr)
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

#!/usr/bin/python3

import tools as t
import numpy as np
from chroms import Chrom
from config import *
import sys
import re

class Mmsample(Chrom):
    """Modern methylation sample class

    This class inherits from Chroms superclass and has attributes name, abbrev, reference,
    method, no_chrs, metadata, chr_names, coord_per_position, and methylation. metadata, chr_names, and methylation are
    lists. no_chrs is determined based on the length of chr_names. An mmsample object is 
    created (with empty defaults): mms = Mmsample(). The attributes can then be populated.
    
    Methods:
    """
    
    def __init__(self, name="unknown", abbrev="unk", species="unknown", reference="", method="", metadata=[], chr_names=[], coord_per_position="", methylation=[], coverage=[]):
        self.name = name
        self.abbrev = abbrev
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
                        coord = fields[1].split(" ")
                        coord = [int(x) for x in coord]
                        self.coord_per_position = coord[0]
                elif i > 6: #7th line is blank, rest are chroms
                    chrom = fields[0]
                    meth = fields[1].split(" ")
                    meth = [float(x) for x in meth] #convert all numbers to float (NaN is a float)
                    self.chr_names.append(chrom)
                    self.methylation.append(meth)
                i += 1
        self.no_chrs = len(self.chr_names) #reassign chrom num based on new info
    
    def bismark_to_mm(self, bisfile, gc_object):
        """Converts Bismark result file (.cov or .bedGraph) into mmsample object
        
        Input: empty Mmsample object, Bismark file name, gcoordinates object (reference)
        Output: populated Mmsample object
        """
        gc = t.load_object(gc_object)
        if bisfile.endswith(".cov"):
            file_type = "cov"
        elif bisfile.endswith(".bedGraph"):
            file_type = "bed"
        else:
            print("File for modern sample must be either .bedGraph or .cov")
            sys.exit(1)
        self.name = mod_name
        self.abbrev = mod_abbrev
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
        with open(bisfile, "rt") as mmfile:
            first = mmfile.readline().rstrip()
            second = mmfile.readline().rstrip()
            third = mmfile.readline().rstrip()
        one = first
        two = second
        if "track" in one:
            one = second
            two = third
        fields1 = one.split("\t")
        fields2 = two.split("\t")
        if int(fields1[1])+1 == int(fields2[1]):  # start in first two consecutive
            coord_per_pos = "2"
        else:
            coord_per_pos = "1"
        self.coord_per_position = coord_per_pos
        meth = []
        cov = []
        if coord_per_pos == "1":
            for chrom in chr_lengths:
                meth.append(np.empty(chrom))  # build list of nan arrays corresponding
                meth[-1].fill(np.nan)         # to chrom lengths for methylation and for coverage
                if file_type == "cov":
                    cov.append(np.empty(chrom))
                    cov[-1].fill(np.nan)
        else:
            for chrom in chr_lengths:
                meth.append(np.empty(chrom*2))  # build list of nan arrays corresponding
                meth[-1].fill(np.nan)         # to twice chrom lengths for methylation and for coverage
                if file_type == "cov":
                    cov.append(np.empty(chrom*2))
                    cov[-1].fill(np.nan)
        with open(bisfile, "rt") as mmfile:
            save_line = None
            unmatched_counter = np.zeros(gc.no_chrs)
            matched_counter = np.zeros(gc.no_chrs)
            if coord_per_pos == "2":
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
                    line1 = line1.rstrip("\n") #remove line feeds
                    fields1 = line1.split("\t")
                    line2 = line2.rstrip("\n") #remove line feeds
                    fields2 = line2.split("\t")
                    chr_name1 = fields1[0]
                    chr_name2 = fields2[0]
                    chr_ind = gc.index([chr_name1])[0]
                    if np.isnan(chr_ind):
                        continue
                    chr_names[chr_name1] = chr_ind
                    start1 = int(fields1[1])
                    start2 = int(fields2[1])
                    if start1 + 1 != start2 or not line2 or chr_name1 != chr_name2:
                        single_flag = True
                        if line2:
                            save_line = line2
                    else:
                        single_flag = False
                        save_line = None
                    if file_type == "bed":  # bedGraph coords are zero-based
                        start1 += 1
                    if single_flag:
                        coord_pos = None
                        try:
                            i = coord_map[chr_ind][start1]  # get index of first start pos in gc
                            coord_pos = "first"
                        except KeyError:
                            try:
                                i = coord_map[chr_ind][start1-1]  # try prev pos--is this 2nd cpg coord?
                                coord_pos = "second"
                            except KeyError:
                                unmatched_counter[chr_ind] += 1
                        except TypeError:
                            print(f"line is {line1}")
                        if coord_pos == "first":
                            matched_counter[chr_ind] += 1
                            meth[chr_ind][i*2] = fields1[3]
                            if file_type == "cov":
                                tot_cov1 = int(fields1[4]) + int(fields1[5])
                                cov[chr_ind][i*2] = tot_cov1
                        elif coord_pos == "second":
                            matched_counter[chr_ind] += 1
                            meth[chr_ind][i*2+1] = fields1[3]
                            if file_type == "cov":
                                tot_cov1 = int(fields1[4]) + int(fields1[5])
                                cov[chr_ind][i*2+1] = tot_cov1
                    else:
                        success = False
                        try:
                            i = coord_map[chr_ind][start1]  # get index of first start pos in gc
                            success = True
                        except:
                            unmatched_counter[chr_ind] += 1
                        if success:
                            matched_counter[chr_ind] += 1
                            meth[chr_ind][i*2] = fields1[3]
                            meth[chr_ind][i*2+1] = fields2[3]
                            if file_type == "cov":
                                tot_cov1 = int(fields1[4]) + int(fields1[5])
                                tot_cov2 = int(fields2[4]) + int(fields2[5])
                                cov[chr_ind][i*2] = tot_cov1
                                cov[chr_ind][i*2+1] = tot_cov2
            else:  # not yet tested
                for line in mmfile:
                    if "track" in line:
                        continue
                    line = line.rstrip("\n") #remove line feeds
                    fields = line.split("\t")
                    chr_name = fields[0]
                    chr_ind = gc.index([chr_name])[0]
                    chr_names[chr_name] = chr_ind
                    start = int(fields[1])
                    if file_type == "bed":  # bedGraph coords are zero-based
                        start += 1
                    try:
                        i = coord_map[chr_ind][start]  # get index of first start pos in gc
                    except:
                        unmatched_counter[chr_ind] += 1
                    matched_counter[chr_ind] += 1
                    meth[chr_ind][i] = fields[3]
                    if file_type == "cov":
                        tot_cov = int(fields[4]) + int(fields[5])
                        cov[chr_ind][i] = tot_cov
                     
        print(f"found {unmatched_counter} unmatched positions and {matched_counter} matched positions") 
        print(f"proportion matched: {matched_counter/chr_lengths}")   
        self.chr_names = list(chr_names.keys())            
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

    def region_methylation(self, region, gc): 
        """Computes methylation in a region as a simple average of the values in all CpGs in the region
        
        Input: genomic region, gcoordinates object of CpG positions
        Output: methylation value in the region
        """
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
    
    def create_mms_from_text_file(self):
        self.parse_infile(modern_infile)
        
    def create_mms_from_bismark_file(self):
        self.bismark_to_mm(bismark_infile, gc_object)
        
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
            fid.write(f"Name: {mname}\nAbbreviation: {self.abbrev}\nSpecies: {self.species}\nReference: {self.reference}\nMethod: {self.method}\n")
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
    mms.parse_infile("/mnt/x/bone_5_short.txt")
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

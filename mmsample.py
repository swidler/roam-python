#!/usr/bin/python3

import tools as t
import numpy

class Mmsample:
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
        return "name: %s\nabbrev: %s\nspecies: %s\nreferece: %s\nmethod: %s\nmetadata: %s\nchr_names: %s\ncoord_per_postion: %s\nmethylation: %s\nno_chrs: %s" % (self.name, self.abbrev, self.species, self.reference, self.method, self.metadata, self.chr_names, self.coord_per_position, self.methylation, self.no_chrs)
    
    def parse_infile(self, infile):
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
                    meth = [int(x) if x.isdigit() else numpy.nan for x in meth] #convert all numbers to int, leave NaN alone
                    self.chr_names.append(chrom)
                    self.methylation.append(meth)
                i += 1
        self.no_chrs = len(self.coord_per_position) #reassign chrom num based on new info

    def indexofchr(self, chr_name): #find index of chrom(s) by name (input and output are lists)
        result = []
        for chrom in chr_name:
            if chrom not in self.chr_names:
                result.append(numpy.nan) #returns NaN as index of nonexistent chrom
            else:
                index = self.chr_names.index(chrom)
                result.append(index)
        return result

    def scale(self):
        for i in range(0,len(self.methylation)):
            biggest = max([x for x in self.methylation[i] if isinstance(x, int)]) #get biggest number in list (ignore string vals)
            if biggest > 1:
                self.methylation[i] = [x/100 if isinstance(x, int) else x for x in self.methylation[i]] #divide numeric vals by 100

    def getmethylation(self, chrom=""):
        if not chrom: #no input, get methylation for all chroms
            result = (self.chr_names, self.methylation)
        elif isinstance(chrom, str): #chrom name entered
            index = self.chr_names.index(chrom) #get index of chrom by name
            result = (chrom, self.methylation[index])
        elif isinstance(chrom, int): #chrom index entered
            result = (self.chr_names[chrom], self.methylation[chrom])
        return result

    def merge(self, report=True):
        for ind in range(0,len(self.chr_names)):
            if self.coord_per_position[ind] == 1: #if coord_per_position has 1, meth vals are merged
                if report:
                    print(f"{self.chr_names[ind]} is already merged")
                continue
            self.methylation[ind] = t.nanconv(self.methylation[ind], "average")
            self.coord_per_position[ind] = 1

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
    meth = mms.getmethylation(chrom)
    print(meth)
    chrom = 2
    meth = mms.getmethylation(chrom)
    print(meth)
    meth = mms.getmethylation()
    print(meth)
    mms.merge()
    print(mms)
    mms.merge()

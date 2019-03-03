#!/usr/bin/python3

from chroms import Chrom
import numpy as np
class Gcoordinates(Chrom):
    def __init__(self, name="", description="", species="unknown", reference="", chr_names=[], coords=[], strand=[], metadata=[]):
        self.name = name
        self.description = description
        self.species = species
        self.reference = reference
        self.chr_names = chr_names
        self.coords = np.array(coords)
        for i in range(0,len(coords)):
            self.coords[i] = np.array(self.coords[i])
        self.strand = strand
        self.metadata = metadata
        self.no_chrs = len(coords)

    def __repr__(self): #defines print of object
        return "name: %s\ndescription: %s\nspecies: %s\nreference: %s\nmetadata: %s\nchr_names: %s\ncoords: %s\nstrand: %s\nno_chrs: %s" % (self.name, self.description, self.species, self.reference, self.metadata, self.chr_names, self.coords, self.strand, self.no_chrs)

    def parse_infile(self, infile):
        i = 0
        temp = []
        with open(infile, "rt") as gcfile:
            for line in gcfile:
                line = line.rstrip("\n") #remove line feeds
                if i == 0:
                    header = line #first line is header
                    self.description = header
                    fields = header.split(" ")
                    ref = fields[2]
                    ref = ref.rstrip(")")
                    ref = ref.lstrip("(")
                    self.reference = ref
                    self.species = fields[0]
                    self.name = fields[1]
                elif i > 1: #2nd line is blank, rest are chroms
                    fields = line.split("\t")
                    self.chr_names.append(fields[0])
                    coords = fields[1].split(",")
                    coords = [int(x) if x.isdigit else np.nan for x in coords] #convert numbers to int, leave nans
                    temp.append(np.asarray(coords))
                i += 1
        self.coords = np.asarray(temp)
        self.no_chrs = len(self.coords) #reassign chrom num based on new info



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

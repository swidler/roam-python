#!/usr/bin/python3

import tools as t
import numpy as np
from chroms import Chrom

class Amsample(Chrom):
    def __init__(self, name="unknown", abbrev="unk", species="unknown", reference="", library="", chr_names=[], coord_per_position=[], no_a = [], no_c = [], no_g = [], no_t = [], g_to_a = [], c_to_t = [], diagnostics = [], p_filters = [], is_filtered = False, is_simulated = False, methylation=[], d_rate = [], metadata=[]):
        self.name = name
        self.abbrev = abbrev
        self.species = species
        self.reference = reference
        self.library = library
        self.chr_names = chr_names
        self.coord_per_position = coord_per_position
        self.no_a = no_a
        self.no_c = no_c
        self.no_g = no_g
        self.no_t = no_t
        self.g_to_a = g_to_a
        self.c_to_t = c_to_t
        self.diagnostics = diagnostics
        self.p_filters = p_filters
        self.is_filtered = is_filtered
        self.is_simulated = is_simulated
        self.methylation = methylation
        self.d_rate = d_rate
        self.metadata = metadata
        self.no_chrs = len(coord_per_position)

    def __repr__(self): #defines print of object
        return "name: %s\nabbrev: %s\nspecies: %s\nreference: %s\nlibrary: %s\nchr_names: %s\ncoord_per_position: %s\nno_a: %s\nno_c: %s\nno_g: %s\nno_t: %s\ng_to_: %s \nc_to_t: %s\ndiagnostics: %s\np_filters: %s\nis_filtered: %s\nis_simulated: %s\nmethylation: %s\nd_rate: %s\nmetadata: %s\nno_chrs: %s" % (self.name, self.abbrev, self.species, self.reference, self.library, self.chr_names, self.coord_per_position, self.no_a, self.no_c, self.no_g, self.no_t, self.g_to_a, self.c_to_t, self.diagnostics, self.p_filters, self.is_filtered, self.is_simulated, self.methylation, self.d_rate, self.metadata, self.no_chrs)

    def get_base_no(self, chrom, base):
        chr_ind = self.indexofchr([chrom])[0]
        if base == "a":
            result = self.no_a
        elif base == "c":
            result = self.no_c
        elif base == "g":
            result = self.no_g
        elif base == "t":
            result = self.no_t
        elif base == "ga":
            result = self.g_to_a
        return result[chr_ind]

if __name__ == "__main__":
    ams = Amsample()
    print(ams)
    ams2 = Amsample(name="First", coord_per_position=[2,2,2,2,2,2,2], no_t=[1,2,3], chr_names=["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7"])
    print(ams2)
    base = ams2.get_base_no("chr2", "t")
    print(f"base: {base}")
    


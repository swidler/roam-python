#!/usr/bin/python3

from chroms import Chrom
import numpy as np

class Gcoordinates(Chrom):
    def __init__(self, name="", description="", species="unknown", reference="", chr_names=[], coords=[], strand=[], metadata=[], metadata_name=[]):
        self.name = name
        self.description = description
        self.species = species
        self.reference = reference
        self.chr_names = chr_names
        self.no_chrs = len(chr_names)
        #for chrom in range(no_chrs)):
        #    self.coords[i] = np.array(self.coords[i])
        #for chrom in range(self.no_chrs):
        strand.extend([[]]*self.no_chrs)
        coords.extend([[]]*self.no_chrs)
        self.strand = strand
        self.coords = coords
        self.metadata = metadata
        self.metadata_name = metadata_name
        

    def __repr__(self): #defines print of object
        return "name: %s\ndescription: %s\nspecies: %s\nreference: %s\nmetadata: %s\nchr_names: %s\ncoords: %s\nstrand: %s\nno_chrs: %s\nmetadata_name: %s" % (self.name, self.description, self.species, self.reference, self.metadata, self.chr_names, self.coords, self.strand, self.no_chrs, self.metadata_name)

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

    def calc_tss(self, genes):
        self.metadata_name.append("UCSC_name")
        for chrom in range(self.no_chrs):
            tss = []
            strands = []
            names = []
            metadata = []
            for ivl in genes:
                if int(ivl.chrom) == chrom+1:  # chrom index val is one less than chrom num
                    start = ivl.start
                    end = ivl.end
                    strand = ivl.strand
                    name = ivl.name
                    UCSC_name = ivl.score
                    if strand == '+':
                        pos = start
                        strand = 1
                    else:
                        pos = end
                        strand = 0
                else:
                    continue    
                tss.append(pos)
                strands.append(strand)
                names.append(name)
                metadata.append(UCSC_name)
            self.coords[chrom] = np.asarray(tss)
            self.strand[chrom] = np.asarray(strands)
            self.metadata.append(metadata)

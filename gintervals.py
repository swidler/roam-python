#!/usr/bin/python3

from chroms import Chrom
import numpy as np
import pybedtools as pbt
import tools as t

class Gintervals(Chrom):
    """Genomic interval class
    
    This class has attributes name, description, species, reference, chr_names, start, end, strand, iname, 
    metadata, metadata_name. The last 7 are lists. no_chrs is determined based on the length of chr_names.
    
    A Gintervals object is created (with empty defaults): g = Gintervals(). The attributes can then be populated.
    """
    def __init__(self, name="", description="", species="unknown", reference="", chr_names=[], start=[], end=[], strand=[], iname=[], metadata=[], metadata_name = []):
        self.name = name
        self.description = description
        self.species = species
        self.reference = reference
        self.chr_names = chr_names
        self.no_chrs = len(chr_names)
        start.extend([[]]*self.no_chrs)
        end.extend([[]]*self.no_chrs)
        iname.extend([[]]*self.no_chrs)
        strand.extend([[]]*self.no_chrs)
        self.start = start
        self.end = end
        self.strand = strand
        self.iname = iname
        self.metadata = metadata
        self.metadata_name = metadata_name

    def __repr__(self): #defines print of object
        return "name: %s\ndescription: %s\nspecies: %s\nreference: %s\nmetadata: %s\nchr_names: %s\nstart: %s\nend: %s\nstrand: %s\niname: %s\nno_chrs: %s" % (self.name, self.description, self.species, self.reference, self.metadata, self.chr_names, self.start, self.end, self.strand, self.iname, self.no_chrs)

    def calc_prom_coords(self, genes, before, after):
        """Given a Ginterval object with data on genes, this function computes promoter intervals, set as 
        <before> nt before the TSS, to <after> nt after the TSS.
        
        Input: genes    gene data from gene bed file
               before   number of nucleotides promoter is extended upstream to the TSS
               after    number of nucleotides promoter is extended downstream to the TSS (into the gene)
        Output: updated Ginterval object
        """
        self.metadata_name.append("UCSC name")
        for chrom in range(self.no_chrs):
            starts = []
            ends = []
            strands = []
            names = []
            metadata = []
            for ivl in genes:
                norm_chrom = t.normalize_chrom(ivl.chrom)
                if norm_chrom is not None and norm_chrom == chrom + 1: # chrom index val is one less than chrom num
                #if int(ivl.chrom) == chrom+1:  # chrom index val is one less than chrom num
                    start = ivl.start
                    end = ivl.end
                    strand = ivl.strand
                    name = ivl.name
                    UCSC_name = ivl.score
                    if strand == '+':
                        prom_start = start - before
                        prom_end = start + after
                        strand = 1
                    else:
                        prom_start = end - after
                        prom_end = end + before
                        strand = 0
                else:
                    next    
                starts.append(prom_start)
                ends.append(prom_end)
                strands.append(strand)
                names.append(name)
                metadata.append(UCSC_name)
            self.start[chrom] = np.asarray(starts)
            self.end[chrom] = np.asarray(ends)
            self.strand[chrom] = np.asarray(strands)
            self.iname[chrom] = np.asarray(names)
            self.metadata.append(metadata)
                
                
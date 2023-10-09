#!/usr/bin/python3

import numpy as np
class Chrom:
    """Superclass with methods common to different sample objects
    """
    def index(self, chr_name): #find index of chrom(s) by name (input and output are lists)
        # TODO: if it receives indices as input, they are returned untouched
        """Gives index of chromosome
        
        Input: chromosome name(s) as a list
        
        Output: index/indices of chromosomes (also as a list)
        """
        result = []
        for chrom in chr_name:
            if "chr" not in chrom and any("chr" in x for x in self.chr_names):  # if chr in self.chr_names, make sure input chrom has "chr"
                chrom = "chr" + chrom
            if chrom not in self.chr_names:
                result.append(np.nan) #returns NaN as index of nonexistent chrom
            else:
                index = self.chr_names.index(chrom)
                result.append(index)
        return result

    def get_methylation(self, chrom=None):
        """gets the methylation vector for a specific chromosome.
        
        Input: chromomosome index or name of chromosome. If absent, methylation
        from all chromosomes is returned.
        Output: chromosome name(s) and corresponding methylation, each as a list (returned in a tuple).
        """
        # todo: add option to take a list of specific chromosomes
        
        if isinstance(chrom, str): #chrom name entered
            index = self.chr_names.index(chrom) #get index of chrom by name
            result = (chrom, self.methylation[index])
        elif isinstance(chrom, int): #chrom index entered
            result = (self.chr_names[chrom], self.methylation[chrom])
        elif not chrom: #no input, get methylation for all chroms
            result = (self.chr_names, self.methylation)
        return result

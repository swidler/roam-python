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
            elif "chr" in chrom and all("chr" not in x for x in self.chr_names):  # if chr_names don't have "chr", remove it from input chrom
                chrom = chrom.replace('chr', '')
            if chrom not in self.chr_names:
                result.append(np.nan) #returns NaN as index of nonexistent chrom
            else:
                index = self.chr_names.index(chrom)
                result.append(index)
        return result


#!/usr/bin/python3

import numpy as np
class Chrom:
    def indexofchr(self, chr_name): #find index of chrom(s) by name (input and output are lists)
        result = []
        for chrom in chr_name:
            if chrom not in self.chr_names:
                result.append(np.nan) #returns NaN as index of nonexistent chrom
            else:
                index = self.chr_names.index(chrom)
                result.append(index)
        return result

    def get_methylation(self, chrom=""):
        if not chrom: #no input, get methylation for all chroms
            result = (self.chr_names, self.methylation)
        elif isinstance(chrom, str): #chrom name entered
            index = self.chr_names.index(chrom) #get index of chrom by name
            result = (chrom, self.methylation[index])
        elif isinstance(chrom, int): #chrom index entered
            result = (self.chr_names[chrom], self.methylation[chrom])
        return result

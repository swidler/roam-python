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

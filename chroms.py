#!/usr/bin/python3

import numpy as np
import tools as t

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

def match_hist(self, chrom, ref):
    """ Wrapper:  histogram-match the reference (ref) to this sample on a given chromosome. 
                  chrom       chrom name (e.g. "chr1")
                  ref         reference (Mmsample)
        Output:   List of methylation values at chrom for ref, after histogram-matching to sample (self) 
    """

    # get chr index
    sampidx = self.index([chrom])[0]
    refidx = ref.index([chrom])[0]
    # Collect methylation values
    if type(self).__name__ == "Amsample":      # If ancient
        meth = self.methylation["methylation"][sampidx]
    elif type(self).__name__ == "Mmsample":    # if modern
        meth = self.methylation[sampidx]
    else:
        raise ValueError(f"Cannot match histogram to {self.name}")
    # smooth reference
    win_size = 11 # hard-coded parameter
    refmeth = t.nansmooth(ref.methylation[refidx], win_size, "same")[0]
    # match ref (query) to sample (target)
    methi = t.basic_match_hist(target=meth, query=refmeth)
    return(methi)

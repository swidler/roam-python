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


    def match_hist(self, chrom, ref):
        """ Matches modern reference (bone5) to sample with histogram-matching. 
        
        Input:    self        Sample (Mmsample/Amsample)
                  chrom       chrom
	          refmeth     Beta-values of the reference at chrom
	          idx_finite  list of finite values in ref at chrom
        
        Output:   List of methylation values at chrom for ref, after histogram-matching to sample (self) 
	"""
        # Collect methylation values from sample
        if type(self).__name__ == "Amsample":      # If ancient
            meth = self.methylation["methylation"][chrom]
        elif type(self).__name__ == "Mmsample":    # if modern 
            meth = self.methylation[chrom]
        else:
            print("Cannot match histogram to ", self.name)
            sys.exit(1)
        # hard-coded parameters
        ref_bins = 1000
        sig_bins = 100
        # use only finite elements
        idx_finite = np.where(np.isfinite(ref.methylation[chrom]))
        refmeth = [ref.methylation[chrom][x] for x in idx_finite[0]] 
        # x-axis for reference signal (sample)
        ref_edges = np.linspace(0,max(refmeth),ref_bins+1)
        sig_edges = np.linspace(0,1,sig_bins+1)
        sig_binwidth = sig_edges[1] - sig_edges[0]
        # smooth the sample
        #vec = t.nansmooth(meth, win_size, "same")[0]
        vec = meth
        vec = [vec[x] for x in idx_finite[0]]
        # generate histograms
        hist = np.histogram(refmeth, bins=ref_edges)[0]
        N = np.sum(hist)
        c = np.cumsum(hist)
        N_ref = c/N
        hist = np.histogram(vec, bins=sig_edges)[0]
        N = np.sum(hist)
        c = np.cumsum(hist)
        N_sig = c/N
        # generate the mapping
        hmap = np.zeros(ref_bins)
        for i in range(ref_bins):
            # find closest value on the CDF
            imin = np.argmin(abs(N_sig - N_ref[i]))
            hmap[i] = sig_edges[imin] + 0.5*sig_binwidth
        # make the range precisely [0,1]
        hmap = (hmap - hmap[0]) / np.ptp(hmap)
        # apply the mapping
        ref1 = np.digitize(refmeth,ref_edges)
        ref1[np.where(ref1 == len(ref_edges))] = len(ref_edges) - 1
        methi = np.empty(len(ref.methylation[chrom]))
        methi[:] = np.nan
        tmp = [hmap[x-1] for x in ref1]
        d = dict(zip(idx_finite[0], tmp))
        for i in idx_finite[0]:
            methi[i] = d[i]
        return(methi)

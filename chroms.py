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


    def basic_match_hist(target, query, qbins=1000, tbins=1000):
        """
        Matches methylation in `query` to the distribution of `target` on a chromosome.

        Input:
            target  vector of methylation values (e.g., sample), assumed in [0, 1]
            query   vector of methylation values to be transformed (e.g., reference)
        Output:
            Vector of methylation values for `query` after histogram-matching to `target`.
            Same shape as `query`.
        """
        target = np.asarray(target, dtype=float)
        query  = np.asarray(query,  dtype=float)

        idx_finite = np.where(np.isfinite(query))
        if idx_finite[0].size == 0:
            # nothing to do
            return np.full_like(query, np.nan, dtype=float)
        qmeth = query[idx_finite]
        tmeth = target[idx_finite]
        qedges = np.linspace(0, np.nanmax(qmeth), qbins + 1)
        tedges = np.linspace(0, 1.0, tbins + 1)
        tbinwidth = tedges[1] - tedges[0]
        # CDF of query
        hist, _ = np.histogram(qmeth, bins=qedges)
        N = hist.sum()
        qN = np.cumsum(hist) / N
        # CDF of target
        hist, _ = np.histogram(tmeth, bins=tedges)
        N = hist.sum()
        tN = np.cumsum(hist) / N
        # build mapping: for each query-quantile, find closest target-quantile
        hmap = np.zeros(qbins)
        for i in range(qbins):
            imin = np.argmin(np.abs(tN - qN[i]))
            hmap[i] = tedges[imin] + 0.5 * tbinwidth
        # normalize hmap to [0, 1]
        if np.ptp(hmap) > 0:
            hmap = (hmap - hmap[0]) / np.ptp(hmap)
        else:
            # degenerate case: all values same
            hmap[:] = hmap[0]

        # apply mapping
        q_bins_idx = np.digitize(qmeth, qedges)
        q_bins_idx[q_bins_idx == len(qedges)] = len(qedges) - 1
        mapped = np.array([hmap[i - 1] for i in q_bins_idx])
        methi = np.full_like(query, np.nan, dtype=float)
        methi[idx_finite] = mapped
        return(methi)


    def match_hist(self, chrom, ref):
        """ Wrapper: histogram-match the reference (ref) to this sample on a given chromosome. 

        Input:    self        Sample (Mmsample/Amsample)
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
        win_size = 11
        refmeth = t.nansmooth(ref.methylation[refidx], win_size, "same")[0]

        # match ref (query) to sample (target)
        methi = basic_match_hist(target=meth, query=refmeth)
        return(methi)

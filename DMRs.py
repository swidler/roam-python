#!/usr/bin/python3

import tools as t
import numpy as np
import itertools
import datetime
import cDMRs as c
import pybedtools as pbt
import gintervals as gint
import gcoordinates as gcoord
import copy
import math
import random
import re
import matplotlib.pyplot as plt

class DMRs:
    """Differentially methylated region class
    
    This class has attributes samples (list), groups (dictionary), species, reference, chromosomes (list),
    cDMRs (list), is_ancient (list), algorithm (list), and no_samples. no_samples is determined based on the 
    length of samples.
    """
    def __init__(self, samples=[], groups={}, species="", reference="", chromosomes=[], cDMRs=[], is_ancient=[], algorithm=[]):
        self.samples = samples
        self.groups = groups
        self.species = species
        self.reference = reference
        self.chromosomes = chromosomes
        self.cDMRs = cDMRs
        self.is_ancient = is_ancient
        self.algorithm = algorithm
        self.no_chromosomes = len(self.chromosomes)
        self.no_samples = len(self.samples) 

    def __repr__(self): #defines print of object
        return "samples: %s\ngroups: %s\nspecies: %s\nreference: %s\nchromosomes: %s\ncDMRs: %s\nis_ancient: %s\nalgorithm: %s\nno_chromosomes: %s\nno_samples: %s" % (self.samples, self.groups, self.species, self.reference, self.chromosomes, self.cDMRs, self.is_ancient, self.algorithm, self.no_chromosomes, self.no_samples)

    @staticmethod
    def findDMRs(idm, iQt, icoord, trunc2clean, imeth, min_bases, min_CpGs, min_Qt):
        """Detects DMRs within the Q signals
        
        Input: idm            cDMR object for 1 chrom
               iQt            the Q-vector
               icoord         coordinates from Gcoordinates object
               trunc2clean    index that maps truncated vectors (e.g. icoord[not_nans]) to the full version (e.g. icoord)
               imeth          estimated methylation
               min_bases      minimum length of a DMR (bases)
               minCpGs        DMR must contain at least min_CpGs CpGs
               min_Qt         DMR must have Qt >= min_Qt
        Output updated idm
        """
        # binarize {Qt}
        bQt = np.array(iQt)
        bQt[iQt>0] = 1
        # find 0->1 transitions and 1->0 transitions
        dbQt = np.diff(bQt, n=1, axis=0)
        idx0to1 = [x for x, y in enumerate(dbQt) if y == 1]
        idx1to0 = [x for x, y in enumerate(dbQt) if y == -1]
        # add a last value to cover the case that bQt ends with a run of 1's
        idx1to0 = np.append(idx1to0, len(bQt))
        # collect all DMRs whose length is at least {minDMRlen}
        for pos in range(len(idx0to1)):
            # beginning and end of a putative DMR (run of Q's)
            start = idx0to1[pos]+1  # why +1?
            end = idx1to0[pos] +1  # nec for python list slicing
            
            # find the precise extent of the DMR
            maxQt = max(iQt[start:end])
            CpGs_inDMR = int(np.where(iQt[start:end]==maxQt)[0][-1]) +1  # to get number of elements, rather than index
            
            # compute the true beginning and end of the DMR in the vectors
            # (remember that iQt has an extra first value and is longer by one)
            tstart = start - 1  #are these indices right?
            tend = tstart + CpGs_inDMR - 1
            # check if putative DMR passes our filters. {dlen} is defined from
            # the first position of the first CpG to the second position of the
            # last CpG.
            dlen = icoord[tend] - icoord[tstart] + 2  #is this right?
            if maxQt >= min_Qt and dlen >= min_bases and CpGs_inDMR >= min_CpGs:
                # compute methylation
                CpG_start = trunc2clean[tstart]
                CpG_end = trunc2clean[tend]
                meth = np.nanmean(imeth[:,CpG_start:CpG_end],axis=1)
                # substitute all in {idm}
                idm.gen_start.append(icoord[tstart]) 
                idm.gen_end.append(icoord[tend]+1)
                idm.no_bases.append(dlen)
                idm.no_CpGs.append(CpGs_inDMR)
                idm.max_Qt.append(maxQt)
                idm.CpG_start.append(CpG_start)
                idm.CpG_end.append(CpG_end)
                idm.grp_methylation_statistic.append(meth)
        idm.no_DMRs = len(idm.gen_start)

    @staticmethod
    def get_region_meth(cdmr, no_samples, samples, samp_meth, coord):
        """Finds methylation for each DMR in each sample
        
        Input: cdmr        cDMR object for 1 chrom
               no_samples  number of samples
               samples     list of sample (Amsample or Mmsample) objects
               samp_meth   2 dimensional zero array of size no_samples by num DMRs in this chrom
               coord       Gcoordinates object
        Output: populated samp_meth array
        """
        for dmr in range(cdmr.no_DMRs):
                #region = f"{cdmr.chromosome}:{cdmr.gen_start[dmr]}-{cdmr.gen_end[dmr]}"
                region = {
                    "chrom":cdmr.chromosome,
                    "start":cdmr.gen_start[dmr],
                    "end":cdmr.gen_end[dmr]
                    }
                #reg_std = t.standardize_region(region)
                #import cProfile
                for samp in range(no_samples):
                    samp_meth[samp,dmr] = samples[samp].region_methylation(region, coord, standardize=False)
                    #cProfile.runctx("samp_meth[samp,dmr] = samples[samp].region_methylation(region, coord)", globals=globals(), locals=locals(), filename="data/logs/regmeth_in_profile")
        return(samp_meth)    
    
    @staticmethod
    def get_regions(cdmr):
        """Gets all regions of DMRs in format usable by pybedtools
        
        Input: cdmr    cDMR object
        Output: list of regons in proper format
        """
        regions = []
        for dmr in range(cdmr.no_DMRs):
                region = f"{cdmr.chromosome} {cdmr.gen_start[dmr]} {cdmr.gen_end[dmr]}"  # format for pybedtools
                #reg_std = t.standardize_region(region)
                regions.append(region)
        return(regions)    

    def groupDMRs(self, samples=[], sample_groups=[], coord=[], d_rate_in=[], chroms=[], winsize_alg={}, fname="DMR_log.txt", win_size="meth", lcf="meth", delta=0.5, min_bases=100, min_Qt=0, min_CpGs=10, max_adj_dist=1000, min_finite=1, max_iterations=20, tol=1e-3, report=True, match_histogram=False, ref=None, win_mod=11, mcpc=3, por=0.667):
        """Detects DMRs between two groups of samples
        
        Input: samples            list of sample (Amsample or Mmsample) objects
               sample_groups      list of sample group names
               coord              Gcoordinates object with CpG coordinates
               d_rate_in          deamination rates for ancient samples (if empty, taken from values in sample)
               chroms             list of chromsome names
               winsize_alg        a dictionary with parameters required to determine window size, see parameters 
                 for determine_shared_winsize.
               fname              log file name
               win_size           window size for smoothing. If 'meth', it is taken as the value used to reconstruct 
                   the methylation in each sample. If 'auto', a recommended value is computed for every chromosome
                   of each sample. Otherwise, it can be a scalar (used for all chromosomes in all samples), a vector 
                   over the samples (same window size is used for all chromosomes of each ancient individual, nan is
                   substituted for each modern individual), or a 2d array with values per individual and chromosome
               lcf                low coverage factor. If 'meth', it is taken as the value used in reconstructing the
                   methylation of each sample
               delta              minimum methylation difference between the two groups
               min_bases          the minimum length of each DMR in bases--shorter DMRs are filtered out
               min_Qt             DMRs with Qt < min_Qt are filtered out
               min_CpGs           DMRs whose number of CpGs is less than min_CpGs are filtered out
               max_adj_dist       max distance between adjacent CpGs within the same DMR. If the distance between
                   consecutive CpG positions is larger than max_adj_dist, the algorithm sets Qt to 0
               min_finite         an array of length no_groups stating the minimum number of ancient samples for 
                   which we require data. If in a position there are not enough samples with data, a NaN is 
                   substituted in this position. It can also be a fraction between 0 and 1, in which case it is 
                   understood as the minimum fraction of the total number of ancient samples in the group
               max_iterations     maximum number of iterations in the Newton-Raphson phase
               tol                tolerance in the Newton-Raphson phase
               report             True if reporting to the display is desired
               match_histogram    whether to perform histogram matching in pooled_methylation
               ref                mmSample reference object (used only for histogram matching).
               win_mod            window size for modern samples
               mcpc               Minimum coverage per CpG in a DMR (on average) in sample to be considered informative
                   If sample in not informative at DMR, reported average methylation will be nan. 
               por                Minimum fraction of informative samples (per DMR) in group of size n to pass 
                   filtering by mcpc. Calculated as ceil(n*por). Will be set to 1 if greater. 
        Output: modified DMR object, Qt_up, Qt_down
        """
        no_samples = len(samples)
        is_ancient = [1 if type(x).__name__ == "Amsample" else 0 for x in samples]
        if all(x==0 for x in is_ancient):
            raise Exception("All samples are modern")
        #to_remove = np.where(np.array(is_ancient) == 0)[0]
        is_meth_lcf = False
        #if len(to_remove) != 0:  # fix--should test length of array
        #    for x in to_remove:
        #        del sample_groups[x]
        #        del samples[x]
        #    no_samples = len(samples)
        #    print(f"In the current version, modern samples are ignored. {len(to_remove)} sample(s) removed.")

        chromosomes = chroms if chroms else coord.chromosomes
        if type(win_size) == str and win_size == "meth":  # futurewarning--fix this and others
            is_auto_win = False
            is_meth_win = True
        elif type(win_size) == str and win_size == "auto":  # fix
            is_auto_win = True
            is_meth_win = False
        else:
            if type(win_size) == str:
                win_size = re.split(",| ", win_size)
            if type(win_size) == list:
                win_size = [int(x) for x in win_size]
            is_auto_win = False
            is_meth_win = False
        if type(lcf) == str and lcf == "meth":  # buggy?
            is_meth_lcf = True
            lcf = np.zeros(no_samples) #should these be nan?
        elif isinstance(lcf, (int, float)):  # pre-methylation
            lcf_in = lcf
            lcf = np.zeros(no_samples)
            for samp in range(no_samples):
                lcf[samp] = lcf_in
        d_rate = np.zeros(no_samples)
        d_rate[:] = np.nan
        if d_rate_in:
            if len(d_rate_in) == no_samples:
                d_rate = d_rate_in  # these are now identical--does this matter?
            else:
                raise Exception(f"Length of d_rate ({len(d_rate_in)}) does not match number of ancient samples ({sum(is_ancient)})")
        #substitute default drates
        for samp in range(no_samples):
            if is_ancient[samp] and np.isnan(d_rate[samp]):
                d_rate[samp] = samples[samp].d_rate["rate"]["global"]
        #process low-coverage-filter
        if is_meth_lcf:
            for samp in range(no_samples):
                if is_ancient[samp]:
                    if isinstance(samples[samp].methylation["lcf"], list):
                        lcf[samp] = samples[samp].methylation["lcf"][0]
                    else:
                        lcf[samp] = samples[samp].methylation["lcf"]
        #process match_histogram
        if ref:
            match_histogram = True
        #process window size
        no_chr = len(chromosomes)
        idx_mod = np.where(np.array(is_ancient) == 0)
        if not is_auto_win:
            if is_meth_win:
                win_size = np.nan * np.ones((no_samples, no_chr))
                for samp in range(no_samples):
                    if is_ancient[samp]:
                        indices = samples[samp].index(chromosomes)
                        smp = samples[samp].methylation["win_size"]
                        win_size[samp,] = [smp[x] for x in indices]
                    else:
                        win_size[samp,] = win_mod
            else:
                #Option 1: same window size for all individuals/chromosomes
                if len(win_size) == 1:
                    if not win_size[0]%2: #win_size is even
                        win_size[0] += 1 #make it odd
                    win_size = win_size * np.ones((no_samples, no_chr))
                    if win_mod:  # prob a silly condition
                        win_size[idx_mod,] = win_mod
                    
                #same W for all chromosomes of an individual  THIS ISN'T RIGHT
                elif (type(win_size) == "list" and (np.array(win_size)).ndim == 1) or win_size.ndim ==1:
                    for samp in range(no_samples):
                        if ~np.isnan(win_size[samp]):
                            if not win_size[samp]%2: #win_size is even
                               win_size[samp] += 1 #make it odd
                    win_size = np.transpose([win_size]) * np.ones(no_chr)
                    if win_mod:  # prob a silly condition
                        win_size[idx_mod,] = win_mod
                    
                #different W for each chromosome and individual
                else:
                    for samp in range(no_samples):
                        for chrom in range(no_chr):
                            if ~np.isnan(win_size[samp, chrom]):
                                if not win_size[samp, chrom]%2: #win_size is even
                                    win_size[samp, chrom] += 1 #make it odd
                    if win_mod:  # prob a silly condition
                        win_size[idx_mod,] = win_mod
                    
        else:
            win_size = np.zeros((no_samples, no_chr))
            coverage = np.nan*np.ones(no_samples)
            drate = np.nan*np.ones(no_samples)
            winsize_alg["coverage"] = coverage
            winsize_alg["drate"] = drate
            for samp in range(no_samples):
                for chrom in range(no_chr):
                    chr_ind = samples[samp].chr_names.index(chromosomes[chrom]) #enough just to use chrom? this doesn't assume same order
                    chr_name = chromosomes[chr_ind] 
                    win_size[samp, chrom] = t.determine_shared_winsize(samples, chr_name, **winsize_alg) 
            if win_mod:  # prob a silly condition
                win_size[idx_mod,] = win_mod
        #group samples by type (eg farmer, hunter-gatherer)
        sample_names = [samples[x].name for x in range(len(samples))]
        types = dict(zip(sample_names, sample_groups))
        #get number of groups
        no_groups = len(set(x for y in types for x in types.values()))  # this must be 2 for algorithm to work
        if no_groups != 2:
                raise Exception(f"Number of groups must be 2 (currently {no_groups})")
        #get number in each group
        group_sizes = [(x, len(list(y))) for x,y in itertools.groupby(sorted(types.values()))]
        group_names = [x[0] for x in group_sizes]
        group_nums = [group_names.index(x)+1 for x in sample_groups]
        positions = []
        grp_ancient = [[],[]]
        self.groups["group_nums"] = group_nums[:]  # this list changes later, so copy here
        self.groups["group_names"] = group_names
        self.groups["no_groups"] = no_groups
        self.groups["positions"] = positions
        #process min_finite
        if np.isscalar(min_finite):
            min_finite = min_finite * np.ones(no_groups)
        for grp in range(no_groups):
            if 0 < min_finite[grp] < 1:
                min_finite[grp] = np.floor(min_finite[grp] * group_sizes[grp][1])
            positions.append([x for x,y in enumerate(samples) if types[y.name] == group_names[grp]])  # get pos of samples with this type
        for grp in range(len(positions)):
            for samp in range(len(positions[grp])):
                grp_ancient[grp].append(is_ancient[positions[grp][samp]])
        for grp in grp_ancient:
            if not all(x==grp[0] for x in grp):
                raise Exception("Groups must be composed of modern or ancient samples, but not both")
        #write to log
        if report:
            fid = open(fname, "w")
            date = datetime.datetime.now()
            date = date.strftime("%c")
            line = "groupDMRs was run on " + date
            sep = "-" * len(line)
            fid.write(f"{line}\n{sep}\n")
            fid.write(f"Comparing {no_samples} samples from {no_groups} groups:\n\n")
            #summarize input
            max_len = np.zeros(no_groups)
            for grp in range(no_groups):
                max_len[grp] = len(group_names[grp])  # this is the string length of the group name
                for pos in positions[grp]:
                    max_len[grp] = max(max_len[grp], len(samples[pos].name))  # increase max_len to longest name in group
                max_len[grp] += 10  # increase further, just to be sure
            #title
            line = "|"
            for grp in range(no_groups):
                line += group_names[grp].center(int(max_len[grp])) + "|"
            sep = "-" * len(line)
            fid.write(f"\t{line}\n\t{sep}\n")
            #line by line
            no_rows = max([x[1] for x in group_sizes])
            for row in range(no_rows):
                line = "|"
                idx = [np.nan for x in group_names]
                for grp in range(no_groups):
                    pos = group_nums.index(grp+1) if grp+1 in group_nums else ""
                    idx[grp] = pos if pos or pos == 0 else np.nan
                    if np.isnan(idx[grp]):
                        word = ""
                    else:
                        word = f"{samples[idx[grp]].name} ({100*d_rate[idx[grp]]:.2f}%)"
                        group_nums[idx[grp]] = np.nan
                    line += word.center(int(max_len[grp])) + "|"
                fid.write(f"\t{line}\n")
            fid.write(f"\t{sep}\n\n")
            #write params of job
            #delta
            fid.write(f"delta = {delta:.2f}\n\t[used to compute the lt statistics]\n")
            #min_bases
            fid.write(f"min_bases = {min_bases}\n\t[minimum length of a DMR (bases). Shorter DMRs are filtered out]\n")
            #min_Qt
            fid.write(f"min_Qt = {min_Qt:.2f}\n\t[a DMR must have Qt >= min_Qt]\n")
            #min_CpGs
            fid.write(f"min_CpGs = {min_CpGs}\n\t[a DMR must contain at least min_CpGs CpGs]\n")
            #max_adj_dist
            fid.write(f"max_adj_dist = {max_adj_dist}\n\t[max distance between adjacent CpGs within the same DMR (bases)]\n")
            #min_finite
            fid.write(f"min_finite = {min_finite}\n\t[minimum number of ancient samples per group for which we require data]\n")
            #lcf
            fid.write(f"lcf = [{lcf}]\n\t[low coverage threshold per sample]\n")
            #max_iterations
            fid.write(f"max_iterations = {max_iterations}\n\t[maximum number of iterations in Newton-Raphson]\n")
            #tol
            fid.write(f"tol = {tol}\n\t[convergence tolerance of Newton-Raphson]\n")
            #match_histogram
            mat = f"on (reference = {ref.name})" if match_histogram else "off"
            fid.write(f"histogram matching = {mat}\n\t[histogram matching of the reconstructed methylation]\n\n")
            
        #initializations
        no_chr = len(chromosomes)
        iS = list(range(no_samples))
        Qt_up = [None]*no_chr
        Qt_down = [None]*no_chr
        for chrom in range(no_chr):
            Qt_up[chrom] = [np.nan] * len(coord.coords[coord.index([chromosomes[chrom]])[0]])
            Qt_down[chrom] = Qt_up[chrom]
        samp_names = [None]*no_samples
        species = [None]*no_samples
        genome_ref = samples[0].reference
        for samp_ind in range(no_samples):
            samp_names[samp_ind] = samples[samp_ind].name
            species[samp_ind] = samples[samp_ind].species
            if samples[samp_ind].reference != genome_ref:
                raise Exception(f"Sample {samples[samp_ind].name} does not use {genome_ref}")
            if not is_ancient[samp_ind]:
                samples[samp_ind].scale()
        self.samples = samp_names
        self.chromosomes = chromosomes
        self.species = species
        self.reference = genome_ref
        #split samples between groups
        if no_groups != 2:
            raise Exception("Currently, the algorithm works only on 2 groups")
        modern_samples = [samples[x] for x in range(len(samples)) if is_ancient[x] == 0]  # are these 2 nec?
        ancient_samples = [samples[x] for x in range(len(samples)) if is_ancient[x] == 1]
        mod_idx = [x for x in range(len(samples)) if is_ancient[x] == 0]
        ancient_idx = [x for x in range(len(samples)) if is_ancient[x] == 1]
        giS = [None]*no_groups
        for grp in range(no_groups):
            giS[grp] = positions[grp]
        # compute reference winsize per chromosome
        ref_winsize = np.round(np.mean(win_size, 0))
        ref_winsize = ref_winsize.astype(int)  # ensure that winsize values are integers
        for chrom in range(no_chr):
            if not ref_winsize[chrom]%2: #win_size is even
                ref_winsize[chrom] += 1 #make it odd
        
        #loop on chroms
        cdm = [c.cDMR() for i in range(no_chr)]
        for chrom in range(no_chr):
            #report
            if report:
                print(f"Processing chromosome {chromosomes[chrom]}")
                fid.write(f"Processing chromosome {chromosomes[chrom]}\n")
            #get number of positions along the chrom
            no_pos = len(coord.coords[coord.index([chromosomes[chrom]])[0]])  # in case of diff chrom order
            #loop on groups
            meth_stat = np.zeros((no_groups, no_pos))  # methylation statistic
            meth_err = np.zeros((no_groups, no_pos))  # standard error in methylation
            for grp in range(no_groups):  
                if grp_ancient[grp][0] == 0:
                    #mij_bar = np.zeros((len(positions[grp]), no_pos))
                    #wij = np.zeros((len(positions[grp]), no_pos))
                    mij_bar = np.zeros((len(mod_idx), no_pos))
                    wij = np.zeros((len(mod_idx), no_pos))
                    for samp in range(len(mod_idx)):
                        idx_chrom = samples[mod_idx[samp]].index([chromosomes[chrom]])[0]
                        [mij_bar[samp], wij[samp]] = samples[mod_idx[samp]].smooth(idx_chrom, [int(x) for x in [win_size[mod_idx[samp], idx_chrom]]])
                        Wj = np.nansum(wij, axis=0)
                        # Calculate mm
                        mm = np.sum(wij * mij_bar, axis=0) / Wj
                        # Calculate dmm
                        dmm = np.sqrt(1 / Wj)
                        # Assign mm and dmm to m and dm respectively
                        meth_stat[grp, :] = mm
                        meth_err[grp, :] = dmm
                else:  
                    [ma, dma] = t.pooled_methylation(np.array(samples)[giS[grp]], [chromosomes[chrom]], win_size=win_size[giS[grp],chrom], lcf=lcf[giS[grp]], min_finite=min_finite[grp], max_iterations=max_iterations, tol=tol, match_histogram=match_histogram, ref=ref, ref_winsize=ref_winsize[chrom])
                    #[ma, dma] = t.pooled_methylation(np.array(samples)[ancient_idx][grp], [chromosomes[chrom]], win_size=win_size[ancient_idx[grp],chrom], lcf=lcf[ancient_idx][grp], min_finite=min_finite[grp], max_iterations=max_iterations, tol=tol, match_histogram=match_histogram, ref=ref, ref_winsize=ref_winsize[chrom])
                    meth_stat[grp,:] = ma[0]  # ma for the first (only, in this case) chrom sent
                    meth_err[grp,:] = dma[0]  # ditto
            # compute the two statistics
            meth_stat[meth_stat>1] = 1
            diffi = meth_stat[0] - meth_stat[1]  # since there must be exactly 2 groups
            idm = np.sqrt(meth_err[0]**2 + meth_err[1]**2)
            #idm = np.sqrt((meth_err[0]*np.sqrt(group_sizes[0][1]))**2 + (meth_err[1]*np.sqrt(group_sizes[1][1]))**2) # Multiply estimator error by sqrt(group_size) to get methylation error in group 
            lt_up = (diffi - delta)/idm
            lt_down = (-diffi - delta)/idm
            not_nans = np.isfinite(lt_up)
            lt_up = lt_up[not_nans]
            lt_down = lt_down[not_nans]
            coordi = coord.coords[coord.index([chromosomes[chrom]])][0]
            coordi = coordi[not_nans]
            #methi = meth[:,not_nans]
            num_finite_pos = len(coordi)
            # create index
            trunc2clean = np.array(range(no_pos))
            trunc2clean = trunc2clean[not_nans]
            # compute the distance between adjacent CpG positions in {coordi}
            # mark all positions whose preceding CpG is at most
            # p_DMRs.max_adj_dist far by 1, and all the others by 0.
            coordi_diff = [coordi[x]-coordi[x-1]-1 for x in range(1,len(coordi))]
            coordi_diff = np.floor(np.array(coordi_diff)/max_adj_dist)
            coordi_diff = [1 if x == 0 else 0 for x in coordi_diff]
            coordi_diff.insert(0,1)  # offset list by 1 to match orig matlab algorithm
            #  initialize {iQt_up}
            iQt_up = np.zeros(num_finite_pos+1)
            # compute {iQt_up} recursively
            for pos in range(num_finite_pos):
                iQt_up[pos+1] = max(0, coordi_diff[pos] * (iQt_up[pos] + lt_up[pos]))
            # initialize {iQt_down}
            iQt_down = np.zeros(num_finite_pos+1)
            # compute {iQt_down} recursively
            for pos in range(num_finite_pos):
                iQt_down[pos+1] = max(0, coordi_diff[pos] * (iQt_down[pos] + lt_down[pos]))
            # make a clean version of {Qt}, that takes into account the many CpG
            # positions we had removed earlier. The function reports this vector,
            # which is useful later for plotting
            Qt_up[chrom] = np.array(Qt_up[chrom])
            Qt_up[chrom][not_nans] = iQt_up[1:]
            Qt_down[chrom] = np.array(Qt_down[chrom])
            Qt_down[chrom][not_nans] = iQt_down[1:]
            # filter and characterize DMRs
            self.findDMRs(cdm[chrom], iQt_up, coordi, trunc2clean, meth_stat, min_bases, min_CpGs, min_Qt)
            self.findDMRs(cdm[chrom], iQt_down, coordi, trunc2clean, meth_stat, min_bases, min_CpGs, min_Qt)
            cdm[chrom].chromosome = chromosomes[chrom]
            # compute methylation in each sample
            samp_meth = np.zeros((len(samples),cdm[chrom].no_DMRs))
            #import cProfile
            #cProfile.runctx("samp_meth = self.get_region_meth(cdm[chrom], no_samples, samples, samp_meth, coord)", globals=globals(), locals=locals(), filename="data/logs/regmeth_profile")
            samp_meth = self.get_region_meth(cdm[chrom], no_samples, samples, samp_meth, coord)
            #meth = np.array(cdm[chrom].methylation)
            #meth = np.transpose(meth)
            #samp_meth = np.insert(samp_meth,0,meth,axis=0)
            cdm[chrom].methylation = samp_meth
            # Filter DMR list by within-group variance
            mpor = 0.7    # hardcoded parameter
            mean_list = []
            gdiff_list = []
            for grp in range(len(giS)):
                mean_list.append([np.nanmean(elements) for elements in zip(*cdm[chrom].methylation[giS[grp]])])    # get group mean
                gdiff_list.append([np.nanmax(elements)-np.nanmin(elements) for elements in zip(*cdm[chrom].methylation[giS[grp]])])    # get within-group diffs
            maxdiff = [abs(a - b)*mpor for a, b in zip(*mean_list)]    # determine max allowed within-group diff by between-group diff
            idx = list(np.where(np.array([np.nanmax(elements) for elements in zip(*gdiff_list)]) <= np.array(maxdiff))[0])    # get index where within-group diffs are no greater than maximum
            del(mpor, mean_list, gdiff_list, maxdiff)
            # Filter DMR list by number of non-informative samples in DMR
            if por > 1:
                por = 1    # minimum fraction of informative samples per group
            mincov = np.array([cc*mcpc for cc in cdm[chrom].no_CpGs])    # min required coverage per CpG
            for grp in range(len(giS)):
                ms = int(np.ceil(por*len(giS[grp])))    # min number of informative samples
                counti = np.zeros(len(cdm[chrom].no_CpGs))    # will count the number of informative samples in group per DMR
                for samp in giS[grp]:
                    scov = np.array([np.nansum(samples[samp].no_t[chrom][start:end+1])+np.nansum(samples[samp].no_c[chrom][start:end+1]) for start,end in zip(np.array(cdm[chrom].CpG_start), np.array(cdm[chrom].CpG_end))])    # get overall cov in each DMR
                    cidx = list(np.where(scov >= mincov))
                    counti[cidx] += 1
                    ncidx = list(np.where(scov < mincov))
                    cdm[chrom].methylation[samp][ncidx] = np.nan # Turn reported methylation in non-informative samples to nan
                gidx = np.where(counti >= ms) # keep only DMRs where at least min number of samples are informative
                idx = np.intersect1d(idx, gidx).tolist() # keep in idx just DMRs that pass this threshold and previous ones
            del(gidx, scov, cidx, ncidx, counti, ms)

            if idx:    # filter DMRs
                cdm[chrom].CpG_start = np.array(cdm[chrom].CpG_start)[idx]
                cdm[chrom].CpG_end = np.array(cdm[chrom].CpG_end)[idx]
                cdm[chrom].gen_start = np.array(cdm[chrom].gen_start)[idx]
                cdm[chrom].gen_end = np.array(cdm[chrom].gen_end)[idx]
                cdm[chrom].no_bases = np.array(cdm[chrom].no_bases)[idx]
                cdm[chrom].no_CpGs = np.array(cdm[chrom].no_CpGs)[idx]
                cdm[chrom].max_Qt = np.array(cdm[chrom].max_Qt)[idx]
                cdm[chrom].methylation = np.array([np.array(cdm[chrom].methylation[x])[idx] for x in range(len(cdm[chrom].methylation))])  
                cdm[chrom].no_DMRs = len(idx)
                cdm[chrom].grp_methylation_statistic = np.array(cdm[chrom].grp_methylation_statistic)[idx]
            if report:
                print(f"\tdetected {cdm[chrom].no_DMRs} DMRs")
                fid.write(f"\tdetected {cdm[chrom].no_DMRs} DMRs\n")
                
        #substitue fields
        alg_props = {}
        alg_props["delta"] = delta
        alg_props["min_bases"] = min_bases
        alg_props["min_Qt"] = min_Qt
        alg_props["min_CpGs"] = min_CpGs
        alg_props["max_adj_dist"] = max_adj_dist
        alg_props["drate"] = d_rate
        alg_props["win_size"] = win_size
        alg_props["min_finite"] = min_finite
        alg_props["lcf"] = lcf
        alg_props["max_iterations"] = max_iterations
        alg_props["tol"] = tol
        alg_props["match_histogram"] = match_histogram
        alg_props["algorithm"] = "groupDMRs"
        alg_props["ref"] = ref
        self.algorithm = alg_props
        self.cDMRs = cdm
        self.no_chromosomes = no_chr
        self.no_samples = no_samples
        
        #close filehandle
        if report:
            fid.close()
        
        return(Qt_up, Qt_down)
    
    def annotate(self, gene_bed, cgi_bed, prom_def=[5000,1000], cust_bed1=None, cust_bed2=None):
        """Retrieves important data about each DMR
        
        Input: gene_bed   bed file with gene data 
               cgi_bed    bed file with CGI data
               prom_def   promoter definition around TSS, a list of 2 values [before, after], where before is the
                   number of nucleotides into the intergenic region, and after is the number of nucleotides 
                   into the gene
        Output: cDMR object with updated annotation
        """
        if gene_bed:
            genes = pbt.BedTool(gene_bed)
            genes_no_dups = genes.groupby(g=[1,2,3,6], c='4,5', o='distinct').cut([0,1,2,4,5,3])  # is stream nec?
            before = int(prom_def[0])
            after = int(prom_def[1])
            proms = gint.Gintervals(chr_names=self.chromosomes)
            proms.calc_prom_coords(genes_no_dups, before, after)
            tss = gcoord.Gcoordinates(chr_names=self.chromosomes, description="TSS positions")
            tss.calc_tss(genes_no_dups)
        else:
            genes = None
        #genes_no_dups = genes.groupby(g=[1,2,3,6], c='4,5', o='distinct').cut([0,1,2,4,5,3], stream=False)
        if cgi_bed:
            cgis = pbt.BedTool(cgi_bed)
        else:
            cgis = None
        if cust_bed1:
            cust1 = pbt.BedTool(cust_bed1)
        else:
            cust1 = None
        if cust_bed2:
            cust2 = pbt.BedTool(cust_bed2)
        else:
            cust2 = None
        #cgis_no_dups = cgis.groupby(g=[1,2,3,6], c='4,5', o='distinct').cut([0,1,2,4,5,3], stream=False)
        
        # loop on chromosomes
        for chrom in range(self.no_chromosomes):
            print(f"chrom {self.chromosomes[chrom]}")
            num_DMRs = self.cDMRs[chrom].no_DMRs
            # add handling for no DMRs?
            # get DMR regions
            regions = self.get_regions(self.cDMRs[chrom])
            dmr_annot = []
            if genes:
                # generate an array of TSSs
                # add handling of no tss?
                itss = tss.coords[chrom]
                itss = itss.astype(float)
                prom_string = ""
                for prom in range(len(proms.start[chrom])):
                    prom_string += f"\n{proms.chr_names[chrom]} {proms.start[chrom][prom]} {proms.end[chrom][prom]} {proms.iname[chrom][prom]} 0 {proms.strand[chrom][prom]}"
                chrom_proms = pbt.BedTool(prom_string, from_string=True)
                #get genes in this chrom
                genes_chrom = genes_no_dups.filter(lambda a: a.chrom == str(chrom+1)).saveas()
            for dmr in range(num_DMRs):
                annot = {}
                region = regions[dmr]
                ivl = pbt.BedTool(region, from_string=True)[0]  # get 1st element to make it an interval object
                if cgis:
                    if cgis.any_hits(ivl):
                        in_CGI = True
                    else:  
                        in_CGI = False
                else:
                    in_CGI = "N/A"
                if cust1:
                    if cust1.any_hits(ivl):
                        in_cust1 = True
                    else:
                        in_cust1 = False
                else:
                    in_cust1 = "N/A"
                if cust2:
                    if cust2.any_hits(ivl):
                        in_cust2 = True
                    else:
                        in_cust2 = False
                else:
                    in_cust2 = "N/A"
                if genes:
                    in_gene = {}
                    region_gene = region.replace("chr", "")
                    ivl_gene = pbt.BedTool(region_gene, from_string=True)[0]  # get 1st element to make it an interval object
                    gene_hits = genes.all_hits(ivl_gene)
                    if gene_hits:
                        in_gene["present"] = True
                        in_gene["name"] = []
                        in_gene["strand"] = []
                        for hit in gene_hits:
                            if hit.name not in in_gene["name"]:
                                in_gene["name"].append(hit.name)
                                if hit.strand == "+":
                                    strand = 1
                                elif hit.strand == "-":
                                    strand = 0
                                else:
                                    strand = np.nan
                                in_gene["strand"].append(strand)
                    else:
                        in_gene["present"] = False
                        in_gene["strand"] = np.nan
                        in_gene["name"] = []
                    in_prom = {}
                    prom_hits = chrom_proms.all_hits(ivl)
                    if prom_hits:
                        in_prom["present"] = True
                        in_prom["name"] = []
                        in_prom["strand"] = []
                        for hit in prom_hits:
                            if hit.name not in in_prom["name"]:
                                in_prom["name"].append(hit.name)
                                in_prom["strand"].append(hit.strand)  # change to int!
                    else:
                        in_prom["present"] = False
                        in_prom["strand"] = np.nan
                        in_prom["name"] = []
                    upstream_TSS = {}
                    up_plus = self.cDMRs[chrom].gen_end[dmr] - itss
                    up_minus = itss - self.cDMRs[chrom].gen_start[dmr]
                    up_plus[np.where(up_plus < 0)] = np.nan
                    up_minus[np.where(up_minus < 0)] = np.nan
                    up_plus[np.where(tss.strand[chrom] == 0)] = np.nan
                    up_minus[np.where(tss.strand[chrom] == 1)] = np.nan
                    closest_plus = np.nanmin(up_plus)
                    idx_plus = np.nan if np.isnan(closest_plus) else np.nanargmin(up_plus)
                    closest_plus = closest_plus - self.cDMRs[chrom].no_bases[dmr] +1
                    if closest_plus < 0:
                        closest_plus = 0
                    closest_minus = np.nanmin(up_minus)
                    idx_minus = np.nan if np.isnan(closest_minus) else np.nanargmin(up_minus)
                    closest_minus = closest_minus - self.cDMRs[chrom].no_bases[dmr] +1
                    if closest_minus < 0:
                        closest_minus = 0
                    if closest_minus == closest_plus:
                        upstream_TSS["dist"] = closest_plus
                        upstream_TSS["name"] = [genes_chrom[int(idx_plus)].name, genes_chrom[int(idx_minus)].name]
                        upstream_TSS["strand"] = [1,0]
                    else:
                        absmin = np.nanmin([closest_minus, closest_plus])
                        if closest_minus == absmin:
                            upstream_TSS["dist"] = closest_minus
                            upstream_TSS["name"] = [genes_chrom[int(idx_minus)].name]
                            upstream_TSS["strand"] = [0]
                        else:
                            upstream_TSS["dist"] = closest_plus
                            upstream_TSS["name"] = [genes_chrom[int(idx_plus)].name]
                            upstream_TSS["strand"] = [1]
                    
                    downstream_TSS = {}
                    down_minus = self.cDMRs[chrom].gen_end[dmr] - itss
                    down_plus = itss - self.cDMRs[chrom].gen_start[dmr]
                    down_plus[np.where(down_plus < 0)] = np.nan
                    down_minus[np.where(down_minus < 0)] = np.nan
                    down_plus[np.where(tss.strand[chrom] == 0)] = np.nan
                    down_minus[np.where(tss.strand[chrom] == 1)] = np.nan
                    closest_plus = np.nanmin(down_plus)
                    idx_plus = np.nan if np.isnan(closest_plus) else np.nanargmin(down_plus)
                    closest_plus = closest_plus - self.cDMRs[chrom].no_bases[dmr] +1
                    if closest_plus < 0:
                        closest_plus = 0
                    closest_minus = np.nanmin(down_minus)
                    idx_minus = np.nan if np.isnan(closest_minus) else np.nanargmin(down_minus)
                    closest_minus = closest_minus - self.cDMRs[chrom].no_bases[dmr] +1
                    if closest_minus < 0:
                        closest_minus = 0
                    if closest_minus == closest_plus:
                        downstream_TSS["dist"] = closest_plus
                        downstream_TSS["name"] = [genes_chrom[int(idx_plus)].name, genes_chrom[int(idx_minus)].name]
                        downstream_TSS["strand"] = [1,0]
                    else:
                        absmin = np.nanmin([closest_minus, closest_plus])
                        if closest_minus == absmin:
                            downstream_TSS["dist"] = closest_minus
                            downstream_TSS["name"] = [genes_chrom[int(idx_minus)].name]
                            downstream_TSS["strand"] = [0]
                        else:
                            downstream_TSS["dist"] = closest_plus
                            downstream_TSS["name"] = [genes_chrom[int(idx_plus)].name]
                            downstream_TSS["strand"] = [1]
                else:
                    in_gene = "N/A"
                    in_prom = "N/A"
                    upstream_TSS = "N/A"
                    downstream_TSS = "N/A"
                    
                annot["in_CGI"] = in_CGI
                annot["in_cust1"] = in_cust1
                annot["in_cust2"] = in_cust2
                annot["in_gene"] = copy.deepcopy(in_gene)
                annot["in_prom"] = copy.deepcopy(in_prom)
                annot["upstream_TSS"] = copy.deepcopy(upstream_TSS)
                annot["downstream_TSS"] = copy.deepcopy(downstream_TSS)
                
                dmr_annot.append(copy.deepcopy(annot))
            
            self.cDMRs[chrom].annotation = (copy.deepcopy(dmr_annot))
            
    def permute(self, no_permutations, samples, coord):
        """Detects permuted DMRs between two groups of samples.
        
        Input: no_permutations    number of permutatios to perform
               samples            the same samples used to generate the original DMRs object
               coord              coordinates of the CpGs in the genome (gcoordinates object)
        Output:     dm_permuted a list of DMRs objects, holding the results of the permutations
        """
        dmp = [DMRs() for i in range(no_permutations)]  # create a list of DMR objects
        no_samples = len(samples)
        grp = self.groups
        groups_base = [grp["group_names"][x-1] for x in grp["group_nums"]]
        alg = self.algorithm
        min_finite_base = alg["min_finite"]
        
        # perform permutations
        num_width = np.ceil(math.log10(no_permutations))
        for permutation in range(no_permutations):
            line = f"permutation #{permutation+1}"
            sep = "-" * len(line)
            print(f"{line}\n{sep}")
            groups = random.sample(groups_base, len(groups_base))
            min_finite = min_finite_base[:]
            num = str(permutation+1).zfill(int(num_width))
            fname = f"p{num}_groupDMRs.txt"
            dmp[permutation].groupDMRs(samples=samples, sample_groups=groups, coord=coord, fname=fname, chroms=self.chromosomes, win_size=alg["win_size"], lcf=alg["lcf"], delta=alg["delta"], min_bases=alg["min_bases"], min_Qt=alg["min_Qt"], min_CpGs=alg["min_CpGs"], max_adj_dist=alg["max_adj_dist"], min_finite=min_finite, max_iterations=alg["max_iterations"], tol=alg["tol"], match_histogram=alg["match_histogram"], ref=alg["ref"])
        return dmp
        
    def permutstat(self, dmp):
        """Computes statistics of permuation test on detected DMRs.
        
        Input: dmp    a list of DMRs objects, holding the results of the permutations
        Output: pstat    dictionary with fields:
                            no_oDMRs: number of observed DMRs in each chromosome.
                            no_pDMRs: number of DMRs in each permutation and each chromosome.
                            chrom_fdr: FDR for each chromosome.
                            tot_fdr: total FDR over the entire genome.
                            chrom_pval: p-value for each chromosome.
                            tot_pval: total p-value over the entire genome.
        """
        # initialize
        pstat = {}

        # get parameters
        no_chrs = self.no_chromosomes
        no_permutations = len(dmp)

        # compute #permuted DMRs
        pstat["no_pDMRs"] = np.full((no_permutations,no_chrs), np.nan)
        tot_pDMRs = np.zeros(no_permutations)
        tot_pDMRs[:] = np.nan
        for permutation in range(no_permutations):
            (tot_pDMRs[permutation], pstat["no_pDMRs"][permutation,:]) = dmp[permutation].noDMRs()
        
        # compute #observed DMRs
        (tot_oDMRs, pstat["no_oDMRs"]) = self.noDMRs()
        
        # calculate FDR
        no_pDMRs_mean = np.mean(pstat["no_pDMRs"], axis=0)
        pstat["chrom_fdr"] = no_pDMRs_mean / pstat["no_oDMRs"]
        pstat["tot_fdr"] = np.mean(tot_pDMRs) / tot_oDMRs
        
        #calculate p-value
        pstat["chrom_pval"] = np.sum(pstat["no_pDMRs"] >= pstat["no_oDMRs"], axis=0) / no_permutations
        pstat["tot_pval"] = sum(tot_pDMRs >= tot_oDMRs) / no_permutations

        return pstat  
            
    def noDMRs(self):
        """ Reports the number of DMRs in each chromosome.
        
        Input: DMR object
        Output: number of DMRs per chromosome, total number of DMRs
        """
        
        no_DMRs = np.zeros(self.no_chromosomes)
        no_DMRs[:] = np.nan
        for chrom in range(self.no_chromosomes):
            no_DMRs[chrom] = self.cDMRs[chrom].no_DMRs
        tot_DMRs = np.sum(no_DMRs)
        return(tot_DMRs, no_DMRs)
    
    @staticmethod   
    def dump_pstat(pstat, fname):
        """Dumps pstat data to text file
        
        Input: pstat dictionary
        Output: text file in format pstat_<time>.txt
        """
        
        with open(fname, "w") as fid:
            fid.write(f"Num pDMRs per chrom: {pstat['no_pDMRs']}\n")
            fid.write(f"Num observed DMRs: {pstat['no_oDMRs']}\n")
            fid.write(f"FDR per chromosome: {pstat['chrom_fdr']}\n")
            fid.write(f"Total FDR: {pstat['tot_fdr']}\n")
            fid.write(f"P-value per chromosome: {pstat['chrom_pval']}\n")
            fid.write(f"Total P-value: {pstat['tot_pval']}\n")
            
            
    
    def dump_DMR(self, fname):
        """Dumps DMR object to text file.
        
        Input: DMR object, iterator
        Output: text file in format DMRs_<time>.txt (directory currently hard-coded).
        """
        fname_gen = fname.replace("DMRs_", "DMR_gen_")
        with open(fname_gen, "w") as fid:
            fid.write(f"Samples: {self.samples}\n")
            fid.write(f"Group Assignment: {self.groups['group_nums']}\nGroup Naming: {self.groups['group_names']}\n")
            fid.write(f"Species: {self.species}\n")
            fid.write(f"Ancient samples: {self.is_ancient}\n")
            fid.write(f"Reference: {self.reference}\n")
            fid.write(f"Chromosomes: {self.chromosomes}\n")
            #fid.write(f"Algorithm: {self.algorithm}\n")
            i = 0
            for key in self.algorithm:
                if key == "win_size":
                    fid.write(f"\twin_size:\n")
                    for sample in self.algorithm[key]:
                        s_name = self.samples[i]
                        fid.write(f"\t\t{s_name}: {sample}\n")
                        i += 1
                elif key == "ref" and self.algorithm[key] != False:
                    fid.write(f"\tref: {self.algorithm['ref'].name}\n")
                else:
                    fid.write(f"\t{key}: {self.algorithm[key]}\n")
        samp_names = ""
        group_names = ""
        for samp in self.samples:
            samp_names += samp + "_average_methylation\t"
        for group in self.groups['group_names']:
            group_names += group + "_meth_statistic\t"
        with open(fname, "w") as fid:
            #fid.write("cDMRs:\n")
            c1 = ""
            c2 = ""
            cg = ""
            g = "\n"
            for chrom in range(self.no_chromosomes):
                if self.cDMRs[chrom].annotation:
                    if self.cDMRs[chrom].annotation[0]['in_cust1'] != "N/A":
                        c1 = "in_cust1\t"
                    if self.cDMRs[chrom].annotation[0]['in_cust2'] != "N/A":
                        c2 = "in_cust2\t"
                    if self.cDMRs[chrom].annotation[0]['in_CGI'] != "N/A":
                        cg = "in_CGI\t"  
                    if self.cDMRs[chrom].annotation[0]['in_gene'] != "N/A":
                        g = "in_gene\tname(s)\tstrand(s)\tin_prom\tname(s)\tstrand(s)\tupstream_TSS\tname(s)\tstrand(s)\tdownstream_TSS\tname(s)\tstrand(s)\n"            
            fid.write(f"Chrom\tDMR#\tout_of\tGenomic_start\tGenomic_end\tCpG_start\tCpG_end\t#CpGs\t#bases\tMax_Qt\t{group_names}{samp_names}{c1}{c2}{cg}{g}")
            for chrom in range(self.no_chromosomes):
                for dmr in range(self.cDMRs[chrom].no_DMRs):
                    
                    fid.write(f"{self.chromosomes[chrom]}\t")
                    fid.write(f"{dmr+1}\t")
                    fid.write(f"{self.cDMRs[chrom].no_DMRs}\t")
                    fid.write(f"{self.cDMRs[chrom].gen_start[dmr]}\t")
                    fid.write(f"{self.cDMRs[chrom].gen_end[dmr]}\t")
                    fid.write(f"{self.cDMRs[chrom].CpG_start[dmr]}\t")
                    fid.write(f"{self.cDMRs[chrom].CpG_end[dmr]}\t")
                    fid.write(f"{self.cDMRs[chrom].no_CpGs[dmr]}\t")
                    fid.write(f"{self.cDMRs[chrom].no_bases[dmr]}\t")
                    fid.write(f"{self.cDMRs[chrom].max_Qt[dmr]}\t")
                    for grp in range(self.groups['no_groups']):
                        fid.write(f"{self.cDMRs[chrom].grp_methylation_statistic[dmr][grp]}\t")
                    for sample in range(self.no_samples):
                        fid.write(f"{self.cDMRs[chrom].methylation[sample][dmr]}\t")
                    if self.cDMRs[chrom].annotation:
                        if self.cDMRs[chrom].annotation[dmr]['in_CGI'] != "N/A":
                            fid.write(f"{self.cDMRs[chrom].annotation[dmr]['in_CGI']}\t")
                        if self.cDMRs[chrom].annotation[dmr]['in_cust1'] != "N/A":
                            fid.write(f"{self.cDMRs[chrom].annotation[dmr]['in_cust1']}\t")
                        if self.cDMRs[chrom].annotation[dmr]['in_cust2'] != "N/A":
                            fid.write(f"{self.cDMRs[chrom].annotation[dmr]['in_cust2']}\t")
                        if self.cDMRs[chrom].annotation[dmr]['in_gene'] != "N/A":
                            fid.write(f"{self.cDMRs[chrom].annotation[dmr]['in_gene']['present']}\t")
                            name = ", ".join(map(str, self.cDMRs[chrom].annotation[dmr]['in_gene']['name']))
                            fid.write(f"{name}\t")
                            strand = self.cDMRs[chrom].annotation[dmr]['in_gene']['strand']
                            if strand != strand:  # only happens when strand is nan
                                clear_strand = strand
                            else:
                                clear_strand = ", ".join(["+" if x == 1 or x == "1" else "-" for x in strand])
                            fid.write(f"{clear_strand}\t")
                            fid.write(f"{self.cDMRs[chrom].annotation[dmr]['in_prom']['present']}\t")
                            name = ", ".join(map(str, self.cDMRs[chrom].annotation[dmr]['in_prom']['name']))
                            fid.write(f"{name}\t")
                            strand = self.cDMRs[chrom].annotation[dmr]['in_prom']['strand']
                            if strand != strand:  # only happens when strand is nan
                                clear_strand = strand
                            else:
                                clear_strand = ", ".join(["+" if x == 1 or x == "1" else "-" for x in strand])
                            fid.write(f"{clear_strand}\t")
                            fid.write(f"{self.cDMRs[chrom].annotation[dmr]['upstream_TSS']['dist']}\t")
                            name = ", ".join(map(str, self.cDMRs[chrom].annotation[dmr]['upstream_TSS']['name']))
                            fid.write(f"{name}\t")
                            strand = self.cDMRs[chrom].annotation[dmr]['upstream_TSS']['strand']
                            if strand != strand:  # only happens when strand is nan
                                clear_strand = strand
                            else:
                                clear_strand = ", ".join(["+" if x == 1 or x == "1" else "-" for x in strand])
                            fid.write(f"{clear_strand}\t")
                            fid.write(f"{self.cDMRs[chrom].annotation[dmr]['downstream_TSS']['dist']}\t")
                            name = ", ".join(map(str, self.cDMRs[chrom].annotation[dmr]['downstream_TSS']['name']))
                            fid.write(f"{name}\t")
                            strand = self.cDMRs[chrom].annotation[dmr]['downstream_TSS']['strand']
                            if strand != strand:  # only happens when strand is nan
                                clear_strand = strand
                            else:
                                clear_strand = ", ".join(["+" if x == 1 or x == "1" else "-" for x in strand])
                            fid.write(f"{clear_strand}\t")
                    fid.write("\n")
                    
                
            
    #def parse_infile(self, infile):
     #   """Populate DMR object from text file
     #   
     #   Input: empty DMR object, file name
     #   Output: populated DMR object
      #  """      
        
     #   with open(infile, "rt") as dmrfile:
     #       for line in dmrfile:
     #           line = line.rstrip("\n") #remove trailing line feeds
      #          line = line.lstrip("\t") #remove leading tabs
      #          fields = line.split(":")  
    def plotmethylation(self, chr_name, DMR_idx, fname):
        """Creates a scatter plot of methylation by group
        
        Input: chromosome name and index of DMR of interest
        Output: png file of scatter plot
        """
        
        # get chrom index by name
        chrom = self.chromosomes.index(chr_name)
        
        # reorder samples by groups
        pos_flat = [x for y in self.groups["positions"] for x in y]
        samp_by_grp = [self.samples[x] for x in pos_flat]
        
        # get methylation values
        meth_groups = self.cDMRs[chrom].methylation[range(self.groups["no_groups"])]
        meth_groups = meth_groups[:, DMR_idx]
        meth_samples = self.cDMRs[chrom].methylation[:, DMR_idx][-self.no_samples:]
        meth_samples = [meth_samples[x] for x in pos_flat]
        
        # scatter plot
        plt.scatter(range(self.no_samples), meth_samples)
        plt.xlim(-0.5, self.no_samples-0.5)
        plt.ylim(0,1)
        plt.xticks(ticks=range(self.no_samples), labels=samp_by_grp, rotation=90)
        plt.ylabel("Methylation")
        
        # beautify
        gs = np.cumsum(np.array([len(x) for x in self.groups["positions"]]))
        gs = -0.5 + np.insert(gs, 0, 0)
        for x in range(1,len(gs)):
            if x < len(gs)-1:
                plt.vlines(gs[x],0,1)
            plt.text(0.5*(gs[x-1]+gs[x]), 1.03, self.groups["group_names"][x-1], ha="center")
            plt.hlines(meth_groups[x-1], gs[x-1], gs[x], linestyles="dotted", linewidths=1)
        #plt.show()
        plt.savefig(fname)
        
    def plot(self, DMR_chrom, DMR_idx, gc, samples, gene_bed, cgi_bed, orderby="groups", widenby=0):
        chr_idx = self.chromosomes.index(DMR_chrom)
        start = self.cDMRs[chr_idx].gen_start[DMR_idx]
        end = self.cDMRs[chr_idx].gen_end[DMR_idx]
        region = {"chrom":DMR_chrom, "start":start, "end":end}
        if orderby == "groups":
           pos_flat = [x for y in self.groups["positions"] for x in y]
           samples = [samples[x] for x in pos_flat]
        t.plot_region(region, gc, samples, gene_bed, cgi_bed, widenby)
        
    def adjust_params(self, sim_dmrs, thresh=0.05, fname="DMR_log.txt", report=True, statfile="fdr_stats.txt"):
        """Finds parameters values that achieve a desired FDR

        Input: observed DMRs
               list of simulated DMR data
               threshold for the FDR (defult:0.05)
               report    True if reporting to the display is desired
               fname     output filename (log file)
        Output:       DMR object where filters on the optimal parameters have been applied.
        """
        most_DMRs = 0 
        thresh_Qt = None
        thresh_CpG = None
        # observed number of DMRs
        obs_noDMRs = self.noDMRs()[0]
        # finding the largest value of the parameters
        with open(statfile, "w") as fid:
            fid.write("cpg\tqt\tcounter\tmean_sim_counter\tratio\n")
            for cpg in range(sim_dmrs[0].algorithm["min_CpGs"], max([max(self.cDMRs[x].no_CpGs, default=0) for x in range(len(self.cDMRs))])+1):  # why +1?
                for qt in range(sim_dmrs[0].algorithm["min_Qt"], int(np.ceil(max([max(self.cDMRs[x].max_Qt, default=0) for x in range(len(self.cDMRs))]))+1)):
                    counter = 0 
                    sim_counter = np.zeros(len(sim_dmrs))
                    # make a grid search, using jumps of 1 (currently fixed default), picking the parameters that obey 
                    # FDR<threshold and providing the largest number of DMRs
                    for chrom in range(self.no_chromosomes):
                        tot = len(np.where((np.array(self.cDMRs[chrom].no_CpGs) >= cpg) & (np.array(self.cDMRs[chrom].max_Qt) >=qt))[0])
                        counter += tot
                        i = 0
                        for dm in sim_dmrs:
                            tot = len(np.where((np.array(dm.cDMRs[chrom].no_CpGs) >= cpg) & (np.array(dm.cDMRs[chrom].max_Qt) >=qt))[0])
                            sim_counter[i] += tot
                            i += 1
                    # evaluate FDR
                    ratio = np.mean(sim_counter)/counter
                    fid.write(f"{cpg}\t{qt}\t{counter}\t{np.mean(sim_counter)}\t{ratio}\n")
                    if np.isnan(ratio):
                        continue
                    elif ratio <= thresh:
                        if counter > most_DMRs:
                            most_DMRs = counter
                            thresh_Qt = qt
                            thresh_CpG = cpg
            print(f"Qt threshold is {thresh_Qt}")
            print(f"CpG threshold is {thresh_CpG}")
            if report:
                with open(fname, "a") as fid:
                    fid.write(f"Qt threshold:{thresh_Qt}, CpG threshold:{thresh_CpG}\n")
        # recompute observed DMRs using chosen parameters    
            adjusted_dm = copy.deepcopy(self)
            #if not thresh_Qt:  # causes errors
            if thresh_Qt == None:
                print("FDR was larger than the threshold for all parameter values")
            else:
                # initalize the cDMRs object
                cdm = [c.cDMR() for i in range(self.no_chromosomes)]
                # populate the object
                for chrom in range(self.no_chromosomes):
                    # find DMRs that pass the threshold
                    chromosome = self.chromosomes[chrom]
                    idx = sorted(list(set(list(np.where(np.array(self.cDMRs[chrom].no_CpGs) >= thresh_CpG)[0])).intersection(list(np.where(np.array(self.cDMRs[chrom].max_Qt) >= thresh_Qt)[0]))))
                    if idx:
                        cdm[chrom].chromosome = chromosome
                        cdm[chrom].CpG_start = np.array(self.cDMRs[chrom].CpG_start)[idx]
                        cdm[chrom].CpG_end = np.array(self.cDMRs[chrom].CpG_end)[idx]
                        cdm[chrom].gen_start = np.array(self.cDMRs[chrom].gen_start)[idx]
                        cdm[chrom].gen_end = np.array(self.cDMRs[chrom].gen_end)[idx]
                        cdm[chrom].no_bases = np.array(self.cDMRs[chrom].no_bases)[idx]
                        cdm[chrom].no_CpGs = np.array(self.cDMRs[chrom].no_CpGs)[idx]
                        cdm[chrom].max_Qt = np.array(self.cDMRs[chrom].max_Qt)[idx]
                        cdm[chrom].methylation = np.array([np.array(self.cDMRs[chrom].methylation[x])[idx] for x in range(len(self.cDMRs[chrom].methylation))])  # will this work?
                        cdm[chrom].no_DMRs = len(idx)
                        cdm[chrom].grp_methylation_statistic = np.array(self.cDMRs[chrom].grp_methylation_statistic)[idx]
                print(f"{most_DMRs} out of {obs_noDMRs} DMRs remain after adjustment to FDR {thresh}")
                adjusted_dm.cDMRs = cdm            
            print("done")
        # add return line! 
        return adjusted_dm
        
        



if __name__ == "__main__":
    dmr = DMRs()
    print(dmr)

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
    def ancient_Newton_Raphson(max_iterations, min_tol, pi, tij, nij):
        """Solves for m when all samples are ancient
        
        Input: max_iterations    maximum number of iterations
               min_tol           convergence tolerance
               pi                deamination rates of all samples (list)
               tij               number of Ts in sample <i>, window <j>
               nij               number of reads in sample <i>, window <j>
        output: m        methylation vector
                dm       standard error of methylation
                m0       initial guess (mainly for debugging purposes)
        """
        #computer useful magnitudes
        Tj = np.nansum(tij,0)
        no_samples = len(pi)
        #change nans to 0s
        tshape = tij.shape  # get shape to reshape vector after replacing nans
        nshape = nij.shape
        tij = tij.flatten("F")  # flatten to find nans
        nij = nij.flatten("F")
        nan_t = [x for x,y in enumerate(tij) if np.isnan(y)]  # find nans
        nan_n = [x for x,y in enumerate(nij) if np.isnan(y)]
        idx = sorted(list(set(nan_t).union(set(nan_n))))  # find union of nans
        tij[idx] = 0  # replace nans
        nij[idx] = 0
        tij = np.reshape(tij, tshape, order="F")  # reshape to original vector structure
        nij = np.reshape(nij, nshape, order="F")
        #compute more params using matrix multiplication
        Tpij = pi.dot(tij)
        Npij = pi.dot(nij)
        #compute inital guess
        m = Tj / (Npij-Tpij)
        m0 = m[:]
        #make iterations
        iteration = 1
        tol = 1
        while iteration<max_iterations and tol>min_tol:
            m_prev = m[:]
            dldm = Tj/m_prev
            d2ldm2 = -Tj/m_prev**2
            for samp in range(no_samples):
                dldm = dldm - (nij[samp] - tij[samp]) * pi[samp] / (1 - m_prev * pi[samp])
                d2ldm2 = d2ldm2 - (nij[samp] - tij[samp]) * pi[samp]**2 / (1 - m_prev * pi[samp])**2
            m = m_prev - dldm/d2ldm2
            tol = max(abs(m-m_prev)/m_prev)
            iteration += 1
        #compute estimaton of the variance
        I = Tj/m**2
        for samp in range(no_samples):
            I = I + (nij[samp] - tij[samp]) * pi[samp]**2 / (1 - m * pi[samp])**2
        dm = np.sqrt(1/I)
        return (m,dm,m0)
    
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
            CpGs_inDMR = int(np.where(iQt[start:end]==maxQt)[0]) +1  # to get number of elements, rather than index
            
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
                idm.methylation.append(meth)
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
                region = f"{cdmr.chromosome}:{cdmr.gen_start[dmr]}-{cdmr.gen_end[dmr]}"
                reg_std = t.standardize_region(region)
                for samp in range(no_samples):
                    samp_meth[samp,dmr] = samples[samp].region_methylation(reg_std, coord)
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

    def groupDMRs(self, samples=[], sample_groups=[], coord=[], d_rate_in=[], chroms=[], fname="groupDMRs.txt", win_size="meth", lcf="meth", delta=0.5, min_bases=100, min_meth=0, max_meth=1, min_Qt=0, min_CpGs=10, max_adj_dist=1000, k=2, min_finite=1, max_iterations=20, tol=1e-3, report=True, real=False, no_permutations=None):
        """Detects DMRs between two groups of samples
        
        Input: samples            list of sample (Amsample or Mmsample) objects
               sample_groups      list of sample group names
               coord              Gcoordinates object with CpG coordinates
               d_rate_in          deamination rates for ancient samples (if empty, taken from values in sample)
               chroms             list of chromsome names
               fname              output file name
               win_size           window size for smoothing. If 'meth', it is taken as the value used to reconstruct 
                   the methylation in each sample. If 'auto', a recommended value is computed for every chromosome
                   of each sample. Otherwise, it can be a scalar (used for all chromosomes in all samples), a vector 
                   over the samples (same window size is used for all chromosomes of each ancient individual, nan is
                   substituted for each modern individual), or a 2d array with values per individual and chromosome
               lcf                low coverage factor. If 'meth', it is taken as the value used in reconstructing the
                   methylation of each sample
               delta              minimum methylation difference between the two groups
               min_bases          the minimum length of each DMR in bases--shorter DMRs are filtered out
               min_meth           sets a lower bound to the methylation levels of the reference, such that in 
                   every position ref_meth = max(ref_meth,ref_min_meth)
               max_meth           sets an upper bound to the methylation levels of the reference, such that in 
                   every position ref_meth = min(ref_meth,ref_max_meth)
               min_Qt             DMRs with Qt < min_Qt are filtered out
               min_CpGs           DMRs whose number of CpGs is less than min_CpGs are filtered out
               max_adj_dist       max distance between adjacent CpGs within the same DMR. If the distance between
                   consecutive CpG positions is larger than max_adj_dist, the algorithm sets Qt to 0
               k                  minimum number of standard errors separating the groups
               min_finite         an array of length no_groups stating the minimum number of ancient samples for 
                   which we require data. If in a position there are not enough samples with data, a NaN is 
                   substituted in this position. It can also be a fraction between 0 and 1, in which case it is 
                   understood as the minimum fraction of the total number of ancient samples in the group
               max_iterations     maximum number of iterations in the Newton-Raphson phase
               tol                tolerance in the Newton-Raphson phase
               report             True if reporting to the display is desired 
               real               False if the analysis is run on simulations
               no_permutations    the number of permutation in case of simulation
        Output: modified DMR object, Qt_up, Qt_down
        """
        no_samples = len(samples)
        is_ancient = [1 if type(x).__name__ == "Amsample" else 0 for x in samples]
        chromosomes = chroms if chroms else coord.chromosomes
        if win_size == "meth":
            is_auto_win = False
            is_meth_win = True
        elif win_size == "auto":
            is_auto_win = True
            is_meth_win = False
        else:
            is_auto_win = False
            is_meth_win = False
        is_meth_lcf = True if lcf == "meth" else False
        lcf = np.zeros(no_samples)
        win_modern = 11
        d_rate = np.zeros(no_samples)
        d_rate[:] = np.nan
        if d_rate_in:
            if len(d_rate_in) == sum(is_ancient):
                y = 0
                for x in range(no_samples):
                    if is_ancient[x] == 1:
                        d_rate[x] = d_rate_in[y]
                        y += 1
            else:
                raise Exception(f"Length of d_rate ({len(d_rate_in)}) does not match number of ancient samples ({sum(is_ancient)})")
        #no_groups = len(groups)
        #min_finite = np.ones(no_groups)
        #self.real = real
        if real == False:
            is_sim = True
            no_simulations = no_permutations
        else:
            is_sim = False
        #guarantee methylation values are in the range [0,1]
        if min_meth > 1:
            min_meth = 0.01 * min_meth
        if max_meth > 1:
            max_meth = 0.01 * max_meth
        #substitute default drates
        for samp in range(no_samples):
            if is_ancient[samp] and np.isnan(d_rate[samp]):
                d_rate[samp] = samples[samp].d_rate["rate"]["global"]
        #process low-coverage-filter
        if is_meth_lcf:
            for samp in range(no_samples):
                if is_ancient[samp]:
                    lcf[samp] = samples[samp].methylation["lcf"][0]
        #process window size
        no_chr = len(chromosomes)
        if not is_auto_win:
            if is_meth_win:
                win_size = np.nan * np.ones((no_samples, no_chr))
                for samp in range(no_samples):
                    if is_ancient[samp]:
                        indices = samples[samp].index(chromosomes)
                        smp = samples[samp].methylation["win_size"]
                        win_size[samp,] = [smp[x] for x in indices]
                    else:
                        win_size[samp,] = win_modern
            else:
                #Option 1: same window size for all individuals/chromosomes
                if len(win_size) == 1:
                    if not win_size%2: #win_size is even
                        win_size += 1 #make it odd
                    win_size = win_size * np.ones((no_samples, no_chr))
                    win_size[np.array(is_ancient)==0,] = win_modern
                #same W for all chromosomes of an individual
                elif win_size.ndim ==1:
                    for samp in range(no_samples):
                        if ~np.isnan(win_size[samp]):
                            if not win_size[samp]%2: #win_size is even
                               win_size[samp] += 1 #make it odd
                    win_size = np.transpose([win_size]) * np.ones(no_chr)
                    win_size[np.array(is_ancient)==0,] = win_modern
                #different W for each chromosome and individual
                else:
                    for samp in range(no_samples):
                        for chrom in range(no_chr):
                            if ~np.isnan(win_size[samp, chrom]):
                                if not win_size[samp, chrom]%2: #win_size is even
                                    win_size[samp, chrom] += 1 #make it odd
                    win_size[np.array(is_ancient)==0,] = win_modern
        else:
            win_size = np.zeros(no_samples, no_chr)
            for samp in range(no_samples):
                if is_ancient[samp]:
                    for chrom in range(no_chr):
                        chr_ind = samples[samp].index(chromosomes[chrom]) #enough just to use chrom?
                        win_size[samp, chrom] = samples[samp].determine_winsize(chr_ind, **winsize_algorithm)
            win_size[np.array(is_ancient)==0,] = win_modern
        #self.win_size = win_size
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
            no_rows = max(group_sizes)[1]
            for row in range(no_rows):
                line = "|"
                idx = [np.nan for x in group_names]
                for grp in range(no_groups):
                    pos = group_nums.index(grp+1) if grp+1 in group_nums else ""
                    idx[grp] = pos if pos or pos == 0 else np.nan
                    if np.isnan(idx[grp]):
                        word = ""
                    else:
                        if is_ancient[idx[grp]]:
                            word = f"{samples[idx[grp]].name} ({100*d_rate[idx[grp]]:.2f}%)"
                        else:
                            word = samples[idx[grp]].name
                        group_nums[idx[grp]] = np.nan
                    line += word.center(int(max_len[grp])) + "|"
                fid.write(f"\t{line}\n")
            fid.write(f"\t{sep}\n\n")
            #write params of job
            #delta
            fid.write(f"delta = {delta:.2f}\n\t[used to compute the lt statistics]\n")
            #k
            fid.write(f"k = {k:.1f}\n\t[used to compute the lt statistics]\n")
            #min_bases
            fid.write(f"min_bases = {min_bases}\n\t[minimum length of a DMR (bases). Shorter DMRs are filtered out]\n")
            #max_adj_dist
            fid.write(f"max_adj_dist = {max_adj_dist}\n\t[max distance between adjacent CpGs within the same DMR (bases)]\n")
            #min_Qt
            fid.write(f"min_Qt = {min_Qt:.2f}\n\t[a DMR must have Qt >= min_Qt]\n")
            #min_CpGs
            fid.write(f"min_CpGs = {min_CpGs}\n\t[a DMR must contain at least min_CpGs CpGs]\n")
            #min_meth and max_meth
            fid.write(f"min_meth = {min_meth:.2f}\n\t[lower bound to methylation levels of the reference, ref_meth = max(ref_meth, min_meth)]\n")
            fid.write(f"max_meth = {max_meth:.2f}\n\t[upper bound to methylation levels of the reference, ref_meth = min(ref_meth, max_meth)]\n")
            #min_finite
            fid.write(f"min_finite = {min_finite}\n\t[minimum number of ancient samples per group for which we require data]\n")
            #max_iterations
            fid.write(f"max_iterations = {max_iterations}\n\t[maximum number of iterations in Newton-Raphson]\n")
            #tol
            fid.write(f"tol = {tol}\n\t[convergence tolerance of Newton-Raphson]\n\n")
            #fid.close()
        #initializations
        no_chr = len(chromosomes)
        iSa = [x for x,y in enumerate(is_ancient) if y == 1]
        iSm = [x for x,y in enumerate(is_ancient) if y == 0]
        Qt_up = [None]*no_chr
        Qt_down = [None]*no_chr
        for chrom in range(no_chr):
            Qt_up[chrom] = [np.nan] * len(coord.coords[coord.index([chromosomes[chrom]])[0]])
            Qt_down[chrom] = Qt_up[chrom]
        samp_names = [None]*no_samples
        species = [None]*no_samples
        reference = samples[0].reference
        for samp_ind in range(no_samples):
            samp_names[samp_ind] = samples[samp_ind].name
            species[samp_ind] = samples[samp_ind].species
            if samples[samp_ind].reference != reference:
                raise Exception(f"Sample {samples[samp_ind].name} does not use {reference}")
            #make sure modern samples are scaled (range 0-1)
            if not is_ancient[samp_ind]:
                samples[samp_ind].scale()
        self.samples = samp_names
        self.is_ancient = is_ancient
        self.chromosomes = chromosomes
        self.species = species
        self.reference = reference
        #split samples between groups
        if no_groups != 2:
            raise Exception("Currently, the algorithm works only on 2 groups")
        giSa = [None]*no_groups
        giSm = [None]*no_groups
        gSa = np.zeros(no_groups)
        gSm = np.zeros(no_groups)
        for grp in range(no_groups):
            giSa[grp] = list(set(iSa) & set(positions[grp]))
            giSm[grp] = list(set(iSm) & set(positions[grp]))
            gSa[grp] = len(giSa[grp])
            gSm[grp] = len(giSm[grp])
        if not all(gSa+gSm):
            raise Exception("Groups may not be empty")
        #loop on chroms
        cdm = [c.cDMR() for i in range(no_chr)]
        for chrom in range(no_chr):
            #report
            if report:
                print(f"Processing chromosome {chromosomes[chrom]}")
                fid.write(f"Processing chromosome {chromosomes[chrom]}\n")
            #find matching chroms in diff samples
            #sample_chrom_idx is the index of the current chrom in each sample
            sample_chrom_idx = np.zeros(no_samples)
            for samp in range(no_samples):
                sample_chrom_idx[samp] = samples[samp].index([chromosomes[chrom]])[0]
            #get number of positions along the chrom
            no_pos = len(coord.coords[chrom])
            #loop on groups
            meth = np.zeros((no_groups, no_pos))  # methylation
            meth_err = np.zeros((no_groups, no_pos))  # standard error in methylation
            for grp in range(no_groups):
                #for ancient--compute tij and nij
                tij = np.zeros((int(gSa[grp]), no_pos))
                nij = np.zeros((int(gSa[grp]), no_pos))
                for samp in range(int(gSa[grp])):
                    samp_id = giSa[grp][samp]
                    [tij[samp], nij[samp]] = samples[samp_id].smooth(int(sample_chrom_idx[samp_id]), int(win_size[samp_id,chrom]))
                    #nij[samp] = win_size[samp_id,chrom]
                    #remove regions with very low coverage
                    lct = t.find_low_coverage_thresh(nij[samp], lcf[samp])
                    nij[samp][nij[samp]<lct] = np.nan
                #for modern--compute mij_bar and wij
                mij_bar = np.zeros((int(gSm[grp]), no_pos))
                wij = np.zeros((int(gSm[grp]), no_pos))
                for samp in range(int(gSm[grp])):
                    samp_id = giSm[grp][samp]
                    [mij_bar[samp], wij[samp]] = samples[samp_id].smooth(int(sample_chrom_idx[samp_id]), int(win_size[samp_id,chrom]))
                    #wij[samp] = win_size[samp_id,chrom]
                #estimate methylation
                if gSa[grp] > 0:
                    #ancient samples
                    [ma, dma, m0] = self.ancient_Newton_Raphson(max_iterations, tol, d_rate[giSa[grp]],tij,nij)
                if gSm[grp] > 0:
                    #modern samples
                    Wj = np.nansum(wij,0)
                    mm = sum(wij * mij_bar, 1)/Wj
                    dmm = np.sqrt(1/Wj)
                #compute final estimator
                if gSm[grp] * gSa[grp] > 0:
                    #both modern and ancient
                    meth[grp] = (gSm[grp]*mm + gSa[grp]*ma) / (gSm[grp] + gSa[grp])
                    meth_err[grp] = 1/(gSm[grp] + gSa[grp]) * np.sqrt(gSm[grp]**2 * dmm**2 + gSa[grp]**2 * dma**2)
                elif gSa[grp] > 0:
                    #only ancient samples
                    meth[grp] = ma
                    meth_err[grp] = dma
                else:
                    #only modern samples
                    meth[grp] = mm
                    meth_err[grp] = dmm
                #substitute nans when there are too many nans in orig data
                finit = np.isfinite(tij) & np.isfinite(nij)
                finit = finit.sum(axis=0)  # this does not work for 1d array
                meth[grp, finit<min_finite[grp]] = np.nan  # any chance this will work?
            # compute the two statistics
            meth[meth>1] = 1
            diffi = meth[0] - meth[1]  # since there must be exactly 2 groups
            idm = np.sqrt(meth_err[0]**2 + meth_err[1]**2)
            lt_up = (diffi - delta)/idm
            lt_down = (-diffi - delta)/idm
            not_nans = np.isfinite(lt_up)
            lt_up = lt_up[not_nans]
            lt_down = lt_down[not_nans]
            coordi = coord.coords[coord.index([chromosomes[chrom]])][0]
            coordi = coordi[not_nans]
            methi = meth[:,not_nans]
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
            # set idx_binary=1 only for positions whose reference methylation
            # in group #2 is below 1-delta (otherwise, clearly there is no DMR)
            idx_binary = np.zeros(num_finite_pos)
            idx_binary[methi[1,:]<1-delta] = 1
            # TODO: test the effect of idx_binary = ones(num_finite_pos,1);
            # compute {iQt_up} recursively
            for pos in range(num_finite_pos):
                iQt_up[pos+1] = max(0, coordi_diff[pos] * idx_binary[pos] * (iQt_up[pos] + lt_up[pos]))
            # initialize {iQt_down}
            iQt_down = np.zeros(num_finite_pos+1)
            # set idx_binary=1 only for positions whose methylation in group #2
            # is above delta
            #TODO: test idx_binary = ones(no_finite_positions,1);
            idx_binary = np.zeros(num_finite_pos)
            idx_binary[methi[1,:]>delta] = 1
            #idx_binary = ones(no_finite_positions,1);
            # compute {iQt_down} recursively
            for pos in range(num_finite_pos):
                iQt_down[pos+1] = max(0, coordi_diff[pos] * idx_binary[pos] * (iQt_down[pos] + lt_down[pos]))
            # make a clean version of {Qt}, that takes into account the many CpG
            # positions we had removed earlier. The function reports this vector,
            # which is useful later for plotting
            Qt_up[chrom] = np.array(Qt_up[chrom])
            Qt_up[chrom][not_nans] = iQt_up[1:]
            Qt_down[chrom] = np.array(Qt_down[chrom])
            Qt_down[chrom][not_nans] = iQt_down[1:]
            # filter and characterize DMRs
            self.findDMRs(cdm[chrom], iQt_up, coordi, trunc2clean, meth, min_bases, min_CpGs, min_Qt)
            self.findDMRs(cdm[chrom], iQt_down, coordi, trunc2clean, meth, min_bases, min_CpGs, min_Qt)
            cdm[chrom].chromosome = chromosomes[chrom]
            # compute methylation in each sample
            samp_meth = np.zeros((len(samples),cdm[chrom].no_DMRs))
            samp_meth = self.get_region_meth(cdm[chrom], no_samples, samples, samp_meth, coord)
            meth = np.array(cdm[chrom].methylation)
            meth = np.transpose(meth)
            samp_meth = np.insert(samp_meth,0,meth,axis=0)
            cdm[chrom].methylation = samp_meth
            if report:
                print(f"\tdetected {cdm[chrom].no_DMRs} DMRs")
                fid.write(f"\tdetected {cdm[chrom].no_DMRs} DMRs\n")
                
        #substitue fields
        self.algorithm = "groupDMRs"
        self.cDMRs = cdm
        self.no_chromosomes = no_chr
        self.no_samples = no_samples
        
        #close filehandle
        if report:
            fid.close()
        
        return(Qt_up, Qt_down)
    
    def annotate(self, gene_bed, cgi_bed, prom_def=[5000,1000]):
        """Retrieves important data about each DMR
        
        Input: gene_bed   bed file with gene data 
               cgi_bed    bed file with CGI data
               prom_def   promoter definition around TSS, a list of 2 values [before after], where before is the
                   number of nucleotides into the intergenic region, and after is the number of nucleotides 
                   into the gene
        Output: cDMR object with updated annotation
        """
        genes = pbt.BedTool(gene_bed)
        #genes_no_dups = genes.groupby(g=[1,2,3,6], c='4,5', o='distinct').cut([0,1,2,4,5,3], stream=False)
        genes_no_dups = genes.groupby(g=[1,2,3,6], c='4,5', o='distinct').cut([0,1,2,4,5,3])  # is stream nec?
        cgis = pbt.BedTool(cgi_bed)
        #cgis_no_dups = cgis.groupby(g=[1,2,3,6], c='4,5', o='distinct').cut([0,1,2,4,5,3], stream=False)
        before = prom_def[0]
        after = prom_def[1]
        proms = gint.Gintervals(chr_names=self.chromosomes)
        proms.calc_prom_coords(genes_no_dups, before, after)
        tss = gcoord.Gcoordinates(chr_names=self.chromosomes, description="TSS positions")
        tss.calc_tss(genes_no_dups)
        # loop on chromosomes
        for chrom in range(self.no_chromosomes):
            print(f"chrom {self.chromosomes[chrom]}")
            num_DMRs = self.cDMRs[chrom].no_DMRs
            # add handling for no DMRs?
            # get DMR regions
            regions = self.get_regions(self.cDMRs[chrom])
            
            # generate an array of TSSs
            # add handling of no tss?
            itss = tss.coords[chrom]
            itss = itss.astype(float)
            prom_string = ""
            for prom in range(len(proms.start[chrom])):
                prom_string += f"\n{proms.chr_names[chrom]} {proms.start[chrom][prom]} {proms.end[chrom][prom]} {proms.iname[chrom][prom]} 0 {proms.strand[chrom][prom]}"
            chrom_proms = pbt.BedTool(prom_string, from_string=True)
            dmr_annot = []
            #get genes in this chrom
            genes_chrom = genes_no_dups.filter(lambda a: a.chrom == str(chrom+1)).saveas()
            for dmr in range(num_DMRs):
                annot = {}
                region = regions[dmr]
                ivl = pbt.BedTool(region, from_string=True)[0]  # get 1st element to make it an interval object
                if cgis.any_hits(ivl):
                    in_CGI = True
                else:  
                    in_CGI = False
                
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
                            in_prom["strand"].append(hit.strand)
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
                    
                annot["in_CGI"] = in_CGI
                annot["in_gene"] = copy.deepcopy(in_gene)
                annot["in_prom"] = copy.deepcopy(in_prom)
                annot["upstream_TSS"] = copy.deepcopy(upstream_TSS)
                annot["downstream_TSS"] = copy.deepcopy(downstream_TSS)
                
                dmr_annot.append(copy.deepcopy(annot))
            
            self.cDMRs[chrom].annotation = (copy.deepcopy(dmr_annot))
       
    def dump_DMR(self):
        """Dumps DMR object to text file.
        
        Input: DMR object
        Output: text file in format DMRs_<time>.txt (directory currently hard-coded).
        """
        time = datetime.datetime.now()
        time = time.strftime("%d-%m-%Y_%H.%M")
        fname = f"data/python_dumps/DMRs_{time}.txt"
        with open(fname, "w") as fid:
            fid.write(f"Samples:\n\t{self.samples}\n")
            fid.write(f"Groups:\n\tAssignment: {self.groups['group_nums']}\n\tNaming: {self.groups['group_names']}\n")
            fid.write(f"Species:\n\t{self.species}\n")
            fid.write(f"Ancient samples:\n\t{self.is_ancient}\n")
            fid.write(f"Reference:\n\t{self.reference}\n")
            fid.write(f"Chromosomes:\n\t{self.chromosomes}\n")
            fid.write(f"Algorithm:\n\t{self.algorithm}\n")
            fid.write("cDMRs:\n")
            for chrom in range(self.no_chromosomes):
                fid.write(f"\tChrom {self.chromosomes[chrom]}:\n")
                fid.write(f"\t\tNumber of DMRs: {self.cDMRs[chrom].no_DMRs}\n")
                fid.write(f"\t\tCpG starts:\n\t\t\t{self.cDMRs[chrom].CpG_start}\n")
                fid.write(f"\t\tCpG ends:\n\t\t\t{self.cDMRs[chrom].CpG_end}\n")
                fid.write(f"\t\tNumber of CpGs:\n\t\t\t{self.cDMRs[chrom].no_CpGs}\n")
                fid.write(f"\t\tGenomic starts:\n\t\t\t{self.cDMRs[chrom].gen_start}\n")
                fid.write(f"\t\tGenomic ends:\n\t\t\t{self.cDMRs[chrom].gen_end}\n")
                fid.write(f"\t\tNumber of bases:\n\t\t\t{self.cDMRs[chrom].no_bases}\n")
                fid.write(f"\t\tMax Qts:\n\t\t\t{self.cDMRs[chrom].max_Qt}\n")
                fid.write("\t\tMethylation:\n")
                for meth in range(len(self.cDMRs[chrom].methylation)):
                    fid.write(f"\t\t\t{self.cDMRs[chrom].methylation[meth]}\n")
                fid.write(f"\t\tAnnotation:\n\t\t\t{self.cDMRs[chrom].annotation}\n")
                
            
            



if __name__ == "__main__":
    dmr = DMRs()
    print(dmr)

#!/usr/bin/python3

import numpy as np
import pickle
import scipy.stats as stats
import re
import sys
import pybedtools as pbt
import matplotlib.pyplot as plt
import matplotlib.axes as ax
import matplotlib.patches as p

from screeninfo import get_monitors

"""This module contains helper functions used in the RoAM process.
"""
def nanmerge(arr, operation):
    """Merges the two adjacent positions of each CpG. The rules of merger are as follows:
        (1) nan + nan = nan
        (2) value + nan = value
        (3) nan + value = value
        (4) value + value = factor*(value+value)
        The variable {factor} is 0.5 if {operation} is 'average', or 1 if {operation} is 'sum'.
        In the case of "max", merge just takes the larger value of the 2 positions.
    
    Input: arr        list with CpG positions for merging
           operation  can be "average", "sum" or "max"
    Output: merged list
    """
    result = []
    if operation == "average":
        factor = 0.5
    elif operation == "sum":
        factor = 1
    elif operation == "max":
        for i in range(0,len(arr),2):
            take_me = max(arr[i],arr[i+1])
            result.append(take_me)
        return result
    for i in range(0,len(arr),2):
        if np.isnan(arr[i]):
            if np.isnan(arr[i+1]):
                result.append(np.nan) #NaN + NaN = NaN
            else:
                result.append(arr[i+1]) #NaN + <num> = <num>
        else:
            if np.isnan(arr[i+1]):
                result.append(arr[i]) #<num> + NaN = <num>
            else:
                merged = (arr[i]+arr[i+1])*factor #num + num = (num + num) * factor. yields float result, with
                result.append(merged)             #sometimes strange precision (eg 0.8200000000000001)
    return result

def standardize_region(region):
    """Takes a region delimited by any combination of spaces, colons, and dashes and standardizes its form into
         a dictionary with keys chrom, start, and end. Commas are also removed from numbers.
         
    Input: region
    Output: dictionary containing start, end, chrom.
    """
    import re
    if isinstance(region, dict):
        return region
    fields = re.split(" +|:+|-+", region)
    chrom = fields[0]
    start = fields[1]
    end = fields[2]
    start = start.replace(",","")
    end = end.replace(",","")
    start = int(start)
    end = int(end)
    std_region = {
            "chrom":chrom,
            "start":start,
            "end":end
            }
    return std_region

def save_object(filename, obj):
    """Pickles object for future use.
    
    Input: filename    path for output file
           obj         object for picking
    Output: pickled object at requested path
    """
    outfile = open(filename, "wb") #open file for writing binaries
    pickle.dump(obj, outfile, -1) #serialize object using highest protocol
    outfile.close()

def load_object(filename):
    """Loads picked object.
    
    Input: filename    path of pickled object file
    output: unpickled object
    """
    infile = open(filename, "rb") #open file for reading binaries
    obj = pickle.load(infile)
    infile.close()
    return obj

    
def nansmooth(vec, winsize, mode):
    """Averages along a moving window ignoring NaNs.
    
    Input: vec       input list (methylation values)
           winsize   window size for smoothing
           mode      mode for convolve function
    Output: res             smoothed array
            zero_for_nan    number of non-Nans in each window
    """
    tpl = np.ones(winsize)
    (nans, zero_for_nan) = get_zeros(vec, tpl, mode)
    res = do_conv(vec, nans, tpl, mode, zero_for_nan)
    return(res, zero_for_nan)

def nanconv(vec, tpl, mode):
    """Convolution ignoring NaNs.
    
    Input: vec    input list (methylation values)
           tpl    template array (ones)
           mode   mode for convolve function
    Output: convoluted array
    """
    (nans, zero_for_nan) = get_zeros(vec, tpl, mode)
    zeros = zero_for_nan>0
    res = do_conv(vec, nans, tpl, mode, zeros)
    return res

def get_zeros(vec, tpl, mode):
    """Produces an array with zeros where NaNs should be present in the output, and ones elsewhere
    
    Input: vec    input list (methylation values)
           tpl    template array (ones)
           mode   mode for convolve function
    Output: nans            array with True where NaNs are in input list
            zero_for_nan    convolution of modified array (with zeros for NaNs) with template
    """
    vec = np.array(vec)
    vec_len = len(vec)
    mod_vec = np.ones(vec_len)
    nans = np.isnan(vec)
    mod_vec[nans] = 0
    zero_for_nan = np.convolve(mod_vec, tpl, mode)
    return (nans, zero_for_nan)

def do_conv(vec, nans, tpl, mode, zero_for_nan):
    """Runs the convolve function
    
    Input: vec            input list (methylation values)
           nans           array with True where NaNs are in input list
           tpl            template array (ones)
           mode           mode for convolve function
           zero_for_nan   convolution of modified array (with zeros for NaNs) with template
    Output: convoluted array
    """
    vec = np.array(vec)
    vec[nans] = 0
    res = np.convolve(vec, tpl, mode)/zero_for_nan #raises divide by zero runtime warnings
    return res

"""The following three functions are used to find the C+T percentiles

Input: arr    array with C+T for each position
       perc   requested percentiles (list)
Output: percentile values (list)
"""
def quantile(x,q):
    n = len(x)
    y = np.sort(x)
    return(np.interp(q, np.linspace(1/(2*n), (2*n-1)/(2*n), n), y))

def prctile(x,p):
    return(quantile(x,np.array(p)/100))

def nan_prctile(arr,perc):
    no_nan = [x for x in arr if ~np.isnan(x)]
    return(prctile(no_nan, perc))

def pop_dist(p,N):
    """ internal function used by bmm"""
    no_dists = len(p)
    B = np.zeros((no_dists, N+1)) 
    vals = range(N+1)
    for x in range(no_dists):
        B[x,] = stats.binom.pmf(vals,N,p[x])
    return B

def log_likelihood(w,B,H):
    """ internal function used by bmm"""
    mat_prod = np.matmul((B+np.spacing(1)).T,w)
    llike = np.matmul(H[0],np.log(mat_prod)) 
    return llike

def bmm(H, p, w, tolerance, p_known):
    """Carries out parameter estimation of binomial mixture model (BMM) using EM algorithm. The model is
        Pr(y_i=k) = sum_k w_k Binom(k|N,p_k) and the algorithm estimates the parameters w_k and p_k.
    
    Input: H            histogram of the observed success counts
           p            initial guess of the vector p (list)
           w            initial guess of the vector w (list)
           tolerance    convergence tolerance (of the loglikelihood)
           p_known      a binary list whose length is the number of mixed distribution, that should be used in the case 
            that some of the p values are known in advance. In this case, the known values should be used in p, and they
            will be identified by those elements in p_known that are true (1)
    Output: p        estimated vector p
            w        estimated vector w
            llike    vector of log-likelihood computed in each of the EM iterations. It can monitor the number of 
            iterations and can be used to verify that the likelihood always increases. Typical values are 1e-3.
    """
    if p_known:
        p_known = np.where(p_known)[0]
        p_val = np.array(p)[p_known]
    else:
        p_known = []
        p_val = np.nan

    #compute basic magnitudes
    N = len(H[0])-1 #N of binomial dist
    n = sum(H[0]) #total number of observations
    G = len(w) #number of mixed distributions
    obs = [x for x in range(N+1)] #list of all poss observations

    #initialize
    B = pop_dist(p,N)
    llike = log_likelihood(w,B,H)
    llike = [llike*2, llike]

    #EM iterations
    while (llike[-2] - llike[-1])/llike[-2] > tolerance:
        #update rgk
        nom = ((np.array([w]).T)*np.ones(N+1))*B #this will work only if dimensions match
        #is above numerator?
        denom = np.array([np.ones(G)]).T*sum(nom)
        rgk = nom/denom
        #update w
        w = 1/n * np.matmul(rgk,H[0].T) 
        #update p
        nom = np.matmul(rgk,(H[0]*obs).T)
        denom = N * np.matmul(rgk,H[0].T)
        p = nom/denom
        p[p_known] = p_val
        #update matrix B
        B = pop_dist(p,N)
        #update log-likelihood
        llike_new = log_likelihood(w,B,H)
        llike = np.append(llike, llike_new) #should be llike+1 elements long
    p = [min(max(0,x),1) for x in p] #make sure p is bet 0 and 1
    return p, w, llike

def nan_eq_2d(list1, list2):
    for x in range(len(list(list1))):
        for y in range(len(list(list1[x]))):
            if list1[x][y] != list2[x][y] and (~np.isnan(list1[x][y]) or ~np.isnan(list2[x][y])):
                print(f"list1({x},{y}) is {list1[x][y]} and list2({x},{y}) is {list2[x][y]}")
                return
    print("lists are the same")

def find_low_coverage_thresh(vec, factor=0.05):
    """Finds a low-coverage threshold. Assumes that the cutoff should be taken as factor*nanmean(vec)
    
    Input: vec    array of values
           factor parameter of the algorithm
    Output: lct    the threshold
    """
    lct= np.ceil(factor * np.nanmean(vec))
    return lct

def build_chr_key(all_chroms):
    """Builds index of chromosomes in bam list, including name variants
    
    Input: list of chromosomes in bam file
    Output: dict keyed by chrom names (and variants) with index as value
    """
    chr_key = {}
    for chrom in all_chroms:
        if "chr" in chrom:
            chr_key[chrom] = all_chroms.index(chrom)
            num = chrom[3:]
            chr_key[num] = all_chroms.index(chrom)
            num = chrom[3:]
            dig = re.search("^\d+", num)
            if not dig:
                chr_key[num.upper()] = all_chroms.index(chrom)
                chr_key[num.lower()] = all_chroms.index(chrom)
                chr_key[num.capitalize()] = all_chroms.index(chrom)
        else:
            chr_key[chrom] = all_chroms.index(chrom)
            chr_key["chr"+chrom] = all_chroms.index(chrom)
            dig = re.search("^\d+", chrom)
            if not dig:
                chr_key[chrom.lower()] = all_chroms.index(chrom)
                chr_key[chrom.upper()] = all_chroms.index(chrom)
                chr_key[chrom.capitalize()] = all_chroms.index(chrom)
                chr_key["chr"+chrom.lower()] = all_chroms.index(chrom)
                chr_key["chr"+chrom.upper()] = all_chroms.index(chrom)
                chr_key["chr"+chrom.capitalize()] = all_chroms.index(chrom)
    return chr_key

def determine_shared_winsize(samples, chrom, coverage=[], drate=[], min_meth=0.2, p0=0.01, max_width=31):
        """Estimates the optimal window size in cases that collecting data from a window is required.
        
        Input: samples      a list of amSamples
               chrom        name of chromosome
               coverage     effective coverage of the chromosome in each of the samples.
               drate        demaination rate of each sample.
               min_meth     minimum methylation level we want to detect.
               p0           the probability to get zero counts in the window if the methylation is min_meth.
               max_width    maximum allowed width. Computed window is not allowed to be larger than 'max_width'.
               
        Output: win_size    recommended window size for each sample, forced to be an odd number.
        """
        no_samples = len(samples)
        if len(coverage) != no_samples:
            print(f"inserted coverage vector (length {len(coverage)}) does not match the number of samples ({no_samples}).")
            print("Ignoring user input.")
            #coverage = np.nan*np.ones(no_samples)
        for samp in range(no_samples):
            chr_num = samples[samp].chr_names.index(chrom)
            coverage[samp] = samples[samp].diagnostics["effective_coverage"][chr_num]
        #else:
        if len(drate) != no_samples:
            print(f"inserted drate vector (length {len(drate)}) does not match the number of samples ({no_samples}).")
            print("Ignoring user input.")
            #drate = np.nan*np.ones(no_samples)
        for samp in range(no_samples):
            drate[samp] = samples[samp].d_rate["rate"]["global"]
        if min_meth>1:
            min_meth = 0.01 * min_meth
        
        #compute the window size
        cvg = sum(coverage)
        drate = max(drate)
        p = min_meth * drate
        win_size = np.ceil(np.log(p0)/np.log(1-p)/cvg)
        #narrow window if too wide
        win_size = min(win_size, max_width)

        #make sure win_size is odd
        if not win_size%2:
            win_size += 1
        
        return win_size
    
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
    no_samples = len(pi)
    #change nans to 0s
    tshape = tij.shape  # get shape to reshape vector after replacing nans
    nshape = nij.shape
    tij = tij.flatten("F")  # flatten to find nans
    nij = nij.flatten("F")
    nan_t = np.where(np.isnan(tij))[0]  # find nans
    nan_n = np.where(np.isnan(nij))[0]
    idx = sorted(list(set(nan_t).union(set(nan_n))))  # find union of nans
    tij[idx] = 0  # replace nans
    nij[idx] = 0
    tij = np.reshape(tij, tshape, order="F")  # reshape to original vector structure
    nij = np.reshape(nij, nshape, order="F")
    #compute more params using matrix multiplication
    Tpij = pi.dot(tij)
    Npij = pi.dot(nij)
    #compute Tj
    Tj = np.nansum(tij,0)
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
        m = m_prev - np.nan_to_num(dldm/d2ldm2)
        tol = np.nanmax(abs(m-m_prev)/m_prev)
        iteration += 1
    #compute estimaton of the variance
    I = Tj/m**2
    for samp in range(no_samples):
        I = I + (nij[samp] - tij[samp]) * pi[samp]**2 / (1 - m * pi[samp])**2
    dm = np.sqrt(1/I)
    return (m,dm,m0)

def pooled_methylation(samples, chroms, win_size="auto", winsize_alg={}, lcf=0.05, drate=[], min_finite=1, match_histogram=True, ref=None, ref_winsize=0, max_iterations=20, tol=0.001):
    """Computes pooled ancient methylation.
    Input: samples            list of amSamples. Methylation is computed based on pooling of information from all these 
                samples.
           chroms             names of chromosomes (list)
           win_size           window size for smoothing. If 'auto', a recommended value is computed for each chromosome.
                Otherwise, it can be a scalar (used for all chromosomes) or a vector over chromosomes
           winsize_algorithm  a dictionary with parameters required to determine window size
           lcf                low coverage factor.
           drate              deamination rates of each of the ancient samples. If provided, its length should match the
                number of ancient samples.
           min_finite         determines the minimum number of ancient samples for which we require data. The methylation
                of a position with fewer informative samples would be NaN. if 'min_finite' is a fraction between 0 
                and 1, it is interpreted as the minimum fraction of ancient samples for which we require data. If
                min_finite is an integer > 1, it is interpreted as the minimum number of ancient samples for which
                we require data
           match_histogram    determines whether to match the histogram to that of a provided reference
           ref                an mmSample representing a reference methylome, used for the match_histogram option.
           ref_winsize        window size used for the reference histogram.
           max_iterations     maximum number of iteration in the Newton-Raphson phase
           tol                tolerance in the Newton-Raphson phase
    Output: m                 pooled methylation. An  array, with a single entry per chromosome.
            dm                standard error of the pooled methylation. An array matching the dimensions of m.
    """
    no_samples = len(samples)
    no_chroms = len(chroms)
    is_auto_win = False
    if no_chroms == 0:
        chroms = samples[0].chr_names
        no_chroms = samples[0].no_chrs
    
    if len(drate) == 0:
        drate = np.nan*np.ones(no_samples)
        for samp in range(no_samples):
            drate[samp] = samples[samp].d_rate["rate"]["global"]
    
    if type(win_size) is str and win_size == "auto":  # error without type checking
        is_auto_win = True
        
    # bring parameters into standard form - win_size    
    if is_auto_win == False:
        if len(win_size) == 1:
            if win_size%2 == 0:  # win size should be odd
                win_size += 1
            win_size = win_size * np.ones(no_chroms)
        else:
            for chrom in range(no_chroms):
                if win_size[chrom]%2 == 0:
                    win_size[chrom] += 1
    else:
        #determine win size here
        win_size = []
        for chrom in chroms:
            win_size.append(determine_shared_winsize(samples, chrom, **winsize_alg))
    # default reference winsize
    if ref_winsize == 0:
        ref_winsize = np.round(np.mean(win_size, 0))
        for chrom in range(no_chr):
            if not ref_winsize[chrom]%2: #win_size is even
                ref_winsize[chrom] += 1 #make it odd
    # bring parameters into standard form - ref
    if ref:
        ref.merge(False)
        ref.scale()
        
    # bring parameters into standard form - drate
    if len(drate) == 1:
        drate = drate * np.ones(no_samples)
        
    # bring parameters into standard form - lcf
    if len(lcf) == 1:
        lcf = lcf * np.ones(no_samples)
        
    # check consistency
    if match_histogram and not ref:
        print("Histogram matching requires a reference")
        sys.exit(1)

    # process min_finite
    if min_finite <= 1:
        min_finite = np.ceil(min_finite*no_samples)
        
    
    m = [[] for x in range(no_chroms)]
    dm = [[] for x in range(no_chroms)]
    
    # loop on chromosomes
    j = 0
    for chrom in chroms:
        # find chromosome index in each sample
        chr_idx = np.nan * np.ones(no_samples)
        for samp in range(no_samples):
            chr_idx[samp] = samples[samp].index([chrom])[0]
        chr_idx = chr_idx.astype(int)  # need to be int to use as indices
        # find number of CpG positions along the chromosome
        no_positions = len(samples[0].no_t[chr_idx[0]])
        # compute tij and nij
        tij = np.zeros((no_samples,no_positions))
        nij = np.zeros((no_samples,no_positions))
        for samp in range(no_samples):
            #[tij(samp,:), nij(samp,:)] = samples{samp}.smooth(c_idx(samp), p_alg.win_size(chr));
            [tij[samp], nij[samp]] = samples[samp].smooth(chr_idx[samp], int(win_size[samp]))
            # remove regions with particularly low coverage
            #lct = findlowcoveragethreshold(nij(samp,:),p_alg.lcf(samp));
            lct = find_low_coverage_thresh(nij[samp], lcf[samp])
            #nij(samp,nij(samp,:)<lct) = nan;
            nij[samp][nij[samp]<lct] = np.nan
        # estimate methylation
        [mi, dmi, m0] = ancient_Newton_Raphson(max_iterations, tol, drate, tij, nij)
    
        # substitute NaNs when there are too many NaNs in the original data
        finit = np.isfinite(tij) & np.isfinite(nij)
        finit = finit.sum(axis=0)  # this does not work for 1d array
        mi[finit<min_finite] = np.nan  # any chance this will work?
        #mi(finit < p_alg.min_finite) = nan;
        
        # modify methylation to match histogram
        idx_finite = np.where(np.isfinite(mi))
        if match_histogram:
            # hard-coded parameters
            ref_bins = 100
            sig_bins = 1000
            # use only finite elements
            sig = [mi[x] for x in idx_finite[0]]
            # x-axis for reference and signal
            ref_edges = np.linspace(0,1,ref_bins+1)
            sig_edges = np.linspace(0,max(sig),sig_bins+1)
            ref_binwidth = ref_edges[1] - ref_edges[0]
            # smooth the reference
            ref_idx = ref.index([chrom])[0]
            #vec = ref.smooth(ref_idx,win_size[j],name=chrom)[0]
            vec = ref.smooth(ref_idx,ref_winsize,name=chrom)[0]
            vec = [vec[x] for x in idx_finite[0]]
            # generate histograms
            hist = np.histogram(vec, bins=ref_edges)[0]
            N = np.sum(hist)
            c = np.cumsum(hist)
            N_ref = c/N
            hist = np.histogram(sig, bins=sig_edges)[0]
            N = np.sum(hist)
            c = np.cumsum(hist)
            N_sig = c/N
            # generate the mapping
            hmap = np.zeros(sig_bins)
            for i in range(sig_bins):
                # find closest value on the CDF
                imin = np.argmin(abs(N_ref - N_sig[i]))
                hmap[i] = ref_edges[imin] + 0.5*ref_binwidth
            # make the range precisely [0,1]
            hmap = (hmap - hmap[0]) / np.ptp(hmap)
            # apply the mapping
            sig1 = np.digitize(sig,sig_edges)
            sig1[np.where(sig1 == len(sig_edges))] = len(sig_edges) - 1
            tmp = [hmap[x-1] for x in sig1]
            d = dict(zip(idx_finite[0], tmp))
            for i in idx_finite[0]:
                mi[i] = d[i]
        else:
            #mi(idx_finite) = min(max(mi(idx_finite),0),1)
            for x in idx_finite[0]:
                mi[x] = min(max(mi[x],0),1)
        
        m[j] = mi
        dm[j] = dmi
        j += 1
    return (m, dm)

def plot_region(region, gc, samples, gene_bed=None, cgi_bed=None, widenby=0):
    panels = {}
    genes = pbt.BedTool(gene_bed)
    #genes_no_dups = genes.groupby(g=[1,2,3,6], c='4,5', o='distinct').cut([0,1,2,4,5,3], stream=False)
    genes_no_dups = genes.groupby(g=[1,2,3,6], c='4,5', o='distinct').cut([0,1,2,4,5,3])  # is stream nec?
    cgis = pbt.BedTool(cgi_bed)
    if gene_bed:
        panels["genes"] = {"type":"interval", "direction":True, "interval_names":True, "data":genes}
    if cgi_bed:
        panels["cgis"] = {"type":"interval", "direction":False, "interval_names":False, "data":cgis}
    no_panels = len(panels)
    region = standardize_region(region)  # nec?
    no_samples = len(samples)
    orig_start = region["start"]
    orig_end = region["end"]
    region["start"] = region["start"] - widenby
    region["end"] = region["end"] + widenby
    gc_chr = gc.chr_names.index(region["chrom"])
    gc_chr_coords = gc.coords[gc_chr]
    start = np.where(gc_chr_coords >= region["start"])[0][0]
    end = np.where(gc_chr_coords <= region["end"])[0][-1]
    orig_start = np.where(gc_chr_coords >= orig_start)[0][0]
    orig_end = np.where(gc_chr_coords <= orig_end)[0][-1]
    width = end - start + 1
    vals = np.zeros((no_samples, width))
    for sample in range(no_samples):
        samp_chr = samples[sample].index([region["chrom"]])[0]
        meth = samples[sample].methylation["methylation"][samp_chr]
        vals[sample] = meth[start:end+1]
    screen_height = get_monitors()[0].height
    # these should go in config
    OPTIMAL_HEIGHT_METH_LANE = 27
    OPTIMAL_SPACE = 20
    OPTIMAL_HEIGHT_INTERVAL_WITH_NAME = 15
    OPTIMAL_HEIGHT_INTERVAL_NO_NAME = 5
    OPTIMAL_HEIGHT_XSCALE = 15
    FIG_TOP_MARGIN = 100
    BOTTOM_MARGIN = 30
    MIN_SHRINKAGE = 0.3
    RIGHT_SHIFT = 0.01
    LEFT_SHIFT = 1
    BACK_COLOR = None
    
    # compute win_height
    meth_height = OPTIMAL_HEIGHT_METH_LANE * no_samples
    win_height = 3*OPTIMAL_SPACE + meth_height + OPTIMAL_HEIGHT_XSCALE + BOTTOM_MARGIN
    panel_order_intervals = []
    panel_height = {}
    for key in panels.keys():
        panel_order_intervals.append(panels[key])
        if panels[key]["interval_names"]:
            panel_height[key] = OPTIMAL_HEIGHT_INTERVAL_WITH_NAME
            if panels[key]["direction"]:
                panel_height[key] *= 2
        else:
            panel_height[key] = OPTIMAL_HEIGHT_INTERVAL_NO_NAME
    win_height += sum(panel_height.values()) + no_panels*OPTIMAL_SPACE
    shrinkage_factor = 1  # in config?
    if win_height > screen_height:
        shrinkage_factor = screen_height / win_height
        if shrinkage_factor > MIN_SHRINKAGE:
            win_height = screen_height
        else:
            shrinkage_factor = MIN_SHRINKAGE
            win_height = MIN_SHRINKAGE * win_height
    """
    # y_min is the bottom y-coordinate of the figure
    y_min = max(screen_height - win_height - FIG_TOP_MARGIN,1)
    # open figure
    fig = plt.subplot()
    bbox = fig.get_window_extent()
    pos = [bbox.xmin, y_min, bbox.width, win_height]
    fig.set_position(pos)
    pos_axes = pos[:]
    pos_axes[0] += pos[2]/10
    pos_axes[2] -= pos[2]/5
    ax.Axes.set_position(fig, pos_axes)  # store in var?
    """  # not sure if any of above relevant
    
    #fig, axes = plt.subplots()
    #fig, [ax3, ax2, ax1] = plt.subplots(3, sharex=True)
    fig, axes = plt.subplots(4, sharex=True)
    plt.ylim([0,1])
    plt.xlim([start, end])
    ax1 = axes[-1]
    pos = ax.Axes.get_position(ax1)
    #pos_axes = [pos.xmin, pos.ymin, 0.8, 0.03]  # currently hardcoded. not sure how to generalize
    pos_ax1 = [pos.xmin, pos.ymin, 0.75, 0.05]
    ax.Axes.set_position(ax1, pos_ax1)
    plt.setp(ax1, xlim=([start, end]), ylim=([0,0.5]))
    ax1.spines['top'].set_visible(False)
    ax1.spines['left'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.get_yaxis().set_visible(False)
    #ax.Axes.set_position(axes, pos_axes)
    #axes.spines['top'].set_visible(False)
    #axes.spines['left'].set_visible(False)
    #axes.spines['right'].set_visible(False)
    #axes.get_yaxis().set_visible(False)
    bar_scale = max(5*round(0.2*0.1*width), 5)
    x_start = start + np.ceil(bar_scale/2)
    plt.axhline(y=0.1, xmin=0.08, xmax=0.18, color="k")
    plt.text(x_start+0.05*bar_scale+4, 0.15, f"{bar_scale} CpGs")
    
    reg_no_chr = region.copy()
    reg_no_chr["chrom"] = reg_no_chr["chrom"].replace("chr", "")  # to match gene file, which has chr nums, not chr<num>
    i = -2
    for key in panels.keys():
        ax2 = axes[i]
        b = panels[key]["data"]
        if key == "genes":
            reg_bed = f'''
            {reg_no_chr["chrom"]} {reg_no_chr["start"]} {reg_no_chr["end"]}
            '''
        elif key == "cgis":
            reg_bed = f'''
            {region["chrom"]} {region["start"]} {region["end"]}
            '''
        a = pbt.BedTool(reg_bed, from_string=True)
        intersection = b.intersect(a, u=True)
        consolidated = intersection.merge(s=True, c="4,5,6", o="distinct")  # as yet untested. used d=-1? 
        #still need to merge intervals?
        interval_span = 1
        subset = consolidated.filter(lambda x: x.strand == "+")  # count lines on plus strand
        tot = len(subset)  # can't do len(subset) more than once (lambda?)
        if panels[key]["direction"] == True:
            if tot != len(consolidated) and tot != 0:  # not all on same strand
                interval_span = 2
        # sort intervals
        starts = [x.start for x in consolidated]
        ends = [x.end for x in consolidated]
        # plot intervals
        for i in range(len(starts)):
        #for st in starts:
            istart = np.where(gc_chr_coords >= starts[i])[0][0]
            truncated_left = False
            if istart < start:
                truncated_left = True
                istart = start
        #for en in ends:
            iend = np.where(gc_chr_coords <= ends[i])[0][-1]
            truncated_right = False
            if iend > end:
                truncated_right = True
                iend = end
            add = 0
            pos = ax.Axes.get_position(ax1)
            pos_ax2 = [pos.xmin, .18, pos.width, 0.1]
            ax.Axes.set_position(ax2, pos_ax2)
            plt.setp(ax2, xlim=([start, end]), ylim=([0,interval_span]))


            if interval_span == 2:
                if consolidated[i].strand == "+":
                    rect = p.Rectangle([istart, 0], iend-istart, 1, facecolor="black")
                    text_y = 0.82 
                else:
                    rect = p.Rectangle([istart, 1], iend-istart, 2, facecolor="black")
                    text_y = 0.32
                    add = 1
            else:
                rect = p.Rectangle([istart, 0], iend-istart, 1, facecolor="black")
                text_y = 1.05
            ax2.add_patch(rect)
            # add white arrows at the ends of truncated intervals
            interval_width = iend - istart + 1
            ar_width = min(width/30,interval_width/4)
            if truncated_left:
                ar_x = [start+2*ar_width, start+ar_width, start, start+ar_width, start+2*ar_width, start+ar_width]
                ar_y = [x + add for x in [0, 0, 0.5, 1, 1, 0.5]]
                coords = list(zip(ar_x, ar_y))
                arrow_left = p.Polygon(coords, facecolor = 'white')
                ax2.add_patch(arrow_left)
            if truncated_right:
                ar_x = [end-2*ar_width, end-ar_width, end, end-ar_width, end-2*ar_width, end-ar_width]
                ar_y = [x + add for x in [0, 0, 0.5, 1, 1, 0.5]]
                coords = list(zip(ar_x, ar_y))
                arrow_right = p.Polygon(coords, facecolor = 'white')
                ax2.add_patch(arrow_right)
            # add name into the interval
            if panels[key]["interval_names"] == True:
                # place in the middle of the interval
                name = consolidated[i].name
                #name = "very long gene name"
                t = plt.text(0.5*(istart+iend), add + text_y, f"{name}", color="w", ha="center")
                #plt.draw()
                plt.savefig("fname")
                # if too big, move to the side
                ext = t.get_window_extent()
                #pos = t.get_position()
                if ext.width > interval_width:
                    if iend-start > end-iend:
                    # place on the left side
                        t.set_position([istart, add+text_y])
                        t.set_ha("right")
                    else:
                        t.set_position([iend, add+text_y])
                        t.set_ha("left")
                    t.set_color("black")
            ax2.spines['top'].set_visible(False)
            ax2.spines['left'].set_visible(False)
            ax2.spines['right'].set_visible(False)
            ax2.spines['bottom'].set_visible(False)
            ax2.tick_params(bottom=False)
            ax2.get_yaxis().set_visible(False) 
            # write interval name
            ax2.text(end+2, 0.5*interval_span, key.capitalize(), ha='left', va='center')  
            # make an arrow designating interval direction
        if starts:
            if panels[key]["direction"] == True:
                if interval_span == 2:
                    ax2.annotate('', xy=(start-1, 0.5), xytext=(start-bar_scale, 0.5), arrowprops=dict(arrowstyle="->", color='k'), annotation_clip=False)
                    ax2.annotate('', xy=(start-bar_scale, 1.5), xytext=(start-1, 1.5), arrowprops=dict(arrowstyle="->", color='k'), annotation_clip=False)
                else:
                    if consolidated[i].strand == "+":
                        ax2.annotate('', xy=(start-1, 0.5), xytext=(start-bar_scale, 0.5), arrowprops=dict(arrowstyle="->", color='k'), annotation_clip=False)
                    else:
                        ax2.annotate('', xy=(start-bar_scale, 0.5), xytext=(start-1, 0.5), arrowprops=dict(arrowstyle="->", color='k'), annotation_clip=False)
                
        plt.show()
        i -= 1
    
    
    
    
    

        

if __name__ == "__main__":
    reg = "chr1:234,543,678-234,567,890"
    new_reg = standardize_region(reg)
    print(new_reg)

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
from matplotlib.colors import ListedColormap
from matplotlib.patches import Rectangle, Polygon
from matplotlib import patheffects as pe
from types import SimpleNamespace
import time

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
            take_me = np.nanmax(np.array([arr[i],arr[i+1]]))  # regular max function will always take nan
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

def normalize_chrom(chrom_str):
	"""Removes  'chr' prefix from chromosome names, if present
	
	Input:  chrom_str	chromosome name
	
	Output: chromomsome number, converted to int
	"""
	# Remove 'chr' prefix if present
	if chrom_str.startswith('chr'):
		chrom_str = chrom_str[3:]
	# Return numeric chromosomes only (skip X/Y/M/etc.)
	return int(chrom_str) if chrom_str.isdigit() else None

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
    """Solves for m when all samples are ancient (vectorized version).

    Input: max_iterations    maximum number of iterations
           min_tol           convergence tolerance
           pi                deamination rates of all samples (1D array-like, length = no_samples)
           tij               number of Ts in sample <i>, window <j>   (shape: no_samples x no_positions)
           nij               number of reads in sample <i>, window <j> (same shape)

    Output: m        methylation vector (1D, length = no_positions)
            dm       standard error of methylation (1D, length = no_positions)
            m0       initial guess (for debugging)
    """
    # Cast to arrays
    pi  = np.asarray(pi,  dtype=float)
    tij = np.asarray(tij, dtype=float)
    nij = np.asarray(nij, dtype=float)

    # Replace NaNs in tij/nij with 0 (same logic as your flatten/union approach)
    mask = np.isnan(tij) | np.isnan(nij)
    if np.any(mask):
        tij = tij.copy()
        nij = nij.copy()
        tij[mask] = 0.0
        nij[mask] = 0.0

    no_samples, no_pos = tij.shape

    # Precompute constants
    Tj   = np.sum(tij, axis=0)           # sum over samples
    Tpij = pi @ tij                      # pi dot tij  -> shape (no_pos,)
    Npij = pi @ nij                      # pi dot nij  -> shape (no_pos,)

    # Initial guess m0 (same formula as before)
    m = Tj / (Npij - Tpij)
    m0 = m.copy()

    # Precompute sample-dependent array A = nij - tij (used in every iteration)
    A = nij - tij                        # shape (no_samples, no_pos)
    pi_col = pi[:, None]                 # shape (no_samples, 1)

    # Newton–Raphson iterations
    iteration = 1
    tol = 1.0
    while iteration < max_iterations and tol > min_tol:
        m_prev = m

        # dℓ/dm = Tj / m_prev - sum_i ( (nij - tij)_ij * π_i / (1 - m_prev_j π_i) )
        # d²ℓ/dm² = -Tj / m_prev² - sum_i ( (nij - tij)_ij * π_i² / (1 - m_prev_j π_i)² )
        m_prev_row = m_prev[None, :]             # shape (1, no_pos)
        den = 1.0 - pi_col * m_prev_row          # shape (no_samples, no_pos)

        # Avoid noisy warnings; behaviour is controlled by nan_to_num later
        with np.errstate(divide="ignore", invalid="ignore"):
            term1 = A * pi_col / den             # contribution to dℓ/dm, shape (no_samples, no_pos)
            term2 = A * (pi_col ** 2) / (den ** 2)

        dldm   = Tj / m_prev - term1.sum(axis=0)
        d2ldm2 = -Tj / (m_prev ** 2) - term2.sum(axis=0)

        step = np.nan_to_num(dldm / d2ldm2)
        m = m_prev - step

        # Same convergence criterion as original: max |Δm / m_prev|
        with np.errstate(divide="ignore", invalid="ignore"):
            rel_change = np.abs((m - m_prev) / m_prev)
        tol = np.nanmax(rel_change)

        iteration += 1

    # Fisher information I and standard error dm = sqrt(1/I)
    m_row = m[None, :]
    den = 1.0 - pi_col * m_row
    with np.errstate(divide="ignore", invalid="ignore"):
        termI = A * (pi_col ** 2) / (den ** 2)

    I = Tj / (m ** 2) + termI.sum(axis=0)
    dm = np.sqrt(1.0 / I)

    return (m, dm, m0)


def pooled_methylation(samples, chroms, win_size="auto", winsize_alg={}, lcf=0.05,
                       drate=[], min_finite=1, match_histogram=True, ref=None,
                       ref_winsize=0, max_iterations=20, tol=0.001):
    """Computes pooled ancient methylation.

    Input:
        samples         list of amSamples. Methylation is computed based on pooling
                        information from all these samples.
        chroms          names of chromosomes (list)
        win_size        window size for smoothing. If 'auto', a recommended value is
                        computed for each chromosome. Otherwise, it can be a scalar
                        (used for all chromosomes) or a vector over samples.
        winsize_alg     dict with parameters required to determine window size
        lcf             low coverage factor
        drate           deamination rates of each ancient sample. If provided, its
                        length should match the number of ancient samples.
        min_finite      minimum number (or fraction) of ancient samples for which
                        we require data in a position; otherwise the position is NaN.
        match_histogram whether to match the histogram to that of a provided reference
        ref             mmSample representing a reference methylome (for matching)
        ref_winsize     window size used for the reference histogram
        max_iterations  maximum number of iterations in Newton–Raphson
        tol             convergence tolerance

    Output:
        m   list (per chromosome) of pooled methylation vectors
        dm  list (per chromosome) of standard errors
    """
    no_samples = len(samples)

    # if chroms not provided, use all chroms from first sample
    no_chroms = len(chroms)
    if no_chroms == 0:
        chroms = samples[0].chr_names
        no_chroms = samples[0].no_chrs

    # bring drate into standard form
    if len(drate) == 0:
        drate = np.nan * np.ones(no_samples)
        for samp in range(no_samples):
            drate[samp] = samples[samp].d_rate["rate"]["global"]
    elif len(drate) == 1:
        drate = drate * np.ones(no_samples)

    # detect auto-win mode
    is_auto_win = isinstance(win_size, str) and (win_size == "auto")

    # bring win_size into standard form
    if not is_auto_win:
        win_size = np.array(win_size, copy=False)
        if win_size.ndim == 0:  # scalar
            w = int(win_size)
            if w % 2 == 0:
                w += 1
            win_size = np.full(no_samples, w, dtype=int)
        else:
            win_size = win_size.astype(float)
            for i in range(win_size.shape[0]):
                if win_size[i] % 2 == 0:
                    win_size[i] += 1
            win_size = win_size.astype(int)
    else:
        # determine per-sample win size from winsize_alg
        tmp = []
        for chrom in chroms:
            tmp.append(determine_shared_winsize(samples, chrom, **winsize_alg))
        # here we use the same window size for all samples on a given chrom
        # but groupDMRs passes per-sample win_size, so in practice you almost
        # never hit this branch in your current pipeline
        win_size = np.array(tmp, dtype=int)

    # default reference winsize (used only when match_histogram is True)
    if ref_winsize == 0:
        # take mean over whatever win_size you gave us and ensure it's odd
        ref_winsize = int(np.round(np.nanmean(win_size)))
        if ref_winsize % 2 == 0:
            ref_winsize += 1

    # bring parameters into standard form – lcf
    if np.isscalar(lcf):
        lcf = lcf * np.ones(no_samples)

    # check consistency of histogram matching
    if match_histogram and not ref:
        print("Histogram matching requires a reference")
        sys.exit(1)

    # prepare reference once (cached on the ref object)
    if ref:
        if not getattr(ref, "_hist_prepared", False):
            ref.merge(False)
            ref.scale()
            setattr(ref, "_hist_prepared", True)

    # process min_finite
    if min_finite <= 1:
        min_finite = int(np.ceil(min_finite * no_samples))

    m = [[] for _ in range(no_chroms)]
    dm = [[] for _ in range(no_chroms)]

    # loop on chromosomes
    j = 0
    for chrom in chroms:
        # find chromosome index in each sample
        chr_idx = np.nan * np.ones(no_samples)
        for samp in range(no_samples):
            chr_idx[samp] = samples[samp].index([chrom])[0]
        chr_idx = chr_idx.astype(int)

        # number of CpG positions along this chromosome
        no_positions = len(samples[0].no_t[chr_idx[0]])

        # compute tij and nij
        tij = np.zeros((no_samples, no_positions))
        nij = np.zeros((no_samples, no_positions))
        for samp in range(no_samples):
            # smooth counts
            [tij[samp], nij[samp]] = samples[samp].smooth(
                chr_idx[samp],
                int(win_size[samp])
            )
            # low coverage threshold
            lct = find_low_coverage_thresh(nij[samp], lcf[samp])
            nij[samp][nij[samp] < lct] = np.nan

        # Newton–Raphson to estimate methylation
        mi, dmi, m0 = ancient_Newton_Raphson(max_iterations, tol, drate, tij, nij)

        # substitute NaNs where there are too few informative samples
        finit = np.isfinite(tij) & np.isfinite(nij)
        finit = finit.sum(axis=0)
        mi[finit < min_finite] = np.nan

        # histogram matching (optional)
        idx_finite = np.where(np.isfinite(mi))
        if match_histogram:
            # hard-coded parameters
            ref_bins = 100
            sig_bins = 1000

            # use only finite elements
            sig = [mi[x] for x in idx_finite[0]]

            # x-axis for reference and signal
            ref_edges = np.linspace(0, 1, ref_bins + 1)
            sig_edges = np.linspace(0, max(sig), sig_bins + 1)
            ref_binwidth = ref_edges[1] - ref_edges[0]

            # smooth the reference
            print("smoothing the reference with window:", ref_winsize)
            ref_idx = ref.index([chrom])[0]
            vec = ref.smooth(ref_idx, ref_winsize, name=chrom)[0]
            vec = [vec[x] for x in idx_finite[0]]

            # reference histogram
            hist = np.histogram(vec, bins=ref_edges)[0]
            N = np.sum(hist)
            c = np.cumsum(hist)
            N_ref = c / N

            # signal histogram
            hist = np.histogram(sig, bins=sig_edges)[0]
            N = np.sum(hist)
            c = np.cumsum(hist)
            N_sig = c / N

            # generate mapping
            hmap = np.zeros(sig_bins)
            for i in range(sig_bins):
                imin = np.argmin(np.abs(N_ref - N_sig[i]))
                hmap[i] = ref_edges[imin] + 0.5 * ref_binwidth

            # make the range precisely [0,1]
            hmap = (hmap - hmap[0]) / np.ptp(hmap)

            # apply the mapping
            sig1 = np.digitize(sig, sig_edges)
            sig1[np.where(sig1 == len(sig_edges))] = len(sig_edges) - 1
            tmp = [hmap[x - 1] for x in sig1]
            d = dict(zip(idx_finite[0], tmp))
            for i in idx_finite[0]:
                mi[i] = d[i]
        else:
            # clamp finite mi to [0,1]
            for x in idx_finite[0]:
                mi[x] = min(max(mi[x], 0), 1)

        m[j] = mi
        dm[j] = dmi
        j += 1

    return (m, dm)

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







# ---------- small helpers for plot region ----------
def _comma(n: int) -> str: return f"{n:,}"

def _cmap(color_blind=False, mode="20and80are0") -> ListedColormap:
    if not color_blind:
        r_ext, g_ext, b_ext = np.array([72,220])/255, np.array([183,59])/255, np.array([48,16])/255
    else:
        r_ext, g_ext, b_ext = np.array([25,220])/255, np.array([0,59])/255,  np.array([51,16])/255
    if mode == "20and80are0":
        r = np.concatenate([np.full(20, r_ext[0]), np.linspace(r_ext[0], r_ext[1], 31), np.full(50, r_ext[1])])
        g = np.concatenate([np.full(50, g_ext[0]), np.linspace(g_ext[0], g_ext[1], 31), np.full(20, g_ext[1])])
        b = np.linspace(b_ext[0], b_ext[1], 101)
    else:
        r = np.linspace(r_ext[0], r_ext[1], 101); g = np.linspace(g_ext[0], g_ext[1], 101); b = np.linspace(b_ext[0], b_ext[1], 101)
    return ListedColormap(np.stack([r,g,b], axis=1), name="methyl_20_80")

def _first_geq_idx(vec, bp): return int(np.searchsorted(vec, bp, side="left"))
def _last_leq_idx(vec, bp):  return int(np.searchsorted(vec, bp, side="right") - 1)
def _region_bed(chrom, start, end): return pbt.BedTool(f"{chrom}\t{start}\t{end}\n", from_string=True)
def _normalize_groups(groups, no_samples, samples=None):
    """
    Returns a list[str] of length no_samples with each sample's group name.
    Supports:
      - your dict: {'group_nums': [...], 'group_names': [...], 'no_groups': N, 'positions': [[...],[...],...]}
      - list[str] already aligned to samples
      - dict[str, list[int]] mapping name -> indices
      - None -> infer from sample.name containing 'farmer'/'hg'
    """
    # your dict format
    if isinstance(groups, dict) and ("positions" in groups or "group_nums" in groups):
        names = groups.get("group_names", [])
        if "positions" in groups:
            labels = [""] * no_samples
            for gi, idxs in enumerate(groups["positions"]):
                name = names[gi] if gi < len(names) else f"group{gi+1}"
                for idx in idxs:
                    if 0 <= idx < no_samples:
                        labels[idx] = name
            return labels
        else:  # group_nums
            nums = groups["group_nums"]
            # preserve first-appearance order
            uniq = []
            for n in nums:
                if n not in uniq:
                    uniq.append(n)
            name_map = {n: (names[i] if i < len(names) else f"group{n}") for i, n in enumerate(uniq)}
            return [name_map.get(n, f"group{n}") for n in nums]

    # dict[name] -> [indices]
    if isinstance(groups, dict):
        labels = [""] * no_samples
        for name, idxs in groups.items():
            for idx in idxs:
                if 0 <= idx < no_samples:
                    labels[idx] = name
        return labels

    # list[str] already aligned
    if isinstance(groups, (list, tuple)) and len(groups) == no_samples:
        return list(groups)

    # infer from sample names
    out = []
    for s in (samples or []):
        nm = getattr(s, "name", "").lower()
        if "farmer" in nm: out.append("Farmers")
        elif "hg" in nm or "hunter" in nm: out.append("HGs")
        else: out.append("Samples")
    return out

# ---------- plot region of DMR ----------
def plot_region(region, gc, samples, gene_bed=None, cgi_bed=None, widenby=0, *,
                groups=None, savepath=None, show=True, close=False):
    """
    groups: can be
      - list[str] of length == len(samples), e.g. ["farmer","farmer","HG",...]
      - dict[str, list[int]] mapping group -> sample indices (0-based)
      - None: will try to infer from sample.name (contains 'farmer' or 'hg')
    """
    chrom = region["chrom"]
    orig_start, orig_end = int(region["start"]), int(region["end"])
    start_bp, end_bp = orig_start - widenby, orig_end + widenby

    # CpG indices
    chr_idx = gc.chr_names.index(chrom)
    gc_chr_coords = gc.coords[chr_idx]
    beg_i = _first_geq_idx(gc_chr_coords, start_bp)
    fin_i = _last_leq_idx(gc_chr_coords, end_bp)
    width = int(fin_i - beg_i + 1)
    
    # DMR (original, pre-widen) indices inside the plotted window
    window = gc_chr_coords[beg_i:fin_i+1]
    beg0_local = int(np.searchsorted(window, orig_start, side='left'))        # first >=
    end0_local = int(np.searchsorted(window, orig_end,   side='right') - 1)   # last <=
    beg0 = beg_i + beg0_local
    end0 = beg_i + end0_local
    dmr_width = int(end0 - beg0 + 1)


    # methylation matrix
    no_samples = len(samples)
    vals = np.zeros((no_samples, width), dtype=float)
    for s_i, s in enumerate(samples):
        samp_chr = s.index([chrom])[0]
        store = s.methylation
        meth = store["methylation"][samp_chr] if isinstance(store, dict) else store[samp_chr]
        vals[s_i, :] = meth[beg_i:fin_i+1]

    have_genes = gene_bed is not None
    have_cgis  = cgi_bed  is not None
    

    # ==== layout: heatmap (top), then CGIs, Genes, then x-axis ====
    height = [2.0 + 0.25*no_samples]    # heatmap
    if have_cgis:  height.append(0.40)
    if have_genes: height.append(1.6)
    height.append(0.8)                   # x-axis
    fig = plt.figure(figsize=(11, 6), constrained_layout=True)
    gs = fig.add_gridspec(len(height), 1, height_ratios=height)
    row = 0

    # ---- HEATMAP ----
    axh = fig.add_subplot(gs[row, 0]); row += 1
    axh.imshow(
        vals, aspect="auto", interpolation="nearest", origin="lower",
        extent=[beg_i, fin_i+1, 1, no_samples+1],
        cmap=_cmap(), vmin=0, vmax=1
    )
    axh.set_xlim([beg_i, fin_i+1])
    axh.set_ylim([1, no_samples+1])
    axh.set_yticks(np.arange(1.5, no_samples+0.6, 1))
    axh.set_yticklabels([s.name for s in samples])
    axh.tick_params(axis="x", bottom=False, labelbottom=False)

    if widenby > 0:
        axh.axvline(beg0, ls=":", color="k")
        axh.axvline(end0, ls=":", color="k")

    # ---- GROUP SEPARATORS + LABELS ----
    sample_groups = _normalize_groups(groups, no_samples, samples)
	
	# Build contiguous blocks in the current sample order
    order, seen = [], {}
    for i, g in enumerate(sample_groups):
        if g not in seen:
            seen[g] = []
            order.append(g)
        seen[g].append(i)
    group_blocks = [(g, seen[g]) for g in order]
	
	# Draw separator lines and right-edge labels
    y_base = 1
    for gname, idxs in group_blocks:
        n = len(idxs)
        cy = y_base + n/2.0
        y_norm = (cy - 1) / no_samples          # map row center to [0..1] in axis coords
        axh.text(1.02, y_norm, gname, transform=axh.transAxes,
         va="center", ha="left", clip_on=False)   # 1.02 pushes it a bit outside

        y_base += n
        if y_base <= no_samples:
            axh.axhline(y_base, color="k", lw=2.2)

    # ---- interval-track helper ----
    def draw_interval_track(ax, bedtool, name, directional, bar_frac=1.0):
    # intersect with region
        rbed = _region_bed(chrom, start_bp, end_bp)
        inter = (bedtool if isinstance(bedtool, pbt.BedTool) else pbt.BedTool(bedtool))
        inter = inter.intersect(rbed, wa=True).sort()
        if len(inter) == 0:
            ax.axis("off"); return
    # one helper for both sides; scales to bar height
        def arrow(side, y0, h):
            arw = max(1, width/30)
            if side == "L":
                xs = [beg_i+2*arw, beg_i+arw, beg_i, beg_i+arw, beg_i+2*arw, beg_i+arw]
            else:
                xs = [fin_i-2*arw, fin_i-arw, fin_i, fin_i-arw, fin_i-2*arw, fin_i-arw]
            ys = (np.array([0,0,.5,1,1,.5]) * h) + y0
            ax.add_patch(Polygon(np.c_[xs, ys], color="w"))

        # --- helpers for genes ---
        def _bp_to_ix(bp):
            xs = np.arange(beg_i, fin_i + 1, dtype=float)
            ys = gc_chr_coords[beg_i:fin_i + 1].astype(float)
            return float(np.clip(np.interp(float(bp), ys, xs, left=xs[0], right=xs[-1]), xs[0], xs[-1]))

        def _pack_lanes(items):
            lanes, ends = [], []
            for f in sorted(items, key=lambda t: (t.start, t.end)):
                placed = False
                for li, e in enumerate(ends):
                    if f.start > e:
                        ends[li] = f.end
                        lanes[li].append(f)
                        placed = True
                        break
                if not placed:
                    lanes.append([f]); ends.append(f.end)
            return lanes

        if directional:
            # collapse transcripts to (name,strand) spans
            has_plus = has_minus = False
            groups = {}
            for f in inter:
                s  = getattr(f, "strand", ".")
                nm = getattr(f, "name", "")
                has_plus  |= (s == "+")
                has_minus |= (s == "-")
                key = (nm, s)
                gbp, ebp = int(f.start), int(f.end)
                if key in groups:
                    a, b = groups[key]
                    groups[key] = (min(a, gbp), max(b, ebp))
                else:
                    groups[key] = (gbp, ebp)

            feats = [SimpleNamespace(start=a, end=b, name=nm, strand=s)
                     for (nm, s), (a, b) in groups.items()]

            plus_lanes  = _pack_lanes([f for f in feats if f.strand == "+"])
            minus_lanes = _pack_lanes([f for f in feats if f.strand == "-"])

            ax.set_xlim([beg_i, fin_i + 1]); ax.set_ylim([0, 2]); ax.axis("off")

            min_w_for_label = max(5, int(width * 0.015))
            min_gap         = max(4, int(width * 0.008))

            def _draw_strand_lanes(lanes, base_y):
                if not lanes:
                    return
                vpad, gap, H = 0.12, 0.08, 1.0
                n = len(lanes)
                lane_h = max(0.02, (H - 2*vpad - (n - 1)*gap) / n)

                for li, lane in enumerate(lanes):
                    y0 = base_y + vpad + li * (lane_h + gap)
                    last_right = -np.inf
                    for feat in lane:
                        gbp, ebp = int(feat.start), int(feat.end)
                        ov_a = max(gbp, start_bp); ov_b = min(ebp, end_bp)
                        has_cols = bool(np.any((gc_chr_coords >= ov_a) & (gc_chr_coords <= ov_b)))
                        label = feat.name
                        if has_cols:
                            ib = _first_geq_idx(gc_chr_coords, ov_a)
                            ie = _last_leq_idx(gc_chr_coords,  ov_b)
                            ibc, iec = max(ib, beg_i), min(ie, fin_i)
                            if ibc > fin_i or iec < beg_i:
                                continue
                            w = iec - ibc + 1

                            ax.add_patch(Rectangle((ibc, y0), w, lane_h, color="k"))
                            if gbp < start_bp: arrow("L", y0, lane_h)
                            if ebp > end_bp:   arrow("R", y0, lane_h)

                            if label:
                                ok_in_bar = (w >= min_w_for_label) and (ibc > last_right + min_gap)
                                if ok_in_bar:
                                    txt = ax.text(ibc + w/2, y0 + 0.5*lane_h, label,
                                                  color="w", ha="center", va="center")
                                    txt.set_path_effects([pe.withStroke(linewidth=3, foreground="k")])
                                    last_right = ibc + w
                                else:
                                # place near the bar (prefer above; if no room, place below), with a white halo
                                    pad_y = 0.08
                                    top_of_block = base_y + 1.0
                                    y_text = y0 + lane_h + pad_y
                                    if y_text > top_of_block - 0.02:      # too close to top of this strand block → put below
                                        y_text = y0 - pad_y

                                    # center horizontally on the bar, clamped to the panel
                                    x_text = float(np.clip(ibc + w/2, beg_i + 2, fin_i - 2))

                                    txt = ax.text(x_text, y_text, label,
                                        ha="center", va="center", color="k", zorder=20, clip_on=False)
                                    txt.set_path_effects([pe.withStroke(linewidth=3, foreground="w")])
                        else:
                            # name-only glyph at bp midpoint with strand arrow
                            x_center = _bp_to_ix(0.5 * (gbp + ebp))
                            txt = ax.text(x_center, y0 + 0.5*lane_h,            # center on the lane
                                label, ha="center", va="center", color="k", zorder=20, clip_on=False)
                            txt.set_path_effects([pe.withStroke(linewidth=3, foreground="w")])

            _draw_strand_lanes(plus_lanes, 0.0)
            _draw_strand_lanes(minus_lanes, 1.0)

        else:
            # CGIs: same as before
            merged = inter.merge()
            ax.set_xlim([beg_i, fin_i + 1]); ax.set_ylim([0, 1]); ax.axis("off")
            for f in merged:
                gbp, ebp = int(f.start), int(f.end)
                ov_a = max(gbp, start_bp); ov_b = min(ebp, end_bp)
                if ov_b < ov_a: continue
                ib = _first_geq_idx(gc_chr_coords, ov_a)
                ie = _last_leq_idx(gc_chr_coords,  ov_b)
                ibc, iec = max(ib, beg_i), min(ie, fin_i)
                if ibc <= iec:
                    w = iec - ibc + 1
                    y0 = 0.5 - float(bar_frac)/2.0
                    ax.add_patch(Rectangle((ibc, y0), w, float(bar_frac), color="k"))
                    if gbp < start_bp: arrow("L", y0, float(bar_frac))
                    if ebp > end_bp:   arrow("R", y0, float(bar_frac))

        # right-edge track tag – vertically center on whichever lanes exist (axes coords)
        if directional:
            has_plus  = len(plus_lanes)  > 0
            has_minus = len(minus_lanes) > 0
            if has_plus and has_minus:
                y_tag = 0.50   # middle of the two-lane block
            elif has_plus:
                y_tag = 0.25   # center of the lower (plus) lane block
            elif has_minus:
                y_tag = 0.75   # center of the upper (minus) lane block
            else:
                y_tag = 0.50
        else:
            y_tag = 0.50

        ax.text(1.02, y_tag, name, transform=ax.transAxes,
            ha="left", va="center", clip_on=False, zorder=10)

        # left-margin arrows (match MATLAB behavior)
        if directional:
            pad = max(2, width/25)
            if has_plus and has_minus:
                ax.text(beg_i - pad, 0.5, r'$\rightarrow$', ha="right", va="center", weight="bold", clip_on=False, zorder=20)
                ax.text(beg_i - pad, 1.5, r'$\leftarrow$', ha="right", va="center", weight="bold", clip_on=False, zorder=20)
            elif has_plus:
                ax.text(beg_i - pad, 0.5, r'$\rightarrow$', ha="right", va="center", weight="bold", clip_on=False, zorder=20)
            elif has_minus:
                ax.text(beg_i - pad, 1.5, r'$\leftarrow$', ha="right", va="center", weight="bold", clip_on=False, zorder=20)

    # ---- CGIs ----
    if cgi_bed is not None:
        ax = fig.add_subplot(gs[row, 0]); row += 1
        draw_interval_track(ax, cgi_bed, "CGIs", directional=False, bar_frac=0.25)

    # ---- Genes ----
    if gene_bed is not None:
        ax = fig.add_subplot(gs[row, 0]); row += 1
        draw_interval_track(ax, gene_bed, "Genes", directional=True)
        
    # ---- x-axis ----
    axx = fig.add_subplot(gs[row, 0])
    axx.set_xlim([beg_i, fin_i+1]); axx.set_ylim([0, 1]); axx.set_yticks([])
    for side in ("top", "right", "left", "bottom"):
        axx.spines[side].set_visible(False)
    # show only the first and last genomic coordinates (bp)
    axx.set_xticks([beg_i, fin_i])
    axx.set_xticklabels([
        _comma(int(gc_chr_coords[beg_i])),
        _comma(int(gc_chr_coords[fin_i]))
    ], fontsize=9)
    bar_scale = max(int(5 * round(0.2 * 0.1 * width)), 5)
    x0 = beg_i + int(np.ceil(bar_scale / 2))
    axx.plot([x0, x0 + bar_scale], [0.5, 0.5], color="k", lw=1.5)
    axx.text(x0 + 0.5 * bar_scale, 0.78, f"{bar_scale} CpGs", ha="center")

    title_str = f"{chrom}:{_comma(orig_start)}-{_comma(orig_end)} ({_comma(dmr_width)} CpG positions)"
    axh.set_title(title_str, loc="center", pad=6, fontweight="bold")

    if savepath: fig.savefig(savepath, dpi=300, bbox_inches="tight")
    if show: plt.show()
    if close: plt.close(fig)
    

        

if __name__ == "__main__":
    reg = "chr1:234,543,678-234,567,890"
    new_reg = standardize_region(reg)
    print(new_reg)

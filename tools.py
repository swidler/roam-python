#!/usr/bin/python3

import numpy as np
import pickle
import scipy.stats as stats
import re
import sys
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
    zero_for_nan = [1 if x>0 else 0 for x in zero_for_nan]
    res = do_conv(vec, nans, tpl, mode, zero_for_nan)
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
            coverage = np.nan*np.ones(no_samples)
            for samp in range(no_samples):
                chr_num = samples[samp].index([chrom])[0]
                coverage[samp] = samples[samp].diagnostics["effective_coverage"][chr_num]
        #else:
        if len(drate) != no_samples:
            print(f"inserted drate vector (length {len(drate)}) does not match the number of samples ({no_samples}).")
            print("Ignoring user input.")
            drate = np.nan*np.ones(no_samples)
            for samp in range(no_samples):
                drate[samp] = samples[samp].d_rate["rate"]["global"]
        if min_meth>1:
            min_meth = 0.01 * min_meth
        
        #compute the window size
        cvg = sum(coverage)  # np.sum?
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
    Tj = np.nansum(tij,0)
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

def pooled_methylation(samples, chroms, win_size="auto", winsize_alg={}, lcf=0.05, drate=[], min_finite=1, match_histogram=True, ref=None, max_iterations=20, tol=0.001):
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
            [tij[samp], nij[samp]] = samples[samp].smooth(chr_idx[samp], int(win_size[j]))
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
            vec = ref.smooth(ref_idx,win_size[j],name=chrom)[0]
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

if __name__ == "__main__":
    reg = "chr1:234,543,678-234,567,890"
    new_reg = standardize_region(reg)
    print(new_reg)

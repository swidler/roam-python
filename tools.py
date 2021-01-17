#!/usr/bin/python3

import numpy as np
import pickle
import scipy.stats as stats
"""This module contains helper functions used in the RoAM process.
"""
def nanmerge(arr, operation):
    """Merges the two adjacent positions of each CpG. The rules of merger are as follows:
        (1) nan + nan = nan
        (2) value + nan = value
        (3) nan + value = value
        (4) value + value = factor*(value+value)
        The variable {factor} is 0.5 if {operation} is 'average', or 1 if {operation} is 'sum'.
    
    Input: arr        list with CpG positions for merging
           operation  either "Average" or "sum
    Output: merged list
    """
    result = []
    if operation == "average":
        factor = 0.5
    elif operation == "sum":
        factor = 1
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


if __name__ == "__main__":
    reg = "chr1:234,543,678-234,567,890"
    new_reg = standardize_region(reg)
    print(new_reg)

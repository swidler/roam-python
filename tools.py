#!/usr/bin/python3

import numpy as np
import pickle
import scipy.stats as stats

def nanmerge(arr, operation):
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
    import re
    if isinstance(region, dict):
        return region
    fields = re.split(" |:|-", region)
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
    outfile = open(filename, "wb") #open file for writing binaries
    pickle.dump(obj, outfile, -1) #serialize object using highest protocol
    outfile.close()

def load_object(filename):
    infile = open(filename, "rb") #open file for reading binaries
    obj = pickle.load(infile)
    infile.close()
    return obj

def nansmooth(vec, winsize, mode):
    tpl = np.ones(winsize)
    (nans, zero_for_nan) = get_zeros(vec, tpl, mode)
    res = do_conv(vec, nans, tpl, mode, zero_for_nan)
    return(res, zero_for_nan)

def nanconv(vec, tpl, mode):
    (nans, zero_for_nan) = get_zeros(vec, tpl, mode)
    zero_for_nan = [1 if x>0 else 0 for x in zero_for_nan]
    res = do_conv(vec, nans, tpl, mode, zero_for_nan)
    return res

def get_zeros(vec, tpl, mode):
    vec = np.array(vec)
    vec_len = len(vec)
    mod_vec = np.ones(vec_len)
    nans = np.isnan(vec)
    mod_vec[nans] = 0
    zero_for_nan = np.convolve(mod_vec, tpl, mode)
    return (nans, zero_for_nan)

def do_conv(vec, nans, tpl, mode, zero_for_nan):
    vec = np.array(vec)
    vec[nans] = 0
    res = np.convolve(vec, tpl, mode)/zero_for_nan #raises divide by zero runtime warnings
    return res

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
    no_dists = len(p)
    B = np.zeros((no_dists, N+1)) 
    vals = range(N+1)
    for x in range(no_dists):
        B[x,] = stats.binom.pmf(vals,N,p[x])
    return B

def log_likelihood(w,B,H):
    mat_prod = np.matmul((B+np.spacing(1)).T,w)
    llike = np.matmul(H[0],np.log(mat_prod)) 
    return llike

def bmm(H, p, w, tolerance, p_known):
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
    return p, w, llike

if __name__ == "__main__":
    reg = "chr1:234,543,678-234,567,890"
    new_reg = standardize_region(reg)
    print(new_reg)

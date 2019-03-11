#!/usr/bin/python3

import numpy as np
import pickle
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
    (nans, zero_for_nan) = get_zeroes(vec, tpl, mode)
    res = do_conv(vec, nans, tpl, mode, zero_for_nan)
    return(res, zero_for_nan)

def nanconv(vec, tpl, mode):
    (nans, zero_for_nan) = get_zeroes(vec, tpl, mode)
    zero_for_nan = [1 if x>0 else 0 for x in zero_for_nan]
    res = do_conv(vec, nans, tpl, mode, zero_for_nan)
    return res

def get_zeroes(vec, tpl, mode):
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
    res = np.convolve(vec, tpl, mode)/zero_for_nan
    return res


if __name__ == "__main__":
    reg = "chr1:234,543,678-234,567,890"
    new_reg = standardize_region(reg)
    print(new_reg)

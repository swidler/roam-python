#!/usr/bin/python3

import numpy
def nanconv(arr, operation):
    result = []
    if operation == "average":
        factor = 0.5
    elif operation == "sum":
        factor = 1
    for i in range(0,len(arr),2):
        if numpy.isnan(arr[i]):
            if numpy.isnan(arr[i+1]):
                result.append(numpy.nan) #NaN + NaN = NaN
            else:
                result.append(arr[i+1]) #NaN + <num> = <num>
        else:
            if numpy.isnan(arr[i+1]):
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
    std_region = {
            "chrom":chrom,
            "start":start,
            "end":end
            }
    return std_region
if __name__ == "__main__":
    reg = "chr1:234,543,678-234,567,890"
    new_reg = standardize_region(reg)
    print(new_reg)

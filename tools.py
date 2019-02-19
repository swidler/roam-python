#!/usr/bin/python3

def nanconv(arr, operation):
    result = []
    if operation == "average":
        factor = 0.5
    elif operation == "sum":
        factor = 1
    for i in range(0,len(arr),2):
        if arr[i] == "NaN":
            if arr[i+1] == "NaN":
                result.append("NaN") #NaN + NaN = NaN
            else:
                result.append(arr[i+1]) #NaN + <num> = <num>
        else:
            if arr[i+1] == "NaN":
                result.append(arr[i]) #<num> + NaN = <num>
            else:
                merged = (arr[i]+arr[i+1])*factor #num + num = (num + num) * factor. yields float result, with
                result.append(merged)             #sometimes strange precision (eg 0.8200000000000001)
    return result

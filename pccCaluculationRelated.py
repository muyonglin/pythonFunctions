###################################################################
#
#calculate pcc and other statistics
#
###################################################################

#calculate the median of a list
def medianManual(list1):
    mylist = sorted(list1)
    length = len(mylist)
    if len(mylist) < 1:
        print("empty list")
        return None
    else:
        if length % 2 == 1:
            return mylist[int(float(length)/2)]
        else:
            return float(mylist[int(float(length)/2)] + mylist[int(float(length)/2)-1])/2
#print medianManual([1,2,3,4])
#print medianManual([1,2,3])

#calculate median using numpy
def medianNumpy(list1):
    import numpy as np
    return np.median(np.array(list1))

#use pearsonr in scipy calculate pcc
def pcc(list1,list2):
    from scipy.stats import pearsonr
    (pcc,pvalue) = pearsonr(list1,list2)
    return (pcc,pvalue)



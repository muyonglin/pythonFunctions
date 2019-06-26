#!/usr/bin/env python
# encoding: utf-8

#read file into pandas dataframe
#polishing columns into the same length
def readFileIntoDataframe(fileName,sep="\t"):
    import pandas as pd
    termList = map(lambda x : x.split("\n")[0].split(sep),open(fileName,"r").readlines())
    termList = list(termList)
    termDataFrame = pd.DataFrame(termList)
    #print(termDataFrame.iloc[0,:])
    return termDataFrame
def combineFpkmsFromCufflinkFiles(fileNameList,sep="\t",fpkmColumn = [0,9]):
    import pandas as pd
    combinedDataFrame = pd.DataFrame()
    for fileName in fileNameList:
        tmpDataFrame = readFileIntoDataframe(fileName,sep)
        combinedDataFrame = pd.concat([combinedDataFrame,tmpDataFrame[fpkmColumn]],axis=1)
    return combinedDataFrame
def getFoldChange(DataFrameIn,fpkmColumnList,startRow=1):
    import pandas as pd
    import numpy as np
    #fpkmColumnList should be [[[ctrl1],[treat1]],[[ctrl2],[treat2]]]
    #when set index of DataFrameIn (in cufflinks fpkm results)
    #to gene name like "WBGene00197333"
    #will raise error like "ValueError: Shape of passed values is (1, 46277), indices imply (1, 46276)"
    #that is caused by duplicated indices
    #don't give index to DataFrameIn (the index will be row numbers)
    #will solve the problem
    foldChange = pd.DataFrame()
    for item in fpkmColumnList:
        for item2 in item:
            for item3 in item2:
                #log2 transformation
                DataFrameIn.iloc[startRow:,item3] = (DataFrameIn.iloc[startRow:,item3].astype(float)+1).apply(lambda x: np.log2(x))
        foldChange = pd.concat([foldChange,DataFrameIn.iloc[startRow:,item[1]].apply(lambda x: np.mean(x),axis=1) - DataFrameIn.iloc[startRow:,item[0]].apply(lambda x: np.mean(x),axis=1)],axis=1)
    return foldChange
def getFoldChangeFromLog2Dataframe(DataFrameLog2In,treatmentList):
    import pandas as pd
    import numpy as np
    #treatmentList should be [[[ctrl1],[treat1]],[[ctrl2],[treat2]]]
    #DataFrameLog2In only contains fpkmValue,
    #column names are sample names
    #indices are unrepeated numbers
    #when set index of DataFrameIn (in cufflinks fpkm results)
    #to gene name like "WBGene00197333"
    #will raise error like "ValueError: Shape of passed values is (1, 46277), indices imply (1, 46276)"
    #that is caused by duplicated indices
    #don't give index to DataFrameIn (the index will be row numbers)
    #will solve the problem
    foldChange = pd.DataFrame()
    fpkmMean = pd.DataFrame()
    for item in treatmentList:
        for item2 in item:
            fpkmMean = pd.concat([fpkmMean,DataFrameLog2In.iloc[:,item2].astype(float).apply(np.mean,axis=1)],axis=1)
        foldChange = pd.concat([foldChange,(fpkmMean.iloc[:,1] - fpkmMean.iloc[:,0])],axis=1)
        fpkmMean = pd.DataFrame()
    return foldChange
def getFoldChangeFromLog2DataframeNoRepeat(DataFrameLog2In,treatmentList):
    import pandas as pd
    import numpy as np
    #treatmentList should be [[[ctrl1],[treat1]],[[ctrl2],[treat2]]]
    #DataFrameLog2In only contains fpkmValue,
    #column names are sample names
    #indices are unrepeated numbers
    #when set index of DataFrameIn (in cufflinks fpkm results)
    #to gene name like "WBGene00197333"
    #will raise error like "ValueError: Shape of passed values is (1, 46277), indices imply (1, 46276)"
    #that is caused by duplicated indices
    #don't give index to DataFrameIn (the index will be row numbers)
    #will solve the problem
    foldChange = pd.DataFrame()
    for item in treatmentList:
        foldChange = pd.concat([foldChange,(DataFrameLog2In.iloc[:,item[1]] - DataFrameLog2In.iloc[:,item[0]])],axis=1)
    return foldChange

#calculate pcc of two list of numbers
def calculatePccFor2ColumnDataFrame(twoColumnPandasDataFrame,columnList=[0,1]):
    from scipy.stats.stats import pearsonr
    if twoColumnPandasDataFrame.shape[1] < 2:
        print("here")
        return [('nan','nan')]
    else:
        return pearsonr(twoColumnPandasDataFrame.iloc[:,columnList[0]].astype(float),twoColumnPandasDataFrame.iloc[:,columnList[1]].astype(float))


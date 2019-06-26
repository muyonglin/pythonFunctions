#####################################################################################
#
#manipulate matrixes in file
#
####################################################################################


#file name could contain path
def readFileIntoList(filename):
    contentList = map(lambda x:x.split(), open(filename).readlines())
    return contentList

#This only return the first column of a filename
#Can be used to read GeneID
def readFileIntoUniqueList(filename):
    contentList = map(lambda x:x.split(), open(filename).readlines())
    return list(set(map(tuple,contentList)))
#read names from 1 file
#extract lines of the other file based on names
#names should in the first column of nameFile
#matrixes should be numbers except row and column names
#rowname column don't have a column name
#no "\t " before col1(the name of second column)
#col1   col2    col3    col4
#row1 num   num     num     num
#row2 num   num     num     num
#nameList should mostly contain 2 column
def extractLinesByNames(nameList,matrixList,erroFileName,withClusterTag="F"):
    #    nameList = readFileIntoList(nameFile)
#    matrix = readFileIntoList(matrixFile)
    emptyDict = {}
    for line in matrixList[1:]:
        emptyDict[line[0]] = [""] * (len(matrixList[0]))
        for i in range(0,len(matrixList[0])):
            emptyDict[line[0]][i] = float(line[i+1])
    emptySubList = []
    erroHandle = open(erroFileName,"w")
    if withClusterTag == "F":
        for entrez in nameList:
            if entrez[0] in emptyDict.keys():
                emptySubList.append((entrez[0],emptyDict[entrez[0]]))
            else:
                erroHandle.write("%s\n" % (str(entrez[0])))

    else:
        for entrez in nameList:
            if entrez[0] in emptyDict.keys():
                emptySubList.append((str(entrez[0])+"_"+str(entrez[1]),emptyDict[entrez[0]]))
            else:
                erroHandle.write("%s\n" % (str(entrez[0])))
    erroHandle.close()
    return emptySubList




def writeSubmatrixByNames(nameFilename,matrixFilename,outFilename):
    nameList = readFileIntoList(nameFilename)
    matrixList = readFileIntoList(matrixFilename)
    subMatrix = extractLinesByNames(nameList,matrixList)
    outHandle = open(outFilename,"w")
    outHandle.write("GeneID\t")
    for i in range(0,len(matrixList[0])):
        outHandle.write("%s\t" % str(matrixList[0][i]))
    for tuples in subMatrix:
        outHandle.write("%s\t" % (str(tuples[0])))
        for i in range(0,len(tuples[1])):
            outHandle.write("%s\t" % (str(tuples[1][i])))
        outHandle.write("\n")
    outHandle.close()



def writeSubmatrixByNamesUnique(nameFilename,matrixFilename,outFilename,erroFileName,columnNum="F"):
    nameList = readFileIntoUniqueList(nameFilename)
    matrixList = readFileIntoList(matrixFilename)
    subMatrix = extractLinesByNames(nameList,matrixList,erroFileName,columnNum)
    outHandle = open(outFilename,"w")
    outHandle.write("GeneID\t")
    for i in range(0,len(matrixList[0])):
        outHandle.write("%s\t" % str(matrixList[0][i]))
    outHandle.write("\n")
    for tuples in subMatrix:
        outHandle.write("%s\t" % (str(tuples[0])))
        for i in range(0,len(tuples[1])):
            outHandle.write("%s\t" % (str(tuples[1][i])))
        outHandle.write("\n")
    outHandle.close()


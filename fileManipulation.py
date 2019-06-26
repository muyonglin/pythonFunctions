#!/usr/bin/env python
# encoding: utf-8
#write list of list to line
def writeListOfListToLine(listName,outFileName='out.txt',sep="\t"):
    outHandle = open(outFileName,'w')
    for item in listName:
        outHandle.write(sep.join(item))
        outHandle.write("\n")
    outHandle.close()


#write list to file, one element per line
def writeElementToLine(listName,outFileName = "values.txt"):
    outHandle = open(outFileName,"w")
    for item in listName:
        outHandle.write(str(item)+"\n")
    outHandle.close()

#parse genebank(.gb) like file
#parse kegg pathway file
#with mutiple sequencial sapces in one line
def parseGenebankLikeFile(fileHandle,fromURL=True):
    #contentList = list(map(lambda x:x.rstrip().split(sep), fileHandle.readlines()))
    contentList = fileHandle.read().decode("utf-8")
    contentList = contentList.split("\n")
    #print(contentList[0:5])
    contentDict = {}
    #for line in fileHandle.readlines():
    for line in contentList:
        #print(line[0])
        if len(line) < 1:
            continue
        elif line[0] != " ":
            itemTmp = line.rstrip().split()
            contentDict[itemTmp[0]] = []
            contentDict[itemTmp[0]].append(itemTmp[1:])
            tmpKey = itemTmp[0]
        elif line[0] == " ":
            contentDict[tmpKey].append(line.rstrip().split())
    #print(contentDict)
    return contentDict

def convertKeggID2EntrezID(ListOfListFromKegg,column=0):
    idFilePath = "/picb/molsysbio/usr/muyonglin/genomeSource/keggPathways/keggId2ncbiGeneIdForCelegans.txt"
    contentList = list(map(lambda x:x.rstrip().split("\t"), open(idFilePath).readlines()))
    contentDict = {}
    for item in contentList:
        contentDict[item[0].split(":")[1]] = item[1].split(":")[1]
    entezIDlist = []
    for  item in ListOfListFromKegg:
        #print(type(item[0]))
        if item[0] in contentDict.keys():
            #print(entezIDlist)
            #print(contentDict)
            entezIDlist.append(contentDict[item[0]])
        else:
            entezIDlist.append(" ")
    return entezIDlist

#returns a list of tuple [(),(),()]
def convertGencodeID2EntrezID(listOfGencodeID,targetIDType="GS",column=0):
    #listOfGencodeID should be [[],[],[]]
    import re
    entrezIdFilePath = "/home/muyonglin/muyonglin/genomeSource/human/douXiaoyangGenome/gencode.v25.metadata.EntrezGene"
    gtfFile = "/home/muyonglin/muyonglin/genomeSource/human/douXiaoyangGenome/gencode.v25.GRCh38.p7.gtf"
    geneSymbolFilePath = "/home/muyonglin/muyonglin/genomeSource/human/douXiaoyangGenome/gencode.v25.metadata.HGNC"
    if targetIDType == "GS":
        T2EcontentList = list(map(lambda x:x.rstrip().split("\t"), open(geneSymbolFilePath).readlines()))
    elif targetIDType == "EI":
        T2EcontentList = list(map(lambda x:x.rstrip().split("\t"), open(entrezIdFilePath).readlines()))
    else:
        print("which ID to convert to ")
        return(0)
    G2TcontentList = [] #  gene_id; transcript_id in 8
    for line in open(gtfFile).readlines():
        if re.search("\ttranscript\t",line):
            G2TcontentList.append(line.rstrip().split("\t"))
    T2EcontentDict = {}
    for item in T2EcontentList:
        T2EcontentDict[item[0]] = item[1]
    G2TcontentDict = {}
    for item in G2TcontentList:
        ENSG = item[8].split(";")[0].split(" ")[1][1:-1]
        ENST = item[8].split(";")[1].split(" ")[2][1:-1]
        G2TcontentDict[ENSG] = ENST
    entrezIDList = []
    #print(contentDict)
    #print(listOfGencodeID)
    #print(contentDict.keys())
    for item in listOfGencodeID:
        #print(item[column])
        if item[column][0] == '"':
            item[column] = item[column][1:-1]
        if item[column] in G2TcontentDict.keys():
            #print(item[column])
            tmpENST = G2TcontentDict[item[column]]
            if tmpENST in T2EcontentDict.keys():
                entrezIDList.append((item[column],tmpENST,T2EcontentDict[tmpENST]))
            else:
                entrezIDList.append((item[column],tmpENST,"NoT2E"))
        else:
            entrezIDList.append((item[column],"NoG2T"))
    return(entrezIDList)

def readGmtFileIntoDict(fileName):
    inHandle = open(fileName)
    contentList = list(map(lambda x:x.split("\t"), inHandle.readlines()))
    inHandle.close()
    #print(len(contentList))
    contentDict = {}
    for item in contentList:
        #print(item[-2:])
        contentDict[item[0]] = item[1:-1]
    return contentDict
def convertKeggID2WormbaseID(listOfKeggIDs):
    idTableFile = "/home/muyonglin/muyonglin/genomeSource/c_elegans/c_elegans.PRJNA13758.WS257.geneOtherIDs.txt"
    contentList = list(map(lambda x:x.rstrip().split("\t"), open(idTableFile).readlines()))
    #print(contentList[1:10])
    contentDict = {}
    returnList = []
    for item in contentList:
        if len(item) < 3:
            continue
        elif item[-1] in contentDict.keys():
            contentDict[item[-1]].append(item[0])
        else:
            contentDict[item[-1]] = item[0]
    for item in listOfKeggIDs:
        if item in contentDict.keys():
            if type(contentDict[item]) == list:
                for item2 in contentDict[item]:
                    returnList.append(item2)
            elif type(contentDict[item]) == str:
                returnList.append(contentDict[item])
        else:
            returnList.append(" ")
    return returnList

def convertWormbaseID2KeggID(listOfWormbaseID):
    idTableFile = "/home/muyonglin/muyonglin/genomeSource/c_elegans/c_elegans.PRJNA13758.WS257.geneOtherIDs.txt"
    contentList = list(map(lambda x:x.rstrip().split("\t"), open(idTableFile).readlines()))
    #print(contentList[1:10])
    contentDict = {}
    returnList = []
    for item in contentList:
        if len(item) < 3:
            continue
        elif item[0] in contentDict.keys():
            contentDict[item[0]].append(item[-1])
        else:
            contentDict[item[0]] = item[-1]
    for item in listOfWormbaseID:
        if item in contentDict.keys():
            if type(contentDict[item]) == list:
                for item2 in contentDict[item]:
                    returnList.append(item2)
            elif type(contentDict[item]) == str:
                returnList.append(contentDict[item])
        else:
            returnList.append(" ")
    return returnList



def writeDictIntoGmt(dictIn,fileName):
    outHandle = open(fileName,"w")
    for key in dictIn.keys():
        outHandle.write("%s\t" %(str(key)))
        for item in dictIn[key]:
            outHandle.write("%s\t" % (str(item)))
        outHandle.write("\n")
    outHandle.close()

#parse gff, extract lines within a genome region
def readGffContentIntoList(fileName):
    inHandle = open(fileName)
    contentList = []
    for line in inHandle.readlines():
        if line[0] == "#":
            continue
        elif line[0] == "":
            break
        else:
            contentList.append(line.rstrip().split("\t"))
    inHandle.close()
    return contentList

def extractLinesWithinARegionFromGff(gffFileName,regionList):
    #regionList should be [chromosome symbol,gene start,gene end]
    contentList = readGffContentIntoList(gffFileName)
    geneInfoList = []
    for item in contentList:
        if item[0] == regionList[0] and int(item[3]) >= regionList[1] and int(item[4]) <= regionList[2]:
            geneInfoList.append(item)
    return geneInfoList

#parse wig files
def readWigIntoList(fileName):
    inHandle = open(fileName)
    contentList = list(map(lambda x:x.rstrip().split("\t"), inHandle.readlines()))
    inHandle.close()
    return contentList

def extractLinesWithinARegionFromWig(wigFileName,regionList):
    contentList = readWigIntoList(wigFileName)
    geneInfoList = []
    for item in contentList[1:]:
        if item[0] == regionList[0] and int(item[1]) >= regionList[1] and int(item[2]) <= regionList[2]:
            geneInfoList.append(item)
    return geneInfoList

###############read column file into a dictionary##################
#read files into dictionary, the first column will be keys
#lines started with # will be omitted
def read_lines_to_dict(handle,sep_string,head="T"):
    i=0
    out_dict = {}
    while True:
        line = handle.readline()
        if line[0] == "#":
            continue
        if line == "":
            i +=1
            if i>100:
                print ("100 empty line")
                break
        if not line:
            break
        if line !="":
            break
    if head == "T":
        out_dict["head"]=line.rstrip().split(sep_string)
    else:
        out_dict[line.rstrip().split(sep_string)[0]]=line.rstrip().split(sep_string)[1:]
    while True:
        line = []
        line = handle.readline()
        if not line:
            return out_dict
            break
        else:
            if line == "":

                continue
            else:
                out_dict[line.rstrip().split(sep_string)[0]]=line.rstrip().split(sep_string)[1:]
    return out_dict

def read_lines_to_dict_SpecKey(fileName,sep_string,keyCol=0,head="T"):
    handle = open(fileName)
    i=0
    out_dict = {}
    while True:
        line = handle.readline()
        if line[0] == "#":
            continue
        if line == "":
            i +=1
            if i>100:
                print ("100 empty line")
                break
        if not line:
            break
        if line !="":
            break
    if head == "T":
        out_dict["head"]=line.rstrip().split(sep_string)
    else:
        out_dict[line.rstrip().split(sep_string)[keyCol]]=line.rstrip().split(sep_string)
    while True:
        line = []
        line = handle.readline()
        if not line:
            return out_dict
            break
        else:
            if line == "":

                continue
            else:
                out_dict[line.rstrip().split(sep_string)[keyCol]]=line.rstrip().split(sep_string)
    handle.close()
    return out_dict


#read multiple files into a dictionary of lists
def readMultipleFileIntoDictOfLists(fileNames,sep_string="\t"):
    Dict = {}
    for item in fileNames:
        contentList = list(map(lambda x:x.rstrip().split(sep_string), open(item).readlines()))
        Dict[item] = contentList
    return Dict

#readlines into list
def read_lines_to_list(fileName,sep_string="No",comment="#"):
    handle = open(fileName)
    i = 0
    out_list = []
    if sep_string=="No":
        while True:
            line = handle.readline()
            if line == "":
                i +=1
                if i>100:
                    print ("100 empty line")
                    break
            if not line:
                break
            if line !="":
                break
        out_list.append(line.rstrip())
        while True:
            if not line:
                return out_list
                break
            else:
                if line == "":
                    return out_list
                    break
                else:
                    line = []
                    line = handle.readline()
                    out_list.append(line.rstrip())
        handle.close()
        return out_list
    else:
        while True:
            line = handle.readline()
            if line == "":
                i +=1
                if i>100:
                    print ("100 empty line")
                    break
            if not line:
                break
            if line !="":
                break
        out_list.append(line.rstrip().split(sep_string))
        while True:
            if not line:
                return out_list
                break
            else:
                if line[0] == comment:
                    continue
                elif line == "":
                    return out_list
                    break
                else:
                    line = []
                    line = handle.readline()
                    out_list.append(line.rstrip().split(sep_string))
        handle.close()
        return out_list


#read lines into lists, omit commented lines
def read_lines_to_list_ommit_comments(fileName,sep_string="No",comment="#"):
    handle = open(fileName)
    contentList = list(map(lambda x: x.rstrip().split(sep_string),
        filter(lambda x: x[0] != comment, handle.readlines())))
    return contentList


def readSoft2Dict(softFileName,index=11):
    """
    if one prob refer to multiple genes, discard this prob
    @return: {prob:[EntrezID,EntrezID],...}
    """
    import gzip
    probe2Entrez = {}
    Flag = False
    if softFileName[-2:] == "gz":
        softHandle = gzip.open(softFileName,"rt")
    else:
        softHandle = open(softFileName,"r")
    softMatrix = softHandle.readlines()
    for line in softMatrix:
        line = line.split("\t")
        #if len(line[0]) <5 :
        #    print(line[0].lower())
        if len(line) <= index:
            continue
        if Flag:
            #print(line)
            if line[0] in probe2Entrez.keys():
                probe2Entrez[line[0]].append(line)
            else:
                probe2Entrez[line[0]] = [line]
        if line[0].lower() == 'id':
            Flag = True
    multipleKeyList = []
    for key in probe2Entrez:    #discard probs refer to multiple genes
        if len(probe2Entrez[key]) > 1:
            multipleKeyList.append(key)
    for key in multipleKeyList:    #can't del keys of dictionary when  iterating it
        del probe2Entrez[key]
    return probe2Entrez

#count reads length in a fastq file
class statReadsLength():
    #@input fastq file name
    #@output text file
    def __init__(self,fastqFileName,seqMaxLen = 51,outFileName = "None"):
        self.fq = fastqFileName
        self.out = outFileName
        self.sqM = seqMaxLen
        self.lL = self.countLength()
        if outFileName != "None":
            self.plotToFile()
    def countLength(self):
        inHandle = open(self.fq,"r")
        lengthList = [0] * self.sqM
        for line in inHandle:
            if line[0] == "@":
                line = inHandle.readline()
                lengthList[len(line)] += 1
        return(lengthList)
    def plotToFile(self):
        import matplotlib
        matplotlib.use("Agg")
        import pylab as plt
        from matplotlib.backends.backend_pdf import PdfPages
        import numpy as np
        #https://pythonspot.com/en/matplotlib-bar-chart/
        pp=PdfPages(self.out)
        x = np.arange(1,len(self.lL) + 1)
        y = self.lL
        plt.bar(x,y,align="center",color="blue")
        plt.xlabel("reads length")
        plt.ylabel("reads number")
        plt.title("reads length distribution")
        #plt.xticks(x)
        pp.savefig()
        pp.close()
#write elements in a dict to lines
def writeDictIntoTxt(dictIn,outFileName="out.txt",sep="\t",head=[]):
    #every key only has one value(not a list dict ...)
    outHandle = open(outFileName,"w")
    if len(head) != 0:
        for item in head:
            outHandle.write(str(item) + sep)
        outHandle.write("\n")
    for key in dictIn.keys():
        outHandle.write(str(key) + sep + str(dictIn[key]) + "\n")
    outHandle.close()
#parse gff files
def extractSpecificGeneTypeFromGTFFile(gtfFile,outName="test.txt",pattern=".*"):
    import os
    cmd = "awk '" + pattern +  "' " + gtfFile + " > " + outName
    os.system(cmd)


##################################################################################
#                       file system
#make directory
def mkdir(folderName):
    import os
    if not os.path.exists(folderName):
        os.makedirs(folderName)
#delete directory
def rmdir(folderName):
    import os
    if os.path.exists(folderName):
        os.system("rm -r " + folderName)

#copy specific files under the same parent folder
def copyThrough(path,pattern,dest="./"):
    import fnmatch, os
    matches = []
    for root, dirnames, filenames in os.walk(path):
        for filename in fnmatch.filter(filenames,pattern):
            matches.append(os.path.join(root,filename))
    for item in matches:
        os.system("cp %s %s" % (item,dest))
#filename could contain path
def readFileIntoList(filename,sep=" "):
    contentList = map(lambda x:x.rstrip().split(sep), open(filename).readlines())
    return list(contentList)



#!/usr/bin/env python2.7
#--coding:utf-8--

#download SRA and GSE
class geoDownload(object):
    """
    @input: GSE number list, destination directory,download type
    @action: Download specified datatype for given GSE number
    @output: None
    """
    import os, re
    from bs4 import BeautifulSoup
    import urllib
    import ssl
    #from functions import mkdir
    #import wgetter
    def __init__(self,GSENumList,destDir,dataType="sra"):
        if dataType == "sra":
            self.downloadSRA(GSENumList,destDir)
    def downloadSRA(self,GSENumList,destDir):
        for item in GSENumList:
            self.os.chdir(self.os.path.abspath(destDir))
            self.os.makedirs(item,exist_ok=True)
            self.os.chdir(item)
            self.queryGSE(item)
            self.os.chdir("../")

        #url = ""
    def queryGSE(self,GSEnum):
        url = 'http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=' + GSEnum
        #context = self.ssl._create_unverified_context()
        while True:
            try:
                context = self.ssl._create_unverified_context()
                response = self.urllib.request.urlopen(url,timeout= 1,context = context)
                html = response.read().decode("utf-8")
    #            print(response)
    #            print(html)
                break
            except:
                print ('try ' + GSEnum + ' again')
                #print(self.urllib.request.urlopen(url,timeout=1))
                #break
        soup = self.BeautifulSoup(html)
        trs = soup.findAll('tr')
        res = []
        index = 0
        for tr in trs:
            tds = tr.findAll('td')
            #f tds[0].text[0:3]
            if len(tds) == 2:
                #print(tds[0].text)
                gsm = tds[0].text
                title = tds[1].text
                if gsm.startswith('GSM'):
                    self.queryGSM(gsm)
                    #description = self.queryGSM(gsm)
                    #mid = [gsm, title, description, '1']
                    #res.append(mid)
                    #index = index + 1
                    #yield mid
        #if outFileName == 'F':
        #yield res
        #else:
        #return trs
    def queryGSM(self,GSM):
        url = 'http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=' + GSM
        print("Trying " + GSM)
        #context = self.ssl._create_unverified_context()
        while True:
            try:
                context = self.ssl._create_unverified_context()

                response = self.urllib.request.urlopen(url, timeout= 1,context = context)
                html = response.read().decode("utf-8")
                break
            except:
                print ('try '+GSM+' again')
                continue
        soup = self.BeautifulSoup(html)
        trs = soup.findAll('tr')
        #DescriptionStatus = "False"
        for tr in trs:
            tds = tr.findAll('td')
            #print(tds)
            #print(len(tds))
            if len(tds) == 4 and tds[3].text == "SRA Experiment":
                ftpPath = tds[2].find('a')['href']
                #print(tds[0].text)
                self.os.makedirs(GSM,exist_ok=True)
                self.os.chdir(GSM)
                self.os.system("wget t -c -r -nd -nH " + ftpPath)
                self.os.chdir("../")





##############################################################
#get gse info from ncbi
class getGSEInfo(object):
    """
    queryGSE will extract the description of all GSMs in a GSE
    queryGSM will download the description of one GSM
    GSEnum should be a list of GSE numbers
    """
    import os, re
    import urllib
    from bs4 import BeautifulSoup
    from Bio import Entrez
    import ssl
    def __init__(self,GSEnum,outFileName="defaultGSE.out",onlyGetGSMNum="F"):
        """TODO: to be defined1. """

        if onlyGetGSMNum == "T":
            GSMNums = []
            for iterms in GSEnum:
                result = self.getGSMNum(iterms)
                #print(result)
                GSMNums.append(result)
            #print(GSMNums)
            GSMOutHandle = open(outFileName,"w")
            for i in range(0,len(GSMNums)):
                GSMOutHandle.write("%s\n" % str(GSEnum[i]))
                for iterm in GSMNums[i]:
                    GSMOutHandle.write("%s\n" % str(iterm))
            GSMOutHandle.close()
        if onlyGetGSMNum == "F":
            outHandle = open (outFileName,'w')
            for iterms in GSEnum:
                result = self.queryGSE(iterms)
                outHandle.write(iterms + '\n')
                for lines in result:
                    #print(lines)
                    #lines1 = [x.encode("utf-8","ignore") for x in lines] #solve 'ascii' codec can't encode character u'\xb0' in position 317: error
                    #outHandle.write('\t'.join(lines1) + '\n')
                    outHandle.write('\t'.join(lines) + '\n')
            outHandle.close()
        #def __iter__(self):
    #    return results
    def getGSMNum(self,GSEnum):
        url = 'http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=' + GSEnum
        while True:
            try:
                #                context = self.ssl._create_unverified_context()
                response = self.urllib.request.urlopen(url,timeout= 1,context=context)
                html = response.read().decode("utf-8")
    #            print(response)
    #            print(html)
                break
            except:
                print ('try ' + GSEnum + ' again')
                #print(self.urllib.request.urlopen(url,timeout=1))
                #break
        soup = self.BeautifulSoup(html)
        trs = soup.findAll('tr')
        res = []
        index = 0
        for tr in trs:
            tds = tr.findAll('td')
            if len(tds) == 2:
                gsm = tds[0].text
                if gsm.startswith("GSM"):
                    res.append(gsm)
        return res
    def queryGSE(self,GSEnum):
        url = 'http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=' + GSEnum
        #context = self.ssl._create_unverified_context()
        while True:
            try:
                context = self.ssl._create_unverified_context()
                response = self.urllib.request.urlopen(url,timeout= 1,context = context)
                html = response.read().decode("utf-8")
    #            print(response)
    #            print(html)
                break
            except:
                print ('try ' + GSEnum + ' again')
                #print(self.urllib.request.urlopen(url,timeout=1))
                #break
        soup = self.BeautifulSoup(html)
        trs = soup.findAll('tr')
        res = []
        index = 0
        for tr in trs:
            tds = tr.findAll('td')
            if len(tds) == 2:
                gsm = tds[0].text
                title = tds[1].text
                if gsm.startswith('GSM'):
                    description = self.queryGSM(gsm)
                    mid = [gsm, title, description, '1']
                    res.append(mid)
                    index = index + 1
                    #yield mid
        #if outFileName == 'F':
        #yield res
        #else:
        return res
    def queryGSM(self,GSM):
        url = 'http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=' + GSM
        print("Trying " + GSM)
        #context = self.ssl._create_unverified_context()
        while True:
            try:
                context = self.ssl._create_unverified_context()

                response = self.urllib.request.urlopen(url, timeout= 1,context = context)
                html = response.read().decode("utf-8")
                break
            except:
                print ('try '+GSM+' again')
                continue
        soup = self.BeautifulSoup(html)
        trs = soup.findAll('tr')
        DescriptionStatus = "False"
        for tr in trs:
            tds = tr.findAll('td')
            #print(len(tds))
            if len(tds) == 2:
                #print(tds[0].text)
                if tds[0].text == 'Description':
                    print (GSM + ' Done!')
                    #print(tds[1].text)
                    DescriptionStatus = "True"
                    if tds[1].text == None:
                        return "No descriptional content"
                    else:
                        return tds[1].text
        if DescriptionStatus == "False":
            return "No descriptional tag"

###############################################################################
#transform probID in cel.gz idat.gz or fastaq files into geneID
class probID2EntrezID(object):
    """
    @input contain family.soft.gz, matrix file, plateformTag
        plateformTag could be "affy" or "idat"
        affy for affymatrix idat for illumina beadchip
    @work flow, read family.soft.gz and matrix file in,
                replace 1st column of matrix with gene id column in soft.gz
                write replaced matrix file out
    @output transformed matrix file, prob ID to EntrezID
    """
    def __init__(self,softFileName,matrixFileName,outFileName,sep=" ",plateformTag="affy"):
        """
        Entrez gene name in 13th column of affymatrix soft file
        Entrez gene ID in 14th column of affymatrix soft file
        Entrez gene name in 10th column of illumina beadchip soft file
        """
        import gzip, os
        import numpy as np
        print("converting %s with ID in %s" %(matrixFileName,softFileName))
        print("write into %s" %(outFileName))
        print("plateform is %s" % (plateformTag))
        index = self.platform2IDcolumn(plateformTag)
        #"rb" is binary mode (default), "rt" is text mode
        softHandle = gzip.open(softFileName,"rt")
        matrixHandle = open(matrixFileName)
        outHandle = open(outFileName,"w")
        self.renameProb(softHandle,index,matrixHandle,outHandle,sep)

    def platform2IDcolumn(self,plateformTag):
        if plateformTag == "idat":
            index = 9
        elif plateformTag == "affy":
            index = 13
        else:
            index = plateformTag
        return index
    def readSoft2Dict(self,softHandle,index=11):
        """
        if one prob refer to multiple genes, discard this prob
        @return: {prob:[EntrezID,EntrezID],...}
        """
        probe2Entrez = {}
        Flag = False
        softMatrix = softHandle.readlines()
        for line in softMatrix:
            line = line.strip().split("\t")
            #if len(line[0]) <5 :
            #    print(line[0].lower())
            if len(line) <= index:
                continue
            if Flag:
                #print(line)
                Entrezs = line[index].split("///")
                for Entrez in Entrezs:
                    tmp = Entrez.strip()
                    if len(tmp) == 0:
                        continue
                    #if probe2Entrez.has_key(line[0]):
                    if line[0] in probe2Entrez.keys():
                        probe2Entrez[line[0]].append(tmp)
                    else:
                        probe2Entrez[line[0]] = [tmp]
            if line[0].lower() == 'id':
                Flag = True
        multipleKeyList = []
        for key in probe2Entrez:    #discard probs refer to multiple genes
            if len(probe2Entrez[key]) > 1:
                multipleKeyList.append(key)
        for key in multipleKeyList:    #can't del keys of dictionary when  iterating it
            del probe2Entrez[key]
        return probe2Entrez
    def renameProb(self,softHandle,index,matrixHandle,outHandle,sep):
        """
        matrix file is seperated by " "
        """
        p2EDict = self.readSoft2Dict(softHandle,index)
        indata = list(map(lambda x:x.rstrip().split(sep), matrixHandle.readlines()))
        outHandle.write("\t")
        #print(indata[1:10])
        for iterm in indata[0]:
            outHandle.write("%s\t" % (str(iterm)))
        outHandle.write("\n")
        for line in indata[1:]:
            #if not p2EDict.has_key(line[0]):
            if line[0] not in p2EDict.keys():
                continue  #skip probs don't represent genes
            if 'NA' in line or 'Inf' in line or '-Inf' in line:
                continue # skip lines contain NA values
            outHandle.write("%s\t" % (str(p2EDict[line[0]][0])))
            for i in range(1,len(line)):
                outHandle.write("%s\t" %(str(line[i])))
            outHandle.write("\n")
##################################################################
class wormBaseID2EntrezID(object):
    def __init__(self,wbFileName,outFileName,sep=",",convertTag="ID"):
        from functions import read_lines_to_dict
        wbDict = self.readID2Dict()
        index = self.Tag2Index(convertTag)
        print("converting %s" %(str(wbFileName)))
        self.renameWBID(wbDict,wbFileName,outFileName,index,sep)
    def readID2Dict(self):
        id2EntrezFile = "/home/muyonglin/muyonglin/genomeSource/c_elegans/Caenorhabditis_elegans.gene_info"
        idHandle = open(id2EntrezFile)
        inList = list(map(lambda x:x.rstrip().split("\t"),idHandle.readlines()))
        idHandle.close()
        wbDict = {}
        #print(inList[1:5])
        #i = 0
        for item in inList[1:]:
            #print(item[5].split(":"))
            if len(item[5].split(":")) == 2:
                wbDict[item[5].split(":")[1]] = item
            else:
                #i += 1
                continue
        #print(i)
        return wbDict
    def Tag2Index(self,convertTag):
        if convertTag == "ID":
            index = 1
        if convertTag == "name":
            index = 2
        return index
    def renameWBID(self,wbDict,wbFileName,outFileName,index,sep):
        """

        """
        inHandle = open(wbFileName)
        outHandle = open(outFileName,"w")
        inData = list(map(lambda x:x.rstrip().split(sep),inHandle.readlines()))
        #print(inData[0:5])
        for item in inData[0]:
            outHandle.write("%s\t" % (str(item)))
        outHandle.write("\n")
        for line in inData[1:]:
            if line[0] not in wbDict.keys():
                #print(line[0])
                for i in line:
                    outHandle.write("%s\t" %(str(i)))
                    #keep ids don't represent genes intacted
                outHandle.write("\n")
                continue
            if 'NA' in line or 'Inf' in line or '-Inf' in line:
                continue # skip lines contain NA values
            #print(index)
            #print(wbDict[line[0]])
            if str(wbDict[line[0]][index]) == "":
                for i in line:
                    outHandle.write("%s\t" %(str(i)))
                    #keep ids don't represent genes intacted
                outHandle.write("\n")
                continue
            outHandle.write("%s\t" % (str(wbDict[line[0]][index])))
            for i in range(1,len(line)):
                outHandle.write("%s\t" % (str(line[i])))
            outHandle.write("\n")
        inHandle.close()
        outHandle.close()

##########################################################
#download supplymentary file in a GSE
def downloadGSEsupplymentary(GSEList,path="./"):
    import os, GEOparse
    from functions import mkdir
    for item in GSEList:
        mkdir(path + item)
        gse = GEOparse.get_GEO(geo = item, destdir = path + item)
        gse.download_supplementary_files(directory = path + item,email="275199027@qq.com")





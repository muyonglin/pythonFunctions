#!/usr/bin/env python
# encoding: utf-8

#get tRNA over view by their transcript ID
class tRNA2AmioAcid():
    '''
    input: files with tRNA transcript ID
    output: tRNA transcript ID and corresponding Amino Acid
    '''
    def __init__(self,tIDFile,outFile,timeout=3):
        self.inF = tIDFile
        self.oF = outFile
        self.to = timeout
    def queryWormBase(self):
        import urllib.request
        import re
        inHandle = open(self.inF,"r")
        resDict = {}
        for line in inHandle:
            tID = line.strip()
            url = "http://www.wormbase.org/rest/widget/transcript/" +\
                    tID + "/overview"
            i = 0
            print(tID + "begin!")
            while True:
                try:
                    response = urllib.request.urlopen(url,timeout = self.to)
                    html = response.read().decode("utf-8")
                    #print(html)
                    #print(tID + " done!")
                    break
                except:
                    print("try " + tID + " again")
                    i += 1
                    if i > 10:
                        print("connection failed 10 times")
                        break
            resDict[tID] = re.findall("(tRNA-.*) ",html)
            #print(re.findall("> (tRNA-.*) <",html))
        outHandle = open(self.oF,"w")
        #print(resDict)
        for key in resDict.keys():
            print(key)
            outHandle.write(str(key) + "\t")
            for item in resDict[key]:
                outHandle.write(str(item) + "\t")
            outHandle.write("\n")
        inHandle.close()
        outHandle.close()


#####################################################
#
#run and parse BICsuperkmeans
#
####################################################


#run BICsuperkmeans in python
#based on bash scripts
#path should be the absolute of a folder, the BIC result will be
#path should ended with '/'
#saved under this folder
#eg/es should be [start,end,step] three numbers
def runBICtest(path,matrixFile,corNum=1,eg=["NA"],es=["NA"],repeat=100,hc=""):
    runBICskmeansMultiCore(path,matrixFile,corNum,eg,es,repeat,hc)
    valueDictionary = parseAppropriateKValue(path,eg,es)
    runPlot(path,valueDictionary)
#multicore version runBICskmeans
def runBICskmeansMultiCore(path,matrixFile,corNum=1,eg=["NA"],es=["NA"],repeat=100,hc=""):
    import os,shutil
    from joblib import Parallel, delayed
    import numpy as np
    os.chdir(path)
    BICPath= '/home/muyonglin/software/BICskmeans/BICSKmeans.py'
    if eg[0] != "NA":
        Parallel(n_jobs=corNum)(delayed(egRun)(egValue,matrixFile,repeat,hc) for egValue in np.arange(eg[0],eg[1],eg[2]))
    if es[0] != "NA":
        Parallel(n_jobs=corNum)(delayed(esRun)(esValue,matrixFile,repeat,hc)for esValue in np.arange(es[0],es[1],es[2]))
#    print(corNum)
def egRun(egValue,matrixFile,repeat,hc):
    #    for i in np.arange(eg[0],eg[1],eg[2]):
    BICPath= '/home/muyonglin/software/BICskmeans/BICSKmeans.py'
    import os, shutil
    dirname = "eg_" + str(egValue)
    os.mkdir(dirname)
    os.chdir(dirname)
    shutil.copyfile("../"+matrixFile,matrixFile)
    os.system("python %s  -f %s  --ks=0 --eg=%d -r %d %s" %(BICPath, matrixFile, egValue,repeat,hc))
    print ("eg%d done!" % (egValue))
    os.chdir("..")
def esRun(esValue,matrixFile,repeat,hc):
    #    for i in np.arange(es[0],es[1],es[2]):
    BICPath= '/home/muyonglin/software/BICskmeans/BICSKmeans.py'
    import os,shutil
    dirname = "es_" + str(esValue)
    os.mkdir(dirname)
    os.chdir(dirname)
    shutil.copyfile("../"+matrixFile,matrixFile)
    os.system("python %s -f %s --kg=0 --es=%d -r %d %s" %(BICPath,matrixFile,esValue,repeat,hc))
    print ("es%d done!" % (esValue))
    os.chdir("..")

#single core version runBICskmeans
def runBICskmeans(path,matrixFile,eg=["NA"],es=["NA"]):
    import os,shutil
    import numpy as np
    os.chdir(path)
    if eg[0] != "NA":
        for i in np.arange(eg[0],eg[1],eg[2]):
            dirname = "eg_" + str(i)
            os.mkdir(dirname)
            os.chdir(dirname)
            shutil.copyfile("../"+matrixFile,matrixFile)
            os.system("python /home/muyonglin/software/BICskmeans/BICSKmeans.py  -f %s --ks=0 --eg=%d" %(matrixFile,i))
            print ("eg%d done!" % (i))
            os.chdir("..")
    if es[0] != "NA":
        for i in np.arange(es[0],es[1],es[2]):
            dirname = "es_" + str(i)
            os.mkdir(dirname)
            os.chdir(dirname)
            shutil.copyfile("../"+matrixFile,matrixFile)
            os.system("python /home/muyonglin/software/BICskmeans/BICSKmeans.py -f %s --kg=0 --es=%d" %(matrixFile,i))
            print ("es%d done!" % (i))
            os.chdir("..")

#parse BICskmeans results
#extract best k value from file name generated from runBICskmeans
#path is the one in runBICskmeans
#the last folder in outPath should not have "/"
#like "out=/home/muyonglin/muyonglin/sin3_microarray/mu_microarray_BICskmeans"
def parseAppropriateKValue(path,eg=["NA"],es=["NA"],egPrefix='eg_',esPrefix='es_'):
    import os,glob
    import numpy as np
    os.chdir(path)
    inlist_eg=[]
    inlist_es=[]
    if eg[0] != "NA":
        for i in np.arange(eg[0],eg[1],eg[2]):
            inlist_eg.append((i,egPrefix + str(i)))
    if es[0] != "NA":
        for i in np.arange(es[0],es[1],es[2]):
            inlist_es.append((i,esPrefix + str(i)))
#    print(inlist_eg)
    row_es=[]
    row_eg=[]
    col_kg=[]
    col_ks=[]
    for value, dirname in inlist_eg:
        filenames = glob.glob(dirname + '/*.kgg')
#        print(filenames)
        if len(filenames) > 1 :
            print ("There are more than one *.kgg file in %s, please check." % (dirname))
        if len(filenames ) == 1:
            num_kg = filenames[0].split('/')[-1].split('_K_G')[-1].split('.kgg')[0]
            row_eg.append(dirname.split('_')[-1])
            col_kg.append(num_kg)
#            print(numkg)
    for value, dirname in inlist_es:
        filenames = glob.glob(dirname + '/*.kag')
        if len(filenames) > 1 :
            print ("There are more than one *.kgg file in %s, please check." % (dirname))
        if len(filenames ) == 1:
            num_ks = filenames[0].split('/')[-1].split('_K_A')[-1].split('.kag')[0]
            row_es.append(dirname.split('_')[-1])
            col_ks.append(num_ks)
    valueDict = {"es" : row_es, "eg" : row_eg , "ks" : col_ks , "kg" : col_kg}
    return valueDict


#when looking for best k valuse
#we should plot eg2kg and es2ks curve and find the inflexion
#use matplot run the scripts imported this function with ipython
#ipython contains scipy and numpy (in anaconda)
#xlist should be es or eg
#ylist should be ks or kg
#xlist and ylist should be a list of numbers (the same length)
def plotTrendCurve(xlist,ylist,pdfNames="kg2eg.pdf",xlabel='EG',ylabel="KG",title="KG-2-EG plot"):
    import matplotlib
    matplotlib.use('Agg')
    import pylab as plt
    from matplotlib.backends.backend_pdf import PdfPages

    #plt.switch_backend('agg') # without this it will show (RuntimeError: Invalid DISPLAY variable)
#    print("plot begin")
    pp = PdfPages(pdfNames)
    #change strings into numbers
    xlist = list(map(float,xlist))
    ylist = list(map(float,ylist))

    plt.figure()
    plt.plot(xlist,ylist,'ro',lw=2.0)
    #p1.plot(xlist,ylist,'ro')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    xlistc = (max(xlist) - min(xlist)) / 10
    ylistc = (max(ylist) - min(ylist)) / 10
    plt.axis([min(xlist) - xlistc, max(xlist) + xlistc, min(ylist) - ylistc, max(ylist) + ylistc])
    pp.savefig()
    pp.close()
#    print(outfile)


#run plot with the results of parseAppropriateKValue
#value input is a dictionary
#path ended with '/'
def runPlot(path,valueDict,outfile=".pdf"):
    import os
    os.chdir(path)
    row_es=valueDict['es']; row_eg=valueDict['eg']
    col_kg=valueDict['kg']; col_ks=valueDict['ks']
#    print(row_es)
#    print(row_eg)
    if row_es != [] and col_ks != []:
        plotTrendCurve(row_es,col_ks,path + "ks2es" + outfile,"es","ks","KS-2-ES plot")
    elif row_eg != [] and col_kg != []:
        plotTrendCurve(row_eg,col_kg,path + "kg2eg" + outfile,"eg","kg","KG-2-EG plot")
    print("plot done!")

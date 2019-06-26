#!/usr/bin/env python
# encoding: utf-8
#GO enrichment with david
def davidEnrichmentAnalysis(geneList,typeIn='ENSEMBL_GENE_ID',annotIn='KEGG_PATHWAY',toolsIn='chartReport'):
    from webRelated import tryUrl
    url = 'http://david.abcc.ncifcrf.gov/api.jsp?type=' + typeIn \
            + '&ids=' + ",".join(geneList)\
            + '&tool=' + toolsIn\
            + '&annot=' + annotIn
    html = tryUrl(url)

#samtools index
def samtoolsIndex(fileName):
    import os
    cmd = "samtools index " + fileName
    print(cmd)
    os.system(cmd)

#uncompress sra, compatible for multiple process
def sra2fastq(fileName,outPath="./",split="--split-3"):
    import os
    #fileName should contains ".sra"
    cmd = "fastq-dump " + split + " " + fileName + " -O " + outPath
    print(cmd)
    os.system(cmd)

def gzip2fastq(fileName):
    import os
    cmd = "gzip -dc " + fileName + " > " + fileName[0:-3]
    print(cmd)
    os.system(cmd)

def fastqc(fileName,outDir="./fastqcOut"):
    import os
    cmd = "fastqc " + fileName +" -o " + outDir
    print(cmd)
    os.system(cmd)

#mapping
def runHisat2(seqFile,spliceFile,indexPath,outDir="./hisat2Out/",threads = 3):
    import os
    #-U for unpaired reads
    #print("begin")
    cmd = "hisat2 --known-splicesite-infile " + spliceFile \
            + " -x " + indexPath + " -U " + seqFile \
            + " -S " + outDir + seqFile.split(".")[0].split("/")[-1] + ".sam" \
            + " -p " + str(threads)
    print(cmd)
    os.system(cmd)

def runTophat2(seqFile,gtfFile,indexPath,outDir="./tophat2Out/" ,threads = 3):
    import os
    #from functions import mkdir
    #mkdir(outDir  + seqFile.split(".")[0].split("/")[-1])
    cmd  = "tophat2 -G " + gtfFile \
            + " -o " + outDir  + seqFile.split(".")[0].split("/")[-1] \
            + " -p " + str(threads)\
            + " " + indexPath + " " + seqFile
    print(cmd)
    os.system(cmd)
def runSTAR(seqFile,gtfFile,indexPath,outDir = "./STARout/" ,threads = 3,pairEnd = "F",threepTrim='0',fivepTrim='3'):
    #for pair end data, seqFile should be ["seqFile1","seqFile2"]
    #for single end data, seqFile should be "seqFile"
    if pairEnd == "T":
        seqFileInput = seqFile[0] + " " + seqFile[1]
    if pairEnd == "F":
        seqFileInput = seqFile
    import os
    """
    --limitOutSJcollapsed 3000000, default 1000000, throw error for 120M reads
    --outSAMstrandField intronMotif, solve error "BAM record error: found spliced alignment without XS attribute
    """
    cmd= "STAR --runMode alignReads --runThreadN " + str(threads) \
            + " --genomeDir " + indexPath + " --readFilesIn "  + seqFileInput \
            + " --clip5pNbases " + fivepTrim \
            + " --outFileNamePrefix " + outDir \
            + seqFileInput.split(".")[0].split("/")[-1] \
            + " --outSAMtype BAM SortedByCoordinate --outSAMreadID Number " \
            +" --outFilterMismatchNmax 10   " \
            + " --sjdbGTFfile " + gtfFile \
            + " --quantMode GeneCounts --outWigType wiggle "\
            +" --outWigStrand Stranded" \
            +" --limitOutSJcollapsed 3000000" \
            +" --outSAMstrandField intronMotif" \
            +" --clip3pNbases " + threepTrim
    print(cmd)
    os.system(cmd)

#run bwa mem
def runBwaMem(seqFile,dbPrefix,,pairEnd='T')
    bwa = '/picb/molsysbio/usr/muyonglin/ts/wesPipe/bwa/bwa-0.7.17/bwa mem'
    import os
    if pairEnd == 'T':
        cmd = bwa + " -k 40 " + dbPrefix + ' ' +seqFile[0] + " " + seqFile[1]
    elif pairEnd == 'F':
        cmd = bwa + " -k 40 " + dbPrefix + ' ' +seqFile[0]
    print(cmd)
    os.system(cmd)

def runCufflinks(alignedFile,gtfFile,outPath="./cufflinksOut/"):
    import os
    cmd = "cufflinks -G " + gtfFile \
            +" -o " + outPath \
            +" " + alignedFile
    print(cmd)
    os.system(cmd)

#merge pair end sequencing results
def runBBmerge(inFile1,inFile2,adapter1="",adapter2="",outDir="",mermory="2000m"):
    import os
    out = inFile1 + inFile2
    outu1 = inFile1 + "_unmerged_1.fq"
    outu2 = inFile2 + "_unmerged_2.fq"
    ihist = inFile1 + inFile2 + ".txt"
    cmd = "/home/muyonglin/muyonglin/software/bbmap/bbmerge.sh in1=" +\
            inFile1 + " in2=" + inFile2 +\
            " out=" + outDir + out +\
            " outu1=" + outDir + outu1 +\
            " outu2=" + outDir + outu2 +\
            " ihist=" + outDir + ihist +\
            " adapter1=" + adapter1 +\
            " adapter2=" + adapter2
    print(cmd)
    os.system(cmd)
#run preRanked GSEA
def runPrerankedGSEA(gmtFile,rnkFile,outDir):
    import os
    cmd = '''java -Xmx1024m xtools.gsea.GseaPreranked \
                -gmx %s \
                -collapse false \
                -mode Median_of_probes \
                -norm meandiv -nperm 1000 \
                -rnk %s \
                -scoring_scheme weighted \
                -rpt_label %s \
                -include_only_symbols true -make_sets true -plot_top_x 20 -rnd_seed timestamp -set_max 5000 -set_min 5 -zip_report false \
                -out %s''' % (gmtFile,rnkFile,rnkFile,outDir)
    os.system(cmd)

#run ShortStack
def runShortStack(fastqFile,outdir,adapter="",corNum=2):
    import os
    if adapter == "":
        cmd = "perl /home/muyonglin/perlFunctions/ShortStack " \
            " --outdir " + outdir +\
            " --bowtie_cores " + str(corNum) + \
            " --sort_mem 2G --readfile " + fastqFile +\
            " --genomefile /home/muyonglin/muyonglin/genomeSource/c_elegans/c_elegans.PRJNA13758.WS257.genomic_softmasked.fa"
    else:
        cmd = "perl /home/muyonglin/perlFunctions/ShortStack --adapter " + adapter + \
            " --outdir " + outdir +\
            " --bowtie_cores " + str(corNum) + \
            " --sort_mem 2G --readfile " + fastqFile +\
            " --genomefile /home/muyonglin/muyonglin/genomeSource/c_elegans/c_elegans.PRJNA13758.WS257.genomic_softmasked.fa"
    os.system(cmd)

#run bowtie2
def runBowtie(fastqFile,outdir,indexPath,seqType='PE',trim5=0,trim3=0,misMatchNum=2,threadsNum=3):
    import os
    if seqType == 'PE':
        cmd = "bowtie -m 1 -v 2 -X 5000 -S -p " + str(threadsNum) \
                + " " + indexPath \
                + " -5 " + str(trim5) \
                + " -3 " + str(trim3) \
                + " -1 " + fastqFile[0] \
                + " -2 " + fastqFile[1] \
                + " > " + outdir +fastqFile[0].split('_')[0] + "_bowtieAlign.sam"
    if seqType == 'SE':
        cmd = "bowtie -m 1 -v 2 -X 5000 -S  -p " + str(threadsNum) \
                + " " + indexPath \
                + " " + fastqFile \
                + " > " + outdir + fastqFile.split('_')[0] + "_bowtieAlign.sam"
    print(cmd)
    os.system(cmd)

#run macs2
def runMACS2Callpeak(bamFile,outpath,species="hs",alignFormat="SAM",pairState='S'):
    import os
    #will throw erro if tagSize is auto determined, 13 in the case
    #macs2 can't auto detect input format unless give it -f
    if pairState == 'P':
        #bamFile in [control.bam,treat.bam] fromat
        cmd = "macs2 callpeak -c " + bamFile[0] + " -t " + bamFile[1]\
                + " -g " + species + " -f " + alignFormat\
                + " --outdir " + outpath + " -B --SPMR -n "\
                + '_'.join([bamFile[0].split("/")[-1],bamFile[1].split("/")[-1]])
    if pairState == 'S':
        cmd = "macs2 callpeak -t " + bamFile + " -g " + species \
                + " -f " + alignFormat \
                + " --outdir " + outpath + " -B --SPMR -n " \
                + bamFile.split("/")[-1]
    print(cmd)
    os.system(cmd)

#run bsmap
def runBsmap(fastqFile,outFile,referenceFasta,seqType='PE',trimLen=144,misMatchNum=3,processsorNum=8):
    #trimming option is not set
    #use unique map
    #outFile should b3 *.sam or *.bam or *.bsp
    #trimLen, -L in bsmap, only map the first N nucleotides
    import os
    bsmap = "/home/muyonglin/muyonglin/seqDataAnaPipeline/bisulfiteSeq/SMAPdigger/bin/bsmap"
    if seqType=='PE':
        cmd = bsmap + ' -a ' + fastqFile[0] \
                +' -b ' + fastqFile[1]\
                +' -d ' + referenceFasta\
                +' -o ' + outFile\
                +' -L ' + str(trimLen)\
                +' -w 10 -v ' + str(misMatchNum) + ' -p ' + str(processsorNum) +' -n 1 -r 0'
    elif seqType == 'SE':
        cmd = bsmap + ' -a ' + fastqFile \
                +' -d ' + referenceFasta\
                +' -o ' + outFile\
                +' -L ' + str(trimLen)\
                +' -w 10 -v ' + str(misMatchNum) + ' -p ' + str(processsorNum) +' -n 1 -r 0'
    print(cmd)
    os.system(cmd)


#################################################################################
#
#                       input preprocess tools
#
##############################################################################
def convertFirstColumnToGmt(fileList,gmtFileName,sep=","):
    #fileList: [["fileName"],["fileName"],...,["fileName"]]
    outDict = {}
    from functions import readFile1thColIntoList
    for item in fileList:
        inList = readFile1thColIntoList(item,sep)
        outDict[item.split("/")[-1]] = inList
    outHandle = open(gmtFileName,"w")
    for key in outDict.keys():
        outHandle.write("%s\t" % (str(key)))
        for item1 in outDict[key]:
            outHandle.write("%s\t"%(str(item1)))
        outHandle.write("\n")

def preRankedrnkFromFirstTwoColumn(dataFileName,rnkFileName,sep="\t",columnList=[0,1],head=True):
    #read  2 column in
    #rank by the second column
    import pandas as pd
    inList = list(map(lambda x:x.rstrip().split(sep), open(dataFileName).readlines()))
    #print(pd.DataFrame(inList).iloc[:5,:])
    if head:
        inList = pd.DataFrame(inList).iloc[1:,columnList]
    else:
        inList = pd.DataFrame(inList).iloc[:,columnList]
    inList.columns = ["geneID","foldChange"]
    inList.iloc[:,1] = inList.iloc[:,1].astype(float)
    #print(inList.iloc[:5,:5])
    inList = inList.sort_values(["foldChange"],ascending=[False])
    inList.to_csv(rnkFileName,sep="\t",index=False,header=False)




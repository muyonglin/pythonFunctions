#!/usr/bin/env python
# encoding: utf-8
###############################################
#
#  dealing with bam and alignment results
#
##############################################

#counts reads in a region
class countReadsNum():
    #@input alignment results(bam/sam),
    #    gtf File (only contains region to be counted)
    #@output dictionary of geneName and counts
    def __init__(self,bamAlignmentFile,gtfFile,geneType="tRNA"):
        self.alF = bamAlignmentFile
        self.gtF = gtfFile
        self.gTy = geneType

    def parseAlignmentGtfAndCount(self):
        import HTSeq
        gtf = HTSeq.GFF_Reader(self.gtF,end_included=True)
        bam = HTSeq.BAM_Reader(self.alF)
        genes = HTSeq.GenomicArrayOfSets( "auto", stranded=False )
        counts = {}
        for feature in gtf:
            if "gene_biotype" in feature.attr.keys():
                if feature.type == "gene" and feature.attr["gene_biotype"] == self.gTy:
                    genes[feature.iv] += feature.name
                    counts[feature.name] = 0
        for aln in bam:
            if aln.aligned:
                iset = None
                for iv2, step_set in genes[aln.iv].steps():
                    if iset is None:
                        iset = step_set.copy()
                    else:
                        iset.intersection_update(step_set)
                if len(iset) == 1:
                    counts[list(iset)[0]] += 1
        return(counts)
def calculateFPKM(bam,gtf,readsNum=1,geneType="tRNA"):
    import HTSeq
    gtf = HTSeq.GFF_Reader(gtf,end_included=True)
    bam = HTSeq.BAM_Reader(bam)
    genes = HTSeq.GenomicArrayOfSets( "auto", stranded=False )
    counts = {}
    lengthOfReads = {}
    for feature in gtf:
        if "gene_biotype" in feature.attr.keys():
            if feature.type == "gene" and feature.attr["gene_biotype"] == geneType:
                genes[feature.iv] += feature.name
                counts[feature.name] = 0
                lengthOfReads[feature.name] = feature.iv.end - feature.iv.start
    for aln in bam:
        if aln.aligned:
            iset = None
            for iv2, step_set in genes[aln.iv].steps():
                if iset is None:
                    iset = step_set.copy()
                else:
                    iset.intersection_update(step_set)
            if len(iset) == 1:
                counts[list(iset)[0]] += 1
    rpkm = {}
    for key in counts.keys():
        #rpkm reads per million per kilobase
        rpkm[key] = counts[key] * 1000 / (lengthOfReads[key] * readsNum)
    return(rpkm)
def calculateFPKMFromtRNAAnnotation(bam,gtf,readsNum=1):
    import HTSeq
    gtf = HTSeq.GFF_Reader(gtf,end_included=True)
    bam = HTSeq.BAM_Reader(bam)
    genes = HTSeq.GenomicArrayOfSets( "auto", stranded=False )
    counts = {}
    lengthOfReads = {}
    for feature in gtf:
        genes[feature.iv] += feature.attr['transcript_id']
        counts[feature.attr['transcript_id']] = 0
        lengthOfReads[feature.attr['transcript_id']] = feature.iv.end - feature.iv.start
    for aln in bam:
        if aln.aligned:
            iset = None
            for iv2, step_set in genes[aln.iv].steps():
                if iset is None:
                    iset = step_set.copy()
                else:
                    iset.intersection_update(step_set)
            if len(iset) == 1:
                counts[list(iset)[0]] += 1
    rpkm = {}
    t2ADict = {}
    with open("/home/muyonglin/muyonglin/genomeSource/c_elegans/WS257_tRNA_tid_AminoAcid.txt") as f:
        for line in f:
            if len(line.strip().split('\t')) < 2:
                t2ADict[line.strip().split('\t')[0]] = line.strip().split('\t')[0]
            else:
                t2ADict[line.strip().split('\t')[0]] = line.strip().split('\t')[1]
    for key in counts.keys():
        #rpkm reads per million per kilobase
        rpkm[key] =( t2ADict[key],counts[key] * 1000 / (lengthOfReads[key] * readsNum))
    return(rpkm)


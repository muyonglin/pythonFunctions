#!/usr/bin/env python
# encoding: utf-8

class runGSEA(object):
    """
    runGSEA
    @input gmtFile rnkFile and output directory should all be list
    @output gsea results
    """
    import os
    from joblib import Parallel, delayed
    def __init__(self,gmtFile,rnkFile,outDir,preRanked = "T",corNum=1):
        if preRanked == "T":
            self.runPrerankedGSEAwithMulticore(gmtFile,rnkFile,outDir,corNum)

    def runPrerankedGSEA(self,gmtFile,rnkFile,outDir):
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
        self.os.system(cmd)
    def runPrerankedGSEAwithMulticore(self,gmtFile,rnkFile,outDir,corNum):
        self.Parallel(n_jobs=corNum)(self.delayed(self.runPrerankedGSEA)(gmtFile[i],rnkFile[i],outDir[i]) for i in range(len(rnkFile)))
def runPrerankedGSEA(gmtFile,rnkFile,outDir,maxSetNum=7000):
    import os
    cmd = '''java -Xmx4g xtools.gsea.GseaPreranked \
                -gmx %s \
                -collapse false \
                -mode Median_of_probes \
                -norm meandiv -nperm 1000 \
                -rnk %s \
                -scoring_scheme weighted \
                -rpt_label %s \
                -include_only_symbols true -make_sets true -plot_top_x 20 -rnd_seed timestamp -set_max %s  -set_min 1 -zip_report false \
                -out %s''' % (gmtFile,rnkFile,rnkFile,str(maxSetNum),outDir)
    os.system(cmd)


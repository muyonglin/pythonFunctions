#!/usr/bin/env python
# encoding: utf-8

#read fpkm from cullinks gene.fpkm_tracking
def extractFPKM(fileName,geneName):
    handle = open(fileName,"r")
    for lines in handle.readlines():
        if lines.split("\t")[0] == geneName:
            return lines.split("\t")[-4]
    handle.close()



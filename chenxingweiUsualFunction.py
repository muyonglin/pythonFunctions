##########################################################################
#
#some reused funcitons in chenxingwei's inflammation2cancer codes
#
########################################################################

__author__="Mu Yonglin"
__date__ = "2016-06-03"
__email__="yonglin.mu@gmail.com"

def readfilelines(filename):
    tmp = map(lambda x:x.split()[0], open(filename).readlines())
    return tmp

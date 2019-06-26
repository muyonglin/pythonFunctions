#!/usr/bin/env python
# encoding: utf-8
#get the source of webpage
def tryUrl(url,outTime=10,sleepTime=10,encoding='utf-8'):
    import urllib.request
    from time import sleep
    res = urllib.request.urlopen(url)
    html = res.read().decode(encoding)
    sleep(sleepTime)
    return(html)


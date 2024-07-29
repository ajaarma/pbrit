#!/usr/bin/python

import re,sys,getopt,os,datetime
sys.path.append("./CONFIG/")
from CONFIG import *


def processArgs(args):

    opts = {}
    xml_file=[]
    try:
        opts,args = getopt.getopt(args,'a:',['help'])

    except getopt.GetoptError as err:
        usage(err)

    for opt,arg in opts:
        if opt=='-h' or opt == '--help':
            usage()
        if opt == '-a':
            xml_file = arg

    return xml_file

def usage(err=[]):

    if err:
        print err
        print "python processInit.py -a <analysis-xml-file> "
    else:
        print "python processInit.py -a <analysis-xml-file> "


if __name__=="__main__":

    xml_file = processArgs(sys.argv[1:])

    obj = CONFIG()
    configDict = obj.getConfigDict(xml_file)
    
    version = configDict["version"]

    gen = configDict["general"]
    dataDir = os.path.abspath(gen['dataDir'])
    annoSourceList = gen['annoSource']['value']
    annoDirList = gen['annoDir']['value']

    annoOnt = ['HPO','DO','GO','MPO']
    for anno in annoSourceList:
        path = dataDir+"/"+version+"/"+anno
        for ele in annoDirList:
            if ele=='OBO':
                if anno in annoOnt:
                    pathD = path+"/"+ele
                    os.system("mkdir -p "+ pathD)
            else:
                #pass
                pathD = path+"/"+ele
                os.system("mkdir -p "+pathD)

    pathMap_tfidf = dataDir+'/'+version+'/MAP/MATRICES/TFIDF'
    pathMap_svd = dataDir+'/'+version+'/MAP/MATRICES/SVD'
    os.system('mkdir -p '+pathMap_tfidf)
    os.system('mkdir -p '+pathMap_svd)


    #goDict = configDict["GO"]


#!/usr/bin/python

import re,sys,getopt,os,datetime
sys.path.append("./PUB/")
sys.path.append("./CONFIG/")
sys.path.append("./MAT/")

from tfidfAbstract import *
from CONFIG import *
from MATRICES import *

def processArgs(args):

    opts = {}
    xml_file=anno_type=[]
    try:
        opts,args = getopt.getopt(args,'a:t:',['help'])

    except getopt.GetoptError as err:
        usage(err)

    for opt,arg in opts:
        if opt=='-h' or opt == '--help':
            usage()
        elif opt == '-a':
            xml_file = arg            
        elif opt=='-t':
            anno_type = arg
    
    return xml_file,anno_type

def usage(err=[]):
    
    if err:
        print "\n\n"
        print err
        print "python processABST.py -a <analysis-xml-file> -t PUB\n\n"
    else:
        print "\nPlease enter a valid command: \n"
        print "python processABST.py -a <analysis-xml-file> -t PUB\n\n"
    sys.exit()


if __name__=="__main__":

    try:
        xml_file,anno_type = processArgs(sys.argv[1:])
        if not xml_file:
            usage()
        if not anno_type:
            usage()

    except TypeError:
        usage()

    obj = CONFIG()
    configDict = obj.getConfigDict(xml_file)

    ontDict = configDict[anno_type]

    for keys in ontDict:
        print keys,"\t",ontDict[keys]


    pathAnno = os.path.abspath(configDict["general"]["dataDir"])+"/"+configDict["version"]+"/"+anno_type

    ontDict = configDict["PUB"]
    annoFile = pathAnno+"/"+ontDict["annoFile"]
    #noMesh = pathAnno+"/"+ontDict["noMesh"]
    #tmp = pathAnno+"/"+ontDict["tmp"]

    abstStem = pathAnno+"/"+ontDict["abstStem"]
    stopWord = pathAnno+"/"+ontDict["stopWord"]
    matRow = pathAnno+"/"+ontDict["matRow"]
    matCol = pathAnno+"/"+ontDict["matCol"]
    matSpr = pathAnno+"/"+ontDict["matSparse"]
    
   
    print "\nCheck if all the input, output files are correctly specified as present in the Analysis.xml file\n"
    print "gene pubmed abstract raw File: ",annoFile
    print "gene pubmed processed stemmed file: ",abstStem
    print "List of Stop Words in the Current Abstract: ",stopWord
    print "gene row file",matRow
    print "gene column file",matCol
    print "gene sparse file",matSpr
    print "\n"

    #sys.exit()
    objAbst = AbstractProcess()
    objMat = MATRICES()

    print "Processsing Raw File and generating Sparse Matrices for: "+anno_type
    objAbst.processRawPubFile(annoFile,abstStem,stopWord,matRow,matCol,matSpr,objMat)

    #print "\nGenerating Sparse Matrices for Annotation source: "+anno_type
    #objMat.getSparseMat(genePP,matRow,matCol,matSpr)

    print "\nGenerating PBS script for running TFIDF score computation: "+anno_type
    pbs_file = obj.writePBSScriptFile(pathAnno,anno_type,configDict)
    print "\nAllDone\n"



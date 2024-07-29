#!/usr/bin/python

import re,sys,getopt,os,datetime
sys.path.append("./PATH/")
sys.path.append("./CONFIG/")
sys.path.append("./MAT/")

from PYMethod import *
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
        print "python processPATH.py -a <analysis-xml-file> -t PATH\n\n"
    else:
        print "\nPlease enter a valid command: \n"
        print "python processPATH.py -a <analysis-xml-file> -t PATH \n\n"
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

    ontDict = configDict["PATH"]
    annoFile = pathAnno+"/"+ontDict["annoFile"]
    #noMesh = pathAnno+"/"+ontDict["noMesh"]
    #tmp = pathAnno+"/"+ontDict["tmp"]


    genePY = pathAnno+"/"+ontDict["genePY"]
    matRow = pathAnno+"/"+ontDict["matRow"]
    matCol = pathAnno+"/"+ontDict["matCol"]
    matSpr = pathAnno+"/"+ontDict["matSparse"]
    
   
    print "\nCheck if all the input, output files are correctly specified as present in the Analysis.xml file\n"
    print "anno raw File: ",annoFile
    print "gene anno Pheno: ",genePY
    print "gene row file",matRow
    print "gene column file",matCol
    print "gene sparse file",matSpr
    print "\n"

    ontMeth = PYMethod()
    objMat = MATRICES()

    print "Processsing Raw File for Annotation source: "+anno_type
    ontMeth.processRawPYFile(annoFile,genePY)
    print "\nGenerating Sparse Matrices for Annotation source: "+anno_type
    objMat.getSparseMat(genePY,matRow,matCol,matSpr)

    print "\nGenerating PBS script for running TFIDF score computation: "+anno_type
    pbs_file = obj.writePBSScriptFile(pathAnno,anno_type,configDict)
    print "\nAllDone\n"


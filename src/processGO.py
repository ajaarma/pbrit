#!/usr/bin/python

import re,sys,getopt,os,datetime
sys.path.append("./ONT/")
sys.path.append("./GO/")
sys.path.append("./CONFIG/")
sys.path.append("./MAT/")
from ONTOLOGY import *
from GOMethod import *
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
        print "python processGO.py -a <analysis-xml-file> -t <anno-source>"
    else:
        print "\nPlease enter a valid command: "
        print "python processGO.py -a <analysis-xml-file> -t <anno-source> \n\n"
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

    goDict = configDict[anno_type]

    for keys in goDict:
        print keys,"\t",goDict[keys]


    pathAnno = os.path.abspath(configDict["general"]["dataDir"])+"/"+configDict["version"]+"/"+anno_type
    annoFile = pathAnno+"/"+goDict["annoFile"]
    obo = pathAnno+"/"+goDict["obo"]
    geneGO = pathAnno+"/"+goDict["geneGO"]
    geneGOPar = pathAnno+"/"+goDict["geneGOPar"]
    geneNoGo = pathAnno+"/"+goDict["geneNoGo"]
    matRow = pathAnno+"/"+goDict["matRow"]
    matCol = pathAnno+"/"+goDict["matCol"]
    matSpr = pathAnno+"/"+goDict["matSparse"]
    
    
    print annoFile
    print obo
    print geneGO
    print geneGOPar
    print geneNoGo
    print matRow
    print matCol
    print matSpr

    oboEn = OBOEngine()
    goMeth = GOMethod()
    objMat = MATRICES()

    #go_hash = oboEn.getGOHashTable(obo)
    #goMeth.processRawGOFile(annoFile,geneGO)
    #goMeth.getFile_GO(go_hash, geneGO, geneGOPar,geneNoGo,"T")
    #objMat.getSparseMat(geneGOPar,matRow,matCol,matSpr)
    pbs_file = obj.writePBSScriptFile(pathAnno,anno_type,configDict)

    #go_hash = obj.getHPOHashTable(inp_file)
    #go_hash = obj.getDOHashTable(inp_file)
    #go_hash = obj.getMPOHashTable(inp_file)


    #for keys in go_hash:
    #    print keys,"\t",go_hash[keys]


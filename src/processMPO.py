#!/usr/bin/python

import re,sys,getopt,os,datetime
sys.path.append("./ONT/")
sys.path.append("./MPO/")
sys.path.append("./CONFIG/")
sys.path.append("./MAT/")
from ONTOLOGY import *
from MPOMethod import *
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
        print "python processMPO.py -a <analysis-xml-file> -t MPO\n\n"
    else:
        print "\nPlease enter a valid command: \n"
        print "python processMPO.py -a <analysis-xml-file> -t MPO \n\n"
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

    annoFile_list = ontDict["annoFile"]["value"]
    annoFile_hmd = pathAnno+"/"+annoFile_list[0]
    annoFile_mgi = pathAnno+"/"+annoFile_list[1]

    obo = pathAnno+"/"+ontDict["obo"]
    geneONT = pathAnno+"/"+ontDict["geneONT"]
    geneONTPar = pathAnno+"/"+ontDict["geneONTPar"]
    geneNoONT = pathAnno+"/"+ontDict["geneNoONT"]
    matRow = pathAnno+"/"+ontDict["matRow"]
    matCol = pathAnno+"/"+ontDict["matCol"]
    matSpr = pathAnno+"/"+ontDict["matSparse"]
    
    
    print annoFile_list
    print annoFile_hmd
    print annoFile_mgi
    print obo
    print geneONT
    print geneONTPar
    print geneNoONT
    print matRow
    print matCol
    print matSpr

    oboEn = OBOEngine()

    ontMeth = MPOMethod()
    objMat = MATRICES()

    ont_hash = oboEn.getMPOHashTable(obo)
    ontMeth.processRawMPOFile(annoFile_hmd,annoFile_mgi,geneONT)
    

    ontMeth.getFile_MPO(ont_hash,geneONT, geneONTPar,geneNoONT,"T")
    objMat.getSparseMat(geneONTPar,matRow,matCol,matSpr)
    pbs_file = obj.writePBSScriptFile(pathAnno,anno_type,configDict)
    
    
    #go_hash = obj.getHPOHashTable(inp_file)
    #go_hash = obj.getDOHashTable(inp_file)
    #go_hash = obj.getMPOHashTable(inp_file)


    #for keys in go_hash:
    #    print keys,"\t",go_hash[keys]


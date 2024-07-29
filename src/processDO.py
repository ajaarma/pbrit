#!/usr/bin/python

import re,sys,getopt,os,datetime
sys.path.append("./ONT/")
sys.path.append("./DO/")
sys.path.append("./CONFIG/")
sys.path.append("./MAT/")
from ONTOLOGY import *
from DOMethod import *
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
        print "python processHPO.py -a <analysis-xml-file> -t <anno-source>"
    else:
        print "\nPlease enter a valid command: "
        print "python processHPO.py -a <analysis-xml-file> -t <anno-source> \n\n"
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
    annoFile_exp = pathAnno+"/"+annoFile_list[0]
    annoFile_txt = pathAnno+"/"+annoFile_list[1]
    annoFile_knl = pathAnno+"/"+annoFile_list[2]

    obo = pathAnno+"/"+ontDict["obo"]
    geneONT = pathAnno+"/"+ontDict["geneONT"]
    geneONTPar = pathAnno+"/"+ontDict["geneONTPar"]
    geneNoONT = pathAnno+"/"+ontDict["geneNoONT"]
    matRow = pathAnno+"/"+ontDict["matRow"]
    matCol = pathAnno+"/"+ontDict["matCol"]
    matSpr = pathAnno+"/"+ontDict["matSparse"]
    
    
    print annoFile_list
    print annoFile_exp
    print annoFile_txt
    print annoFile_knl
    print obo
    print geneONT
    print geneONTPar
    print geneNoONT
    print matRow
    print matCol
    print matSpr

    oboEn = OBOEngine()
    ontMeth = DOMethod()
    objMat = MATRICES()

    ont_hash = oboEn.getDOHashTable(obo)
    ontMeth.processRawDOFile(annoFile_exp,annoFile_txt,annoFile_knl,geneONT)
    ontMeth.getFile_DO(ont_hash,geneONT, geneONTPar,geneNoONT,"T")
    objMat.getSparseMat(geneONTPar,matRow,matCol,matSpr)
    pbs_file = obj.writePBSScriptFile(pathAnno,anno_type,configDict)

    #go_hash = obj.getHPOHashTable(inp_file)
    #go_hash = obj.getDOHashTable(inp_file)
    #go_hash = obj.getMPOHashTable(inp_file)


    #for keys in go_hash:
    #    print keys,"\t",go_hash[keys]


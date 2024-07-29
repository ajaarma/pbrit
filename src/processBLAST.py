#!/usr/bin/python

import re,sys,getopt,os,datetime

script_path = os.path.dirname(os.path.abspath( __file__ ))

sys.path.append(script_path+"/SEQ/")
sys.path.append(script_path+"/CONFIG/")
sys.path.append(script_path+"/MAT/")
sys.path.append(script_path+"/HGNC/")

from CONFIG import *
from MATRICES import *
from SEQ import *
from HGNC import *

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
        elif opt=='-a':
            xml_file = arg
        elif opt=='-t':
            anno_type = arg

    return xml_file,anno_type

def usage(err=[]):

    if err:
        print '\n\n'
        print err
        print 'python processBLAST.py -a <analysis-xml-file> -t SEQ\n\n'
    else:
        print '\nPlease enter a valid command: \n'
        print 'python processBLAST.py -a <analysis-xml-file> -t SEQ\n\n'
    
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
        print keys,'\t',ontDict[keys]

    pathAnno = os.path.abspath(configDict['general']['dataDir'])+'/'+\
                                            configDict['version']+'/'+anno_type

    ontDict = configDict[anno_type]
    seqOutBlast = pathAnno+'/'+ontDict['blastOut']
    seqOutBit = pathAnno+'/'+ontDict['blastOutBit']
    seqOutNorm = pathAnno+'/'+ontDict['blastOutNorm']
    matScript = script_path+'/'+ontDict['rscript']
    matOut = pathAnno+'/'+ontDict['matOutFile']
    
    print "\nCheck if all the input, output files are correctly specified as "+\
                                            "present in the Analysis.xml file\n"
    print "BLAST results file: ",seqOutBlast
    print "BLAST results parsed bit score output: ",seqOutBit
    print "BLAST results parsed bit score normalized output: ",seqOutNorm
    print "BLAST results Matrix script file: ",matScript
    print "BLAST results Matrix out file: ",matOut
    print "\n"

    ontMeth = SEQMethod()
    objMat = MATRICES()
    objH = HGNC()

    ensg_genes_db = ontMeth.getDBUniqueGeneNames(configDict,objH)
    print len(ensg_genes_db)
    print 'Processing Raw File to generate QUERY TARGET Normalized Bit score: '+\
                                                                        anno_type

    #ontMeth.getQueryTargetBitScore(seqOutBlast,seqOutBit,ensg_genes_db)
    #ontMeth.getQueryTargetNormBitScore(seqOutBit,seqOutNorm,ensg_genes_db)

    print '\nGenerating PBS script for running BLAST Sequence similarity\n'
    pathAnno = os.path.dirname(os.path.dirname(matOut))
    pbs_file = obj.writePBSScriptFile(pathAnno,anno_type,configDict)
    print '\nAll Done\n'





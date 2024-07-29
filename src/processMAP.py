#!/usr/bin/python

import re,sys,os

script_path = os.path.dirname(os.path.abspath(__file__))

sys.path.append(script_path+'/CONFIG/')
sys.path.append(script_path+'/UTIL/')
sys.path.append(script_path+'/UNIPROT/')
sys.path.append(script_path+'/HGNC/')

from CONFIG import *
from UTIL import *
from UNIPROT import *
from HGNC import *

if __name__=='__main__':

    objC = CONFIG()
    objT = UTIL()
    objH = HGNC()
    objU = UNIPROT()

    # Get command line arguments
    cmdDict = objT.getMAPargs()
    xml_file = cmdDict['config']
    anno_type = cmdDict['anno']
    launch_type = objT.str2bool(cmdDict['launch'])
    
    # Get components of CONFG-XML file for given anno-type
    configDict = objC.getConfigDict(xml_file)
    pathAnno = configDict['general']['dataDir']+'/'+configDict['version']+\
               '/'+anno_type

    inp_hugo_map_file = os.path.dirname(pathAnno)+'/'+configDict[anno_type]['i']
    out_hgnc_ensg_file = pathAnno+'/'+configDict[anno_type]['o']

    uniprot_id_map_file = os.path.dirname(pathAnno)+'/UNIPROT/'+\
                          configDict['UNIPROT']['idMap']
    ppi_hgnc_file = os.path.dirname(pathAnno)+'/PPI/'+configDict['PPI']['matRow']
    ppi_ensg_map_file = pathAnno+'/'+configDict[anno_type]['p']
  
    #Get all the HGNC ids inbuilt in internal pBRIT database
    pbrit_hgnc_list = objT.getAllHgncPbritDb(configDict,os.path.dirname(pathAnno))
    
    #HUGO-approved-Ensembl Mapping
    objH.processRawHGNCMapFile(inp_hugo_map_file,out_hgnc_ensg_file,pbrit_hgnc_list)
    
    #Map Uniprot-accession Id of PPI to ENSEMBL
    objU.mapPPIUniprotEnsembl(uniprot_id_map_file,ppi_ensg_map_file,ppi_hgnc_file)

    #Write PBS script for launching MAP job
    pbs_file = objC.writePBSScriptFile(pathAnno,anno_type,configDict)
    if launch_type:
        os.system('qsub '+pbs_file)
    

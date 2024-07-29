#!/usr/bin/python

import re,sys,os

script_path = os.path.dirname(os.path.abspath(__file__))

sys.path.append(script_path+'/CONFIG/')
sys.path.append(script_path+'/UTIL/')

from CONFIG import *
from UTIL import *

if __name__=='__main__':

    objC = CONFIG()
    objT = UTIL()

    # Get command line arguments
    cmdDict = objT.getCMDargs()
    xml_file = cmdDict['config']
    anno_type = cmdDict['anno']
    launch_type = objT.str2bool(cmdDict['launch'])
    
    # Get components of CONFG-XML file for given anno-type
    configDict = objC.getConfigDict(xml_file)
    r_version = configDict['general']['tools']['R']
    svd_script = configDict['general']['tools']['svd']
    data_dir = configDict['general']['dataDir']
    pathAnno = configDict['general']['dataDir']+'/'+configDict['version']+\
               '/MAP/'

    #Write PBS script for launching MAP job
    pbs_file = objC.writePBSScriptFile(pathAnno,anno_type,configDict)
    if launch_type:
        os.system('qsub '+pbs_file)
     

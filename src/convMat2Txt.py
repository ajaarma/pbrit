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
    task_type = cmdDict['anno']
    launch_type = objT.str2bool(cmdDict['launch'])
    
    # Get components of CONFG-XML file for given anno-type
    configDict = objC.getConfigDict(xml_file)
    pathAnno = '/'.join([configDict['general']['dataDir'],configDict['version'],
                         'MAP'])

    # Path directories
    matrix_dir = '/'.join([pathAnno,'MATRICES'])
    tfidf_dir = '/'.join([pathAnno,'MATRICES/TFIDF'])
    svd_dir = '/'.join([pathAnno,'MATRICES/SVD'])
    tfidf_comp_dir = '/'.join([pathAnno,'MATRICES/COMP/TFIDF'])
    svd_comp_dir = '/'.join([pathAnno,'MATRICES/COMP/SVD'])

    # PBS script for Composite Matrices for TFIDF & SVD
    anno_source_prefix = ['PB','GO','HP','DO','MP','GD','HD','PY','PP','BL']

    pbs_file = objC.writePBSScriptFile(pathAnno,task_type,configDict,'SVD')
    pbs_file = objC.writePBSScriptFile(pathAnno,task_type,configDict,'TFIDF')
    pbs_file = objC.writePBSScriptFile(pathAnno,task_type,configDict,'MAT_X_TFIDF')
    pbs_file = objC.writePBSScriptFile(pathAnno,task_type,configDict,'MAT_Y_TFIDF')
    pbs_file = objC.writePBSScriptFile(pathAnno,task_type,configDict,'MAT_X_SVD')
    pbs_file = objC.writePBSScriptFile(pathAnno,task_type,configDict,'MAT_Y_SVD')

    for a_p in anno_source_prefix:
        pbs_file = objC.writePBSScriptFile(pathAnno,task_type,configDict,
                                                                'MAT_'+a_p+'_ANNO')


    if launch_type:
        os.system('qsub '+pbs_file)
    

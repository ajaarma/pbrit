#!/usr/bin/python

########################################################################
#
# Description: Generic class for creating XML based Configuration files; Specific functionalities to create directories, specifying path to input files depending upon the annotation sources.
#
#
#
#
########################################################################

from xml.etree import ElementTree
from xml.dom import minidom
from xml.etree.ElementTree import Element, SubElement, Comment, tostring
import xmltodict
import os
import re

class CONFIG:

    def __init__(self,elements=[]):
        self.__elements={}
        for e in elements:
            self.__elements[e]=1


    def display(self):
        print "Inside CONFIG class"

    def prettify(elem):
       """Return a pretty-printed XML string for the Element. """
       rough_string = ElementTree.tostring(elem, 'utf-8')
       reparsed = minidom.parseString(rough_string)
       return reparsed.toprettyxml(indent="  ")
    

    def getConfigTop(self, version_type):

        top = Element('config')
        comment = Comment('Configuration file for processing pBRIT annotation sources')
        top.append(comment)
        version = SubElement(top,"version")
        version.text = version_type

        return(top)
        

    def genConfigOnt(self,top, anno_type,data_dir,raw_anno_file,out_pre_file, 
                     out_gont_file, out_gNont_file, out_row_file, out_col_file, 
                                                        out_spr_file, obo_file):
        
        #top = getConfig()
        
        anno = SubElement(top,'annotation')
        anno.text = anno_type

        data = SubElement(anno,'data')
        raw_anno = SubElement(data,'RawAnnoFile')
        raw_anno.text = raw_anno_file
        out_pre   = SubElement(data,'preProcess')
        out_pre.text = out_pre_file
        out_gont = SubElement(data,'geneOnt')
        out_gont.text = out_gont_file
        out_gNont = SubElement(data,'NoOnt')
        out_gNont.text = out_gNont_file
        mat = SubElement(anno,'matrix')
        out_row = SubElement(mat,'matRow')
        out_row.text = out_row_file
        out_col = SubElement(mat,'matCol')
        out_col.text = out_col_file
        out_spr = SubElement(mat,'matSpr')
        out_spr.text = out_spr_file
         
        obo = SubElement(top,'OboFile')
        obo.text = obo_file

        return(top)

    def getConfigDict(self,analysis_file,module_name=[]):

        with open(analysis_file) as fd:
            doc = xmltodict.parse(fd.read())
        
        #print "The item classes present in the Analysis Configuration Files are: "+",".join(doc.keys())
        doc = doc["config"]
        if module_name:
            config_dict = doc[module_name]
        else:
            config_dict = doc


        return config_dict

    def writePBSScriptFile(self,data_dir,anno_type,configDict,method=[]):
        ''' Subroutine to generate PBS script for processing TFIDF-Score 
        and BLAST Sequence similarity matrix '''
    
        if anno_type=='SEQ':
            #pbs_file = os.path.dirname(os.path.dirname(data_dir))+\
            pbs_file = os.path.abspath(data_dir)+"/QSUB"+"/SIM_Score"+anno_type+".sh"
        
        elif anno_type=='MAP':
            pbs_file = os.path.abspath(data_dir)+"/QSUB"+"/MAP_Score"+".sh"
        
        elif anno_type=='COMP':
            if method in ['TFIDF','SVD']:
                pbs_file = os.path.abspath(data_dir)+"/QSUB"+"/compositeMatrix_"+method+".sh"
            elif re.search('MAT',method):
                pbs_file = os.path.abspath(data_dir)+"/QSUB"+"/convert_"+method+".sh"
        else:
            pbs_file = os.path.abspath(data_dir)+"/QSUB"+"/TFIDF_Score_"+\
                                                                anno_type+".sh"

        wh = open(pbs_file,"w")
        
        qsubDict = configDict["qsub"]
        if configDict.has_key(anno_type):
            annoDict = configDict[anno_type]

        dash = "-"
       
        anno_type_orig = anno_type

        print >>wh,"#!/bin/bash"
        print >>wh,"\n"
        for key in qsubDict:
            if key==dash:
                pass
            elif key=="N":
                print >>wh,"#PBS "+dash+key+" "+qsubDict[key]+"_"+anno_type
            elif key=="d":
                print >>wh,"#PBS "+dash+key+" "+data_dir+"/"
            elif key=="o":

                if anno_type=='COMP':
                    anno_type = method
                else:
                    anno_type = anno_type

                print >>wh,"#PBS "+dash+key+" "+data_dir+"/QOUT/"+anno_type+\
                                                                    "_qsub_o.txt"
            elif key=="e":
                if anno_type=='COMP':
                    anno_type = method
                else:
                    anno_type = anno_type
                
                print >>wh,"#PBS "+dash+key+" "+data_dir+"/QOUT/"+anno_type+\
                                                                    "_qsub_e.txt"
            elif key=="l":
                l_list = qsubDict["l"]["value"]
                if anno_type=='SEQ':
                    l_list[2]='mem=80gb'
                print >>wh,"#PBS "+dash+key+" "+l_list[0]+":"+l_list[1]+","+\
                                                                        l_list[2] 
            else:
                print >>wh,"#PBS "+dash+key+" "+qsubDict[key]

        print >>wh,"\n\n\n"

        script_path = os.path.dirname(os.path.dirname(
                                            os.path.abspath(__file__))
                                     )
        
        anno_type = anno_type_orig
        # Commands for each different type of Anno sources       
        if anno_type=='SEQ':
            cmd = configDict['general']['tools']['R']+' '+script_path+'/'+\
                  annoDict['rscript']+' -i '+\
                  data_dir+'/'+annoDict['blastOutNorm']+' -o '+data_dir+' -m '+\
                                                data_dir+'/'+annoDict['matOutFile']
            print >>wh,cmd
            wh.close()
        
        elif anno_type=='MAP':
            
            rscript = configDict['general']['tools']['R']
            rscript_map = script_path+'/'+configDict[anno_type]['rscript']
            hgnc_ensg_map_file = configDict['general']['dataDir']+'/'+\
                       configDict['version']+'/MAP/'+configDict[anno_type]['o']
            
            ppi_ensg_map_file = configDict['general']['dataDir']+'/'+\
                       configDict['version']+'/MAP/'+configDict[anno_type]['p']
           
            data_dir = os.path.abspath(configDict['general']['dataDir']+\
                                            '/'+configDict['version'])+'/'            
            
            cmd = ' '.join([rscript,rscript_map,'-i',hgnc_ensg_map_file,
                                                '-p',ppi_ensg_map_file,
                                                '-o',data_dir
                           ])
            
            print >>wh,cmd
            wh.close()

        elif anno_type == 'COMP':
            rscript = configDict['general']['tools']['R']
            compBin = '/'.join([script_path,configDict['COMP']['compBin'] ])
            convBin = '/'.join([script_path,configDict['COMP']['convBin'] ])
            mat_dir = '/'.join([configDict['general']['dataDir'],
                                            configDict['version'],
                                                    'MAP/MATRICES'
                               ])
            pathAnno = '/'.join([configDict['general']['dataDir'],
                                    configDict['version'], 'MAP'
                                 ])
  
            if method in ['TFIDF','SVD']:

                cmd_comp = ' '.join([rscript,compBin,'-d',mat_dir,'-o',mat_dir,
                                                                    '-m',method
                               ])

                print >>wh,cmd_comp
                #print >>wh,cmd_conv

                wh.close()

            elif re.search('MAT',method):
                
                if re.search('X_TFIDF',method):
                    matComp = '/'.join([mat_dir,'COMP/TFIDF/composite_X'])
                    tmpOut = '/tmp/tfidf_composite_X.txt'
                    matOut = '/'.join([mat_dir,'COMP/TFIDF/composite_X.txt.gz'])
                    
                elif re.search('Y_TFIDF',method):
                    matComp = '/'.join([mat_dir,'COMP/TFIDF/composite_Y'])
                    tmpOut = '/tmp/tfidf_composite_Y.txt'
                    matOut = '/'.join([mat_dir,'COMP/TFIDF/composite_Y.txt.gz'])

                elif re.search('X_SVD',method):
                    matComp = '/'.join([mat_dir,'COMP/SVD/composite_X'])
                    tmpOut = '/tmp/svd_composite_X.txt'
                    matOut = '/'.join([mat_dir,'COMP/SVD/composite_X.txt.gz'])

                elif re.search('Y_SVD',method):
                    matComp = '/'.join([mat_dir,'COMP/SVD/composite_Y'])
                    tmpOut = '/tmp/svd_composite_Y.txt'
                    matOut = '/'.join([mat_dir,'COMP/SVD/composite_Y.txt.gz'])
                else:
                    an_strs = re.split('\_',method)
                    if an_strs[1] !='BL':
                        
                        tfidf_dir = '/'.join([pathAnno,'MATRICES/TFIDF'])
                        matComp = '/'.join([tfidf_dir,'ensg.'+an_strs[1]+'.TFIDF'])
                        tmpOut = '/'.join(['/tmp','ensg.'+an_strs[1]+'.TFIDF.txt'])
                        matOut = '/'.join([tfidf_dir,'ensg.'+an_strs[1]+'.TFIDF.txt.gz'])
                    else:
                        pathSeq = '/'.join([os.path.dirname(pathAnno),'SEQ'])
                        matComp = '/'.join([pathSeq,'MATRICES/BL.ENSG.sim.mat'])
                        tmpOut = '/'.join(['/tmp','BL.ENSG.sim.mat.txt'])
                        matOut = '/'.join([pathSeq,'MATRICES/BL.ENSG.sim.mat.txt.gz'])
            
                cmd_conv = ' '.join([rscript,convBin,'-d',matComp,'-o',tmpOut,
                                                                  '-m',method
                               ])
                cmd_conv_mv = ' '.join(['mv ',tmpOut+'.gz',matOut])

                print >>wh,cmd_conv
                print >>wh,cmd_conv_mv
                wh.close()

        elif anno_type=='SVD':
            
            r_version = configDict['general']['tools']['R']
            svd_script = script_path+'/'+configDict['general']['tools']['svd']
            anno_list = configDict['general']['annoSource']
            data_dir = configDict['general']['dataDir']+'/'+configDict['version']

            tfidf_mat_anno = ['HP','DO','HD','GD','GO','MP','PB','PY','PP']

            for anno in tfidf_mat_anno:
                print >>wh,'echo \"Computing SVD for: '+anno+'\"'
                if anno != 'PB':
                    tfidf_mat_path = data_dir+'/MAP/MATRICES/TFIDF/ensg.'+anno+'.TFIDF'
                elif anno == 'PB':
                    tfidf_mat_path = data_dir+'/PUB/MATRICES/mat.PB.TFIDF'
                svd_mat_path = os.path.abspath(data_dir+'/MAP/MATRICES/SVD/')
                cmd_svd = [r_version,svd_script,'-d',tfidf_mat_path,'-a',anno,
                                                            '-o',svd_mat_path]
                print >>wh,' '.join(cmd_svd)+'\n'

        else:
            #anno_list = configDict[]
            print annoDict
            print anno_type
            cmd = configDict['general']['tools']['R']+' '+script_path+'/'+\
                  annoDict["rscript"]+" "+annoDict["matSparse"]+"  "+\
                  annoDict["matRow"]+"  "+annoDict["matCol"]+"  "+data_dir+\
                                                    "/MATRICES"+" "+anno_type
            print >>wh,cmd
            wh.close()

        return pbs_file


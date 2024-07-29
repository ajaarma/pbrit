#!/usr/bin/python

import sys
from xml.etree import ElementTree
from xml.etree.ElementTree import Element, SubElement, Comment
#from ElementTree_pretty import prettify
from xml.dom import minidom
sys.path.append("/home/aakumar/pbrit_scripts/CONFIG/")
from CONFIG import *

def prettify(elem):
    """Return a pretty-printed XML string for the Element. """
    rough_string = ElementTree.tostring(elem, 'utf-8')
    reparsed = minidom.parseString(rough_string)
    return reparsed.toprettyxml(indent="  ")



if __name__=="__main__":


    data_dir = "/home/shared_medgen_aorta/pbrit_data/"
    version  = "V3"
    #raw_anno_file_list = ["HPO_process.txt","GO_process.txt","DO_process.txt","MPO_process.txt"]
    
    objC = CONFIG()
    top = objC.getConfigTop("V3")

    anno_list = ["HPO", "GO","DO","MPO"]
    for ele in anno_list:
        
        raw_anno_file = ele+"_raw.txt"
        out_pre_file  = ele+"_preprocess.txt"
        out_gont_file = "gene_"+ele+".txt"
        out_gNont_file = "gene_no_"+ele+".txt"
        out_row_file = "gene_"+ele+"_row.txt"
        out_col_file = "gene_"+ele+"_col.txt"
        out_spr_file = "gene_"+ele+"_spr.txt"
        obo_file = "OBO_"+ele+".obo"
        
        top = objC.genConfigOnt(top, ele,data_dir,raw_anno_file, out_pre_file, out_gont_file, out_gNont_file, out_row_file, out_col_file, out_spr_file, obo_file)


    print prettify(top)

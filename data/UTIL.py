#!/us/bin/python

import re,sys,os
import argparse

class UTIL:

    def __init__(self,elements=[]):
        self.__elements={}
        for e in elements:
            self.__elements[e]=1

    def str2bool(self,v):
        return v.lower() in ('yes','true','t','1')


    def getCMDargs(self):
        ''' Processing Command line arguments for Mapping HGNC-ENSMBL '''

        cmdDict = {}
        parser = argparse.ArgumentParser(
                   description='python <process***.py> -a <xml-file> -p <proj-name>\n'
                 )

        parser.add_argument('-a','--analXML',help='Analysis XML file: Analysis.xml',
                             action='store',dest='xmlFile',required=True
                           )
        parser.add_argument('-t','--annoType',help='Annotation Source',
                             action='store',dest='annoType',required=True
                           )
        parser.add_argument('-l','--launch',help='Boolean:Submit PBS script',
                             action='store',dest='launch',required=True
                           )
        self.results = parser.parse_args()
        self.xml_file = os.path.abspath(self.results.xmlFile)
        self.anno_type = self.results.annoType
        self.launch_job = self.results.launch

        cmdDict = {'config':self.xml_file,'anno':self.anno_type,
                   'launch':self.launch_job
                  }

        return cmdDict

    def getAllHgncPbritDb(self,configDict,pathAnno):
        
        ''' Function to extract all the HGNC gene ids present in the internal 
        database of pBRIT. The datbase all excluding: Abstract (mapped to Ensembl), 
        Sequence (Mapped to Ensembl), PPI (mapped to Uniprot-ID, MAP,HUGO '''

        anno_list = configDict['general']['annoSource']['value']
        
        hgnc_dict = {}

        for anno_type in anno_list:
            if anno_type in ['MAP','UNIPROT','PUB','SEQ','HUGO','PPI']:
                pass
            else:
                hgnc_file = pathAnno+'/'+anno_type+'/INPUT/HGNC_rownames.txt'
                fh = open(hgnc_file)

                for lines in fh:
                    lines = lines.strip()
                    strs = re.split('\t',lines);strs = [x.strip().upper() for x in strs]
                    hgnc_id = strs[1]
                    hgnc_dict[hgnc_id] = 1

                fh.close()

        return hgnc_dict.keys()

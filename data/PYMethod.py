#!/usr/bin/python

import re,os
from collections import OrderedDict


class PYMethod:

    def __init__(self, elements=[]):
        self.__elements={}
        for e in elements:
            self.__elements[e]=1


    def processRawPYFile(self,annoFile,genePY):

        py_hash = OrderedDict()

        fh = open(annoFile)
        wh = open(genePY,"w")

        for lines in fh:
            lines = lines.strip()
            if not re.search("external|hgnc|source",lines):
                strs = re.split("\t",lines)
                #print strs
                py_term = strs[0].strip()
                gene_list = re.split("\,",strs[3].strip())
                
                for gene_id in gene_list:
                    gene_id = gene_id.strip().upper()
                    if py_hash.has_key(gene_id):
                        tmp = py_hash[gene_id]
                        tmp.append(py_term)
                        py_hash[gene_id] = tmp
                        tmp = []

                    else:
                        py_hash[gene_id] = [py_term]
                

        for keys in py_hash:
            print >>wh,keys+"\t"+",".join(py_hash[keys])

        wh.close()

        


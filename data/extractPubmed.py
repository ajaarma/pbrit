#!/usr/bin/python

##################################################
# Title: Program to download abstract from PUBMED
#
# Usage: python extractPubmed.py  Inputfile.txt
#
# Author: Ajay Anand Kumar ajay.kumar@ua.ac.be
#
#################################################




import sys
sys.path.append("/home/aakumar/aakumar/WORK/PUBMED_TEXT/")
from Bio import Entrez,Medline
import re



def getAbstract(pmid):
   try:
      handle = Entrez.efetch(db="pubmed",id=pmid,rettype="medline",retmode="text")
      records = Medline.parse(handle)
      records = list(records)
      #print records
      abst_list = []
      title = records[0]['TI']
      abst_list.append(title)
      abst = records[0]['AB']
      abst_list.append(abst)
   except:
      pass
   
   return abst_list

if __name__ == "__main__":

	Entrez.email = "aakumar1707@gmail.com"

	data_dir = "/home/shared_data_medgen_aorta/pbrit_data/PUB/"
	
	inp_file = data_dir+"PUB_UNIQUE.txt"
	out_file = data_dir+"PUB_ABSTRACT.txt"
	no_abst = data_dir+"NO_PUBMED_ABSTRACT.txt"
	fh = open(inp_file)
   
   #gene_hash = {}
	op = open(out_file,"w")
	wh = open(no_abst,"w")
   
	for lines in fh:
		lines = lines.strip()
      		strs = re.split("\t",lines)
      		strs1 = re.split("\,",strs[0])
      		abst_list = []
      		print strs[0]
      		for ele in strs1:
         		ele = ele.strip()
         		abst = getAbstract(ele)
         		abst_list = abst
      			if abst_list:
         			print >>op,strs[0],"\t"," ".join(abst_list)
      			else:
	
	 			print >>wh,strs[0].strip()
	
  	wh.close()
   	op.close()

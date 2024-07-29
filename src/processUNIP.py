#!/usr/bin/python

import re,sys
sys.path.append("UNIPROT/")
sys.path.append("HGNC/")
from UNIPROT import *
from HGNC import *

def execUNIPROT(inp_file,out_file, ensg_hgnc, ensg_go_pr,out_pub):

	obj = UNIPROT()
	obj.display()
	col_ind = [19,7]
	obj.parseUNIP(inp_file, out_file, col_ind,ensg_hgnc)
	obj.processENSG_GO(out_file, ensg_go_pr)
	obj.processENSG_PUB(inp_file,out_pub, ensg_hgnc)


def execHGNC(inp_file, out_file):
	
	obj = HGNC()
	obj.display()
	ensg_hgnc = obj.parseHGNC(inp_file, out_file)
	
	return ensg_hgnc


if __name__=="__main__":

	data_dir = "/home/shared_data_medgen_aorta/pbrit_data/"
	ensg_hugo = data_dir+"ENSG_HUGO/Ensg_Hugo.txt"
	ensg_hgnc_map = data_dir+"ENSG_HUGO/HGNC_ENSG_MAP.txt"

	id_map_file = data_dir+"UNIPROT/HUMAN_9606_idmapping_selected.tab"
	ensg_go = data_dir+"GO/ENSG_GO_preprocess.txt"
	ensg_go_pr = data_dir+"GO/ENSG_GO_processed.txt"

	out_ensg_pub = data_dir+"PUB/ENSG_PUB_tmp.txt"
	
		
	ensg_hgnc = execHGNC(ensg_hugo,ensg_hgnc_map) #inp_file: Ensg HGNC map file obtained from biomart; out_file: Name of the output file;
	#for keys in ensg_hgnc:
	#	print keys, "\t", "||".join(ensg_hgnc[keys])
	execUNIPROT(id_map_file, ensg_go, ensg_hgnc,ensg_go_pr, out_ensg_pub) #inp_file: 
	#execUNIPROT(ensg_go,ensg_go_pr)
	

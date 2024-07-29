#!/usr/bin/python

import re,sys

sys.path.append("/home/aakumar/pbrit_scripts/PUB/")
from ABSTRACT import *


if __name__=="__main__":

	data_dir = "/home/shared_data_medgen_aorta/pbrit_data/PUB/"
	ensg_pub = "ENSG_PUB.txt"
	pub_abst = "PUB_ABSTRACT.txt"
	out_file = "ENSG_PUBMED_ABSTRACT.txt"
	obj = ABSTRACT()
	
	pub_hash = obj.getPubHash(data_dir+pub_abst)
	ensg_pub_hash = obj.getEnsgPubHash(data_dir+ensg_pub)
	obj.mapEnsgAbstract(ensg_pub_hash,pub_hash,data_dir+out_file)
	
	

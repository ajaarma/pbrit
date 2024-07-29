#!/usr/bin/python

######################################
#
# Usage: python extractingIDS.py /home/shared_data_medgen_aorta/pbrit_data/SEQ/Ensg_Biomart_fasta_HUMAN.txt /home/shared_data_medgen_aorta/pbrit_data/ENSG_HUGO/HGNC_ENSG_MAP.txt > /home/shared_data_medgen_aorta/pbrit_data/SEQ/ENSG_MAP_PROT_SEQ.txt 
#
#
#
#
######################################


import re
import sys
file = sys.argv[1]

fh = open(file)
fh1 = open(sys.argv[2])

ensg_hash = {}

for lines in fh1:
	lines = lines.strip()
	ensg_id = re.split("\t",lines)[1].strip()
	ensg_hash[ensg_id] = 1

fh1.close()

prot = []
gene_id = []
data_hash = {}
data_hash2= {}
count = 0
tmp = []


for lines in fh:
	lines = lines.strip()
	
	if not re.search(">",lines):
		prot.append(lines)

	else:
		tmp = ''.join(prot)
		strs = re.split("\>|\|",lines)
		data_hash2[strs[1].strip()] = 1
		ensg_id = lines#strs[0].strip()
		gene_id.append(ensg_id)
		if count >=1:
			data_hash[gene_id[count-1]] = tmp
		#else:
		#	data_hash[gene_id] = tmp
		prot = []
		count = count+1

tmp = ''.join(prot)
data_hash[gene_id[count-1]] = tmp


for keys in data_hash:
	if not re.search("Sequence",data_hash[keys]):
		strs = re.split("\>|\|",keys)
		#print strs
		#sys.exit()
		if ensg_hash.has_key(strs[1].strip()):
			#print ">"+strs[1]
			print keys
		#print keys
			print data_hash[keys]
			#pass
		#else:
		#	print keys

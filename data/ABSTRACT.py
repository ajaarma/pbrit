#!/usr/bin/python

import re,sys

class ABSTRACT:

    def __init__(self, elements=[]):
        self.__elements={}
        for e in elements:
            self.__elements[e] = 1

	def getPubHash(self, abst_file):
		ab = open(abst_file)
        abs_hash = {}

        for lines in ab:
            lines = lines.strip()
            strs = re.split("\t",lines)
            abs_hash[strs[0].strip()] = strs[1]

        ab.close()

		return abs_hash

	def getEnsgPubHash(self, ensg_pub_file):
		
		ep = open(ensg_pub_file)
		ensg_hash = {}
		
		for lines in ep:
		    lines = lines.strip()
			strs = re.split("\t",lines)
			try:
				ensg_id = strs[0].strip()
				pub_list = re.split("\;",strs[1].strip())
				if ensg_hash.has_key(ensg_id):
					tmp = ensg_hash[ensg_id]
					tmp = tmp+pub_list
					ensg_hash[ensg_id] = list(set(tmp))
					tmp = []

				else:
					ensg_hash[ensg_id] = pub_list
			except:
				pass
			
		ep.close()
		
		return ensg_hash

	def mapEnsgAbstract(self, ensg_hash, pub_hash,out_file):


		wh = open(out_file,"w")
		
		for ensg_id in ensg_hash:

		    pub_strs = ensg_hash[ensg_id]
			abst_list = []
			for e in pub_strs:
				e = e.strip()
				if pub_hash.has_key(e):
					abst_list.append(pub_hash[e])

			print >>wh,ensg_id+"\t"+''.join(abst_list)

		wh.close()


#!/usr/bin/python

import re

##############################################################################################
#
# Description: Generic class for processing Gene Ontology DAG structre. Generating parent/ancestor terms for the GO classes
#
#
# Methods: GOMethod: Method for retrieving parent/ancestor terms for any given GO Classes
#
#
#
##############################################################################################

class GOMethod:
    parents = []
    sorted_parents={}
    alt_id={}
    no_go={}
    altid_flag_val='t'

    def __init__(self, elements=[]):
        self.__elements={}
        for e in elements:
            self.__elements[e]=1


    def getGOParents(self, GO_hash, go_id):

        self.GO_hash = GO_hash
        self.go_id = go_id

        if self.go_id=='0008150' or self.go_id == '0003674' or self.go_id == '0005575':
            if not go_id in GOMethod.parents:
                GOMethod.parents.append(self.go_id)

        else:
            if self.GO_hash.has_key(self.go_id):
                keys = self.go_id
                if self.GO_hash[keys]['is_a']:
                    isA_go = self.GO_hash[keys]['is_a']
                    for e in isA_go:
                        GOMethod.parents.append(e)
                        self.getGOParents(self.GO_hash,e)
                
                if self.GO_hash[keys]['part_of']:
                    isA_part = self.GO_hash[keys]['part_of']
                    for e in isA_part:
                        GOMethod.parents.append(e)
                        self.getGOParents(self.GO_hash,e)

                if self.GO_hash[keys]['consider']:
                    isA_consider = self.GO_hash[keys]['consider']
                    for e in isA_consider:
                        GOMethod.parents.append(e)
                        self.getGOParents(self.GO_hash,e)
            else:
                print "Key not founds, No Parents available for GO: "+self.go_id

                if GOMethod.altid_flag_val == 't':
                    print "Looking for the alternate ids.............."
                    flag=0

                    for keys in self.GO_hash:
                        alt_id_tmp = self.GO_hash[keys]['alt_id']
                        for ele in alt_id_tmp:
                            if ele == self.go_id:
                                print "Alternate id for GO:"+self.go_id,"is GO:"+keys,"\n"
                                flag=1
                                GOMethod.alt_id[self.go_id] = keys
                                if self.GO_hash[keys]['is_a']:
                                    isA_go = self.GO_hash[keys]['is_a']
                                    for e in isA_go:
                                        GOMethod.parents.append(e)
                                        self.getGOParents(self.GO_hash,e)

                                if self.GO_hash[keys]['part_of']:
                                    isA_part = self.GO_hash[keys]['part_of']
                                    for e in isA_part:
                                        GOMethod.parents.append(e)
                                        self.getGOParents(self.GO_hash,e)

                                if self.GO_hash[keys]['consider']:
                                    isA_consider = self.GO_hash[keys]['part_of']
                                    for e in isA_consider:
                                        GOMethod.parents.append(e)
                                        self.getGOParents(self.GO_hash,e)
                                        
                                if self.GO_hash[keys]['replaced_by']:
                                    isA_replaced = self.GO_hash[keys]['replaced_by']
                                    for e in isA_replaced:
                                        GOMethod.parents.append(e)
                                        self.getGOParents(self.GO_hash,e)

                    if flag==0:
                        GOMethod.no_go[self.go_id]=1
                        print "No alternate id found for GO:"+self.go_id+"\n"

        for outer in GOMethod.parents:
            GOMethod.sorted_parents[outer]=1


        return GOMethod.parents



    def getFile_GO(self, GO_hash, inp_file, out_gene_go, out_no_go, altid_val):

        self.go_lib = {}
        GOMethod.altid_flag_val = altid_val

        op = open(out_gene_go,'w')
        ng = open(out_no_go,'w')

        gene_id_list = []
        self.gene_go_hash = {}
        parents_list = []
        no_go = []

        fh = open(inp_file)

        lines = fh.readline()

        while lines:
            lines = lines.strip()
            exp = re.split('\t',lines)
            self.comb_parents={}
            if exp:
                gene_id = exp[0].strip()
                gene_id_list.append(gene_id)
                gene_go_strs = re.split("\s|\,|\;",exp[1])
                gene_go_strs = [x for x in gene_go_strs if x]
            
                go_list = []
                for ele_go in gene_go_strs:
                    ele_go = ele_go.strip()
                    if ele_go:
                        exp_strs = re.split("\:",ele_go)
                        go_num = exp_strs[1]
                        go_list.append(go_num)
                        if self.go_lib.has_key(go_num):
                            parents_list = self.go_lib[go_num]
                            for ele in parents_list:
                                self.comb_parents[ele] = 1
                            parents_list = []
                        else:
                            if not GOMethod.no_go.has_key(go_num):
                                try:
                                    parents_list = self.getGOParents(GO_hash,go_num)
                                    if parents_list:
                                        self.go_lib[go_num] = parents_list
                                        for ele in parents_list:
                                            self.comb_parents[ele] = 1
                                    else:
                                        no_go.append(go_num)
                                except IndexError:
                                    print "Error occurred:",go_num
                                GOMethod.parents = []
                                GOMethod.sorted_parents = {}
                                parents_list = []

                if self.comb_parents:
                    print >>op, gene_id+"\t"+",".join(go_list)+"\t"+','.join(self.comb_parents.keys())
                else:
                    print >>op,gene_id+"\t"+','.join(go_list)

            lines = fh.readline()

        if GOMethod.no_go:
            for keys in GOMethod.no_go:
                print >>ng,keys
        self.comb_parents = {}
        GOMethod.no_go = {}
        op.close()
        ng.close()
        fh.close()

    def processRawGOFile(self,raw_file, geneGO):
        
        gene_hash = {}
        go_id = []

        fh = open(raw_file)

        for lines in fh:
            lines = lines.strip()

            if re.search("\!",lines):
                pass
            else:
                strs = re.split("\t",lines)
                #print strs
                go_id = []

                try:
                    gene_id = strs[2].strip().upper()
                    go_id = strs[4].strip()

                except IndexError:
                    gene_id = []
                    go_id = []

                if gene_hash.has_key(gene_id):
                    tmp = gene_hash[gene_id]
                    tmp.append(go_id)
                    gene_hash[gene_id] = list(set(tmp))
                    tmp = []
                else:
                    if gene_id !='' and go_id !='':
                        gene_hash[gene_id] = [go_id]

        fh.close()
        wh = open(geneGO,"w")

        for keys in gene_hash:
            print >>wh,keys+"\t"+",".join(gene_hash[keys])

        wh.close()


#!/usr/bin/python

import re,sys
from collections import OrderedDict

class DOMethod:

    parents =[]
    sorted_parents = OrderedDict()
    children = []#Added for incorporating getChildren function
    sorted_children = OrderedDict()
    alt_id = OrderedDict()
    no_do = OrderedDict()
    altid_flag_val = 't'

    def getDOParents(self, DO_hash, do_id):

        self.DO_hash = DO_hash
        self.do_id = do_id

        if self.do_id == '4':
            if not do_id in DOMethod.parents:
                DOMethod.parents.append(self.do_id)

        else:
            if self.DO_hash.has_key(self.do_id):
                keys = self.do_id
                if self.DO_hash[keys]['is_a']:
                    isA_do = self.DO_hash[keys]['is_a']
                    for e in isA_do:
                        DOMethod.parents.append(e)
                        self.getDOParents(self.DO_hash,e)
                else:
                    if DOMethod.altid_flag_val =='t':
                        print "Looking for Alternate ids....."
                        flag=0
                        for keys in self.DO_hash:
                            alt_id_tmp = self.DO_hash[keys]['alt_id']
                            for ele in alt_id_tmp:
                                if ele==self.do_id:
                                    flag=1
                                    DOMethod.alt_id[self.do_id] = keys
                                    if self.DO_hash[keys]['is_a']:
                                        isA_do = self.DO_hash[keys]['is_a']
                                        for e in isA_do:
                                            self.getDOParents(self.DO_hash,e)
                        if flag==0:
                            DOMethod.no_do[self.do_id]=1
                            print "No alternate ID found... for DO:"+self.do_id

        for outer in DOMethod.parents:
            DOMethod.sorted_parents[outer]=1

        return DOMethod.parents


    def getFile_DO(self, DO_hash, inp_file, out_gene_do, out_no_do, altid_val):

        self.do_lib = OrderedDict()
        DOMethod.altid_flag_val = altid_val

        op = open(out_gene_do,'w')
        ng = open(out_no_do,'w')

        gene_id_list = []
        #self.gene_hpo_hash = {}
        parents_list = []
        no_do = []

        fh = open(inp_file)

        lines = fh.readline()

        while lines:
            lines = lines.strip()
            exp = re.split("\t",lines)
            self.comb_parents = OrderedDict()
            gene_do_list = []

            if exp:
                gene_id = exp[0].strip()
                gene_id_list.append(gene_id)
                gene_do_strs = re.split("\s|\,|\;",exp[1])
                gene_do_strs = [x for x in gene_do_strs if x]

                for ele_do in gene_do_strs:
                    ele_do = ele_do.strip()
                    if ele_do:
                        exp_strs = re.split("\:",ele_do)
                        do_num = exp_strs[1]
                        gene_do_list.append(do_num)

                        if self.do_lib.has_key(do_num):
                            parents_list = self.do_lib[do_num]
                            for ele in parents_list:
                                self.comb_parents[ele] = 1
                            parents_list = []
                        else:
                            if not DOMethod.no_do.has_key(do_num):
                                try:
                                    parents_list = self.getDOParents(DO_hash,do_num)
                                    if parents_list:
                                        self.do_lib[do_num] = parents_list
                                        for ele in parents_list:
                                            self.comb_parents[ele] = 1
                                    else:
                                        no_do.append(do_num)

                                except IndexError:
                                    print "Error occurred:",do_num

                                DOMethod.parents = []
                                DOMethod.sorted_parents = OrderedDict()
                                parents_list = []

                if self.comb_parents:
                    print >>op, gene_id+"\t"+",".join(gene_do_list)+"\t"+','.join(self.comb_parents.keys())
                else:
                    print "Gene Ids with HPO Ids having NO Parent HPO terms: " + lines
                    print >>op, gene_id+"\t"+",".join(gene_do_list)

            lines = fh.readline()

        if DOMethod.no_do:
            for keys in DOMethod.no_do:
                print >>ng,keys
        self.comb_parents = OrderedDict()
        DOMethod.no_do = OrderedDict()

        op.close()
        fh.close()
        ng.close()


    def processRawDOFile(self, inp_exp, inp_text, inp_knl, out_file):

        gene_hash = OrderedDict()

        try:
            fh1 = open(inp_exp)
            fh2 = open(inp_text)
            fh3 = open(inp_knl)

            wh = open(out_file,"w")

            for lines in fh1:
                lines = lines.strip()
                strs = re.split("\t",lines)
                gene_id = strs[1].strip().upper()
                do_id = strs[2].strip()
                if re.search("^DO",do_id):
                    tmp = []
                    if gene_hash.has_key(gene_id):
                        
                        #strs1 = re.split("\,",gene_hash[gene_id])
                        tmp = gene_hash[gene_id]
                        tmp.append(do_id)
                        #gene_hash[gene_id] = ",".join(list(set(do_id)))
                        gene_hash[gene_id] = list(set(tmp))
                        tmp = []
                    else:
                        gene_hash[gene_id] = [do_id]


            for lines in fh2:
                lines = lines.strip()
                strs = re.split("\t",lines)
                gene_id = strs[1].strip().upper()
                do_id = strs[2].strip()
                
                if re.search("^DO",do_id):
                    tmp = []
                    if gene_hash.has_key(gene_id):
                        tmp = gene_hash[gene_id]
                        tmp.append(do_id)
                        gene_hash[gene_id] = list(set(tmp))
                        tmp = []
                    else:
                        gene_hash[gene_id] = [do_id]

            for lines in fh3:
                lines = lines.strip()
                strs = re.split("\t",lines)
                gene_id = strs[1].strip().upper()
                do_id = strs[2].strip()

                if re.search("^DO",do_id):
                    tmp = []
                    if gene_hash.has_key(gene_id):
                        tmp = gene_hash[gene_id]
                        tmp.append(do_id)
                        gene_hash[gene_id] = list(set(tmp))
                        tmp = []
                    else:
                        gene_hash[gene_id] = [do_id]


            for keys in gene_hash:
                #if hgnc_ensg.has_key(keys):
                    #ensg_id = hgnc_ensg[keys][0]
                hgnc_id = keys
                print >>wh,hgnc_id+"\t"+",".join(gene_hash[keys])
            wh.close()
            fh1.close()
            fh2.close()
            fh3.close()

        except IOError:

            print "Not a valid path/file. Please enter correct path or file"
            sys.exit()
        #except:
        #    sys.exit()


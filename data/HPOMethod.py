#!/usr/bin/python

import re,sys
from collections import OrderedDict

class HPOMethod:

    parents =[]
    sorted_parents = OrderedDict()
    children = []#Added for incorporating getChildren function
    sorted_children = OrderedDict()
    alt_id = OrderedDict()
    no_hp = OrderedDict()
    altid_flag_val = 't'

    def getHPOParents(self, HP_hash, hp_id):

        self.HP_hash = HP_hash
        self.hp_id = hp_id

        if self.hp_id == '0000001':
            if not hp_id in HPOMethod.parents:
                HPOMethod.parents.append(self.hp_id)

        else:
            if self.HP_hash.has_key(self.hp_id):
                keys = self.hp_id
                if self.HP_hash[keys]['is_a']:
                    isA_hp = self.HP_hash[keys]['is_a']
                    for e in isA_hp:
                        HPOMethod.parents.append(e)
                        self.getHPOParents(self.HP_hash,e)
                else:
                    if HPOMethod.altid_flag_val =='t':
                        print "Looking for Alternate ids....."
                        flag=0
                        for keys in self.HP_hash:
                            alt_id_tmp = self.HP_hash[keys]['alt_id']
                            for ele in alt_id_tmp:
                                if ele==self.hp_id:
                                    flag=1
                                    HPOMethod.alt_id[self.hp_id] = keys
                                    if self.HP_hash[keys]['is_a']:
                                        isA_hp = self.HP_hash[keys]['is_a']
                                        for e in isA_hp:
                                            self.getHPOParents(self.HP_hash,e)
                        if flag==0:
                            HPOMethod.no_hp[self.hp_id]=1
                            print "No alternate ID found... for HP:"+self.hp_id

        for outer in HPOMethod.parents:
            HPOMethod.sorted_parents[outer]=1

        return HPOMethod.parents


    def getFile_HPO(self, HPO_hash, inp_file, out_gene_hpo, out_no_hpo, altid_val):

        self.hpo_lib = OrderedDict()
        HPOMethod.altid_flag_val = altid_val

        op = open(out_gene_hpo,'w')
        ng = open(out_no_hpo,'w')

        gene_id_list = []
        #self.gene_hpo_hash = {}
        parents_list = []
        no_hp = []

        fh = open(inp_file)

        lines = fh.readline()

        while lines:
            lines = lines.strip()
            exp = re.split("\t",lines)
            self.comb_parents = OrderedDict()
            gene_hp_list = []

            if exp:
                gene_id = exp[0].strip()
                gene_id_list.append(gene_id)
                gene_hp_strs = re.split("\s|\,|\;",exp[1])
                gene_hp_strs = [x for x in gene_hp_strs if x]

                for ele_hp in gene_hp_strs:
                    ele_hp = ele_hp.strip()
                    if ele_hp:
                        exp_strs = re.split("\:",ele_hp)
                        hp_num = exp_strs[1]
                        gene_hp_list.append(hp_num)

                        if self.hpo_lib.has_key(hp_num):
                            parents_list = self.hpo_lib[hp_num]
                            for ele in parents_list:
                                self.comb_parents[ele] = 1
                            parents_list = []
                        else:
                            if not HPOMethod.no_hp.has_key(hp_num):
                                try:
                                    parents_list = self.getHPOParents(HPO_hash,hp_num)
                                    if parents_list:
                                        self.hpo_lib[hp_num] = parents_list
                                        for ele in parents_list:
                                            self.comb_parents[ele] = 1
                                    else:
                                        no_hp.append(hp_num)

                                except IndexError:
                                    print "Error occurred:",hp_num

                                HPOMethod.parents = []
                                HPOMethod.sorted_parents = OrderedDict()
                                parents_list = []

                if self.comb_parents:
                    print >>op, gene_id+"\t"+",".join(gene_hp_list)+"\t"+','.join(self.comb_parents.keys())
                else:
                    print "Gene Ids with HPO Ids having NO Parent HPO terms: " + lines
                    print >>op, gene_id+"\t"+",".join(gene_hp_list)

            lines = fh.readline()

        if HPOMethod.no_hp:
            for keys in HPOMethod.no_hp:
                print >>ng,keys
        self.comb_parents = OrderedDict()
        HPOMethod.no_hp = OrderedDict()

        op.close()
        fh.close()
        ng.close()


    def processRawHPOFile(self, inp_file, out_file):

        hgnc_hpo = OrderedDict()
        try:
            fh = open(inp_file)
            wh = open(out_file,"w")

            lines = fh.readline()
            if re.search("#",lines):
                lines = fh.readline()

            while lines:
                lines = lines.strip()
                strs = re.split('\t',lines)
                hgnc_id = strs[1].strip().upper()
                hpo_id = strs[3].strip()

                #if hgnc_ensg.has_key(hgnc_id):
                if hgnc_hpo.has_key(hgnc_id):
                    tmp = hgnc_hpo[hgnc_id]
                    tmp.append(hpo_id)
                    hgnc_hpo[hgnc_id] = tmp
                else:
                    hgnc_hpo[hgnc_id] = [hpo_id]

                lines = fh.readline()

            for keys in hgnc_hpo:
                #if hgnc_ensg.has_key(keys):
                    #ensg_id = hgnc_ensg[keys][0]
                hgnc_id = keys
                print >>wh,hgnc_id+"\t"+",".join(hgnc_hpo[keys])
            wh.close()

        except IOError:

            print "Not a valid path/file. Please enter correct path or file"
            sys.exit()
        except:
            sys.exit()


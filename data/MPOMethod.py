#!/usr/bin/python

import re,sys,os
from collections import OrderedDict

class MPOMethod:

    parents =[]
    sorted_parents = OrderedDict()
    children = []#Added for incorporating getChildren function
    sorted_children = OrderedDict()
    alt_id = OrderedDict()
    no_mp = OrderedDict()
    altid_flag_val = 't'

    def getMPOParents(self, MP_hash, mp_id):

        self.MP_hash = MP_hash
        self.mp_id = mp_id

        if self.mp_id == '0000001':
            if not mp_id in MPOMethod.parents:
                MPOMethod.parents.append(self.mp_id)

        else:
            if self.MP_hash.has_key(self.mp_id):
                keys = self.mp_id
                if self.MP_hash[keys]['is_a']:
                    isA_mp = self.MP_hash[keys]['is_a']
                    for e in isA_mp:
                        MPOMethod.parents.append(e)
                        self.getMPOParents(self.MP_hash,e)
                else:
                    if MPOMethod.altid_flag_val =='t':
                        print "Looking for Alternate ids....."
                        flag=0
                        for keys in self.MP_hash:
                            alt_id_tmp = self.MP_hash[keys]['alt_id']
                            for ele in alt_id_tmp:
                                if ele==self.mp_id:
                                    flag=1
                                    MPOMethod.alt_id[self.mp_id] = keys
                                    if self.MP_hash[keys]['is_a']:
                                        isA_mp = self.MPO_hash[keys]['is_a']
                                        for e in isA_mp:
                                            self.getMPOParents(self.MP_hash,e)
                        if flag==0:
                            MPOMethod.no_mp[self.mp_id]=1
                            print "No alternate ID found... for DO:"+self.mp_id

        for outer in MPOMethod.parents:
            MPOMethod.sorted_parents[outer]=1

        return MPOMethod.parents


    def getFile_MPO(self, MP_hash, inp_file, out_gene_mp, out_no_mp, altid_val):

        self.mp_lib = OrderedDict()
        MPOMethod.altid_flag_val = altid_val

        op = open(out_gene_mp,'w')
        ng = open(out_no_mp,'w')

        gene_id_list = []
        #self.gene_hpo_hash = {}
        parents_list = []
        no_mp = []

        fh = open(inp_file)

        lines = fh.readline()

        while lines:
            lines = lines.strip()
            exp = re.split("\t",lines)
            self.comb_parents = OrderedDict()
            gene_mp_list = []

            if exp:
                gene_id = exp[0].strip()
                gene_id_list.append(gene_id)
                gene_mp_strs = re.split("\s|\,|\;",exp[1])
                gene_mp_strs = [x for x in gene_mp_strs if x]

                for ele_mp in gene_mp_strs:
                    ele_mp = ele_mp.strip()
                    if ele_mp:
                        exp_strs = re.split("\:",ele_mp)
                        mp_num = exp_strs[1]
                        gene_mp_list.append(mp_num)

                        if self.mp_lib.has_key(mp_num):
                            parents_list = self.mp_lib[mp_num]
                            for ele in parents_list:
                                self.comb_parents[ele] = 1
                            parents_list = []
                        else:
                            if not MPOMethod.no_mp.has_key(mp_num):
                                try:
                                    parents_list = self.getMPOParents(MP_hash,mp_num)
                                    if parents_list:
                                        self.mp_lib[mp_num] = parents_list
                                        for ele in parents_list:
                                            self.comb_parents[ele] = 1
                                    else:
                                        no_mp.append(mp_num)

                                except IndexError:
                                    print "Error occurred:",mp_num

                                MPOMethod.parents = []
                                MPOMethod.sorted_parents = OrderedDict()
                                parents_list = []

                if self.comb_parents:
                    print >>op, gene_id+"\t"+",".join(gene_mp_list)+"\t"+','.join(self.comb_parents.keys())
                else:
                    print "Gene Ids with HPO Ids having NO Parent HPO terms: " + lines
                    print >>op, gene_id+"\t"+",".join(gene_mp_list)

            lines = fh.readline()

        if MPOMethod.no_mp:
            for keys in MPOMethod.no_mp:
                print >>ng,keys
        self.comb_parents = OrderedDict()
        MPOMethod.no_mp = OrderedDict()

        op.close()
        fh.close()
        ng.close()


    def processRawMPOFile(self,inp_hmd,inp_mgi,out_file):

        hmd_hash = OrderedDict()
        mgi_hash = OrderedDict()

        try:
            hmd = open(inp_hmd)
            mgi = open(inp_mgi)

            tmp_file = os.path.dirname(out_file)+'/tmp_preprocess.txt'
            wh1 = open(tmp_file,"w")

            for lines in hmd:
                lines = lines.strip()
                strs = re.split("\t",lines)
                try:
                    mgi_id = strs[5].strip()
                except:
                    mgi_id = []

                try: 
                    hgnc_id = strs[0].strip().upper()
                except:
                    hgnc_id = []

                try:
                    mpo_id = strs[6].strip()
                    mpo_strs = re.split("\s",mpo_id)
                except:
                    mpo_strs = []
                    mpo_id = []

                if hmd_hash.has_key(mgi_id):

                    tmp = hmd_hash[mgi_id]['mpo']
                    for ele in mpo_strs:
                        tmp.append(ele)

                    hmd_hash[mgi_id] = {'hgnc':hgnc_id,'mpo':list(set(tmp))}
                    tmp = []

                else:
                    hmd_hash[mgi_id] = {'hgnc':hgnc_id,'mpo':mpo_strs}

            

            for lines in mgi:

                lines = lines.strip()
                strs = re.split("\t",lines)
                
                try:
                    mgi_id = strs[5].strip()
                    mgi_strs = re.split("\,",mgi_id)
                except:
                    mgi_strs = []
                    mgi_id = []

                try:
                    mpo_id = strs[3].strip()
                except:
                    mpo_id = []

                for mgi_id in mgi_strs:
                    if mgi_hash.has_key(mgi_id):
                        tmp = mgi_hash[mgi_id]
                        tmp.append(mpo_id)
                        mgi_hash[mgi_id] = tmp
                        tmp = []
                    else:
                        mgi_hash[mgi_id] = [mpo_id]

            hmd_mgi_diff = set(hmd_hash.keys()).difference(set(mgi_hash.keys()))
            mgi_hmd_diff = set(mgi_hash.keys()).difference(set(hmd_hash.keys()))

            for key in hmd_hash:
                hgnc_id = hmd_hash[key]["hgnc"]
                mpo_hmd = hmd_hash[key]["mpo"]
                mpo_comb = []
                if mgi_hash.has_key(key):
                    mpo_mgi = mgi_hash[key]
                    mpo_comb = mpo_mgi+mpo_hmd
                    print >>wh1,hgnc_id.upper()+"\t"+",".join(list(set(mpo_comb)))
                else:
                    if mpo_hmd:
                        #mpo_comb = list(set(mpo_hmd))
                        print >>wh1,hgnc_id.upper()+"\t"+",".join(list(set(mpo_hmd)))

            wh1.close()
            hmd.close()
            mgi.close()

            #Process the tmp file and combine all MPO terms for input HGNC genes
            fh = open(tmp_file) 
            wh = open(out_file,'w')

            hgnc_hash = {}

            for lines in fh:
                lines = lines.strip()
                strs = re.split('\t',lines)
                strs = [x.strip() for x in strs]
                hgnc_id = strs[0].upper()
                mpo_strs = re.split('\,',strs[1])

                if hgnc_hash.has_key(hgnc_id):
                    tmp = hgnc_hash[hgnc_id]
                    tmp = tmp+mpo_strs
                    hgnc_hash[hgnc_id] = list(set(tmp))
                else:
                    hgnc_hash[hgnc_id] = mpo_strs

            fh.close()
            
            for keys in hgnc_hash:
                print >>wh,keys+'\t'+','.join(hgnc_hash[keys])
            wh.close()

        except IOError:
            print "Not a valid path/file. Please enter correct path or file"
            sys.exit()
        #except:
        #    sys.exit()


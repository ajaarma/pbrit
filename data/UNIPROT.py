#!/usr/bin/python


##################################################################################
# UNIPROT: General Class of UNIPROT. Consist of mehtods for parsing the UNIPROT ID Mapping file for HUMAN.
#
# Consist of methods for extracting respective columns of the file which is accroding to the annotaiton sources (GO, PUBMED Ids, UNIPROT-KB)
#
#
#
###################################################################################

import re,sys



class UNIPROT:

    def __init__(self,elements=[]):
        self.__elements={}
        for e in elements:
            self.__elements[e]=1
        
    def display(self):
        print "Inside the UNIPROT Class"

    

    def parseUNIP(self, inp_file, out_file, col_ind, ensg_hgnc):

        try:
            self.fh = open(inp_file)
            self.wh = open(out_file,"w")
        except:
            print "The provided file could not be found. Please enter the correct path or file name"
            sys.exit()  

        self.go_hash = {}
        
        for lines in self.fh:
            lines = lines.strip()
            strs = re.split("\t",lines) 
    
            self.go_list = []
            self.ensg_list = []
            
            for e in strs[0:8]:
                if re.search("^GO\:",e):
                    self.go_list = e

            for e in strs[len(strs)-7:len(strs)]:
                if re.search("^ENSG",e):
                    self.ensg_list = e

            if len(self.go_list) !=0 and len(self.ensg_list) !=0:
                
                strs_ensg = re.split("\;",self.ensg_list)
                for ele in strs_ensg:
                    if ensg_hgnc.has_key(ele):
                        print >>self.wh, ele+"\t"+self.go_list
                        break
        
        self.fh.close()
        self.wh.close()


    def processENSG_GO (self, inp_file, out_file):

        ensg_go = {}
        self.fh = open(inp_file)
        print "Processing Ensg GO file"
        
        for lines in self.fh:
            lines = lines.strip()
            strs = re.split("\t",lines)
            ensg_id = strs[0].strip()
            go_strs = re.split("\;",strs[1].strip())
            
            if ensg_go.has_key(ensg_id):
                tmp = ensg_go[ensg_id]
                tmp = tmp+go_strs
                ensg_go[ensg_id] = list(set(tmp))
            else:
                ensg_go[ensg_id] = go_strs

        self.fh.close()
        self.wh = open(out_file,"w")

        for keys in ensg_go:
            print >>self.wh,keys,"\t",";".join(ensg_go[keys])
        
        self.wh.close()

    def processENSG_PUB(self, inp_file,out_pub, ensg_hgnc):

        self.fh = open(inp_file)
        self.wh = open(out_pub,"w")
        ensg_pub = {}
        pub_hash = {}

        for lines in self.fh:
            lines = lines.strip()
            strs = re.split("\t",lines)
            #print strs
            self.pub_list = strs[15].strip()
            #try:
            #   self.add_pub_list = strs[21].strip()
            #except:
            #   self.add_pub_list = []
            #   pass
            self.add_pub_list = []
            self.ensg_list = []
            
            for e in strs[len(strs)-7:len(strs)]:
                if re.search("^ENSG",e):
                    self.ensg_list = e  
            
            if (len(self.pub_list) !=0 and len(self.ensg_list) !=0) :

                pub_strs = re.split("\;",self.pub_list)
                if self.add_pub_list:
                    add_pub_strs = re.split("\;",self.add_pub_list)
                    pub_strs = pub_strs+add_pub_strs
                else:
                    pub_strs = pub_strs
                
                for e in pub_strs:
                    e = e.strip()
                    pub_hash[e] = 1

                strs_ensg = re.split("\;",self.ensg_list)
                for ele in strs_ensg:
                    if ensg_hgnc.has_key(ele):
                        print >>self.wh, ele+"\t"+";".join(list(set(pub_strs)))
                        break
        self.fh.close()
        self.wh.close()

        self.wh = open(re.split("ENSG_|tmp_|.txt",out_pub)[0]+"PUB_UNIQUE.txt","w")

        for keys in pub_hash:
            print >>self.wh,keys

        self.wh.close()

    def mapPPIUniprotEnsembl(self,inp_file,out_file,ppi_file):
        ''' Function to process Input UNIPROT File and extract Uniprot Accession ID 
            and Ensembl ID '''

        fh = open(inp_file)
        wh = open(out_file,'w')

        #Processing Uniprot-Accession PPI file of pBRIT internal database
        ph = open(ppi_file)
        ppi_hash = {}

        for lines in ph:
            lines = lines.strip()
            strs = re.split('\t',lines)
            ppi_id = strs[1].strip()
            ppi_hash[ppi_id] = 1

        ph.close()

        #Processing Uniprot database to extract Uniprot Accession and Ensembl ID
        unip_ensg = {}

        for lines in fh:
            lines = lines.strip()
            strs = re.split('\t',lines) 
            strs = [x.strip() for x in strs]
            unip_id = re.split('\_',strs[1])[0].upper()

            ensg_list = []
            for e in strs[len(strs)-7:len(strs)]:
                if re.search('^ENSG',e):
                    ensg_list = e

            #Create unip-ensg Dictionary
            if unip_ensg.has_key(unip_id):
                tmp = unip_ensg[unip_id]
                tmp.append(ensg_id)
                unip_ensg[unip_id] = list(set(tmp))
            else:
                if ensg_list:
                    ensg_id = re.split('\;|\,',ensg_list)[0]
                    unip_ensg[unip_id] = [ensg_id]



        fh.close()

        #Output PPI-ID and corresponding Ensembl map
        for ppi_id in ppi_hash:
            ppi_id = ppi_id.upper()
            if unip_ensg.has_key(ppi_id):
                print >>wh,ppi_id+'\t'+unip_ensg[ppi_id][0]

        wh.close()

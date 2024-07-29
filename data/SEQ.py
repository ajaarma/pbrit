#!/usr/bin/python

import re,sys,os

class SEQMethod:

    def __init__(self,elements=[]):
        self.__elements={}
        for e in elements:
            self.__elements[e]=1

    def getDBUniqueGeneNames(self,configDict,objH):
        ''' Retrive all the unique gene names from all 
            the internal databases of pBRIT '''

        version = configDict['version']
        data_dir = configDict['general']['dataDir']+'/'+version

        anno_src_list = configDict['general']['annoSource']['value']
        ensg_hgnc_file = data_dir+'/HUGO/'+configDict['HUGO']['ensgHugo']
        hgnc_ensg_file = data_dir+'/HUGO/'+configDict['HUGO']['hugoEnsg']

        hgnc_dict = {}
        ensg_dict = {}

        for anno in anno_src_list:
            try:
                anno_dir = os.path.abspath(data_dir+'/'+anno+'/INPUT')
                anno_row_file = anno_dir+'/HGNC_rownames.txt'
                fh = open(anno_row_file)
                for lines in fh:
                    lines = lines.strip()
                    strs = re.split('\t',lines);
                    strs = [x.strip() for x in strs]
                    if not anno in ['SEQ','HUGO','UNIPROT','PUB']:
                        hgnc_id = strs[1]
                        hgnc_dict[hgnc_id] = 1
                    elif anno in ['PUB']:
                        ensg_id = strs[0]
                        ensg_dict[ensg_id] = 1
                fh.close()
            except:
                pass

        ensg_hugo_dict = objH.getEnsgHgncDict(ensg_hgnc_file)
        hugo_ensg_dict = objH.getEnsgHgncDict(hgnc_ensg_file)

        ensg_uniq_genes_db = []

        for hgnc_id in hgnc_dict:
            if hugo_ensg_dict.has_key(hgnc_id):
                ensg_id = hugo_ensg_dict[hgnc_id][0]
                ensg_uniq_genes_db.append(ensg_id)

        for engs_id in ensg_dict.keys():
            ensg_uniq_genes_db.append(ensg_id)

        return list(set(ensg_uniq_genes_db)) 
            

    def getQueryTargetBitScore(self,seqOutFile,seqOutBit,ensg_genes_db):
        ''' Subroutine to process BLAST Output sequences '''

        fh = open(seqOutFile)
        wh = open(seqOutBit,'w')

        print '-- Parsing the BLAST Output results'
        count = 0
        for lines in fh:
            lines = lines.strip()
            if re.search('^#',lines):
                continue
            else:
                strs = re.split('\t',lines);
                strs=[x.strip() for x in strs if x]
                query_id = re.split('\|',strs[0])
                target_id = re.split('\|',strs[1])
                bit_score = str(strs[11])
                print >>wh,'\t'.join([query_id[0],target_id[0],bit_score])
                count = count+1

        fh.close()
        wh.close()

    def getQueryTargetNormBitScore(self,seqOutBit,seqOutNorm,ensg_genes_db):

        fh = open(seqOutBit)

        query_hash = {}
        
        print ' -- Reading the parsed BLAST output and BIT score file'
        for lines in fh:
            lines = lines.strip()
            strs = re.split('\t',lines);
            strs = [x.strip() for x in strs]
            query_id = strs[0]
            target_id = strs[1]
            strs_jt = query_id+'_'+target_id

            if query_hash.has_key(strs_jt):
                exp = re.split('\,',query_hash[strs_jt])
                exp.append(strs[2])
                query_hash[strs_jt] = ','.join(exp)
            else:
                query_hash[strs_jt] = strs[2]

        fh.close()

        wh = open(seqOutNorm,'w')

        print ' -- Normalizing the BIT scores'
        for keys in query_hash:
            exp1 = re.split('\_',keys)
            exp2 = re.split('\,',query_hash[keys])
            exp2 = [float(x.strip()) for x in exp2]
            #print >>wh,exp1,'\t',exp2
            #if exp1[0] in ensg_genes_db:
            print >>wh,'\t'.join([keys,exp1[0],exp1[1],
                                        str(sum(exp2)/len(exp2))
                                 ]
                                )

        wh.close()

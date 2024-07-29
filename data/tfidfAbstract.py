#!/usr/bin/python

################################################################
# Program to generate sparse binary matrix format for Abstract
#
# Usage: python tfidfExtract.py  pubmed_abstrsct_copy.txt
#
# Output: Generates two files
#         a) Stemmed pubmed abstract
#         b) Word count hash table
#
# Author: Ajay Anand Kumar. Faculty of medical Genetics
# mail:        ajay.kumar@ua.ac.be
#  
################################################################


import sys
import re
from gensim import *
from nltk.stem import SnowballStemmer	

class AbstractProcess:

    master_word_hash= {}

    def __init__(self, elements=[]):
        #pass
        self.__elements={}
        for e in elements:
            self.__elements[e] = 1

    def processAbstract(self,abstr_str_arg):
        words_hash = {}
        words_list = []
        abstr = abstr_str_arg
        abstr_str = re.split("[\s|\,|\W|\_]",abstr)
        abstr_str = [x for x in abstr_str if x]
        for ele in abstr_str:
            ele = ele.strip()
            if not ele.isdigit():
                words_list.append(ele)
        abstr_processed = ' '.join(words_list)
   
        return abstr_processed

    def getStopWordsFile(self,stopWordFile,abstr_str_arg):

        stop_word = {}

        fh = open(abstr_str_arg)
        for lines in fh:
            lines = lines.strip()
            strs = re.split("\t",lines)
            abst = strs[1].strip()
            abst_processed = self.processAbstract(abst)
            abst_strs = re.split("\s",abst_processed)
            for e in abst_strs:
                e = e.lower()
                if stop_word.has_key(e):
                    count = stop_word[e]
                    count = count+1
                    stop_word[e] = count
                else:
                    stop_word[e] = 1


        wh = open(stopWordFile,"w")
        for keys in sorted(stop_word,key=stop_word.get,reverse=True):
            print >>wh,keys,"\t",stop_word[keys]
        
        wh.close() 
            
 
    def processStopWords(self,abstr_str_arg,stopWordFile):
        #fh = open("stop_words.txt")
        fh = open(stopWordFile)
        stop_hash = {}
        for lines in fh:
            lines = lines.strip()
            stop_hash[lines] = 1
       
        fh.close()
        stoplist1 = 'the of in and they to that is by with for we was are as from this an were these which on not or have at be but also has been it into more their most they there a b c d e f g h i j k l m n o p q r s t u v w x y z'
        stoplist1_strs = re.split("\s",stoplist1)
        for ele in stoplist1_strs:
            stop_hash[ele] = 1
        stoplist = " ".join(stop_hash.keys())
        abstr_processed = re.split("\s",abstr_str_arg)
        documents = abstr_processed
        texts = [[word for word in document.lower().split() if word not in stoplist]  for document in documents]
        text_buffer = []
        for word in texts:
            if word:
                text_buffer.append(word[0])
                AbstractProcess.master_word_hash[word[0]] = 1
     
        return(' '.join(text_buffer))
  
    def processStemming(self,abstr_str_arg):
 
        abstr_processed = re.split("\s",abstr_str_arg)
        stemmer = SnowballStemmer("english")
 
        stem_words = []
        for word in abstr_processed:
            word = word.strip()
            word = stemmer.stem(word)     
            stem_words.append(str(word))

        final_text = ' '.join(stem_words)#stem_words.keys()
      
        return final_text

    def getPubSparseMat(self,stemmed_file,row_file,col_file,sparse_file):

        fh = open(stemmed_file)
        cfh = open(stemmed_file)
      
        #wrh = open("pubmed_abstract_rownames.txt","w")
        wrh = open(row_file,"w")
        wch = open(col_file,"w")
        wsh = open(sparse_file,"w")

        words_index_hash = {}
        gene_word_freq_hash = {}
      
        index_count = 1 

        for lines in fh:
            word_freq_hash = {}
            lines = lines.strip()
            try:
                strs = re.split("\t",lines)
                gene_id = strs[0].strip()
                gene_id_count = 1
         
                abstr_str = re.split("\s",strs[1].strip())
            except IndexError:
                abstr_str = []
            count = []
            for ele in abstr_str:
                ele = ele.strip()
                if word_freq_hash.has_key(ele):
                    count = word_freq_hash[ele]
                    count  = count+1
                    word_freq_hash[ele] = count
                else:
                    word_freq_hash[ele] = 1
                    words_index_hash[ele] = 1
                    count = [] 
         
            gene_word_freq_hash[gene_id] = word_freq_hash         	


        master_word_list = words_index_hash.keys()
        for ele in master_word_list:
            ele = ele.strip()
            words_index_hash[ele] = index_count
            index_count = index_count+1
      
        for keys in sorted(words_index_hash,key=words_index_hash.get,reverse=False):
            print >>wch,keys,"\t",words_index_hash[keys]
         
         
        gene_row_ind = 1
        for lines in cfh:
            lines = lines.strip()
            strs = re.split("\t",lines)
            try:
                gene_id = strs[0].strip()
                print >>wrh,gene_id         
                abstr_str = re.split("\s",strs[1].strip()) 
                for ele in abstr_str:
                    ele = ele.strip()
                    col_ind = words_index_hash[ele]
                    word_freq_hash = gene_word_freq_hash[gene_id]
                    freq_count = word_freq_hash[ele]
                    print >>wsh,gene_row_ind,"\t",col_ind,"\t",freq_count	
                    gene_row_ind = gene_row_ind+1     
            except IndexError:
                gene_row_ind = gene_row_ind+1           
        
        wrh.close()
        wch.close()
        wsh.close()

    def processRawPubFile(self,rawFile,stem_file,stopWord,matRow,matCol,matSpr,objMat):

        fh = open(rawFile)
        wh = open(stem_file,"w")

        obj = AbstractProcess() 
        #obj.getStopWordsFile(stopWord,rawFile)
        #sys.exit()
        print "\nStemming the Abstract file\n"
        for lines in fh:
            lines = lines.strip()
            strs = re.split("\t",lines)
            gene_id = strs[0].strip()
            abst = strs[1].strip()
            abstr_processed = obj.processAbstract(abst)
            text_processed  = obj.processStopWords(abstr_processed,stopWord)
            text_stemm = obj.processStemming(text_processed)
            print >>wh,gene_id+"\t"+text_stemm
        
        
        fh.close()
        wh.close()
        
        print "\n  Generating Sparse Matrices"
        #objMat = MATRICES()
        objMat.getPubSparseMat(stem_file,matRow,matCol,matSpr)


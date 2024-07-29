#!/usr/bin/python

import re,os
from collections import OrderedDict
import itertools


class PPMethod:

    def __init__(self, elements=[]):
        self.__elements={}
        for e in elements:
            self.__elements[e]=1


    def processRawPPFile(self,annoFile,genePP,matRow,matCol,matSpr):

        pp_hash = OrderedDict()
        col_hash = {}

        fh = open(annoFile)
        wh = open(genePP,"w")
        wr = open(matRow,"w")
        wc = open(matCol,"w")
        ws = open(matSpr,"w")

        for lines in fh:
            lines = lines.strip()
            if not re.search("^#",lines):
                strs = re.split("\t",lines)
                int_ppt = strs[2].strip()
                #print strs
                try:
                    strs1 = re.split("\,|\_HUMAN|\,",int_ppt)
                    #print strs1
                    a = [x.upper() for x in strs1 if x]
                    tmp = list(itertools.permutations(a,2))
                    for e in tmp:
                        print >>wh,"\t".join(list(e)),"\t",strs[3].strip()
                    for ele in a:
                        ele = ele.strip()
                        if (ele !='interaction_participants'):
                            if re.search("\.",ele):
                                strs2 = re.split("\.",ele)
                                pp_hash[strs2[1].strip()] = 1
                            else:
                                pp_hash[ele] = 1
                except IndexError:
                    pass

        fh.close()
        wh.close()

        count = 1
        for keys in pp_hash:
            col_hash[keys] = count
            count = count+1

        for keys in sorted(col_hash,key=col_hash.get, reverse=False):
            print >>wr,col_hash[keys],'\t',keys
            print >>wc,keys,"\t",col_hash[keys]
        wr.close()
        wc.close()

        fh = open(genePP)

        for lines in fh:
            lines = lines.strip()
            strs = re.split("\t",lines)
            a1 = strs[0].strip()
            a2 = strs[1].strip()
            if re.search("\.",a1):
                strs1 = re.split("\.",a1)
                a1 = strs1[1].strip()
            else:
                a1 = a1.strip()

            if re.search("\.",a2):
                strs1 = re.split("\.",a2)
                a2 = strs1[1].strip()
            else:
                a2 = a2.strip()
            
            if re.search("NA",strs[2]):
                #print lines
                print >>ws,col_hash[a1],"\t",col_hash[a2],"\t",0,"\t",strs[2]
            else:
                print >>ws,col_hash[a1],"\t",col_hash[a2],"\t",1,"\t",strs[2]


        fh.close()
        ws.close()

        


#!/usr/bin/python

import re,sys,os

if __name__=="__main__":

    inp_file = os.path.abspath(sys.argv[1])
    ont_row = os.path.abspath(sys.argv[2])

    fh = open(inp_file)

    hgnc_ensg = {}
    ensg_hgnc = {}

    for lines in fh:
        lines = lines.strip()

        strs = re.split('\t',lines);
        approv_symb = re.split('\,',strs[1].strip())
        
        #try:
        #    prev_symb = re.split('\,',strs[4].strip())
        #except:
        #    prev_symb = []
        try:
            alias_symb = re.split('\,',strs[5].strip())
        except:
            alias_symb = []

        #symbol_strs = approv_symb+prev_symb+alias_symb
        symbol_strs = approv_symb+alias_symb
        symbol_strs = [x.strip().upper() for x in symbol_strs if x]
        
        try:
            ensg_id = re.split('\,',strs[9].strip())
            
        except:
            ensg_id = ['NA']

        #print 'outside: ',symbol_strs,'\t',ensg_id
        if ensg_id[0] !='' and ensg_id[0] !='NA':
            #print 'Inside: ',symbol_strs,'\t',ensg_id
            for symbol_id in symbol_strs:
                symbol_id = symbol_id.upper()

                if hgnc_ensg.has_key(symbol_id):
                    tmp = hgnc_ensg[symbol_id]
                    tmp.append(ensg_id[0])
                    hgnc_ensg[symbol_id] = list(set(tmp))
                else:
                    hgnc_ensg[symbol_id] = list(set(ensg_id))
            
            for ele in ensg_id:
                if ensg_hgnc.has_key(ele):
                    tmp = ensg_hgnc[ele]
                    tmp += symbol_strs
                    tmp = [x.upper() for x in tmp]
                    ensg_hgnc[ele] = list(set(tmp))
                else:
                    ensg_hgnc[ele] = symbol_strs

    fh = open(ont_row)

    for lines in fh:
        lines = lines.strip()
        strs = re.split('\t',lines)
        hgnc_id = strs[1].strip()
        if hgnc_ensg.has_key(hgnc_id):
            print hgnc_id,'\t',hgnc_ensg[hgnc_id][0]


    #for hgnc_id,ensg_id in hgnc_ensg.items():
    #for hgnc_id,ensg_id in ensg_hgnc.items():
        #if len(ensg_id)>1:
        #    print hgnc_id,'\t',ensg_id
        #print hgnc_id,'\t',ensg_id



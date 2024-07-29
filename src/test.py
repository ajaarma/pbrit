#!/usr/bin/python

import re,sys,os

script_path = os.path.dirname(os.path.abspath(__file__))
sys.path.append(script_path+'/UNIPROT/')

from UNIPROT import *

if __name__=='__main__':

    inp_file = sys.argv[1]

    out_file = sys.argv[2]

    ppi_file = sys.argv[3]

    objU = UNIPROT()

    objU.mapPPIUniprotEnsembl(inp_file,out_file,ppi_file)


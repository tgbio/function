#!/usr/bin/env python

import os
import re
import subprocess

from Bio.Blast import NCBIXML
from Bio import Entrez
from Bio import SeqIO

#bogus stuff here 

genotypes = ["Cvi", "Ler"]
roi = {1:[[174952, 3112850], [23698501, 24905939]]} #hardwired, change to accept a variable from another script

Entrez.email = ""

print "starting with:"
os.getcwd()
for g in genotypes: # g = geno type
    geno = g 
    print geno
    geno_path = "./" + geno
    os.chdir(geno_path)
    os.getcwd()
    for ch in roi: # roi is a dictionary, with key = chrom, value = list of positions 
        print ch        
        for reg in roi[ch]:
            print "Processing results for ", geno, ch, reg
            reg_descrip = "Chr_" + str(ch) + "_Reg_" + str(reg[0]) + "_" + str(reg[1])
            reg_dir_name = "./" + reg_descrip
            
            os.chdir(reg_dir_name)
            print os.getcwd()
            file_list = [f for f in os.listdir(".") if os.path.isfile(os.path.join(".", f))]
            print "file_list", file_list
            tree_dict = {}
            msa_dict = {}
            for fl in file_list:
                if fl[-12:] == "_semphy.tree":
                    tree_dict[fl[:-12]] = fl 
                elif fl[-15:] == "_probcons.fasta":
                    msa_dict[fl[:-15]] = fl
            assert len(tree_dict) == len(msa_dict)
            for e in tree_dict:
                tree_file = tree_dict[e]
                msa_file = msa_dict[e]
                mapp_command = "-jar ~/MAPP.jar -f " + "./" + msa_file + " -t " + "./" + tree_file + " -o " + e + "_MAPP_results.txt"                    
                print mapp_command                
                subprocess.call(mapp_command, shell=True)
            os.chdir("..")
    os.chdir("..")

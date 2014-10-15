#parsing the xml output
import os
import re
import subprocess

from Bio.Blast import NCBIXML
from Bio import Entrez
from Bio import SeqIO

genotypes = ["Cvi", "Ler"]
roi = {1:[[174952, 3112850], [23698501, 24905939]]}

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
            for fl in file_list:
                if fl[-15:] == "_probcons.fasta":
                    print "mfa file:", fl
                    probcons_msa = fl
                    print "Sending probcons multiple sequence alignment to semphy"
                    semphy_out = fl[:-15] + "_semphy.txt"
                    semphy_tree = fl[:-15] + "_semphy.tree"
                    semphy_log = fl[:-15] + "_semphy_log.txt"
                    semphy_command = "semphy -s "+ probcons_msa +" -o " + semphy_out + " -T " + semphy_tree + " -l " + semphy_log + " -a 20 --jtt --posteriorDTME -O -v 5"
                    print semphy_command
                    subprocess.call(semphy_command, shell=True)
            os.chdir("..")
    os.chdir("..")# <codecell>
#sending the psi blast results to probcons 
#change directory to probcons

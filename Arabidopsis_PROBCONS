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
            print "Parsing psi-blast results for ", geno, ch, reg
            reg_descrip = "Chr_" + str(ch) + "_Reg_" + str(reg[0]) + "_" + str(reg[1])
            reg_dir_name = "./" + reg_descrip
            
            os.chdir(reg_dir_name)
            os.getcwd()
            file_list = [f for f in os.listdir(".") if os.path.isfile(os.path.join(".", f))]
            print "file_list", file_list     
            for fl in file_list:
                if fl[-4:] == ".xml":
                    print "xml file", fl
                    seq_out_name = fl[:-4] + ".fasta"
                    org_fasta_file = fl[:-8] + ".fasta"
                    result_handle = open(fl)
                    blast_records = NCBIXML.parse(result_handle)
                    evalue = 0.05
                    id_list = []
                    evalue_dict = {}
                    for blast_record in blast_records:
                        for alignment in blast_record.alignments:
                            for hsp in alignment.hsps:
                                if hsp.expect < evalue:
                                    print "e-value", hsp.expect
                                    #print "Alignment title", alignment.title
                                    start_gi = re.search("gi\|", alignment.title) #within alignment.title retrieve the gene ids of the sequences gene id follows "gi|" and precedes "|ref|"
                                    end_gi = re.search("\|ref\|", alignment.title)
                                    #print "start_gi.end()", start_gi.end() 
                                    #print "end_gi.start()", end_gi.start() 
                                    
                                    gi_id = str(alignment.title[start_gi.end():end_gi.start()])
                                    if evalue in evalue_dict:
                                        evalue_dict[hsp.expect].append(gi_id)
                                    else:
                                        evalue_dict[hsp.expect] = [gi_id]
                    sort_evals = [k for k in sorted(evalue_dict.keys())] #sort the keys (evalues) in ascending order
                    id_list = []
                    for se in sort_evals:
                        for gi in evalue_dict[se]: #in case there is more than one gene id for a single evalue
                            if gi in id_list:#keep only unique ids in the id list
                                continue
                            else:
                                id_list.append(gi)
                    if len(id_list) >99:# capture only the 99 ids with the lowest e-values
                        id_list = id_list[:99]
                    
                    print len(id_list)                    
                    
                    #return the full sequences corresponding to the ids of the lowest 9 values in the psi-blast hits (id_list) plus the query sequence
                    gi_str = ",".join(id_list) # based on Biopython tutorial 9.14.3  Searching, downloading, and parsing GenBank records
                    handle = Entrez.efetch(db="protein", id=gi_str, rettype="gb", retmode="text")  #look up those ids at ncbi, retrieve their sequences and write the fasta seqs to file
                    records = SeqIO.parse(handle, "gb")
                    count_records = 0
                    
                    with open(seq_out_name, "w") as psi_seqs: 
                        for seq_record in SeqIO.parse(org_fasta_file, "fasta"):
                            org_header = ">" + str(seq_record.id) + "\n"
                            psi_seqs.write(org_header)
                            org_seq = str(seq_record.seq) + "\n"
                            psi_seqs.write(org_seq)
                            count_records += 1
                        for record in records:
                            hit_header = ">"+ str(record.id) + "\n"
                            psi_seqs.write(hit_header)
                            hit_seq = str(record.seq) + "\n"
                            psi_seqs.write(hit_seq)
                            count_records += 1
                    handle.close()
                    print "Wrote", count_records, "sequences from psi-blast results"
                    print "Sending", count_records, "sequences to probcons"
                    probcons_msa = fl[:-8] + "_probcons.mfa" 
                    probcons_command = "probcons " + seq_out_name +  " > " + probcons_msa
                    subprocess.call(probcons_command, shell=True) # submit the evalue sorted psi-blast results to multiple sequence alignment with probcons
                    
                    #print "Sending probcons multiple sequence alignment to semphy"
                    #semphy_out = fl[:-8] + "_semphy.txt"
                    #semphy_tree = fl[:-8] + "_semphy.tree"
                    #semphy_command = "semphy -s probcons_msa -o semphy_results_test.txt -T sempthy_test.tree -l semphy_log.txt -a 20 --jtt --posteriorDTME -O -v 5"
                    #subprocess.call(semphy_command, shell=True)
            os.chdir("..")
    os.chdir("..")# <codecell>
#sending the psi blast results to probcons 
#change directory to probcons
#

############################
from Bio import SeqIO
#from Bio.Seq import Seq
#from BCBio.GFF import GFFExaminer
#from BCBio import GFF
#from Bio.Alphabet import IUPAC
#import pprint
import csv
import re
#import sys
#import string
#import os

#from Bio.SeqRecord import SeqRecord
#from Bio.Alphabet import generic_dna
from Bio import SeqFeature
#from Bio.SeqFeature import FeatureLocation
#from Bio.SeqFeature import SeqFeature
#from Bio import Entrez
#from httplib import HTTPException
#import subprocess
#import csv
#import fasta
#import string
#import sys
#import re
#import random
#import os
#import math
##from sets import Set
#from scipy import stats
#from collections import Counter
#from Bio.Alphabet import IUPAC
#from Bio import Entrez
#from Bio.Seq import Seq
#from Bio import SeqFeature
#from Bio import SeqRecord
#from Bio.SeqFeature import FeatureLocation
#from Bio import SeqIO
#from httplib import HTTPException

# <codecell>
##this takes the file ordered_SNP_pos.csv (which is NOT ordered -by chrom, then pos- coming out of SHM_Arabidopsis_Good_SNPs_to_csv_05-06-14.py, but csv ordered by hand)

ordered_snp_list = [[],[],[],[],[]]
with open("ordered_SNP_pos.csv", "rU") as o_snp:
    o_snp_reader = csv.reader(o_snp)
    for r in o_snp_reader:
        if r[0] == "Chr":
            continue
        ordered_snp_list[int(r[0])-1].append(r[1])
        

# <codecell>
gid_dict = {1:"240254421", 2:"240254678", 3:"240255695", 4:"240256243", 5:"240256493"} #GI ids for TAIR10

# <codecell>

#index_genbank_features adapted from http://www2.warwick.ac.uk/fac/sci/moac/people/students/peter_cock/python/genbank/#indexing_features
#
##http://www.insdc.org/documents/feature_table.html
##
def index_genbank_features(gb_record, feature_type, qualifier) :
    answer = dict()
    for (index, feature) in enumerate(gb_record.features) :
        if feature.type==feature_type :
            if qualifier in feature.qualifiers :
                for value in feature.qualifiers[qualifier] :
                    if value in answer :
                        answer[value].append(index)
                    else :
                        answer[value] = [index]
    return answer

cds_index = []
print "MAKE CDS FEATURE INDEX"
for c in gid_dict: # for each of the original chromomsome sequences
    print "Processing chrom", c
    filename = "GI" + str(gid_dict[c]) + ".gbk"
    print filename
    for gb_record in SeqIO.parse(open(filename, "r"), "gb"): #based on http://www2.warwick.ac.uk/fac/sci/moac/people/students/peter_cock/python/genbank/
        loc_tag_cds_dict = index_genbank_features(gb_record, "CDS", "locus_tag") #based on http://www2.warwick.ac.uk/fac/sci/moac/people/students/peter_cock/python/genbank/
        #print "loc_tag_cds_dict", loc_tag_cds_dict        
        cds_index.append(loc_tag_cds_dict)
        
print "len(cds_index)", len(cds_index)
chr_counter = 1
for i in cds_index:
    print len(i), "number of cds features on Chr", chr_counter
    chr_counter +=1

## <codecell>
####there are a few cds' in the GenBank files that produce errors when translated, catch them and exclude those cds' from further analysis 
print "LOOK FOR TRANSLATION ERRORS"

with open("Translation_errors_Col_Ref.csv", "wb") as error_file:
    error_writer = csv.writer(error_file)
    error_writer.writerow(["Chr", "GBfeature_index"])
    for c in gid_dict:
        count_temp = 0
        print "Processing chrom", c
        gb_file = "GI" + str(gid_dict[c]) + ".gbk"
        error_counter = 0
        for gb_record in SeqIO.parse(open(gb_file, "r"), "gb"):
            print "len(cds_index[c-1])", len(cds_index[c-1])  
            for key in cds_index[c-1]: #a list of dictionaries, with one dictionary per chromosome, key = AT_id, value = index 
                val = cds_index[c-1][key]
                for index in val:
                    count_temp +=1
                    if count_temp %1000 ==0:
                        print "count_temp",  count_temp
                        #print index
                    
                    
                    gb_feature = gb_record.features[index] 
                    straight_gb = gb_feature.extract(gb_record.seq)
                    try:
                        trans_gb = straight_gb.translate(cds=True)#Ok to force cds = True because we're testing wether the GB feature annotated as a cds, can be translated without error
                        str(gb_feature.qualifiers['translation'][0]) == str(trans_gb)
                    except:
                        error_counter +=1
                        error_writer.writerow([c, index]) #write the chromosome and index to file
                        continue
            print "Chromosome", c, ", number translation errors:", error_counter
                    
                                    
                    
# <codecell>
error_list = [[], [], [], [], []]

with open("Translation_errors_Col_Ref.csv", "rU") as errors:
    error_reader = csv.reader(errors)
    for r in error_reader:
        if r[0] == "Chr":
            continue
        error_list[int(r[0])-1].append(r[1])

print "error list", error_list
chr_counter = 1
for e in error_list:
    print len(e), "tranlation errors on Chr", chr_counter
    chr_counter +=1

###### <codecell>
######look up cds features in Cvi genbank file, translate cds and ask if it is the same as the "official" translation, if so locate the position of the mutation
print "LOOK FOR NON-SYN SNPs"
def trans_cds_feats(c, gb_file): 
    count_non_sym = 0
    for gb_record in SeqIO.parse(open(gb_file, "r"), "gb"):
        print "Processing ", gb_file
        find_file_prefix = re.search(".gbk", gb_file)
        file_prefix = gb_file[:find_file_prefix.start()]
        print file_prefix
        seq_out_file = file_prefix + "_trans_results.csv"
        with open(seq_out_file, "w") as out:
            out_writer = csv.writer(out)
            out_writer.writerow(["AT_id", "chr", "feat_start", "feat_end", "feat_strand", "protein_id", "product", "function", "translation", "gb_translation", "mutated_aa"])
            
            print "Number of CDS features ", len(cds_index[int(c)-1]), " Col Genome"
            
            new_geno_tag_cds_dict = index_genbank_features(gb_record, "CDS", "locus_tag") #based on http://www2.warwick.ac.uk/fac/sci/moac/people/students/peter_cock/python/genbank/
            #print "loc_tag_cds_dict", loc_tag_cds_dict        
            
            #new_geno_cds = 0
            #for f in gb_record.features:
            #    if f.type == "CDS":
            #        new_geno_cds += 1
            
            #print "Number of CDS features ", new_geno_cds, " SNP Genome"
            #assert len(cds_index[int(c)-1])  == new_geno_cds # make sure that there are the same number of features in the orginal (Col) Genbank record as there are in the new Cvi or Ler Genbank record

            print "Number of CDS features ", len(new_geno_tag_cds_dict), " SNP Genome"
            assert len(cds_index[int(c)-1])  == len(new_geno_tag_cds_dict) # make sure that there are the same number of features in the orginal (Col) Genbank record as there are in the new Cvi or Ler Genbank record


            for key in cds_index[int(c)-1]: #a list of dictionaries, with one dictionary per chromosome, key = AT_id, value = [index]               
                val = cds_index[int(c)-1][key]
                for index in val: #sometimes there is more than one record for AT id
                    if index in error_list:
                        continue
                    gb_feature = gb_record.features[index] 
                    straight_gb = gb_feature.extract(gb_record.seq)
                    
                    startposition = gb_feature.location.start.position
                    endposition = gb_feature.location.end.position
                    
                    strand = gb_feature.location.strand
                    #feat_position = (startposition, endposition, strand)

                    temp_snp_list = []#list of positions of SNPs in a feature
                    snp_count = 0 
                    for k in ordered_snp_list[c-1]:#ordered list of SNP pos
                        if int(k) in gb_feature:
                            temp_snp_list.append(int(k))
                            
                            snp_count +=1
                            mut_list = ""
                            trans_gb = straight_gb.translate()
                            find_stop = re.search("\*", str(trans_gb))
                            if find_stop == None:
                                mut_list = "NoStop"
                                out_writer.writerow([key, c, startposition, endposition, strand, str(gb_feature.qualifiers.get('protein_id', "")), str(gb_feature.qualifiers.get('product', "")), str(gb_feature.qualifiers.get('function', "")), str(trans_gb), str(gb_feature.qualifiers['translation'][0]), mut_list])            
                                break # we could try to track down the next stop codon and send that sequence into function prediction - but it would take some doing - TODO?
                            trans_seq = trans_gb[:find_stop.start()] #the translated sequence without the "*", so that we can ask if the sequence itself differs from the Genbank translation
                            print trans_seq   
                            print type(trans_seq)
                            if len(trans_seq) == 0:
                                break
                            if str(trans_seq)[0] !="M":
                                #print "NoStart"                                            
                                mut_list = "NoStart"
                                out_writer.writerow([key, c, startposition, endposition, strand, str(gb_feature.qualifiers.get('protein_id', "")), str(gb_feature.qualifiers.get('product', "")), str(gb_feature.qualifiers.get('function', "")), str(trans_seq), str(gb_feature.qualifiers['translation'][0]), mut_list])            
                                break #Potentially catostrophic SNP, any remaining SNPs in this featureare irrelevant
                            if len(str(trans_seq)) < len(str(gb_feature.qualifiers['translation'][0])):
                                #print "Pre-mature stop"
                                stop_at = len(str(trans_seq))-1
                                org_aa = str(gb_feature.qualifiers['translation'][0][stop_at])
                                mut_note = org_aa + str(stop_at) + "*" #Pre-mature stop codon                                   
                                mut_list = mut_note
                                out_writer.writerow([key, c, startposition, endposition, strand, str(gb_feature.qualifiers.get('protein_id', "")), str(gb_feature.qualifiers.get('product', "")), str(gb_feature.qualifiers.get('function', "")), str(trans_seq), str(gb_feature.qualifiers['translation'][0]), mut_list])            
                                break #Potentially catostrophic SNP, any remaining SNPs in this featureare irrelevant                                
                            if str(trans_seq) != str(gb_feature.qualifiers['translation'][0]):
                                for n in range(min(len(str(trans_seq)), len(str(gb_feature.qualifiers['translation'][0])))):
                                    print "n", n
                                    print "len(str(trans_seq))", len(str(trans_seq))
                                    print "len(str(gb_feature.qualifiers['translation'][0])", len(str(gb_feature.qualifiers['translation'][0]))
                                    if str(trans_seq[n]) != str(gb_feature.qualifiers['translation'][0][n]):
                                        aa_pos = n+1
                                        org_aa = str(gb_feature.qualifiers['translation'][0][n])
                                        cvi_aa = str(trans_seq[n])
                                        if aa_pos == 1: #if the first aa is not the same then the Cvi SNP interupted the start codon
                                            #print "NoStart"                                            
                                            mut_list = "NoStart"
                                            out_writer.writerow([key, c, startposition, endposition, strand, str(gb_feature.qualifiers.get('protein_id', "")), str(gb_feature.qualifiers.get('product', "")), str(gb_feature.qualifiers.get('function', "")), str(trans_seq), str(gb_feature.qualifiers['translation'][0]), mut_list])            
                                            break #Potentially catostrophic SNP, any remaining SNPs in this featureare irrelevant
                                        else:
                                            count_non_sym +=1
                                            mut_note = org_aa + str(aa_pos) + cvi_aa
                                            #print "mut_note", mut_note
                                            if len(mut_list) >1:
                                                mut_list = mut_list + "," + mut_note
                                            else:
                                                mut_list = mut_note

                                out_writer.writerow([key, c, startposition, endposition, strand,  str(gb_feature.qualifiers.get('protein_id', "")), str(gb_feature.qualifiers.get('product', "")), str(gb_feature.qualifiers.get('function', "")), str(trans_seq), str(gb_feature.qualifiers['translation'][0]), mut_list])            
                                break #break if it finds one snp in record - only need to know that *a* snp is in feature because sequence is already modified for ALL SNPs and looping over pos in dict
                            
    print "Number of non-synonymous SNPs for " + file_prefix + ": ", count_non_sym
    
for c in gid_dict:
    for prefix in ["Cvi_SNPs_Chr", "Ler_SNPs_Chr"]:
        gb_file = prefix + str(c) + ".gbk"
        trans_cds_feats(c, gb_file)


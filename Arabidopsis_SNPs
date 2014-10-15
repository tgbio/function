
############################
from Bio import SeqIO
#from Bio.Seq import Seq
from BCBio.GFF import GFFExaminer
from BCBio import GFF
#from Bio.Alphabet import IUPAC
import pprint
import csv
import re
import sys
import string
import os

from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
from Bio import SeqFeature
from Bio.SeqFeature import FeatureLocation
from Bio.SeqFeature import SeqFeature
from Bio import Entrez
from httplib import HTTPException
import subprocess
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
#####Retreive genbank records for each chromosome of TAIR10
gid_dict = {1:"240254421", 2:"240254678", 3:"240255695", 4:"240256243", 5:"240256493"} #GI ids for TAIR10
#gid_dict = {1: "238621221"}# , 2:"238621222", 3:"238621223", 4:"238621224", 5:"238481659"} #GI ids for TAIR9
#gid_dict = {2:"240254678"} #temporary for debugging

##for c in gid_dict:
###    #print "Processing records for chromosome", c
###    Entrez.email = "trudi.gulick001@umb.edu"
###    filename = "GI" + str(gid_dict[c]) + ".gbk"
###    try: #adapted from http://biopython.org/pipermail/biopython/2011-October/007555.html
###        if not os.path.isfile(filename):
###            net_handle = Entrez.efetch(db="nucleotide", rettype="gb", retmode="text", id=gid_dict[c])
###            out_handle = open(filename, "w")
###            out_handle.write(net_handle.read())
###            out_handle.close()
###            net_handle.close()
###            print "Saved gb record for chromosome", c
###                        
###                
###    except HTTPException, e:
###        print "Network problem: %s" % e
###        print "Second (and final) attempt..."
###        if not os.path.isfile(filename):
###            net_handle = Entrez.efetch(db="nucleotide", rettype="gb", retmode="text", id=gid_dict[c])
###            out_handle = open(filename, "w")
###            out_handle.write(net_handle.read())
###            out_handle.close()
###            net_handle.close()
###            print "Saved gb record for chromosome", c
##
##
#### <codecell>
#in_file = "./SNP/tair9_EckerCvi_LerHomozygous_snp tiny.gff"
#### <codecell>
##examiner = GFFExaminer()
##in_handle = open(in_file)
##pprint.pprint(examiner.available_limits(in_handle))
###ex = examiner.available_limits(in_handle)
###for g in ex['gff_id']:
###    print g, ex['gff_id'][g]
##
###in_handle.close()
###
#### <codecell>
###read the snps from the gff and separate them into a cvi and a ler dictionary 

tair_ver = "TAIR_9"
tair_ver = "TAIR_10"

if tair_ver == "TAIR_9":

#cvi_snps = [{}, {}, {}, {}, {}] #cvi_snps[chrom-1][(begin_pos, end_pos)] = [col_nt, cvi_nt]
#ler_snps = [{}, {}, {}, {}, {}] #ler_snps[chrom-1][(begin_pos, end_pos)] = [col_nt, ler_nt]
##chr_dict = [dict(gff_id=('Chr5',))]#[dict(gff_id=('Chr1',)), dict(gff_id=('Chr2',)), dict(gff_id=('Chr3',)), dict(gff_id=('Chr4',)), dict(gff_id=('Chr5',))]
#chr_dict = [dict(gff_id=["Chr1"]), dict(gff_id=["Chr2"]), dict(gff_id=["Chr3"]), dict(gff_id=["Chr4"]), dict(gff_id=["Chr5"])]
#
#
#for limit_info in chr_dict:
#    #print limit_info
#    for f in limit_info:
#        chrom = int(limit_info[f][0][-1])#will only work with single digit chromosomes
#        print "Reading the SNPs for Chr:", chrom, "from gff"
#    in_handle = open(in_file)
#    for rec in GFF.parse(in_handle, limit_info=limit_info):
#        for feature in rec.features:
#            
#            snp_begin = int(feature.location.start)
#            snp_end = int(feature.location.end)
#            assert snp_end-snp_begin==1 #check that the polymorphisms actually are actually *single* nt 
#            #snp_pos = snp_begin
#            #print snp_pos
#        
#           # print "feature.qualifiers[Note]", feature.qualifiers["Note"]
#        
#            col_nt = feature.qualifiers["Note"][2]
#            geno_nt = feature.qualifiers["Note"][3]
#            #print col_nt, geno_nt
#        
#            if feature.type == "Ecker_Cvi_snps":
#                cvi_snps[chrom-1][(snp_begin, snp_end)] = [col_nt, geno_nt]
#            elif feature.type == "Ler_homozygous_snp":
#                ler_snps[chrom-1][(snp_begin, snp_end)] = [col_nt, geno_nt]
#    in_handle.close()
###
###   # <codecell>
#####take the raw list of snps and keep only the positions where the col nt listed in the gff is the same as the nt at that position in Col genome
###    
#new_cvi_snps = [{}, {}, {}, {}, {}]    # new_cvi_snps = [{}, {}, {}, {}, {}] cvi_snps[chrom-1][pos] = [col_nt, cvi_nt]
#new_ler_snps = [{}, {}, {}, {}, {}]
#for c in gid_dict:
#    print "Processing chrom", c
#    filename = "GI" + str(gid_dict[c]) + ".gbk"
#    print filename
#    same_count1 = 0 
#    diff_count1 = 0
#
#    same_count2 = 0 
#    diff_count2 = 0
#
#    for gb_record in SeqIO.parse(open(filename, "r"), "gb"): #
#        chr_seq = gb_record.seq        
#        for p in cvi_snps[c-1]: #cvi_snps = [{}, {}, {}, {}, {}] cvi_snps[chrom-1][pos] = [col_nt, cvi_nt]
#            cvi_col_nt = cvi_snps[c-1][p][0]
#            if chr_seq[p[0]] == cvi_col_nt:
#                new_cvi_snps[c-1][p[0]] = cvi_snps[c-1][p]
#                same_count1 +=1
#            else:
#                diff_count1 +=1
#                
# 
#        print "Number positions same, cvi/reference", same_count1
#        print "Number positions different, cvi/reference", diff_count1
#
#        for ps in ler_snps[c-1]: #cvi_snps = [{}, {}, {}, {}, {}] cvi_snps[chrom-1][pos] = [col_nt, cvi_nt]           
#            ler_col_nt = ler_snps[c-1][ps][0]
#            if chr_seq[ps[0]] == ler_col_nt:
#                new_ler_snps[c-1][ps[0]] = ler_snps[c-1][ps]
#                same_count2 +=1
#            else:
#                diff_count2 +=1
#        print "Number positions same, ler/reference", same_count2
#        print "Number positions different, ler/reference", diff_count2
#        
#count_h = 1
#for h in new_cvi_snps:
#    print "New number of Cvi SNPs, Chr", count_h, ":", len(h)
#    count_h += 1 
#    
#count_k = 1
#for k in new_ler_snps:
#    print "New number of Ler SNPs, Chr", count_k, ":", len(k)
#    count_k += 1 
##
#####    #break
##### <codecell>            
##### take the ler and cvi snps (where Col nt at that position matches the nt in the genome, done above) and add only those postiions where ler and cvi are polymorpic relative to each other
#all_snps = [{}, {}, {}, {}, {}] #one dict per chromosome, keys = begining position if begin+1 = end, value = [col nt, cvi nt, ler nt]
#col_discrep = 0
#same_snp = 0
#for g in range(5):#TEMP for testing limit to 1, range should = 5
#    for p in new_ler_snps[g]: #ler/Cvi_snps[chrom-1][snp_pos] = [col_nt, geno_nt] where snp_pos is a the begining position of the SNP
#        if p in new_cvi_snps[g]: #SNP in both Ler and Cvi dicts
#            if new_ler_snps[g][p][1] != new_cvi_snps[g][p][1]: #if the Cvi and Ler nt are not the same
#                if new_ler_snps[g][p][0] == new_cvi_snps[g][p][0]: #is a position occurs is polymorphic in both genotypes, make sure the reported Col nt is the same in both dictionaries
#                    col_nt = new_ler_snps[g][p][0]
#                    all_snps[g][p] = [col_nt, new_cvi_snps[g][p][1], new_ler_snps[g][p][1]]
#                else:
#                    col_discrep +=1
#                    continue
#            else: #if Cvi and Ler nt are the same
#                print "Ler and Cvi have same SNP relative to Col at:", g+1, p
#                same_snp +=1
#        else: #SNP in Ler but not Cvi dict
#            col_nt = new_ler_snps[g][p][0]
#            all_snps[g][p] = [col_nt, col_nt, new_ler_snps[g][p][1]] #assuming that if snp not in a genotype the nt at that position is the reference nt (Col)
#    for q in new_cvi_snps[g]: #SNP in Cvi but not Ler dict
#        if q not in new_ler_snps[g]:
#            col_nt = new_cvi_snps[g][q][0]
#            all_snps[g][q] = [col_nt, new_cvi_snps[g][q][1], col_nt] #assuming that if snp not in a genotype the nt at that position is the reference nt (Col)
#
#count_useful_snps = 0
#for z in range(len(all_snps)):
#    count_useful_snps = count_useful_snps + len(all_snps[z])
#    
#
#print count_useful_snps, "total number of useful SNPs."  
#print same_snp, "sites are polymorphic relative to Columbia but not relative to each other."                
#print col_discrep, "positions in Col were reported as diff nts by the Ler and Cvi dicts."
#
# 
#
## <codecell>
###this takes a riduiculously long time to run, so am just appending positions then opening the resulting csv file and sorting
##generate an ordered list of the dictionary ids
#
##
#print "Writing useable SNPs to file, file will need to be ordered before continuing."
#with open("ordered_SNP_pos.csv", "w") as ordered_SNP:
#    o_snp_writer = csv.writer(ordered_SNP)
#    o_snp_writer.writerow(["Chr", "snp_pos", "Col nt", "Cvi nt", "Ler nt"])
#    for n in range(len(all_snps)): ##len(all_snps) = the number of chromsomes, all_snps[g][p] = [col_nt, new_cvi_snps[g][p][1], new_ler_snps[g][p][1]]
#        for j in all_snps[n]: # for key in dict
#            o_snp_writer.writerow([n+1, j, all_snps[n][j][0], all_snps[n][j][1], all_snps[n][j][2]]) #chromosome, position, col_nt, cvi_nt, ler_nt
#            
#print count_useful_snps, "writen to file."


# <codecell>


ordered_snp_list = [[],[],[],[],[]]
with open("ordered_SNP_pos tiny2.csv", "rU") as o_snp:
    o_snp_reader = csv.reader(o_snp)
    for r in o_snp_reader:
        if r[0] == "Chr":
            continue
        #print r[0], r[1]
        ordered_snp_list[int(r[0])-1].append(r[1])
        

   # <codecell>

#creating new genbank record adapted from http://www.biostars.org/p/57549/, #2-4 The addition of features to Genbank records wasn't working, features of Ler and Cvi GB file identical, only the sequences differ
#
#def make_gb_file(seq, gb_record, name, c):
#    new_record = SeqRecord(seq)
#    new_record.seq.alphabet = generic_dna
#    new_record.features = gb_record.features
#        
#    out_name = name + str(c) + ".gbk"
#    output_handle = open(out_name, "w")
#    SeqIO.write(new_record, output_handle, "gb")
#    print "File ", out_name, " created"
#    output_handle.close()
#
#for c in gid_dict:
#    print "Processing chrom", c
#    filename = "GI" + str(gid_dict[c]) + ".gbk"
#    count_ler_snps = 0
#    count_cvi_snps = 0
#    for gb_record in SeqIO.parse(open(filename, "r"), "gb"): 
#        print "original len(gb_record)", len(gb_record)
#        print "original len(gb_record.features)", len(gb_record.features)
#        print "original len(gb_record.annotations))", len(gb_record.annotations)
#        
#        with open("ordered_SNP_pos.csv", "rU") as o_snp:
#
#            o_snp_reader = csv.reader(o_snp)
#            prev_snp_pos = 0
#            for r in o_snp_reader:
#                if r[0] == "Chr":
#                    continue
#                chrom = int(r[0])
#                if chrom != c: #if SNP not on current chromosome "c"
#                    continue
#
#                pos = int(r[1])
#                #print pos
#                #print prev_snp_pos
#                #print gb_record.seq[prev_snp_pos:pos]
#                
#                #exit()
#                
#                col_nt = r[2]
#                assert len(col_nt) == 1  
#                cvi_nt = r[3]
#                assert len(cvi_nt) == 1
#                ler_nt = r[4]
#                assert len(ler_nt) == 1
#
#                assert gb_record.seq[pos] == col_nt #assert that the nucleotide in the Col genome matches what the SNP list says it should be
#                temp_seq = gb_record.seq[prev_snp_pos:pos] #this should be the SNP free segment of the genome between SNP positions                 
#                #print "len(temp_seq)", len(temp_seq)
#                if prev_snp_pos == 0:
#                    if col_nt != cvi_nt and col_nt != ler_nt: #cvi and ler are both polymorphic relative to col
#                        cvi_seq = temp_seq + cvi_nt
#                        ler_seq = temp_seq + ler_nt
#                        count_cvi_snps +=1
#                        count_ler_snps +=1
#                    elif col_nt != cvi_nt and col_nt == ler_nt: #cvi is polymorphic relative to col, but not ler
#                        cvi_seq = temp_seq + cvi_nt
#                        ler_seq = temp_seq + col_nt
#                        count_cvi_snps +=1
#                    elif col_nt == cvi_nt and col_nt != ler_nt: #ler is polymorphic relative to col, but not cvi
#                        cvi_seq = temp_seq + col_nt
#                        ler_seq = temp_seq + ler_nt 
#                        count_ler_snps +=1
#                    prev_snp_pos = pos+1
#
#                else:
#                    if col_nt != cvi_nt and col_nt != ler_nt: #cvi and ler are both polymorphic relative to col
#                        cvi_seq = cvi_seq + temp_seq + cvi_nt
#                        ler_seq = ler_seq + temp_seq + ler_nt
#                        count_cvi_snps +=1
#                        count_ler_snps +=1
#                    elif col_nt != cvi_nt and col_nt == ler_nt: #cvi is polymorphic relative to col, but not ler
#                        cvi_seq = cvi_seq + temp_seq + cvi_nt
#                        ler_seq = ler_seq + temp_seq + col_nt
#                        count_cvi_snps +=1
#                    elif col_nt == cvi_nt and col_nt != ler_nt: #ler is polymorphic relative to col, but not cvi
#                        cvi_seq = cvi_seq + temp_seq + col_nt
#                        ler_seq = ler_seq + temp_seq + ler_nt 
#                        count_ler_snps +=1
#                    prev_snp_pos = pos+1
#                
#                #print "ler_seq", ler_seq
#                print len(ler_seq)
#                
#                #print "cvi_seq", cvi_seq
#                print len(cvi_seq)
#                
#                #print "str(gb_record.seq[:pos])", str(gb_record.seq[:pos])
#                
#                assert len(ler_seq) == len(str(gb_record.seq[:pos]))+1
#                assert len(cvi_seq) == len(str(gb_record.seq[:pos]))+1 
#
#            cvi_seq = cvi_seq + str(gb_record.seq[prev_snp_pos:])
#            ler_seq = ler_seq + str(gb_record.seq[prev_snp_pos:])
#
#            assert len(cvi_seq) == len(ler_seq) == len(str(gb_record.seq))
#         
#        print "count_cvi_snps", count_cvi_snps
#        print "count_ler_snps", count_ler_snps         
#        
#        assert str(cvi_seq) != str(gb_record.seq) 
#        assert str(ler_seq) != str(gb_record.seq)
#        assert str(cvi_seq) != str(ler_seq)
#        
#        rec_dict = {"Cvi_SNPs_Chr":cvi_seq, "Ler_SNPs_Chr":ler_seq}
#        
#        for name in rec_dict:
#            seq = rec_dict[name]
#            make_gb_file(seq, gb_record, name, c)

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

for c in gid_dict: 
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

# <codecell>
###there are a few cds' in the GenBank files that produce errors when translated, catch them and exclude those cds' from further analysis 
#
#count_temp = 0
#with open("Translation_errors_Col_Ref.csv", "w") as error_file:
#    error_writer = csv.writer(error_file)
#    error_writer.writerow(["Chr", "GBfeature_index"])
#    for c in gid_dict:
#        print "Processing chrom", c
#        gb_file = "GI" + str(gid_dict[c]) + ".gbk"
#        error_counter = 0
#        for gb_record in SeqIO.parse(open(gb_file, "r"), "gb"):
#            print "len(cds_index[c-1])", len(cds_index[c-1])  
#            for key in cds_index[c-1]: #a list of dictionaries, with one dictionary per chromosome, key = AT_id, value = index 
#                val = cds_index[c-1][key]
#                for index in val:
#                    count_temp +=1
#                    if count_temp %1000 ==0:
#                        print "count_temp",  count_temp
#                        #print index
#                    
#                    
#                    gb_feature = gb_record.features[index] 
#                    straight_gb = gb_feature.extract(gb_record.seq)
#                    try:
#                        trans_gb = straight_gb.translate(cds=True)#Ok to force cds = True because we're testing wether the GB feature annotated as a cds, can be translated without error
#                        str(gb_feature.qualifiers['translation'][0]) == str(trans_gb)
#                    except:
#                        error_counter +=1
#                        error_writer.writerow([c, index]) #write the chromosome and index to file
#                        continue
#            print "Chromosome", c, ", number translation errors:", error_counter
#                    
                                    
                    
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
#
##### <codecell>
#####look up cds features in Cvi genbank file, translate cds and ask if it is the same as the "official" translation, if so locate the position of the mutation
####
#####count_temp = 0
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
            
            print "Number of CDS features", len(cds_index[int(c)-1])  

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
                                break
                            trans_seq = trans_gb[:find_stop.start()] #the translated sequence without the "*", so that we can ask if the sequence itself differs from the Genbank translation
                            if trans_seq[0] !="M":
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
#
##
## <codecell>
#####define regions
#roi = {1:[[174952, 3112850], [23698501, 24905939]], 4:[[11438290, 16924197]]}
#
## <codecell>
#genotype_trans = {1:["Cvi_SNPs_Chr1_trans_results.csv", "Ler_SNPs_Chr1_trans_results.csv"]}
#for i in genotype_trans:
#    print "Chr", i
#    for j in genotype_trans[i]:
#        
#        with open(j, "r") as trans_file:
#            trans_reader = csv.reader(trans_file)
#            find_geno = re.search("_", j)
#            geno = j[:find_geno.start()]
#            print geno
#            out_filename = "ROI_seqs_" + geno + "_Chr" + str(i) + ".csv"
#            with open(out_filename, "w") as result_file:
#                result_writer = csv.writer(result_file)
#                result_writer.writerow(["Chr", "QTL_region_begin", "QTL_region_end", "AT_id", "feat_start", "feat_end", "feat_strand", "protein_id", "product", "function", "translation", "gb_translation", "mutated_aa"])
#                
#                for r in trans_reader:
#                    if r[0] == "AT_id":
#                        continue
#                    #print "r[1]",r[1]
#                    if int(r[1]) in roi: #If chromosome from translation file is a key in the dict roi, then there is a QTL region on the chrom
#                        for reg in roi[int(r[1])]:
#
#                            if ((int(r[2]) or int(r[3])) <= reg[1]) and ((int(r[2]) or int(r[3])) >= reg[0]): #if the begining of the feature(r[2]) is greater than the begining of the QTL regionnd less than the end or the end is greater than 
#                                result_writer.writerow([str(i), str(reg[0]), str(reg[1]), r[0], r[2], r[3], r[4], r[5], r[6], r[7], r[8], r[9], r[10]])
## <codecell>
##this section takes the roi results split by genotype and creates a dictionary with value Ler, Cvi and value = dictionary with keys = regions values = [(NP id, gb_translation sequence, mutations), (...), ...]
##this section is only useful if the SIFT feature of comma separated RefSeq Id, mutations were working, its not
#roi_files = {1:["ROI_seqs_Cvi_Chr1.csv", "ROI_seqs_Ler_Chr1.csv"]}
#for chr in roi_files:
#    geno_dict = {}
#    for f in roi_files[chr]:
#        find_start_geno = re.search("ROI_seqs_", f) # look up genotype
#        find_end_geno = re.search("_Chr", f)
#        geno = f[find_start_geno.end():find_end_geno.start()]
#        print geno
#        
#        if geno not in geno_dict:
#            geno_dict[geno] = {}
#        with open(f, "r") as roi_results:
#            roi_reader = csv.reader(roi_results)
#            for l in roi_reader:
#                if l[0] == "Chr":
#                    continue
#                if (l[0], l[1], l[2]) in geno_dict[geno]:
#                    np_name = string.strip(l[7], "['']")
#                    snp_list = (np_name, l[11], l[12])
#                    #print snp_list
#                    geno_dict[geno][(l[0], l[1], l[2])].append(snp_list)
#                else:                                                   ####2-7 check that its looping/setting correctly through here
#                    np_name = string.strip(l[7], "['']")
#                    snp_list = (np_name, l[11], l[12])
#                    #print snp_list
#                    geno_dict[geno][(l[0], l[1], l[2])] = [snp_list]
## <codecell>
#print geno_dict
## <codecell>
##make sure there are no duplicate entries
#for g in geno_dict:
#    for r in geno_dict[g]:
#        
#        print len(geno_dict[g][r])
#        s = [f for f in set(geno_dict[g][r])]
#        print len(s)
## <codecell>       
#os.chdir("/Users/trudigulick/Dropbox/Stuart/Arabidopsis/Working")
## <codecell>
#for g in geno_dict: # g = geno type
#    #print "begining loop"
#    #os.getcwd()
#    geno = g 
#    print geno
#    geno_path = "./" + geno
#    if os.path.exists(geno_path) == False:
#        os.mkdir(geno_path)
#    os.chdir(geno_path)
#    #print "changed into geno dir"
#    #os.getcwd()
#    for reg in geno_dict[g]: # geno_dict[g] is a dictionary, with key = region 
#        reg_descrip = "Chr_" + reg[0] + "_Reg_" + reg[1] + "_" + reg[2]
#        reg_dir_name = "./" + reg_descrip
#        if os.path.exists(reg_dir_name) == False:
#            os.mkdir(reg_dir_name)
#        os.chdir(reg_dir_name)
#        #print "changed into region dir"
#        #os.getcwd()
#        for np in geno_dict[g][reg]: #at is a list where the first element is the RefSeq id and the other elements are amino acid mutaions             
#            if  "*" in np[2]:
#                print "**** ", geno, reg, np[0], " premature stop ", np[2], " ****"
#            elif "NoStart" in np[2]:
#                print "**** ", geno, reg, np[0], " no start ", np[2], " ****"
#            elif "NoStop" in np[2]:
#                print "**** ", geno, reg, np[0], " no stop ", np[2], " ****"
#            else:
#                fasta_file = geno + "_" + str(reg[0]) + "_" + str(reg[1])+ "_" + str(reg[2]) + "_" + np[0] + ".fasta"
#                mut_file = geno + "_" + str(reg[0]) + "_" + str(reg[1])+ "_" + str(reg[2]) + "_" + np[0] + "_aa_subs.txt"
#                with open(fasta_file, "w") as fasta:
#                    header = ">" + np[0] + '\n'          
#                    fasta.write(header)
#                    fasta.write(np[1])
#                    with open(mut_file, "w") as subs:
#                        subs.write(np[2])
#        #print "end region loop"
#        os.chdir("..")
#        #os.getcwd()
#    #print "end geno loop"
#    os.chdir("/Users/trudigulick/Dropbox/Stuart/Arabidopsis/Working")
#    #os.getcwd()
 
## <codecell>
##this section takes the roi results split by genotype and creates a dictionary with value Ler, Cvi and value = dictionary with keys = regions values = [(NP id, mutations, mutations, .. ), (...), ...]
##this section is only useful if the SIFT feature of comma separated RefSeq Id, mutations were working, its not
#roi_files = {1:["ROI_seqs_Cvi_Chr1.csv", "ROI_seqs_Ler_Chr1.csv"]}
#for chr in roi_files:
#    geno_dict = {}
#    for f in roi_files[chr]:
#        find_start_geno = re.search("ROI_seqs_", f) # look up genotype
#        find_end_geno = re.search("_Chr", f)
#        geno = f[find_start_geno.end():find_end_geno.start()]
#        print geno
#        if geno not in geno_dict:
#            geno_dict[geno] = {}
#        with open(f, "r") as roi_results:
#            roi_reader = csv.reader(roi_results)
#            for l in roi_reader:
#                if l[0] == "Chr":
#                    continue
#                if (l[0], l[1], l[2]) in geno_dict[geno]:
#                    #snp_list = [string.strip(sp, "['']") for sp in l[12]]
#                    snp_list = l[12].split(",")
#                    #print snp_list
#                    #print type(snp_list)
#                    snp_list.insert(0, string.strip(l[7], "['']"))
#                    snp_list = tuple(snp_list)
#                    geno_dict[geno][(l[0], l[1], l[2])].append(snp_list)
#                else:                                                   ####2-7 check that its looping/setting correctly through here
#                    snp_list = l[12].split(",")
#                    snp_list.insert(0, string.strip(l[7], "['']"))
#                    snp_list = tuple(snp_list)
#                    geno_dict[geno][(l[0], l[1], l[2])] = [snp_list]
#
##make sure there are no duplicate entries
#for g in geno_dict:
#    for r in geno_dict[g]:
#        
#        print len(geno_dict[g][r])
#        s = [f for f in set(geno_dict[g][r])]
#        print len(s)
#
#for g in geno_dict: # g = geno type
#    geno = g 
#    print geno
#    for reg in geno_dict[g]: # geno_dict[g] is a dictionary, with key = region 
#        print reg        
#        AT_id_file = geno + "_" + str(reg[0]) + "_" + str(reg[1])+ "_" + str(reg[2]) + "_ATids_muts.csv"
#        with open(AT_id_file, "w") as AT_by_region:
#            AT_reg_writer = csv.writer(AT_by_region)            
#            for at in geno_dict[g][reg]: #at is a list where the first element is the RefSeq id and the other elements are amino acid mutaions  
#                if  "*" in at[1]:
#                    print at
#                elif "NoStart" in at[1]:
#                    print at
#                elif "NoStop" in at[1]:
#                    print at
#                else:
#                    AT_reg_writer.writerow(at)
#
## <codecell>

####loop over the trans files for each genotype, if feature is in range add it to a dictionary with key = ATid
####print number of AT ids in range
####print number of items in AT values
###collapse the records with non-synomouse SNPs if appropriate (ie where there are alternative mRNA's associated one AT seq, where the alternative slicing doesn't affect the protein encoded)
##read the records into a dictionary with key = AT id, value = list of lists containing chr, position, NP id,product, translation of genotype seq, official translation, mutations
##for items in dictionary, for items in value list, retain only the unique values
##do the same for the other genotype
##write the unique records of each to file
#
## <codecell>
##define regions
##pull out records within that region, counting them 
##write the regions to file

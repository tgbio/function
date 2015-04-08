############################
from Bio import SeqIO
from BCBio import GFF
import csv
import os
from Bio import Entrez
from httplib import HTTPException


# <codecell>
#####Retreive genbank records for each chromosome of TAIR10

gid_dict = {1:"240254421", 2:"240254678", 3:"240255695", 4:"240256243", 5:"240256493"} #GI ids for TAIR10


for c in gid_dict:
    #print "Processing records for chromosome", c
    Entrez.email = "trudi.gulick001@umb.edu"
    filename = "GI" + str(gid_dict[c]) + ".gbk"
    try: #adapted from http://biopython.org/pipermail/biopython/2011-October/007555.html
        if not os.path.isfile(filename):
            net_handle = Entrez.efetch(db="nucleotide", rettype="gb", retmode="text", id=gid_dict[c])
            out_handle = open(filename, "w")
            out_handle.write(net_handle.read())
            out_handle.close()
            net_handle.close()
            print "Saved gb record for chromosome", c
                        
                
    except HTTPException, e:
        print "Network problem: %s" % e
        print "Second (and final) attempt..."
        if not os.path.isfile(filename):
            net_handle = Entrez.efetch(db="nucleotide", rettype="gb", retmode="text", id=gid_dict[c])
            out_handle = open(filename, "w")
            out_handle.write(net_handle.read())
            out_handle.close()
            net_handle.close()
            print "Saved gb record for chromosome", c

# <codecell>
#in_file = "tair9_EckerCvi_LerHomozygous_snp.gff"
#
#examiner = GFFExaminer()
#in_handle = open(in_file)
#ex = examiner.available_limits(in_handle)
#for g in ex['gff_id']:
##    print g, ex['gff_id'][g]
#
#in_handle.close()
# <codecell>

##read the snps from the gff and separate them into a Cvi and a Ler dictionary 
in_file = "tair9_EckerCvi_LerHomozygous_snp.gff"
cvi_snps = [{}, {}, {}, {}, {}] #cvi_snps[chrom-1][(begin_pos, end_pos)] = [col_nt, cvi_nt]
ler_snps = [{}, {}, {}, {}, {}] #ler_snps[chrom-1][(begin_pos, end_pos)] = [col_nt, ler_nt]
chr_dict = [dict(gff_id=["Chr1"]), dict(gff_id=["Chr2"]), dict(gff_id=["Chr3"]), dict(gff_id=["Chr4"]), dict(gff_id=["Chr5"])]


for limit_info in chr_dict:
    for f in limit_info:
        chrom = int(limit_info[f][0][-1])#will only work with single digit chromosomes
        print "Reading the SNPs for Chr:", chrom, "from gff"
    in_handle = open(in_file)
    for rec in GFF.parse(in_handle, limit_info=limit_info):
        for feature in rec.features: #feature.qualifiers formatted: {'Note': ['Cvi', '1', 'T', 'A'], 'source': ['EckerSNPs'], 'type': ['substitution'], 'Name': ['Ecker_325451']}
            snp_begin = int(feature.location.start) # zero based counting
            snp_end = int(feature.location.end) # zero based counting
            assert snp_end-snp_begin==1 #check that the polymorphisms actually are actually *single* nt 
        
            col_nt = feature.qualifiers["Note"][2]
            geno_nt = feature.qualifiers["Note"][3]
        
            if feature.type == "Ecker_Cvi_snps":
                cvi_snps[chrom-1][(snp_begin, snp_end)] = [col_nt, geno_nt]
            elif feature.type == "Ler_homozygous_snp":
                ler_snps[chrom-1][(snp_begin, snp_end)] = [col_nt, geno_nt]
    in_handle.close()

# <codecell>
######take the raw list of snps and keep only the positions where the col nt listed in the gff is the same as the nt at that position in Col genome
   
new_cvi_snps = [{}, {}, {}, {}, {}]         
new_ler_snps = [{}, {}, {}, {}, {}]
for c in gid_dict:
    print "Processing chrom", c
    filename = "GI" + str(gid_dict[c]) + ".gbk"
    print filename
    same_count1 = 0 
    diff_count1 = 0
    
    same_count2 = 0 
    diff_count2 = 0
        
    for gb_record in SeqIO.parse(open(filename, "r"), "gb"): #
        chr_seq = gb_record.seq        
        for p in cvi_snps[c-1]: #cvi_snps = [{}, {}, {}, {}, {}] cvi_snps[chrom-1][pos] = [col_nt, cvi_nt]
            cvi_col_nt = cvi_snps[c-1][p][0] #Columbia nt at position of Cvi SNP
            if chr_seq[p[0]] == cvi_col_nt: #if the genome sequence for the same chromosome is the same as the Col nt
                new_cvi_snps[c-1][p[0]] = cvi_snps[c-1][p]
                same_count1 +=1
            else:
                diff_count1 +=1
                

        print "Number positions same, cvi/reference", same_count1
        print "Number positions different, cvi/reference", diff_count1

        for ps in ler_snps[c-1]: #cvi_snps = [{}, {}, {}, {}, {}] cvi_snps[chrom-1][pos] = [col_nt, cvi_nt]           
            ler_col_nt = ler_snps[c-1][ps][0]
            if chr_seq[ps[0]] == ler_col_nt:
                new_ler_snps[c-1][ps[0]] = ler_snps[c-1][ps]
                same_count2 +=1
            else:
                diff_count2 +=1
                

        print "Number positions same, ler/reference", same_count2
        print "Number positions different, ler/reference", diff_count2
    
   
count_h = 1
for h in new_cvi_snps:
    print "New number of Cvi SNPs, Chr", count_h, ":", len(h)
    count_h += 1 

count_k = 1
for k in new_ler_snps:
    print "New number of Ler SNPs, Chr", count_k, ":", len(k)
    count_k += 1 
# <codecell>            
###### take the ler and cvi snps (where Col nt at that position matches the nt in the genome, done above) and add only those postiions where ler and cvi are polymorpic relative to each other
all_snps = [{}, {}, {}, {}, {}] #one dict per chromosome, keys = begining position if begin+1 = end, value = [col nt, cvi nt, ler nt]
col_discrep = 0
same_snp = 0
for g in range(5):#TEMP for testing limit to 1, range should = 5
    for p in new_ler_snps[g]: #ler/Cvi_snps[chrom-1][snp_pos] = [col_nt, geno_nt] where snp_pos is a the begining position of the SNP
        if p in new_cvi_snps[g]: #SNP (position) in both Ler and Cvi dicts
            if new_ler_snps[g][p][1] != new_cvi_snps[g][p][1]: #if the Cvi and Ler nt are not the same
                if new_ler_snps[g][p][0] == new_cvi_snps[g][p][0]: #is a position occurs is polymorphic in both genotypes, make sure the reported Col nt is the same in both dictionaries
                    col_nt = new_ler_snps[g][p][0]
                    all_snps[g][p] = [col_nt, new_cvi_snps[g][p][1], new_ler_snps[g][p][1]]
                else:
                    col_discrep +=1
                    continue
            else: #if Cvi and Ler nt are the same
                print "Ler and Cvi have same SNP relative to Col at:", g+1, p
                same_snp +=1
        else: #SNP in Ler but not Cvi dict (Cvi nt is same as Col)
            col_nt = new_ler_snps[g][p][0]
            all_snps[g][p] = [col_nt, col_nt, new_ler_snps[g][p][1]] #assuming that if SNP not in a genotype the nt at that position is the reference nt (Col)
    for q in new_cvi_snps[g]: #SNP in Cvi but not Ler dict
        if q not in new_ler_snps[g]:
            col_nt = new_cvi_snps[g][q][0]
            all_snps[g][q] = [col_nt, new_cvi_snps[g][q][1], col_nt] #assuming that if SNP not in a genotype the nt at that position is the reference nt (Col)

count_useful_snps = 0
for z in range(len(all_snps)):
    count_useful_snps = count_useful_snps + len(all_snps[z])
    

print count_useful_snps, "total number of useful SNPs."  
print same_snp, "sites are polymorphic relative to Columbia but not relative to each other."                
print col_discrep, "positions in Col were reported as diff nts by the Ler and Cvi dicts."

 

# <codecell>
#write the useable SNP locations to a csv, entries need to be sorted by Chromosome, with SNP positions in assending order

print "Writing useable SNPs to file, file will need to be ordered before continuing."
with open("ordered_SNP_pos.csv", "w") as ordered_SNP:
    o_snp_writer = csv.writer(ordered_SNP)
    o_snp_writer.writerow(["Chr", "snp_pos", "Col nt", "Cvi nt", "Ler nt"])
    for n in range(len(all_snps)): ##len(all_snps) = the number of chromsomes, all_snps[g][p] = [col_nt, new_cvi_snps[g][p][1], new_ler_snps[g][p][1]]
        for j in all_snps[n]: # for key in dict
            o_snp_writer.writerow([n+1, j, all_snps[n][j][0], all_snps[n][j][1], all_snps[n][j][2]]) #chromosome, position, col_nt, cvi_nt, ler_nt
            
print count_useful_snps, "written to file."

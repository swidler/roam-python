#!/usr/bin/python3

import amsample as a
import mmsample as m
import tools as t

#create Mmsample object
mms = m.Mmsample()
mms.parse_infile("bone5.txt")
outfile = "objects/bone5"
t.save_object(outfile, mms)
#create object
ams = a.Amsample(name="I1116", abbrev="1116")                                                                      
#ams = a.Amsample(name="Ust_Ishim", abbrev="Ust")                                                                   
#ams = a.Amsample(name="Altai_Neanderthal", abbrev="Alt")
#mat = a.Amsample()                                                                                                 
#ams2 = a.Amsample(name="First", coord_per_position=[2,2,2,2,2,2,2], no_t=[1,2,3], chr_names=["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7"])                                                                           

#vars for bam_to_am
chr_lengths = [249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,141213431,135534747,135006516,133851895,115169878,107349540,102531392,90354753,81195210,78077248,59128983,63025520,48129895,51304566,155270560,59373566]
#ust_chr_lengths = [249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022, 141213431, 135534747, 135006516, 133851895, 115169878, 107349540, 102531392, 90354753, 81195210, 78077248, 59128983, 63025520, 48129895, 51304566, 155270560, 59373566, 16569]

#create object from bam file
ams.bam_to_am(filename="/slow/nomi/I1116_chr8.bam", library="double", chr_lengths=chr_lengths, genome_seq="/slow/nomi/hg19.fa.gz", species="Homo sapiens", trim_ends=True, chroms=[7])
#ams.bam_to_am(filename="../../I1116.bam", library="double", chr_lengths=chr_lengths, genome_seq="../../hg19.fa.gz", species="Homo sapiens", trim_ends=True)
#ams.bam_to_am(filename="../../ust_ishim.bam", library="single", chr_lengths=ust_chr_lengths, genome_seq="../../hg19.fa.gz", species="Homo sapiens")
#ams.bam_to_am(filename="../../Alt_chr22.bam", library="single", chr_lengths=chr_lengths, genome_seq="../../hg19.fa.gz", species="Homo sapiens", chroms=[21], tot_chroms=[21])
#ams.bam_to_am(filedir="../../altai", library="single", chr_lengths=chr_lengths, genome_seq="../../hg19.fa.gz", species="Homo sapiens", chroms=[20, 21], tot_chroms=[20, 21])
#ams.bam_to_am(filedir="../../altai", file_per_chrom=True, library="single", chr_lengths=chr_lengths, genome_seq="../../hg19.fa.gz", species="Homo sapiens")                                                                              

#get object info from text file
#ams.parse_infile("data/python_dumps/Altai_Neanderthal_bam.txt")                                       
#ams.parse_infile("data/python_dumps/Altai_Neanderthal_diag2_from_matlab.txt")                                       
#ams.parse_infile("data/matlab_dumps/altai.txt")                                       
#ams.parse_infile("../../u_1116.txt")                                                                                
#ams.parse_infile("data/python_dumps/I1116_drate.txt")
#ams.parse_infile("data/python_dumps/I1116_meth.txt")                                                       

#input file
#infile = "objects/U1116"
#infile = "objects/U1116_diag"
#infile = "objects/U1116_filtered"                                                                                   
#infile = "objects/bone_5"

#load object from (pickled) input file
#ams = t.load_object(infile)                                                                                         
#mms = t.load_object(infile)

#run with profiling
#import cProfile                                                                                    
#cProfile.run("ams = t.load_object(infile)", "data/logs/load_profile")
#cProfile.run("ams.diagnose()", "data/logs/amsample_profile")                                                        
#cProfile.run("t.save_object(outfile, ams)", "data/logs/save_profile")
#cProfile.run("ams.bam_to_am(library='double', chr_lengths=chr_lengths, genome_seq='../../hg19.fa.gz', species='Homo sapiens')", "data/logs/bam_profile")

#run various methods
#ams.diagnose()                                                                                                      
#ams.filter()
#ams.estimate_drate(ref=mms)                                                                                         
#ams.reconstruct_methylation()
#ams.simulate(0.018598, mms) #rate comes from drate global in I1116_meth.txt
#ams.simulate(0, mms) #rate for testing

#get and print attributes
#name = ams.name
#print(f"name: {name}")
#num = ams.no_chrs
#print(f"num of chroms: {num}")
#meth_py = ams.methylation["methylation"]  
#meth_mat = mat.methylation["methylation"]
#for chrom in range(ams.no_chrs):                                                                                    
#    py = np.array(meth_py[chrom])
#    diffs = py - meth_mat[chrom]                                                                                    
#    diffs_nonan = [abs(x) for x in diffs if ~np.isnan(x)]
#    res = max(diffs_nonan)                                                                                          
#    print(f"Max methylation diff for {chrom}: {res}")
#base = ams.get_base_no("chr2", "c")
#print(f"base: {base[0:20]}")                                                                                        
#(no_t, no_ct) = ams.smooth("chr5", 17)                                                                              
#print(f"no_t: {no_t[0:20]}")
#print(f"no_ct: {no_ct[0:25]}")


#print object
#print(ams)
#print(ams2)

#output file
outfile = "objects/U1116_8_2"
#outfile = "objects/U1116_filtered_drate"                                                                            

#save (pickle) object to output file
t.save_object(outfile, ams)                                                                            
#t.save_object(outfile, ams)                                                                                         

#roam stage
#stage = "bam"
#stage = "drate"                                                                                                     
#stage = "meth"                                                                                                      
#stage = "sim_0"
#stage = "filt2_from_matlab"

#dump object to text file
#ams.dump(stage)

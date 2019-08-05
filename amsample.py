#!/usr/bin/python3

import tools as t
import numpy as np
import math
#import matplotlib.pyplot as plt
#import matplotlib.mlab as mlab
import scipy.stats as stats
import copy
from chroms import Chrom
import pysam
from Bio import SeqIO
import gzip
import datetime

class Amsample(Chrom):
    def __init__(self, name="unknown", abbrev="unk", species="unknown", reference="", library="", chr_names=[], coord_per_position=[], no_a = [], no_c = [], no_g = [], no_t = [], g_to_a = [], c_to_t = [], diagnostics = {}, p_filters = {}, is_filtered = False, is_simulated = False, methylation={}, d_rate = {}, metadata=[]):
        self.name = name
        self.abbrev = abbrev
        self.species = species
        self.reference = reference
        self.library = library
        self.chr_names = chr_names
        self.coord_per_position = coord_per_position
        self.no_a = no_a
        self.no_c = no_c
        self.no_g = no_g
        self.no_t = no_t
        self.g_to_a = g_to_a
        self.c_to_t = c_to_t
        self.diagnostics = diagnostics
        self.p_filters = p_filters
        self.is_filtered = is_filtered
        self.is_simulated = is_simulated
        self.methylation = methylation
        self.d_rate = d_rate
        self.metadata = metadata
        self.no_chrs = len(coord_per_position)

    def __repr__(self): #defines print of object
        return "name: %s\nabbrev: %s\nspecies: %s\nreference: %s\nlibrary: %s\nchr_names: %s\ncoord_per_position: %s\nno_a: %s\nno_c: %s\nno_g: %s\nno_t: %s\ng_to_a: %s \nc_to_t: %s\ndiagnostics: %s\np_filters: %s\nis_filtered: %s\nis_simulated: %s\nmethylation: %s\nd_rate: %s\nmetadata: %s\nno_chrs: %s" % (self.name, self.abbrev, self.species, self.reference, self.library, self.chr_names, self.coord_per_position, self.no_a, self.no_c, self.no_g, self.no_t, self.g_to_a, self.c_to_t, self.diagnostics, self.p_filters, self.is_filtered, self.is_simulated, self.methylation, self.d_rate, self.metadata, self.no_chrs)

    def diag_eq(self, other):
        if len(self.diagnostics) == 0 and len(other.diagnostics) == 0:
            return True
        else:
            return(self.diagnostics["effective_coverage"] == other.diagnostics["effective_coverage"])
    
    def filt_eq(self, other):
        if len(self.p_filters) == 0 and len(other.p_filters) == 0:
            return True
        else:
            return(self.p_filters["method"] == other.p_filters["method"] and list(self.p_filters["max_coverage"]) == list(other.p_filters["max_coverage"]) and self.p_filters["max_TsPerCoverage"] == other.p_filters["max_TsPerCoverage"] and list(self.p_filters["max_g_to_a"]) == list(other.p_filters["max_g_to_a"]) and list(self.p_filters["max_a"]) == list(other.p_filters["max_a"]))
    
    def meth_eq(self, other):
        if len(self.methylation) == 0 and len(other.methylation) == 0:
            return True
        else:
            return(self.methylation == other.methylation) #need to add keys once populated

    def drate_eq(self, other):
        if len(self.d_rate) == 0 and len(other.d_rate == 0):
            return True
        else:
            return(self.d_rate["global"] == other.d_rate["global"] and self.d_rate["dglobal"] == other.d_rate["dglobal"] and list(self.d_rate["local"]) == list(other.d_rate["local"]) and list(self.d_rate["dlocal"]) == list(other.d_rate["dlocal"]) and list(self.d_rate["no_positions"]) == list(other.d_rate["no_positions"]))



    def __eq__(self, other): #overload == operator
        return(self.name == other.name and self.abbrev == other.abbrev and self.species == other.species and
                self.reference == other.reference and self.library == other.library and 
                self.chr_names == other.chr_names and self.coord_per_position == other.coord_per_position and
                self.no_a == other.no_a and self.no_c == other.no_c and self.no_g == other.no_g and
                self.no_t == other.no_t and self.g_to_a == other.g_to_a and self.c_to_t == other.c_to_t and
                self.diag_eq(other) and self.filt_eq(other) and self.meth_eq(other) and self.drate_eq(other) and
                self.is_filtered == other.is_filtered and self.is_simulated == other.is_simulated and
                self.metadata == other.metadata and self.no_chrs == other.no_chrs)

    @staticmethod
    def find_indel(cig, seq, qual):
        types = {0:"M", 1:"I", 2:"D", 3:"N", 4:"S"}
        start = 0
        cig_dict = {}
        for element in cig:
            cat = types[element[0]]
            end = start + element[1]
            cig_dict[start] = str(end) + cat
            start = end
        keys = list(cig_dict.keys())
        keys.sort(reverse=True)
        for key in keys:
            val = cig_dict[key]
            cat = val[-1]
            val = int(val[:-1])
            if cat == "M" or cat == "S":  # what should be done with soft clipped?
                continue
            elif cat == "I":
                seq = seq[:key]+seq[val:]
                del qual[key:val]
            elif cat == "D" or cat == "N":
                size = val - key
                seq = seq[:key] + "0"*size + seq[key:]
                qual[key:key] = [0 for x in range(size)]
        return(seq, qual)


    def bam_to_am(self, filename="../../I1116.bam", library=None, chr_lengths=None, genome_seq=None, org=None, species=None, chr_format="num", chroms=list(range(23)), save_aligned=False, folder="", tot_chroms=list(range(23)), trim_ends=False, mapq_thresh=20, qual_thresh=20, asian_african=""):  # genome_seq is filename, orig qual_thresh=53, subtract 33 for 0 start
        corrupt_reads = 0
        low_mapq_reads = 0
        reads_long_cigar = 0
        long_cigar_not_i_d = 0
        xflag = 0
        if len(tot_chroms) == len(chroms):
            self.no_t = [[]]*(max(chroms)+1)
            self.no_c = [[]]*(max(chroms)+1)
            self.no_g = [[]]*(max(chroms)+1)
            self.no_a = [[]]*(max(chroms)+1)
            self.c_to_t = [[]]*(max(chroms)+1)
            self.g_to_a = [[]]*(max(chroms)+1)

        bam = pysam.AlignmentFile(filename, "rb")
        chrom_names = bam.references[0:max(tot_chroms)+1]  # lists names of all chroms up to max present in num order
        with gzip.open(genome_seq, "rt") as fas:
            records = list(SeqIO.parse(fas, "fasta"))
                #for record in SeqIO.parse(fas, "fasta"):
        for chrom in tot_chroms:
            chrom_name = chrom_names[chrom]
            date = datetime.datetime.now()
            date = date.strftime("%c")
            print(f"Starting chromosome {chrom_name}: {date}")
            chrom_bam = bam.fetch(chrom_name)
            tot_a_plus = np.zeros(chr_lengths[chrom])
            tot_a_minus = np.zeros(chr_lengths[chrom])
            tot_g_plus = np.zeros(chr_lengths[chrom])
            tot_g_minus = np.zeros(chr_lengths[chrom])
            tot_c_plus = np.zeros(chr_lengths[chrom])
            tot_c_minus = np.zeros(chr_lengths[chrom])
            tot_t_plus = np.zeros(chr_lengths[chrom])
            tot_t_minus = np.zeros(chr_lengths[chrom])
            #write aligned reads to file?
            for record in records:  # these are from the genome fasta file
                num, seq = record.id, str(record.seq)
                if num == "chr"+chrom_name:
                    seq = seq.upper()  # change all letters to uppercase
                    c_in_genome = [1 if x == "C" else 0 for x in seq]
                    c_in_genome[-1] = 0  # c in last pos can't be cpg
                    g_in_genome = [1 if x == "G" else 0 for x in seq]
            c_in_cpg = np.zeros(chr_lengths[chrom])
            g_in_cpg = np.zeros(chr_lengths[chrom])
            c_idx = []
            g_idx = []
            for i in range(len(c_in_genome)):
                if c_in_genome[i] and g_in_genome[i+1]:
                    c_in_cpg[i] = 1
                    g_in_cpg[i+1] = 1
                    c_idx.append(i)
                    g_idx.append(i+1)
            if asian_african:
                mutations = [x[2] for x in asian_african if x[0] == "chr"+chrom_name]  # start coord, = 0-based point mutation pos
                to_remove = []
                for i in c_idx:
                    if i in mutations:
                        to_remove.append(i)
                        to_remove.append(i+1)
                for i in g_idx:
                    if i in mutations:
                        to_remove.append(i)
                        to_remove.append(i-1)
                for i in to_remove:
                    c_in_cpg[i] = 2
                    g_in_cpg[i] = 2
            cpg_plus = [x for x,y in enumerate(c_in_cpg) if y == 1]
            cpg_minus = [x for x,y in enumerate(g_in_cpg) if y == 1]
            for read in chrom_bam:
                if read.is_unmapped:
                    print(f"Read {read.qname} on chrom {chrom_name} is unampped. Skipping.")
                    continue
                seq = read.query_sequence
                cig = read.cigar
                flag = read.flag
                mapq = read.mapping_quality
                qual = read.query_qualities[:]  # copy, so as not to change in object
                qual = [1 if x>qual_thresh else 0 for x in qual]
                pos = read.reference_start
                strand = "-" if read.is_reverse else "+"
                if mapq < mapq_thresh:
                    low_mapq_reads += 1
                    continue
                if len(seq) != len(qual):
                    corrupt_reads += 1
                    continue
                #if read.is_duplicate:
                #   continue
                if len(cig) != 1:
                    (seq, qual) = self.find_indel(cig, seq, qual)
                if pos + len(seq) > chr_lengths[chrom]:
                    print(f"Chrom {chrom_name} is {chr_lengths[chrom]} bp. Current read ({read.qname}) starts at {pos} and is {len(seq)} bp")
                    continue

                if trim_ends:
                    if strand == "+":
                        #seq = seq[1:-2]
                        #del qual[0]
                        #del qual[-2:]
                        qual[0] = 0
                        qual[-2:] = [0,0]
                    else:
                        #seq = seq[2:-1]
                        #del qual[:2]
                        #del qual[-1]
                        qual[:2] = [0,0]
                        qual[-1] = 0
                read_a = [1 if x=="A" else 0 for x in seq]
                read_a = np.array(read_a)*qual
                read_g = [1 if x=="G" else 0 for x in seq]
                read_g = np.array(read_g)*qual
                read_c = [1 if x=="C" else 0 for x in seq]
                read_c = np.array(read_c)*qual
                read_t = [1 if x=="T" else 0 for x in seq]
                read_t = np.array(read_t)*qual
                if strand == "+":
                    tot_a_plus[pos:pos+len(read_a)] += read_a
                    tot_g_plus[pos:pos+len(read_g)] += read_g
                    tot_c_plus[pos:pos+len(read_c)] += read_c
                    tot_t_plus[pos:pos+len(read_t)] += read_t
                else:
                    tot_a_minus[pos:pos+len(read_t)] += read_t
                    tot_g_minus[pos:pos+len(read_c)] += read_c
                    tot_c_minus[pos:pos+len(read_g)] += read_g
                    tot_t_minus[pos:pos+len(read_a)] += read_a
            chr_a = np.zeros(len(cpg_plus)*2)
            chr_g = np.zeros(len(cpg_plus)*2)
            chr_c = np.zeros(len(cpg_plus)*2)
            chr_t = np.zeros(len(cpg_plus)*2)
            if library == "single":
                chr_a[0::2] = tot_a_plus[cpg_plus]
                chr_a[1::2] = tot_a_minus[cpg_minus]
                chr_g[0::2] = tot_g_plus[cpg_plus]
                chr_g[1::2] = tot_g_minus[cpg_minus]
                chr_c[0::2] = tot_c_plus[cpg_plus]
                chr_c[1::2] = tot_c_minus[cpg_minus]
                chr_t[0::2] = tot_t_plus[cpg_plus]
                chr_t[1::2] = tot_t_minus[cpg_minus]
            elif library == "double":
                chr_a[0::2] = tot_a_plus[cpg_plus]+tot_t_minus[cpg_plus]
                chr_a[1::2] = tot_a_minus[cpg_minus]+tot_t_plus[cpg_minus]
                chr_g[0::2] = tot_g_plus[cpg_plus]+tot_c_minus[cpg_plus]
                chr_g[1::2] = tot_g_minus[cpg_minus]+tot_c_plus[cpg_minus]
                chr_c[0::2] = tot_c_plus[cpg_plus]+tot_g_minus[cpg_plus]
                chr_c[1::2] = tot_c_minus[cpg_minus]+tot_g_plus[cpg_minus]
                chr_t[0::2] = tot_t_plus[cpg_plus]+tot_a_minus[cpg_plus]
                chr_t[1::2] = tot_t_minus[cpg_minus]+tot_a_plus[cpg_minus]
            c_to_t = chr_t/(chr_c+chr_t)
            if library == "single":
                g_to_a = chr_a/(chr_a+chr_g)
            else:
                g_to_a = []
            self.no_a[chrom] = chr_a
            self.no_g[chrom] = chr_g
            self.no_c[chrom] = chr_c
            self.no_t[chrom] = chr_t
            self.c_to_t[chrom] = c_to_t
            self.g_to_a[chrom] = g_to_a
        self.library = library
        self.coord_per_position = [2]*len(tot_chroms)
        self.no_chrs = len(tot_chroms)
        self.chr_names = ["chr"+bam.references[x] for x in chroms]  # lists names of chroms requested in order sent
        #add object name
        bam.close()        

    def parse_infile(self, infile):
        i = 0
        flag = 0 #to keep track of which set of chr lines up to in file
        j = 0
        self.chr_names = []
        with open(infile, "rt") as amfile:
            for line in amfile:
                line = line.rstrip("\n") #remove trailing line feeds
                line = line.lstrip("\t") #remove leading tabs
                fields = line.split(":")
                if len(fields) > 1:
                    fields[1] = fields[1].lstrip() #remove leading space
                if i < 5: #first 5 lines are headers
                    if fields[0] == "Name":
                        self.name = fields[1]
                    elif fields[0] == "Abbreviation":
                        self.abbrev = fields[1]
                    elif fields[0] == "Species":
                        self.species = fields[1]
                    elif fields[0] == "Reference":
                        self.reference = fields[1]
                    elif fields[0] == "Library":
                        self.library = fields[1]
                elif i == 6: #filtered line
                    if "not" in line:
                        self.is_filtered = False
                        flag = 1 #no first set of chr lines for filtered
                    else:
                        self.is_filtered = True
                elif i > 6:
                    if fields[0] == "Method":
                        self.p_filters["method"] = fields[1]
                    elif fields[0] == "max_coverage":
                        cov = fields[1].split()
                        cov = [int(x) if x.isdigit() else np.nan for x in cov] #convert all numbers to int, leave NaN alone
                        self.p_filters["max_coverage"] = cov
                    elif fields[0] == "max_TsPerCoverage":
                        self.p_filters["max_TsPerCoverage"] = []
                    elif "chr" in fields[0] and not flag: #no flag till first set of chr lines done
                        ts = fields[1].split()
                        ts = [int(x) if x.isdigit() else np.nan for x in ts] #convert all numbers to int, leave NaN alone
                        self.p_filters["max_TsPerCoverage"].append(ts)
                        #self.chr_names.append(fields[0])
                    elif fields[0] == "max_aOga" or fields[0] == "max_g_to_a":
                        flag = 1 #after first set of chr lines
                        aga = fields[1].split()
                        aga = [float(x) for x in aga] #convert all numbers to float (NaN is a float)
                        self.p_filters["max_g_to_a"] = aga
                    elif fields[0] == "max_No_As":
                        max_a = fields[1].split()
                        max_a = [int(x) if x.isdigit() else np.nan for x in max_a] #convert all numbers to int, leave NaN alone
                        self.p_filters["max_a"] = max_a
                    elif "simulated" in fields[0]:
                        self.is_simulated = False if "not" in fields[0] else True
                    elif fields[0] == "Coordinates per position":
                        coord = fields[1].split()
                        coord = [int(x) if x.isdigit() else np.nan for x in coord] #convert all numbers to int, leave NaN alone
                        self.coord_per_position = coord
                    elif fields[0] == "Deamination rate":
                        self.d_rate["rate"] = {"local":[], "dlocal":[], "no_positions":[]}
                    elif fields[0] == "Reference":
                        self.d_rate["ref"] = fields[1]
                    elif fields[0] == "min_beta":
                        self.d_rate["min_beta"] = float(fields[1])
                    elif fields[0] == "min_coverage":
                        self.d_rate["min_coverage"] = float(fields[1])
                    elif fields[0] == "global":
                        self.d_rate["rate"]["global"] = float(fields[1])
                    elif fields[0] == "dglobal":
                        self.d_rate["rate"]["dglobal"] = float(fields[1])
                    elif fields[0] == "local":
                        local = fields[1].split()
                        local = [float(x) for x in local] #convert all numbers to float (NaN is a float)
                        self.d_rate["rate"]["local"] = local
                    elif fields[0] == "dlocal":
                        dlocal = fields[1].split()
                        dlocal = [float(x) for x in dlocal] #convert all numbers to float (NaN is a float)
                        self.d_rate["rate"]["dlocal"] = dlocal
                    elif fields[0] == "no_positions":
                        no_positions = fields[1].split()
                        no_positions = [float(x) for x in no_positions] #convert all numbers to float (NaN is a float)
                        self.d_rate["rate"]["no_positions"] = no_positions
                    elif fields[0] == "Effective coverage":
                        eff_cov = fields[1].split()
                        eff_cov = [float(x) for x in eff_cov] #convert all numbers to float (NaN is float)
                        self.diagnostics["effective_coverage"] = eff_cov
                    elif "chr" in fields[0] and flag == 1: #chr lines in all files
                        self.chr_names.append(fields[0])
                    elif fields[0] == "No_As":
                        no_a = fields[1].split()
                        no_a = [int(x) if x.isdigit() else np.nan for x in no_a] #convert all numbers to int, leave NaN alone
                        self.no_a.append(no_a)
                    elif fields[0] == "No_Cs":
                        no_c = fields[1].split()
                        no_c = [int(x) if x.isdigit() else np.nan for x in no_c] #convert all numbers to int, leave NaN alone
                        self.no_c.append(no_c)
                    elif fields[0] == "No_Gs":
                        no_g = fields[1].split()
                        no_g = [int(x) if x.isdigit() else np.nan for x in no_g] #convert all numbers to int, leave NaN alone
                        self.no_g.append(no_g)
                    elif fields[0] == "No_Ts":
                        no_t = fields[1].split()
                        no_t = [int(x) if x.isdigit() else np.nan for x in no_t] #convert all numbers to int, leave NaN alone
                        self.no_t.append(no_t)
                    elif fields[0] == "aOga" or fields[0] == "g_to_a":
                        g_to_a = fields[1].split()
                        g_to_a = [int(x) if x.isdigit() else np.nan for x in g_to_a] #convert all numbers to int, leave NaN alone
                        self.g_to_a.append(g_to_a)
                    elif fields[0] == "tOct" or fields[0] == "c_to_t":
                        c_to_t = fields[1].split()
                        c_to_t = [int(x) if x.isdigit() else np.nan for x in c_to_t] #convert all numbers to int, leave NaN alone
                        self.c_to_t.append(c_to_t)
                    elif fields[0] == "Reconstructed methylation":
                        self.methylation["methylation"] = []
                        flag = 2 #done with second set of chr lines
                    elif fields[0] == "win_size":
                        win = fields[1].split()
                        win = [int(x) if x.isdigit() else np.nan for x in win] #convert all numbers to int, leave NaN alone
                        self.methylation["win_size"] = win
                    elif fields[0] == "slope":
                        slope = fields[1].split()
                        slope = [float(x) for x in slope] #convert all numbers to float (NaN is float)
                        self.methylation["slope"] = slope
                    elif fields[0] == "intercept":
                        intercept = fields[1].split()
                        intercept = [float(x) for x in intercept] #convert all numbers to float (NaN is float)
                        self.methylation["intercept"] = intercept
                    elif fields[0] == "lcf":
                        lcf = fields[1].split()
                        lcf = [float(x) for x in lcf] #convert all numbers to float (NaN is float)
                        self.methylation["lcf"] = lcf
                    elif "chr" in fields[0] and flag == 2: #third set of chr lines
                        meth = fields[1].split()
                        meth = [float(x) for x in meth] #convert all numbers to float (NaN is float)
                        self.methylation["methylation"].append(meth)
                i += 1
        self.no_chrs = len(self.coord_per_position)

    def get_base_no(self, chrom, base):
        if isinstance(chrom, str):
            chr_ind = self.indexofchr([chrom])[0]
        else:
            chr_ind = chrom
        if base == "a":
            result = self.no_a
        elif base == "c":
            result = self.no_c
        elif base == "g":
            result = self.no_g
        elif base == "t":
            result = self.no_t
        elif base == "ga":
            result = self.g_to_a
        return result[chr_ind]

    def smooth(self, chrom, winsize):
        no_t = self.get_base_no(chrom, "t")
        no_c = self.get_base_no(chrom, "c")
        no_t = np.array(no_t) #convert to numpy array in order to add elementwise
        no_c = np.array(no_c)
        no_ct = no_t + no_c
        tpl = np.ones(winsize)
        no_t = t.nanconv(no_t, tpl, "same")
        no_ct = t.nanconv(no_ct, tpl, "same")
        return(no_t, no_ct)

    def region_methylation(self, region, gc):
        region = t.standardize_region(region) #take region input and change to dictionary
        chrom = region["chrom"] #get chrom from region dict
        chr_ind = gc.indexofchr([chrom])[0] #find index of region chrom in gc object

        cpg_start = np.where(gc.coords[chr_ind] >= region["start"])[0][0] #get index of first
        cpg_end = np.where(gc.coords[chr_ind] <= region["end"])[0][-1]    #and last coords in region
        chr_ind = self.indexofchr([chrom])[0] #find index of chrom in ams object
        no_t = np.nansum(self.no_t[chr_ind][cpg_start:cpg_end+1])
        no_ct = no_t + np.nansum(self.no_c[chr_ind][cpg_start:cpg_end+1])
        meth = self.methylation["slope"][chr_ind] * no_t / no_ct + self.methylation["intercept"][chr_ind]
        if not np.isnan(meth):
            meth = min(max(meth,0),1)
        return meth

    def diagnose(self, fname=None, span=5, strict=True, tolerance=1e-3, compare=True, max_c_to_t=0.25, max_g_to_a=0.25, max_coverage=100, low_coverage=1):
        if fname is None:
            fname = "data/logs/"+self.name+"_diagnostics.txt"
        
        #initialize
        no_chrs = self.no_chrs
        eff_coverage = np.zeros(no_chrs)
        coverage_threshold = np.zeros(no_chrs)
        c_to_t_threshold = [None] * no_chrs #set number of elements to add later by chrom index (use append instead?)
        fid = open(fname, "w")

        #compute number of rows and columns in output figures
        no_rows = math.floor(math.sqrt(no_chrs))
        if no_rows:
            no_cols = math.ceil(no_chrs/no_rows)
        #fig_coverage = plt.figure("coverage") #if we decide to plot figures, find faster utility
        #fig_coverage.suptitle(self.name)
        #fig_perc_removed_per_coverage = plt.figure("removed_per_coverage");
        #fig_perc_removed_per_coverage.suptitle(self.name)
        #fig_perc_removed = plt.figure("removed")
        #fig_perc_removed.suptitle(self.name)
        
        #Loop on chromosomes
        for chrom in range(no_chrs):
            #report
            if self.chr_names:
                print(f"Diagnosing {self.chr_names[chrom]}")
                fid.write(f"Diagnosing {self.chr_names[chrom]}\n")
            else:
                print(f"Diagnosing chrom #{chrom+1}")
                fid.write(f"Diagnosing chrom #{chrom+1}\n")
            
            #vectors of current chromosome
            no_t = self.no_t[chrom]
            no_t = np.array(no_t, dtype="float") #convert to numpy array in order to add elementwise
            no_c = self.no_c[chrom]
            no_c = np.array(no_c)
            no_ct = no_t + no_c
            if self.library == "SS": #test this function on single stranded data
                no_a = self.no_a[chrom]
                g_to_a = self.g_to_a[chrom]
                if not g_to_a:
                    no_g = self.no_g[chrom]
                    no_g = np.array(no_g, dtype="float") #convert to numpy array in order to add elementwise
                    no_a = np.array(no_a)
                    if no_g:
                        g_to_a = no_a / (no_a + no_g)
            fid.write(f"\tNumber of CpG positions: {len(no_ct):,d}\n")
            fid.write(f"\tInitial coverage: {np.nanmean(no_ct):.1f}\n")
            
            #compute the upper threshold of {nt} in a few steps
            #crude outlier removal using very loose criterion
            percentiles = t.nan_prctile(no_ct, [25,50,75])
            fid.write(f"\t[C+T]: median = {percentiles[1]:.0f}\n")
            fid.write(f"\t     [25th 75th] percentiles = {percentiles[[0,2]]}")
            iqrange = percentiles[2] - percentiles[0]
            thresh = percentiles[2] + span*iqrange
            fid.write(f"=> initial threshold was set to {thresh:.0f}\n")
            to_remove = np.where(no_ct > thresh)[0]
            fid.write(f"\t     {len(to_remove):,d} positions ({100*len(to_remove)/len(no_ct):,.2f}%) removed as no_ct > {thresh:.0f}\n")
            for x in to_remove:
                no_ct[x] = np.nan
            
            #evaluate the parameters of the normal distribution of {nt}
            nt_zeros = np.where(no_ct == 0)[0]
            nz_ct = np.copy(no_ct) #copy array
            for x in nt_zeros:
                nz_ct[x] = np.nan
            #N = plt.hist(nz_ct,bins=range(int(thresh)))
            more_to_remove = []
            if strict:
                print("in strict loop") #handle later
            coverage_threshold[chrom] = thresh
            
            #plot the coverage
            #plt.figure("coverage")
            #plt.subplot(no_rows,no_cols,chrom+1) #chrom num, not index num
            #plt.bar(len(N[0]), N[0])
            #plt.plot([thresh, thresh], list(plt.gca().get_ylim()), "k")
            if self.chr_names:
                xlabel = self.chr_names[chrom]
            else:
                xlabel = f"Chromosome #{chrom+1}"
            #plt.xlabel(xlabel)
            #plt.show()

            #remove outliers
            to_remove = [to_remove, more_to_remove] #more_to_remove gets vals in skipped strict loop
            for x in to_remove:
                for y in x:
                    no_ct[y] = np.nan
            for x in to_remove:
                for y in x:
                    no_t[y] = np.nan
            
            #effective coverage
            eff_coverage[chrom] = np.nanmean(no_ct)
            fid.write(f"\t     Effective coverage: {eff_coverage[chrom]:.1f}\n")
            max_c = float(max_coverage)
            
            #compare to previous PCR-duplicate threshold
            prev_to_remove = np.where(no_ct>max_c)[0]
            fid.write(f"\t     COMPARISON: {len(prev_to_remove):,d} positions ({100*len(prev_to_remove)/len(no_ct):.2f}%) ")
            fid.write(f"would have been removed if the PCR-duplicate threshold were {max_coverage}\n")
            
            #report on the zero bin
            fid.write(f"\t{len(nt_zeros):,d} positions ({100*len(nt_zeros)/len(no_ct):.2f}%) have no_ct = 0\n")

            #analyze each coverage level independently
            th_c2g = np.zeros(int(thresh))
            chr_cov = [x for x in range(low_coverage, int(thresh+2))] 
            no_pos_tot = np.zeros(int(thresh+1))
            no_pos_removed = np.zeros(int(thresh+1))
            for cover in range(int(thresh-1),low_coverage-2, -1): 
                fid.write(f"\tAnalyzing [xt] coming from coverage {cover+1}:\n")
                idx = np.where(no_ct==cover+1)[0]
                no_pos_cover = len(idx)
                no_pos_tot[cover] = no_pos_cover 
                fid.write(f"\t\tTotal number of positions: {no_pos_cover:,d}\n")
                H = np.histogram(no_t[idx],bins=range(cover+3)) #+1 to compensate for 0 base, +1 for slice
                #ie in matlab, 0:5 = [0,1,2,3,4,5], but in python, array of range(5) = [0,1,2,3,4]
                p, w, llike = t.bmm(H, [0.01, 0.5, 0.99], [0.9, 0.05, 0.05], tolerance, [0, 1, 0])
                fid.write("\t\tEstimated homozygous mutations: ")
                fid.write(f"Pr(meth) = {p[2]:.2f}; Weight = {w[2]:.3f}\n")
                fid.write("\t\tEstimated heterozygous mutations: ")
                fid.write(f"Pr(meth) = {p[1]:.2f}; Weight = {w[1]:.3f}\n")
                fid.write("\t\tEstimated deaminated positions: ")
                fid.write(f"Pr(meth) = {p[0]:.2f}; Weight = {w[0]:.3f}\n")
                thresh = np.ceil((np.log(w[1]/w[0]) + (cover+1) * np.log(1/2/ (1-p[0]))) / np.log(p[0]/(1-p[0]))) -1
                th_c2g[cover] = thresh
                more_to_remove = idx[no_t[idx]>thresh] 
                no_pos_removed[cover] = len(more_to_remove)
                fid.write(f"\t\t{len(more_to_remove):,d} positions ({100*len(more_to_remove)/len(no_ct):.2f}%) removed ")
                fid.write(f"as xt > {thresh:.0f} ")
                fid.write(f"(corresponding to C2T-ratio of {(thresh+1)/(cover+1):.3f})\n")
                fid.write(f"\t\tThe threshold ({thresh:.0f}) is expected to remove ")
                fid.write(f"{w[0]*no_pos_cover*stats.binom.sf(thresh,cover+1,p[0]):.1f} true positives\n")
                fid.write(f"\t\tand to include {w[1]*no_pos_cover*stats.binom.cdf(thresh-1,cover+1,0.5)+w[2]*no_pos_cover*stats.binom.cdf(thresh-1,cover+1,p[2]):.1f} false positives\n")
                if thresh == cover+1:
                    #what happens with no threshold?
                    fid.write("\t\tIf no threshold is used (no positions removed), ")
                    fid.write(f"then {(w[1]+w[2])*no_pos_cover:.1f} false positives will be retained\n")
                #compute positions removed using G->A and C->T
                if compare:
                    if np.isfinite(max_c_to_t):
                        c2t_crit = [x for x in idx if no_t[x]>1 and no_t[x]/no_ct[x] >= max_c_to_t]
                        fid.write("\t\tCOMPARISON: Using C->T ")
                        fid.write(f"{len(c2t_crit):,d} positions are removed.\n")
                    if self.library == "SS": #test sections with this condition
                        g2a_crit = [x for x in idx if (no_a[x]==1 and a_to_g[x]>=max_a_to_g) or no_a[x]>1]
                        fid.write("\t\tCOMPARISON: Using G->A ")
                        fid.write(f"{len(g2a_crit):,d} are removed.")
                        if np.isfinite(max_c_to_t):
                            g2a_additional = np.setdiff1d(g2a_crit, c2t_crit)
                            fid.write(f" Of them, {len(g2a_additional):,d} positions were not removed ")
                            fid.write("by the C->T ration threshold.\n")
                            flat_remove = [x for y in more_to_remove for x in y] #flatten 2d list to do intersect
                            crits = c2t_crit+g2a_crit
                            in_both = list(set(crits).intersection(set(flat_remove)))
                            fid.write("\t\tCOMPARISON: ")
                            if len(in_both) == len(more_to_remove):
                                fid.write("All positions removed only by looking at [xt] ")
                                fid.write(f"({len(more_to_remove):,d}) also removed when looking at ")
                                fid.write("both C->T and G->A.\n")
                            else:
                                fid.write(f"Of the {len(more_to_remove):,d} positions removed only by looking at ")
                                fid.write(f"[xt], {len(in_both):,d} ({100*len(in_both)/len(more_to_remove):.1f}%)")
                                fid.write(" are also removed when looking at both C->T and G->A.\n")
                        else:
                            fid.write("\n")
                        in_crit = list(set(g2a_crit).intersection(set(flat_remove)))
                        fid.write("\t\tCOMPARISON: ")
                        if len(in_crit) == len(more_to_remove):
                            fid.write("All positions removed only by looking at [xt] ")
                            fid.write(f"({len(more_to_remove):,d}) also removed when looking at G->A.\n")
                        else:
                            fid.write(f"Of the {len(more_to_remove):,d} positions removed only by looking at ")
                            fid.write(f"[xt], {len(in_crit):,d} ({100*len(in_crit)/len(more_to_remove):.1f}%)")
                            fid.write(" are also removed when looking at G->A.\n")
                #remove the positions
                for x in more_to_remove:
                    no_ct[x] = np.nan #assumes 1d array
                for x in more_to_remove:
                    no_t[x] = np.nan
            c_to_t_threshold[chrom] = th_c2g

            #plot ratio of removed to total per coverage
            quotient = no_pos_removed/no_pos_tot*100
            #plt.figure("removed_per_coverage")
            #plt.subplot(no_rows,no_cols,chrom+1) #chrom num, not index num
            #plt.plot(chr_cov,quotient,"bp")
            if self.chr_names:
                xlabel = self.chr_names[chrom]
            else:
                xlabel = f"Chromosome #{chrom+1}"
            #plt.xlabel(xlabel)
            ylabel = "%removed(coverage)/total(coverage)"
            #plt.ylabel(ylabel)
            #plt.grid(True)
            #plt.show()

            #plot ratio of removed to total
            tot = sum(no_pos_tot)
            #plt.figure("removed")
            #plt.subplot(no_rows,no_cols,chrom+1)
            #plt.plot(chr_cov,no_pos_removed/tot*100,"bp")
            if self.chr_names:
                xlabel = self.chr_names[chrom]
            else:
                xlabel = f"Chromosome #{chrom+1}"
            #plt.xlabel(xlabel)
            ylabel = "%removed(coverage)/total(coverage)"
            #plt.ylabel(ylabel)
            #plt.grid(True)
        #plt.show()

        fid.close()
        self.diagnostics["effective_coverage"] = eff_coverage
        self.p_filters["method"] = ""
        self.p_filters["max_coverage"] = coverage_threshold
        self.p_filters["max_TsPerCoverage"] = c_to_t_threshold
        self.p_filters["max_g_to_a"] = []
        self.p_filters["max_a"] = []

    @staticmethod
    def extend_removed(to_remove):
        odds = np.where(to_remove%2)
        evens = np.where(~to_remove%2)
        to_remove = list(to_remove) + list(to_remove[odds]-1) + list(to_remove[evens]+1)
        to_remove = np.unique(to_remove)
        return to_remove

    def filter(self, max_c_to_t = None, min_t = 1, merge = True, fname = None, max_g_to_a = .25, max_a = 1, method = None, max_coverage = None, max_TsPerCoverage = None):
        #initialize
        no_chr = self.no_chrs
        if fname == None:
            fname = "data/logs/"+self.name+"_filter.txt"
        is_c_to_t = False
        is_ts_per_cov = False
        if max_coverage == None:
            max_coverage = self.p_filters["max_coverage"] if self.p_filters["max_coverage"].any() else 100
        if method == None:
            if self.library == "SS": #test this function on single stranded data
                method = "both"
            else:
                method = "c_to_t"
        if max_TsPerCoverage is not None:
            is_ts_per_cov = True
        else:
            if max_c_to_t is not None:
                max_TsPerCoverage = max_c_to_t
            else:
                max_TsPerCoverage = self.p_filters["max_TsPerCoverage"] if self.p_filters["max_TsPerCoverage"] else .25
        if max_c_to_t is not None:
            is_c_to_t = True
        #'max_tOct' and 'max_TsPerCoverage' cannot be used at the same time
        if is_c_to_t and is_ts_per_cov:
            print("Both 'max_TsPerCoverage' and 'max_c_to_t' were used in the input")
        #bring input parameters into standard form - max_coverage
        if np.isscalar(max_coverage):
            max_coverage = max_coverage * np.ones(no_chr)
        #bring input parameters into standard form - max_TsPerCoverage
        if np.isscalar(max_TsPerCoverage):
            tmp_max = max_TsPerCoverage
        else:
            tmp_max = max_TsPerCoverage.copy() 
        if is_c_to_t: #max_c_to_t
            if np.isscalar(max_TsPerCoverage): #single ratio
                for chrom in range(no_chr):
                    max_TsPerCoverage[chrom] = math.ceil(tmp_max * range(max_coverage[chrom])-1)
            else: #ratio per chrom or per coverage per chrom
                for chrom in range(no_chr):
                    max_TsPerCoverage[chrom] = math.ceil(tmp_max[chrom] * range(max_coverage[chrom])-1)
        else: #max_TsPerCoverage or default
            if np.isscalar(max_TsPerCoverage): #single number
                for chrom in range(no_chr):
                    max_TsPerCoverage[chrom] = tmp_max * np.ones(int(max_coverage[chrom]))
            else: #number per chrom
                for chrom in range(no_chr):
                    max_TsPerCoverage[chrom] = tmp_max[chrom] * np.ones(int(max_coverage[chrom]))
        #bring input parameters into standard form - max_g_to_a
        if np.isscalar(max_g_to_a):
            max_g_to_a = max_g_to_a * np.ones(no_chr)
        #bring input parameters into standard form - max_a
        if np.isscalar(max_a):
            max_a = max_a * np.ones(no_chr)
        #bring input parameters into standard form - min_t
        for chrom in range(no_chr):
            vec = max_TsPerCoverage[chrom] 
            vec = [min_t if x < min_t else x for x in vec]
            vec[0:min_t-1] = range(1,min_t) #min_t for a given coverage can't be more than the coverage
            max_TsPerCoverage[chrom] = vec
        #f_ams = copy.deepcopy(self) #copy object. copy incredible slow. nec? ie does it matter if ams is changed?
        f_ams = self #don't need to worry abt overwriting, since output has diff filename
        f_ams.p_filters["method"] = method
        f_ams.p_filters["max_coverage"] = max_coverage
        f_ams.p_filters["max_TsPerCoverage"] = max_TsPerCoverage
        f_ams.p_filters["max_g_to_a"] = max_g_to_a
        f_ams.p_filters["max_a"] = max_a
        fid = open(fname, "w")
        #loop on chromosomes
        tot_removed = 0
        for chrom in range(no_chr):
            #report
            if self.chr_names:
                print(f"Filtering {self.chr_names[chrom]}")
                fid.write(f"Filtering {self.chr_names[chrom]}\n")
            else:
                print(f"Filtering chrom #{chrom+1}")
                fid.write(f"Filtering chrom #{chrom+1}\n")
            #vectors of current chromosome
            no_t = self.no_t[chrom]
            no_t = np.array(no_t, dtype="float") #convert to numpy array in order to add elementwise
            no_c = self.no_c[chrom]
            no_c = np.array(no_c)
            no_ct = no_t + no_c
            if method != "c_to_t":
                no_a = self.no_a[chrom]
                g_to_a = self.g_to_a[chrom]
            no_pos = len(no_t)
            print("\tNumber of CpG positions ", end="")
            print(f"({self.coord_per_position[chrom]} coordinates per position): {no_pos:,d}")
            fid.write("\tNumber of CpG positions ")
            fid.write(f"({self.coord_per_position[chrom]} coordinates per position): {no_pos:,d}\n")
            #remove positions whose coverage is too high
            to_remove = np.where(no_ct>max_coverage[chrom])[0]
            no_removed = len(to_remove)
            fid.write(f"\t{no_removed:,d} positions ({100*no_removed/no_pos:.2f}%) removed as No_CTs > {int(max_coverage[chrom])}\n")
            #loop on each coverage level and remove positions with high No_Ts
            if method != "g_to_a":
                #loop over all coverage levels
                for cover in range(int(max_coverage[chrom]),0,-1):
                    idx = np.where(no_ct==cover)[0]
                    fid.write(f"\tcoverage {cover} ({len(idx):,d} positions):")
                    more_to_remove = idx[no_t[idx]>max_TsPerCoverage[chrom][cover-1]] 
                    more_to_remove = self.extend_removed(more_to_remove)
                    no_removed = len(to_remove) #redefined for every coverage level
                    to_remove = np.unique(list(to_remove)+list(more_to_remove))
                    no_removed = len(to_remove) - no_removed
                    fid.write(f"\t{no_removed:,d} (extended) positions ")
                    fid.write(f"were removed as No_Ts > {int(max_TsPerCoverage[chrom][cover-1])}\n")
            #remove more positions if data on A's and G's is available
            if method != "c_to_t":
                more_to_remove = np.where((no_a<=max_no_a[chrom] and g_to_a>=max_g_to_a[chrom]) or no_a>max_a[chrom])[0]
                more_to_remove = self.extend_removed(more_to_remove)
                no_removed = len(to_remove)
                to_remove = np.unique(list(to_remove)+list(more_to_remove))
                no_removed = len(to_remove) - no_removed
                fid.write(f"\t{no_removed:,d} additional positions ")
                fid.write(f"removed as (1) No_As > {int(max_a[chrom])} or (2) No_As <= {int(max_a[chrom])} ")
                fid.write(f"and aOga >= {max_g_to_a[chrom]:.2f}\n")
            #remove positions and keep only relevant vectors in p_CpGs
            no_removed = len(to_remove)
            print(f"Overall {no_removed:,d} positions ({100*no_removed/no_pos:.2f}%) were removed")
            fid.write(f"Overall {no_removed:,d} positions ({100*no_removed/no_pos:.2f}%) were removed\n")
            for x in to_remove:
                no_t[x] = np.nan
                no_ct[x] = np.nan
            np.array(no_ct) #for elementwise subtraction
            np.array(no_t)
            no_c = no_ct - no_t #does this behave? (should be elementwise subtraction)
            #merge positions
            if merge:
                if f_ams.coord_per_position[chrom] == 1: #test on single stranded
                    if f_ams.chr_names:
                        print(f"{f_ams.chr_names[chrom]} had already gone through merger")
                    else:
                        print(f"Chromosome #{chrom} had already gone through merger")
                else:
                    f_ams.coord_per_position[chrom] = 1
                    no_t = t.nanmerge(no_t, "sum")
                    no_c = t.nanmerge(no_c, "sum")
            #substitute into f_ams
            f_ams.no_t[chrom] = no_t
            f_ams.no_c[chrom] = no_c
            f_ams.diagnostics["effective_coverage"][chrom] = (np.nansum(no_t) + np.nansum(no_c))/len(no_t)
            if f_ams.no_a: f_ams.no_a[chrom] = []
            if f_ams.no_g: f_ams.no_g[chrom] = []
            if f_ams.g_to_a: f_ams.g_to_a[chrom] = []
            if f_ams.c_to_t: f_ams.c_to_t[chrom] = []
            tot_removed = tot_removed + no_removed
        f_ams.is_filtered = True
        print(f"In total {tot_removed:,d} positions were removed")
        #close file
        fid.close()

    def estimate_drate(self, method="reference", global_meth=np.nan, min_cov=1, ref=[], min_beta=1):
        ref_params = {}
        global_params = {}
        if global_meth > 1: #can be dec or %
            global_meth = 0.01*global_meth
        global_params["global_methylation"] = global_meth
        ref_params["ref"] = ref
        if min_beta > 1:
            min_beta = 0.01*min_beta
        ref_params["min_beta"] = min_beta
        ref_params["min_coverage"] = min_cov
        global_params["min_coverage"] = min_cov
        if method == "global":
            meth_params = global_params
            if np.isnan(meth_params["global_methylation"]):
                raise Exception("No global_methylation param provided when using 'global' method")
        elif method == "reference":
            meth_params = ref_params
            if not meth_params["ref"]:
                raise Exception("No reference provided when using 'reference' method")
        elif method != "estref":
            raise Exception(f"Unknown method '{method}'")
        print(f"Estimating deamination rate using the '{method}' method")
        if method == "global":
            print(f"\tglobal_methylation: {meth_params['global_methylation']:.2f}")
        elif method == "reference":
            print(f"\tmin_beta: {meth_params['min_beta']:.2f}")
        if method == "global" or method == "reference":
            print(f"\tmin_coverage: {meth_params['min_coverage']:d}")
        
        #sanity check
        if not self.is_filtered:
            raise Exception(f"{self.name} is not filtered")
        if not all(x == 1 for x in self.coord_per_position):
            raise Exception(f"{self.name} is not merged")
        
        #verify reference is merged and scaled
        if method == "reference":
            meth_params["ref"].merge() #ref should be type mms, so this should work
            meth_params["ref"].scale()
        
        #initialize
        drate = {"global":np.nan, "dglobal":np.nan, "local":np.full((self.no_chrs),np.nan), "dlocal":np.full((self.no_chrs),np.nan), "no_positions":np.full((self.no_chrs),np.nan)}
        factor = 1
        if method == "global":
            factor = 1/meth_params["global_methylation"]
        tot_t = 0
        tot_ct = 0
        
        #loop on chromosomes
        for chrom in range(self.no_chrs):
            
            #report
            print(f"Estimating deamination rate in {self.chr_names[chrom]}")
            
            #vectors of current chrom
            no_t = self.no_t[chrom]
            no_t = np.array(no_t) #convert to numpy array in order to add elementwise
            no_c = self.no_c[chrom]
            no_c = np.array(no_c)
            no_ct = no_t + no_c

            #remove positions for which the reference has low beta-values
            if method == "reference":
                #Take only positions where {ref}>=min_beta
                include = np.where(np.array(meth_params["ref"].get_methylation(self.chr_names[chrom])[1])>=meth_params["min_beta"])[0]
                no_t = no_t[include]
                no_ct = no_ct[include]

            #remove positions that are not covered high enough
            include = np.where(no_ct >= meth_params["min_coverage"])[0]
            no_t = no_t[include]
            no_ct = no_ct[include]

            #compute estimates based on the current chromosome
            drate["local"][chrom] = factor * np.nansum(no_t)/np.nansum(no_ct)
            drate["no_positions"][chrom] = sum(np.isfinite(no_t))
            drate["dlocal"][chrom] = factor * math.sqrt(drate["local"][chrom]*(1-drate["local"][chrom])/drate["no_positions"][chrom])

            #accumulate sums
            tot_t = tot_t + np.nansum(no_t)
            tot_ct = tot_ct + np.nansum(no_ct)
        
        #evaluate global degradation rate
        drate["global"] = factor * tot_t/tot_ct
        drate["dglobal"] = factor * math.sqrt(drate["global"] * (1-drate["global"])/sum(drate["no_positions"]))

        #plug vals into object
        if method == "reference":
            self.d_rate = {"method":"reference", "rate":drate, "ref":meth_params["ref"].name, "min_beta":meth_params["min_beta"], "min_coverage":meth_params["min_coverage"]}
        elif method == "global":
            self.d_rate = {"method":"global", "rate":drate, "global_methylation":meth_params["global_methylation"], "min_coverage":meth_params["min_coverage"]}
        elif method == "estref":
            self.d_rate = {"method":"estref", "rate":drate}
    
    def determine_winsize(self, chrom, method="prob", coverage=None, drate=None, min_meth=0.2, p0=0.01, k_inv=1/2.5, max_width=31):
        if coverage == None:
            coverage = self.diagnostics["effective_coverage"][chrom]
        if drate == None:
            drate=self.d_rate["rate"]["global"]
        if not method=="prob" and not method=="relerror":
            raise Exception(f"Unknown method '{method}'")
        if min_meth>1:
            min_meth = 0.01 * min_meth
        if method == "prob":
            param = p0
        else:
            param = k_inv

        #compute the window size
        p = min_meth * drate
        if method == "prob":
            win_size = np.ceil(np.log(param)/np.log(1-p)/coverage)
        else:
            win_size = ceil((1-p)/coverage/p/param**2)

        #narrow window if too wide
        win_size = min(win_size, max_width)

        #make sure win_size is odd
        if not win_size%2:
            win_size += 1
        
        return win_size
    
    def reconstruct_methylation(self, win_size="auto", winsize_alg={}, function="log", slope=None, intercept=[0], lcf=0.05, report=True):
        no_chr = self.no_chrs
        if slope == None:
            slope=[1/self.d_rate["rate"]["global"]]
        if isinstance(win_size, str):
            auto_win = True
        else:
            auto_win = False
        
        #bring parameters into standard format - win_size
        if not auto_win:
            if len(win_size) == 1:
                if not win_size%2:
                    win_size += 1 #win_size should be odd
                win_size = win_size * np.ones(no_chr)
            else:
                for chrom in range(no_chr):
                    if not win_size[chrom]%2:
                        win_size[chrom] += 1 #win_size should be odd
        else:
            win_size = np.zeros(no_chr)
            for chrom in range(no_chr):
                win_size[chrom] = self.determine_winsize(chrom, **winsize_alg)
        win_size = [int(x) for x in win_size]

        #bring parameters into standard formt - slope
        if len(slope) == 1:
            slope = slope * np.ones(no_chr)

        #bring parameters into standard formt - intercept
        if len(intercept) == 1:
            intercept = intercept * np.ones(no_chr)
        intercept = [int(x) for x in intercept]

        meth = []
        for chrom in range(no_chr):
            if self.chr_names:
                print(f"Computing methylation in {self.chr_names[chrom]}")
            else:
                print(f"Computing methylation in chrom #{chrom+1}")
            #get smoothed No_Ts and No_CTs
            (no_t, no_ct) = self.smooth(chrom, int(win_size[chrom]))
            #remove regions with particularly low coverage
            lct = t.find_low_coverage_thresh(no_ct, lcf)
            no_ct = [np.nan if x < lct else x for x in no_ct]
            #compute methylation
            if function == "log":
                tmp_exp = np.exp(2*slope[chrom]*no_t/no_ct)
                methi = 2*tmp_exp/(tmp_exp+1)-1
            else:
                methi = slope[chrom]*no_t/no_ct+intercept[chrom]
                methi = min(max(methi,0),1) #keep between 0 and 1 (nans untouched)
            if report:
                print(f"Average methylation: {np.nanmean(methi):.2f}")
            meth.append(methi)
        self.methylation = {"methylation":meth, "win_size":win_size, "slope": slope, "intercept":intercept, "lcf":lcf}

    def simulate(self, degrad_rate, mms, report=True):
        
        #initialize
        self.name = self.name + "_(sim)"
        self.is_simulated = True

        #make sure methylation in {mms} is scaled to [0,1]
        mms.scale()

        #parameters needed for the simulation
        meth_map = mms.get_methylation()[1]

        for chrom in range(self.no_chrs):
            if report:
                print(f"Processing {self.chr_names[chrom]} ... ")
            no_c = np.array(self.no_c[chrom])
            tot_ct = no_c + self.no_t[chrom]
            self.no_t[chrom] = np.random.binomial(tot_ct, degrad_rate*np.array(meth_map[chrom]))
            self.no_t[chrom] = [max(0,x) for x in self.no_t[chrom]] #change bad (neg) vals resulting from nans to 0
            self.no_c[chrom] = tot_ct - self.no_t[chrom]
            if report:
                print("done")

    def dump(self, stage):  # currently, dump works only for sorted list of sequential chroms
        aname = self.name
        fname = "data/python_dumps/" + aname + "_" + stage + ".txt"
        with open(fname, "w") as fid:
            fid.write(f"Name: {aname}\nAbbreviation: {self.abbrev}\nSpecies: {self.species}\nReference: {self.reference}\n")
            fid.write(f"Library: {self.library}\n\n")
            filt = "" if self.is_filtered else "not "
            fid.write(f"This sample is {filt}filtered\n")
            for key in self.p_filters.keys():
                if key == "method":
                    fid.write(f"\tMethod: {self.p_filters['method']}\n")
                elif key == "max_coverage":
                    max_cov = [int(x) if ~np.isnan(x) else "NaN" for x in self.p_filters['max_coverage']]
                    fid.write(f"\tmax_coverage: {' '.join(map(str, max_cov))}\n")
                elif key == "max_TsPerCoverage":
                    fid.write("\tmax_TsPerCoverage:\n")
                    for chrom in range(self.no_chrs):
                        max_t = [int(x) if ~np.isnan(x) else "NaN" for x in self.p_filters['max_TsPerCoverage'][chrom]]
                        fid.write(f"\t\t{self.chr_names[chrom]}: {' '.join(map(str, max_t))}\n")
                elif key == "max_g_to_a":
                    fid.write(f"\tmax_g_to_a: {' '.join(map(str, self.p_filters['max_g_to_a']))}\n")
                elif key == "max_a":
                    max_a = [int(x) if ~np.isnan(x) else "NaN" for x in self.p_filters['max_a']]
                    fid.write(f"\tmax_No_As: {' '.join(map(str, max_a))}\n")
            sim = "" if self.is_simulated else "not "
            fid.write(f"This sample is {sim}simulated\n")
            fid.write(f"Coordinates per position: {' '.join(map(str, self.coord_per_position))}\n")
            fid.write("Deamination rate:\n")
            for key in self.d_rate.keys():
                if key == "ref":
                    fid.write(f"\tReference: {self.d_rate['ref']}\n")
                elif key == "min_beta":
                    fid.write(f"\tmin_beta: {int(self.d_rate['min_beta'])}\n")
                elif key == "min_coverage":
                    fid.write(f"\tmin_coverage: {int(self.d_rate['min_coverage'])}\n")
                elif key == "rate":
                    fid.write(f"\tglobal: {self.d_rate['rate']['global']}\n\tdglobal: {self.d_rate['rate']['dglobal']}\n")
                    fid.write(f"\tlocal: {' '.join(map(str, self.d_rate['rate']['local']))}\n")
                    fid.write(f"\tdlocal: {' '.join(map(str, self.d_rate['rate']['dlocal']))}\n")
                    fid.write(f"\tno_positions: {' '.join(map(str, self.d_rate['rate']['no_positions']))}\n")
            for key in self.diagnostics.keys():
                if key == "effective_coverage":
                    fid.write(f"Effective coverage: {' '.join(map(str, self.diagnostics['effective_coverage']))}\n")
            for chrom in range(self.no_chrs):
                fid.write(f"\n{self.chr_names[chrom]}:\n")
                no_a = [int(x) if ~np.isnan(x) else "NaN" for x in self.no_a[chrom]]
                fid.write(f"No_As: {' '.join(map(str, no_a))}\n")
                no_c = [int(x) if ~np.isnan(x) else "NaN" for x in self.no_c[chrom]]
                fid.write(f"No_Cs: {' '.join(map(str, no_c))}\n")
                no_g = [int(x) if ~np.isnan(x) else "NaN" for x in self.no_g[chrom]]
                fid.write(f"No_Gs: {' '.join(map(str, no_g))}\n")
                no_t = [int(x) if ~np.isnan(x) else "NaN" for x in self.no_t[chrom]]
                fid.write(f"No_Ts: {' '.join(map(str, no_t))}\n")
                if len(self.g_to_a) != 0:
                    g_to_a = [int(x) if ~np.isnan(x) else "NaN" for x in self.g_to_a[chrom]]
                    fid.write(f"g_to_a: {' '.join(map(str, g_to_a))}\n")
                c_to_t = [int(x) if ~np.isnan(x) else "NaN" for x in self.c_to_t[chrom]]
                fid.write(f"c_to_t: {' '.join(map(str, c_to_t))}\n")
            fid.write("Reconstructed methylation:\n")
            for key in self.methylation.keys():
                if key == "win_size":
                    fid.write(f"\twin_size: {' '.join(map(str, self.methylation['win_size']))}\n")
                elif key == "slope":
                    fid.write(f"\tslope: {' '.join(map(str, self.methylation['slope']))}\n")
                elif key == "intercept":
                    fid.write(f"\tintercept: {' '.join(map(str, self.methylation['intercept']))}\n")
                elif key == "lcf":
                    fid.write(f"\tlcf: {self.methylation['lcf']}\n")
                elif key == "methylation":
                    for chrom in range(self.no_chrs):
                #meth = [round(x, 6) if ~np.isnan(x) else "NaN" for x in self.methylation["methylation"][chrom]]
                        meth = self.methylation["methylation"][chrom]
                        fid.write(f"\t{self.chr_names[chrom]}: {' '.join(map(str, meth))}\n")
            




if __name__ == "__main__":
    ams = Amsample(name="I1116", abbrev="1116")
    #mat = Amsample()
    #print(ams)
    #ams.diagnose()
    #ams2 = Amsample(name="First", coord_per_position=[2,2,2,2,2,2,2], no_t=[1,2,3], chr_names=["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7"])
    #print(ams2)
    #ams.parse_infile("../../u_1116.txt")
    #ams.parse_infile("data/python_dumps/I1116_drate.txt")
    #ams.parse_infile("data/python_dumps/I1116_meth.txt")
    #outfile = "objects/U1116"
    #t.save_object(outfile, ams)
    #infile = "objects/U1116"
    import cProfile
    #cProfile.run("ams = t.load_object(infile)", "data/logs/load_profile")
    
    #infile = "objects/U1116_diag"
    #infile = "objects/U1116_filtered"
    #ams = t.load_object(infile)
    #infile = "objects/bone_5"
    #mms = t.load_object(infile)
    #name = ams.name
    #print(f"name: {name}")
    #num = ams.no_chrs
    #print(f"num of chroms: {num}")
    #cProfile.run("ams.diagnose()", "data/logs/amsample_profile")
    #ams.diagnose()
    #ams.filter()
    #ams.estimate_drate(ref=mms)
    #stage = "drate"
    #ams.reconstruct_methylation()
    #meth_py = ams.methylation["methylation"]
    #meth_mat = mat.methylation["methylation"]
    #for chrom in range(ams.no_chrs):
    #    py = np.array(meth_py[chrom])
    #    diffs = py - meth_mat[chrom]
    #    diffs_nonan = [abs(x) for x in diffs if ~np.isnan(x)]
    #    res = max(diffs_nonan)
    #    print(f"Max methylation diff for {chrom}: {res}")
    #stage = "meth"
    #ams.dump(stage)
    #outfile = "objects/U1116_filtered_drate"
    #cProfile.run("t.save_object(outfile, ams)", "data/logs/save_profile")
    #t.save_object(outfile, ams)
    #base = ams.get_base_no("chr2", "c")
    #print(f"base: {base[0:20]}")
    #(no_t, no_ct) = ams.smooth("chr5", 17)
    #print(f"no_t: {no_t[0:20]}")
    #print(f"no_ct: {no_ct[0:25]}")
    #ams.simulate(0.018598, mms) #rate comes from drate global in I1116_meth.txt
    #ams.simulate(0, mms) #rate for testing
    #stage = "sim_0"
    #ams.dump(stage)
    chr_lengths = [249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,141213431,135534747,135006516,133851895,115169878,107349540,102531392,90354753,81195210,78077248,59128983,63025520,48129895,51304566,155270560,59373566]
    #cProfile.run("ams.bam_to_am(library='double', chr_lengths=chr_lengths, genome_seq='../../hg19.fa.gz', species='Homo sapiens')", "data/logs/bam_profile")
    #ams.bam_to_am(library="double", chr_lengths=chr_lengths, genome_seq="../../hg19.fa.gz", species="Homo sapiens", chroms=[15], tot_chroms=[15])
    ams.bam_to_am(library="double", chr_lengths=chr_lengths, genome_seq="../../hg19.fa.gz", species="Homo sapiens", trim_ends=True)
    stage = "bam"
    ams.dump(stage)

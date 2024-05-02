#!/usr/bin/python3

import tools as t
import numpy as np
import math
import matplotlib.pyplot as plt
import scipy.stats as stats
import copy
from chroms import Chrom
import pysam
import datetime
import glob
import re
#from config import *
import sys


class Amsample(Chrom):
    """Ancient methylation sample class

    This class inherits from Chroms superclass and has attributes name, species, reference,
    library, is_simulated, no_chrs, metadata, chr_names, coord_per_position, no_a, no_g, no_c, no_t,
    g_to_a, c_to_t, diagnostics, p_filters, methylation, and d_rate. The last 5 are
    dictionaries, while 8 of the previous 9 (metadata through c_to_t, excluding coord_per_position) are lists. 
    no_chrs is determined based on the length of chr_names. 
        
    An amsample object is created (with empty defaults): ams = Amsample(). The attributes can then be populated.
    """
    def __init__(self, name="unknown", species="unknown", reference="", library="", chr_names=[], coord_per_position="", no_a = [], no_c = [], no_g = [], no_t = [], g_to_a = [], c_to_t = [], diagnostics = {}, p_filters = {}, is_filtered = False, is_simulated = False, methylation={}, d_rate = {}, metadata=[]):
        self.name = name
        self.species = species
        self.reference = reference
        self.library = library
        self.chr_names = chr_names
        self.coord_per_position = coord_per_position
        self.no_a = copy.deepcopy(no_a)  # if these are just assigned (not copied), when creating 2 objects, these
        self.no_c = copy.deepcopy(no_c)  # elements will be identical and changing one will change the other
        self.no_g = copy.deepcopy(no_g)
        self.no_t = copy.deepcopy(no_t)
        self.g_to_a = copy.deepcopy(g_to_a)
        self.c_to_t = copy.deepcopy(c_to_t)
        self.diagnostics = copy.deepcopy(diagnostics)
        self.p_filters = copy.deepcopy(p_filters)
        self.is_filtered = is_filtered
        self.is_simulated = is_simulated
        self.methylation = copy.deepcopy(methylation)
        self.d_rate = copy.deepcopy(d_rate)
        self.metadata = copy.deepcopy(metadata)
        self.no_chrs = len(chr_names)

    def __repr__(self): #defines print of object
        return "name: %s\nspecies: %s\nreference: %s\nlibrary: %s\nchr_names: %s\ncoord_per_position: %s\nno_a: %s\nno_c: %s\nno_g: %s\nno_t: %s\ng_to_a: %s \nc_to_t: %s\ndiagnostics: %s\np_filters: %s\nis_filtered: %s\nis_simulated: %s\nmethylation: %s\nd_rate: %s\nmetadata: %s\nno_chrs: %s" % (self.name, self.species, self.reference, self.library, self.chr_names, self.coord_per_position, self.no_a, self.no_c, self.no_g, self.no_t, self.g_to_a, self.c_to_t, self.diagnostics, self.p_filters, self.is_filtered, self.is_simulated, self.methylation, self.d_rate, self.metadata, self.no_chrs)

    @staticmethod
    def find_indel(cig, seq, qual):
        """Adjusts sequence and quality values based on CIGAR string
        
        Input: CIGAR, query sequence, query quality (all lists)
        
        Output: new values for query sequence and query quality (lists)
        """
        types = {0:"M", 1:"I", 2:"D", 3:"N", 4:"S"}
        start = 0
        cig_vals = []
        seq_new = ""
        seq_old = seq
        qual_new = []
        for element in cig:
            cat = types[element[0]]
            val = element[1]
            cig_vals.append(str(val)+cat)
        cig_string = "".join(cig_vals)
        for val in cig_vals:
            cat = val[-1]
            val = int(val[:-1])
            if cat == "M":
                seq_new += seq[0:val]
                seq = seq[val:]
                qual_new.append(qual[:val])
                qual = qual[val:]
            elif cat == "S":
                seq_new = seq_old
                qual_new = [0 for x in range(len(seq_new))]  # throw out soft-clipped reads by making qual = 0 for each pos
                return(seq_new, qual_new)
            elif cat == "I":
                seq = seq[val:]
                qual = qual[val:]  # commented till temp bug removed
                continue
            elif cat == "D" or cat == "N":
                size = val
                seq_new += "0"*size
                qual_new.append([0 for x in range(size)])
        #if "I" in cig_string:
            #qual_new = qual_new[0:len(seq_new)]  # temporarily adding a bug to match matlab code
        qual_new = [x for y in qual_new for x in y]  # flatten list
        return(seq_new, qual_new)

    def process_bam(self, bam, chrom_name, chrom, gc, library, trim_ends, mapq_thresh, qual_thresh, chr_lengths, chrom_names, bam_chrom):
        """Goes over bam files to extract data
        
        """
        corrupt_reads = 0
        low_mapq_reads = 0
        reads_long_cigar = 0
        long_cigar_not_i_d = 0
        xflag = 0
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
        cpg_chrom = chrom_name
        cpg_plus = gc.coords[chrom]-1
        cpg_minus = gc.coords[chrom]
        for read in chrom_bam:
            if read.is_unmapped:
                print(f"Read {read.qname} on chrom {chrom_name} is unmapped. Skipping.")
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
            if read.is_duplicate:  # new in python script--matlab didn't skip duplicates
               continue
            if len(cig) != 1:
                (seq, qual) = self.find_indel(cig, seq, qual)
            if pos + len(seq) > chr_lengths[chrom]:
                print(f"Chrom {chrom_name} is {chr_lengths[chrom]} bp. Current read ({read.qname}) starts at {pos} and is {len(seq)} bp")
                continue
            if trim_ends:
                if strand == "+":
                    qual[0] = 0
                    qual[-2:] = [0,0]
                else:
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
            chr_a[0::2] = tot_a_minus[cpg_plus]
            chr_a[1::2] = tot_a_plus[cpg_minus]
            chr_g[0::2] = tot_g_minus[cpg_plus]
            chr_g[1::2] = tot_g_plus[cpg_minus]
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
        
        self.chr_names[chrom] = "chr"+chrom_names[bam_chrom]

    def bam_to_am(self, filename="", filedir=None, file_per_chrom=False, library=None, chr_lengths=None, species=None, chroms=list(range(23)), trim_ends=False, mapq_thresh=20, qual_thresh=20, gc_object=""):  # genome_seq is filename, orig qual_thresh=53, subtract 33 for 0 start
        """Converts bam files to Amsample objects
        
        Input: 
            filename           name and full path of file (if there's only 1)
            filedir            name of dir (for multiple files)
            file_per_chrom     True if there is exactly one file for each chromosome
            library            single or double
            chr_lengths        list of chromosome lengths in the same order as the chromosomes are given
            species
            chroms             list of chromosomes (defaults to chroms 1-22, X, Y)
            trim_ends          True to trim ends during processing, False if this has already been done
            mapq_thresh        threshold for read quality
            qual_thresh        threshold for read nucleotide quality
            gc_object          path for CpG file
        
        Output: Amsample object, populated with data from bam file
        """
        self.library = library
        self.species = species
        self.chr_names = [None]*(len(chroms))
        self.no_t = [[]]*(len(chroms))
        self.no_c = [[]]*(len(chroms))
        self.no_g = [[]]*(len(chroms))
        self.no_a = [[]]*(len(chroms))
        self.c_to_t = [[]]*(len(chroms))
        self.g_to_a = [[]]*(len(chroms))
        
        gc = t.load_object(gc_object)  # load CpG gcoordinates object
        
        if filedir:
            filenames = glob.glob(filedir+"/*.bam")
            if file_per_chrom:
                chrom_names = [x.split("_")[-1].split(".")[0] for x in filenames]  # filename format: <your label>_chr<chrom>.bam
                file_ref = {}
        if file_per_chrom:
            for i in range(len(chrom_names)):
                file_ref[chrom_names[i]] = filenames[i]
            i = 0
            for chrom in chroms:
                try:
                    bamfile = file_ref[chrom]  # take files in chrom order, regardless of filename order
                except:
                    print(f"No bam file for chromosome {chrom}. Please make sure chromosomes in list match filenames.")
                    sys.exit(1)
                bam = pysam.AlignmentFile(bamfile, "rb")
                all_chroms = bam.references
                chrom_key = t.build_chr_key(all_chroms)
                try:
                    chrom_num = chrom_key[chrom]
                except KeyError:
                    print(f"Chromosome {chrom} cannot be found in the bam file. Please remove it (and its length) from the config file.")
                    sys.exit(1)
                chrom_names = bam.references[0:chrom_num+1]  # lists names of all chroms up to max present in num order
                if self.chr_names[i]:
                    continue
                chrom_name = chrom_names[chrom_num]
                if bam.count(chrom_name) == 0:
                    continue
                chrom_pos = i
                self.process_bam(bam, chrom_name, chrom_pos, gc, library, trim_ends, mapq_thresh, qual_thresh, chr_lengths, chrom_names, chrom_num)
                bam.close()
                i += 1
        else:
            bam = pysam.AlignmentFile(filename, "rb")
            all_chroms = bam.references
            chrom_key = t.build_chr_key(all_chroms)
            chrom_index = []
            input_index = {}
            i = 0
            for chrom in chroms:
                try:
                    chrom_index.append(chrom_key[chrom])
                    input_index[chrom_key[chrom]] = i
                except KeyError:
                    print(f"Chromosome {chrom} cannot be found in the bam file. Please remove it (and its length) from the config file.")
                    sys.exit(1)
                i += 1
            chrom_names = bam.references[0:max(chrom_index)+1]  # lists names of all chroms up to max present in num order
            for chrom in chrom_index:
                chrom_name = chrom_names[chrom]
                chrom_pos = input_index[chrom]
                self.process_bam(bam, chrom_name, chrom_pos, gc, library, trim_ends, mapq_thresh, qual_thresh, chr_lengths, chrom_names, chrom)
            bam.close()  
        if len(self.no_t[0])/len(gc.coords[0]) == 2:      
            self.coord_per_position = "2"
        elif len(self.no_t[0])/len(gc.coords[0]) == 1:
            self.coord_per_position = "1"
        else:
            print(f"ratio of chrom coords to cpg coords is {len(self.no_t[1])/len(gc.coords[0])}")
        self.no_chrs = len(chroms)
        #add object name

    def parse_infile(self, infile):
        """Populate Amsample object from text file
        
        Input: empty Amsample object, file name
        Output: populated Amsample object
        """
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
                    if fields[0] == "Name" and (self.name == None or self.name == "unknown"):
                        self.name = fields[1]
                    elif fields[0] == "Species":
                        self.species = fields[1]
                    elif fields[0] == "Reference":
                        self.reference = fields[1]
                    elif fields[0] == "Library":
                        self.library = fields[1]
                    elif fields[0] == "Chromosomes":
                        self.chr_names = re.sub("[\[\]\']", "", fields[1]).split(", ")
                elif i == 5: #filtered line
                    if "False" in line:
                        self.is_filtered = False
                        #flag = 1 #no first set of chr lines for unfiltered except in new diagnose, fixed below
                    else:
                        self.is_filtered = True
                elif i == 6: #p filters line
                    if "False" in line:
                        flag = 1 #no first set of chr lines for unfiltered except in new diagnose, fixed below
                    #else:
                        #self.is_filtered = True
                elif i > 6:
                    if fields[0] == "Method" and i == 7:
                        self.p_filters["method"] = fields[1]
                        flag = 0  # fix for newer files with diagnostics for unfiltered
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
                        #self.is_simulated = False if "False" in fields[1] else True
                        self.is_simulated = False if "False" in line else True
                    elif fields[0] == "Coordinates per position":
                        coord = fields[1].split()
                        coord = [int(x) if x.isdigit() else np.nan for x in coord] #convert all numbers to int, leave NaN alone
                        self.coord_per_position = coord[0]
                    elif fields[0] == "Deamination rate":
                        if "True" in fields[1]:
                            self.d_rate["rate"] = {}
                        else:
                            continue
                    elif fields[0] == "Reference":
                        self.d_rate["ref"] = fields[1]
                    elif fields[0] == "Method":
                        self.d_rate["method"] = fields[1]
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
                        #self.chr_names.append(fields[0])
                        continue
                    elif fields[0] == "No_As" and "True" in fields[1]:
                        line = next(amfile)
                        no_a = line.split()
                        no_a = [int(x) if x.isdigit() else np.nan for x in no_a] #convert all numbers to int, leave NaN alone
                        self.no_a.append(no_a)
                    elif fields[0] == "No_Cs" and "True" in fields[1]:
                        line = next(amfile)
                        no_c = line.split()
                        no_c = [int(x) if x.isdigit() else np.nan for x in no_c] #convert all numbers to int, leave NaN alone
                        self.no_c.append(no_c)
                    elif fields[0] == "No_Gs" and "True" in fields[1]:
                        line = next(amfile)
                        no_g = line.split()
                        no_g = [int(x) if x.isdigit() else np.nan for x in no_g] #convert all numbers to int, leave NaN alone
                        self.no_g.append(no_g)
                    elif fields[0] == "No_Ts" and "True" in fields[1]:
                        line = next(amfile)
                        no_t = line.split()
                        no_t = [int(x) if x.isdigit() else np.nan for x in no_t] #convert all numbers to int, leave NaN alone
                        self.no_t.append(no_t)
                    elif (fields[0] == "aOga" or fields[0] == "g_to_a") and "True" in fields[1]:
                        line = next(amfile)
                        g_to_a = line.split()
                        g_to_a = [float(x) for x in g_to_a] #convert all numbers to float
                        g_to_a = np.array(g_to_a)
                        self.g_to_a.append(g_to_a)
                    elif (fields[0] == "tOct" or fields[0] == "c_to_t") and "True" in fields[1]:
                        line = next(amfile)
                        c_to_t = line.split()
                        c_to_t = [float(x) for x in c_to_t] #convert all numbers to float
                        self.c_to_t.append(c_to_t)
                    elif fields[0] == "Reconstructed methylation" and fields[1] == "True":
                        self.methylation["methylation"] = []
                        flag = 2 #done with second set of chr lines
                    elif fields[0] == "win_size":
                        win = fields[1].split()
                        win = [int(x) if x.isdigit() else np.nan for x in win] #convert all numbers to int, leave NaN alone
                        self.methylation["win_size"] = win
                    elif fields[0] == "algorithm":
                        self.methylation["algorithm"] = fields[1]
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
        self.no_chrs = len(self.chr_names)

    def get_base_no(self, chrom, base):
        """Get relevant list from object
        
        Input: chromosome (can be name or index), base (a, g, c, t, ga)
        Output: requested list for input chromosome
        """
        if isinstance(chrom, str):
            chr_ind = self.index([chrom])[0]
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
        """Smooths C and CT vectors
        
        Input: chromosome, window size
        Output: smoothed vectors
        """
        no_t = self.get_base_no(chrom, "t")
        no_c = self.get_base_no(chrom, "c")
        no_t = np.array(no_t) #convert to numpy array in order to add elementwise
        no_c = np.array(no_c)
        no_ct = no_t + no_c
        tpl = np.ones(winsize)
        no_t = t.nanconv(no_t, tpl, "same")
        no_ct = t.nanconv(no_ct, tpl, "same")
        return(no_t, no_ct)

    def region_methylation(self, region, gc, standardize=True):
        """Computes methylation in a specific region
        
        Input:  region   Genomic coordinates (chrom, start, end, delimited by :, -, or space(s))
                gc       Gcoordinates object of CpG positions
        Output: meth     methylation value in the region
        """
        if standardize:
            region = t.standardize_region(region) #take region input and change to dictionary
        chrom = region["chrom"] #get chrom from region dict
        chr_ind = gc.index([chrom])[0] #find index of region chrom in gc object

        cpg_start = np.where(gc.coords[chr_ind] >= region["start"])[0][0] #get index of first
        cpg_end = np.where(gc.coords[chr_ind] <= region["end"])[0][-1]    #and last coords in region
        chr_ind = self.index([chrom])[0] #find index of chrom in ams object
        no_t = np.nansum(self.no_t[chr_ind][cpg_start:cpg_end+1])
        no_ct = no_t + np.nansum(self.no_c[chr_ind][cpg_start:cpg_end+1])
        if self.methylation:
            meth = self.methylation["slope"][chr_ind] * no_t / no_ct + self.methylation["intercept"][chr_ind]
        else:  # used for DMRs (pooled_methylation) for samples that are pre-reconstruct_methylation
            meth = no_t/(no_ct * self.d_rate["rate"]["global"])  # slope = 1/global drate, ignore intercept
        if not np.isnan(meth):
            meth = min(max(meth,0),1)
        return meth

    def diagnose(self, fname=None, span=5, tolerance=1e-3, low_coverage=1, logdir=None, picdir=None):
        """Computes basic statistics on each input chromosome, and recommends what thresholds to use when excluding 
            PCR duplicates and true mutations.
            
        Input:  fname        name and path of output file
                span         parameter for outlier removal in the array no_ct
                tolerance    tolerance of BMM model convergence
                
        Output: Amsample object with updated diagnostics and p_filters fields
        """
        if fname is None:
            fname = logdir+self.name+"_diagnostics.txt"
        
        #initialize
        no_chrs = self.no_chrs
        #no_chrs = 3
        eff_coverage = np.zeros(no_chrs)
        coverage_threshold = np.zeros(no_chrs)
        c_to_t_threshold = [None] * no_chrs #set number of elements to add later by chrom index (use append instead?)
        fid = open(fname, "w")

        #compute number of rows and columns in output figures
        no_rows = math.floor(math.sqrt(no_chrs))
        if no_rows:
            no_cols = math.ceil(no_chrs/no_rows)
        #fig_coverage = plt.figure("coverage") 
        #fig_coverage.suptitle(self.name)
        #fig_perc_removed_per_coverage = plt.figure("removed_per_coverage");
        #fig_perc_removed_per_coverage.suptitle(self.name)
        #fig_perc_removed = plt.figure("removed")
        #fig_perc_removed.suptitle(self.name)
        
        #Loop on chromosomes
        time = datetime.datetime.now()
        time = time.strftime("%d-%m-%Y_%H_%M")
            
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
            if self.library == "single": #test this function on single stranded data
                no_a = self.no_a[chrom]
                g_to_a = self.g_to_a[chrom]
                if len(g_to_a) == 0:  
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
            N = np.histogram(nz_ct,bins=range(int(thresh)+2))
            more_to_remove = []
            
            imx = N[0].argmax()  # get index of N's max val
            points = [*range(imx-1,imx+2)]
            if points[0] == 0:
                points = [x+1 for x in points]
            p2 = np.polyfit(points,N[0][points],2)
            mu = -.5*p2[1]/p2[0]
            N0 = -.25*p2[1]**2/p2[0] + p2[2]
            imu = np.ceil(mu)
            idx = imu-1 + np.where(N[0][int(imu):-1]<.1*N0)[0][1]  # index 0 is 0 in python, so skip
            delta = idx - mu
            f = N[0][int(idx)]/N0
            sig = delta/np.sqrt(2*np.log(1/f))
            fid.write(f"\t     the distribution of [nt] is approximated as N(mu,sig)=N({mu},{sig})\n")
            #take the threshold as the first value where the expectation is to get less than one count
            total_counts = sum(N[0])
            thresh = np.ceil(mu + sig*np.sqrt(2*np.log(total_counts/sig/np.sqrt(2*np.pi))))
            fid.write(f"\t     Final threshold was set to {thresh}\n")
            more_to_remove = np.where(no_ct>thresh)[0]
            fid.write(f"\t     {len(more_to_remove):,d} positions ({100*len(more_to_remove)/len(no_ct):,.2f}%) removed as nt > {thresh:.0f}\n")
            
            coverage_threshold[chrom] = thresh
            
            #plot the coverage
            plt.figure("coverage", figsize=[15, 10]) #plots chr1 in sep window
            plt.suptitle(self.name)
            plt.subplot(no_rows,no_cols,chrom+1) #chrom num, not index num ##plots rest of chroms in 1 window, chr labeling off by 1
            
            if self.chr_names:
                xlabel = self.chr_names[chrom]
            else:
                xlabel = f"Chromosome #{chrom+1}"
            
            plt.xlabel(xlabel)
            N = plt.hist(nz_ct,bins=range(int(thresh)), label="coverage") #creates plots for coverage
            #plt.bar(len(N[0]), N[0])
            plt.plot([thresh, thresh], list(plt.gca().get_ylim()), "k")
            plt.subplots_adjust(wspace=.5, hspace=.5)
            pic_file = picdir + self.name + "_cov_" + time
            plt.savefig(pic_file)
        
            #remove outliers
            to_remove = [to_remove, more_to_remove]
            for x in to_remove:
                for y in x:
                    no_ct[y] = np.nan
            for x in to_remove:
                for y in x:
                    no_t[y] = np.nan
            
            #effective coverage
            eff_coverage[chrom] = np.nanmean(no_ct)
            fid.write(f"\t     Effective coverage: {eff_coverage[chrom]:.1f}\n")
            
            #max_c = float(max_coverage)
            
            #compare to previous PCR-duplicate threshold     REMOVE THESE LINES?
            #prev_to_remove = np.where(no_ct>max_c)[0]
            #fid.write(f"\t     COMPARISON: {len(prev_to_remove):,d} positions ({100*len(prev_to_remove)/len(no_ct):.2f}%) ")
            #fid.write(f"would have been removed if the PCR-duplicate threshold were {max_coverage}\n")
            
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
                p[0] = max(p[0],0.001)
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
                fid.write(f"\t\t\tand to include {w[1]*no_pos_cover*stats.binom.cdf(thresh-1,cover+1,0.5)+w[2]*no_pos_cover*stats.binom.cdf(thresh-1,cover+1,p[2]):.1f} false positives\n")
                if thresh == cover+1:
                    #what happens with no threshold?
                    fid.write("\t\tIf no threshold is used (no positions removed), ")
                    fid.write(f"then {(w[1]+w[2])*no_pos_cover:.1f} false positives will be retained\n")
                
                #remove the positions
                for x in more_to_remove:
                    no_ct[x] = np.nan #assumes 1d array
                for x in more_to_remove:
                    no_t[x] = np.nan
            
            c_to_t_threshold[chrom] = th_c2g

            #plot ratio of removed to total per coverage
            quotient = no_pos_removed/no_pos_tot*100
            plt.figure("removed_per_coverage", figsize=[15, 10])
            plt.subplot(no_rows,no_cols,chrom+1) #chrom num, not index num
            plt.plot(chr_cov,quotient,"b.")
            if self.chr_names:
                xlabel = self.chr_names[chrom]
            else:
                xlabel = f"Chromosome #{chrom+1}"
            plt.xlabel(xlabel, fontsize=9)
            ylabel = "%removed(coverage)/\ntotal(coverage)"
            plt.ylabel(ylabel, fontsize=9)
            plt.grid(True)
            plt.subplots_adjust(wspace=.5, hspace=.5)
            pic_file = picdir + self.name + "_rpc_" + time
            plt.savefig(pic_file)

            #plot ratio of removed to total
            tot = sum(no_pos_tot)
            plt.figure("removed", figsize=[15, 10])
            plt.subplot(no_rows,no_cols,chrom+1)
            plt.plot(chr_cov,no_pos_removed/tot*100,"b.")
            if self.chr_names:
                xlabel = self.chr_names[chrom]
            else:
                xlabel = f"Chromosome #{chrom+1}"
            plt.xlabel(xlabel)
            ylabel = "%removed(coverage)/total"
            plt.ylabel(ylabel)
            plt.grid(True)
            plt.subplots_adjust(wspace=.5, hspace=.5)
            pic_file = picdir + self.name + "_rem_" + time
            plt.savefig(pic_file, bbox_inches="tight")
            
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
        """If a CpG position is marked to be removed, both consecutive CpG positions are marked to be removed.
        
        Input:  positions to be removed
        Output: modified list of positions to be removed
        """
        odds = np.where(to_remove%2)
        evens = np.where(~to_remove%2)
        to_remove = list(to_remove) + list(to_remove[odds]-1) + list(to_remove[evens]+1)
        to_remove = np.unique(to_remove)
        return to_remove

    def filter(self, max_c_to_t = 0.25, merge = True, max_g_to_a = .25, method = None, use_max_TsPerCoverage = True, max_coverage = None,  fname = None, min_t = 1, max_a = 1, logdir=None):
        """Removes information from CpG sites that did not pass various quality control tests
        
        Input:    max_c_to_t         threshold used to identify sites with a true C->T mutation. All positions 
                    where c_to_t >= max_c_to_t are removed. This parameter can be one of the following:
                    (1) scalar, in which case it is used for all chromosomes, and all coverage levels.
                    (2) vector over chromosomes, with a different threshold for each chromosome (and possibly for 
                    each coverage level).
                  min_t              the filter 'max_c_to_t' is applied only to positions where no_t > min_t.
                  merge              whether to merge the two consecutive coordinates of every CpG position.
                  fname              output (log) file name.
                  max_g_to_a         a threshold used to identify sites with a true C->T mutation. Applicable only
                    when the library field is 'single'. All positions where g_to_a >= max_g_to_a are removed. 
                    A rule of thumb to determine it is to set max_g_to_a = 2 * (1/(coverage/2)), otherwise, 
                    take the same value as 'max_c_to_t'. This parameter can be one of the following:
                    (1) scalar, in which case it is used for all chromosomes, and all coverage levels.
                    (2)vector over chromosomes, with a different threshold for each chromosome.
                  max_a              all positions where no_a > max_a are removed. Applicable only when the 
                    library field is 'single'.
                  method             method used to remove true C->T mutations. It can be one of the following:
                    (1) 'c_to_t' uses information from no_t and no_c only. Uses the parameter 'max_c_to_t' above.
                    (2) 'both' combines c_to_t option, above, with g_to_a, using information from no_t, no_c, no_a, 
                    and no_g. Uses the parameters 'max_c_to_t', 'max_g_to_a' and 'max_a' above. Applicable only when the
                    library field is 'single'.
                  max_coverage       a coverage threshold used to remove PCR duplicates. All positions for which
                    no_c > max_coverage are removed. This parameter can be one of the following:
                    (1) scalar, in which case it is used for all chromosomes.
                    (2)vector over chromosomes, with a different threshold for each chromosome.
                  use_max_TsPerCoverage  max_TsPerCoverage is a threshold used to identify sites with a true C->T 
                    mutation. All positions with a certain coverage no_t+no_c=C, where no_t > max_TsPerCoverage(C) 
                    are removed. This parameter is an array over chromosomes, with an array of max_TsPerCoverage for 
                    each coverage level in each chromosome. Default is to use this over user-entered max_c_to_t.
        Output: Amsample object with removed sites replaced with NaNs
        """
        #initialize
        no_chr = self.no_chrs
        if fname == None:
            fname = logdir+self.name+"_filter_log.txt"
        if max_coverage == None:
            max_coverage = self.p_filters["max_coverage"] if any(self.p_filters["max_coverage"]) else 100
        if method == None or not method:
            if self.library == "single": #test this function on single stranded data
                method = "both"
            else:
                method = "c_to_t"
        if use_max_TsPerCoverage:
            max_TsPerCoverage = self.p_filters["max_TsPerCoverage"] if self.p_filters["max_TsPerCoverage"] else max_c_to_t
        else:
            max_TsPerCoverage = max_c_to_t
        #bring input parameters into standard form - max_coverage
        if np.isscalar(max_coverage):
            max_coverage = max_coverage * np.ones(no_chr)
        #bring input parameters into standard form - max_TsPerCoverage
        if np.isscalar(max_TsPerCoverage):
            tmp_max = max_TsPerCoverage
        else:
            tmp_max = max_TsPerCoverage.copy() 
        if not use_max_TsPerCoverage: #max_c_to_t
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
            print(f"({self.coord_per_position} coordinates per position): {no_pos:,d}")
            fid.write("\tNumber of CpG positions ")
            fid.write(f"({self.coord_per_position} coordinates per position): {no_pos:,d}\n")
            #remove positions whose coverage is too high
            to_remove = np.where(no_ct>max_coverage[chrom])[0]
            no_removed = len(to_remove)
            fid.write(f"\t{no_removed:,d} positions ({100*no_removed/no_pos:.2f}%) removed as No_CTs > {int(max_coverage[chrom])}\n")
            #loop on each coverage level and remove positions with high No_Ts
            if method != "g_to_a":
                #loop over all coverage levels
                max_t_int = [int(x) if ~np.isnan(x) else np.nan for x in max_TsPerCoverage[chrom]]
                for cover in range(int(max_coverage[chrom]),0,-1):
                    idx = np.where(no_ct==cover)[0]
                    fid.write(f"\tcoverage {cover} ({len(idx):,d} positions):")
                    more_to_remove = idx[no_t[idx]>max_TsPerCoverage[chrom][cover-1]] 
                    more_to_remove = self.extend_removed(more_to_remove)
                    no_removed = len(to_remove) #redefined for every coverage level
                    to_remove = np.unique(list(to_remove)+list(more_to_remove))
                    no_removed = len(to_remove) - no_removed
                    fid.write(f"\t{no_removed:,d} (extended) positions ")
                    fid.write(f"were removed as No_Ts > {max_t_int[cover-1]}\n")
            #remove more positions if data on A's and G's is available
            if method != "c_to_t":
                more_to_remove = np.union1d(np.intersect1d(np.where(no_a<=max_a[chrom]), np.where(g_to_a>=max_g_to_a[chrom])), np.where(no_a>max_a[chrom]))
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
                if f_ams.coord_per_position == "1": #test on single coord per pos 
                    if f_ams.chr_names:
                        print(f"{f_ams.chr_names[chrom]} had already gone through merger")
                    else:
                        print(f"Chromosome #{chrom} had already gone through merger")
                else:
                    #f_ams.coord_per_position = "1"
                    if f_ams.library == "double":
                        operation = "max"
                    else:
                        operation = "sum"
                    no_t = t.nanmerge(no_t, operation)
                    no_c = t.nanmerge(no_c, operation)
            #substitute into f_ams
            f_ams.no_t[chrom] = no_t
            f_ams.no_c[chrom] = no_c
            f_ams.diagnostics["effective_coverage"][chrom] = (np.nansum(no_t) + np.nansum(no_c))/len(no_t)
            if f_ams.no_a: f_ams.no_a[chrom] = []
            if f_ams.no_g: f_ams.no_g[chrom] = []
            if f_ams.g_to_a: f_ams.g_to_a[chrom] = []
            if f_ams.c_to_t: f_ams.c_to_t[chrom] = []
            tot_removed = tot_removed + no_removed
        if merge and f_ams.coord_per_position != "1":
            f_ams.coord_per_position = "1"
        f_ams.is_filtered = True
        print(f"In total {tot_removed:,d} positions were removed")
        #close file
        fid.close()

    def estimate_drate(self, method="reference", global_meth=np.nan, min_cov=1, ref=[], min_beta=1):
        """Estimates deamination rate
        
        Input: method        the method used to perform the estimation. Can be one of:
                 'reference' the preferred and most accurate method, which should be used when we have a vector of 
                  beta-values as a reference.
                 'global' should be used in cases we have no measured reference, and the estimation is based
                  on the pre-assumption of the total level of methylation in the genome.
               global_meth   the estimated value of global genomic methylation, either as a fraction or as a
                  percentage. Applicable for 'method'='global'.
               min_cov       minimum coverage of sites that are used for the estimation.
               ref           Mmsample object containing the beta-values of the reference. Applicable for 
                  'method'='reference'.
               min_beta      minimum beta value to take for the estimation. Applicable for 'method'='reference'.
        Output: Amsample object with udpated 'd_rate' field. This field is a dictionary with keys:
                  'rate', a dictionary that contains:
                    'global' global deamination rate (computed from all chromosomes)
                    'dglobal' STD of the global deamination rate (assuming binomial distribution)
                    'local' local deamination rate, computed for each chromosome separately
                    'dlocal' STD of the local deamination rate (assuming binomial distribution)
                    'no_positions' in each chromosome, number of positions upon which deamination rate was computed
                  and, based in the input, 'method', 'global_methylation', 'ref' (the name of the Mmsample object),
                  'min_beta', and 'min'coverage'.
        """
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
        else:
            raise Exception(f"Unknown method '{method}'")
        print(f"Estimating deamination rate using the '{method}' method")
        if method == "global":
            print(f"\tglobal_methylation: {meth_params['global_methylation']:.2f}")
        elif method == "reference":
            print(f"\tmin_beta: {meth_params['min_beta']:.2f}")
        #if method == "global" or method == "reference":
        print(f"\tmin_coverage: {meth_params['min_coverage']:d}")
        
        #sanity check
        #if not self.is_filtered:
        #    raise Exception(f"{self.name} is not filtered")
        if self.coord_per_position != "1" and self.coord_per_position != 1 and 1 not in self.coord_per_position:  # handle string and int (to cover both old and new objects)
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
            if re.search("[Mm]",self.chr_names[chrom]):
                continue
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
                if self.chr_names[chrom] not in ref.chr_names:
                    print(f"Chromosome {self.chr_names[chrom]} is not in reference sample")
                    continue
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
        drate["dglobal"] = factor * math.sqrt(drate["global"] * (1-drate["global"])/np.nansum(drate["no_positions"]))

        #plug vals into object
        if method == "reference":
            self.d_rate = {"method":"reference", "rate":drate, "ref":meth_params["ref"].name, "min_beta":meth_params["min_beta"], "min_coverage":meth_params["min_coverage"]}
        elif method == "global":
            self.d_rate = {"method":"global", "rate":drate, "global_methylation":meth_params["global_methylation"], "min_coverage":meth_params["min_coverage"]}
    
    def determine_winsize(self, chrom, method="prob", min_meth=0.2, p0=0.01, k_recip=1/2.5, max_width=31):
        """Estimates the optimal window size in cases that collecting data from a window is required.
        
        Input: chrom        index of chromosome
               method       either 'prob' (for probability) or 'relerror' (for relative error).
               min_meth     minimum methylation level we want to detect.
               p0           applicable if 'method'='prob'. It is the probability to get zero counts in the 
                 window if the methylation is min_meth.
               k_recip        applicable if 'method'='relerror'. It is one over the maximum relative error in
                 estimating the methylation in a window whose true methylation is min_meth. It also means that the
                 mean is far (k standard deviations) from zero.
               max_width    maximum allowed width. Computed window is not allowed to be larger than 'max_width'.
        Output: win_size    recommended window size, forced to be an odd number.
        """
        coverage = self.diagnostics["effective_coverage"][chrom]
        drate=self.d_rate["rate"]["global"]
        if not method=="prob" and not method=="relerror":
            raise Exception(f"Unknown method '{method}'")
        if min_meth>1:
            min_meth = 0.01 * min_meth
        
        #compute the window size
        p = min_meth * drate
        if method == "prob":
            win_size = np.ceil(np.log(p0)/np.log(1-p)/coverage)
        else:
            win_size = np.ceil((1-p)/coverage/p/k_recip**2)

        #narrow window if too wide
        win_size = min(win_size, max_width)

        #make sure win_size is odd
        if not win_size%2:
            win_size += 1
        
        return win_size
    
    def reconstruct_methylation(self, win_size="auto", winsize_alg={}, function="histogram", slope=None, intercept=[0], ref=[], lcf=0.05):
        """Computes methylation from c_to_t data, based on some function of the C->T ratio (no_t/no_ct).
        
        Input: win_size        window size for smoothing. If 'auto', a recommended value is computed for each 
                 chromosome. Otherwise, it can be a list with a single value (where the value is used for all 
                 chromosomes) or a value for each chromosome.
               winsize_alg     a dictionary with parameters required to determine window size, see parameters 
                 for determine_winsize.
               function        to compute methylation as a function of the C->T ratio. Options are:
                 'histogram', where the function is computed by histogram matching to a reference methylome.
                 'linear', where the function is: meth = slope * no_t / no_ct + intercept
                 'logistic', where the function is: meth = tanh(slope * no_t / no_ct).
               slope           a parameter used for the 'linear' and the 'logistic' functions. It can be a list with a 
               single value (where the value is used for all chromosomes) or a value for each chromosome.
               intercept       a parameter for used for the 'linear' function, and determines the intercept of the
                 linear transformation. Can be a list with a single value (where the value is used for 
                 all chromosomes) or a value for each chromosome.
               ref             Mmsample object containing the beta-values of the reference.
               lcf             low coverage factor.
               
        Output: Amsample object with udpated methylation field.
        """
        no_chr = self.no_chrs
        if slope == None or slope == "":
            slope=[1/self.d_rate["rate"]["global"]]
        if win_size == "auto":
            auto_win = True
        elif isinstance(win_size, list):
            auto_win = False
        else:
            win_size = re.split(",| ", win_size)
            win_size = [int(x) for x in win_size]
            auto_win = False
        
        #bring parameters into standard format - win_size
        if not auto_win:
            if len(win_size) == 1:
                if not win_size[0]%2:
                    win_size[0] += 1 #win_size should be odd
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

        #bring parameters into standard format - slope
        if len(slope) == 1:
            slope = float(slope[0]) * np.ones(no_chr)

        #bring parameters into standard format - intercept
        if len(intercept) == 1:
            intercept = float(intercept[0]) * np.ones(no_chr)
        intercept = [float(x) for x in intercept]
        
            
        meth = []
        for chrom in range(no_chr):
            if self.chr_names:
                print(f"Computing methylation in {self.chr_names[chrom]}")
            else:
                print(f"Computing methylation in chrom #{chrom+1}")  # does this work with partial chrom list?
            #get smoothed No_Ts and No_CTs
            (no_t, no_ct) = self.smooth(chrom, int(win_size[chrom]))
            #remove regions with particularly low coverage
            lct = t.find_low_coverage_thresh(no_ct, lcf)
            no_ct = [np.nan if x < lct else x for x in no_ct]
            #compute methylation
            c_to_t = no_t/no_ct
            if function == "histogram":
                ref.merge()
                ref.scale()
                # hard-coded parameters
                ref_bins = 100
                sig_bins = 1000
                # use only finite elements
                idx_finite = np.where(np.isfinite(c_to_t))
                sig = [c_to_t[x] for x in idx_finite[0]]
                #sig = c_to_t[idx_finite]  # !!
                # x-axis for reference and signal
                ref_edges = np.linspace(0,1,ref_bins+1)
                sig_edges = np.linspace(0,max(sig),sig_bins+1)
                ref_binwidth = ref_edges[1] - ref_edges[0]
                # smooth the reference
                vec = ref.smooth(None,win_size[chrom],name=self.chr_names[chrom])[0]
                vec = [vec[x] for x in idx_finite[0]]
                # generate histograms
                hist = np.histogram(vec, bins=ref_edges)[0]
                N = np.sum(hist)
                c = np.cumsum(hist)
                N_ref = c/N
                hist = np.histogram(sig, bins=sig_edges)[0]
                N = np.sum(hist)
                c = np.cumsum(hist)
                N_sig = c/N
                # generate the mapping
                hmap = np.zeros(sig_bins)
                for i in range(sig_bins):
                    # find closest value on the CDF
                    imin = np.argmin(abs(N_ref - N_sig[i]))
                    hmap[i] = ref_edges[imin] + 0.5*ref_binwidth
                # make the range precisely [0,1]
                hmap = (hmap - hmap[0]) / np.ptp(hmap)
                # apply the mapping
                sig1 = np.digitize(sig,sig_edges)
                sig1[np.where(sig1 == len(sig_edges))] = len(sig_edges) - 1
                methi = np.empty(len(c_to_t))
                methi[:] = np.nan
                tmp = [hmap[x-1] for x in sig1]
                d = dict(zip(idx_finite[0], tmp))
                for i in idx_finite[0]:
                    methi[i] = d[i]
                            
            elif function == "logistic" or function == "log":
                methi = np.tanh(slope[chrom]*c_to_t)
            elif function == "linear" or function == "lin":
                methi = slope[chrom]*c_to_t+intercept[chrom]
                methi = np.minimum(np.maximum(methi,0),1) #keep between 0 and 1 (nans untouched)
            print(f"Average methylation: {np.nanmean(methi):.2f}")
            meth.append(methi)
        self.methylation = {"methylation":meth, "algorithm":function, "win_size":win_size, "slope": slope, "intercept":intercept, "lcf":lcf}

    def simulate(self, mms, chrom=None):
        """Simulates Cs and Ts of ancient DNA based on the degradation rate and coverage of Amsample object,
            assuming a methylation given by Mmsample object.
        
        Input: mms            Mmsample object that contains the reference methylation.
               
        Output: Amsample object with modified no_t and no_c values, based on simulation.
        """
        #initialize
        # if "sim" not in self.name:  # not sure why this was done
            # self.name += "__sim_"
        self.is_simulated = True

        #make sure methylation in {mms} is scaled to [0,1]
        mms.scale()
        mms.merge()
        degrad_rate = self.d_rate["rate"]["global"]
        #parameters needed for the simulation
        meth_map = mms.get_methylation(chrom=chrom)[1]

        for chrom in range(self.no_chrs):
            print(f"Processing {self.chr_names[chrom]} ... ")
            no_c = np.array(self.no_c[chrom])
            tot_ct = no_c + self.no_t[chrom]
            ct_nan_idx = np.argwhere(np.isnan(tot_ct))
            meth_nan_idx = np.argwhere(np.isnan(meth_map[chrom]))
            tot_nan_idx = np.concatenate((ct_nan_idx, meth_nan_idx))
            uniq_nan_idx = set(tot_nan_idx.flatten())
            ct_float = [0 if i in uniq_nan_idx else tot_ct[i] for i in range(len(tot_ct))]
            #ct_int = [int(x) for x in ct_float] #nec?
            meth_float = [0 if i in uniq_nan_idx else meth_map[chrom][i] for i in range(len(meth_map[chrom]))]
            #no_t = np.random.binomial(ct_int, degrad_rate*np.array(meth_float))  # np.array?
            no_t = np.random.binomial(ct_float, degrad_rate*np.array(meth_float))  # np.array?
            no_t_fl = no_t.astype(float)
            no_t_fl[list(uniq_nan_idx)] = np.nan
            self.no_t[chrom] = no_t_fl
            self.no_c[chrom] = tot_ct - self.no_t[chrom]
            print("done")

    def dump(self, stage, chroms=None, dir="", bed=True, gc_object="", object_dir=""):
        """Dumps Amsample object to text file.
        
        Input: stage    name indicating which part of the process has been dumped, for use in file name.
               chroms   indices of chromsomes to dump. 
               dir      output directory, specified in input params or config file
        Output: text file in format <object_name>_<stage>.txt (directory currently hard-coded).
        """
        aname = self.name
        aname = re.sub("\s", "_", aname)
        fname = dir + aname + "_" + stage + ".txt"
        with open(fname, "w") as fid:
            fid.write(f"Name: {aname}\nSpecies: {self.species}\nReference: {self.reference}\n")
            fid.write(f"Library: {self.library}\n")
            fid.write(f"Chromosomes: {self.chr_names}\n")
            #filt = "" if self.is_filtered else "not "
            #fid.write(f"This sample is {filt}filtered\n")
            fid.write(f"Is filtered: {self.is_filtered}\n")
            fid.write("P filters: ")
            if self.p_filters:
                fid.write("True\n")
            else:
                fid.write("False\n")
            for key in self.p_filters.keys():
                if key == "method":
                    fid.write(f"\tMethod: {self.p_filters['method']}\n")
                elif key == "max_coverage":
                    max_cov = [int(x) if ~np.isnan(x) else "NaN" for x in self.p_filters['max_coverage']]
                    fid.write(f"\tmax_coverage: {' '.join(map(str, max_cov))}\n")
                elif key == "max_TsPerCoverage":
                    fid.write("\tmax_TsPerCoverage:\n")
                    for chrom in range(self.no_chrs):
                        max_t = [int(x) if ~np.isnan(x) and ~np.isinf(x) else "NaN" for x in self.p_filters['max_TsPerCoverage'][chrom]]
                        fid.write(f"\t\t{self.chr_names[chrom]}: {' '.join(map(str, max_t))}\n")
                elif key == "max_g_to_a":
                    fid.write(f"\tmax_g_to_a: {' '.join(map(str, self.p_filters['max_g_to_a']))}\n")
                elif key == "max_a":
                    max_a = [int(x) if ~np.isnan(x) else "NaN" for x in self.p_filters['max_a']]
                    fid.write(f"\tmax_No_As: {' '.join(map(str, max_a))}\n")
            #sim = "" if self.is_simulated else "not "
            #fid.write(f"This sample is {sim}simulated\n")
            fid.write(f"Is simulated: {self.is_simulated}\n")
            fid.write(f"Coordinates per position: {self.coord_per_position}\n")
            fid.write("Deamination rate: ")
            if self.d_rate:
                fid.write("True\n")
            else:
                fid.write("False\n")
            for key in self.d_rate.keys():
                if key == "ref":
                    fid.write(f"\tReference: {self.d_rate['ref']}\n")
                elif key == "method":
                    fid.write(f"\tMethod: {self.d_rate['method']}\n")
                elif key == "min_beta":
                    fid.write(f"\tmin_beta: {int(self.d_rate['min_beta'])}\n")
                elif key == "min_coverage":
                    fid.write(f"\tmin_coverage: {int(self.d_rate['min_coverage'])}\n")
                elif key == "rate":
                    for subkey in self.d_rate["rate"].keys():
                        if subkey == "global":
                            fid.write(f"\tglobal: {self.d_rate['rate']['global']}\n")
                        elif subkey == "dglobal":
                            fid.write(f"\tdglobal: {self.d_rate['rate']['dglobal']}\n")
                        elif subkey == "local":
                            fid.write(f"\tlocal: {' '.join(map(str, self.d_rate['rate']['local']))}\n")
                        elif subkey == "dlocal":
                            fid.write(f"\tdlocal: {' '.join(map(str, self.d_rate['rate']['dlocal']))}\n")
                        elif subkey == "no_positions":
                            fid.write(f"\tno_positions: {' '.join(map(str, self.d_rate['rate']['no_positions']))}\n")
            fid.write("Diagnostics: ")
            if self.diagnostics:
                fid.write("True\n")
            else:
                fid.write("False\n")
            for key in self.diagnostics.keys():
                if key == "effective_coverage":
                    fid.write(f"\tEffective coverage: {' '.join(map(str, self.diagnostics['effective_coverage']))}\n")
            if not chroms:
                chroms = range(self.no_chrs)
            for chrom in chroms:
                fid.write(f"{self.chr_names[chrom]}:\n")
                no_a = [int(x) if ~np.isnan(x) else "NaN" for x in self.no_a[chrom]] if self.no_a else 0
                val = f"True\n{' '.join(map(str, no_a))}" if no_a else "False"
                fid.write(f"No_As: {val}\n")
                no_c = [int(x) if ~np.isnan(x) else "NaN" for x in self.no_c[chrom]] if self.no_c else 0
                val = f"True\n{' '.join(map(str, no_c))}" if no_c else "False"
                fid.write(f"No_Cs: {val}\n")
                #if self.no_g:  # single stranded from matlab doesn't have
                no_g = [int(x) if ~np.isnan(x) else "NaN" for x in self.no_g[chrom]] if self.no_g else 0
                val = f"True\n{' '.join(map(str, no_g))}" if no_g else "False"
                fid.write(f"No_Gs: {val}\n")
                no_t = [int(x) if ~np.isnan(x) else "NaN" for x in self.no_t[chrom]] if self.no_t else 0
                val = f"True\n{' '.join(map(str, no_t))}" if no_t else "False"
                fid.write(f"No_Ts: {val}\n")
                #if len(self.g_to_a) != 0:
                g_to_a = [float(x) for x in self.g_to_a[chrom]] if self.g_to_a else 0
                val = f"True\n{' '.join(map(str, g_to_a))}" if g_to_a else "False"
                fid.write(f"g_to_a: {val}\n")
                #if self.c_to_t:  # single stranded from matlab doesn't have
                c_to_t = [float(x) for x in self.c_to_t[chrom]] if self.c_to_t else 0
                val = f"True\n{' '.join(map(str, c_to_t))}" if c_to_t else "False"
                fid.write(f"c_to_t: {val}\n")
            fid.write("Reconstructed methylation: ")
            if self.methylation:
                fid.write("True\n")
            else:
                fid.write("False\n")
            for key in self.methylation.keys():
                if key == "win_size":
                    fid.write(f"\twin_size: {' '.join(map(str, self.methylation['win_size']))}\n")
                elif key == "algorithm":
                    fid.write(f"\talgorithm: {self.methylation['algorithm']}\n")
                elif key == "slope":
                    fid.write(f"\tslope: {' '.join(map(str, self.methylation['slope']))}\n")
                elif key == "intercept":
                    fid.write(f"\tintercept: {' '.join(map(str, self.methylation['intercept']))}\n")
                elif key == "lcf":
                    fid.write(f"\tlcf: {self.methylation['lcf']}\n")
                elif key == "methylation":
                    if bed:
                        gc = t.load_object(gc_object)
                        meth_bed = {}
                    #for chrom in range(self.no_chrs):
                    for chrom in range(len(self.methylation["methylation"])):
                        if bed:
                            meth_bed[chrom] = []
                            for i in range(len(self.methylation["methylation"][chrom])):
                                meth_bed[chrom].append(f"{self.chr_names[chrom]}\t{gc.coords[chrom][i]}\t{gc.coords[chrom][i]+1}\t{self.methylation['methylation'][chrom][i]}\n")
                        meth = self.methylation["methylation"][chrom]
                        fid.write(f"\t{self.chr_names[chrom]}: {' '.join(map(str, meth))}\n")
                    if bed:
                        bed_fname = dir + aname + "_" + stage + ".bed"
                        with open(bed_fname, "w") as bid:
                            for chrom in range(len(meth_bed)):
                                for line in range(len(meth_bed[chrom])):
                                    bid.write(meth_bed[chrom][line])
        if stage == "meth":
            outfile = object_dir + aname
            t.save_object(outfile, self)




if __name__ == "__main__":
    #ams = Amsample(name="I1116", abbrev="1116")
    #ams = Amsample(name="Ust_Ishim", abbrev="Ust")
    ams = Amsample(name="Altai_Neanderthal")

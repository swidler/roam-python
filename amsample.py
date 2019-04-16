#!/usr/bin/python3

import tools as t
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from chroms import Chrom

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
                        self.chr_names.append(fields[0])
                    elif fields[0] == "max_aOga":
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
                        self.d_rate["local"] = []
                        self.d_rate["dlocal"] = []
                        self.d_rate["no_positions"] = []
                    elif fields[0] == "Reference":
                        self.d_rate["ref"] = fields[1]
                    elif fields[0] == "min_beta":
                        self.d_rate["min_beta"] = float(fields[1])
                    elif fields[0] == "min_coverage":
                        self.d_rate["min_coverage"] = float(fields[1])
                    elif fields[0] == "global":
                        self.d_rate["global"] = float(fields[1])
                    elif fields[0] == "dglobal":
                        self.d_rate["dglobal"] = float(fields[1])
                    elif fields[0] == "local":
                        local = fields[1].split()
                        local = [float(x) for x in local] #convert all numbers to float (NaN is a float)
                        self.d_rate["local"] = local
                    elif fields[0] == "dlocal":
                        dlocal = fields[1].split()
                        dlocal = [float(x) for x in dlocal] #convert all numbers to float (NaN is a float)
                        self.d_rate["dlocal"] = dlocal
                    elif fields[0] == "no_positions":
                        no_positions = fields[1].split()
                        no_positions = [float(x) for x in no_positions] #convert all numbers to float (NaN is a float)
                        self.d_rate["no_positions"] = no_positions
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
                    elif fields[0] == "aOga":
                        g_to_a = fields[1].split()
                        g_to_a = [int(x) if x.isdigit() else np.nan for x in g_to_a] #convert all numbers to int, leave NaN alone
                        self.g_to_a.append(g_to_a)
                    elif fields[0] == "tOct":
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
        chr_ind = self.indexofchr([chrom])[0]
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
        c_to_t_threshold = []
        fid = open(fname, "w")

        #compute number of rows and columns in output figures
        no_rows = math.floor(math.sqrt(no_chrs))
        if no_rows:
            no_cols = math.ceil(no_chrs/no_rows)
        fig_coverage = plt.figure("coverage")
        fig_coverage.suptitle(self.name)
        fig_perc_removed_per_coverage = plt.figure("removed_per_coverage");
        fig_perc_removed_per_coverage.suptitle(self.name)
        fig_perc_removed = plt.figure("removed")
        fig_perc_removed.suptitle(self.name)
        
        #Loop on chromosomes
        for chrom in range(no_chrs):
            #report
            if self.chr_names:
                print(f"Diagnosing chrom {self.chr_names[chrom]}")
                fid.write(f"Diagnosing chrom {self.chr_names[chrom]}\n")
            else:
                print(f"Diagnosing chrom #{chrom+1}")
                fid.write(f"Diagnosing chrom #{chrom+1}\n")
            
            #vectors of current chromosome
            no_t = self.no_t[chrom]
            no_t = np.array(no_t) #convert to numpy array in order to add elementwise
            no_c = self.no_c[chrom]
            no_c = np.array(no_c)
            no_ct = no_t + no_c
            if self.library == "SS":
                no_a = self.no_a[chrom]
                g_to_a = self.g_to_a[chrom]
                if not g_to_a:
                    no_g = self.no_g[chrom]
                    no_g = np.array(no_g) #convert to numpy array in order to add elementwise
                    no_a = np.array(no_a)
                    if no_g:
                        g_to_a = no_a / (no_a + no_g)
            fid.write(f"\tNumber of CpG positions: {len(no_ct):,d}\n")
            fid.write(f"\tInitial coverage: {np.nanmean(no_ct)}\n")
            
            #compute the upper threshold of {nt} in a few steps
            #crude outlier removal using very loose criterion
            percentiles = t.nan_prctile(no_ct, [25,50,75])
            fid.write(f"\t[C+T]: median = {percentiles[1]}\n")
            fid.write(f"\t     [25th 75th] percentiles = {percentiles[[0,2]]}")
            iqrange = percentiles[2] - percentiles[0]
            thresh = percentiles[2] + span*iqrange
            fid.write(f"=> initial threshold was set to {thresh}\n")
            to_remove = np.where(no_ct > thresh)[0]
            fid.write(f"\t     {len(to_remove):,d} positions ({100*len(to_remove)/len(no_ct):,.2f}%) removed as no_ct > {thresh}\n")
            no_ct = [float(x) for x in no_ct]
            for x in to_remove:
                no_ct[x] = np.nan
            no_ct = np.array(no_ct)
            
            #evaluate the parameters of the normal distribution of {nt}
            nt_zeros = np.where(no_ct == 0)[0]
            nz_ct = np.copy(no_ct) #copy array
            for x in nt_zeros:
                nz_ct[x] = np.nan
            N = plt.hist(nz_ct,bins=range(int(thresh)))
            more_to_remove = []
            if strict:
                print("in strict loop") #handle later
            coverage_threshold[chrom] = thresh
            
            #plot the coverage
            plt.figure("coverage")
            plt.subplot(no_rows,no_cols,chrom+1) #chrom num, not index num
            plt.bar(len(N[0]), N[0])
            plt.plot([thresh, thresh], list(plt.gca().get_ylim()), "k")
            if self.chr_names:
                xlabel = self.chr_names[chrom]
            else:
                xlabel = f"Chromosome #{chrom+1}"
            plt.xlabel(xlabel)
            
            #remove outliers
            to_remove = [to_remove, more_to_remove] #more_to_remove gets vals in skipped strict loop
            for x in to_remove:
                for y in x:
                    no_ct[y] = np.nan
            no_t = [float(x) for x in no_t]
            for x in to_remove:
                for y in x:
                    no_t[y] = np.nan
            
            #effective coverage
            eff_coverage[chrom] = np.nanmean(no_ct)
            fid.write(f"\t     Effective coverage: {eff_coverage[chrom]}\n")
            max_c = float(max_coverage)
            
            #compare to previous PCR-duplicate threshold
            prev_to_remove = np.where(no_ct>max_c)[0]
            fid.write(f"\t     COMPARISON: {len(prev_to_remove):,d} positions ({100*len(prev_to_remove)/len(no_ct)}%) ")
            fid.write(f"would have been removed if the PCR-duplicate threshold were {max_coverage}\n")
            
            #report on the zero bin
            fid.write(f"\t{len(nt_zeros):,d} positions ({100*len(nt_zeros)/len(no_ct)}%) have no_ct = 0\n")

            #analyze each coverage level independently
            th_c2g = np.zeros(int(thresh))
            chr_cov = [x for x in range(low_coverage, int(thresh+1))]
            no_pos_tot = np.zeros(int(thresh+1))
            no_pos_removed = np.zeros(int(thresh))
            for cover in range(int(thresh),low_coverage, -1):
                fid.write(f"\tAnalyzing [xt] coming from coverage {cover}\n")
                idx = np.where(no_ct==cover)[0]
                no_pos_cover = len(idx)
                no_pos_tot[cover] = no_pos_cover
                fid.write(f"\t\tTotal number of positions: {no_pos_cover:,d}\n")
                H = plt.hist(np.array(no_t)[idx],bins=range(cover))
                p, w, llike = t.bmm(H, [0.01, 0.5, 0.99], [0.9, 0.05, 0.05], tolerance, [0, 1, 0])
                fid.write("\t\tEstimated homozygous mutations: ")
                fid.write(f"Pr(meth) = {p[2]}; Weight = {w[2]}\n")
                fid.write("\t\tEstimated heterozygous mutations: ")
                fid.write(f"Pr(meth) = {p[1]}; Weight = {w[1]}\n")
                fid.write("\t\tEstimated deaminated positions: ")
                fid.write(f"Pr(meth) = {p[0]}; Weight = {w[0]}\n")




if __name__ == "__main__":
    ams = Amsample()
    print(ams)
    #ams.diagnose()
    #ams2 = Amsample(name="First", coord_per_position=[2,2,2,2,2,2,2], no_t=[1,2,3], chr_names=["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7"])
    #print(ams2)
    #ams.parse_infile("../../u_1116.txt")
    #outfile = "objects/U1116"
    #t.save_object(outfile, ams)
    infile = "objects/U1116"
    ams = t.load_object(infile)
    name = ams.name
    print(f"name: {name}")
    num = ams.no_chrs
    print(f"num of chroms: {num}")
    ams.diagnose()
    base = ams.get_base_no("chr2", "c")
    print(f"base: {base[0:20]}")
    (no_t, no_ct) = ams.smooth("chr5", 17)
    print(f"no_t: {no_t[0:20]}")
    print(f"no_ct: {no_ct[0:25]}")
    


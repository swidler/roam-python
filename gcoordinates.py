#!/usr/bin/python3

from chroms import Chrom
import numpy as np
import tools as t
import gzip
from Bio import SeqIO

class Gcoordinates(Chrom):
    def __init__(self, name="", description="", species="unknown", reference="", chr_names=[], coords=[], strand=[], metadata=[], metadata_name=[]):
        self.name = name
        self.description = description
        self.species = species
        self.reference = reference
        self.chr_names = chr_names
        self.no_chrs = len(chr_names)
        #for chrom in range(no_chrs)):
        #    self.coords[i] = np.array(self.coords[i])
        #for chrom in range(self.no_chrs):
        strand.extend([[]]*self.no_chrs)
        coords.extend([[]]*self.no_chrs)
        self.strand = strand
        self.coords = coords
        self.metadata = metadata
        self.metadata_name = metadata_name
        

    def __repr__(self): #defines print of object
        return "name: %s\ndescription: %s\nspecies: %s\nreference: %s\nmetadata: %s\nchr_names: %s\ncoords: %s\nstrand: %s\nno_chrs: %s\nmetadata_name: %s" % (self.name, self.description, self.species, self.reference, self.metadata, self.chr_names, self.coords, self.strand, self.no_chrs, self.metadata_name)

    def parse_infile(self, infile):
        """Parses a Gcoordinates object from a text file.
        
        Input: self    empty Gcoordinates object
               infile  input filename and path
        Output: populated Gcoordinates object
        """
        i = 0
        temp = []
        with open(infile, "rt") as gcfile:
            for line in gcfile:
                line = line.rstrip("\n") #remove line feeds
                if i == 0:
                    header = line #first line is header
                    self.description = header
                    fields = header.split(" ")
                    ref = fields[2]
                    ref = ref.rstrip(")")
                    ref = ref.lstrip("(")
                    self.reference = ref
                    self.species = fields[0]
                    self.name = fields[1]
                elif i > 1: #2nd line is blank, rest are chroms
                    fields = line.split("\t")
                    self.chr_names.append(fields[0])
                    coords = fields[1].split(",")
                    coords = [int(x) if x.isdigit else np.nan for x in coords] #convert numbers to int, leave nans
                    temp.append(np.asarray(coords))
                i += 1
        self.coords = np.asarray(temp)
        self.no_chrs = len(self.coords) #reassign chrom num based on new info
        
    def create_cpg_file(self, genome_file, outfile):
        """Builds pickled CpG for reference genome
        
        Output: pickled file in format <species_name>_cpg_coords.P in object dir
        """
        with gzip.open(genome_file, "rt") as fas:
            records = list(SeqIO.parse(fas, "fasta"))
            
        # Map fasta record ids to sequences once (avoid scanning records repeatedly)
        rec_by_id = {r.id: str(r.seq).upper() for r in records}
        for chrom, chrom_name in enumerate(self.chr_names):

            # Try both "<name>" and "chr<name>" conventions, plus chrM/chrMT weirdness
            candidates = [chrom_name, "chr" + chrom_name]

            # Handle common mitochondrial naming mismatch
            # If user asked for chrMT, accept fasta's chrM
            if chrom_name in ("MT", "chrMT"):
                candidates += ["chrM", "MT", "chrMT"]
                
            seq = None
            for cid in candidates:
                if cid in rec_by_id:
                    seq = rec_by_id[cid]
                    break

            if seq is None:
                raise ValueError(
                    f"Chromosome '{chrom_name}' not found in FASTA. Tried ids: {candidates}"
                )

            # Find CpG: positions where seq[i]=='C' and seq[i+1]=='G'
            # Store 1-based positions of the C (like your original code)
            cpg_plus = [i + 1 for i in range(len(seq) - 1) if seq[i] == "C" and seq[i + 1] == "G"]

            self.coords[chrom] = np.array(cpg_plus, dtype=np.int32)

        t.save_object(outfile, self)


    def calc_tss(self, genes):
        """Computes Transcription Start Sites (TSSs)
        
        Input: genes    bed file with all genes (no duplicates)
        Output: Gcoordinates object updated with TSS info
        """
        self.metadata_name.append("UCSC_name")
        for chrom in range(self.no_chrs):
            tss = []
            strands = []
            names = []
            metadata = []
            for ivl in genes:
                norm_chrom = t.normalize_chrom(ivl.chrom)
                if norm_chrom is not None and norm_chrom == chrom + 1: # chrom index val is one less than chrom num
                #if int(ivl.chrom) == chrom+1:  # chrom index val is one less than chrom num
                    start = ivl.start
                    end = ivl.end
                    strand = ivl.strand
                    name = ivl.name
                    UCSC_name = ivl.score
                    if strand == '+':
                        pos = start
                        strand = 1
                    else:
                        pos = end
                        strand = 0
                else:
                    continue    
                tss.append(pos)
                strands.append(strand)
                names.append(name)
                metadata.append(UCSC_name)
            self.coords[chrom] = np.asarray(tss)
            self.strand[chrom] = np.asarray(strands)
            self.metadata.append(metadata)

#!/usr/bin/python3

import tools as t
import numpy as np
import copy

class cDMR:
    """Chromosome-specific DMRs
    
    This class has attributes chromosome (name), CpG_start, CpG_end, gen_start, gen_end, no_bases, no_CpGs, max_Qt, 
    methylation, annotation, all of which are lists, other than the chromosome. no_DMRs is determined based on the length
    of max_Qt.
    
    A cDMR object is created (with empty defaults): cdm = cDMR(). The attributes can then be populated.
    """
    def __init__(
        self,
        chromosome="",
        CpG_start=None,
        CpG_end=None,
        gen_start=None,
        gen_end=None,
        no_bases=None,
        no_CpGs=None,
        max_Qt=None,
        grp_methylation_statistic=None,
        methylation=None,
        annotation=None,
    ):
        self.chromosome = chromosome
        CpG_start = [] if CpG_start is None else CpG_start
        CpG_end = [] if CpG_end is None else CpG_end
        gen_start = [] if gen_start is None else gen_start
        gen_end = [] if gen_end is None else gen_end
        no_bases = [] if no_bases is None else no_bases
        no_CpGs = [] if no_CpGs is None else no_CpGs
        max_Qt = [] if max_Qt is None else max_Qt
        grp_methylation_statistic = (
            [] if grp_methylation_statistic is None else grp_methylation_statistic
        )
        methylation = [] if methylation is None else methylation
        annotation = [] if annotation is None else annotation
        self.CpG_start = copy.deepcopy(CpG_start)  # if these are just assigned (not copied), when creating 2 objects, these
        self.CpG_end = copy.deepcopy(CpG_end)  # elements will be identical and changing one will change the other
        self.gen_start = copy.deepcopy(gen_start)
        self.gen_end = copy.deepcopy(gen_end)
        self.no_bases = copy.deepcopy(no_bases)
        self.no_CpGs = copy.deepcopy(no_CpGs)
        self.max_Qt = copy.deepcopy(max_Qt)
        self.methylation = copy.deepcopy(methylation)
        self.grp_methylation_statistic = copy.deepcopy(grp_methylation_statistic)
        self.annotation = copy.deepcopy(annotation)
        self.no_DMRs = len(self.max_Qt)

    def __repr__(self): #defines print of object
        return "chromosome: %s\nCpG_start: %s\nCpG_end: %s\ngen_start: %s\ngen_end: %s\nno_bases: %s\nno_CpGs: %s\nmax_Qt: %s\nmethylation:  %s\ngrp_methylation_statistic: %s\nannotation: %s\nno_DMRs: %s" % (self.chromosome, self.CpG_start, self.CpG_end, self.gen_start, self.gen_end, self.no_bases, self.no_CpGs, self.max_Qt, self.methylation, self.grp_methylation_statistic, self.annotation, self.no_DMRs)

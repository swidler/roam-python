#!/usr/bin/python3

from chroms import Chrom
import numpy as np
import tools as t
import gcoordinates as gcoord
import perscoords as pers
import argparse
import os

parser = argparse.ArgumentParser(description="Create a personalized CpG coordinate file for a sample based on lists of positions to add and remove from an existing coords file")

parser.add_argument("name", help="sample name")
parser.add_argument("object_dir", help="directory for saved (pickled) object files (include final /)")
parser.add_argument("toadd_file", help="path to file with positions to add (chrom \\t position)")
parser.add_argument("torem_file", help="path to file with positions to remove (chrom \\t position)")
parser.add_argument("gc_file", help="CpG file")

args = parser.parse_args()

name=args.name
print("Creating a personalized gc file for ", name)

def create_personalize_gc(gc_file, name, object_dir, toadd_file, torem_file):
    # Ref gc file
    gc = t.load_object(gc_file)
    
    # Position to add
    toadd = pers.pers_coords()
    toadd.parse_posinfile(toadd_file)
    
    # Position to remove
    torem = pers.pers_coords()
    torem.parse_posinfile(torem_file)
    
    gc.personalize_gc(name=name, object_dir=object_dir, toadd=toadd, torem=torem)


create_personalize_gc(gc_file=args.gc_file, name=name, object_dir=args.object_dir, toadd_file=args.toadd_file, torem_file=args.torem_file)

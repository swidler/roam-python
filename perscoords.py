#!/usr/bin/python3

from chroms import Chrom
import numpy as np
import tools as t


class pers_coords(Chrom):
    def __init__(self, name="", chr_names=[], coords=[]):
        self.name = name
        self.chr_names = chr_names
        self.coords = coords

    def __repr__(self): #defines print of object
        return "name: %s\nchr_names: %s\ncoords: %s" % (self.name, self.chr_names, self.coords)

    def parse_posinfile(self, poslist):
        """Parses a pers_coords object from a text files.
        
        Input: self   empty Gcoordinates object
               poslist  input pos_to_add/remove filename and path
        Output: populated pers_coords object
        """
        # Initialize an empty dictionary to store positions to add
        posdict = {}
        with open(poslist, "rt") as file:
            for line in file:
                line = line.strip()
                chrom, position = line.split("\t")
                # If the chromosome is already in the dictionary, append the position
                if chrom in posdict:
                    posdict[chrom].append(int(position))  # positions should be integers
                else:
                    # Otherwise, create a new list with the current position
                    posdict[chrom] = [int(position)]

        self.chr_names = list(posdict.keys())
        self.coords = np.asarray([np.asarray(x) for x in list(posdict.values())], dtype=object)
        del posdict 

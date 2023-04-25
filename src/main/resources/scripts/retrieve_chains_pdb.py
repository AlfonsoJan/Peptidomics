#!/usr/bin/env python3

import MDAnalysis as mda
import numpy as np
import json
import sys

__author__ = "Wouter Zeevat"

def get_chains(path):
    u = mda.Universe(path)
    selection = 'protein and (name N or name CA or name C or name O)'
    proteins = u.select_atoms(selection)
    prev = count = {protein:0 for protein in np.unique(proteins.segids)}


    for protein in proteins:
        chain = protein.segid
        if prev[chain] != protein.resid:
            prev[chain] = protein.resid
            count[chain] += 1
    print(json.dumps({p:str(count[p]) for p in count}))

def main(args):
    get_chains(args[1])


if __name__ == "__main__":
    sys.exit(main(sys.argv))
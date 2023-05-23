#!/usr/bin/env python3

"""
Module that retrieves the chains from a pdb file AS JSON data
"""

import json
import sys
import MDAnalysis as mda
import numpy as np

__author__ = "Wouter Zeevat"

def get_chains(path):
    """
    This function returns the chains of the pdb file
    """

    protein = mda.Universe(path)
    selection = 'protein and (name N or name CA or name C or name O)'
    proteins = protein.select_atoms(selection)
    prev = count = {protein:0 for protein in np.unique(proteins.segids)}

    # Counts all the unique residues count
    for protein in proteins:
        chain = protein.segid
        if prev[chain] != protein.resid:
            prev[chain] = protein.resid
            count[chain] += 1
    print(json.dumps({p:str(count[p]) for p in count}))

def main(args):
    """
    Main function that calls everything
    """
    get_chains(args[1])


if __name__ == "__main__":
    sys.exit(main(sys.argv))

#!/usr/bin/env python3

import MDAnalysis as mda
import numpy as np
import sys

def read_pdb(file, split_number):
    """
    Reads pdb and spits it out!
    """

    # Read a protein in
    u = mda.Universe(file)
    #u = mda.Universe('../Proteins/5fq5.pdb')

    # Select only (protein) backbone
    selection = 'protein and (name N or name CA or name C or name O)'
    protein = u.select_atoms(selection)
    coordinates = protein.positions

    # Split on breaks. Both segment/chain breaks and
    # within-chain breaks need to be accounted for.
    # The following assumes the atom order N-CA-C-O
    parts = np.split(
        coordinates, #             N[1:]               C[:-1]
        4 * np.where(((coordinates[4::4] - coordinates[2:-3:4]) ** 2).sum(axis=1) > 3)[0] + 4
    )
    n = int(split_number) * 4 # aminoacids * backbone atoms
    P = np.array([p[i:i+n] for p in parts for i in range(0, len(p)-n, 4)])
    P -= P.mean(axis=1)[:, None, :]
    return P

def main(args):
    np.save(f'{args[2]}.npy', read_pdb(args[1], args[3]))


if __name__ == "__main__":
    sys.exit(main(sys.argv))

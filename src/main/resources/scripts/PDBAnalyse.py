#!/usr/bin/env python3

import sys
import numpy as np
from collections import Counter
import json

__author__ = "Wouter Zeevat"

AA1 = 'A C D E F G H I K L M N P Q R S T V W Y'.split()
AA3 = 'ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR'.split()
AA321 = dict(zip(AA3, AA1))

def read_pdb_backbone(pdbfile):
    with open(pdbfile) as pdb:
        backbone = [
            line.strip() for line in pdb
            if line.startswith('ATOM') and line[12:16] in (' N  ', ' CA ', ' C  ')
        ]

    chainsnos = [ atom[20:22].strip() for atom in backbone ]
    residues = [ atom[17:20] for atom in backbone ]
    atomnames = [ atom[12:16].strip() for atom in backbone ]
    coordinates = np.loadtxt([ atom[30:54] for atom in backbone ])
    atomnos = [ atom[4:12].strip() for atom in backbone ]
    return residues, coordinates, atomnos, atomnames, chainsnos

def peptidize(pdbfile, pepsize):
    '''
    Read peptides from given size from PDB file
    '''
    sequence, coordinates, atomnos, atomnames, chainsnos = read_pdb_backbone(pdbfile)

    # Because each residue has three atoms
    pepsize *= 3

    # Separate parts of chains
    breaks = np.where(((coordinates[1:] - coordinates[:-1]) ** 2).sum(axis=1) > 4)[0] + 1
    # For the coordinates
    parts = np.split(coordinates, breaks)
    # For the atom numbers, peptides and chain letters
    breaks = [0, *breaks.tolist(), None]
    partatomnos = [ atomnos[a:b] for a, b in zip(breaks[:-1], breaks[1:]) ]
    partseq = [ sequence[a:b] for a, b in zip(breaks[:-1], breaks[1:]) ]
    partchains = [ chainsnos[a:b] for a, b in zip(breaks[:-1], breaks[1:]) ]
    # Coordinates per part
    pepcoords = np.array([
        part[i:i + pepsize]
        for part in parts
        for i in range(len(part) - pepsize)
    ])
    pepcoords -= pepcoords.mean(axis=1)[:,None,:]
    # Peptides per part
    peptides = [
        Counter(part[i:i+pepsize:3]).most_common()[0][0]
        for part in partseq
        for i in range(len(part) - pepsize)
    ]
    # Atom numbers per part
    atomnos = [
        [min(part[i:i + pepsize], key=int), max(part[i:i + pepsize], key=int)]
        for part in partatomnos
        for i in range(len(part) - pepsize)
    ]
    # Chain letters per part
    chainsnos = [
        Counter(part[i:i+pepsize:3]).most_common()[0][0]
        for part in partchains
        for i in range(len(part) - pepsize)
    ]
    pep_information = list(zip(peptides, atomnos, chainsnos))
    return pep_information, pepcoords

def princana(coordinates):
    distances = dist_function(coordinates)
    values, vectors = np.linalg.eigh(np.cov(distances.T))
    # These are the characteristic structure values for each peptide
    scores = distances @ vectors[:, ::-1][:, :3]

    return scores

def dist_function(coordinates):
    distances = ((coordinates[:, :, None] - coordinates[:, None, :]) ** 2).sum(axis=3)
    rows, columns = np.triu_indices_from(distances[0], 1)
    distances = distances[:, rows, columns].reshape((len(distances), -1))
    distances -= distances.mean(axis=0)
    return distances

def parse_to_json(scores, peptide_information):
    print(json.dumps({
        idx: {
            "peptide": val[0],
            "atomnos": {"min": val[1][0], "max": val[1][1]},
            "chain": val[2][0],
            "x": scores[:, 0][idx],
            "y": scores[:, 1][idx],
            "z": scores[:, 2][idx],
        }
        for idx, val in enumerate(peptide_information)
    }))

"""
Main function that calls everything
"""
def main(args):
    filename = args[1]
    oligo_length = int(args[2])
    peptide_information, pepcoords = peptidize(filename, oligo_length)
    scores = princana(pepcoords)
    parse_to_json(scores, peptide_information)


if __name__ == "__main__":
    sys.exit(main(sys.argv))

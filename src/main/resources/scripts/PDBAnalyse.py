#!/usr/bin/env python3

import sys
import numpy as np
from collections import Counter
import json
import requests
import time

__author__ = "Wouter Zeevat"

AA1 = 'A C D E F G H I K L M N P Q R S T V W Y'.split()
AA3 = 'ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR'.split()
AA321 = dict(zip(AA3, AA1))
LETTERS = 'H B E G I T S P None'.split()
MEANING = 'ALPHA-HELIX BETA-BRIDGE BETA-LADDER 3-HELIX 5-HELIX HYDROGEN-BONDED-TURN BEND PLOY-PROLINE-HELICES LOOP/IRREGULAR'.split()
LETTERS_DICT = dict(zip(LETTERS, MEANING))

def start_job(pdb_code):
    form_data = {'data': pdb_code}
    url_create_code = f"https://www3.cmbi.umcn.nl/xssp/api/create/pdb_id/dssp/"
    r = requests.post(url_create_code, data=form_data)
    r.raise_for_status()
    job_id = json.loads(r.text)['id']
    return job_id

def get_secondary_structure(pdb_code):
    job_id = start_job(pdb_code)
    ready = False
    error = False
    while not ready:
        url_check_status = f"https://www3.cmbi.umcn.nl/xssp/api/status/pdb_id/dssp/{job_id}"
        r = requests.get(url_check_status)
        r.raise_for_status()
        status = json.loads(r.text)['status']
        if status == 'SUCCESS':
            ready = True
        elif status in ['FAILURE', 'REVOKED']:
            error = True
            ready = True
        else:
            time.sleep(5)
    if error:
        return []
    url_result = f"https://www3.cmbi.umcn.nl/xssp/api/result/pdb_id/dssp/{job_id}"
    r = requests.get(url_result)
    r.raise_for_status()
    result = json.loads(r.text)['result']
    result = result.split("\n")[:-1]
    chains = []
    res_numbers = []
    codes = []
    header = False
    for line in result:
        if header:
            chains.append(line[10:12].strip())
            res_numbers.append(line[5:10].strip())
            code = line[16]
            if code == " ":
                code = "None"
            codes.append(LETTERS_DICT[code])
        if line.lstrip().startswith("#"):
            header = True
    return list(zip(chains, res_numbers, codes))

def structure_parser(chainsnos, res_number, result_dssp):
    if result_dssp.size == 0:
        return ["undefined"] * len(res_number)
    structures = []
    t = np.array(list(zip(chainsnos, res_number)))
    for index in range(0, t.shape[0]):
        chain = t[index][0]
        res_num = int(t[index][1])
        result_index = np.where((result_dssp[:,:-1] == t[index]).all(axis=1))[0][0]
        structures.append(result_dssp[result_index][2])
    return np.array(structures)

def read_pdb_backbone(pdbfile, pdb_code):
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
    residues_number = [ atom[22:27].strip() for atom in backbone ]
    structures = structure_parser(chainsnos, residues_number, np.array(get_secondary_structure(pdb_code)))
    return residues, coordinates, atomnos, atomnames, chainsnos, structures

def peptidize(pdbfile, pepsize, pdb_code):
    '''
    Read peptides from given size from PDB file
    '''
    sequence, coordinates, atomnos, atomnames, chainsnos, structures = read_pdb_backbone(pdbfile, pdb_code)

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
    partstructure = [ structures[a:b] for a, b in zip(breaks[:-1], breaks[1:]) ]
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
    # Structure per part
    pep_structures = np.array([
        Counter(part[i:i+pepsize]).most_common()[0][0]
        for part in partstructure
        for i in range(len(part) - pepsize)
    ])
    pep_information = list(zip(peptides, atomnos, chainsnos, pep_structures))
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
            "structure": val[3],
        }
        for idx, val in enumerate(peptide_information)
    }))

"""
Main function that calls everything
"""
def main(args):
    filename = args[1]
    oligo_length = int(args[2])
    pdb_code = args[3]
    peptide_information, pepcoords = peptidize(filename, oligo_length, pdb_code)
    scores = princana(pepcoords)
    parse_to_json(scores, peptide_information)


if __name__ == "__main__":
    sys.exit(main(sys.argv))

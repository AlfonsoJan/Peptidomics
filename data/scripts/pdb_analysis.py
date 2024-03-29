#!/usr/bin/env python3

"""Performs the analysis creating the JSON coordinates together with the metadata"""

import json
import time
from pathlib import Path
import sys
import numpy as np
import requests

__author__ = "John Busker, Tsjerk Wassenaar & Wouter Zeevat"

LETTERS = 'H B E G I T S P None'.split()
MEANING = 'ALPHA-HELIX BETA-BRIDGE BETA-LADDER 3-HELIX 5-HELIX HYDROGEN-BONDED-TURN ' \
          'BEND POLY-PROLINE-HELICES LOOP/IRREGULAR'.split()
LETTERS_DICT = dict(zip(LETTERS, MEANING))
API_CALL_COUNTER_MAX = 2

# https://stackoverflow.com/questions/26646362/numpy-array-is-not-json-serializable
class NumpyEncoder(json.JSONEncoder):
    """
    Encoder for JSON data
    """
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)

def start_job(pdb_code):
    """
    Starts the job and creates a json object
    :param pdb_code:
    :return:
    """
    form_data = {'data': pdb_code}
    url_create_code = "https://www3.cmbi.umcn.nl/xssp/api/create/pdb_id/dssp/"
    request = requests.post(url_create_code, data=form_data)
    request.raise_for_status()
    job_id = json.loads(request.text)['id']
    return job_id

def get_secondary_structure(pdb_code):
    """
    Gets the pdb file from a code and saves the important data
    :param pdb_code:
    :return:
    """
    api_call_counter = 0
    job_id = start_job(pdb_code)
    ready = False
    error = False
    while not ready:
        url_check_status = f"https://www3.cmbi.umcn.nl/xssp/api/status/pdb_id/dssp/{job_id}"
        request = requests.get(url_check_status)
        request.raise_for_status()
        status = json.loads(request.text)['status']
        if status == 'SUCCESS':
            ready = True
        elif status in ['FAILURE', 'REVOKED'] or api_call_counter == API_CALL_COUNTER_MAX:
            error = True
            ready = True
        else:
            api_call_counter += 1
            time.sleep(5)
    if error:
        return []
    url_result = f"https://www3.cmbi.umcn.nl/xssp/api/result/pdb_id/dssp/{job_id}"
    request = requests.get(url_result)
    request.raise_for_status()
    result = json.loads(request.text)['result']
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

def get_middle_peptide(cords, pepsize):
    """
    Returns the middle 3 items of a list, and gets the middle 2 on an even number length list
    :param cords:
    :param pepsize:
    :return:
    """
    index = int((len(cords) - 1)/2)
    if pepsize / 3 == 1:
        return cords
    if ((pepsize / 3) % 2) == 0:
        return cords[index: index+2]
    return cords[index-1:index+2]

def structure_parser(chainsnos, res_number, result_dssp):
    """
    Parses the structure by using the chain number and the res number.

    :param chainsnos:
    :param res_number:
    :param result_dssp:
    :return:
    """
    if result_dssp.size == 0:
        return ["undefined"] * len(res_number)
    structures = []
    pep_info = np.array(list(zip(chainsnos, res_number)))
    for index in range(0, pep_info.shape[0]):
        result_index = np.where((result_dssp[:,:-1] == pep_info[index]).all(axis=1))[0][0]
        structures.append(result_dssp[result_index][2])
    return np.array(structures)

def read_pdb_backbone(pdbfile, pdb_code):
    """
    Reads the backbone from the pdb structure
    :param pdbfile:
    :param pdb_code:
    :return:
    """
    with open(pdbfile) as pdb:
        backbone = [
            line.strip() for line in pdb
            if line.startswith('ATOM') and line[12:16] in (' N  ', ' CA ', ' C  ')
        ]

    chainsnos = [ atom[20:22].strip() for atom in backbone ]
    residues = [ atom[17:20] for atom in backbone ]
    coordinates = np.loadtxt([ atom[30:54] for atom in backbone ])
    atomnos = [ atom[4:12].strip() for atom in backbone ]
    residues_number = [ atom[22:27].strip() for atom in backbone ]
    structures = structure_parser(chainsnos, residues_number,
                                  np.array(get_secondary_structure(pdb_code)))
    return residues, coordinates, atomnos, chainsnos, structures

def peptidize(pdbfile, pepsize, pdb_code):
    """
    Read peptides from given size from PDB file
    """
    residues, coordinates, atomnos, chainsnos, structures = \
        read_pdb_backbone(pdbfile, pdb_code)

    # Because each residue has three atoms
    pepsize *= 3

    # Separate parts of chains
    breaks = np.where(((coordinates[1:] - coordinates[:-1]) ** 2)
                      .sum(axis=1) > 4)[0] + 1
    # For the coordinates
    parts = np.split(coordinates, breaks)
    # For the atom numbers, peptides and chain letters
    breaks = [0, *breaks.tolist(), None]
    partatomnos = [ atomnos[a:b] for a, b in zip(breaks[:-1], breaks[1:]) ]
    partres = [ residues[a:b] for a, b in zip(breaks[:-1], breaks[1:]) ]
    partchains = [ chainsnos[a:b] for a, b in zip(breaks[:-1], breaks[1:]) ]
    partstructure = [ structures[a:b] for a, b in zip(breaks[:-1], breaks[1:]) ]
    # Coordinates per part
    pepcoords = np.array([
        part[i:i + pepsize]
        for part in parts
        for i in range(0, len(part) - pepsize, 3)
    ])
    pepcoords -= pepcoords.mean(axis=1)[:,None,:]
    # Peptides per part
    pepresidues = [
        get_middle_peptide(part[i:i + pepsize:3], pepsize)
        for part in partres
        for i in range(0, len(part) - pepsize, 3)
    ]
    # Atom numbers per part
    atomnos = [
        [min(part[i:i + pepsize], key=int), max(part[i:i + pepsize], key=int)]
        for part in partatomnos
        for i in range(0, len(part) - pepsize, 3)
    ]
    # Chain letters per part
    chainsnos = np.array([
        get_middle_peptide(part[i:i + pepsize:3], pepsize)
        for part in partchains
        for i in range(0, len(part) - pepsize, 3)
    ])
    # Structure per part
    pep_structures = np.array([
        get_middle_peptide(part[i:i + pepsize:3], pepsize)
        for part in partstructure
        for i in range(0, len(part) - pepsize, 3)
    ])
    pep_information = list(zip(pepresidues, atomnos, chainsnos, pep_structures))
    return pep_information, pepcoords

def princana(coordinates, vectors):
    """
    Performs principal component analysis!
    :param coordinates:
    :param vectors:
    :return:
    """
    distances = dist_function(coordinates)
    # These are the characteristic structure values for each peptide
    projected_data = distances @ vectors
    return projected_data

def dist_function(coordinates):
    """
    Makes distance matrices
    :param coordinates:
    :return:
    """
    distances = ((coordinates[:, :, None] - coordinates[:, None, :]) ** 2).sum(axis=3)
    rows, columns = np.triu_indices_from(distances[0], 1)
    distances = distances[:, rows, columns].reshape((len(distances), -1))
    #distances -= distances.mean(axis=0)
    return distances

def parse_to_json(projected_data, peptide_information, scores):
    """
    Parses all the data to JSON for exporting
    :param projected_data:
    :param peptide_information:
    :param scores:
    :return:
    """
    data = {
        idx: {
            "peptide": val[0],
            "atomnos": {"min": val[1][0], "max": val[1][1]},
            "chain": val[2],
            "x": projected_data[:, 0][idx],
            "y": projected_data[:, 1][idx],
            "z": projected_data[:, 2][idx],
            "structure": val[3]
        }
        for idx, val in enumerate(peptide_information)
    }
    data["scores"] = {
        "x": scores[:, 0],
        "y": scores[:, 1],
        "z": scores[:, 2]
    }
    print(json.dumps(data, cls=NumpyEncoder))

def main(args):
    """
    Main function that calls everything
    """
    filename = args[1]
    oligo_length = int(args[2])
    pdb_code = args[3]

    try:
        peptide_information, pepcoords = peptidize(filename, oligo_length, pdb_code)
        vectors = np.genfromtxt(
            f'{Path(__file__).parent.parent}/vectors/vectors_{oligo_length}.csv',
            delimiter=',', comments='#'
        )[1:]
        scores = np.genfromtxt(
            f'{Path(__file__).parent.parent}/scores/scores_downsampled_{oligo_length}.csv',
            delimiter=',', comments='#'
        )[1:]
        projected_data = princana(pepcoords, vectors)
        parse_to_json(projected_data, peptide_information, scores)
    except Exception as exception:
        print(json.dumps({"error": "Untracked error: " + str(exception).replace("`", "\"")}))


if __name__ == "__main__":
    sys.exit(main(sys.argv))

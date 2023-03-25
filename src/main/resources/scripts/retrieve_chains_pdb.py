#!/usr/bin/env python3

import MDAnalysis as mda
import numpy as np
import json
import requests
import tempfile
import os
import sys

def get_chains(pdb_id):
    fd, path = tempfile.mkstemp(suffix=".pdb")
    try:
        with os.fdopen(fd, 'w') as tmp:
            url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
            with requests.get(url, stream=True) as r:
                r.raise_for_status()
                for line in r.iter_lines():
                    if line:
                        if line.decode('utf-8').startswith("ATOM"):
                            tmp.write(f"{line.decode('utf-8')}\n")
        u = mda.Universe(path)
        selection = 'protein'
        proteins = u.select_atoms(selection)

        prev = count = {protein:0 for protein in np.unique(proteins.segids)}


        for protein in proteins:
            chain = protein.segid
            if prev[chain] != protein.resid:
                prev[chain] = protein.resid
                count[chain] += 1
        print(json.dumps({p:str(count[p]) for p in count}))
    finally:
        os.close(fd)
        os.remove(path)

def main(args):
    get_chains(args[1])


if __name__ == "__main__":
    sys.exit(main(sys.argv))

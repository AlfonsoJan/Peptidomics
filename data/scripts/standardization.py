import glob
import numpy as np
import sys
import csv
import re

def get_files():
    files_list = glob.glob("pdb\*.pdb")
    np.random.seed(14)
    np.random.shuffle(files_list)
    return files_list

def get_pdb_codes(file_list):
    pdb_codes = np.array(list(map(lambda v: re.sub(r'pdb\\','', v) , np.array(file_list))))
    pdb_codes = np.array(list(map(lambda v: re.sub(r'.pdb','', v) , pdb_codes)))
    pdb_codes = np.array_split(pdb_codes, round(len(pdb_codes) / 10))
    return pdb_codes

def read_pdb_backbone(pdbfile):
    with open(pdbfile) as pdb:
        backbone = [
            line.strip() for line in pdb
            if line.startswith('ATOM') and line[12:16] in (' N  ', ' CA ', ' C  ')
        ]
    coordinates = np.loadtxt([ atom[30:54] for atom in backbone ])
    return coordinates


def dist_function(coordinates):
    distances = ((coordinates[:, :, None] - coordinates[:, None, :]) ** 2).sum(axis=3)
    rows, columns = np.triu_indices_from(distances[0], 1)
    distances = distances[:, rows, columns].reshape((len(distances), -1))
    distances -= distances.mean(axis=0)
    return distances

def peptidize(pdbfile, pepsize):
    # Because each residue has three atoms
    pepsize *= 3
    # The coordinates
    coordinates = read_pdb_backbone(pdbfile)
    # Separate parts of chains
    breaks = np.where(((coordinates[1:] - coordinates[:-1]) ** 2).sum(axis=1) > 4)[0] + 1
    parts = np.split(coordinates, breaks)
    # Split the coordinates
    pepcoords = np.array([
        part[i:i + pepsize]
        for part in parts
        for i in range(len(part) - pepsize)
    ])
    return pepcoords

def peptidize_all(pdb_list, pepsize):
    coordinates = []
    for pdb in pdb_list:
        crd = peptidize(pdb, pepsize)
        coordinates.extend(crd)
    coordinates = np.array(coordinates)
    coordinates -= coordinates.mean(axis=1)[:, None, :]
    return coordinates

def perform_pca(distances, pdb_codes, pepsize, total_files):
    _, vectors = np.linalg.eigh(np.cov(distances.T))
    vectors = vectors[:, ::-1][:, :3]
    with open(f"vectors/vectors_{pepsize}.csv", "w", newline='') as my_csv:
        my_csv.write(f"# eigen vectors based on {distances.shape[0]} coordinates\n")
        my_csv.write(f"# eigen vectors based on {total_files} pdb files\n")
        for sub_list in pdb_codes:
            my_csv.write(f"# {', '.join(sub_list)}\n")
        csvWriter = csv.writer(my_csv, delimiter=',')
        csvWriter.writerow(["X", "Y", "Z"])
        csvWriter.writerows(vectors)
    return vectors

def scores(distances, vectors, pepsize, pdb_codes, total_files):
    scores = distances @ vectors
    with open(f"scores/scores_{pepsize}.csv", "w", newline='') as my_csv:
        my_csv.write(f"# scores based on {distances.shape[0]} coordinates\n")
        my_csv.write(f"# scores based on {total_files} pdb files\n")
        for sub_list in pdb_codes:
            my_csv.write(f"# {', '.join(sub_list)}\n")
        csvWriter = csv.writer(my_csv, delimiter=',')
        csvWriter.writerow(["X", "Y", "Z"])
        csvWriter.writerows(scores)

def main(args):
    total_files = 15
    file_list = get_files()[:total_files]
    pdb_codes = get_pdb_codes(file_list)
    for pepsize in range(30, 0, -1):
        print(f"{pepsize=}")
        coordinates = peptidize_all(file_list, pepsize)
        distances = dist_function(coordinates)
        vectors = perform_pca(distances, pdb_codes, pepsize, total_files)
        scores(distances, vectors, pepsize, pdb_codes, total_files)

if __name__ == "__main__":
    sys.exit(main(sys.argv))

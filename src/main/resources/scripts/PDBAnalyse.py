#!/usr/bin/env python3

import MDAnalysis as mda
import numpy as np
import json
import sys

__author__ = "Wouter Zeevat"

class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)


def getAtomPDB(path, param):
    u = mda.Universe(path)
    selection = 'protein and (name N or name CA or name C or name O)'
    protein = u.select_atoms(selection)
    coordinates = protein.positions
    parts = np.split(coordinates,
                     4 * np.where(((coordinates[4::4] - coordinates[2:-3:4]) ** 2).sum(axis=1) > 3)[0] + 4)
    n = int(param) * 4
    P = np.array([p[i:i + n] for p in parts for i in range(0, len(p) - n, 4)])
    P -= P.mean(axis=1)[:, None, :]
    return P


def dimPlot(P):
    D = ((P[:, :, None, :] - P[:, None, :, :]) ** 2).sum(axis=3).reshape((len(P), -1)) ** 0.5
    D -= D.mean(axis=0)
    vals, vecs = np.linalg.eigh(D.T @ (D / len(D)))
    scores = np.arange(len(vals)) + 1, 100 * vals[::-1] / sum(vals)
    return scores


def scatter(P):
    D = ((P[:, :, None, :] - P[:, None, :, :]) ** 2).sum(axis=3).reshape((len(P), -1)) ** 0.5
    D -= D.mean(axis=0)
    vals, vecs = np.linalg.eigh(D.T @ (D / len(D)))
    scores = D @ vecs[:, [-1, -2, -3]]
    return scores.T


def main(args):

    P = getAtomPDB(args[1], args[2])
    P_compare = getAtomPDB(args[3], args[2])
    dim_plot = dimPlot(P)
    scores_scatter = scatter(P)
    scores_scatter_compare = scatter(P_compare)
    print(json.dumps(
        {
            'dim': {
                'x': dim_plot[0],
                'y': dim_plot[1],
            },
            'scatter': {
                'x': scores_scatter[0],
                'y': scores_scatter[1],
                'z': scores_scatter[2],
            }, "compare": {
            'x': scores_scatter_compare[0],
            'y': scores_scatter_compare[1],
            'z': scores_scatter_compare[2],
        }}
        , cls=NumpyEncoder)
    )


if __name__ == "__main__":
    sys.exit(main(sys.argv))
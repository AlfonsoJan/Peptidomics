#!/usr/bin/env python3

import numpy as np
import sys
import json

class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)


def get_plot(P):
    D = ((P[:, :, None, :] - P[:, None, :, :]) ** 2).sum(axis=3).reshape((len(P), -1)) ** 0.5
    D -= D.mean(axis=0)
    vals, vecs = np.linalg.eigh(D.T @ (D / len(D)))
    scores = D @ vecs[:, [-1, -2, -3]]
    print(json.dumps({'x': scores.T[0], 'y': scores.T[1], 'z': scores.T[2]}, cls=NumpyEncoder))

def main(args):
    p = np.load(args[1])
    get_plot(p)

if __name__ == "__main__":
    sys.exit(main(sys.argv))
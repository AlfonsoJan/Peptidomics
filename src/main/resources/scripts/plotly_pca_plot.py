#!/usr/bin/env python3

import numpy as np
import sys
import json


class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)


def get_plot(p):
    d = ((p[:, :, None, :] - p[:, None, :, :]) ** 2).sum(axis=3).reshape((len(p), -1)) ** 0.5
    d -= d.mean(axis=0)
    vals, vecs = np.linalg.eigh(d.T @ (d / len(d)))
    scores = d @ vecs[:, [-1, -2, -3]]
    print(json.dumps({'x': scores.T[0], 'y': scores.T[1], 'z': scores.T[2]}, cls=NumpyEncoder))


def main(args):
    p = np.load(args[1])
    get_plot(p)


if __name__ == "__main__":
    sys.exit(main(sys.argv))

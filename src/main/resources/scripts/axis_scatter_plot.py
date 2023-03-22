#!/usr/bin/env python3

"""
Creates the first plot from the numpy matrix
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
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

    px = 1/plt.rcParams['figure.dpi']

    fig, ax = plt.subplots(figsize=(1200*px, 500*px))
    plt.box(False)
    ax.grid(color='#AAAAAA', linestyle=(5, (10, 3)), axis='y')
    ax.set_axisbelow(True)

    ax.tick_params(axis='x', colors='#999999')
    ax.tick_params(axis='y', colors='#999999')


    scores = D @ vecs[:, [-1, -2, -3]]
    scores = scores[:, :2].T
    print(json.dumps({'x': scores[0], 'y': scores[1]}, cls=NumpyEncoder))

def main(args):
    p = np.load(args[1])
    get_plot(p)

if __name__ == "__main__":
    sys.exit(main(sys.argv))

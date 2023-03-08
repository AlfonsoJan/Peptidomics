#!/usr/bin/env python3

"""
Creates the first plot from the numpy matrix
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
from io import BytesIO
import base64

def get_plot(P):
    D = ((P[:, :, None, :] - P[:, None, :, :]) ** 2).sum(axis=3).reshape((len(P), -1)) ** 0.5
    D -= D.mean(axis=0)

    vals, vecs = np.linalg.eigh(D.T @ (D / len(D)))

    fig = plt.figure(figsize= (5,5))

    plt.scatter(np.arange(len(vals)) + 1, 100 * vals[::-1] / sum(vals))
    plt.xlim(0, 10)
    fig_file = BytesIO()
    fig.savefig(fig_file, format="png")
    fig_file.seek(0)
    base_figure = base64.b64encode(fig_file.getvalue()).decode("ascii")
    return base_figure

def main(args):
    p = np.load('./output.npy')
    get_plot(p)

if __name__ == "__main__":
    sys.exit(main(sys.argv))
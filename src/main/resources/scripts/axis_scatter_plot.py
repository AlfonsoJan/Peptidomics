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

    fig, ax = plt.subplots(figsize=(12,5))
    plt.box(False)
    ax.grid(color='#AAAAAA', linestyle=(5, (10, 3)), axis='y')
    ax.set_axisbelow(True)

    ax.tick_params(axis='x', colors='#999999')
    ax.tick_params(axis='y', colors='#999999')


    scores = D @ vecs[:, [-1, -2, -3]]
    plt.scatter(*scores[:, :2].T, s=50, alpha=0.3)
    #plt.axes().set_aspect('equal')
    fig_file = BytesIO()
    fig.savefig(fig_file, format="png", dpi=75)
    fig_file.seek(0)
    base_figure = base64.b64encode(fig_file.getvalue()).decode("ascii")
    return base_figure

def main(args):
    p = np.load(args[1])
    print(get_plot(p))

if __name__ == "__main__":
    sys.exit(main(sys.argv))

#!/usr/bin/env python3

"""Module that downsamples plot data!"""

import csv
import numpy as np
import pandas as pd


__author__ = "Jan Alfonso Busker"

NUM_CLUSTERS = 2000

def main():
    """
    Main function
    Downsamples the data to num_clusters
    :return:
    """
    for pepsize in range(30, 0, -1):
        print(f"{pepsize=}")
        scores = pd.read_csv(f'./scores/scores_{pepsize}.csv', sep=",",
                             header=0, comment='#').to_numpy()

        # Get the number of points in the original array
        num_points_original = scores.shape[0]

        # Randomly select indices of the downsampled points
        downsampled_indices = np.random.choice(num_points_original, size=NUM_CLUSTERS, replace=False)

        # Select the corresponding coordinates based on the indices
        downsampled_coordinates = scores[downsampled_indices]

        with open(f"scores/scores_downsampled_{pepsize}.csv", "w", newline='') as my_csv:
            my_csv.write(f"# Down sampled to {NUM_CLUSTERS} coordinates\n")
            csv_writer = csv.writer(my_csv, delimiter=',')
            csv_writer.writerow(["X", "Y", "Z"])
            csv_writer.writerows(downsampled_coordinates)

if __name__ == "__main__":
    main()

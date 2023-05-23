#!/usr/bin/env python3

"""Module that downsamples plot data!"""

import csv
import numpy as np
import pandas as pd
from sklearn.cluster import KMeans


__author__ = "Jan Alfonso Busker"

NUM_CLUSTERS = 1000

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

        # Apply K-means clustering
        kmeans = KMeans(n_clusters=NUM_CLUSTERS)
        kmeans.fit(scores)

        # Get the cluster labels for each data point
        cluster_labels = kmeans.labels_

        # # Initialize an array to store the downsized coordinates
        downsampled_coordinates = np.zeros((NUM_CLUSTERS, 3))

        for cluster in range(NUM_CLUSTERS):
            cluster_indices = np.where(cluster_labels == cluster)[0]
            cluster_coordinates = scores[cluster_indices]
            representative_coordinate = np.mean(cluster_coordinates, axis=0)
            downsampled_coordinates[cluster] = representative_coordinate

        with open(f"scores/scores_downsampled_{pepsize}.csv", "w", newline='') as my_csv:
            my_csv.write(f"# Down sampled to {NUM_CLUSTERS} coordinates\n")
            csv_writer = csv.writer(my_csv, delimiter=',')
            csv_writer.writerow(["X", "Y", "Z"])
            csv_writer.writerows(downsampled_coordinates)

if __name__ == "__main__":
    main()

import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
import csv

num_clusters = 1000
for pepsize in range(30, 0, -1):
    print(f"{pepsize=}")
    scores = pd.read_csv(f'./scores/scores_{pepsize}.csv', sep=",", header=0, comment='#').to_numpy()

    # Apply K-means clustering
    kmeans = KMeans(n_clusters=num_clusters)
    kmeans.fit(scores)

    # Get the cluster labels for each data point
    cluster_labels = kmeans.labels_

    # # Initialize an array to store the downsized coordinates
    downsampled_coordinates = np.zeros((num_clusters, 3))

    for cluster in range(num_clusters):
        cluster_indices = np.where(cluster_labels == cluster)[0]
        cluster_coordinates = scores[cluster_indices]
        representative_coordinate = np.mean(cluster_coordinates, axis=0)
        downsampled_coordinates[cluster] = representative_coordinate

    with open(f"scores/scores_downsampled_{pepsize}.csv", "w", newline='') as my_csv:
        my_csv.write(f"# Down sampled to {num_clusters} coordinates\n")
        csvWriter = csv.writer(my_csv, delimiter=',')
        csvWriter.writerow(["X", "Y", "Z"])
        csvWriter.writerows(downsampled_coordinates)

####################################################################
# Given a csv file of kmer frequencies for reads within each barcode,
# create a UMAP that will cluster the reads based on the similarity to
# each other. UMAP is a general purpose manifold learning and dimension
# reduction algorithm. 
# Output: A dataframe of reads coordinates and a dataframe of reads and 
# the groups they were clustered into.
# (Following https://umap-learn.readthedocs.io/en/latest/basic_usage.html)
#####################################################################

import argparse
import hdbscan
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import adjusted_rand_score, adjusted_mutual_info_score
import sys
import umap.umap_ as umap


parser = argparse.ArgumentParser()

parser.add_argument("-kmer", "--kmer_freq_matrix", help="path to kmer frequency matrix file")
parser.add_argument("-out", "--output_path", help="path for the output file")
parser.add_argument("-res_out", "--results_out", help="path to the results output file")
parser.add_argument("-nc", "--n_components", help="number of components used for UMAP clustering")
parser.add_argument("-nn", "--umap_n_neighbors", help="number of neighbours for UMAP clustering (local versus global structure)")
parser.add_argument("-min_cs", "--hdbscan_min_cluster_size", help="min number of reads that will create a cluster")
parser.add_argument("-csm", "--hdbscan_cluster_selection_method", help="method used to select clusters in hdbscan")

args = parser.parse_args()

# Read in the kmer frequency file
file_name = args.kmer_freq_matrix
out_path = args.output_path
name = file_name.split('/')[-1].split('.')[0]

df = pd.read_csv(file_name, index_col=0)
species = df.index.tolist()

sys.stdout.write('\tUMAP: Construct a UMAP object. \n')
# Construct a UMAP object
clusterable_reducer = umap.UMAP(
    n_neighbors=int(args.umap_n_neighbors),
    min_dist=0.0,
    n_components=int(args.n_components),
    random_state=42)

standard_reducer = umap.UMAP(random_state=42)

# Convert each feature into z-scores (number of standard deviations from the mean) for comparability.
read_data = df[:].values
scaled_read_data = StandardScaler().fit_transform(read_data)

# Use the fit_transform method which first calls fit and then returns the transformed data as a numpy array.
# UMAP reduces down to a 2D matrix. 
# Each row of the array is a 2-dimensional representation of the corresponding barcode read.
standard_embedding = standard_reducer.fit_transform(scaled_read_data)

clusterable_embedding = clusterable_reducer.fit_transform(scaled_read_data)

sys.stdout.write('\tUMAP: cluster clusters using HDBSCAN. \n')
labels = hdbscan.HDBSCAN(
    min_cluster_size=int(args.hdbscan_min_cluster_size),
    cluster_selection_epsilon=0.05,
    cluster_selection_method=args.hdbscan_cluster_selection_method
).fit_predict(clusterable_embedding)

clustered = (labels >= 0)
plt.figure(figsize=(14, 8))
plt.rcParams['axes.facecolor'] = 'black'

# Plot the points that weren't clustered (as gray dots)
plt.scatter(standard_embedding[~clustered, 0],
            standard_embedding[~clustered, 1],
            facecolors='none', edgecolors='gray',
            s=5.0)

# Plot the points that were clustered (as colored dots)
plt.scatter(standard_embedding[clustered, 0],
            standard_embedding[clustered, 1],
            c=labels[clustered],
            s=5.0,
            cmap='plasma')

plt.title('UMAP projection of \n' + name + ' \nreads kmer frequency matrix and hdbscan clustering \n(' 
            + str(max(labels)+1) + ' clusters found)')

# Label each cluster on the UMAP projection.
# Calculate the middle of the cluster and print the label there.
label_coordinates = {}
label_middle = {}
# First gather all coordinates based on the label
for i in range(0, len(labels)):
    if labels[i] not in label_coordinates.keys():
        label_coordinates[labels[i]] = [(standard_embedding[i][0], standard_embedding[i][1])]
    else:
        label_coordinates[labels[i]].append((standard_embedding[i][0], standard_embedding[i][1]))

# Next, calculate the mean for each coordinate.
for k in label_coordinates.keys():
    x_total, y_total = 0,0
    for j in label_coordinates[k]:
        x_total += j[0]
        y_total += j[1]
    x_mean = x_total / float(len(label_coordinates[k]))
    y_mean = y_total / float(len(label_coordinates[k]))
    label_middle[str(k)] = [x_mean, y_mean]

for i in range (0, max(labels)+1):
    plt.text(label_middle[str(i)][0], label_middle[str(i)][1], str(i), color='white')
    

# Make a table of reads and which cluster they were clustered into
sys.stdout.write('\tUMAP: Write output files. \n')
out_file = out_path + '/' + name + '_clustering.csv'

plt.savefig( args.results_out + '/' + name + '_umap.png')

with open(out_file, 'w') as out:
    for i in range(0, len(labels)):
        out.write(species[i] + ',' + str(labels[i]) + ',' +  str(standard_embedding[i][0]) + ',' + str(standard_embedding[i][1]) + '\n') 

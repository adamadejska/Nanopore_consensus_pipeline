#####################################################################
# After the initial cllustering using UMAP (umap_clustering.py), we 
# need to refine some of the clusters to rescue some reads that 
# accidentally might have grouped with some other species.
# To divide them, make an alignment on full reads and create a tree
# based on the alignment. if there are any deep branches, we can divide
# the subtrees on those deep branches closer to the root.
# Output: a file containing a list of reads and the new, refined clusters.
#####################################################################

import argparse
from Bio import SeqIO
from Bio.Cluster import kcluster
from ete3 import PhyloTree
import numpy as np
import os




output_dir = '/home/ada/Desktop/16S_alignments/scripts/16S_alignment/reads_clustering_pipeline/main_pipeline/src/tmp/mock6_1000_varying_tmp'
job_name = 'mock6_1000_varyingaa'
fasta_file = '/home/ada/Desktop/16S_alignments/scripts/16S_alignment/reads_clustering_pipeline/main_pipeline/src/tmp/mock6_1000_varying_tmp/fastq_files/mock6_1000_varyingaa'
out_file = output_dir + '/'+ job_name + '_refined_clusters.txt'

final_out = open(out_file, 'w')

# Make a folder for storing reads of each cluster as separate fasta files
cluster_fasta_dir = output_dir + '/cluster_fasta/'

if not os.path.isdir(cluster_fasta_dir):
    os.mkdir(cluster_fasta_dir)

os.system('rm  ' + cluster_fasta_dir + '/*')

# Read in the fastq file and save which read has what sequence
name_to_seq = {}
for record in SeqIO.parse(fasta_file, "fasta"):
    name = str(record.id).split(' ')[0].strip()
    seq = str(record.seq)
    name_to_seq[name] = seq

cluster_file = '/home/ada/Desktop/16S_alignments/scripts/16S_alignment/reads_clustering_pipeline/main_pipeline/src/tmp/mock6_1000_varying_tmp/mock6_1000_varyingaa_clustering.csv'
# For each cluster, make a separate fasta file for each read.

clusters = []
# Check the number of clusters in the file
with open(cluster_file, 'r') as f:
    for line in f:
        c = line.split(',')[1]
        clusters.append(int(c))


for c in range(4, 5):
    print('\t Refinement: processing cluster ' + str(c))

    read_counter = 0
    guide_tree_file = output_dir + '/' + job_name + '_guide_tree_cluster_' + str(c) + '.dnd'
    read_names = []

    with open(cluster_file, 'r') as f:
        for line in f:
            line = line.split(',')
            read_name = line[0].split('/')[-1].strip()
            cluster = int(line[1])
            if cluster == c:
                read_names.append(read_name)
                tmp_fasta_file = cluster_fasta_dir + '/' + read_name + '.fasta'
                with open(tmp_fasta_file, 'w') as out:
                    out.write('>' + read_name + '\n')
                    out.write(name_to_seq[read_name] + '\n')
                    read_counter += 1

    print('\t Refinement: creating a tree for ' + str(read_counter) + ' reads.')
    os.system('mashtree '+ cluster_fasta_dir + '/*   > ' + guide_tree_file)

    guide_tree_file = '/home/ada/Desktop/16S_alignments/scripts/16S_alignment/reads_clustering_pipeline/main_pipeline/src/tmp/mock6_1000_varying_tmp/cluster_fasta/guide_tree_test.txt'

    t = PhyloTree(guide_tree_file)
    std = np.std(np.array([t.get_distance(i) for i in t.get_leaf_names()]))
    print(std)


    if std > 0.01:
        print('\t Refinement: cluster ' + str(c) + ' passed the refinement checkpoint. Proceeding with splitting.')
        
        leafs = list(t.get_leaf_names())
        data = np.zeros((len(leafs), len(leafs)))
        for row in range(0, len(leafs)):
            for col in range(0,len(leafs)):
                if row == col:
                    data[row][col] = 0
                else:
                    data[row][col] = t.get_distance(leafs[row], leafs[col])
                    data[col][row] = t.get_distance(leafs[row], leafs[col])

        std = np.std(data)
        print(std)
        print(len(read_names))
        # k means clustering, determine the number of clusters (k) by the number of reads we work with
        k = int(read_counter / 250)  # roughly, each cluster should have 250 reads. More than enough for consensus.
        clusterid, error, nfound = kcluster(data, nclusters=k)
        clusterid = np.array(list(clusterid))

        print(len(clusterid))

        read_names = np.array(read_names)
        # divide the read names into each group
        subgroups = {}
        for k in range(0, max(clusterid)+1):
            subgroups[k] = read_names[clusterid==k]

        for k,v in subgroups.items():
            for read in v:
                read = read.split('.')[0]
                cluster_num = str(c) + '_' + str(k)
                final_out.write(read + ',' + cluster_num + '\n')

    else:
        # If the std of branch length is not enough, just copy all the reads of the cluster. We won't
        # try to split it further.
        for read in read_names:
            final_out.write(read + ',' + str(c) + '\n')
   
   
    os.system('rm  ' + cluster_fasta_dir + '/*')

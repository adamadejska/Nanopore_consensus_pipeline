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


output_dir = '/home/ada/Desktop/16S_alignments/scripts/16S_alignment/reads_clustering_pipeline/main_pipeline/src/tmp/mock5_1000_varying_tmp'
job_name = 'mock5_1000_varyingaa'
fasta_file = '/home/ada/Desktop/16S_alignments/scripts/16S_alignment/reads_clustering_pipeline/main_pipeline/src/tmp/mock5_1000_varying_tmp/fastq_files/mock5_1000_varyingaa'
out_file = output_dir + '/'+ job_name + '_refined_clusters.txt'
dependencies_dir = '/home/ada/Desktop/16S_alignments/scripts/16S_alignment/reads_clustering_pipeline/main_pipeline/src/dependencies'
current_dir = os.getcwd()

final_out = open(out_file, 'w')

# Read in the fastq file and save which read has what sequence
name_to_seq = {}
for record in SeqIO.parse(fasta_file, "fasta"):
    name = str(record.id).split(' ')[0].strip()
    seq = str(record.seq)
    name_to_seq[name] = seq

cluster_file = '/home/ada/Desktop/16S_alignments/scripts/16S_alignment/reads_clustering_pipeline/main_pipeline/src/tmp/mock5_1000_varying_tmp/mock5_1000_varyingaa_clustering.csv'
# For each cluster, make a separate fasta file for each read.

clusters = []
# Check the number of clusters in the file
with open(cluster_file, 'r') as f:
    for line in f:
        c = line.split(',')[1]
        clusters.append(int(c))



for c in range(0, max(clusters)+1):
    print('\t Refinement: processing cluster ' + str(c))

    read_counter = 0
    guide_tree_file = output_dir + '/' + job_name + '_guide_tree_cluster_' + str(c) + '.dnd'
    read_names = []
    indiv_cluster_file = output_dir + '/' + job_name + '_tmp_ind_cluster_' + str(c) + '.fasta'

    with open(indiv_cluster_file, 'w') as out:
        with open(cluster_file, 'r') as f:
            for line in f:
                line = line.split(',')
                read_name = line[0].split('/')[-1].strip()
                cluster = int(line[1])
                if cluster == c:
                    out.write('>' + read_name + '\n')
                    out.write(name_to_seq[read_name] + '\n')
                    read_counter += 1
   
    print('\t Refinement: creating a tree for ' + str(read_counter) + ' reads.')

    os.system(dependencies_dir + '/mbed/./mBed -infile ' + indiv_cluster_file + ' -method SparseMap -useSeedsOnly true ')

    #t = PhyloTree(guide_tree_file)
    #std = np.std(np.array([t.get_distance(i) for i in t.get_leaf_names()]))
    #print(std)
        
    dist_file = current_dir + '/distMat.out'

    reads_counter = 0 
    with open(dist_file, 'r') as f:
        for line in f:
            reads_counter += 1

    print('-------')
    print(reads_counter)
    data = np.zeros((reads_counter, reads_counter))
    read_names = []
    row_n = 0
    with open(dist_file, 'r') as f:
        for line in f:
            read_names.append(line.split()[0])
            row = line.split()[1:]
            col_n = 0
            for col in row:
                if row_n == col_n:
                    data[row_n][col_n] = 0
                else:
                    data[row_n][col_n] = col
                    data[col_n][row_n] = col
                col_n += 1
            row_n += 1

    std = np.std(data)
    final_out.write(str(std))
    
    if std > 0.02:
        print('\t Refinement: cluster ' + str(c) + ' passed the refinement checkpoint. Proceeding with splitting.')
        

        # k means clustering, determine the number of clusters (k) by the number of reads we work with
        k = int(reads_counter / 150)  # roughly, each cluster should have 250 reads. More than enough for consensus.
        print(k)
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
   
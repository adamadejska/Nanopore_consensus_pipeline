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
from ete3 import PhyloTree
import numpy as np
import os


def recursive_subclustering(root, t, parent):
    """
    This function divides MSA guide trees into subclades based on the differences in branch lengths.
    The function is called recursively, traversing the tree from bottom to the top and if it detects
    a change in branch length that is above a set threshold, it will call the leafs attached to that 
    point as two separate subclusters. 
    Returns a list of nodes that indicate the start of a new subcluster.
    """
    if root.get_children() == []:
        return([], t.get_distance(root, parent))
    
    c1, c2 = root.get_children()

    val0, len0 = recursive_subclustering(c1, t, root)
    val1, len1 = recursive_subclustering(c2, t, root)

    if np.std(np.array([len0, len1])) > 0.0001:
        ret0, ret1 = [c1], [c2]
    else:
        ret0, ret1 = [], []

    ret = []

    if val0 == []:
        ret= ret0 + ret
    else: 
        ret = val0 + ret

    if val1 == []:
        ret = ret + ret1
    else:
        ret = ret + val1
    
    return(ret, max(len0, len1) + t.get_distance(root, parent))


parser = argparse.ArgumentParser()

parser.add_argument("-fasta", "--fasta_file", help="path to FASTQ file")
parser.add_argument("-clusters", "--cluster_file", help="path to clusters file")
parser.add_argument("-out", "--output_directory", help="path to the output directory")
parser.add_argument("-name", "--job_name", help="the name of the job / what the outfiles will use as part of their name")

args = parser.parse_args()

# Read in the fastq and clusters files
fasta_file = args.fasta_file
cluster_file = args.cluster_file

name = cluster_file.split('/')[-1].replace('_kmer_matrix_clustering.csv', '')
out_file = args.output_directory + '/'+ args.job_name + '_refined_clusters.txt'

# Read in the fastq file and save which read has what sequence
name_to_seq = {}
for record in SeqIO.parse(fasta_file, "fasta"):
    name = str(record.id).split(' ')[0].strip()
    seq = str(record.seq)
    name_to_seq[name] = seq

clusters = []
# Check the number of clusters in the file
with open(cluster_file, 'r') as f:
    for line in f:
        c = line.split(',')[1]
        clusters.append(int(c))


# Read the cluster file. For each cluster, use clustal omega to create alignment and guide tree
final_out = open(out_file, 'w')
indiv_cluster_file = args.output_directory + '/tmp_ind_cluster.fasta'
for c in range(0, max(clusters)):
    print('\t Refinement: processing cluster ' + str(c))
    with open(indiv_cluster_file, 'w') as out:
        with open(cluster_file, 'r') as f:
            for line in f:
                line = line.split(',')
                read_name = line[0].split('/')[-1].strip()
                cluster = int(line[1])
                if cluster == c:
                    out.write('>' + read_name + '\n')
                    out.write(name_to_seq[read_name] + '\n')

    # run clustal omega on the tmp individual cluster file. Use the alignment and resulting guide tree 
    # to refine our clustering even more.
    guide_tree_file = args.output_directory + '/tmp_guide_tree.xml'
    clustal_out_file = args.output_directory + '/tmp_clustal_alignments.txt'
    os.system('/home/ada/Desktop/16S_alignments/scripts/16S_alignment/reads_clustering_pipeline/main_pipeline/src/dependencies/./clustalo -i ' 
              + indiv_cluster_file + '  --guidetree-out=' + guide_tree_file + ' --out ' + clustal_out_file + ' --force')


    # Check if the tree is worth refining. If the standard deviation of the leaf - root distances are 
    # long enough, we should probably split the tree.
    t = PhyloTree(guide_tree_file)
    std = np.std(np.array([t.get_distance(i) for i in t.get_leaf_names()]))

    if std > 0.0001:
        print('\t Refinement: cluster ' + str(c) + ' passed the refinement checkpoint. Proceeding with splitting.')
        
        leafs = list(t.get_leaf_names())
        print('\t Refinement: checking  ' + str(len(leafs)) + ' leafs.')
        
        # Run the recursive function for determining subclusters.
        nodes, dist = (recursive_subclustering(t, t, t))

        counter = 0
        for i in nodes:
            cluster = i.get_leaf_names()
            if len(cluster) > 3:
                counter += 1
                cluster_name = str(c) + '_' + str(counter)
                for j in cluster:
                    final_out.write(j + ',' + cluster_name + '\n')
    
    else:
        # If the std of branch length is not enough, just copy all the reads of the cluster. We won't
        # try to split it further.
        with open(indiv_cluster_file, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    name = line[1:].strip()
                    final_out.write(name + ',' + str(c) + '\n')

    os.system('rm ' + guide_tree_file)  # clustalO will not overwrite the guidetree file.
    os.system('rm ' + clustal_out_file)

final_out.close()
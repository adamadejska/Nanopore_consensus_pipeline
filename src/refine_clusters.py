#####################################################################
# After the initial cllustering using UMAP (umap_clustering.py), we 
# need to refine some of the clusters to rescue some reads that 
# accidentally might have grouped with some other species.
# To divide them, make an alignment on full reads and create a tree
# based on the alignment. if there are any deep branches, we can divide
# the subtrees on those deep branches closer to the root.
# Output: a file containing a list of reads and the new, refined clusters.
#####################################################################

from Bio import Phylo
from Bio import SeqIO
import argparse
import os


# Code from: https://stackoverflow.com/questions/3067529/a-set-union-find-algorithm
class DisjointSet(object):

    def __init__(self):
        self.leader = {} # maps a member to the group's leader
        self.group = {} # maps a group leader to the group (which is a set)

    def add(self, a, b):
        leadera = self.leader.get(a)
        leaderb = self.leader.get(b)
        if leadera is not None:
            if leaderb is not None:
                if leadera == leaderb: return # nothing to do
                groupa = self.group[leadera]
                groupb = self.group[leaderb]
                if len(groupa) < len(groupb):
                    a, leadera, groupa, b, leaderb, groupb = b, leaderb, groupb, a, leadera, groupa
                groupa |= groupb
                del self.group[leaderb]
                for k in groupb:
                    self.leader[k] = leadera
            else:
                self.group[leadera].add(b)
                self.leader[b] = leadera
        else:
            if leaderb is not None:
                self.group[leaderb].add(a)
                self.leader[a] = leaderb
            else:
                self.leader[a] = self.leader[b] = a
                self.group[a] = set([a, b])


parser = argparse.ArgumentParser()

parser.add_argument("-fastq", "--fastq_file", help="path to FASTQ file")
parser.add_argument("-clusters", "--cluster_file", help="path to clusters file")
parser.add_argument("-out", "--output_directory", help="path to the output directory")

args = parser.parse_args()

# Read in the fastq and clusters files
fastq_file = args.fastq_file
cluster_file = args.cluster_file

name = cluster_file.split('/')[-1].replace('_kmer_matrix_clustering.csv', '')
out_file = args.output_directory + '/'+ name + '_refined_clusters.txt'

# Read in the fastq file and save which read has what sequence
name_to_seq = {}
for record in SeqIO.parse(fastq_file, "fasta"):
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
for c in range(-1, max(clusters)):
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
              + indiv_cluster_file + '  --guidetree-out=' + guide_tree_file + ' --out ' + clustal_out_file)

    # Parse the guide tree
    tree = Phylo.read(guide_tree_file, 'newick')

    # Refine the clusters based on the guide tree of Clustal Omega. If the distance between 2 leafs is
    # greater than a threshold, consider them not related.
    threshold = 0.35
    # Get the names of the leafs
    leafs = list(tree.get_terminals())
    names = [x.name for x in leafs]
    
    ds = DisjointSet()
    for i in names:
        for j in names:
            dist = tree.distance(i, j)
            if dist <= threshold:
                ds.add(i,j)

    counter = 0
    for i in ds.group.keys():
        cluster = list(ds.group[i])
        if len(cluster) > 3:
            counter += 1
            cluster_name = str(c) + '_' + str(counter)
            for j in cluster:
                final_out.write(j + ',' + cluster_name + '\n')

    os.system('rm ' + guide_tree_file)  # clustalO will not overwrite the guidetree file.
    os.system('rm ' + clustal_out_file)

final_out.close()
#####################################################################
# Given a file of read names and cluster numbers, make consensus
# sequences for all clusters. Output a file with cluster number and
# the consensus sequence.
#####################################################################

import argparse
from Bio import SeqIO
import os

# Use the function for consensus creation to create consesus on reads based on the groups created by MSA
def reverse_strand(strand, orientation):
    # If the orientation is '-' then reverse the strand.

    d = {'A':'T','T':'A','C':'G','G':'C'}
    strand = strand.upper()
    if orientation == '+':
        return(strand)
    else:
        strand = strand[::-1]
        new_strand = ''
        for l in strand:
            new_strand += d[l] 
        return(new_strand)

        
def make_consensus_seq(bacteria_list, full_sequences, out_name, tmp_path, dep_path):
    # This function creates a consensus sequence for the list of bacteria provided
    # and the dictionary of full sequences from the 16S barcode files. 

    # First, make sure all sequences are in the same orientation.
    # Use minimap2 to put the sequences in the right orientation then read the output and reorient the sequences as needed.
    target_file = tmp_path + '/target_tax_' + out_name + '.fasta'
    with open(target_file, 'w') as f:
        bacteria = bacteria_list[0]
        f.write('>' + bacteria + '\n')     # Arbitrarily use the first sequence. This doesn't matter.
        f.write(full_sequences[bacteria]+'\n')

    fasta_file = tmp_path + '/cluster_seq_tax_' + out_name + '.fasta'
    # Make a file to store the sequences and names that need to be aligned.
    with open(fasta_file, 'w') as f:
        for bacteria in bacteria_list:
            f.write('>' + bacteria + '\n')
            f.write(full_sequences[bacteria]+'\n')

    minimap_out = tmp_path + '/minimap_out_tax_' + out_name + '.txt'
    os.system(dep_path + '/minimap2/./minimap2 ' + target_file + ' ' +  fasta_file + ' > ' +  minimap_out + ' 2> /dev/null')

    # Reorient the sequences based on minimap2 output before further processing.
    corrected_fasta = tmp_path + '/cluster_seq_tax_corrected_' + out_name + '.fasta'
    fasta = open(corrected_fasta,'w')
    with open(minimap_out, 'r') as f:
        for line in f:
            line = line.split('\t')
            name, orientation = line[0], line[4]
            strand = reverse_strand(full_sequences[name], orientation)
            fasta.write('>' + name + '\n')
            fasta.write(strand+'\n')
    fasta.close()

    # Now that the sequences are corrected, align them.

    msa_file = tmp_path + '/cluster_aligned_tax_' + out_name + '.msf'
    # Use MAFFT to align the sequences and create a temporary MSA file
    os.system('mafft ' + corrected_fasta + ' > ' + msa_file + ' 2> /dev/null')

    consensus_file = tmp_path + '/cluster_consensus_tax_' + out_name + '.fasta'
    # Use em_cons for creating a consensus sequence from the multiple alignment. 
    os.system('em_cons ' + msa_file + ' ' + consensus_file + ' 2> /dev/null')

    # Create a consensus without the 'N's
    no_n_consensus = ''
    with open(consensus_file, 'r') as f:
        f.readline()    # header
        for line in f:
            no_n_consensus += line.replace("n", "").strip()
    
    return(no_n_consensus)


parser = argparse.ArgumentParser()

parser.add_argument("-fasta", "--fasta_file", help="path to FASTA file")
parser.add_argument("-clusters", "--cluster_file", help="path to clusters file")
parser.add_argument('-tmp_out', "--tmp_output_path", help="path to the tmp folder")
parser.add_argument("-out", "--output_path", help="path to the output folder")
parser.add_argument("-dep", "--dependencies_path", 
                    help="path to the dependenceis folder (folder with blast, minimap, etc.)")
parser.add_argument("-name", "--job_name", help="the name of the job / what the outfiles will use as part of their name")
                   
args = parser.parse_args()

fasta_file = args.fasta_file
cluster_file = args.cluster_file
tmp_path = args.tmp_output_path



# Make a dictionary of full sequences from the barcode's fastq file. 
full_sequences = {}
for record in SeqIO.parse(fasta_file, "fasta"):
    name = str(record.id)
    full_sequences[name] = str(record.seq)


cluster_to_reads = {}
with open(cluster_file, 'r') as f:
    for line in f:
        line = line.split(',')
        id, c = line[0].split('/')[-1], line[1]
        c = c.strip()
        
        if c not in cluster_to_reads.keys():
            cluster_to_reads[c] = [id]
        else:
            cluster_to_reads[c].append(id)

# Make a consensus sequence for each cluster.
out = args.output_path + '/' + args.job_name + '_consensus.txt'

with open(out, 'w') as f:
    for k, v in cluster_to_reads.items():
        #if k == '-1':   # Skip the unclassified reads
        #    continue

        consensus = make_consensus_seq(v, full_sequences, args.job_name, tmp_path, args.dependencies_path)
        if len(consensus) > 0:
            f.write('>cluster ' + k + '\n')
            f.write(consensus + '\n')


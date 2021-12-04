####################################################################
# Given a file with a list of noise reads and a list of consensus sequences,
# combine both of those files into one fasta file.
####################################################################

import argparse
from Bio import SeqIO
import os

parser = argparse.ArgumentParser()

parser.add_argument("-fastq", "--fastq_file", help="path to FASTQ file")
parser.add_argument("-consensus", "--consensus_folder", help="path to the folder with consensus files")
parser.add_argument("-noise", "--noise_file", help="path to the noise file")
parser.add_argument("-out", "--output_directory", help="path to the output directory")
parser.add_argument("-name", "--job_name", help="the name of the job / what the outfiles will use as part of their name")

args = parser.parse_args()

out_file = args.output_directory + '/' + args.job_name + '_final.fasta'

# Find all the consensus file paths.
consensus_files = [i for i in os.listdir(args.consensus_folder) if i.endswith('_consensus.txt')]

# Write in the consensus sequences into one out file.
with open(out_file, 'w') as f:
    read_counter = 0
    for i in consensus_files:
        with open(args.consensus_folder + i, 'r') as c:
            for line in c:
                if not line.startswith('>'):
                    read_name = 'consensus_' + str(read_counter)
                    f.write('>' + read_name + '\n')
                    f.write(line)
                    read_counter += 1

# Make a list of all the names of the noise reads. We'll need to find their sequence in the OG fastq file.
noise_reads = []

with open(args.noise_file, 'r') as f:
    for line in f:
        noise_reads.append(line.split('/')[-1].strip())

# Go through the fastq file and add the noise sequences to the final out file.
with open(out_file, 'a') as f:
    for record in SeqIO.parse(args.fastq_file, "fastq"):
        name = str(record.id)
        seq = str(record.seq)
        if name in noise_reads:
            f.write('>' + name + '\n')
            f.write(seq + '\n')


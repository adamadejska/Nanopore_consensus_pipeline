############################################################
# This script joins consensus from noise reads and consensus from 
# each smaller job.
############################################################

import argparse
import os

parser = argparse.ArgumentParser()

parser.add_argument("-cons", "--consensus_directory", help="path to the consensus directory")
parser.add_argument("-out", "--output_directory", help="path to the output directory")
parser.add_argument("-name", "--job_name", help="the name of the job / what the outfiles will use as part of their name")

args = parser.parse_args()

out_file = args.output_directory + '/' + args.job_name + '_final_consensus.fasta'

# Find all the consensus file paths.
consensus_files = [i for i in os.listdir(args.consensus_directory) if i.endswith('_consensus.txt')]

# Write in the consensus sequences into one out file.
with open(out_file, 'w') as f:
    read_counter = 0
    for i in consensus_files:
        with open(args.consensus_directory + i, 'r') as c:
            for line in c:
                if not line.startswith('>'):
                    read_name = 'consensus_' + str(read_counter)
                    f.write('>' + read_name + '\n')
                    f.write(line)
                    read_counter += 1
#####################################################################
# Given a file with the consensus sequences, run it through BLASTN.
#####################################################################

import argparse
from Bio import SeqIO
import os

parser = argparse.ArgumentParser()

parser.add_argument("-consensus", "--consensus_file", help="path to the cluster consensus file")
parser.add_argument("-out", "--output_path", help="path to the output folder")
parser.add_argument("-dep", "--dependencies_path", 
                    help="path to the dependenceis folder (folder with blast, minimap, etc.)")
parser.add_argument("-name", "--job_name", help="the name of the job / what the outfiles will use as part of their name")


args = parser.parse_args()
consensus = args.consensus_file

# Using the consensus file, run each sequence in BLASTN
blast_out =  args.output_path + '/' + args.job_name + '_blastn_result.fa'
print('consensus path: ' + consensus)
print('out_folder: ' + blast_out)
blast_path = args.dependencies_path + '/ncbi-blast-2.12.0+/blastdb/'
os.system('blastn -db ' + blast_path + '16S_ribosomal_RNA -query ' + args.consensus_file + ' -out ' + blast_out + ' -max_target_seqs 5')

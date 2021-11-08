############################################################
# Before we start the kmer frequency and clustering analysis
# on the 16S sequences, we need to check the quality of the reads.
############################################################

import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser()

parser.add_argument("-fastq", "--fastq_file", help="path to FASTQ file")
parser.add_argument("-out", "-output_name", help="name for the output file")

args = parser.parse_args()

# Make sure the reads are long enough to be 16S sequences (around 1.5 k bp in length)
name = args.output_name
out_file = '/home/ada/Desktop/16S_alignments/scripts/16S_alignment/reads_clustering_pipeline/tmp_data/'+ name + 'QC.fasta'

with open(out_file, 'w') as out:
    for record in SeqIO.parse(args.fastq_file, "fastq"):
        name = str(record.id)
        seq = str(record.seq)
        if len(seq) > 1200 and len(seq) < 1700:
            out.write('>' + name + '\n')
            out.write(seq + '\n')


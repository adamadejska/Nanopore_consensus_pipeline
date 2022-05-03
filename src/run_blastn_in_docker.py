#####################################################################################
# Run BLASTN remotely. Do so in the background and run each consensus read separately
# to save time. 
#####################################################################################


import argparse
import os

# Parse all the arguments from the flags.
parser = argparse.ArgumentParser()

parser.add_argument("-consensus", "--consensus_file", help="path to the cluster consensus file")
parser.add_argument("-out", "--output_path", help="path to the output folder")
parser.add_argument("-dbpath", "--database_path", 
                    help="path to the database files folder (folder with files downloaded from NCBI)")
parser.add_argument("-db_name", "--database_name", help="name of the database used for BLASTN")


args = parser.parse_args()

consensus_file = args.consensus_file
current_dir = '/'.join(__file__.split('/')[:-1])

consensus_file_name = consensus_file.split('/')[-1]
blast_db = args.database_name

# Run BLASTN from a docker image
os.system(current_dir + '/./run_blast.sh ' + args.database_path + ' ' + args.output_path + ' ' + consensus_file_name + ' ' + blast_db)

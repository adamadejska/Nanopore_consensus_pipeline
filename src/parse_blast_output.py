#####################################################################
# Given a BLASTN ouput file, parse it and organize the IDs for each cluster
# in a more organized way. Include the metadata for each BLAST hit. 
#####################################################################

"""
Sample BLASTN output:

Query= cluster 11

Length=1567
                                                                      Score        E
Sequences producing significant alignments:                          (Bits)     Value

NR_164943.1 Wenzhouxiangella salilacus strain 15181 16S ribosomal...  2228       0.0  
NR_152707.1 Wenzhouxiangella sediminis strain XDB06 16S ribosomal...  2217       0.0  
NR_136878.1 Wenzhouxiangella marina strain Ma-11 16S ribosomal RN...  2137       0.0  
NR_074775.1 Alkalilimnicola ehrlichii strain MLHE-1 16S ribosomal...  1960       0.0  
NR_148757.1 Thiolapillus brandeum strain Hiromi 1 16S ribosomal R...  1893       0.0  


>NR_164943.1 Wenzhouxiangella salilacus strain 15181 16S ribosomal RNA, partial 
sequence
Length=1524

 Score = 2228 bits (1206),  Expect = 0.0
 Identities = 1424/1525 (93%), Gaps = 31/1525 (2%)
 Strand=Plus/Plus

Sample parsed identity file:

cluster#, ID, full_name, score, e_value, identity_%, gaps_%
cluster 11, NR_164943.1, Wenzhouxiangella salilacus strain 15181 16S ribosomal RNA, 1206, 0.0, 93, 2
...
"""

import argparse


# Add command line arguments for the in and out files.
parser = argparse.ArgumentParser()

parser.add_argument("-blast", "--blast_file", help="path to BLASTN output file")
parser.add_argument("-out", "--output_path", help="path to parsed output file")

args = parser.parse_args()

blast_file = args.blast_file
out_path = args.output_path

file_name = blast_file.split('/')[-1].split('_')[0]
out_file = out_path + '/' + file_name + '_final_cluster_identities.csv'


# Read the BLAST output file line by line and extract the info on the top hit for each cluster.
# Save the results to the appropriate output file.
with open(out_file, 'w') as out:
    out.write('cluster#, ID, full_name, score, e_value, identity_%, gaps_%\n')
    with open(blast_file, 'r') as f:
        new_query = False
        for line in f:
            if line.startswith('Query=') and not new_query:
                cluster = line.split('=')[1].strip()
                new_query = True
            if line.startswith('>') and new_query:
                line = line.split(',')[0]
                id = line.split(' ')[0][1:]
                full_name = ' '.join(line.split(' ')[1:]).strip()

            if line.startswith(' Score') and new_query:
                score = line.split(',')[0].split('(')[-1][:-1]
                e_value = line.split(',')[1].split('=')[1].strip()
            if line.startswith(' Identities') and new_query:
                identity = line.split(',')[0].split('(')[-1][:-2]
                gaps = line.split(',')[1].split('(')[-1][:-3]

                out.write(cluster + ',' + id + ',' + full_name + ',' + score + ',' + e_value + 
                        ',' + identity + ',' + gaps + '\n')
                new_query = False
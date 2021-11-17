#######################################################################
# Main function of the pipeline that runs all the scripts. 
#######################################################################

import argparse
import os
import sys 
import yaml


def parse_config_file(config_path):
    """
    This function reads in the config yaml from the config flag.
    Returns a dict of parameters to be used in the pipeline.
    """

    sys.stdout.write('Main: Checking the config file for errors. \n')
    # Read in the input yaml
    with open(config_path, 'r') as conf:
        input_config = yaml.load(conf.read(), yaml.FullLoader)

    # Make sure that all required fields are filled in and exist.
    required_fields = ['fastq', 'tmp', 'results', 'dependencies', 'parameters']
    required_parameters = ['kmer_len', 'umap_n_components', 'umap_n_neighbors', 'hdbscan_min_cluster_size',
                            'hdbscan_cluster_selection_method', 'refine_threshold', 'ram']

    for field in required_fields:
        if field == 'fastq':
            try:
                fastq = input_config[field]
                file = open(fastq, 'r')
            except FileNotFoundError:
                sys.stderr.write("The FASTQ file does not exist. Please check the path to your file.\n")
                sys.exit()
        elif field != 'parameters':
            if not os.path.isdir(input_config[field]):
                sys.stderr.write('The path for ' + field + ' does not exist. Please check the path to the folder.\n')
                sys.exit()
        elif field == 'refinement':
            if not isinstance(input_config[field], bool):
                sys.stderr.write('The refinement field requires a boolean value (True / False).\n')
                sys.exit()
        else:
            for param in required_parameters:
                if param == 'hdbscan_cluster_selection_method':
                    if input_config[field][param] not in ['leaf', 'eom']:
                        sys.stderr.write('The hdbscan_cluster_selection_method provided is not supported by hdbscan.'+
                                            ' Please choose either leaf or eom or check the spelling \n')
                        sys.exit()
                elif param == 'refine_threshold':
                    if not isinstance(input_config[field][param], float) or input_config[field][param] > 1:
                        sys.stderr.write('The threshold provided is not a float number or is greater than 1.\n')
                        sys.exit()
                else:
                    if not isinstance(input_config[field][param], int):
                        sys.stderr.write('The provided number for ' + param + ' is not a whole number.\n')
                        sys.exit()

    sys.stdout.write('Main: The config file looks correct. \n')
    return(input_config)        


def launch_pipeline(input_parameters):
    """
    This function launches all the scripts in the correct order.  The DAG can be viewed in the README.
    """

    sys.stdout.write('Main: Launching the pipeline.\n')
    current_dir = '/'.join(__file__.split('/')[:-1])
    fastq_file = input_parameters['fastq'].split('/')[-1]

    ## Run QC script.
    sys.stdout.write('Main: Launch Quality Control script.\n')
    os.system('python3 ' + current_dir + '/src/quality_control.py -fastq ' + input_parameters['fastq'] +
                     ' -out ' + input_parameters['tmp'])

    qc_file_path = input_parameters['tmp'] + fastq_file.split('.')[0] + '.fasta'
    # Check if the QC script outputted a correct file to correct directory.
    if not os.path.exists(qc_file_path):
        sys.stderr.write('Main: the QC script did not produce expected fasta file in the tmp directory.\n')
        sys.exit()

    ## Run kmer freq matrix script.
    fasta_file = fastq_file.split('.')[0] + '.fasta'
    sys.stdout.write('Main: Launch kmer frequency matrix construction.\n')
    os.system('julia ' + current_dir + '/src/main.jl ' + input_parameters['tmp'] + ' ' + fasta_file + ' ' +
                     input_parameters['tmp'] + ' ' + str(input_parameters['parameters']['kmer_len']))

    kmer_file_path = input_parameters['tmp'] + fastq_file.split('.')[0] + '_kmer_matrix.csv'
    # Check if the kmer freq script outputted a correct file to correct directory.
    if not os.path.exists(kmer_file_path):
        sys.stderr.write('Main: the kmer freq matrix script did not produce expected csv file in the tmp directory.\n')
        sys.exit()

    ## Run UMAP clustering script.
    sys.stdout.write('Main: Launch UMAP clustering script.\n')
    os.system('python3 ' + current_dir + '/src/umap_clustering.py -kmer ' + kmer_file_path + 
                ' -out ' + input_parameters['tmp'] + ' -nc ' + str(input_parameters['parameters']['umap_n_components']) +
                ' -nn ' + str(input_parameters['parameters']['umap_n_neighbors']) + ' -min_cs ' +
                str(input_parameters['parameters']['hdbscan_min_cluster_size']) + ' -csm ' + 
                input_parameters['parameters']['hdbscan_cluster_selection_method'] + ' -res_out ' +
                input_parameters['results'])

    clustering_file_path = input_parameters['tmp'] + fastq_file.split('.')[0] + '_kmer_matrix_clustering.csv'
    # Check if the UMAP clustering outputted a correct file to correct directory.
    if not os.path.exists(clustering_file_path):
        sys.stderr.write('Main: the UMAP clustering script did not produce expected csv file in the tmp directory.\n')
        sys.exit()

    if input_parameters['refinement']:
        ## Run refine clusters script.
        sys.stdout.write('Main: Launch refine clusters script.\n')
        os.system('python3 ' + current_dir + '/src/refine_clusters.py -fasta ' + qc_file_path + 
                    ' -clusters ' + clustering_file_path + ' -out ' + input_parameters['tmp'] + 
                    ' -ref_thr ' + str(input_parameters['parameters']['refine_threshold']))

        refined_file_path = input_parameters['tmp'] + fastq_file.split('.')[0] + '_refined_clusters.txt'
        # Check if the UMAP clustering outputted a correct file to correct directory.
        if not os.path.exists(refined_file_path):
            sys.stderr.write('Main: the refinement script did not produce expected txt file in the tmp directory.\n')
            sys.exit()

        ## Run consensus script.
        sys.stdout.write('Main: Launch define consensus script.\n')
        os.system('python3 ' + current_dir + '/src/make_consensus.py -fasta ' + qc_file_path + 
                    ' -clusters ' + refined_file_path + ' -tmp_out ' + input_parameters['tmp'] + 
                    ' -out ' + input_parameters['results'] + ' -dep ' + input_parameters['dependencies'])

        consensus_file_path = input_parameters['results'] + fastq_file.split('.')[0] + '_consensus.txt'
        # Check if the consensus outputted a correct file to correct directory.
        if not os.path.exists(consensus_file_path):
            sys.stderr.write('Main: the consensus script did not produce expected txt file in the tmp directory.\n')
            sys.exit()
    else:
        ## Run consensus script.
        sys.stdout.write('Main: Launch define consensus script.\n')
        os.system('python3 ' + current_dir + '/src/make_consensus.py -fasta ' + qc_file_path + 
                    ' -clusters ' + clustering_file_path + ' -tmp_out ' + input_parameters['tmp'] + 
                    ' -out ' + input_parameters['results'] + ' -dep ' + input_parameters['dependencies'])

        consensus_file_path = input_parameters['results'] + fastq_file.split('.')[0] + '_consensus.txt'
        # Check if the consensus outputted a correct file to correct directory.
        if not os.path.exists(consensus_file_path):
            sys.stderr.write('Main: the consensus script did not produce expected txt file in the tmp directory.\n')
            sys.exit()

    ## Run BLASTN parsing script.
    sys.stdout.write('Main: Launch BLASTN output parsing script.\n')
    blast_file_path = input_parameters['results'] + fastq_file.split('.')[0] + '_blastn_result.fa'
    os.system('python3 ' + current_dir + '/src/parse_blast_output.py -blast ' + blast_file_path +  
                ' -out ' + input_parameters['results'])

    parsed_result_file_path = input_parameters['results'] + fastq_file.split('.')[0] + '_final_cluster_identities.csv'
    # Check if the UMAP clustering outputted a correct file to correct directory.
    if not os.path.exists(parsed_result_file_path):
        sys.stderr.write('Main: the BLAST parsing script did not produce expected txt file in the tmp directory.\n')
        sys.exit()


def main():
    """
    This is the main function of the consensus generating pipeline
    """

    sys.stdout.write('Main: Start the main function. \n')
    parser = argparse.ArgumentParser(prog='consensus generating pipeline',
                                     description='Create a consensus reads for Nanopore generated reads',
                                     epilog='Contact Ada Madejska (amaddejska@ucsb.edu) if you encounter '
                                     'any problems while running this pipeline')
    inputs = parser.add_mutually_exclusive_group(required=True)
    inputs.add_argument('--config_file', dest='config_file', help='Config file to be used in the '
                        'run.', type=str, default=None)

    parser.add_argument('--max-cores-per-job', dest='max_cores', help='Maximum cores to use per '
                        'job. This value should be set to the number of cpus on the '
                        'smallest node in a cluster.',
                        type=int, required=False, default=None)

    params = parser.parse_args()
    params.config_file = os.path.abspath(params.config_file)
    if params.max_cores:
        max_cores = params.max_cores
    else:
        max_cores = 1
    
    input_parameters = parse_config_file(params.config_file)

    launch_pipeline(input_parameters)

if __name__ == '__main__':
    main()
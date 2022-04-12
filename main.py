#######################################################################
# Main function of the pipeline that runs all the scripts. 
#######################################################################

import argparse
import multiprocessing as mp
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
                            'hdbscan_cluster_selection_method', 'reads_per_job']

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
                else:
                    if not isinstance(input_config[field][param], int):
                        sys.stderr.write('The provided number for ' + param + ' is not a whole number.\n')
                        sys.exit()

    # Set up the name of the data that we will use for naming all the files.
    input_config['name'] = input_config['fastq'].split('/')[-1].split('.')[0]

    sys.stdout.write('Main: The config file looks correct. \n')
    return(input_config)        


def setup_environment(input_parameters):
    """
    This function creates a separate folder for temporary files in the tmp folder and results folder.
         Will make multiprocessing much easier to manage. 
    """

    sys.stdout.write('Main: Creating directories. \n')
    # Create a separate folder in the tmp directory.
    tmpr_folder = input_parameters['tmp'] + '/' + input_parameters['name'] + '_tmp'
    if not os.path.isdir(tmpr_folder):
        os.mkdir(tmpr_folder)

        input_parameters['tmp'] = tmpr_folder + '/'
    else:
        sys.stdout.write('Main: Cowardly refuses to overwrite the existing folder ' + tmpr_folder + '\n')
        sys.stdout.write('Please, delete or move the folder to a different location \n')
        sys.exit()

    # Create a separate folder in the results directory
    results_folder = input_parameters['results'] + '/' + input_parameters['name'] + '_results'
    if not os.path.isdir(results_folder):
        os.mkdir(results_folder)

        input_parameters['results'] = results_folder + '/'
    else:
        sys.stdout.write('Main: Cowardly refuses to overwrite the existing folder ' + results_folder + '\n')
        sys.stdout.write('Please, delete or move the folder to a different location \n')
        sys.exit()

    # Create a folder for subfiles of the FASTQ file. 
    fastq_folder = tmpr_folder + '/fastq_files'
    if not os.path.isdir(fastq_folder):
        os.mkdir(fastq_folder)

        input_parameters['fastq_folder'] = fastq_folder + '/'
    else:
        sys.stdout.write('Main: Cowardly refuses to overwrite the existing folder ' + fastq_folder + '\n')
        sys.stdout.write('Please, delete or move the folder to a different location \n')
        sys.exit()

    sys.stdout.write('Main: Directories created sucessfully. \n')
    return(input_parameters)


def divide_fastq_file(input_parameters, qc_fastq_path):
    """
    This function divides the original FASTQ file into a smaller files with # of reads indicated by the user.
    Each file will be analyzed separately during multiprocessing step.
    """
    sys.stdout.write('Main: Divide fastq file into smaller jobs.\n')
    # split the files using linux split function
    os.system('split -l ' + str(input_parameters['parameters']['reads_per_job']*2) + ' ' + qc_fastq_path +
              ' ' + input_parameters['fastq_folder'] + input_parameters['name'])

    fastq_list = [input_parameters['fastq_folder'] + i for i in os.listdir(input_parameters['fastq_folder'])]
    return(fastq_list)


def run_QC(input_parameters):
    """
    This function is a first step of the pipeline. It runs the QC step on the whole fastq file before
    it gets split into smaller files.
    """
    sys.stdout.write('Main: Launching the pipeline.\n')
    current_dir = '/'.join(__file__.split('/')[:-1])
    
    sys.stdout.write('Main: Launch quality control script. \n')
    ## Run QC script.
    os.system('python3 ' + current_dir + '/src/quality_control.py -fastq ' + input_parameters['fastq'] +
                     ' -out ' + input_parameters['tmp'])

    qc_file_path = input_parameters['tmp'] + input_parameters['name'] + '.fasta'
    # Check if the QC script outputted a correct file to correct directory.
    if not os.path.exists(qc_file_path):
        sys.stderr.write('Main: the QC script did not produce expected fasta file in the tmp directory.\n')
        sys.exit()

    return(qc_file_path)


def launch_pipeline(input_parameters, qc_file_path, job_name):
    """
    This function launches all the scripts in the correct order.  The DAG can be viewed in the README.
    """
    current_dir = '/'.join(__file__.split('/')[:-1])
    input_parameters['name'] = qc_file_path.split('/')[-1]

    sys.stdout.write('Main: Working on file ' + qc_file_path.split('/')[-1] + '\n')

    ## Run kmer freq matrix script.
    fasta_file = qc_file_path.split('/')[-1]
    sys.stdout.write('Main: ' + input_parameters['name'] + ' Launch kmer frequency matrix construction.\n')
    os.system('julia ' + current_dir + '/src/main.jl ' + input_parameters['fastq_folder'] + ' ' + 
                     fasta_file + ' ' + input_parameters['tmp'] + ' ' + 
                     str(input_parameters['parameters']['kmer_len']))

    kmer_file_path = input_parameters['tmp'] + input_parameters['name'] + '_kmer_matrix.csv'
    # Check if the kmer freq script outputted a correct file to correct directory.
    if not os.path.exists(kmer_file_path):
        sys.stderr.write('Main: ' + input_parameters['name'] + ' the kmer freq matrix script did not produce expected csv file in the tmp directory.\n')
        sys.exit()

    ## Run UMAP clustering script.
    sys.stdout.write('Main: ' + input_parameters['name'] + ' Launch UMAP clustering script.\n')
    os.system('python3 ' + current_dir + '/src/umap_clustering.py -kmer ' + kmer_file_path + 
                ' -out ' + input_parameters['tmp'] + ' -nc ' + str(input_parameters['parameters']['umap_n_components']) +
                ' -nn ' + str(input_parameters['parameters']['umap_n_neighbors']) + ' -min_cs ' +
                str(input_parameters['parameters']['hdbscan_min_cluster_size']) + ' -csm ' + 
                input_parameters['parameters']['hdbscan_cluster_selection_method'] + ' -res_out ' +
                input_parameters['results'] + ' -fn ' + input_parameters['name'] + ' -jn ' + 
                job_name)

    clustering_file_path = input_parameters['tmp'] + input_parameters['name'] + '_clustering.csv'
    # Check if the UMAP clustering outputted a correct file to correct directory.
    if not os.path.exists(clustering_file_path):
        sys.stderr.write('Main: ' + input_parameters['name'] + ' the UMAP clustering script did not produce expected csv file in the tmp directory.\n')
        sys.exit()
    
    if input_parameters['refinement']:
        
        ## Run refine clusters script.
        sys.stdout.write('Main: ' + input_parameters['name'] + ' Launch refine clusters script.\n')
        os.system('python3 ' + current_dir + '/src/refine_clusters.py -fasta ' + qc_file_path + 
                    ' -clusters ' + clustering_file_path + ' -out ' + input_parameters['tmp'] + 
                    ' -name ' + input_parameters['name'] + ' -dep ' + input_parameters['dependencies'])
        
        
        refined_file_path = input_parameters['tmp'] + input_parameters['name'] + '_refined_clusters.txt'
        # Check if the UMAP clustering outputted a correct file to correct directory.
        if not os.path.exists(refined_file_path):
            sys.stderr.write('Main: ' + input_parameters['name'] + ' the refinement script did not produce expected txt file in the tmp directory.\n')
            sys.exit()
        
        ## Run consensus script.
        sys.stdout.write('Main: ' + input_parameters['name'] + ' Launch define consensus script.\n')
        os.system('python3 ' + current_dir + '/src/make_consensus.py -fasta ' + qc_file_path + 
                    ' -clusters ' + refined_file_path + ' -tmp_out ' + input_parameters['tmp'] + 
                    ' -out ' + input_parameters['tmp'] + ' -dep ' + input_parameters['dependencies'] + 
                    ' -name ' + input_parameters['name'])

        consensus_file_path = input_parameters['tmp'] + input_parameters['name'] + '_consensus.txt'
        # Check if the consensus outputted a correct file to correct directory.
        if not os.path.exists(consensus_file_path):
            sys.stderr.write('Main: ' + input_parameters['name'] + ' the consensus script did not produce expected txt file in the tmp directory.\n')
            sys.exit()
    else:
        ## Run consensus script.
        sys.stdout.write('Main: ' + input_parameters['name'] + ' Launch define consensus script.\n')
        os.system('python3 ' + current_dir + '/src/make_consensus.py -fasta ' + qc_file_path + 
                    ' -clusters ' + clustering_file_path + ' -tmp_out ' + input_parameters['tmp'] + 
                    ' -out ' + input_parameters['tmp'] + ' -dep ' + input_parameters['dependencies'] +
                    ' -name ' + input_parameters['name'])

        consensus_file_path = input_parameters['tmp'] + input_parameters['name'] + '_consensus.txt'
        # Check if the consensus outputted a correct file to correct directory.
        if not os.path.exists(consensus_file_path):
            sys.stderr.write('Main: ' + input_parameters['name'] + ' the consensus script did not produce expected txt file in the tmp directory.\n')
            sys.exit()

    return()


def final_clustering(input_parameters, qc_fasta_path, job_name):
    """
    After analyzing the small files separately, we need to do one last cummulative analysis of all the consensus 
        we found and noise reads we've collected.
    This will result in one comprehensive file with all cluster identities.  
    """
    current_dir = '/'.join(__file__.split('/')[:-1])

    sys.stdout.write('Main: Initiating the final analysis of all jobs. \n')

    sys.stdout.write('Main: Concatenate consensus sequences and noise reads. \n')
    os.system('python3 ' + current_dir + '/src/make_noise_fastq.py -fastq ' + input_parameters['fastq'] +
              ' -noise ' + input_parameters['tmp'] + '/' + job_name + '_noise.csv' + 
              ' -out ' + input_parameters['tmp'] + ' -name ' + job_name)

    noise_reads_file_path = input_parameters['tmp'] + job_name + '_noise.fasta'
    # Check if the make_noise_fastq outputted a correct file to correct directory.
    if not os.path.exists(noise_reads_file_path):
        sys.stderr.write(noise_reads_file_path)
        sys.stderr.write('Main: the make_noise_fastq script did not produce expected txt file in the tmp directory.\n')
        sys.exit()

    ## Run kmer freq matrix script.
    sys.stdout.write('Main: Launch final kmer frequency matrix construction.\n')
    os.system('julia ' + current_dir + '/src/main.jl ' + input_parameters['tmp'] + ' ' + 
                     job_name + '_noise.fasta' + ' ' + input_parameters['tmp'] + ' ' + 
                     str(input_parameters['parameters']['kmer_len']))

    kmer_file_path = input_parameters['tmp'] + job_name + '_noise_kmer_matrix.csv'
    # Check if the kmer freq script outputted a correct file to correct directory.
    if not os.path.exists(kmer_file_path):
        sys.stderr.write('Main: the kmer freq matrix script did not produce expected csv file in the tmp directory.\n')
        sys.exit()

    ## Run UMAP clustering script.
    sys.stdout.write('Main: Launch final UMAP clustering script.\n')
    os.system('python3 ' + current_dir + '/src/umap_clustering.py -kmer ' + kmer_file_path + 
                ' -out ' + input_parameters['tmp'] + ' -nc ' + str(input_parameters['parameters']['umap_n_components']) +
                ' -nn ' + str(input_parameters['parameters']['umap_n_neighbors']) + ' -min_cs ' +
                str(input_parameters['parameters']['hdbscan_min_cluster_size']) + ' -csm ' + 
                input_parameters['parameters']['hdbscan_cluster_selection_method'] + ' -res_out ' +
                input_parameters['results'] + ' -fn noise -jn noise ')

    clustering_file_path = input_parameters['tmp'] + 'noise_clustering.csv'
    # Check if the UMAP clustering outputted a correct file to correct directory.
    if not os.path.exists(clustering_file_path):
        sys.stderr.write('Main: the UMAP clustering script did not produce expected csv file in the tmp directory.\n')
        sys.exit()

    ## Run consensus script.
    sys.stdout.write('Main: Launch define consensus script.\n')
    os.system('python3 ' + current_dir + '/src/make_consensus.py -fasta ' + qc_fasta_path + 
                ' -clusters ' + clustering_file_path + ' -tmp_out ' + input_parameters['tmp'] + 
                ' -out ' + input_parameters['tmp'] + ' -dep ' + input_parameters['dependencies'] +
                ' -name noise')

    consensus_file_path = input_parameters['tmp'] + 'noise_consensus.txt'
    # Check if the consensus outputted a correct file to correct directory.
    if not os.path.exists(consensus_file_path):
        sys.stderr.write('Main: the consensus script did not produce expected txt file in the tmp directory.\n')
        sys.exit()

    sys.stdout.write('Main: Launch concatenate consensus files script.\n')
    os.system('python3 ' + current_dir + '/src/concatenate_consensus.py -cons ' + input_parameters['tmp'] + 
                ' -out ' + input_parameters['results'] + ' -name ' + job_name)

    final_consensus_file_path = input_parameters['results'] + job_name + '_final_consensus.fasta'
    # Check if the concatenation outputted a correct file to correct directory.
    if not os.path.exists(final_consensus_file_path):
        sys.stderr.write('Main: the concatenation script did not produce expected txt file in the results directory.\n')
        sys.exit()

    ## Run BLASTN on the consensus sequences.
    sys.stdout.write('\nMain: Launch BLASTN script. \n')
    os.system('python3 ' + current_dir + '/src/run_blastn.py -consensus ' + final_consensus_file_path +
             ' -out ' + input_parameters['tmp'] + ' -dep ' + input_parameters['dependencies'] + 
             ' -name ' + job_name)

    # Check if the BLASTN outputted a correct file to correct directory.
    blastn_result_file_path = input_parameters['tmp'] + job_name + '_blastn_result.fa'
    if not os.path.exists(blastn_result_file_path):
        sys.stderr.write(blastn_result_file_path)
        sys.stderr.write('Main: the BLASTN script did not produce expected txt file in the tmp directory.\n')
        sys.exit()
    
    ## Run BLASTN parsing script.
    sys.stdout.write('Main: Launch BLASTN output parsing script.\n')
    os.system('python3 ' + current_dir + '/src/parse_blast_output.py -blast ' + blastn_result_file_path +  
                ' -out ' + input_parameters['results'] + ' -name ' + job_name)

    parsed_result_file_path = input_parameters['results'] + job_name + '_final_cluster_identities.csv'
    # Check if the BLASTN parsing outputted a correct file to correct directory.
    if not os.path.exists(parsed_result_file_path):
        sys.stderr.write(parsed_result_file_path)
        sys.stderr.write('Main: the BLAST parsing script did not produce expected txt file in the results directory.\n')
        sys.exit()

    return()


def collect_result(result):
    global results
    results.append(result)


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

    params = parser.parse_args()
    params.config_file = os.path.abspath(params.config_file)
    
    input_parameters = parse_config_file(params.config_file)

    setup_environment(input_parameters)

    qc_fasta_path = run_QC(input_parameters)

    fastq_file_list = divide_fastq_file(input_parameters, qc_fasta_path)

    job_name = input_parameters['name']

    # Mulitprocessing step of the pipeline.
    pool = mp.Pool(mp.cpu_count())
    for i, file in enumerate(fastq_file_list):
        pool.apply_async(launch_pipeline, args=(input_parameters, file, job_name), callback=collect_result)
    
    pool.close() 
    pool.join()

    final_clustering(input_parameters, qc_fasta_path, job_name)

    sys.stdout.write('Main: Finished sucessfully.\n')


if __name__ == '__main__':
    results = []
    main()
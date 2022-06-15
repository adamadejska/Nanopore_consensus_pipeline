
# Parts from https://stackoverflow.com/questions/65492490/how-to-conda-install-cuda-enabled-pytorch-in-a-docker-container 
# Use nvidia/cuda image
FROM nvidia/cuda:10.2-cudnn7-devel-ubuntu18.04

# set bash as current shell
#RUN chsh -s /bin/bash
#SHELL ["/bin/bash", "-c"]

# install anaconda
RUN apt-key adv --fetch-keys http://developer.download.nvidia.com/compute/cuda/repos/ubuntu1604/x86_64/3bf863cc.pub

RUN apt-get update
RUN apt-get install -y wget bzip2 ca-certificates libglib2.0-0 libxext6 libsm6 libxrender1 git mercurial subversion && \
        apt-get clean
RUN wget --quiet https://repo.anaconda.com/archive/Anaconda3-2020.02-Linux-x86_64.sh -O ~/anaconda.sh && \
        /bin/bash ~/anaconda.sh -b -p /opt/conda && \
        rm ~/anaconda.sh && \
        ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
        echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
        find /opt/conda/ -follow -type f -name '*.a' -delete && \
        find /opt/conda/ -follow -type f -name '*.js.map' -delete && \
        /opt/conda/bin/conda clean -afy

# set path to conda
ENV PATH /opt/conda/bin:$PATH


# setup conda virtual environment
COPY ./requirements.yml /tmp/requirements.yml
RUN conda update conda 
RUN conda env create -f /tmp/requirements.yml

RUN echo "conda activate deepdpm-env" >> ~/.bashrc
ENV PATH /opt/conda/envs/deepdpm-env/bin:$PATH
ENV CONDA_DEFAULT_ENV $deepdpm-env

#RUN conda --version

# Clone deepDPM code from github
RUN git clone https://github.com/BGU-CS-VIL/DeepDPM.git

#COPY ./rundeepdpm.sh DeepDPM/rundeepdpm.sh


# Run image
# docker build -t deepdpm .
# docker run -it --rm --name deepdpm deepdpm
# ctrl + p, ctrl + q <- exit interactive mode
# docker exec deepdpm python DeepDPM/DeepDPM.py --dataset synthetic --latent_dim 2 --log_emb every_n_epochs --log_emb_every 1

# Run deepdpm
# 1. conda activate deepdpm-env    # no need, the image already activates the env
# 2. python DeepDPM/DeepDPM.py --dataset synthetic --latent_dim 2 --log_emb every_n_epochs --log_emb_every 1
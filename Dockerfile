FROM tutum/curl:latest

# Install miniconda to /miniconda
RUN curl -LO http://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    bash Miniconda3-latest-Linux-x86_64.sh -p /miniconda -b && \
    rm Miniconda3-latest-Linux-x86_64.sh

# set the environment
ENV PATH="/miniconda/bin:$PATH" \
 LC_ALL=C

# Install miniconda environment
RUN conda update --all && \
    conda config --add channels bioconda && \
    conda config --add channels conda-forge && \
    conda config --add channels anaconda && \
    conda create -n WGS_COVID19 python=3.6.8 pysam=0.15.2 nextflow=19.10.0 \
    simplejson=3.17.0 pandas=1.0.1 matplotlib=3.1.2 graphviz=2.42.3

## add pangolin
## commit a0f2dbab550daaf3c870fc5ec5a32ded948c06f4
#RUN conda install git && \
#    cd / && \
#    git clone https://github.com/cov-lineages/pangolin.git && \
#    cd /pangolin && \
#    conda env create -f environment.yml
#
#RUN echo "source activate pangolin" > /etc/.bashrc && \
#    cd /pangolin && \
#    /miniconda/envs/pangolin/bin/python setup.py install && \
#    /miniconda/envs/pangolin/bin/pangolin -v

## activate by bash mycoprofiler env
RUN echo "source activate WGS_COVID19" > ~/.bashrc
ENV PATH=/miniconda/envs/WGS_COVID19/bin:${PATH}

# set shell
SHELL ["/bin/bash", "-c"]

# add codebase to docker
ADD ./ /workflow
WORKDIR /workflow
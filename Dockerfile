FROM tutum/curl:latest

# Install miniconda to /miniconda
RUN curl -LO http://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    bash Miniconda3-latest-Linux-x86_64.sh -p /miniconda -b && \
    rm Miniconda3-latest-Linux-x86_64.sh

# set the environment
ENV PATH="/miniconda/bin:$PATH" LC_ALL=C

# Install miniconda environment
RUN conda update --all && \
    conda config --add channels bioconda && \
    conda config --add channels conda-forge && \
    conda config --add channels anaconda && \
    conda create -n WGS_COVID19 python=3.6.8 pysam=0.15.2 nextflow=20.10.0 \
    simplejson=3.17.0 pandas=1.0.1 matplotlib=3.1.2 graphviz=2.42.3

## activate by bash mycoprofiler env
RUN echo "source activate WGS_COVID19" > ~/.bashrc
ENV PATH=/miniconda/envs/WGS_COVID19/bin:${PATH}

# set shell
SHELL ["/bin/bash", "-c"]

# add codebase to docker
# .dockerignore arranges which will be excluded in the image build
ADD ./ /workflow
WORKDIR /workflow
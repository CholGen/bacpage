FROM mambaorg/micromamba:1.4.9
LABEL authors="Nate Matteson"
MAINTAINER Nate M <natem@scripps.edu>

WORKDIR /home/mambauser/
RUN micromamba install -y -n base -c conda-forge git
ARG MAMBA_DOCKERFILE_ACTIVATE=1
RUN micromamba install -y -n base -c bioconda -c condaforge bacpage && \
    micromamba clean --all --yes

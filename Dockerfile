FROM mambaorg/micromamba:1.4.9
LABEL authors="Nate Matteson"
MAINTAINER Nate M <natem@scripps.edu>

WORKDIR /home/mambauser/
COPY --chown=$MAMBA_USER:$MAMBA_USER . ./bacpage/

WORKDIR /home/mambauser/bacpage/
RUN micromamba install -y -n base -f environment_docker.yaml \
    && micromamba clean --all --yes
ARG MAMBA_DOCKERFILE_ACTIVATE=1
RUN pip install .

N# Specify path to gain access to conda environment without terminal model
ENV PATH="/opt/conda/bin:/opt/conda/condabin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin"

WORKDIR /home/mambauser/

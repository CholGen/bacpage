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

WORKDIR /home/mambauser/

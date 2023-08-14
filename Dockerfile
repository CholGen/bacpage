FROM mambaorg/micromamba:1.4.9
LABEL authors="Nate Matteson"
MAINTAINER Nate M <natem@scripps.edu>

RUN micromamba install -y -n base -c conda-forge git
RUN git clone --depth=1 https://github.com/watronfire/eureka.git
RUN micromamba install -y -n base -f eureka/environment.yaml && \
    micromamba clean --all --yes

ENTRYPOINT ["top", "-b"]
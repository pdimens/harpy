# HARPY qc workflow

FROM mambaorg/micromamba

#USER root

#RUN apt-get update && apt-get install -y perl \
#    && rm -rf /var/lib/apt/lists/*

#USER $MAMBA_USER

#COPY --chown=$MAMBA_USER:$MAMBA_USER qc.yaml /tmp/env.yaml

#RUN micromamba install --name base --yes --file /tmp/env.yaml \
#    && micromamba clean --all --yes

RUN micromamba install --yes --name base -c bioconda -c conda-forge \
    falco \ 
    fastp \ 
    multiqc \ 
    r-base \ 
    r-dplyr \ 
    r-dt \ 
    r-flexdashboard \ 
    r-ggplot2 \ 
    r-magrittr \ 
    r-plotly \ 
    r-tidyr \ 
    r-viridislite

#ARG MAMBA_DOCKERFILE_ACTIVATE=1
RUN which perl

COPY --chown=$MAMBA_USER:$MAMBA_USER ../../workflow/scripts/* /home/$MAMBA_USER/scripts

RUN chmod +x /home/$MAMBA_USER/scripts/* && cp home/scripts/* $CONDA_PREFIX/bin/
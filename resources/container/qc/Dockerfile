# HARPY QC workflow

FROM mambaorg/micromamba

USER root

RUN apt-get update && apt-get install -y perl \
    && rm -rf /var/lib/apt/lists/*

USER $MAMBA_USER

#COPY --chown=$MAMBA_USER:$MAMBA_USER qc.lock /tmp/env.lock

#RUN micromamba install --name base --yes --file /tmp/env.lock \
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
    r-knitr \ 
    r-magrittr \ 
    r-plotly \ 
    r-rmarkdown \ 
    r-tidyr \ 
    r-viridislite

ARG MAMBA_DOCKERFILE_ACTIVATE=1

COPY /workdir/pdimens/harpy/workflow/scripts/ /home/

RUN chmod +x /home/scripts/* && cp home/scripts/* $CONDA_PREFIX/bin/
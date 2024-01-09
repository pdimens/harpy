# HARPY qc workflow

FROM mambaorg/micromamba

RUN micromamba install --yes --name base -c bioconda -c conda-forge \
    falco \ 
    fastp \ 
    multiqc \ 
    r-base \ 
    r-dt \ 
    r-flexdashboard \ 
    r-ggplot2 \ 
    r-plotly \ 
    r-tidyr

#ARG MAMBA_DOCKERFILE_ACTIVATE=1

ADD --chown=$MAMBA_USER:$MAMBA_USER workflow/scripts/ $CONDA_PREFIX/bin/
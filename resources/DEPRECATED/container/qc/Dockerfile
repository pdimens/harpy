# HARPY workflows

FROM mambaorg/micromamba AS qc
RUN micromamba install --yes --name base -c bioconda -c conda-forge \
    falco \ 
    fastp \ 
    multiqc \ 
    r-base \ 
    r-dt \ 
    r-flexdashboard \ 
    r-ggplot2 \ 
    r-plotly \ 
    r-tidyr && \
    micromamba clean --all --yes
ARG MAMBA_DOCKERFILE_ACTIVATE=1
COPY --chown=$MAMBA_USER:$MAMBA_USER workflow/scripts/ workflow/report/ $CONDA_PREFIX/bin/

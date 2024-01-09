FROM mambaorg/micromamba AS impute
RUN micromamba install --yes --name base -c bioconda -c conda-forge \
    bcftools \ 
    bioconductor-complexheatmap \ 
    llvm-openmp \ 
    multiqc \ 
    r-base \ 
    r-circlize \ 
    r-dt \ 
    r-flexdashboard \ 
    r-ggplot2 \ 
    r-ggridges \ 
    r-plotly \ 
    r-stitch \ 
    r-tidyr && \
    micromamba clean --all --yes
ARG MAMBA_DOCKERFILE_ACTIVATE=1
COPY --chown=$MAMBA_USER:$MAMBA_USER workflow/scripts/ workflow/report/ $CONDA_PREFIX/bin/

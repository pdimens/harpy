FROM mambaorg/micromamba AS snp
RUN micromamba install --yes --name base -c bioconda -c conda-forge \
    bcftools \ 
    bioconductor-complexheatmap \ 
    freebayes \ 
    pysam \ 
    r-base \ 
    r-circlize \ 
    r-dt \ 
    r-flexdashboard \ 
    r-ggplot2 \ 
    r-ggridges \ 
    r-plotly \ 
    r-tidyr \ 
    sambamba \ 
    samtools \ 
    seqtk \ 
    tabix && \
    micromamba clean --all --yes
ARG MAMBA_DOCKERFILE_ACTIVATE=1
COPY --chown=$MAMBA_USER:$MAMBA_USER src/harpy/scripts/ src/harpy/reports/ $CONDA_PREFIX/bin/

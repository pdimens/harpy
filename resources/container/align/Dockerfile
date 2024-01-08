# HARPY align workflow

FROM mambaorg/micromamba

#COPY --chown=$MAMBA_USER:$MAMBA_USER align.lock /tmp/env.lock

#RUN micromamba install --name base --yes --file /tmp/env.lock \
#    && micromamba clean --all --yes


RUN micromamba install --yes --name base -c bioconda -c conda-forge \
    bwa \ 
    ema \
    icu \
    libzlib \
    multiqc \
    llvm-openmp \ 
    pysam \
    r-base \
    r-circlize \
    r-dplyr \
    r-dt \ 
    r-flexdashboard \ 
    r-ggplot2 \ 
    r-ggridges \ 
    r-magrittr \ 
    r-plotly \ 
    r-tidyr \ 
    r-viridislite \ 
    sambamba \ 
    samtools \ 
    seqtk \ 
    tabix \ 
    xz

RUN micromamba clean --all --yes



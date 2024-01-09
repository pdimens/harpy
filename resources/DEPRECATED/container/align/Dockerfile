
FROM mambaorg/micromamba AS align
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
    r-dt \ 
    r-flexdashboard \ 
    r-ggplot2 \ 
    r-ggridges \ 
    r-plotly \ 
    r-tidyr \ 
    sambamba \ 
    samtools \ 
    seqtk \ 
    tabix \ 
    xz && \
    micromamba clean --all --yes
ARG MAMBA_DOCKERFILE_ACTIVATE=1
COPY --chown=$MAMBA_USER:$MAMBA_USER workflow/scripts/ workflow/report/ $CONDA_PREFIX/bin/

FROM mambaorg/micromamba
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
    r-knitr \ 
    r-magrittr \ 
    r-plotly \ 
    r-rmarkdown \ 
    r-tidyr \ 
    r-viridislite \ 
    sambamba \ 
    samtools \ 
    seqtk \ 
    tabix \ 
    xz

RUN micromamba clean --all --yes



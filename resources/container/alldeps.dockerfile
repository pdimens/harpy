FROM mambaorg/micromamba
RUN micromamba install --yes --name base -c bioconda -c conda-forge \
    bcftools \ 
    bioconductor-complexheatmap \ 
    bwa \ 
    ema \ 
    falco \ 
    fastp \ 
    freebayes \ 
    hapcut2 \ 
    icu \ 
    libzlib \ 
    leviathan \ 
    llvm-openmp \ 
    multiqc \ 
    naibr-plus \ 
    python =3.10 \ 
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
    r-stitch \ 
    r-tidyr \ 
    r-viridislite \ 
    sambamba \ 
    samtools \ 
    seqtk \ 
    tabix \ 
    whatshap \ 
    xz

RUN micromamba clean --all --yes



FROM mambaorg/micromamba
RUN micromamba install --yes --name base -c bioconda -c conda-forge \
    bcftools \ 
    bioconductor-complexheatmap \ 
    freebayes \ 
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
    tabix

RUN micromamba clean --all --yes



# HARPY impute workflow

FROM mambaorg/micromamba

#COPY --chown=$MAMBA_USER:$MAMBA_USER impute.lock /tmp/env.lock

#RUN micromamba install --name base --yes --file /tmp/env.lock \
#    && micromamba clean --all --yes

RUN micromamba install --yes --name base -c bioconda -c conda-forge \
    bcftools \ 
    bioconductor-complexheatmap \ 
    llvm-openmp \ 
    multiqc \ 
    r-base \ 
    r-circlize \ 
    r-dplyr \ 
    r-dt \ 
    r-ggplot2 \ 
    r-ggridges \ 
    r-magrittr \ 
    r-plotly \ 
    r-rmarkdown \ 
    r-stitch \ 
    r-tidyr \ 
    r-viridislite

RUN micromamba clean --all --yes



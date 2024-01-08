# HARPY varaints sv workflow

FROM mambaorg/micromamba

#COPY --chown=$MAMBA_USER:$MAMBA_USER sv.lock /tmp/env.lock

#RUN micromamba install --name base --yes --file /tmp/env.lock \
#    && micromamba clean --all --yes

RUN micromamba install --yes --name base -c bioconda -c conda-forge \
    bcftools \ 
    bioconductor-complexheatmap \ 
    leviathan \ 
    naibr-plus \ 
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
    whatshap \ 
    xz

RUN micromamba clean --all --yes



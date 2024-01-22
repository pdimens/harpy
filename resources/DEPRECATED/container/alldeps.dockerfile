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
COPY --chown=$MAMBA_USER:$MAMBA_USER workflow/scripts/ workflow/report/ $CONDA_PREFIX/bin/

FROM mambaorg/micromamba AS sv
RUN micromamba install --yes --name base -c bioconda -c conda-forge \
    bcftools \ 
    bioconductor-complexheatmap \ 
    leviathan \ 
    naibr-plus \ 
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
    whatshap \ 
    xz && \
    micromamba clean --all --yes
ARG MAMBA_DOCKERFILE_ACTIVATE=1
COPY --chown=$MAMBA_USER:$MAMBA_USER workflow/scripts/ workflow/report/ $CONDA_PREFIX/bin/

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

FROM mambaorg/micromamba AS phase
RUN micromamba install --yes --name base -c bioconda -c conda-forge \
    bcftools \ 
    bioconductor-complexheatmap \ 
    hapcut2 \ 
    multiqc \ 
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
    samtools && \
    micromamba clean --all --yes
ARG MAMBA_DOCKERFILE_ACTIVATE=1
COPY --chown=$MAMBA_USER:$MAMBA_USER workflow/scripts/ workflow/report/ $CONDA_PREFIX/bin/

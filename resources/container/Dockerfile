FROM condaforge/mambaforge:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="d1f47fc29da7b5a95fe4d2add88c3a7dd6977254e24239fe934fd3f681c0a6fd"

# Step 1: Retrieve conda environments

# Conda environment:
#   source: .harpy_envs/align.yaml
#   prefix: /conda-envs/4538a63544ffaad3a6e49e8d25bfc6cb
#   name: align
#   channels:
#     - bioconda
#     - conda-forge
#     - defaults
#   dependencies:
#     - bioconda::bwa
#     - bioconda::ema
#     - conda-forge::icu
#     - conda-forge::libzlib
#     - bioconda::minimap2
#     - bioconda::samtools=1.20
#     - bioconda::seqtk
#     - bioconda::tabix
#     - conda-forge::xz
RUN mkdir -p /conda-envs/4538a63544ffaad3a6e49e8d25bfc6cb
COPY .harpy_envs/align.yaml /conda-envs/4538a63544ffaad3a6e49e8d25bfc6cb/environment.yaml

# Conda environment:
#   source: .harpy_envs/phase.yaml
#   prefix: /conda-envs/4914663d751c6d6b4e3b0b349fc0b652
#   name: phase
#   channels:
#     - bioconda
#     - conda-forge
#     - defaults
#   dependencies:
#     - bioconda::hapcut2
#     - bioconda::whatshap
RUN mkdir -p /conda-envs/4914663d751c6d6b4e3b0b349fc0b652
COPY .harpy_envs/phase.yaml /conda-envs/4914663d751c6d6b4e3b0b349fc0b652/environment.yaml

# Conda environment:
#   source: .harpy_envs/qc.yaml
#   prefix: /conda-envs/6e3a7dc7152f45301056861193dd3db5
#   name: qc
#   channels:
#     - bioconda
#     - conda-forge
#     - defaults
#   dependencies:
#     - bioconda::falco
#     - bioconda::fastp
#     - bioconda::multiqc
#     - bioconda::pysam=0.22
RUN mkdir -p /conda-envs/6e3a7dc7152f45301056861193dd3db5
COPY .harpy_envs/qc.yaml /conda-envs/6e3a7dc7152f45301056861193dd3db5/environment.yaml

# Conda environment:
#   source: .harpy_envs/r.yaml
#   prefix: /conda-envs/f769d2b509f06bc58137776df3f21936
#   name: r
#   channels:
#     - bioconda
#     - conda-forge
#     - defaults
#   dependencies:
#     - conda-forge::r-xml2
#     - bioconda::bioconductor-complexheatmap
#     - conda-forge::r-highcharter
#     - conda-forge::r-circlize
#     - r::r-biocircos
#     - conda-forge::r-dt
#     - conda-forge::r-flexdashboard
#     - conda-forge::r-ggplot2
#     - conda-forge::r-ggridges
#     - conda-forge::r-plotly
#     - conda-forge::r-tidyr
#     - bioconda::r-stitch
RUN mkdir -p /conda-envs/f769d2b509f06bc58137776df3f21936
COPY .harpy_envs/r.yaml /conda-envs/f769d2b509f06bc58137776df3f21936/environment.yaml

# Conda environment:
#   source: .harpy_envs/simulations.yaml
#   prefix: /conda-envs/6b5e3ed520a8389d84dfb0a12d5f07d4
#   name: simulations
#   channels:
#     - bioconda
#     - conda-forge
#     - defaults
#   dependencies:
#     - conda-forge::perl
#     - bioconda::perl-math-random
#     - bioconda::perl-inline-c
#     - bioconda::perl-parse-recdescent
#     - conda-forge::numpy
#     - bioconda::dwgsim
#     - alienzj::msort
RUN mkdir -p /conda-envs/6b5e3ed520a8389d84dfb0a12d5f07d4
COPY .harpy_envs/simulations.yaml /conda-envs/6b5e3ed520a8389d84dfb0a12d5f07d4/environment.yaml

# Conda environment:
#   source: .harpy_envs/snp.yaml
#   prefix: /conda-envs/b6f3d64eb4a83d06bc8a2731a1f32414
#   name: snp
#   channels:
#     - bioconda
#     - conda-forge
#     - defaults
#   dependencies:
#     - bioconda::bcftools=1.20
#     - bioconda::freebayes=1.3.6
RUN mkdir -p /conda-envs/b6f3d64eb4a83d06bc8a2731a1f32414
COPY .harpy_envs/snp.yaml /conda-envs/b6f3d64eb4a83d06bc8a2731a1f32414/environment.yaml

# Conda environment:
#   source: .harpy_envs/sv.yaml
#   prefix: /conda-envs/6f7bf3772ef42cc35bb95f8f02843984
#   name: sv
#   channels:
#     - bioconda
#     - conda-forge
#     - defaults
#   dependencies:
#     - bioconda::leviathan
#     - bioconda::naibr-plus
RUN mkdir -p /conda-envs/6f7bf3772ef42cc35bb95f8f02843984
COPY .harpy_envs/sv.yaml /conda-envs/6f7bf3772ef42cc35bb95f8f02843984/environment.yaml

# Step 2: Generate conda environments

RUN mamba env create --prefix /conda-envs/4538a63544ffaad3a6e49e8d25bfc6cb --file /conda-envs/4538a63544ffaad3a6e49e8d25bfc6cb/environment.yaml && \
    mamba env create --prefix /conda-envs/4914663d751c6d6b4e3b0b349fc0b652 --file /conda-envs/4914663d751c6d6b4e3b0b349fc0b652/environment.yaml && \
    mamba env create --prefix /conda-envs/6e3a7dc7152f45301056861193dd3db5 --file /conda-envs/6e3a7dc7152f45301056861193dd3db5/environment.yaml && \
    mamba env create --prefix /conda-envs/f769d2b509f06bc58137776df3f21936 --file /conda-envs/f769d2b509f06bc58137776df3f21936/environment.yaml && \
    mamba env create --prefix /conda-envs/6b5e3ed520a8389d84dfb0a12d5f07d4 --file /conda-envs/6b5e3ed520a8389d84dfb0a12d5f07d4/environment.yaml && \
    mamba env create --prefix /conda-envs/b6f3d64eb4a83d06bc8a2731a1f32414 --file /conda-envs/b6f3d64eb4a83d06bc8a2731a1f32414/environment.yaml && \
    mamba env create --prefix /conda-envs/6f7bf3772ef42cc35bb95f8f02843984 --file /conda-envs/6f7bf3772ef42cc35bb95f8f02843984/environment.yaml && \
    mamba clean --all -y

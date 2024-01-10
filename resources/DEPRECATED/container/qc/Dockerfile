# HARPY workflows

FROM mambaorg/micromamba:alpine AS qc
RUN micromamba install --yes --name base -c bioconda -c conda-forge \
    falco \ 
    fastp \ 
    multiqc && \ 
    micromamba clean --all --yes && \
    find -name '*.a' -delete && \
    rm -rf /env/conda-meta && \
    rm -rf /env/include && \
    find -name '__pycache__' -type d -exec rm -rf '{}' '+' && \
    rm -rf /env/lib/python3.9/site-packages/pip /env/lib/python3.9/idlelib /env/lib/python3.9/ensurepip \
        /env/lib/libasan.so.5.0.0 \
        /env/lib/libtsan.so.0.0.0 \
        /env/lib/liblsan.so.0.0.0 \
        /env/lib/libubsan.so.1.0.0 \
        /env/bin/x86_64-conda-linux-gnu-ld \
        /env/bin/sqlite3 \
        /env/bin/openssl \
        /env/share/terminfo

#FROM mambaorg/micromamba AS qc
#RUN micromamba install --yes --name base -c bioconda -c conda-forge \
#    falco \ 
#    fastp \ 
#    multiqc \ 
#    r-base \ 
#    r-dt \ 
#    r-flexdashboard \ 
#    r-ggplot2 \ 
#    r-plotly \ 
#    r-tidyr && \
#    micromamba clean --all --yes
#ARG MAMBA_DOCKERFILE_ACTIVATE=1
#COPY --chown=$MAMBA_USER:$MAMBA_USER workflow/scripts/ workflow/report/ $CONDA_PREFIX/bin/

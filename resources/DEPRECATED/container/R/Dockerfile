FROM mambaorg/micromamba:alpine
RUN micromamba install --yes --name base -c conda-forge \
    r-flexdashboard \ 
    r-ggplot2       \ 
    r-plotly        \ 
    r-tidyr         \
    r-dt &&         \ 
    micromamba clean --all --yes
#ARG MAMBA_DOCKERFILE_ACTIVATE=1
#COPY --chown=$MAMBA_USER:$MAMBA_USER workflow/scripts/ workflow/report/ $CONDA_PREFIX/bin/

from harpy.common.conda import CONDA_ENVS
from harpy import __version__

localrules: all, conda_env

rule all:
    input:
        collect("{conda}.env", conda = config.get("envs", list(CONDA_ENVS.keys())))

rule conda_env:
    output: "{conda}.env"
    conda: "envs/{conda}.yaml"
    container: "docker://pdimens/harpy:{conda}_" + __version__
    shell: "touch {output}"

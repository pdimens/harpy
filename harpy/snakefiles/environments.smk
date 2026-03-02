from harpy.common.conda import CONDA_ENVS

rule all:
    input:
        collect("{conda}.env", conda = config.get("envs", list(CONDA_ENVS.keys())))

rule conda_env:
    output: "{conda}.env"
    conda: "envs/{conda}.yaml"
    container: "docker://pdimens/harpy:{conda}_" + f"{VERSION}"
    shell: "touch {output}"

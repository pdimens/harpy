
if config.get("envs", None):
    out_envs = config["envs"]
else:
    out_envs =  ["align", "assembly", "metassembly", "phase", "qc", "r", "simulations", "stitch", "variants"] 

if config.get("spades", None):
    out_envs.append("spades")

rule all:
    input:
        collect("{conda}.env", conda = out_envs)

rule conda_env:
    output: "{conda}.env"
    container: "docker://pdimens/harpy:latest"
    conda: "envs/{conda}.yaml"
    shell: "touch {output}"

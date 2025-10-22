
out_envs = config.get("envs", ["align", "assembly", "metassembly", "phase", "qc", "report", "simulations", "stitch", "variants"])

rule all:
    input:
        collect("{conda}.env", conda = out_envs)

rule conda_env:
    output: "{conda}.env"
    conda: "envs/{conda}.yaml"
    shell: "touch {output}"

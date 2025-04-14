import os

out_envs =  ["align", "assembly", "metassembly", "phase", "qc", "r", "simulations", "stitch", "variants"]
if config.get("spades", None):
    out_envs.append("spades")
    outdir = os.path.join(os.getcwd(), "localenv/")
    else:
    # spades isn't added b/c it has a post-setup script, i.e. incompatible with containerization
        outdir = ""

rule all:
    input:
        collect(outdir + "{conda}.env", conda = out_envs)

rule conda_env:
    output: outdir + "{conda}.env"
    conda: "envs/{{conda}}.yaml"
    shell: "touch {output}"

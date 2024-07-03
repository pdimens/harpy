import os
import shutil

onsuccess:
    shutil.rmtree(f'.harpy_envs', ignore_errors=True)

rule all:
    input:
        collect("{conda}.env", conda = ["align", "phase", "qc", "r","simulations", "snp","sv"])

rule qc:
    output: "qc.env"
    conda: os.getcwd() + "/.harpy_envs/qc.yaml"
    shell: "touch {output}"

rule align:
    output: "align.env"
    conda: os.getcwd() + "/.harpy_envs/align.yaml"
    shell: "touch {output}"

rule snp:
    output: "snp.env"
    conda: os.getcwd() + "/.harpy_envs/snp.yaml"
    shell: "touch {output}"

rule sv:
    output: "sv.env"
    conda: os.getcwd() + "/.harpy_envs/sv.yaml"
    shell: "touch {output}"

rule phase:
    output: "phase.env"
    conda: os.getcwd() + "/.harpy_envs/phase.yaml"
    shell: "touch {output}"

rule r:
    output: "r.env"
    conda: os.getcwd() + "/.harpy_envs/r.yaml"
    shell: "touch {output}"

rule stitch:
    output: "stitch.env"
    conda: os.getcwd() + "/.harpy_envs/stitch.yaml"
    shell: "touch {output}"

rule simulations:
    output: "simulations.env"
    conda: os.getcwd() + "/.harpy_envs/simulations.yaml"
    shell: "touch {output}"
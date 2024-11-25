import os
import shutil

envdir = os.path.join(os.getcwd(), "container/workflow/envs")
# spades isn't added b/c it has a post-setup script
rule all:
    input:
        collect("{conda}.env", conda = ["align", "assembly", "metassembly", "phase", "qc", "r", "simulations", "stitch", "variants"])

rule qc:
    output: "qc.env"
    conda: f"{envdir}/qc.yaml"
    shell: "touch {output}"

rule align:
    output: "align.env"
    conda: f"{envdir}/align.yaml"
    shell: "touch {output}"

rule variants:
    output: "variants.env"
    conda: f"{envdir}/variants.yaml"
    shell: "touch {output}"

rule phase:
    output: "phase.env"
    conda: f"{envdir}/phase.yaml"
    shell: "touch {output}"

rule r:
    output: "r.env"
    conda: f"{envdir}/r.yaml"
    shell: "touch {output}"

rule stitch:
    output: "stitch.env"
    conda: f"{envdir}/stitch.yaml"
    shell: "touch {output}"

rule simulations:
    output: "simulations.env"
    conda: f"{envdir}/simulations.yaml"
    shell: "touch {output}"

rule assembly:
    output: "assembly.env"
    conda: f"{envdir}/assembly.yaml"
    shell: "touch {output}"

rule metassembly:
    output: "metassembly.env"
    conda: f"{envdir}/metassembly.yaml"
    shell: "touch {output}"
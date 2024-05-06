import os
import shutil

condachannels = ["bioconda","conda-forge","defaults"]
environ = {
    "qc" : ["bioconda::falco", "bioconda::fastp", "bioconda::multiqc", "bioconda::pysam=0.22"],
    "align": ["bioconda::bwa", "bioconda::ema","conda-forge::icu","conda-forge::libzlib", "bioconda::minimap2", "bioconda::samtools=1.20", "bioconda::seqtk", "bioconda::tabix", "conda-forge::xz"],
    "snp": ["bioconda::bcftools=1.20", "bioconda::freebayes=1.3.6"],
    "sv": ["bioconda::leviathan", "bioconda::naibr-plus"],
    "phase" : ["bioconda::hapcut2", "bioconda::whatshap"],
    "simulations" : ["conda-forge::perl", "bioconda::perl-math-random", "bioconda::perl-inline-c", "bioconda::perl-parse-recdescent", "conda-forge::numpy", "bioconda::dwgsim", "alienzj::msort"],
    "r" : ["conda-forge::r-xml2", "bioconda::bioconductor-complexheatmap", "conda-forge::r-highcharter", "conda-forge::r-circlize", "r::r-biocircos", "conda-forge::r-dt", "conda-forge::r-flexdashboard", "conda-forge::r-ggplot2", "conda-forge::r-ggridges", "conda-forge::r-plotly", "conda-forge::r-tidyr", "bioconda::r-stitch"]
}
os.makedirs(os.getcwd() + ".harpy_envs", exist_ok = True)
for i in environ:
    # overwrites existing
    with open(os.getcwd() + f"/.harpy_envs/{i}.yaml", "w") as yml:
        yml.write(f"name: {i}\n")
        yml.write("channels:\n  - ")
        yml.write("\n  - ".join(condachannels))
        yml.write("\ndependencies:\n  - ")
        yml.write("\n  - ".join(environ[i]) + "\n")

onsuccess:
    shutil.rmtree(f'.harpy_envs', ignore_errors=True)

rule all:
    input:
        collect("{conda}.env", conda = ["qc","align","snp","sv","phase","r","simulations"])

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

rule r_env:
    output: "r.env"
    conda: os.getcwd() + "/.harpy_envs/r.yaml"
    shell: "touch {output}"

rule simulations:
    output: "simulations.env"
    conda: os.getcwd() + "/.harpy_envs/simulations.yaml"
    shell: "touch {output}"
maxlen = config["maxlen"]
fpextra = config["extra"]
seq_dir = config["seq_directory"]
fqext = config["fqext"]
samplenames = config["samplenames"]

rule trim_fastp:
    input:
        fw = seq_dir + "/{sample}" + fqext[0],
        rv = seq_dir + "/{sample}" + fqext[1],
    output:
        fw = "Trimming/{sample}.R1.fq.gz",
        rv = "Trimming/{sample}.R2.fq.gz"
    log:
        html = "Trimming/logs/{sample}.html",
        json = "Trimming/logs/json/{sample}.json",
        err = "Trimming/logs/err/{sample}.log"
    threads: 2
    message: "Removing adapters + quality trimming: {wildcards.sample}"
    wildcard_constraints:
		sample = "[a-zA-Z0-9_-]*"
    params:
        maxlen = f"--max_len1 {maxlen}",
        extra = fpextra
    shell:
		"""
        fastp --trim_poly_g --cut_right --detect_adapter_for_pe {params} --thread {threads} -i {input.fw} -I {input.rv} -o {output.fw} -O {output.rv} -h {log.html} -j {log.json} 2> {log.err}
        """

rule all:
    input: expand("Trimming/{sample}{ext}", sample = samplenames, ext = fqext)
    shell: "multiqc Trimming/logs/"
    default_target: True
maxlen = config["maxlen"]
fpextra = config["extra"]
seq_dir = config["seq_directory"]
fqext = config["fqext"]
samplenames = config["samplenames"]

rule trim_fastp:
    input:
        fw = seq_dir + "/{sample}" + fqext[0],
        rv = seq_dir + "/{sample}" + fqext[1]
    output:
        fw = "Trimming/{sample}.R1.fq.gz",
        rv = "Trimming/{sample}.R2.fq.gz"
    log:
        html = "Trimming/logs/{sample}.html",
        json = "Trimming/logs/json/{sample}.json",
        err = "Trimming/logs/err/{sample}.log"
    threads: 2
	benchmark: "Benchmark/Trimming/{sample}.txt"
    wildcard_constraints: 
        sample = "[a-zA-Z0-9_-]*"
    message: "Removing adapters + quality trimming: {wildcards.sample}"
    params:
        maxlen = f"--max_len1 {maxlen}",
        extra = fpextra
    shell: "fastp --trim_poly_g --cut_right --detect_adapter_for_pe {params} --thread {threads} -i {input.fw} -I {input.rv} -o {output.fw} -O {output.rv} -h {log.html} -j {log.json} 2> {log.err}"

rule reports:
    input: expand("Trimming/{sample}{ext}", sample = samplenames, ext = [".R1.fq.gz", ".R2.fq.gz"])
	benchmark: "Benchmark/Trimming/report.txt"
    default_target: True
    shell: "multiqc Trimming/logs/"
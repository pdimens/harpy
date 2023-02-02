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
		serr = "Trimming/logs/err/{sample}.log"
	benchmark: "Benchmark/Trimming/{sample}.txt"
	message: "Removing adapters + quality trimming: {wildcards.sample}"
	wildcard_constraints: 
		sample = "[a-zA-Z0-9_-]*"
	threads: 2
	params:
		maxlen = f"--max_len1 {maxlen}",
		extra = fpextra
	shell: "fastp --trim_poly_g --cut_right --detect_adapter_for_pe {params} --thread {threads} -i {input.fw} -I {input.rv} -o {output.fw} -O {output.rv} -h {log.html} -j {log.json} 2> {log.serr}"

rule reports:
	input: expand("Trimming/{sample}{ext}", sample = samplenames, ext = [".R1.fq.gz", ".R2.fq.gz"])
	benchmark: "Benchmark/Trimming/report.txt"
	default_target: True
	shell: "multiqc Trimming/logs/"
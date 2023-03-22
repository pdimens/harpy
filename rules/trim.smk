maxlen = config["maxlen"]
extra = config.get("extra", "") 
seq_dir = config["seq_directory"]
fqext = config["fqext"]
samplenames = config["samplenames"]

rule trimFastp:
	input:
		fw = seq_dir + "/{sample}" + fqext[0],
		rv = seq_dir + "/{sample}" + fqext[1]
	output:
		fw = "Trimming/{sample}.R1.fq.gz",
		rv = "Trimming/{sample}.R2.fq.gz",
		json = "Trimming/logs/json/{sample}.fastp.json"
	log:
		html = "Trimming/logs/html/{sample}.html",
		serr = "Trimming/logs/err/{sample}.log"
	benchmark: "Benchmark/Trimming/{sample}.txt"
	message: "Removing adapters + quality trimming: {wildcards.sample}"
	wildcard_constraints: 
		sample = "[a-zA-Z0-9_-]*"
	threads: 2
	params:
		maxlen = f"--max_len1 {maxlen}",
		extra = extra
	shell: 
		"fastp --trim_poly_g --cut_right --detect_adapter_for_pe {params} --thread {threads} -i {input.fw} -I {input.rv} -o {output.fw} -O {output.rv} -h {log.html} -j {output.json} 2> {log.serr}"

rule createReport:
	input: expand("Trimming/logs/json/{sample}.fastp.json", sample = samplenames)
	output: "Trimming/logs/trim.report.html"
	shell: 
		"multiqc Trimming/logs/json -m fastp --force --filename {output} --quiet --no-data-dir 2>/dev/null"

rule trimCheck:
	input: 
		expand("Trimming/{sample}{ext}", sample = samplenames, ext = [".R1.fq.gz", ".R2.fq.gz"]),
		"Trimming/logs/trim.report.html"
	default_target: True
	message: "Trimming is complete!"

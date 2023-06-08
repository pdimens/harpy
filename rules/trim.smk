maxlen 		= config["maxlen"]
extra 		= config.get("extra", "") 
seq_dir 	= config["seq_directory"]
fqext 		= config["fqext"]
samplenames = config["samplenames"]

rule trimFastp:
	input:
		fw   = seq_dir + "/{sample}" + fqext[0],
		rv   = seq_dir + "/{sample}" + fqext[1]
	output:
		fw   = "Trim/{sample}.R1.fq.gz",
		rv   = "Trim/{sample}.R2.fq.gz",
		json = "Trim/logs/json/{sample}.fastp.json"
	log:
		html = "Trim/logs/html/{sample}.html",
		serr = "Trim/logs/err/{sample}.log"
	benchmark:
		"Benchmark/Trim/{sample}.txt"
	message:
		"Removing adapters + quality trimming: {wildcards.sample}"
	wildcard_constraints: 
		sample = "[a-zA-Z0-9_-]*"
	threads: 2
	params:
		maxlen = f"--max_len1 {maxlen}",
		extra = extra
	shell: 
		"fastp --trim_poly_g --cut_right --detect_adapter_for_pe {params} --thread {threads} -i {input.fw} -I {input.rv} -o {output.fw} -O {output.rv} -h {log.html} -j {output.json} 2> {log.serr}"

rule createReport:
	input: 
		json = expand("Trim/logs/json/{sample}.fastp.json", sample = samplenames),
		fr   = expand("Trim/{sample}.R1.fq.gz", sample = samplenames),
		rv   = expand("Trim/{sample}.R2.fq.gz", sample = samplenames)
	output:
		"Trim/logs/trim.report.html"
	message:
		"Sequencing quality filtering and trimming is complete!"
	default_target: True
	shell: 
		"multiqc Trim/logs/json -m fastp --force --filename {output} --quiet --no-data-dir 2>/dev/null"
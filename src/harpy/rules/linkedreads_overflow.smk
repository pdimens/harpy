rule genome_faidx:
    input:
        outdir + "/workflow/input/hap.{hap}.fasta"
    output: 
        outdir + "/workflow/input/hap.{hap}.fasta.fai"
    message:
        "Indexing haplotype {wildcards.hap}"
    shell:
        "samtools faidx --fai-idx {output} {input}"

if not barcodes:
    rule download_barcodes:
        output:
            barcodefile
        message:
            "Downloading list of standard 10X barcodes"
        run:
            from urllib.request import urlretrieve
            _ = urlretrieve("https://raw.githubusercontent.com/aquaskyline/LRSIM/master/4M-with-alts-february-2016.txt", output[0])

rule create_molecules_hap:
    input:
        outdir + "/workflow/input/hap.{hap}.fasta"
    output:
        temp(multiext(outdir + "/dwgsim.{hap}.12", ".bwa.read1.fastq.gz" ,".bwa.read2.fastq.gz", ".mutations.txt", ".mutations.vcf"))
    log:
        outdir + "/logs/dwgsim.hap.{hap}.log"
    params:
        readpairs = config["read_pairs"] * 500000,
        outerdist = config["outer_distance"],
        distsd = config["distance_sd"],
        prefix = lambda wc: outdir + "/dwgsim." + wc.get("hap") + ".12"
    conda:
        os.getcwd() + "/.harpy_envs/simulations.yaml"
    message:
        "Creating reads from {input}"
    shell:
        """
        dwgsim -N {params.readpairs} -e 0.0001,0.0016 -E 0.0001,0.0016 -d {params.outerdist} -s {params.distsd} -1 135 -2 151 -H -y 0 -S 0 -c 0 -R 0 -r 0 -F 0 -o 1 -m /dev/null {input} {params.prefix} 2> {log}
        """

rule interleave_dwgsim_output:
    default_target: True
    input:
        collect(outdir + "/dwgsim.{{hap}}.12.bwa.read{rd}.fastq.gz", rd = [1,2]) 
    output:
        outdir + "/dwgsim.{hap}.12.fastq"
    message:
        "Decompressing DWGSIM haplotype {wildcards.hap} output"
    shell:
        "seqtk mergepe {input} > {output}"

rule lrsim:
    default_target: True
    input:
        hap1 = f"{outdir}/dwgsim.0.12.fastq",
        hap2 = f"{outdir}/dwgsim.1.12.fastq",
        fai1 = outdir + "/workflow/input/hap.0.fasta.fai",
        fai2 = outdir + "/workflow/input/hap.1.fasta.fai",
        barcodes = barcodefile
    output:
        collect(outdir + "/sim.{hap}.{ext}", hap = [0,1], ext = ["fp", "manifest"])
    log:
        f"{outdir}/logs/LRSIM.log"
    params:
        lrsim = f"{outdir}/src/harpy/scripts/LRSIMharpy.pl",
        static = "-o 1 -d 2 -u 4",
        proj_dir = f"-r {outdir}",
        prefix = f"-p {outdir}/sim",
        outdist  = f"""-i {config["outer_distance"]}""",
        dist_sd  = f"""-s {config["distance_sd"]}""",
        n_pairs  = f"""-x {config["read_pairs"]}""",
        mol_len  = f"""-f {config["molecule_length"]}""",
        parts    = f"""-t {config["partitions"]}""",
        mols_per = f"""-m {config["molecules_per_partition"]}"""
    threads:
        workflow.cores
    conda:
        os.getcwd() + "/.harpy_envs/simulations.yaml"
    message:
        "Running LRSIM to generate linked reads from\nhaplotype 1: {input.hap1}\nhaplotype 2: {input.hap2}" 
    shell: 
        "perl {params} -g {input.hap1},{input.hap2} -b {input.barcodes} -z {threads} 2> {log}"

rule sort_manifest:
    input:
        outdir + "/sim.{hap}.manifest"
    output:
        outdir + "/sim.{hap}.sort.manifest"
    shell:
        "msort -kn1 {input} > {output}"

rule extract_reads:
    input:
        manifest = outdir + "/sim.{hap}.sort.manifest",
        dwg_hap = outdir + "/dwgsim/hap{hap}.12.fastq"
    output:
        outdir + "/sim_hap{hap}_R1_001.fastq.gz",
        outdir + "/sim_hap{hap}_R2_001.fastq.gz"
    params:
        lambda wc: outdir + "/sim_hap" + wc.get("hap")
    shell:
        "extractReads {input} {output}"

rule convert_haplotag:
    input:
        fw = outdir + "/sim_hap{hap}_R1_001.fastq.gz",
        rv = outdir + "/sim_hap{hap}_R2_001.fastq.gz",
        barcodes = barcodefile
    output:
        fw = outdir + "/sim_hap{hap}_haplotag.R1.fq.gz",
        rv = outdir + "/sim_hap{hap}_haplotag.R2.fq.gz"
    log:
        outdir + "/workflow/10XtoHaplotag_{hap}.txt" 
    message:
        "Converting 10X barcodes to haplotag format"
    shell:
        "10xtoHaplotag.py -f {input.fw} -r {input.rv} -b {input.barcodes} -l {log}"

rule log_workflow:
    input:
        collect(outdir + "/hap{hap}_haplotag.R{fw}.fq.gz", hap = [1,2], fw = [1,2])
    output:
        outdir + "/workflow/simulate.reads.workflow.summary"
    params:
        static = "-o 1 -d 2 -u 4",
        proj_dir = f"-r {outdir}",
        prefix = f"-p {outdir}/sim",
        outdist  = f"""-i {config["outer_distance"]}""",
        dist_sd  = f"""-s {config["distance_sd"]}""",
        n_pairs  = f"""-x {config["read_pairs"]}""",
        mol_len  = f"""-f {config["molecule_length"]}""",
        parts    = f"""-t {config["partitions"]}""",
        mols_per = f"""-m {config["molecules_per_partition"]}"""
    message:
        "Summarizing the workflow: {output}"
    run:
        with open(output[0], "w") as f:
            _ = f.write("The harpy simulate reads module ran using these parameters:\n\n")
            _ = f.write(f"Genome haplotype 1: {gen_hap1}\n")
            _ = f.write(f"Genome haplotype 2: {gen_hap2}\n")
            _ = f.write(f"Barcode file: {barcodefile}\n")
            _ = f.write("LRSIM was started from step 3 (-u 3) with these parameters:\n")
            _ = f.write("    " + f"LRSIMharpy.pl -g genome1,genome2 " + " ".join(params) + "\n")
            _ = f.write("10X style barcodes were converted in haplotag BX:Z tags using:\n")
            _ = f.write("    " + f"10xtoHaplotag.py")
            _ = f.write("\nThe Snakemake workflow was called via command line:\n")
            _ = f.write("    " + str(config["workflow_call"]) + "\n")
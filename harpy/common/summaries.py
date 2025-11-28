"""Basic functions to write the workflow summaries"""
import os
import re

class Summary:
    def __init__(self, config: dict):
        self.config: dict = config

    def get_text(self) -> str:
        return self.__getattribute__(self.config["Workflow"]["name"])()

    def align_bwa(self) -> str:
        ignore_bx = self.config["Workflow"]["linkedreads"]["type"] == "none"
        is_standardized = self.config["Workflow"]["linkedreads"]["standardized"]
        keep_unmapped = self.config["Parameters"]["keep-unmapped"]
        extra 		= self.config["Parameters"].get("extra", "") 
        genomefile 	= self.config["Inputs"]["reference"]
        quality = self.config["Parameters"]["min-map-quality"]

        unmapped = "" if keep_unmapped else "-F 4"
        bx_mode = "--barcode-tag BX" if not ignore_bx else ""
        bwa_static = "-C -v 2" if is_standardized else "-v 2"
        extra   = extra

        summary = ["The harpy align bwa workflow ran using these parameters:"]
        summary.append(f"The provided genome: {genomefile}")
        align = "Sequences were aligned with BWA using:\n"
        align += f'\tbwa mem {bwa_static} {extra} -R "@RG\\tID:SAMPLE\\tSM:SAMPLE" genome forward_reads reverse_reads |\n'
        align += f"\tsamtools view -h {unmapped} -q {quality}"
        summary.append(align)
        standardization = "Barcodes were standardized in the aligments using:\n"
        standardization += "\tstandardize_barcodes_sam > {output} < {input}"
        summary.append(standardization)
        duplicates = "Duplicates in the alignments were marked following:\n"
        duplicates += "\tsamtools collate |\n"
        duplicates += "\tsamtools fixmate |\n"
        duplicates += f"\tsamtools sort -T SAMPLE --reference {genomefile} -m 2000M |\n"
        duplicates += f"\tsamtools markdup -S {bx_mode} -d 100 (2500 for novaseq)"
        summary.append(duplicates)
        sm = "The Snakemake command invoked:\n"
        sm += f"\t{self.config["Workflow"]['snakemake']['relative']}"
        summary.append(sm)
        return "\n\n".join(summary)

    def align_strobe(self) -> str:
        genomefile 	= self.config["Inputs"]["reference"]
        ignore_bx = self.config["Workflow"]["linkedreads"]["type"] == "none"
        is_standardized = self.config["Workflow"]["linkedreads"]["standardized"]
        keep_unmapped = self.config["Parameters"]["keep-unmapped"]

        quality = self.config["Parameters"]["min-map-quality"]
        unmapped_strobe = "" if keep_unmapped else "-U"
        unmapped = "" if keep_unmapped else "-F 4"
        bx_mode = "--barcode-tag BX" if not ignore_bx else ""
        static = "-C" if is_standardized else ""
        extra = self.config["Parameters"].get("extra", "") 

        summary = ["The harpy align strobe workflow ran using these parameters:"]
        summary.append(f"The provided genome: {genomefile}")
        align = "Sequences were aligned with strobealign using:\n"
        align += f"\tstrobealign {unmapped_strobe} {static} --rg-id=SAMPLE --rg=SM:SAMPLE {extra} genome reads.F.fq reads.R.fq |\n"
        align += f"\t\tsamtools view -h {unmapped} -q {quality}"
        summary.append(align)
        standardization = "Barcodes were standardized in the aligments using:\n"
        standardization += "\tstandardize_barcodes_sam > {output} < {input}"
        summary.append(standardization)
        duplicates = "Duplicates in the alignments were marked following:\n"
        duplicates += "\tsamtools collate |\n"
        duplicates += "\tsamtools fixmate |\n"
        duplicates += f"\tsamtools sort -T SAMPLE --reference {genomefile} -m 2000M |\n"
        duplicates += f"\tsamtools markdup -S {bx_mode} -d 100 (2500 for novaseq)"
        summary.append(duplicates)
        sm = "The Snakemake command invoked:\n"
        sm += f"\t{self.config["Workflow"]['snakemake']['relative']}"
        summary.append(sm)
        return "\n\n".join(summary)

    def assembly(self) -> str:
        # SPADES
        max_mem      = self.config["Parameters"]["spades"]["max-memory"]
        k_param      = self.config["Parameters"]["spades"]["k"]
        spades_extra = self.config["Parameters"]["spades"].get("extra", "")
        # ARCS
        mapq       = self.config["Parameters"]["tigmint"]["minimum-mapping-quality"]
        mismatch   = self.config["Parameters"]["tigmint"]["mismatch"]
        mol_dist   = self.config["Parameters"]["tigmint"]["molecule-distance"]
        mol_len    = self.config["Parameters"]["tigmint"]["molecule-length"]
        span       = self.config["Parameters"]["tigmint"]["span"]
        min_align  = self.config["Parameters"]["arcs"]["minimum-aligned-reads"]
        min_contig = self.config["Parameters"]["arcs"]["minimum-contig-length"]
        seq_id     = self.config["Parameters"]["arcs"]["minimum-sequence-identity"]
        arcs_extra = self.config["Parameters"]["arcs"].get("extra", "")
        links      = self.config["Parameters"]["links"]["minimum-links"]
        k_param = k_param
        max_mem = max_mem // 1000
        spades_extra = spades_extra
        params = [
            "-C scaffold",
            "-j THREADS",
            "draft=spades",
            "reads=interleaved",
            "t=THREADS",
            f"mapq={mapq}",
            f"nm={mismatch}",
            f"dist={mol_dist}",
            f"minsize={mol_len}",
            f"span={span}",
            f"c={min_align}",
            f"z={min_contig}",
            f"s={seq_id}",
            f"l={links}",
            arcs_extra
        ]
        summary = ["The harpy assemble workflow ran using these parameters:"]
        spades = "Reads were assembled using cloudspades:\n"
        spades += f"\tspades.py -t THREADS -m {max_mem} --gemcode1-1 FQ1 --gemcode1-2 FQ2 --isolate -k {k_param} {spades_extra}"
        summary.append(spades)
        arcs = "The draft assembly was error corrected and scaffolded with Tigmint/ARCS/LINKS:\n"
        arcs += f"\tarcs-make arcs-tigmint {' '.join(params[3:])}"
        summary.append(arcs)
        sm = "The Snakemake command invoked:\n"
        sm += f"\t{self.config["Workflow"]['snakemake']['relative']}"
        summary.append(sm)
        return "\n\n".join(summary)

    def deconvolve(self) -> str:
        kmer_length = self.config["Parameters"]["kmer-length"]
        window_size = self.config["Parameters"]["window-size"]
        density 	= self.config["Parameters"]["density"] 
        dropout     = self.config["Parameters"]["dropout"]

        summary = ["The harpy deconvolve workflow ran using these parameters:"]
        interleave = "fastq files were interleaved with seqtk:\n"
        interleave += "\tseqtk mergepe forward.fq reverse.fq"
        summary.append(interleave)
        deconv = "Deconvolution occurred using QuickDeconvolution:\n"
        deconv += f"\tQuickDeconvolution -t threads -i infile.fq -o output.fq -k {kmer_length} -w {window_size} -d {density} -a {dropout}"
        summary.append(deconv)
        recover = "The interleaved output was split back into forward and reverse reads with seqtk:\n"
        recover += "\tseqtk seq -1 interleaved.fq | gzip > file.R1.fq.gz\n"
        recover += "\tseqtk seq -2 interleaved.fq | gzip > file.R2.fq.gz"
        summary.append(recover)
        sm = "The Snakemake command invoked:\n"
        sm += f"\t{self.config["Workflow"]['snakemake']['relative']}"
        summary.append(sm)
        return "\n\n".join(summary)

    def demultiplex_meier2021(self) -> str:
        schemafile = self.config["Inputs"]["demultiplex_schema"]
        qxrx = self.config["Parameters"]["qx-rx"]
        unknown_samples = self.config["Parameters"]["samples"]
        unknown_barcodes = self.config["Parameters"]["barcodes"]

        R1 = self.config["Inputs"]["R1"],
        R2 = self.config["Inputs"]["R2"],
        I1 = self.config["Inputs"]["I1"],
        I2 = self.config["Inputs"]["I2"],
        outdir = f"--samples {os.getcwd()}",
        qxrx = "--rx --qx" if qxrx else "",
        unknown_barcodes = "--undetermined-barcodes _unknown_barcodes" if unknown_barcodes else "",
        unknown_samples = "--undetermined-samples _unknown_samples" if unknown_samples else ""

        summary = ["The harpy demultiplex workflow ran using these parameters:"]
        summary.append("Linked Read Barcode Design: Meier et al. 2021")
        inputs = "The multiplexed input files:\n"
        inputs += f"\tread 1: {R1}\n"
        inputs += f"\tread 2: {R2}\n"
        inputs += f"\tindex 1: {I1}\n"
        inputs += f"\tindex 2: {I2}\n"
        inputs += f"Sample demultiplexing schema: {schemafile}"
        summary.append(inputs)
        demux = "Samples were demultiplexed using:\n"
        demux += f"\tdmox --R1 --R2 --I1 --I2 {outdir} {qxrx} {unknown_barcodes} {unknown_samples}"
        summary.append(demux)
        qc = "QC checks were performed on demultiplexed FASTQ files using:\n"
        qc += "\tfalco -skip-report -skip-summary -data-filename output input.fq.gz"
        summary.append(qc)
        sm = "The Snakemake command invoked:\n"
        sm += f"\t{self.config["Workflow"]['snakemake']['relative']}"
        summary.append(sm)
        return "\n\n".join(summary)

    def impute(self) -> str:
        region = self.config["Parameters"].get("region", None)
        if region:
            _,positions = region.split(":")
            startpos,endpos,buffer = [int(i) for i in positions.split("-")]
            regiontext = f"\t\tregionStart = {startpos},\n"
            regiontext += f"\t\tregionEnd = {endpos},\n"
            regiontext += f"\t\tbuffer = {buffer},\n"
        else:
            regiontext = ""
        paramfiletext = "\t".join(open(self.config["Inputs"]["parameters"], "r").readlines())
        summary = ["The harpy impute workflow ran using these parameters:"]
        summary.append(f"The provided variant file: {self.config['Inputs']['vcf']}")
        preproc = "Preprocessing was performed with:\n"
        preproc += "\tbcftools view -M2 -v snps --regions CONTIG INFILE |\n"
        preproc += """\tbcftools query -i '(STRLEN(REF)==1) & (STRLEN(ALT[0])==1) & (REF!="N")' -f '%CHROM\\t%POS\\t%REF\\t%ALT\\n'"""
        summary.append(preproc)
        stitchparam = f"The STITCH parameter file: {self.config['Inputs']['parameters']}\n"
        stitchparam += f"\t{paramfiletext}"
        summary.append(stitchparam)
        stitch = "Within R, STITCH was invoked with the following parameters:\n"
        stitch += "\tSTITCH(\n"
        stitch += "\t\tmethod = model,\n"
        stitch += "\t\tposfile = posfile,\n"
        stitch += "\t\tbamlist = bamlist,\n"
        stitch += "\t\tnCores = ncores,\n"
        stitch += "\t\tnGen = ngen,\n"
        stitch += "\t\tchr = chr,\n"
        stitch += regiontext
        stitch += "\t\tK = k,\n"
        stitch += "\t\tS = s,\n"
        stitch += "\t\tuse_bx_tag = usebx,\n"
        stitch += "\t\tbxTagUpperLimit = bxlimit,\n"
        stitch += "\t\tniterations = 40,\n"
        stitch += "\t\tswitchModelIteration = 39,\n"
        stitch += "\t\tsplitReadIterations = NA,\n"
        if self.config["Parameters"]["grid-size"] > 1:
            stitch += f"\t\tgridWindowSize = {self.config["Parameters"]['grid-size']}\n"
        stitch += "\t\toutputdir = outdir,\n"
        stitch += "\t\toutput_filename = outfile\n\t)"
        stitchextra = "Additional STITCH parameters provided (overrides existing values above):\n"
        stitchextra += "\t" + self.config["Parameters"].get("extra", "None")
        summary.append(stitchextra)
        sm = "The Snakemake command invoked:\n"
        sm += f"\t{self.config["Workflow"]['snakemake']['relative']}"
        summary.append(sm)
        return "\n\n".join(summary)

    def metassembly(self) -> str:
        bx = self.config["Workflow"]["linkedreads"]["barcode-tag"]
        max_mem = self.config["Parameters"]["spades"]["max-memory"]
        k_param = self.config["Parameters"]["spades"]["k"]
        ignore_bx = self.config["Parameters"]["spades"]["ignore-barcodes"]
        extra = self.config["Parameters"]["spades"].get("extra", "")
        spadesdir = f"{'cloudspades' if not ignore_bx else 'spades'}_assembly"

        summary = ["The harpy metassembly workflow ran using these parameters:"]  
        bxsort = "FASTQ inputs were sorted by their linked-read barcodes:\n"
        bxsort += "\tsamtools import -T \"*\" FQ1 FQ2 |\n"
        bxsort += f"\tsamtools sort -O SAM -t {bx} |\n"  
        bxsort += "\tsamtools fastq -T \"*\" -1 FQ_out1 -2 FQ_out2"  
        summary.append(bxsort)
        bxappend = "Barcoded-sorted FASTQ files had \"-1\" appended to the barcode to make them Athena-compliant:\n"  
        bxappend += f"\tsed 's/{bx}:Z:[^[:space:]]*/&-1/g' FASTQ | bgzip > FASTQ_OUT"  
        summary.append(bxappend)
        if not ignore_bx:
            spades = "Reads were assembled using cloudspades:\n"
            spades += f"\tspades.py -t THREADS -m {max_mem} --gemcode1-1 FQ1 --gemcode1-2 FQ2 --meta -k {k_param} {extra}"
        else:
            spades = "Reads were assembled using spades:\n"
            spades += f"\tmetaspades.py -t THREADS -m {max_mem} -k {k_param} {extra} -1 FQ_1 -2 FQ2 -o {spadesdir}"
        summary.append(spades)
        align = "Original input FASTQ files were aligned to the metagenome using BWA:\n"
        align += "\tbwa mem -C -p spades.contigs FQ1 FQ2 | samtools sort -O bam -"
        summary.append(align)
        interleaved = "Barcode-sorted Athena-compliant sequences were interleaved with seqtk:\n"
        interleaved += "\tseqtk mergepe FQ1 FQ2 > INTERLEAVED.FQ"
        summary.append(interleaved)
        athena = "Athena ran with the config file Harpy built from the files created from the previous steps:\n"
        athena += "\tathena-meta --config athena.config"
        summary.append(athena)
        sm = "The Snakemake command invoked:\n"
        sm += f"\t{self.config["Workflow"]['snakemake']['relative']}"
        summary.append(sm)
        return "\n\n".join(summary)

    def phase_snp(self) -> str:
        bc_type           = self.config["Workflow"]["linkedreads"]["type"]
        pruning           = self.config["Parameters"]["prune"]
        map_qual          = self.config["Parameters"]["min-map-quality"]
        base_qual         = self.config["Parameters"]["min-base-quality"]
        molecule_distance = self.config["Parameters"]["distance-threshold"]
        extra             = self.config["Parameters"].get("extra", "") 
        variantfile       = self.config["Inputs"]["vcf"]
        invalid_regex = {
            "haplotagging" : "'$4 !~ /[ABCD]00/'",
            "stlfr" : "'$4 !~ /^0_|_0_|_0$/'",
            "tellseq": "'$4 !~ /N/'"
        }
        linkarg = "--10x 0" if bc_type == "none" else "--10x 1"
        indelarg   = "--indels 1 --ref reference.fasta" if self.config["Inputs"].get("reference", None) else ""
        hairs_params = f"{indelarg} {linkarg} --mmq {map_qual} --mbq {base_qual} --nf 1 --maxfragments 1500000"
        prune = f"--threshold {pruning}" if pruning > 0 else "--no_prune 1"

        summary = ["The harpy phase bam workflow ran using these parameters:"]
        summary.append(f"The provided variant file: {variantfile}")
        hetsplit = "The variant file was split by sample and filtered for heterozygous sites using:\n"
        hetsplit += "\tbcftools view -s SAMPLE | bcftools view -m 2 -M 2 -i \'GT=\"het\"\'"
        summary.append(hetsplit)
        phase = "Phasing was performed using the components of HapCut2:\n"
        phase += f"\textractHAIRS {hairs_params} --bam sample.bam --VCF sample.vcf --out sample.unlinked.frags\n"
        if bc_type != "none":
            phase += "\t awk " + invalid_regex.get(bc_type, "'$4 !~ /N/'") + " sample.unlinked.frags > sample.frags.filt"
            phase += f"\tLinkFragments.py --bam sample.bam --VCF sample.vcf --fragments sample.frags.filt --out sample.linked.frags -d {molecule_distance}\n"
            phase += f"\tHAPCUT2 --fragments sample.linked.frags --vcf sample.vcf --out sample.blocks --nf 1 --error_analysis_mode 1 --call_homozygous 1 --outvcf 1 {prune} {extra}\n"
        else:
            phase += f"\tHAPCUT2 --fragments sample.unlinked.frags --vcf sample.vcf --out sample.blocks --nf 1 --error_analysis_mode 1 --call_homozygous 1 --outvcf 1 {prune} {extra}\n"
        summary.append(phase)
        annot = "Variant annotation was performed using:\n"
        annot += "\tbcftools annotate -a sample.phased.vcf -c CHROM,POS,FMT/GT,FMT/PS,FMT/PQ,FMT/PD -m +HAPCUT\n"
        annot += "\tbcftools merge --output-type b samples.annot.bcf"
        summary.append(annot)
        sm = "The Snakemake command invoked:\n"
        sm += f"\t{self.config["Workflow"]['snakemake']['relative']}"
        summary.append(sm)
        return "\n\n".join(summary)

    def phase_bam(self) -> str:
        extra             = self.config["Parameters"].get("extra", "") 
        variantfile       = self.config["Inputs"]["vcf"]

        summary = ["The harpy phase snp workflow ran using these parameters:"]
        summary.append(f"The provided variant file: {variantfile}")
        validsplit = "The input alignments had their records filtered for valid barcodes:\n"
        validsplit += "\tdjinn filter-invalid --invalid filtered/sample sample.bam"
        summary.append(validsplit)
        phase = "Phasing was performed using whatshap:\n"
        phaseparam = "--linked-read-distance-cutoff {moldist} --tag-supplementary copy-primary --no-supplementary-strand-match --supplementary-distance {moldist} --ignore-read-groups --skip-missing-contigs"
        phase += f"\twhatshap haplotag --sample name --reference input.ref {phaseparam} {extra} input.vcf input.bam"
        summary.append(phase)
        concataln = "Invalid-barcode alignments were added back to the phased alignments using:\n"
        concataln += "\tsamtools merge sample.phased.bam sample.invalid.bam | samtools sort -"
        summary.append(concataln)
        sm = "The Snakemake command invoked:\n"
        sm += f"\t{self.config["Workflow"]['snakemake']['relative']}"
        summary.append(sm)
        return "\n\n".join(summary)

    def qc(self) -> str:
        minlen = f"--length_required {self.config["Parameters"]['min-len']}"
        maxlen = f"--max_len1 {self.config["Parameters"]['max-len']}"
        extra = self.config["Parameters"].get("extra", "") 
        trim_adapters = self.config["Parameters"].get("trim_adapters", None)
        if trim_adapters:
            trim_arg = "--detect_adapter_for_pe" if trim_adapters == "auto" else f"--adapter_fasta {trim_adapters}"
        else:
            trim_arg = "--disable_adapter_trimming"
        dedup = "-D" if self.config["Parameters"]["deduplicate"] else ""

        summary = ["The harpy qc workflow ran using these parameters:"]
        fastp = "fastp ran using:\n"
        fastp += "\tfastp --trim_poly_g --cut_right " + " ".join([minlen,maxlen,trim_arg,dedup,extra])
        summary.append(fastp)
        sm = "The Snakemake command invoked:\n"
        sm += f"\t{self.config["Workflow"]['snakemake']['relative']}"
        summary.append(sm)
        return "\n\n".join(summary)
    
    def snp_freebayes(self) -> str:
        ploidy 		= self.config["Parameters"]["ploidy"]
        extra 	    = self.config["Parameters"].get("extra", "") 
        regions_input = self.config["Inputs"]["regions"]
        genomefile 	= os.path.basename(self.config["Inputs"]["reference"])
        groupings 	= self.config["Inputs"].get("groupings", None)

        params = f"-p {ploidy} "
        params += f"--populations {groupings} " if groupings else ''
        params += extra

        summary = ["The harpy snp freebayes workflow ran using these parameters:"]
        summary.append(f"The provided reference genome: {genomefile}")
        summary.append(f"Genomic positions for which variants were called: {regions_input}")
        varcall = "The freebayes parameters:\n"
        varcall += f"\tfreebayes -f REFERENCE -L samples.list -r REGION {params} |\n"
        varcall += "\tbcftools sort -"
        summary.append(varcall)
        merged = "The variants identified in the intervals were merged into the final variant file using:\n"
        merged += "\tbcftools concat -f bcf.files -a --remove-duplicates"
        summary.append(merged)
        normalize = "The variants were normalized using:\n"
        normalize += "\tbcftools norm -m -both -d both -c w"
        summary.append(normalize)
        sm = "The Snakemake command invoked:\n"
        sm += f"\t{self.config["Workflow"]['snakemake']['relative']}"
        summary.append(sm)
        return "\n\n".join(summary)
    
    def snp_mpileup(self) -> str:
        mp_extra = self.config["Parameters"].get("extra", "")
        genomefile = self.config["Inputs"]["reference"]
        groupings = self.config["Inputs"].get("groupings", [])
        region_input = self.config["Inputs"]["regions"]
        ploidy = self.config["Parameters"]["ploidy"]
        params = f"--ploidy {ploidy} --populations "
        params += f"{groupings}" if groupings else "-"

        summary = ["The harpy snp freebayes workflow ran using these parameters:"]
        summary.append(f"The provided reference genome: {genomefile}")
        summary.append(f"Genomic positions for which variants were called: {region_input}")
        mpileup = "The mpileup parameters:\n"
        mpileup += f"\tbcftools mpileup --fasta-ref REFERENCE --region REGION --bam-list BAMS --annotate AD --output-type b {mp_extra}"
        summary.append(mpileup)
        bcfcall = "The bcftools call parameters:\n"
        bcfcall += f"\tbcftools call --multiallelic-caller {params} --variants-only --output-type b |\n"
        bcfcall += "\tbcftools sort -"
        summary.append(bcfcall)
        merged = "The variants identified in the intervals were merged into the final variant file using:\n"
        merged += "\tbcftools concat -f bcf.files -a --remove-duplicates"
        summary.append(merged)
        normalize = "The variants were normalized using:\n"
        normalize += "\tbcftools norm -m -both -d both -c w"
        summary.append(normalize)
        sm = "The Snakemake command invoked:\n"
        sm += f"\t{self.config["Workflow"]['snakemake']['relative']}"
        summary.append(sm)
        return "\n\n".join(summary)

    def validate_bam(self) -> str:
        lr_platform = self.config["Workflow"]["linkedreads"]["type"]

        summary = ["The harpy validate bam workflow ran using these parameters:"]
        valids = "Validations were performed with:\n"
        valids += f"\tcheck_bam {lr_platform} sample.bam > sample.txt"
        summary.append(valids)
        sm = "The Snakemake command invoked:\n"
        sm += f"\t{self.config["Workflow"]['snakemake']['relative']}"
        summary.append(sm)
        return "\n\n".join(summary)

    def validate_fastq(self) -> str:
        lr_platform = self.config["Workflow"]["linkedreads"]["type"]
        
        summary = ["The harpy validate fastq workflow ran using these parameters:"]
        valids = "Validations were performed with:\n"
        valids += f"\tcheck_fastq {lr_platform} sample.fastq > sample.txt"
        summary.append(valids)
        sm = "The Snakemake command invoked:\n"
        sm += f"\t{self.config["Workflow"]['snakemake']['relative']}"
        summary.append(sm)
        return "\n\n".join(summary)

    def sv_leviathan(self) -> str:
        if "groupings" in self.config["Inputs"]:
            return self.leviathan_pop()
        else:
            return self.leviathan()

    def leviathan(self) -> str:
        genomefile = os.path.basename(self.config["Inputs"]["reference"])
        min_size = self.config["Parameters"]["min-size"]
        min_bc = self.config["Parameters"]["min-barcodes"]
        iterations = self.config["Parameters"]["iterations"]
        small_thresh = self.config["Parameters"]["variant-thresholds"]["small"]
        medium_thresh = self.config["Parameters"]["variant-thresholds"]["medium"]
        large_thresh = self.config["Parameters"]["variant-thresholds"]["large"]
        duplcates_thresh = self.config["Parameters"]["variant-thresholds"]["duplicates"]
        extra = self.config["Parameters"].get("extra", "") 
        params = " ".join([
            f"-v {min_size}",
            f"-c {min_bc}",
            f"-B {iterations}",
            f"-s {small_thresh}",
            f"-m {medium_thresh}",
            f"-l {large_thresh}",
            f"-d {duplcates_thresh}",
            extra
        ])

        summary = ["The harpy sv leviathan workflow ran using these parameters:"]
        summary.append(f"The provided reference genome: {genomefile}")
        bc_idx = "The barcodes were indexed using:\n"
        bc_idx += "\tLRez index bam -p -b INPUT"
        summary.append(bc_idx)
        svcall = "Leviathan was called using:\n"
        svcall += f"\tLEVIATHAN -b INPUT -i INPUT.BCI -g GENOME {params}"
        summary.append(svcall)
        sm = "The Snakemake command invoked:\n"
        sm += f"\t{self.config["Workflow"]['snakemake']['relative']}"
        summary.append(sm)
        return "\n\n".join(summary)

    def leviathan_pop(self) -> str:
        groupfile 	= self.config["Inputs"]["groupings"]
        genomefile 	= os.path.basename(self.config["Inputs"]["reference"])
        extra 		= self.config["Parameters"].get("extra", "") 
        min_size      = self.config["Parameters"]["min-size"]
        min_bc      = self.config["min-barcodes"]
        iterations  = self.config["Parameters"]["iterations"]
        small_thresh = self.config["Parameters"]["variant-thresholds"]["small"]
        medium_thresh = self.config["Parameters"]["variant-thresholds"]["medium"]
        large_thresh = self.config["Parameters"]["variant-thresholds"]["large"]
        duplcates_thresh = self.config["Parameters"]["variant-thresholds"]["duplicates"]
        params = " ".join([
            f"-v {min_size}",
            f"-c {min_bc}",
            f"-B {iterations}",
            f"-s {small_thresh}",
            f"-m {medium_thresh}",
            f"-l {large_thresh}",
            f"-d {duplcates_thresh}",
            extra
        ])
        summary = ["The harpy sv leviathan workflow ran using these parameters:"]
        summary.append(f"The provided reference genome: {genomefile}")
        summary.append(f"The provided populations grouping file: {groupfile}")
        concat = "The alignments were concatenated using:\n"
        concat += "\tconcatenate_bam --bx -b samples.list > groupname.bam"
        summary.append(concat)
        bc_idx = "The barcodes were indexed using:\n"
        bc_idx += "\tLRez index bam -p -b INPUT"
        summary.append(bc_idx)
        svcall = "Leviathan was called using:\n"
        svcall += f"\tLEVIATHAN -b INPUT -i INPUT.BCI -g GENOME {params}"
        summary.append(svcall)
        sm = "The Snakemake command invoked:\n"
        sm += f"\t{self.config["Workflow"]['snakemake']['relative']}"
        summary.append(sm)
        return "\n\n".join(summary)

    def sv_naibr(self) -> str:
        if "groupings" in self.config["Inputs"]:
            return self.naibr_pop()
        else:
            return self.naibr()

    def naibr(self) -> str:
        genomefile = os.path.basename(self.config["Inputs"]["reference"])
        extra = self.config["Parameters"].get("extra", None) 
        min_size = self.config["Parameters"]["min-size"]
        min_barcodes = self.config["Parameters"]["min-barcodes"]
        min_quality  = self.config["Parameters"]["min-map-quality"]
        mol_dist    = self.config["Parameters"]["molecule-distance"]
        argdict = {
        "min_mapq" : min_quality,
        "d"        : mol_dist,
        "min_sv"   : min_size,
        "k"        : min_barcodes
        }
        if extra:
            words = [i for i in re.split(r"\s|=", extra) if len(i) > 0]
            for i in zip(words[::2], words[1::2], strict = True):
                if "blacklist" in i or "candidates" in i:
                    argdict[i[0].lstrip("-")] = i[1]
        
        summary = ["The harpy sv naibr workflow ran using these parameters:"]
        summary.append(f"The provided reference genome: {genomefile}")
        naibr = "naibr variant calling ran using these configurations:\n"
        naibr += "\tbam_file=BAMFILE\n"
        naibr += "\tprefix=PREFIX\n"
        naibr += "\toutdir=Variants/naibr/PREFIX\n"
        naibr += "\n\t".join([f"{k}={v}" for k,v in argdict.items()])
        summary.append(naibr)
        sm = "The Snakemake command invoked:\n"
        sm += f"\t{self.config["Workflow"]['snakemake']['relative']}"
        summary.append(sm)
        return "\n\n".join(summary)

    def naibr_pop(self) -> str:
        genomefile   = os.path.basename(self.config["Inputs"]["reference"])
        groupfile    = self.config["Inputs"]["groupings"]
        extra        = self.config["Parameters"].get("extra", None) 
        min_size     = self.config["Parameters"]["min-size"]
        min_barcodes = self.config["Parameters"]["min-barcodes"]
        min_quality  = self.config["Parameters"]["min-map-quality"]
        mol_dist     = self.config["Parameters"]["molecule-distance"]
        argdict = {
        "min_mapq" : min_quality,
        "d"        : mol_dist,
        "min_sv"   : min_size,
        "k"        : min_barcodes
        }
        if extra:
            words = [i for i in re.split(r"\s|=", extra) if len(i) > 0]
            for i in zip(words[::2], words[1::2], strict = True):
                if "blacklist" in i or "candidates" in i:
                    argdict[i[0].lstrip("-")] = i[1]

        summary = ["The harpy sv naibr workflow ran using these parameters:"]
        summary.append(f"The provided reference genome: {genomefile}")
        summary.append(f"The provided populations grouping file: {groupfile}")
        concat = "The alignments were concatenated using:\n"
        concat += "\tconcatenate_bam -b samples.list > groupname.bam"
        summary.append(concat)
        naibr = "naibr variant calling ran using these configurations:\n"
        naibr += "\tbam_file=BAMFILE\n"
        naibr += "\tprefix=PREFIX\n"
        naibr += "\toutdir=Variants/naibr/PREFIX\n"
        naibr += "\n\t".join([f"{k}={v}" for k,v in argdict.items()])
        summary.append(naibr)
        sm = "The Snakemake command invoked:\n"
        sm += f"\t{self.config["Workflow"]['snakemake']['relative']}"
        summary.append(sm)
        return "\n\n".join(summary)

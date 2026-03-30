"""Basic functions to write the workflow summaries"""
import os
from harpy.common.file_ops import naibr_extra

class Summary:
    def __init__(self, version, config: dict):
        self.summary: list[str] = [f"Harpy Version: {version}"]
        self.WORKFLOW    = config.get('Workflow') or {}
        self.PARAMETERS = config.get('Parameters') or {}
        self.INPUTS     = config['Inputs']

    def get(self) -> str:
        self.__getattribute__(self.WORKFLOW["name"])()
        self.summary.append("The Snakemake command invoked:\n\t" + self.WORKFLOW['snakemake']['relative'])
        return "\n\n".join(self.summary)

    def align_bwa(self):
        ignore_bx = self.WORKFLOW.get("linkedreads", {}).get("type", 'none') == "none"
        bx_tag = self.WORKFLOW.get("linkedreads", {}).get("standardized", {}).get("BX", False)
        vx_tag = self.WORKFLOW.get("linkedreads", {}).get("standardized", {}).get("VX", False)

        is_standardized = bx_tag and vx_tag
        keep_unmapped = self.PARAMETERS.get("keep-unmapped", False)
        extra 		= self.PARAMETERS.get("extra", "") 
        genomefile 	= self.INPUTS["reference"]
        quality = self.PARAMETERS.get("min-map-quality", 30)

        unmapped = "" if keep_unmapped else "-F 4"
        bx_mode = "--barcode-tag BX" if not ignore_bx else ""
        bwa_static = "-C -v 2" if is_standardized else "-v 2"
        extra   = extra

        align = "Sequences were aligned with BWA using:\n"
        align += f'\tbwa mem {bwa_static} {extra} -R "@RG\\tID:SAMPLE\\tSM:SAMPLE" genome forward_reads reverse_reads |\n'
        align += f"\tsamtools view -h {unmapped} -q {quality}"
        duplicates = "Duplicates in the alignments were marked following:\n"
        duplicates += "\tsamtools collate |\n"
        duplicates += "\tsamtools fixmate |\n"
        duplicates += f"\tsamtools sort -T SAMPLE -m 2000M |\n"
        duplicates += f"\tsamtools markdup -S {bx_mode} -d 100 (2500 for novaseq)"
        standardization = "Barcodes were standardized to BX + VX format in the aligments using:\n"
        standardization += "\tdjinn-standardize {input.bam} > {output.bam}"
        self.summary.append("The harpy align bwa workflow ran using these parameters:")
        self.summary.append(f"The provided genome: {genomefile}")
        self.summary.append(align)
        if not ignore_bx:
            self.summary.append(standardization)
        self.summary.append(duplicates)

    def align_strobe(self):
        ignore_bx = self.WORKFLOW.get("linkedreads", {}).get("type", 'none') == "none"
        bx_tag = self.WORKFLOW.get("linkedreads", {}).get("standardized", {}).get("BX", False)
        vx_tag = self.WORKFLOW.get("linkedreads", {}).get("standardized", {}).get("VX", False)

        is_standardized = bx_tag and vx_tag
        keep_unmapped = self.PARAMETERS.get("keep-unmapped", False)
        extra 		= self.PARAMETERS.get("extra", "") 
        genomefile 	= self.INPUTS["reference"]
        quality = self.PARAMETERS.get("min-map-quality", 30)

        unmapped_strobe = "" if keep_unmapped else "-U"
        unmapped = "" if keep_unmapped else "-F 4"
        bx_mode = "--barcode-tag BX" if not ignore_bx else ""
        static = "-C" if is_standardized else ""
        extra = self.PARAMETERS.get("extra", "") 

        align = "Sequences were aligned with strobealign using:\n"
        align += f"\tstrobealign {unmapped_strobe} {static} --rg-id=SAMPLE --rg=SM:SAMPLE {extra} genome reads.F.fq reads.R.fq |\n"
        align += f"\t\tsamtools view -h {unmapped} -q {quality}"
        duplicates = "Duplicates in the alignments were marked following:\n"
        duplicates += "\tsamtools collate |\n"
        duplicates += "\tsamtools fixmate |\n"
        duplicates += f"\tsamtools sort -T SAMPLE --reference {genomefile} -m 2000M |\n"
        duplicates += f"\tsamtools markdup -S {bx_mode} -d 100 (2500 for novaseq)"
        standardization = "Barcodes were standardized in the aligments using:\n"
        standardization += "\tstandardize-barcodes-sam > {output} < {input}"
        self.summary.append("The harpy align strobe workflow ran using these parameters:")
        self.summary.append(f"The provided genome: {genomefile}")
        self.summary.append(align)
        if not ignore_bx:
            self.summary.append(standardization)
        self.summary.append(duplicates)

    def assembly(self):
        # SPADES
        max_mem      = self.PARAMETERS.get("spades", {}).get("max-memory", 'auto')
        k_param      = self.PARAMETERS.get("spades", {}).get("k", 10000)
        spades_extra = self.PARAMETERS.get("spades", {}).get("extra", "")
        # ARCS
        mapq       = self.PARAMETERS.get("tigmint", {}).get("min-mapping-quality", 0)
        mismatch   = self.PARAMETERS.get("tigmint", {}).get("mismatch", 5)
        mol_dist   = self.PARAMETERS.get("tigmint", {}).get("molecule-distance", 50000)
        mol_len    = self.PARAMETERS.get("tigmint", {}).get("molecule-length", 2000)
        span       = self.PARAMETERS.get("tigmint", {}).get("span", 20)
        min_align  = self.PARAMETERS.get("arcs", {}).get("min-aligned-reads", 5)
        min_contig = self.PARAMETERS.get("arcs", {}).get("min-contig-length", 500)
        seq_id     = self.PARAMETERS.get("arcs", {}).get("min-sequence-identity", 98)
        arcs_extra = self.PARAMETERS.get("arcs", {}).get("extra", "")
        links      = self.PARAMETERS.get("links", {}).get("min-links", 5)

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
        spades = "Reads were assembled using cloudspades:\n"
        spades += f"\tspades.py -t THREADS -m {max_mem} --gemcode1-1 FQ1 --gemcode1-2 FQ2 --isolate -k {k_param} {spades_extra}"
        arcs = "The draft assembly was error corrected and scaffolded with Tigmint/ARCS/LINKS:\n"
        arcs += f"\tarcs-make arcs-tigmint {' '.join(params[3:])}"
        self.summary.append("The harpy assemble workflow ran using these parameters:")
        self.summary.append(spades)
        self.summary.append(arcs)

    def deconvolve(self):
        kmer_length = self.PARAMETERS.get("kmer-length", 21)
        window_size = self.PARAMETERS.get("window-size", 40)
        density 	= self.PARAMETERS.get("density", 3) 
        dropout     = self.PARAMETERS.get("dropout", 0)

        interleave = "fastq files were interleaved with seqtk:\n"
        interleave += "\tseqtk mergepe forward.fq reverse.fq"
        deconv = "Deconvolution occurred using QuickDeconvolution:\n"
        deconv += f"\tQuickDeconvolution -t threads -i infile.fq -o output.fq -k {kmer_length} -w {window_size} -d {density} -a {dropout}"
        recover = "The interleaved output was split back into forward and reverse reads with seqtk:\n"
        recover += "\tseqtk seq -1 interleaved.fq | gzip > file.R1.fq.gz\n"
        recover += "\tseqtk seq -2 interleaved.fq | gzip > file.R2.fq.gz"
        self.summary.append("The harpy deconvolve workflow ran using these parameters:")
        self.summary.append(interleave)
        self.summary.append(deconv)
        self.summary.append(recover)

    def impute(self):
        region = self.PARAMETERS.get("region", None)
        if region:
            _,positions = region.split(":")
            startpos,endpos,buffer = [int(i) for i in positions.split("-")]
            regiontext = f"\t\tregionStart = {startpos},\n"
            regiontext += f"\t\tregionEnd = {endpos},\n"
            regiontext += f"\t\tbuffer = {buffer},\n"
        else:
            regiontext = ""
        paramfiletext = "\t".join(open(self.INPUTS["parameters"], "r").readlines())
        preproc = "Preprocessing was performed with:\n"
        preproc += "\tbcftools view -M2 -v snps --regions CONTIG INFILE |\n"
        preproc += """\tbcftools query -i '(STRLEN(REF)==1) & (STRLEN(ALT[0])==1) & (REF!="N")' -f '%CHROM\\t%POS\\t%REF\\t%ALT\\n'"""
        stitchparam = f"The STITCH parameter file: {self.INPUTS['parameters']}\n"
        stitchparam += f"\t{paramfiletext}"
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
        if self.PARAMETERS.get("grid-size", 1) > 1:
            stitch += f"\t\tgridWindowSize = {self.PARAMETERS.get("grid-size", 1)}\n"
        stitch += "\t\toutputdir = outdir,\n"
        stitch += "\t\toutput_filename = outfile\n\t)"
        stitchextra = "Additional STITCH parameters provided (overrides existing values above):\n"
        stitchextra += "\t" + self.PARAMETERS.get("extra", "None")
        self.summary.append("The harpy impute workflow ran using these parameters:")
        self.summary.append(f"The provided variant file: {self.INPUTS['vcf']}")
        self.summary.append(preproc)
        self.summary.append(stitchparam)
        self.summary.append(stitchextra)

    def metassembly(self):
        BX_TAG       = self.WORKFLOW.get("linkedreads", {})["barcode_tag"]
        max_mem      = self.PARAMETERS.get("spades", {}).get("max_memory", 10000)
        k_param      = self.PARAMETERS.get("spades", {}).get("k", 'auto')
        ignore_bx    = self.PARAMETERS.get("spades", {}).get("ignore_barcodes", False)
        extra        = self.PARAMETERS.get("spades", {}).get("extra", "")
        force_athena = self.PARAMETERS.get("athena", {}).get("force", False)
        force = "--force_reads" if force_athena else ""
        extra = self.PARAMETERS["spades"].get("extra", "")
        spadesdir = f"{'cloudspades' if not ignore_bx else 'spades'}_assembly"

        bxsort = "FASTQ inputs were sorted by their linked-read barcodes:\n"
        bxsort += "\tsamtools import -T \"*\" FQ1 FQ2 |\n"
        bxsort += f"\tsamtools sort -O SAM -t {BX_TAG} |\n"  
        bxsort += "\tsamtools fastq -T \"*\" -1 FQ_out1 -2 FQ_out2"  
        bxappend = "Barcoded-sorted FASTQ files had \"-1\" appended to the barcode to make them Athena-compliant:\n"  
        bxappend += f"\tsed 's/{BX_TAG}:Z:[^[:space:]]*/&-1/g' FASTQ | bgzip > FASTQ_OUT"  
        if not ignore_bx:
            spades = "Reads were assembled using cloudspades:\n"
            spades += f"\tspades.py -t THREADS -m {max_mem} --gemcode1-1 FQ1 --gemcode1-2 FQ2 --meta -k {k_param} {extra}"
        else:
            spades = "Reads were assembled using spades:\n"
            spades += f"\tmetaspades.py -t THREADS -m {max_mem} -k {k_param} {extra} -1 FQ_1 -2 FQ2 -o {spadesdir}"
        align = "Original input FASTQ files were aligned to the metagenome using BWA:\n"
        align += "\tbwa mem -C -p spades.contigs FQ1 FQ2 | samtools sort -O bam -"
        interleaved = "Barcode-sorted Athena-compliant sequences were interleaved with seqtk:\n"
        interleaved += "\tseqtk mergepe FQ1 FQ2 > INTERLEAVED.FQ"
        athena = "Athena ran with the config file Harpy built from the files created from the previous steps:\n"
        athena += f"\tathena-meta {force} --config athena.config"
        self.summary.append("The harpy metassembly workflow ran using these parameters:")
        self.summary.append(bxsort)
        self.summary.append(bxappend)
        self.summary.append(spades)
        self.summary.append(align)
        self.summary.append(interleaved)
        self.summary.append(athena)

    def phase_snp(self):
        bc_type           = self.WORKFLOW.get("linkedreads", {}).get("type", 'none')
        pruning           = self.PARAMETERS.get("prune", 30)
        map_qual          = self.PARAMETERS.get("min-map-quality", 20)
        base_qual         = self.PARAMETERS.get("min-base-quality", 13)
        molecule_distance = self.PARAMETERS.get("distance-threshold", 100000)

        extra             = self.PARAMETERS.get("extra", "") 
        variantfile       = self.INPUTS["vcf"]
        invalid_regex = {
            "haplotagging" : "'$4 !~ /[ABCD]00/'",
            "stlfr" : "'$4 !~ /^0_|_0_|_0$/'",
            "tellseq": "'$4 !~ /N/'"
        }
        linkarg = "--10x 0" if bc_type == "none" else "--10x 1"
        indelarg   = "--indels 1 --ref reference.fasta" if self.INPUTS.get("reference", None) else ""
        hairs_params = f"{indelarg} {linkarg} --mmq {map_qual} --mbq {base_qual} --nf 1 --maxfragments 1500000"
        prune = f"--threshold {pruning}" if pruning > 0 else "--no_prune 1"

        hetsplit = "The variant file was split by sample and filtered for heterozygous sites using:\n"
        hetsplit += "\tbcftools view -s SAMPLE | bcftools view -m 2 -M 2 -i \'GT=\"het\"\'"
        phase = "Phasing was performed using the components of HapCut2:\n"
        phase += f"\textractHAIRS {hairs_params} --bam sample.bam --VCF sample.vcf --out sample.unlinked.frags\n"
        if bc_type != "none":
            phase += "\t awk " + invalid_regex.get(bc_type, "'$4 !~ /N/'") + " sample.unlinked.frags > sample.frags.filt"
            phase += f"\tLinkFragments.py --bam sample.bam --VCF sample.vcf --fragments sample.frags.filt --out sample.linked.frags -d {molecule_distance}\n"
            phase += f"\tHAPCUT2 --fragments sample.linked.frags --vcf sample.vcf --out sample.blocks --nf 1 --error_analysis_mode 1 --call_homozygous 1 --outvcf 1 {prune} {extra}\n"
        else:
            phase += f"\tHAPCUT2 --fragments sample.unlinked.frags --vcf sample.vcf --out sample.blocks --nf 1 --error_analysis_mode 1 --call_homozygous 1 --outvcf 1 {prune} {extra}\n"
        annot = "Variant annotation was performed using:\n"
        annot += "\tbcftools annotate -a sample.phased.vcf -c CHROM,POS,FMT/GT,FMT/PS,FMT/PQ,FMT/PD -m +HAPCUT\n"
        annot += "\tbcftools merge --output-type b samples.annot.bcf"
        self.summary.append("The harpy phase bam workflow ran using these parameters:")
        self.summary.append(f"The provided variant file: {variantfile}")
        self.summary.append(hetsplit)
        self.summary.append(phase)
        self.summary.append(annot)

    def phase_bam(self):
        mol_dist    = self.PARAMETERS.get("distance-threshold", 100000)
        extra       = self.PARAMETERS.get("extra", "") 
        ploidy      = self.PARAMETERS.get("ploidy", 2) 
        variantfile = self.INPUTS["vcf"]
        params = [
            f"--ploidy {ploidy}",
            f"-d {mol_dist}",
            "--tag-supplementary copy-primary",
            "--no-supplementary-strand-match",
            f"--supplementary-distance {3 * mol_dist}",
            "--ignore-read-groups",
            "--skip-missing-contigs",
        ]
        params = " ".join(params)
        validsplit = "The input alignments had their records filtered for valid barcodes:\n"
        validsplit += "\tdjinn sam filter-invalid --invalid sample.invalid.bam sample.bam"
        phase = "Phasing was performed using whatshap:\n"
        phase += f"\twhatshap haplotag --sample name --reference input.ref {params} {extra} input.vcf input.bam"
        concataln = "Invalid-barcode alignments were added back to the phased alignments using:\n"
        concataln += "\tsamtools merge sample.phased.bam sample.invalid.bam | samtools sort -"
        self.summary.append("The harpy phase snp workflow ran using these parameters:")
        self.summary.append(f"The provided variant file: {variantfile}")
        self.summary.append(validsplit)
        self.summary.append(phase)
        self.summary.append(concataln)


    def preprocess_meier2021(self):
        schemafile = self.INPUTS["schema"]
        qxrx             = self.PARAMETERS.get("qx-rx", False)
        unknown_samples  = self.PARAMETERS.get("samples", False)
        unknown_barcodes = self.PARAMETERS.get("barcodes", False)

        R1 = self.INPUTS["R1"],
        R2 = self.INPUTS["R2"],
        I1 = self.INPUTS["I1"],
        I2 = self.INPUTS["I2"],
        outdir = f"--samples {os.getcwd()}",
        qxrx = "--rx --qx" if qxrx else "",
        unknown_barcodes = "--undetermined-barcodes _unknown_barcodes" if unknown_barcodes else "",
        unknown_samples = "--undetermined-samples _unknown_samples" if unknown_samples else ""

        inputs = "The multiplexed input files:\n"
        inputs += f"\tread 1: {R1}\n"
        inputs += f"\tread 2: {R2}\n"
        inputs += f"\tindex 1: {I1}\n"
        inputs += f"\tindex 2: {I2}\n"
        inputs += f"Sample demultiplexing schema: {schemafile}"
        demux = "Samples were demultiplexed using:\n"
        demux += f"\tdmox --R1 --R2 --I1 --I2 {outdir} {qxrx} {unknown_barcodes} {unknown_samples}"
        qc = "QC checks were performed on demultiplexed FASTQ files using:\n"
        qc += "\tfalco -skip-report -skip-summary -data-filename output input.fq.gz"
        self.summary.append("The harpy preprocess workflow ran using these parameters:")
        self.summary.append("Linked Read Barcode Design: Meier et al. 2021")
        self.summary.append(inputs)
        self.summary.append(demux)
        self.summary.append(qc)

    def preprocess_gih(self):
        me_seq   = self.PARAMETERS.get("ME-sequence", "AGATGTGTATAAGAGACAG")
        mismatch = self.PARAMETERS.get("ME-mismatch", 1) 
        minlen   = self.PARAMETERS.get("min-length", 10) 

        stagger = "Input FASTQs had the ME sequence identified and removed, then provided a nucleotide padding sequence (if necessary):\n"
        stagger += f"\tgih-stagger --me {me_seq} --max-mismatch {mismatch} --min-len {minlen} --stats output.stats FQ1 FQ2 > output.sam"
        pheniqs = "The resulting interleaved unaligned SAM file was then piped into the Pheniqs for barcode demultiplexing:\n"
        pheniqs += "\tpheniqs mux --output output.sam --quality -c pheniqs.config.json --report output.json < input.sam"
        recode = "Finally, the nucleotides were converted into standardized haplotagging ACBD format and FASTQ format:\n"
        recode += "\tgih-convert input.sam fq1 fq2 > output.stats"
        qc = "QC checks were performed on demultiplexed FASTQ files using:\n"
        qc += "\tfalco -skip-report -skip-summary -data-filename output input.fq.gz"
        self.summary.append("The harpy preprocess workflow ran using these parameters:")
        self.summary.append("Linked Read Barcode Design: Iqbal et al. 2026")
        self.summary.append(stagger)
        self.summary.append(pheniqs)
        self.summary.append(recode)
        self.summary.append(qc)

    def qc(self):
        min_len 	  = self.PARAMETERS.get("min-len", 30)
        max_len 	  = self.PARAMETERS.get("max-len", 150)
        extra 	      = self.PARAMETERS.get("extra", "") 
        trim_adapters = self.PARAMETERS.get("trim_adapters", None)
        dedup         = self.PARAMETERS.get("deduplicate", False)

        if trim_adapters:
            trim_arg = "--detect_adapter_for_pe" if trim_adapters == "auto" else f"--adapter_fasta {trim_adapters}"
        else:
            trim_arg = "--disable_adapter_trimming"

        params = [
            "--trim_poly_g",
            "--cut_right",
            f"--length_required {min_len}",
            f"--max_len1 {max_len}",
            trim_arg,
            "-D" if dedup else "",
            extra
        ]

        fastp = "fastp ran using:\n"
        fastp += "\tfastp  " + " ".join(params)
        self.summary.append("The harpy qc workflow ran using these parameters:")
        self.summary.append(fastp)
    
    def snp_freebayes(self):
        ploidy 		  = self.PARAMETERS.get("ploidy", 2)
        extra 	      = self.PARAMETERS.get("extra", "") 
        genomefile 	  = self.INPUTS["reference"]
        regions_input = self.INPUTS["regions"]
        groupings 	  = self.INPUTS.get("groupings", [])

        params = f"-p {ploidy} "
        params += f"--populations {groupings} " if groupings else ''
        params += extra

        varcall = "The freebayes parameters:\n"
        varcall += f"\tfreebayes -f REFERENCE -L samples.list -r REGION {params} |\n"
        varcall += "\tbcftools sort -"
        merged = "The variants identified in the intervals were merged into the final variant file using:\n"
        merged += "\tbcftools concat -f bcf.files -a --remove-duplicates"
        normalize = "The variants were normalized using:\n"
        normalize += "\tbcftools norm -m -both -d both -c w"
        self.summary.append("The harpy snp freebayes workflow ran using these parameters:")
        self.summary.append(f"The provided reference genome: {genomefile}")
        self.summary.append(f"Genomic positions for which variants were called: {regions_input}")
        self.summary.append(varcall)
        self.summary.append(merged)
        self.summary.append(normalize)
    
    def snp_mpileup(self):
        ploidy 		 = self.PARAMETERS.get("ploidy", 2)
        mp_extra 	 = self.PARAMETERS.get("extra", "")
        genomefile 	 = self.INPUTS["reference"]
        groupings 	 = self.INPUTS.get("groupings", [])
        region_input = self.INPUTS["regions"]

        params = f"--ploidy {ploidy} --populations "
        params += f"{groupings}" if groupings else "-"

        mpileup = "The mpileup parameters:\n"
        mpileup += f"\tbcftools mpileup --fasta-ref REFERENCE --region REGION --bam-list BAMS --annotate AD --output-type b {mp_extra}"
        bcfcall = "The bcftools call parameters:\n"
        bcfcall += f"\tbcftools call --multiallelic-caller {params} --variants-only --output-type b |\n"
        bcfcall += "\tbcftools sort -"
        merged = "The variants identified in the intervals were merged into the final variant file using:\n"
        merged += "\tbcftools concat -f bcf.files -a --remove-duplicates"
        normalize = "The variants were normalized using:\n"
        normalize += "\tbcftools norm -m -both -d both -c w"
        self.summary.append("The harpy snp freebayes workflow ran using these parameters:")
        self.summary.append(f"The provided reference genome: {genomefile}")
        self.summary.append(f"Genomic positions for which variants were called: {region_input}")
        self.summary.append(mpileup)
        self.summary.append(bcfcall)
        self.summary.append(merged)
        self.summary.append(normalize)


    def sv_leviathan(self):
        genomefile = os.path.basename(self.INPUTS["reference"])
        groupfile 	= self.INPUTS.get("groupings", None)
        extra         = self.PARAMETERS.get("extra", "")
        min_size      = self.PARAMETERS.get("min-size", 1000)
        min_bc        = self.PARAMETERS.get("min-barcodes", 2)
        iterations    = self.PARAMETERS.get("iterations", 50)
        small_thresh  = self.PARAMETERS.get("variant-thresholds", {}).get("small", 95)
        medium_thresh = self.PARAMETERS.get("variant-thresholds", {}).get("medium", 95)
        large_thresh  = self.PARAMETERS.get("variant-thresholds", {}).get("large", 95)
        dup_thresh    = self.PARAMETERS.get("variant-thresholds", {}).get("duplicates", 10)

        params = " ".join([
            f"-v {min_size}",
            f"-c {min_bc}",
            f"-B {iterations}",
            f"-s {small_thresh}",
            f"-m {medium_thresh}",
            f"-l {large_thresh}",
            f"-d {dup_thresh}",
            extra
        ])

        bc_idx = "The barcodes were indexed using:\n"
        bc_idx += "\tLRez index bam -p -b INPUT"
        svcall = "Leviathan was called using:\n"
        svcall += f"\tLEVIATHAN -b INPUT -i INPUT.BCI -g GENOME {params}"
        self.summary.append("The harpy sv leviathan workflow ran using these parameters:")
        self.summary.append(f"The provided reference genome: {genomefile}")
        if groupfile:
            self.summary.append(f"The provided populations grouping file: {groupfile}")
            self.summary.append("The alignments were concatenated using:\n\tdjinn sam concat samples.bam... > group.bam")
        self.summary.append(bc_idx)
        self.summary.append(svcall)

    def sv_naibr(self):
        genomefile   = os.path.basename(self.INPUTS["reference"])
        groupfile    = self.INPUTS.get("groupings", None)
        extra        = self.PARAMETERS.get("extra", None) 
        min_size     = self.PARAMETERS.get("min-size", 1000)
        min_barcodes = self.PARAMETERS.get("min-barcodes", 2)
        min_quality  = self.PARAMETERS.get("min-map-quality", 30)
        mol_dist     = self.PARAMETERS.get("molecule-distance", 100000)

        argdict = naibr_extra(
            {"min_mapq" : min_quality, "d" : mol_dist, "min_sv" : min_size, "k": min_barcodes},
            extra
        )
        naibr = "naibr variant calling ran using these configurations:\n"
        naibr += "\tbam_file=BAMFILE\n"
        naibr += "\tprefix=PREFIX\n"
        naibr += "\toutdir=Variants/naibr/PREFIX\n"
        naibr += "\n\t".join([f"{k}={v}" for k,v in argdict.items()])
        self.summary.append("The harpy sv naibr workflow ran using these parameters:")
        self.summary.append(f"The provided reference genome: {genomefile}")
        if groupfile:
            self.summary.append(f"The provided populations grouping file: {groupfile}")
            self.summary.append("The alignments were concatenated using:\n\tdjinn sam concat samples.bam... > group.bam")
        self.summary.append(naibr)


    def validate_bam(self) -> str:
        lr_platform = self.WORKFLOW.get("linkedreads", {}).get("type", 'none')
        self.summary.append("The harpy validate bam workflow ran using these parameters:")
        valids = "Validations were performed with:\n"
        valids += f"\tharpy-utils check-bam {lr_platform} sample.bam > sample.txt"
        self.summary.append(valids)
        sm = "The Snakemake command invoked:\n"
        sm += f"\t{self.WORKFLOW['snakemake']['relative']}"
        self.summary.append(sm)
        return "\n\n".join(self.summary)

    def validate_fastq(self):
        lr_platform = self.WORKFLOW.get("linkedreads", {}).get("type", 'none')
        valids = "Validations were performed with:\n"
        valids += f"\tharpy-utils check-fastq {lr_platform} sample.fastq > sample.txt"
        self.summary.append("The harpy validate fastq workflow ran using these parameters:")
        self.summary.append(valids)

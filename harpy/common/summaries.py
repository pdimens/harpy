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
        return self.summary + "\nThe Snakemake command invoked:\n\t" + self.WORKFLOW['snakemake']['relative']

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
        bwa_static = "-y -x sr" if is_standardized else "-x sr"
        extra   = extra

        self.summary = f'''The harpy align bwa workflow ran using these parameters:
The provided genome: {genomefile}

Sequences were aligned with minibwa using:
    minibwa map {bwa_static} {extra} -R "@RG\\tID:SAMPLE\\tSM:SAMPLE" genome forward_reads reverse_reads |
    samtools view -h {unmapped} -q {quality}"

Duplicates in the alignments were marked following:
    samtools collate |
    samtools fixmate |
    samtools sort -T SAMPLE -m 2000M |
    samtools markdup -S {bx_mode} -d 100   #(2500 for novaseq)

If linked reads, barcodes were standardized to BX + VX format in the aligments using:
    djinn-standardize {{input.bam}} > {{output.bam}}
'''

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

        self.summary = f'''The harpy align strobe workflow ran using these parameters:
The provided genome: {genomefile}

Sequences were aligned with strobealign using:
    strobealign {unmapped_strobe} {static} --rg-id=SAMPLE --rg=SM:SAMPLE {extra} genome reads.F.fq reads.R.fq |
    samtools view -h {unmapped} -q {quality}

Duplicates in the alignments were marked following:
    samtools collate |
    samtools fixmate |
    samtools sort -T SAMPLE --reference {genomefile} -m 2000M |
    samtools markdup -S {bx_mode} -d 100   #(2500 for novaseq)

If linked reads, barcodes were standardized to BX + VX format in the aligments using:
    djinn-standardize {{input.bam}} > {{output.bam}}
'''

    def align_minimap(self):
        ignore_bx = self.WORKFLOW.get("linkedreads", {}).get("type", 'none') == "none"
        bx_tag = self.WORKFLOW.get("linkedreads", {}).get("standardized", {}).get("BX", False)
        vx_tag = self.WORKFLOW.get("linkedreads", {}).get("standardized", {}).get("VX", False)
        tech = self.PARAMETERS.get("aligner-technology", "sr") 
        is_standardized = bx_tag and vx_tag
        keep_unmapped = self.PARAMETERS.get("keep-unmapped", False)
        extra 		= self.PARAMETERS.get("extra", "")
        genomefile 	= self.INPUTS["reference"]
        quality = self.PARAMETERS.get("min-map-quality", 30)

        unmapped = "" if keep_unmapped else "-F 4"
        bx_mode = "--barcode-tag BX" if not ignore_bx else ""
        tech = f"-ax map-{tech}" if tech != "sr" else "-ax sr"
        static = "-y --MD" if is_standardized else "--MD"
        extra = self.PARAMETERS.get("extra", "")

        self.summary = f'''The harpy align minimap workflow ran using these parameters:
The provided genome: {genomefile}

Sequences were aligned with minimap2 using:
    minimap2 {tech} {static} -R \"@RG\\tID:SAMPLE\\tSM:SAMPLE\" {extra} genome reads.F.fq reads.R.fq |
    samtools view -h {unmapped} -q {quality}

Duplicates in the alignments were marked following:
    samtools collate |
    samtools fixmate |
    samtools sort -T SAMPLE --reference {genomefile} -m 2000M |
    samtools markdup -S {bx_mode} -d 100   #(2500 for novaseq)

If linked reads, barcodes were standardized to BX + VX format in the aligments using:
    djinn-standardize {{input.bam}} > {{output.bam}}
'''

    def align_arachne(self):
        keep_unmapped = self.PARAMETERS.get("keep-unmapped", False)
        extra 		= self.PARAMETERS.get("extra", "")
        genomefile 	= self.INPUTS["reference"]
        quality = self.PARAMETERS.get("min-map-quality", 30)

        unmapped = "" if keep_unmapped else "-F 4"
        extra = self.PARAMETERS.get("extra", "")

        self.summary = f'''The harpy align arachne workflow ran using these parameters:
The provided genome: {genomefile}

FASTQ files were standardized, sorted by barcode, and filtered for valid barcodes:
    arachne prep prefix reads.F.fq reads.R.fq

Valid-barcoded FASTQ files were aligned with arachne:
    arachne align -s samplename {genomefile} reads.valid.F.fq reads.valid.R.fq |
    samtools sort -u > arachne.bam

Invalid-barcoded FASTQ files were aligned with minibwa:
    minibwa map {extra} -R "@RG\\tID:SAMPLE\\tSM:SAMPLE" genome reads.invalid.F.fq reads.invalid.R.fq  |
    samtools view -h {unmapped} -q {quality}

Duplicates in the minibwa alignments were marked following:
    samtools collate |
    samtools fixmate |
    samtools sort -T SAMPLE --reference {genomefile} -m 2000M |
    samtools markdup -S -d 100   #(2500 for novaseq)

Processed arachne and minibwa alignments were concatenated and sorted:
    samtools cat --output-fmt-option level=0 arachne.bam minibwa.bam |
    samtools sort -O BAM > final.bam
'''

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
        self.summary =  f'''The harpy assemble workflow ran using these parameters
Reads were assembled using cloudspades:
    spades.py -t THREADS -m {max_mem} --gemcode1-1 FQ1 --gemcode1-2 FQ2 --isolate -k {k_param} {spades_extra}

The draft assembly was error corrected and scaffolded with Tigmint/ARCS/LINKS:
    arcs-make arcs-tigmint {' '.join(params[3:])}
'''

    def deconvolve(self):
        kmer_length = self.PARAMETERS.get("kmer-length", 21)
        window_size = self.PARAMETERS.get("window-size", 40)
        density 	= self.PARAMETERS.get("density", 3)
        dropout     = self.PARAMETERS.get("dropout", 0)

        self.summary = f'''The harpy deconvolve workflow ran using these parameters:
fastq files were interleaved with seqtk:
    seqtk mergepe forward.fq reverse.fq

Deconvolution occurred using QuickDeconvolution:
    QuickDeconvolution -t threads -i infile.fq -o output.fq -k {kmer_length} -w {window_size} -d {density} -a {dropout}

The interleaved output was split back into forward and reverse reads with seqtk:
    seqtk seq -1 interleaved.fq | gzip > file.R1.fq.gz
    seqtk seq -2 interleaved.fq | gzip > file.R2.fq.gz
'''

    def impute(self):
        region = self.PARAMETERS.get("region", None)
        window = self.PARAMETERS.get("window-size", None)
        buffer = self.PARAMETERS.get("buffer", None)
        regiontext = ""
        if region:
            _,positions = region.split(":")
            startpos,endpos = map(int, positions.split("-"))
            regiontext += f"\t\tregionStart = {startpos},\n"
            regiontext += f"\t\tregionEnd = {endpos},\n"
            regiontext += f"\t\tbuffer = {buffer},\n"
        elif window:
            regiontext += f"\t\tbuffer = {buffer}"

        if self.PARAMETERS.get("grid-size", 1) > 1:
            gridparam = f"\n\t\tgridWindowSize = {self.PARAMETERS.get("grid-size", 1)}\n"
        paramfiletext = "\t".join(open(self.INPUTS["parameters"], "r").readlines())
        self.summary = f'''The harpy impute workflow ran using these parameters:

The provided variant file: {self.INPUTS['vcf']}

Preprocessing was performed with:
    bcftools view -M2 -v snps --regions CONTIG INFILE |
    bcftools query -i '(STRLEN(REF)==1) & (STRLEN(ALT[0])==1) & (REF!="N")' -f '%CHROM\\t%POS\\t%REF\\t%ALT\\n'

The STITCH parameter file: {self.INPUTS['parameters']}
    {paramfiletext}"
    
Within R, STITCH was invoked with the following parameters:
    STITCH(
        method = model,
        posfile = posfile,
        bamlist = bamlist,
        nCores = ncores,
        nGen = ngen,
        chr = chr,
{regiontext}
        K = k,
        S = s,
        use_bx_tag = usebx,
        bxTagUpperLimit = bxlimit,
        niterations = 40,
        switchModelIteration = 39,{gridparam}
        splitReadIterations = NA,
        outputdir = outdir
        output_filename = outfile
    )

Additional STITCH parameters provided (overrides existing values above):
        {self.PARAMETERS.get("extra", "None")}
'''

    def metassembly(self):
        BX_TAG       = self.WORKFLOW.get("linkedreads", {}).get("barcode-tag", "BX")
        max_mem      = self.PARAMETERS.get("spades", {}).get("max-memory", 10000)
        k_param      = self.PARAMETERS.get("spades", {}).get("k", 'auto')
        ignore_bx    = self.PARAMETERS.get("spades", {}).get("ignore-barcodes", False)
        extra        = self.PARAMETERS.get("spades", {}).get("extra", "")
        force_athena = self.PARAMETERS.get("athena", {}).get("force", False)
        force = "--force_reads" if force_athena else ""
        extra = self.PARAMETERS["spades"].get("extra", "")
        spadesdir = f"{'cloudspades' if not ignore_bx else 'spades'}_assembly"

        self.summary = f'''The harpy metassembly workflow ran using these parameters:
FASTQ inputs were sorted by their linked-read barcodes and had '-1' appended to the barcode to make them Athena-compliant:
    samtools import -T "*" FQ1 FQ2 |
    sed s/{BX_TAG}:Z:[^[:space:]]*/&-1/g |
    samtools sort -O SAM -t {BX_TAG} |
    samtools fastq -T "*" -1 FQ_out1 -2 FQ_out2

If the data were linked-reads, reads were assembled using cloudspades:
    spades.py -t THREADS -m {max_mem} --gemcode1-1 FQ1 --gemcode1-2 FQ2 --meta -k {k_param} {extra}
    
Otheriwse, they were assembled using spades:
    metaspades.py -t THREADS -m {max_mem} -k {k_param} {extra} -1 FQ_1 -2 FQ2 -o {spadesdir}

Original input FASTQ files were aligned to the metagenome using BWA:
    bwa mem -C -p spades.contigs FQ1 FQ2 | samtools sort -O bam -

Barcode-sorted Athena-compliant sequences were interleaved with seqtk:
    seqtk mergepe FQ1 FQ2 > INTERLEAVED.FQ
    
Athena ran with the config file Harpy built from the files created from the previous steps:
    athena-meta {force} --config athena.config
'''

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

        self.summary =  f'''The harpy phase bam workflow ran using these parameters:
The provided variant file: {variantfile}

The variant file was split by sample and filtered for heterozygous sites using:
    bcftools view -s SAMPLE | bcftools view -m 2 -M 2 -i 'GT="het"'

Phasing was performed using the components of HapCut2:
    extractHAIRS {hairs_params} --bam sample.bam --VCF sample.vcf --out sample.unlinked.frags
    
If the data were linked-reads, invalid barcodes were filtered out:
    awk {invalid_regex.get(bc_type, "'$4 !~ /N/'")} sample.unlinked.frags > sample.frags.filt
    LinkFragments.py --bam sample.bam --VCF sample.vcf --fragments sample.frags.filt --out sample.linked.frags -d {molecule_distance}
    HAPCUT2 --fragments sample.linked.frags --vcf sample.vcf --out sample.blocks --nf 1 --error_analysis_mode 1 --call_homozygous 1 --outvcf 1 {prune} {extra}
    
If the data were not linked reads, HapCut2 ran ignoring linked-read specific information:
    HAPCUT2 --fragments sample.unlinked.frags --vcf sample.vcf --out sample.blocks --nf 1 --error_analysis_mode 1 --call_homozygous 1 --outvcf 1 {prune} {extra}
    
Variant annotation was performed using:
    bcftools annotate -a sample.phased.vcf -c CHROM,POS,FMT/GT,FMT/PS,FMT/PQ,FMT/PD -m +HAPCUT
    bcftools merge --output-type b samples.annot.bcf
'''

    def phase_bam(self):
        mol_dist    = self.PARAMETERS.get("distance-threshold", 100000)
        extra       = self.PARAMETERS.get("extra", "")
        ploidy      = self.PARAMETERS.get("ploidy", 2)
        linked      = self.PARAMETERS.get("use-linked-info", True)
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
        self.summary = f'''The harpy phase bam workflow ran using these parameters:
The provided variant file: {variantfile}

If inputs were linked-read data,the input alignments had their records filtered for valid barcodes:
    djinn sam filter-invalid --invalid sample.invalid.bam sample.bam

Phasing was performed using whatshap:
    whatshap haplotag --sample name --reference input.ref {params} {extra} input.vcf input.bam

If linked-reads, invalid-barcode alignments were added back to the phased alignments using:
    samtools merge sample.phased.bam sample.invalid.bam | samtools sort -

Alignments were sorted using:
    samtools sort sample.phased.bam
'''

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

        self.summary = f'''The harpy preprocess meier2021 workflow ran using these parameters:
Linked Read Barcode Design: Meier et al. 2021
The multiplexed input files:
    read 1: {R1}
    read 2: {R2}
    index 1: {I1}
    index 2: {I2}

Sample demultiplexing schema: {schemafile}

Samples were demultiplexed using:
    dmox --R1 --R2 --I1 --I2 {outdir} {qxrx} {unknown_barcodes} {unknown_samples}
 
QC checks were performed on demultiplexed FASTQ files using:
    falco -skip-report -skip-summary -data-filename output input.fq.gz
'''
    def preprocess_gih(self):
        me_seq   = self.PARAMETERS.get("ME-sequence", "AGATGTGTATAAGAGACAG")
        mismatch = self.PARAMETERS.get("ME-mismatch", 1)
        minlen   = self.PARAMETERS.get("min-length", 10)

        self.summary = f'''The harpy preprocess gih workflow ran using these parameters:
Linked Read Barcode Design: Iqbal et al. 2026

Input FASTQs had the ME sequence identified and removed, then provided a nucleotide padding sequence (if necessary):
    gih-stagger --me {me_seq} --max-mismatch {mismatch} --min-len {minlen} --stats output.stats FQ1 FQ2 > output.sam

The resulting interleaved unaligned SAM file was then piped into the Pheniqs for barcode demultiplexing:
    pheniqs mux --output output.sam --quality -c pheniqs.config.json --report output.json < input.sam

Finally, the nucleotides were converted into standardized haplotagging ACBD format and FASTQ format:
    gih-convert input.sam fq1 fq2 > output.stats
    
QC checks were performed on demultiplexed FASTQ files using:
    falco -skip-report -skip-summary -data-filename output input.fq.gz

Use the output 'adapters.fasta' file as input into harpy qc
'''

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

        self.summary = f'''The harpy qc workflow ran using these parameters:
fastp ran using:
    fastp  {" ".join(params)}
'''

    def snp_freebayes(self):
        ploidy 		  = self.PARAMETERS.get("ploidy", 2)
        extra 	      = self.PARAMETERS.get("extra", "")
        genomefile 	  = self.INPUTS["reference"]
        regions_input = self.INPUTS["regions"]
        groupings 	  = self.INPUTS.get("groupings", [])

        params = f"-p {ploidy} "
        params += f"--populations {groupings} " if groupings else ''
        params += extra

        self.summary = f'''The harpy snp freebayes workflow ran using these parameters:
        The provided reference genome: {genomefile}

        Genomic positions for which variants were called: {regions_input}

The freebayes parameters:
    freebayes -f REFERENCE -L samples.list -r REGION {params} |
    bcftools sort -

The variants identified in the intervals were merged into the final variant file using:
    bcftools concat -f bcf.files -a --remove-duplicates

The variants were normalized using:
    bcftools norm -m -both -d both -c w
'''
    def snp_mpileup(self):
        ploidy 		 = self.PARAMETERS.get("ploidy", 2)
        mp_extra 	 = self.PARAMETERS.get("extra", "")
        genomefile 	 = self.INPUTS["reference"]
        groupings 	 = self.INPUTS.get("groupings", [])
        region_input = self.INPUTS["regions"]

        params = f"--ploidy {ploidy} --populations "
        params += f"{groupings}" if groupings else "-"

        self.summary = f'''The harpy snp mpileup workflow ran using these parameters:
The provided reference genome: {genomefile}

Genomic positions for which variants were called: {region_input}
   
The mpileup parameters:
    bcftools mpileup --fasta-ref REFERENCE --region REGION --bam-list BAMS --annotate AD --output-type b {mp_extra}
    
The bcftools call parameters:
    bcftools call --multiallelic-caller {params} --variants-only --output-type b |
    bcftools sort -

The variants identified in the intervals were merged into the final variant file using:
    bcftools concat -f bcf.files -a --remove-duplicates"

The variants were normalized using:
    bcftools norm -m -both -d both -c w
'''

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
        if groupfile:
            grp = f"The provided populations grouping file: {groupfile}\n\n" "The alignments were concatenated using:\n\tdjinn sam concat samples.bam... > group.bam\n\n"
        else:
            grp = ""

        self.summary = f'''The harpy sv leviathan workflow ran using these parameters
The provided reference genome: {genomefile}

{grp}The barcodes were indexed using LRez
    LRez index bam -p -b INPUT

Leviathan was called using:
    LEVIATHAN -b INPUT -i INPUT.BCI -g GENOME {params}
'''

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

        self.summary = f'''The harpy sv naibr workflow ran using these parameters:
The provided reference genome: {genomefile}

naibr variant calling ran using these configurations:
    bam_file=BAMFILE
    prefix=PREFIX
    outdir=Variants/naibr/PREFIX
    {"\n\t".join([f"{k}={v}" for k,v in argdict.items()])}
'''

    def validate_bam(self) -> str:
        lr_platform = self.WORKFLOW.get("linkedreads", {}).get("type", 'none')
        self.summary = f'''The harpy validate bam workflow ran using these parameters:
Validations were performed with:
    harpy-utils check-bam {lr_platform} sample.bam > sample.txt
'''

    def validate_fastq(self):
        lr_platform = self.WORKFLOW.get("linkedreads", {}).get("type", 'none')
        self.summary = f'''The harpy validate fastq workflow ran using these parameters:
Validations were performed with:
harpy-utils check-fastq {lr_platform} sample.fastq > sample.txt
'''

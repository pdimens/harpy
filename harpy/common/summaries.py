"""Basic functions to write the workflow summaries"""
import os

class Summary:
    def __init__(self, config: dict):
        self.config: dict = config

    def get_text(self) -> str:
        return self.__getattribute__(self.config["workflow"])()

    def align_bwa(self) -> str:
        ignore_bx = self.config["linkedreads"]["type"] == "none"
        is_standardized = self.config["linkedreads"]["standardized"]
        keep_unmapped = self.config["keep_unmapped"]
        extra 		= self.config.get("extra", "") 
        genomefile 	= self.config["inputs"]["reference"]
        quality = self.config["alignment_quality"]

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
        sm = "The Snakemake workflow was called via command line:\n"
        sm += f"\t{self.config['snakemake']['relative']}"
        summary.append(sm)
        return "\n\n".join(summary)

    def align_strobe(self) -> str:
        genomefile 	= self.config["inputs"]["reference"]
        ignore_bx = self.config["linkedreads"]["type"] == "none"
        is_standardized = self.config["linkedreads"]["standardized"]
        keep_unmapped = self.config["keep_unmapped"]

        quality = self.config["alignment_quality"]
        unmapped_strobe = "" if keep_unmapped else "-U"
        unmapped = "" if keep_unmapped else "-F 4"
        bx_mode = "--barcode-tag BX" if not ignore_bx else ""
        static = "-C" if is_standardized else ""
        extra = self.config.get("extra", "") 

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
        sm += f"\t{self.config['snakemake']['relative']}"
        summary.append(sm)
        return "\n\n".join(summary)

    def assembly(self) -> str:
        # SPADES
        max_mem      = self.config["spades"]["max_memory"]
        k_param      = self.config["spades"]["k"]
        spades_extra = self.config["spades"].get("extra", "")
        # ARCS
        mapq       = self.config["tigmint"]["minimum_mapping_quality"]
        mismatch   = self.config["tigmint"]["mismatch"]
        mol_dist   = self.config["tigmint"]["molecule_distance"]
        mol_len    = self.config["tigmint"]["molecule_length"]
        span       = self.config["tigmint"]["span"]
        min_align  = self.config["arcs"]["minimum_aligned_reads"]
        min_contig = self.config["arcs"]["minimum_contig_length"]
        seq_id     = self.config["arcs"]["minimum_sequence_identity"]
        arcs_extra = self.config["arcs"].get("extra", "")
        links      = self.config["links"]["minimum_links"]
        k_param = k_param
        max_mem = max_mem // 1000
        spades_extra = spades_extra
        params = [
            f"-C scaffold",
            "-j THREADS",
            f"draft=spades",
            f"reads=interleaved",
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
        arcs += f"\tarcs-make arcs-tigmint {" ".join(params[3:])}"
        summary.append(arcs)
        sm = "The Snakemake command invoked:\n"
        sm += f"\t{self.config['snakemake']['relative']}"
        summary.append(sm)
        return "\n\n".join(summary)

    def deconvolve(self) -> str:
        kmer_length = self.config["kmer_length"]
        window_size = self.config["window_size"]
        density 	= self.config["density"] 
        dropout     = self.config["dropout"]

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
        sm += f"\t{self.config['snakemake']['relative']}"
        summary.append(sm)
        return "\n\n".join(summary)

    def demultiplex_meier2021(self) -> str:
        schemafile = self.config["inputs"]["demultiplex_schema"]
        qxrx = self.config["retain"]["qx_rx"]
        unknown_samples = self.config["retain"]["samples"]
        unknown_barcodes = self.config["retain"]["barcodes"]

        R1 = self.config["inputs"]["R1"],
        R2 = self.config["inputs"]["R2"],
        I1 = self.config["inputs"]["I1"],
        I2 = self.config["inputs"]["I2"],
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
        sm += f"\t{self.config['snakemake']['relative']}"
        summary.append(sm)
        return "\n\n".join(summary)

    def impute(self) -> str:
        print(self.config.keys())
        region = self.config.get("region", None)
        if region:
            _,positions = region.split(":")
            startpos,endpos,buffer = [int(i) for i in positions.split("-")]
        paramfiletext = "\t".join(open(self.config["inputs"]["parameters"], "r").readlines())
        summary = ["The harpy impute workflow ran using these parameters:"]
        summary.append(f"The provided variant file: {self.config["inputs"]["vcf"]}")
        preproc = "Preprocessing was performed with:\n"
        preproc += "\tbcftools view -M2 -v snps --regions CONTIG INFILE |\n"
        preproc += """\tbcftools query -i '(STRLEN(REF)==1) & (STRLEN(ALT[0])==1) & (REF!="N")' -f '%CHROM\\t%POS\\t%REF\\t%ALT\\n'"""
        summary.append(preproc)
        stitchparam = f"The STITCH parameter file: {self.config["inputs"]["parameters"]}\n"
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
        if region:
            stitch += f"\t\tregionStart = {startpos},\n"
            stitch += f"\t\tregionEnd = {endpos},\n"
            stitch += f"\t\tbuffer = {buffer},\n"
        stitch += "\t\tK = k,\n"
        stitch += "\t\tS = s,\n"
        stitch += "\t\tuse_bx_tag = usebx,\n"
        stitch += "\t\tbxTagUpperLimit = bxlimit,\n"
        stitch += "\t\tniterations = 40,\n"
        stitch += "\t\tswitchModelIteration = 39,\n"
        stitch += "\t\tsplitReadIterations = NA,\n"
        if self.config["grid_size"] > 1:
            stitch += f"\t\tgridWindowSize = {self.config["grid_size"]}\n"
        stitch += "\t\toutputdir = outdir,\n"
        stitch += "\t\toutput_filename = outfile\n\t)"
        stitchextra = "Additional STITCH parameters provided (overrides existing values above):\n"
        stitchextra += "\t" + self.config.get("stitch_extra", "None")
        summary.append(stitchextra)
        sm = "The Snakemake command invoked:\n"
        sm += f"\t{self.config['snakemake']['relative']}"
        summary.append(sm)
        return "\n\n".join(summary)

    def metassembly(self) -> str:
        bx = self.config["linkedreads"]["barcode_tag"]
        max_mem = self.config["spades"]["max_memory"]
        k_param = self.config["spades"]["k"]
        ignore_bx = self.config["spades"]["ignore_barcodes"]
        extra = self.config["spades"].get("extra", "")
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
        sm += f"\t{self.config['snakemake']['relative']}"
        summary.append(sm)
        return "\n\n".join(summary)

    def phase(self) -> str:
        bc_type           = self.config["linkedreads"]["type"]
        pruning           = self.config["phasing"]["prune"]
        map_qual          = self.config["phasing"]["min_map_quality"]
        base_qual         = self.config["phasing"]["min_base_quality"]
        molecule_distance = self.config["linkedreads"]["distance_threshold"]
        extra             = self.config.get("extra", "") 
        variantfile       = self.config["inputs"]["vcf"]["file"]
        invalid_regex = {
            "haplotagging" : "'$4 !~ /[ABCD]00/'",
            "stlfr" : "'$4 !~ /^0_|_0_|_0$/'",
            "tellseq": "'$4 !~ /N/'"
        }
        linkarg = "--10x 0" if bc_type == "none" else "--10x 1"
        indelarg   = "--indels 1 --ref reference.fasta" if self.config["inputs"].get("reference", None) else ""
        hairs_params = f"{indelarg} {linkarg} --mmq {map_qual} --mbq {base_qual} --nf 1 --maxfragments 1500000"
        prune = f"--threshold {pruning}" if pruning > 0 else "--no_prune 1"
        extra = extra

        summary = ["The harpy phase workflow ran using these parameters:"]
        summary.append(f"The provided variant file: {variantfile}")
        hetsplit = "The variant file was split by sample and filtered for heterozygous sites using:\n"
        hetsplit += "\tbcftools view -s SAMPLE | bcftools view -m 2 -M 2 -i \'GT=\"het\"\'"
        summary.append(hetsplit)
        phase = "Phasing was performed using the components of HapCut2:\n"
        phase += f"\textractHAIRS {hairs_params} --bam sample.bam --VCF sample.vcf --out sample.unlinked.frags\n"
        if bc_type != "none":
            phase += f"\t awk " + invalid_regex.get(bc_type, "'$4 !~ /N/'") + " sample.unlinked.frags > sample.frags.filt"
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
        sm = f"\t{self.config['snakemake']['relative']}"
        summary.append(sm)
        return "\n\n".join(summary)

    def qc(self) -> str:
        minlen = f"--length_required {self.config["min_len"]}"
        maxlen = f"--max_len1 {self.config["max_len"]}"
        extra = self.config.get("extra", "") 
        trim_adapters = self.config.get("trim_adapters", None)
        if trim_adapters:
            trim_arg = "--detect_adapter_for_pe" if trim_adapters == "auto" else f"--adapter_fasta {trim_adapters}"
        else:
            trim_arg = "--disable_adapter_trimming"
        dedup = "-D" if self.config["deduplicate"] else ""

        summary = ["The harpy qc workflow ran using these parameters:"]
        fastp = "fastp ran using:\n"
        fastp += f"\tfastp --trim_poly_g --cut_right " + " ".join([minlen,maxlen,trim_arg,dedup,extra])
        summary.append(fastp)
        sm = "The Snakemake command invoked:\n"
        sm += f"\t{self.config['snakemake']['relative']}"
        summary.append(sm)
        return "\n\n".join(summary)

    def simulate_snpindel(self) -> str:
        genome = self.config["inputs"]["genome"]
        snp_vcf = self.config["snp"].get("vcf", None)
        indel_vcf = self.config["indel"].get("vcf", None)
        heterozygosity = float(self.config["heterozygosity"]["ratio"])
        only_vcf = self.config["heterozygosity"]["only_vcf"]
        outprefix = self.config["prefix"]
        randomseed = self.config.get("random_seed", None)

        snp = False 
        indel = False
        variant_params = ""
        if snp_vcf or indel_vcf:
            if snp_vcf:
                snp = True
                snp_vcf_correct = snp_vcf[:-4] + ".vcf.gz" if snp_vcf.lower().endswith("bcf") else snp_vcf
                variant_params += f" -snp_vcf {snp_vcf_correct}"
            if indel_vcf:
                indel = True
                indel_vcf_correct = indel_vcf[:-4] + ".vcf.gz" if indel_vcf.lower().endswith("bcf") else indel_vcf
                variant_params += f" -indel_vcf {indel_vcf_correct}"
        else:
            snp_count = self.config["snp"].get("count", None)
            indel_count =  self.config["indel"].get("count", None)
            if snp_count:
                snp = True
                variant_params += f" -snp_count {snp_count}"
                snp_constraint = self.config["snp"].get("gene_constraints", None)
                variant_params += f" -coding_partition_for_snp_simulation {snp_constraint}" if snp_constraint else ""
                ratio = self.config["snp"].get("titv_ratio", None)
                variant_params += f" -titv_ratio {ratio}" if ratio else ""
            if indel_count:
                indel = True
                variant_params += f" -indel_count {indel_count}"
                ratio = self.config["indel"].get("indel_ratio", None)
                variant_params += f" -ins_del_ratio {ratio}" if ratio else ""
                size_alpha = self.config["indel"].get("size_alpha", None)
                variant_params += f" -indel_size_powerlaw_alpha {size_alpha}" if size_alpha else ""        
                size_constant = self.config["indel"].get("size_constant", None)
                variant_params += f" -indel_size_powerlaw_constant {size_constant}" if size_constant else ""        

            centromeres = self.config["inputs"].get("centromeres", None)
            variant_params += f" -centromere_gff {centromeres}" if centromeres else ""
            genes = self.config["inputs"].get("genes", None)
            variant_params += f" -gene_gff {genes}" if genes else ""
            exclude = self.config["inputs"].get("excluded_chromosomes", None)
            variant_params += f" -excluded_chr_list {exclude}" if exclude else ""
            variant_params += f" -seed {randomseed}" if randomseed else ""

        snp = f"-snp_vcf haplotype_X/{outprefix}.snp.hapX.vcf" if snp else ""
        indel = f"-indel_vcf haplotype_X/{outprefix}.indel.hapX.vcf" if indel else ""

        summary = ["The harpy simulate snpindel workflow ran using these parameters:"]
        summary.append(f"The provided genome: {genome}")
        summary.append(f"Heterozygosity specified: {heterozygosity}")
        haploid = "Haploid variants were simulated using simuG:\n"    
        haploid += f"\tsimuG -refseq {genome} -prefix {outprefix} {variant_params}"
        summary.append(haploid)
        if heterozygosity > 0 and not only_vcf:
            diploid = "Diploid variants were simulated after splitting by the heterozygosity ratio:\n"
            diploid += f"\tsimuG -refseq {genome} -prefix hapX {snp} {indel}"
            summary.append(diploid)
        sm = "The Snakemake command invoked:\n"
        sm += f"\t{self.config['snakemake']['relative']}"
        summary.append(sm)
        return "\n\n".join(summary)
    
    def simulate_inversion(self) -> str:
        return self.simulate_variants()

    def simulate_cnv(self) -> str:
        return self.simulate_variants()      

    def simulate_translocation(self) -> str:
        return self.simulate_variants()

    def simulate_variants(self) -> str:
        variant = self.config["workflow"].split("_")[-1]
        outprefix = self.config["prefix"]
        genome = self.config["inputs"]["genome"]
        vcf = self.config[variant].get("vcf", None)
        heterozygosity = float(self.config["heterozygosity"]["ratio"])
        only_vcf = self.config["heterozygosity"]["only_vcf"]
        randomseed = self.config.get("random_seed", None)

        if vcf:
            variant_params = f"-{variant}_vcf {variant}.vcf"
        else:
            variant_params = f"-{variant}_count " + str(self.config[variant]["count"])
            centromeres = self.config["inputs"].get("centromeres", None)
            variant_params += f" -centromere_gff {centromeres}" if centromeres else ""
            genes = self.config["inputs"].get("genes", None)
            variant_params += f" -gene_gff {genes}" if genes else ""
            exclude = self.config["inputs"].get("excluded_chromosomes", None)
            variant_params += f" -excluded_chr_list {exclude}" if exclude else ""
            variant_params += f" -seed {randomseed}" if randomseed else ""
            if variant in ["inversion", "cnv"]:  
                variant_params += f" -{variant}_min_size " +  str(self.config[variant]["min_size"])
                variant_params += f" -{variant}_max_size " +  str(self.config[variant]["max_size"])
            if variant == "cnv":
                variant_params += f" -duplication_tandem_dispersed_ratio " +  str(self.config[variant]["duplication_ratio"])
                variant_params += f" --cnv_max_copy_number " +  str(self.config[variant]["max_copy"])
                variant_params += f" --cnv_gain_loss_ratio " +  str(self.config[variant]["gain_ratio"])

        summary = [f"The harpy simulate {variant} workflow ran using these parameters:"]
        summary.append(f"The provided genome: {genome}")
        summary.append(f"Heterozygosity specified: {heterozygosity}")
        haploid = "Haploid variants were simulated using simuG:\n"    
        haploid += f"\tsimuG -refseq {genome} -prefix {outprefix} {variant_params}"
        summary.append(haploid)
        if heterozygosity > 0 and not only_vcf:
            diploid = f"Diploid variants were simulated after splitting by the heterozygosity ratio:\n"
            diploid += f"\tsimuG -refseq {genome} -prefix HAP_PREFIX {variant}.vcf hapX.vcf"
            summary.append(diploid)
        sm = "The Snakemake command invoked:\n"
        sm += f"\t{self.config['snakemake']['relative']}"
        summary.append(sm)
        return "\n\n".join(summary)
    
    def snp_freebayes(self) -> str:
        ploidy 		= self.config["ploidy"]
        extra 	    = self.config.get("extra", "") 
        regions_input = self.config["inputs"]["regions"]
        genomefile 	= os.path.basename(self.config["inputs"]["reference"])
        groupings 	= self.config["inputs"].get("groupings", None)

        params = f"-p {ploidy} "
        params += f"--populations {groupings} " if groupings else ''
        params += extra

        summary = ["The harpy snp freebayes workflow ran using these parameters:"]
        summary.append(f"The provided reference genome: {genomefile}")
        summary.append(f"Genomic positions for which variants were called: {regions_input}")
        varcall = "The freebayes parameters:\n"
        varcall += f"\tfreebayes -f REFERENCE -L samples.list -r REGION {params} |\n"
        varcall += f"\tbcftools sort -"
        summary.append(varcall)
        merged = "The variants identified in the intervals were merged into the final variant file using:\n"
        merged += "\tbcftools concat -f bcf.files -a --remove-duplicates"
        summary.append(merged)
        normalize = "The variants were normalized using:\n"
        normalize += "\tbcftools norm -m -both -d both -c w"
        summary.append(normalize)
        sm = "The Snakemake command invoked:\n"
        sm += f"\t{self.config['snakemake']['relative']}"
        summary.append(sm)
        return "\n\n".join(summary)
    
    def snp_mpileup(self) -> str:
        mp_extra = self.config.get("extra", "")
        genomefile = self.config["inputs"]["reference"]
        groupings = self.config["inputs"].get("groupings", [])
        region_input = self.config["inputs"]["regions"]
        ploidy = self.config["ploidy"]
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
        sm += f"\t{self.config['snakemake']['relative']}"
        summary.append(sm)
        return "\n\n".join(summary)

    def validate_bam(self) -> str:
        lr_platform = self.config["linkedreads"]["type"]

        summary = ["The harpy validate bam workflow ran using these parameters:"]
        valids = "Validations were performed with:\n"
        valids += f"\tcheck_bam {lr_platform} sample.bam > sample.txt"
        summary.append(valids)
        sm = "The Snakemake command invoked:\n"
        sm += f"\t{self.config['snakemake']['relative']}"
        summary.append(sm)
        return "\n\n".join(summary)

    def validate_fastq(self) -> str:
        lr_platform = self.config["linkedreads"]["type"]
        
        summary = ["The harpy validate fastq workflow ran using these parameters:"]
        valids = "Validations were performed with:\n"
        valids += f"\tcheck_fastq {lr_platform} sample.fastq > sample.txt"
        summary.append(valids)
        sm = "The Snakemake command invoked:\n"
        sm += f"\t{self.config['snakemake']['relative']}"
        summary.append(sm)
        return "\n\n".join(summary)
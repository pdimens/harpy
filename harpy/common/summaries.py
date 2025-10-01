"""Basic functions to write the workflow summaries"""

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
        sm = "The Snakemake command that was used:\n"
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
        sm = "The Snakemake command that was used:\n"
        sm += f"\t{self.config['snakemake']['relative']}"
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
        sm = "The Snakemake command that was used:\n"
        sm += f"\t{self.config['snakemake']['relative']}"
        summary.append(sm)
        return "\n\n".join(summary)

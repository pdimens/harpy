
@click.command(no_args_is_help = True, context_settings=dict(allow_interspersed_args=False), epilog = "Documentation: https://pdimens.github.io/harpy/workflows/align/ema")
@click.option('-x', '--extra-params', type = EmaParams(), help = 'Additional ema align parameters, in quotes')
@click.option('-d', '--fragment-density',  is_flag = True, show_default = True, default = False, help = 'Perform read fragment density optimization')
@click.option('-w', '--depth-window', default = 50000, show_default = True, type = click.IntRange(min = 50), help = 'Interval size (in bp) for depth stats')
@click.option('-b', '--ema-bins', default = 500, show_default = True, type = click.IntRange(1,1000, clamp = True), help="Number of barcode bins")
@click.option('-u', '--keep-unmapped',  is_flag = True, default = False, help = 'Retain unmapped sequences in the output')
@click.option('-q', '--min-quality', default = 30, show_default = True, type = click.IntRange(0, 40, clamp = True), help = 'Minimum mapping quality to pass filtering')
@click.option('-o', '--output-dir', type = click.Path(exists = False, resolve_path = True), default = "Align/ema", show_default=True,  help = 'Output directory name')
@click.option('-p', '--platform', type = click.Choice(['haplotagging', '10x'], case_sensitive=False), default = "haplotagging", show_default=True, help = "Linked read type\n[haplotagging, 10x]")
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(4,999, clamp = True), help = 'Number of threads to use')
@click.option('-l', '--barcode-list', type = click.Path(exists=True, dir_okay=False, resolve_path=True), help = "File of known barcodes for 10x linked reads")
@click.option('--setup-only',  is_flag = True, hidden = True, default = False, help = 'Setup the workflow and exit')
@click.option('--contigs',  type = ContigList(), help = 'File or list of contigs to plot')
@click.option('--container',  is_flag = True, default = False, help = 'Use a container instead of conda', callback=container_ok)
@click.option('--hpc',  type = HPCProfile(), help = 'HPC submission YAML configuration file')
@click.option('--quiet', show_default = True, default = 0, type = click.Choice([0, 1, 2]), help = '`0` all output, `1` show one progress bar, `2` no output')
@click.option('--skip-reports',  is_flag = True, show_default = True, default = False, help = 'Don\'t generate HTML reports')
@click.option('--snakemake', type = SnakemakeParams(), help = 'Additional Snakemake parameters, in quotes')
@click.argument('reference', type=InputFile("fasta", gzip_ok = True), required = True, nargs = 1)
@click.argument('inputs', required=True, type=click.Path(exists=True, readable=True, resolve_path=True), nargs=-1)
def ema(reference, inputs, output_dir, platform, barcode_list, fragment_density, depth_window, keep_unmapped, threads, ema_bins, skip_reports, extra_params, min_quality, snakemake, quiet, hpc, container, contigs, setup_only):
    """
    Align sequences to reference genome using EMA

    Provide the reference fasta followed by the fastq files and/or directories at the end of the
    command as individual files/folders, using shell wildcards
    (e.g. `data/axolotl*.fastq.gz`), or both.

    EMA may improve mapping, but it also marks split reads as secondary
    reads, making it less useful for variant calling with leviathan. The barcode
    list is a file of known barcodes (in nucleotide format, one per line) that lets EMA know what
    sequences at the beginning of the forward reads are known barcodes.
    """
    workflow = Workflow("align_ema", "align_ema.smk", output_dir, quiet)
    workflow.setup_snakemake(container, threads, hpc, snakemake)
    workflow.reports = ["align_stats.qmd", "align_bxstats.qmd"]
    workflow.conda = ["align", "r", "qc"]

    ## checks and validations ##
    platform = platform.lower()
    # the tellseq stuff isn't impremented yet, but this is a placeholder for that (wishful thinking)
    if platform in ["tellseq", "10x"] and not barcode_list:
        print_error("missing barcode list", f"{platform} technology requires a list of known barcodes.")
        if platform == "10x":
            print_solution("Running EMA requires 10X barcodes provided to [green]--barcode-list[/]. A standard 10X barcode list can be downloaded from [dim]https://github.com/10XGenomics/cellranger/tree/master/lib/python/cellranger/barcodes[/dim]")
        else:
            print_solution("Running EMA requires TELLseq barcodes provided to [green]--barcode-list[/]. They can be acquired from the TELL-read software [dim]https://www.illumina.com/products/by-type/informatics-products/basespace-sequence-hub/apps/universal-sequencing-tell-seq-data-analysis-pipeline.html[/dim]")
        sys.exit(1)
    if platform == "haplotagging" and barcode_list and not quiet:
        print_notice("Haplotagging data does not require a barcode list and the file provided to [green]--barcode-list[/] will be ignored.")
    fqlist, sample_count = parse_fastq_inputs(inputs, "INPUTS")
    check_fasta(reference)
    if contigs:
        fasta_contig_match(contigs, reference)
    if barcode_list:
        validate_barcodefile(barcode_list, False, quiet, gzip_ok=False, haplotag_only=True)

    workflow.config = {
        "workflow" : workflow.name,
        "alignment_quality" : min_quality,
        "keep_unmapped" : keep_unmapped,
        "fragment_density_optimization": fragment_density,
        "depth_windowsize" : depth_window,
        "platform" : platform,
        "EMA_bins" : ema_bins,
        **({'extra': extra_params} if extra_params else {}),
        "snakemake" : {
            "log" : workflow.snakemake_log,
            "absolute": workflow.snakemake_cmd_absolute,
            "relative": workflow.snakemake_cmd_relative,
        },
        "conda_environments" : workflow.conda,
        "reports" : {
            "skip": skip_reports,
            **({'plot_contigs': contigs} if contigs else {'plot_contigs': "default"}),
        },
        "inputs" : {
            "reference": reference,
            **({'barcode_list': barcode_list} if barcode_list else {}),
            "fastq": fqlist
        }
    }

    workflow.start_text = workflow_info(
        ("Samples:",sample_count),
        ("Reference:", os.path.basename(reference)),
        ("Platform:", platform),
        ("Output Folder:", os.path.basename(output_dir) + "/")
    )

    workflow.initialize(setup_only)


read_with_bc = config["read_with_barcode"]

rule pheniqs_demux:
    input:
        FQ1 = if read_with_bc == 1 else ,
        FQ2 = if read_with_bc == 1 else ,
        pheniqs_conf = "workflow/pheniqs_config.json"
    output:
        FQ1 = "pheniqs/{sample}_extract_R1.bam",
        FQ2 = "pheniqs/{sample}_extract_R2.bam"
    log:
        "logs/pheniqs/{sample}.json"
    shell:
        """
        # If barcodes in R1 header, then use R1 in extract call, else use R2
        if [[ $params.file_containing_barcodes == "R1" ]]; then
            input_fastq_file=`echo "$reads" | cut -d " " -f 1`
            read2_in=`echo "$reads" | cut -d " " -f 2`

            output_fastq_file="\${prefix}"_extract_R1.bam
            read2_out="\${prefix}"_extract_R2.bam
        else
            input_fastq_file=`echo "$reads" | cut -d " " -f 2`
            read2_in=`echo "$reads" | cut -d " " -f 1`
        fi

        pheniqs mux \\
            --input {input.FQ1} \\
            --input {input.FQ2} \\
            --output {output.FQ1} \\
            --output {output.FQ2} \\
            -c {input.pheniqs_conf} \\
            --quality \\
            --report {log}
        """

rule determine_stagger:
    input:
        FQ1 = ,
        FQ2 = ,
        summary = 
    shell:
        """
        prefix=`echo $reads | cut -d "$params.delimiter" -f$params.field`

        input_R1_file=`echo $reads | cut -d " " -f1`
        input_R2_file=`echo $reads | cut -d " " -f2`
        info_summary_file=`echo "$reads" | cut -d " " -f 3`

        # Get the info summary file and determine whether stagger sequence is being used
        # ln -s $params.findMESequence_output_dir/\${prefix}_info_summary.txt .
        # include_stagger=\$(python $baseDir/stagger_check.py "\${prefix}_info_summary.txt")
        include_stagger=\$(python $baseDir/stagger_check.py \$info_summary_file)

        # Determine bead complexity
        ln -s $baseDir/ .
        samtools view {input.FQ1} | awk 'NR>100000 && NR<=2000000 {n=split(\$0,arr,"CB:Z:"); if(n>1){split(arr[2],result,"\t"); print result[1];}}' > temp_barcodes.txt
        complexity=\$(python $baseDir/determine_bead_complexity.py --include_stagger "\$include_stagger" --temp_barcodes temp_barcodes.txt)
        grep "\$complexity" 12nt-barcodes-with-stagger-and-bead-complexity.txt > 12nt-barcodes-subset.txt

        # Convert 43nt barcode sequences to ‘haplotagging’ format AxxCxxBxxDxx for both R1 and R2 files.
        awk -v stagger="\$include_stagger" -f $baseDir/convert_bam_to_ACBD_format.awk 12nt-barcodes-subset.txt <(samtools view {input.FQ1}) | pigz -p 8 > \${prefix}_converted.R1.fastq.gz
        mv barcode_log.log \${prefix}_converted_R1_barcode_log.log
        awk -v stagger="\$include_stagger" -f $baseDir/convert_bam_to_ACBD_format.awk 12nt-barcodes-subset.txt <(samtools view {input.FQ1}) | pigz -p 8 > \${prefix}_converted.R2.fastq.gz
        mv barcode_log.log \${prefix}_converted_R2_barcode_log.log

        # Produce QC report
        python $baseDir/knee_and_segment_count_plots.py --input "\${prefix}_converted_R1_barcode_log.log" --add_cutoff_line $add_cutoff_line
        python $baseDir/knee_and_segment_count_plots.py --input "\${prefix}_converted_R2_barcode_log.log" --add_cutoff_line $add_cutoff_line

        multiqc . -f --filename \${prefix}.multiqc.report
        """
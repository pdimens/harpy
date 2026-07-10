#! /usr/bin/env bash

PREFIX=$1
FQ1=$2
FQ2=$3

cutadapt -a "CTGTCTCTTATACACATCT" -A "AGATCGGAAGAGC" \
  --nextseq-trim=20 \
  --minimum-length 50 \
  --too-short-output ${PREFIX}_short.R1.fq \
  --too-short-paired-output ${PREFIX}_short.R2.fq \
  --cores 10 --match-read-wildcards --action trim \
  --output ${PREFIX}_trim.R1.fq.gz \
  --paired-output ${PREFIX}_trim.R2.fq.gz \
  ${FQ1} ${FQ2} > ${PREFIX}_precheck-inserts.log



### part 2
cutadapt -a CTGTCTCTTATACACATCT \\
  -A CTGTCTCTTATACACATCT \\
  --nextseq-trim 20 \\
  -m 50 \\
  --cores 10 \\
  --match-read-wildcards \\
  -o \${prefix}_trimmed_R1.fastq.gz \\
  -p \${prefix}_trimmed_R2.fastq.gz \\
  \$hUMI_R1_fastq_file \\
  \$hUMI_R2_fastq_file \\
  > \${prefix}_trim3pAd.log

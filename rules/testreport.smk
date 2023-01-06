
rule all:
    input: "Variants/mpileup/variants.raw.stats"
    output: "Variants/report.html"
    message: "Generating bcftools report"
    script: "../utilities/bcftoolsreport.Rmd"

rule all:
    input: "VariantCall/variants.raw.stats"
    output: "VariantCall/report.html"
    message: "Generating bcftools report"
    script: "../utilities/bcftoolsreport.rmd"
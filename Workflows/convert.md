---
label: Convert
description: Convert between linked-read data formats
icon: arrow-switch
order: 10
---

# :icon-arrow-switch: Convert between data formats
[!badge variant="secondary" text="linked reads"]

Regrettably, the bright minds who developed various linked-read technologies cannot seem to agree on a unified data format.
That's annoying at best and hinders the field of linked-read analysis at worst, as there are pieces of very clever software
that are specific to a narrow set of linked-read data formats. Until such a day where there is a concensus, Harpy provides
the means to convert between the various popular linked-read data formats. 

==- Creating order in chaos :warning:
You will notice one of the formats is called
`standard`, and that is our attempt to _encourage_ linked-read practioners to consider a unified and practical data format
in which a barcode **of any format** is encoded is in the `BX:Z` tag of a FASTQ/BAM file (a standard SAM-compliant tag) and the validation for the barcode
(whether it is valid or not according to the technology) is encoded in the `VX:i` tag as either `0` (invalid) or `1` (valid).

As an example, if you have stLFR barcoded data, whose barcodes take the form `1_2_3`, the barcode `54_0_1123` would be considered
invalid because stLFR barcodes with a `0` as one of the segments are invalid (missing/ambiguous segment). The `standard` data format,
regardless of FASTQ or BAM, would have the barcoded as `BX:Z:54_0_1123` and the validation as `VX:i:0`.
===


:::info Useless trivia
This module was written while Pavel was waiting at a mechanic shop for his car to be repaired. During development,
it was called `lr-switcheroo`.
:::
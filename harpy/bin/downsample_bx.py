import pysam
import random

input_bam = "input.bam"
output_bam = "downsampled.bam"
tag_key = "XT"  # Replace with your tag key
tag_value = "desired_value"  # Replace with your tag value
downsample_fraction = 0.1  # Fraction of reads to retain

with pysam.AlignmentFile(input_bam, "rb") as infile, pysam.AlignmentFile(output_bam, "wb", template=infile) as outfile:
    for read in infile:
        if read.has_tag(tag_key) and read.get_tag(tag_key) == tag_value:
            if random.random() < downsample_fraction:
                outfile.write(read)
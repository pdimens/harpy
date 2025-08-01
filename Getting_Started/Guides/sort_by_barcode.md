---
author: 
  - name: Pavel Dimens
    avatar: https://cals.cornell.edu/sites/default/files/styles/faculty/public/2024-09/afs-headshot-high-res-2cropped_0.jpg
date: 2024-11-05
category: guides
description: Sorting data by linked-read barcode
icon: list-ordered
image: https://t4.ftcdn.net/jpg/07/40/71/29/360_F_740712956_dyX3F3ehpbdqjCZT0RdqkJ8Aniu7fecl.jpg
---

# :icon-list-ordered: Sort data by barcode
You would think sorting data would be a no-brainer, and in most cases it is.
You can use `seqtk` or `seqkit` to sort FASTQ/A files by their IDs, `samtools` to sort
SAM/BAM/CRAM files by name or coordinates. However, in the world of linked-read
data, sometimes you may need to sort your FASTQ (or BAM) files by the
linked-read barcode. The way to do that wasn't initially obvious to the Harpy/haplotag
team, so this article serves to make this knowledge widely available to linked-read
adopters.

## Sorting Alignments
Let's start with BAM (or SAM/CRAM) files because the process is much simpler.
Since the linked-read barcode is stored in a `BX:Z` tag (or less often as `BC:Z:`),
we can use a little feature of `samtools sort` to guide the sort by the barcode:
> -t TAG     
> Sort first by the value in the alignment tag TAG, then by position or name (if using -n or -N)

The `-t` option then makes it pretty trivial to sort an alignment file by barcode:
```bash
samtools sort -t BX file.bam > sorted.bam
```

The above command will accomplish sorting by whatever kind of barcode is listed in
the `BX:Z` tag. If your barcode was in the `BC:Z` tag, you would use `-t BC`.

!!! BX and not BX:Z
Notice that you only specify `BX` rather than `BX:Z`, because the `:Z` part
of the tag is a SAM specification to indicate what type of data is in the tag (i.e. `TAG:TYPE:VALUE`).
A tag is named by its first two letters and you cannot have clashing names in a single record
(e.g. having both `BX:Z` and `BX:i` would be invalid because they are both named `BX`).
!!!

## Sorting FASTQ
Sorting FASTQ files by barcode is trickier, only because there aren't (to our knowledge!)
any existing convenience methods to do it. Like any bioinformatics puzzle, you could
probably solve it with a sophisticated AWK command, but HTSlib tools are so much more
efficient and built for these exact purposes. The process to accomplish this includes
3 steps that will be shown at the end as a single pipe.

### 1. convert FASTQ to SAM
Yep, we're solving our problem by doing a simple file conversion to SAM/BAM. That's
the easiest way to do it, surprisingly. FASTQ files can be converted to unmapped BAM
files using `samtools import`, which would also interleave the forward and reverse
reads into a single file. The `-T "*"` argument preserves all the tags between file formats.
```bash
samtools import -T "*" sample_01.R1.fq sample_01.R2.fq > sample_01.sam 
```

### 2. sort the SAM by barcode
Exactly like shown above to sort a SAM/BAM file with `samtools sort`, we're going
to do the same on the unmapped SAM file we just created:
```bash
samtools sort -O SAM -t BX sample_01.sam > sample_01.sort.sam
```

### 3. convert SAM back to FASTQ
Now that the data have been sorted, we need to convert it back into forward and reverse
FASTQ files using `samtools fastq`. The `-T "*"` argument once again preserves all the tags
between file formats. The `-1` and `-2` arguments are the forward and reverse output FASTQ files,
respectively.
```bash
samtools fastq -T "*" -1 sample_01.sort.R1.fq -2 sample_01.sort.R2.fq sample_01.sort.sam
```

### as a single pipe
Rather than splitting out these three processess, you can stream/pipe them in a single workflow:
```bash
samtools import -T "*" sample_01.R1.fq sample_01.R2.fq |
samtools sort -O SAM -t BX |
samtools fastq -T "*" -1 sample_01.sort.R1.fq -2 sample_01.sort.R2.fq
```
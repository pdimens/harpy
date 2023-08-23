---
label: Common Issues
icon: alert
order: 1
---

Lots of stuff can go wrong during an analysis. The intent of this page is to highlight
common issues you may experience during analysis and ways to address these issues.

## Failures during imputation or phasing
If you use `bamutils clipOverlap` on alignments that are used for the `impute` or
`phase` modules, they will cause both programs to error. We don't know why, but they do.

**Solution**: Do not clip overlapping alignments for bam files you intend to use for
the `impute` or `phase` modules. Harpy does not clip overlapping alignments, so
alignments produced by Harpy should work just fine.

## Alignment file name and ID: tag mismatch
Aligning a sample to a genome via Harpy will insert the sample name (based on the file name)
into the alignment header (the `@RG ID:name SM:name` tag). It likewise expects, through various steps,
that the sample names in resulting vcf files match the filenames of associated bam files. This creates 
problems when manually renaming alignment files after the creation of any vcf files. If you rename the 
bam file, the alignments will still have the original sample name hardcoded into the file header. 
Harpy will check for this and will preemtively warn you of a mismatch between file name and encoded
sample name. Due to certain expectations of the workflow, this mismatch will absolutely cause things
to fail, hence the pre-flight check.

**Solution**: If you need to rename a bam file, do so using the `renamebam` script bundled with Harpy, which is a just a thin veneer over `samtools addreplacerg` with some extra validations.
```bash
renamebam input.bam newname
```
Call the script with no arguments to see the full usage instructions.


**More cases will be added here as they become apparent to us**
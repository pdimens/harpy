# Fixes
## Reports
- preprocess GIH uses log-scaled binning to reduce data significantly and prevent altair crashes
- align aggregate linked read Total vs. plot uses reversed tubro palette and a dropdown selector instead of radio
- samtools stats stats-boxes sets more realistic cutoffs for %mapped, %properly paired, %optical dupes 
- coerce bins in bcftools report into numeric type
- force columns pertaining to contigs into String types
    - fixes altair compatability for nominal data
    - fixes table parsing errors that arise when contigs are both named as a number and as strings (e.g. a contig named '1' and another named 'X')
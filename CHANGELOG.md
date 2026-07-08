# Fixes
## Reports
- preprocess GIH uses log-scaled binning to reduce data significantly and prevent altair crashes
- align aggregate linked read Total vs. plot uses reversed tubro palette and a dropdown selector instead of radio
- samtools stats stats-boxes sets more realistic cutoffs for %mapped, %properly paired, %optical dupes 
- coerce bins in bcftools report into numeric type
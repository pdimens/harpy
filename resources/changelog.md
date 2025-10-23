# new
- `diagnose` now has 2 subcommands:
  - `stall`: same as previous `diagnose` behavior, where it runs snakemake with `--dry-run --debug-dag`
  - `rule`: attempt to directly run the failing rule of a workflow as identified in the snakemake log, will attempt to run snakemake to generate missing inputs if necessary
- `harpy resume` has new `--direct` option to call Snakemake directly without harpy intervention

# deprecations
- harpy convert
- harpy downsample
- harpy simulate linkedreads

# changes
- simplified the rich-click theming
- bwa-mem2 replaces bwa in align bwa
- impute needs a minimum of 5 biallelic snps per contig
- alignment during metassembly no longer outputs unmapped reads or alignments with mapq < 10
- added hidden `--force` option to metassembly to force athena to run even if fq/bam don't pass its internal qc
- options error borders are now yellow, making it consistent with other errors

# fixes
- impute workflow has explicit output plot filename declarations to catch errors better

# In progress but not ready
- replace manually hacked progressbars with executor plugin
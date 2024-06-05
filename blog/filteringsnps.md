---
author: Pavel Dimens
date: 2024-06-05
category: guides
description: A gentle introduction to the wild world of filtering SNPs
icon: filter
image: https://vre.eucanshare.bsc.es/vre/tools/SNP_Filtering/assets/home/logo.png
---

# :icon-filter: Filtering Variants
The discussion around filtering SNPs and indels is _massive_ and many researchers go about it differently, each very
opinionated as to why their method is the best. As a starting point, have a look at how the authors of `HTSlib` give [a
technical overview of variant filtering](http://www.htslib.org/workflow/filter.html). It's a dense read, but does offer
insights and considerations for SNP/indel filtering. Here are some of the basic things to be mindful of for variant filtering:

==- using bcftools to filter
The best and fastest way to filter variants is to use [bcftools](https://samtools.github.io/bcftools/bcftools.html#expressions),
which has a bit of a learning curve, but its power is unmatched. Filtering can be achieved using either `bcftools view` or `bcftools filter`
and the filtering expression can either be `-i` to **include** sites or `-e` to **exclude** sites matching the expression: 
```bash
# bcftools view approach
bcftools view -i 'EXPRESSION' input.vcf > output.vcf
bcftools view -e 'EXPRESSION' input.vcf > output.vcf

# bcftools filter approach
bcftools filter -i 'EXPRESSION' input.vcf > output.vcf
bcftools filter -e 'EXPRESSION' input.vcf > output.vcf
```
In either case, you can add `-Ob` to output a compressed `bcf` (recommended) file instead of an uncompressed `vcf` file (default). The
[EXPRESSION](https://samtools.github.io/bcftools/bcftools.html#expressions) is extremely flexible and multiple expressions can be chained
with `||` (OR) and `&&` (AND).
```bash filtering expression examples
# -e to EXCLUDE
bcftools view -Ob -e 'QUAL <= 10 || DP > 35 || MQBZ < -3 || RPBZ < -3 || RPBZ > 3 || FORMAT/SP > 32 || SCBZ > 3' in.vcf > out.bcf

# -i to INCLUDE, this example would result in the same output as the -e example
bcftools filter -Ob -i 'QUAL > 10 || DP <= 35 || MQBZ >= -3 || RPBZ >= -3 || RPBZ <= 3 || FORMAT/SP <= 32 || SCBZ <= 3' in.vcf > out.bcf
```
===

#### genotype quality (QUAL)
You will obviously want higher quality genotype calls to remove false positives. The HTSlib guide suggests at least `50` (e.g. `-i 'QUAL>=50'`),
but we typically filter much higher at `90` or more (e.g. `-i 'QUAL>=90'`).

#### read depth (DP)
Variant sites with too few reads backing up the genotype might be false positives, although this may not hold true for very
low-coverage data. Conversely, a maximum cut off is important because sites with very high read depths (relative to the distribution of read depth)
are likely repetitive ones mapping to multiple parts of the genome. You could used fixed values for these thresholds that make sense for your data.
One scalable approach is to define the thresholds as quantiles, such as the `0.01` and `0.99` quantiles of read depths, which would remove the
sites with the lowest 1% and highest 1% read depths. These are example quantiles and they don't need to be symmetric. It would be best to
plot the distribution of site depths to assess what makes sense for your data. Unfortunately, `bcftools` does not have internal routines to calculate
quantiles, but you can do it all from the command line using a combination of `bcftools query` and `awk` (separated onto 3 lines here for demonstration purposes):
```bash # find a specific depth quantile
bcftools query -f "%DP\n" input.bcf |\
  sort -n |\
  awk '{all[NR] = $0} END{print all[int(NR*0.95 - 0.5)]}'
```
- **line 1**: extract the depth for every site in the vcf file, one per line
- **line 2**: sort the values numerically
- **line 3**: find the 95th quantile
    - change the `0.95` in `NR*0.95` to whatever quantile you want
    - the `- 0.5` part rounds down and may need to be adjusted for your quantile

#### minor allele frequency (MAF)
It's usually advisable to set a minor allele frequency threshold with which to remove sites below that threshold. The reasoning
is that if a MAF is too low, it might be because of incorrectly called genotypes in a very small handful of individuals (e.g. one or two). The MAF threshold is again dependent on your data, although it's
fairly common to use `0.05` (e.g. `-i 'MAF>0.05'`) to `0.10` (e.g. `-i 'MAF>0.10'`).

#### missing data (F_MISSING)
Missing data is, frankly, not terribly useful. The amount of missing data you're willing to tolerate will depend on your study, but
it's common to remove sites with >20% missing data (e.g. `-e 'F_MISSING>0.2'`). This can be as strict (or lenient) as you want; it's not uncommon to see very
conservative filtering at 10% or 5% missing data. **However**, you can impute missing genotypes to recover
missing data! Harpy can leverage linked-read information to impute genotypes with the [!badge corners="pill" text="impute"](../Modules/impute.md)
module. You should try to impute genotypes first before filtering out sites based on missingness.

---
author: Pavel Dimens
date: 2024-09-30
category: guides
description: A realistic workflow to simulate variants
icon: git-merge-queue
image: https://neupsykey.com/wp-content/uploads/2017/03/ch02_f01.jpg
---

# :icon-git-merge-queue: Simulating variants
You may want to (and are encouraged to) simulate data before investing in the
costs associated with linked-read sample preparation and subsequent sequencing. 
Harpy provides both a variant and linked-read simulators and this tutorial serves to
show a real-world workflow starting with a haploid genome and creating a diploid genome
with the variants we want in it, which will then be fed into linked-read simulation. The
process might seem a little roundabout due to the limitations of the underlying software,
but it shouldn't be too bad to wrap your head around it! Ultimately, you will create
linked-reads from the resulting genome and then aligning those reads onto your **original**
genome to identify those variants.

===  :icon-checklist: You will need
- a genome assembly in FASTA format: [!badge variant="success" text=".fasta"] [!badge variant="success" text=".fa"] [!badge variant="success" text=".fasta.gz"] [!badge variant="success" text=".fa.gz"]
!!! Picking a genome
For simplicity and shorter runtimes, this tutorial can be followed with a [Drosophila melanogaster](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001215.4/)
genome. Otherwise use your favorite genome! For learning purposes, the fewer contigs the better.
!!!
===

## 1. Add snps and indels
To keep this tutorial simple, but with a certain amount of real-world complexity, let's say you are interested in **inversions** and
whether your linked-reads will be able to identify inversions in your system. To simulate inversions, we will first simulate SNPs and indels,
then simulate inversions onto that. The SNPs and indels serve to create typical variants in the downstream linked-reads because it's highly
unlikely that the _only_ variants in your data are a few inversions. 

!!! heterozygosity
By specifying a heterozygosity value for `harpy simulate ...`, we can make sure that our diploid haplotypes aren't exclusively homozygous for the alternative allele of the variants we are introducing.
!!!

For demonstrative purposes, let's say I wanted to simulate SNPs and indels like so:
- 500 indels
- 75k snps
- introduce a heterozygosity of 10% for the variants
  - based on https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1203202/pdf/255.pdf

The Harpy command to accomplish this is:
```bash
harpy simulate snpindel -m 500 -n 75000 -z 0.1 --conda -o sim_random_snpindel GENOME.fa
```
==- :icon-terminal: code explation
- `-m` is the number of ranbom indels (`500`)
- `-n` is the number of random snps (`500`)
- `-z` is the level of heterozygosity (`0.1` = 10%)
- `-o` is the name of the output directory
    - specifying this so subsequent runs don't overwrite each other
- `--conda` is optional and a matter of runtime preference
- `GENOME.fa` is your genome file
===

### 📝 output to keep
The two VCFs of snps and indels

## 2. Add random inversions
Next, we will need to do pretty much the same thing but for inversions. If you wanted to manually
create inversions in specific areas or with specific lengths, this would be a good starting point too since
you could manually modify the resulting VCF to create the specific inversions you want. We won't be covering
that here, but you should hopefully be able to intuit how to do that by the end of this tutorial.

```bash
harpy simulate inversion --conda -n 20 -z 0.1 --min-size 30000 dmel.nosex.fa
```
==- :icon-terminal: code explation
- `-n` is the number of random inversions (`20`)
- `--min-size` is the minimum inversion size (`30000`)
    - default is 1000 bp, which was arbitrarily made bigger in this example
- `-z` is the level of heterozygosity (`0.1` = 10%)
- `-o` is the name of the output directory
- `--conda` is optional and a matter of runtime preference
- `GENOME.fa` is your genome file
===

### 📝 output to keep
The two inversion VCFs

---
## ⏸️ Checkpoint ⏸️
So it seems like we did a bunch of simulating already, but we aren't done just yet. What we have done so far is create
random variants (SNPs, indels, inversions), and now we need to create a diploid genome from them. Armed with the VCF files
output from [**Step 1**](#1-add-snps-and-indels) and [**Step 2**](#2-add-random-inversions), we can now build our genome using "known" variants.

---

## 3. Create diploid genome from "known" snps and indels
Using the VCFs from [**Step 1**](#1-add-snps-and-indels), we will run Harpy twice, once for each haplotype, using the corresponding VCFs:

### haplotype 1
```bash
harpy simulate snpindel --conda --snp-vcf SNP.hap1.vcf --indel-vcf indel.hap1.vcf -o sim_snp_hap1 GENOME.fa
```

==- :icon-terminal: code explation
- `--snp-vcf` is the vcf of snps for haplotype 1 from [**Step 1**](#1-add-snps-and-indels)
- `--indel-vcf` is the vcf of indels for haplotype 1 from [**Step 1**](#1-add-snps-and-indels)
- `-o` is the name of the output directory
- `--conda` is optional and a matter of runtime preference
- `GENOME.fa` is your genome file
===

### haplotype 2
```bash
harpy simulate snpindel --conda --snp-vcf SNP.hap2.vcf --indel-vcf indel.hap2.vcf -o sim_snp_hap2 GENOME.fa
```
==- :icon-terminal: code explation
- `--snp-vcf` is the vcf of snps for haplotype 2 from [**Step 1**](#1-add-snps-and-indels)
- `--indel-vcf` is the vcf of indels for haplotype 2 from [**Step 1**](#1-add-snps-and-indels)
- `-o` is the name of the output directory
- `--conda` is optional and a matter of runtime preference
- `GENOME.fa` is your genome file
===

### 📝 output to keep
The resulting genomes for both haplotype 1 and haplotype 2

## 4. Add "known" inversions to the two haplotypes
Now, all that's left is to repeat this process with the inversion VCFs from [**Step 2**](#2-add-random-inversions):

### haplotype 1
```bash
harpy simulate inversion --conda --vcf VCF_inversions -o sim_inversions_hap1 geno.hap1.fa
```
==- :icon-terminal: code explation
- `--vcf` is the vcf of inversions for haplotype 1 from [**Step 2**](#2-add-random-inversions)
- `-o` is the name of the output directory
- `--conda` is optional and a matter of runtime preference
- `geno.hap1.fa` is the resulting genome from [**Step 3: haplotype 1**](#haplotype-1)
 - i.e. haplotype 1 with the simulated snps, and indels
===

### haplotype 2
```bash
harpy simulate inversion --conda --vcf VCF_inversions -o inversions_hap2 geno.hap2.fa
```
==- :icon-terminal: code explation
- `--vcf` is the vcf of inversions for haplotype 1 from [**Step 2**](#2-add-random-inversions)
- `-o` is the name of the output directory
- `--conda` is optional and a matter of runtime preference
- `geno.hap1.fa` is the resulting genome from [**Step 3: haplotype 2**](#haplotype-2)
 - i.e. haplotype 2 with the simulated snps and indels
===

### 📝 output to keep
The resulting genomes for both haplotype 1 and haplotype 2. These two fasta files make up your diploid genome! 🎉🎉

-------

## 5. Simulating linked-reads
Now that you have heterozygous haplotypes created from your starting genome, you can simulate linked-reads from it using
`harpy simulate linkedreads`. A simple implementation of that could look like:
```bash
harpy simulate linkedreads --conda -t 4 HAP1.fa HAP2.fa
```
==- :icon-terminal: code explation
- `--conda` is optional and a matter of runtime preference
- `-t` is the number of threads to use (`4`)
- `HAP1.fa` is the resulting genome from [**Step 4: haplotype 1**](#haplotype-1-1)
 - i.e. haplotype 1 with the simulated snps, indels, and inversions
- `HAP2.fa` is the resulting genome from [**Step 4: haplotype 2**](#haplotype-2-1)
 - i.e. haplotype 2 with the simulated snps, indels, and inversions
===

"""Module with python-click types for command-line level validations of program-specific extra-params inputs"""

import os
import click
from shlex import split as shellsplit
from shlex import join as shelljoin

def sanitize_shell(sh):
    return shelljoin(shellsplit(sh))

class FastpParams(click.ParamType):
    """A class for a click type that validates fastp extra-params."""
    name = "fastp_params"
    def convert(self, value, param, ctx):
        harpy_options = "--length_required -l --trim_poly_g --max_len1 -b --detect_adapter_for_pe --disable_adapter_trimming -A --dedup -D".split() 
        valid_options = "--failed_out -6 --phred64 --reads_to_process --fix_mgi_id -a --adapter_sequence --adapter_sequence_r2 --adapter_fasta -f --trim_front1 -t --trim_tail1 -b --max_len1 -F --trim_front2 -T --trim_tail2 -B --max_len2 --dup_calc_accuracy --dont_eval_duplication --poly_g_min_len -x --trim_poly_x --poly_x_min_len -5 --cut_front -3 --cut_tail -W --cut_window_size -M --cut_mean_quality --cut_front_window_size --cut_front_mean_quality --cut_tail_window_size --cut_tail_mean_quality --cut_right_window_size --cut_right_mean_quality -Q --disable_quality_filtering -q --qualified_quality_phred -u --unqualified_percent_limit -n --n_base_limit -e --average_qual -L --disable_length_filtering -l --length_required --length_limit -y --low_complexity_filter -Y --complexity_threshold --filter_by_index1 --filter_by_index2 --filter_by_index_threshold -c --correction --overlap_len_require --overlap_diff_limit --overlap_diff_percent_limit -U --umi --umi_loc --umi_len --umi_prefix --umi_skip -p --overrepresentation_analysis -P --overrepresentation_sampling -s --split -S --split_by_lines -d --split_prefix_digits".split()
        opts = 0
        docs = "https://github.com/OpenGene/fastp?tab=readme-ov-file#all-options"
        for i in shellsplit(value):
            if i.startswith("-"):
                opts += 1
                if i in harpy_options:
                    self.fail(f"{i} is already used by Harpy when calling fastp.", param, ctx)
                if i not in valid_options:
                    self.fail(f"{i} is not a valid fastp option. See the fastp documentation for a list of available options: {docs}.", param, ctx)
        if opts < 1:
            self.fail(f"No valid options recognized. Available fastp options begin with one or two dashes (e.g. --phred64 or -a). See the fastp documentation for a list of available options: {docs}.", param, ctx)
        return sanitize_shell(value)

class BwaParams(click.ParamType):
    """A class for a click type that validates bwa extra-params."""
    name = "bwa_params"
    def convert(self, value, param, ctx):
        harpy_options = "-C -v -t -R".split() 
        valid_options = "-k -w -d -r -c -P -A -B -O -E -L -U -p -T -a -H -M ".split()
        opts = 0
        docs = "https://bio-bwa.sourceforge.net/bwa.shtml"
        for i in shellsplit(value):
            if i.startswith("-"):
                opts += 1
                if i in harpy_options:
                    self.fail(f"{i} is already used by Harpy when calling bwa mem.", param, ctx)
                if i not in valid_options:
                    self.fail(f"{i} is not a valid bwa mem option. See the bwa documentation for a list of available options: {docs}.", param, ctx)
        if opts < 1:
            self.fail(f"No valid options recognized. Available bwa options begin with one dash (e.g. -M). See the bwa documentation for a list of available options: {docs}.", param, ctx)
        return sanitize_shell(value)

class EmaParams(click.ParamType):
    """A class for a click type that validates ema extra-params."""
    name = "ema_params"
    def convert(self, value, param, ctx):
        harpy_options = "-t -p -d -r -R -x".split() 
        valid_options = "-i".split()
        opts = 0
        docs = "https://github.com/arshajii/ema"
        for i in shellsplit(value):
            if i.startswith("-"):
                opts += 1
                if i in harpy_options:
                    self.fail(f"{i} is already used by Harpy when calling ema.", param, ctx)
                if i not in valid_options:
                    self.fail(f"{i} is not a valid ema option. See the ema documentation for a list of available options: {docs}.", param, ctx)
        if opts < 1:
            self.fail(f"No valid options recognized. Available ema options begin with one dash (e.g. -i). See the ema documentation for a list of available options: {docs}.", param, ctx)
        return sanitize_shell(value)

class StrobeAlignParams(click.ParamType):
    """A class for a click type that validates strobealign extra-params."""
    name = "strobealign_params"
    def convert(self, value, param, ctx):
        harpy_options = "--use-index -i --create-index -r -N -t -U -C --rg-id --rg".split() 
        valid_options = "-x -A -B -O -E -L -f -S -M -R -m -k -l -u -c -s -b --aux-len --aemb --eqx --no-PG --details".split()
        opts = 0
        docs = "https://github.com/ksahlin/strobealign?tab=readme-ov-file#command-line-options"
        for i in shellsplit(value):
            if i.startswith("-"):
                opts += 1
                if i in harpy_options:
                    self.fail(f"{i} is already used by Harpy when calling strobealign.", param, ctx)
                if i not in valid_options:
                    self.fail(f"{i} is not a valid strobealign option. See the strobealign documentation for a list of available options: {docs}.", param, ctx)
        if opts < 1:
            self.fail(f"No valid options recognized. Available strobealign options begin with one or two dashes (e.g. --eqx or -L). See the strobealign documentation for a list of available options: {docs}.", param, ctx)
        return sanitize_shell(value)

class SpadesParams(click.ParamType):
    """A class for a click type that validates spades extra-params."""
    name = "spades_params"
    def convert(self, value, param, ctx):
        harpy_options = "-t -m -k --gemcode1-1 --gemcode1-2 -o --isolate --pe1-1 --pe1-2".split() 
        valid_options = "--dataset --pacbio --nanopore --sanger --trusted-contigs --untrusted-contigs --assembly-graph --cov-cutoff --phred-offset --custom-hmms --gfa11".split()
        valid_options += [f"--mp{x}-{orient}" for x in range(1,10) for orient in ["1","2","12", "fr", "rf", "ff"]]
        valid_options += [f"--hqmp{x}-{orient}" for x in range(1,10) for orient in ["1","2","12", "s", "fr", "rf", "ff"]]
        opts = 0
        docs = "http://ablab.github.io/spades/running.html"
        for i in shellsplit(value):
            if i.startswith("-"):
                opts += 1
                if i in harpy_options:
                    self.fail(f"{i} is already used by Harpy when calling spades.", param, ctx)
                if i not in valid_options:
                    self.fail(f"{i} is not a valid spades option. See the spades documentation for a list of available options: {docs}.", param, ctx)
        if opts < 1:
            self.fail(f"No valid options recognized. Available spades options begin with two dashes (e.g. --cov-cutoff). See the spades documentation for a list of available options: {docs}.", param, ctx)
        return sanitize_shell(value)

class ArcsParams(click.ParamType):
    """A class for a click type that validates ARCS extra-params."""
    name = "arcs_params"
    def convert(self, value, param, ctx):
        harpy_options = "draft reads t mapq nm dist minsize span c z s l base_name".split() 
        valid_options = "G cut longmap window as trim ac u multfile g graph gap tsv barcodecounts m index_multiplicity d max_degree e end_length r error_percent k k_value j j_index B bin_sizeN D dist_est no_dist_est dist_median dist_upper dist_tsv samples_tsv P pair f d k o e a b r p x".split()
        opts = 0
        docs = "\nTigmint: https://github.com/bcgsc/tigmint\nARCS: https://github.com/bcgsc/arcs\nLINKS: https://github.com/bcgsc/links"
        for i in shellsplit(value):
            if i.startswith("-"):
                self.fail(f"{i} begins with a dash, which would be interpreted as an argument to arcs-make rather than arcs. To avoid unexpected errors, arguments to arcs-make are disallowed. If this was inteded to be an argument to arcs, try using " + i.lstrip("-") + "=VAL instead", param, ctx)
            if "=" in i:
                opts += 1
                argsplit = i.split("=")
                if len(argsplit) != 2:
                     self.fail(f"{i} is not a valid arcs option. Valid arcs options begin without dashes and must be in the form ARG=VAL, without spaces (e.g. k=15). See the documentation for a list of available options.{docs}", param, ctx)
                arg = argsplit[0].strip()
                if arg in harpy_options:
                    self.fail(f"{arg} is already used by Harpy when calling arcs.", param, ctx)
                if arg not in valid_options:
                    self.fail(f"{arg} is not a valid arcs option. See the documentation for a list of available options.{docs}", param, ctx)
        if opts < 1:
            self.fail(f"No valid options recognized. Valid arcs options begin without dashes and must be in the form ARG=VAL, without spaces (e.g. k=15). See the documentation for a list of available options.{docs}", param, ctx)
        return sanitize_shell(value)

class StitchParams(click.ParamType):
    """A class for a click type that validates stitch extra-params. Sanitizes and corrects different input styles to work with STITCH cli"""
    name = "stitch_params"
    def convert(self, value, param, ctx):
        harpy_options = "--method --posfile --bamlist --nCores --nGen --chr --buffer --regionStart --regionEnd --K --S --use_bx_tag --bxTagUpperLimit --outputdir --output_filename --tempdir".split() 
        valid_options = "--nStarts --genfile --B_bit_prob --outputInputInVCFFormat --downsampleToCov --downsampleFraction --readAware --chrStart --chrEnd --maxDifferenceBetweenReads --maxEmissionMatrixDifference --alphaMatThreshold --emissionThreshold --iSizeUpperLimit --bqFilter --niterations --shuffleHaplotypeIterations --splitReadIterations --expRate --maxRate --minRate --Jmax --regenerateInput --originalRegionName --keepInterimFiles --keepTempDir --switchModelIteration --generateInputOnly --restartIterations --refillIterations --downsampleSamples --downsampleSamplesKeepList --subsetSNPsfile --useSoftClippedBases --outputBlockSize --outputSNPBlockSize --inputBundleBlockSize --genetic_map_file --reference_haplotype_file --reference_legend_file --reference_sample_file --reference_populations --reference_phred --reference_iterations --reference_shuffleHaplotypeIterations --initial_min_hapProb --initial_max_hapProb --regenerateInputWithDefaultValues --plot_shuffle_haplotype_attempts --save_sampleReadsInfo --gridWindowSize --shuffle_bin_nSNPs --shuffle_bin_radius --keepSampleReadsInRAM --useTempdirWhileWriting --output_haplotype_dosages".split()
        opts = 0
        docs = "https://github.com/rwdavies/STITCH/blob/master/Options.md"
        clean_args = []
        for i in shellsplit(value):
            if not i.startswith("-"):
                i = "--" + i.lstrip("-")
            if "=" in i:
                opts += 1
                argsplit = [j.strip() for j in i.split("=")]
                if len(argsplit) != 2:
                    self.fail(f"{i} is not in the proper format for STITCH. STITCH options must be in the form ARG=VAL (e.g. --downsampleFraction=0.5). See the stitch documentation for a list of available options: {docs}", param, ctx)
                arg = argsplit[0]
                if arg in harpy_options:
                    self.fail(f"{arg} is already used by Harpy when calling STITCH.", param, ctx)
                if arg not in valid_options:
                    self.fail(f"{arg} is not a valid STITCH option. See the STITCH documentation for a list of available options: {docs}", param, ctx)
                clean_args.append("=".join(argsplit))
            else:
                self.fail(f"{i} is not in the proper format for STITCH. STITCH options must be in the form ARG=VAL (e.g. --downsampleFraction=0.5). See the stitch documentation for a list of available options: {docs}", param, ctx)
        if opts < 1:
            self.fail(f"No valid options recognized. STITCH options begin with a double-dash and must be in the form --ARG=VAL (e.g. --downsampleFraction=0.5). See the stitch documentation for a list of available options: {docs}.", param, ctx)
        return sanitize_shell(" ".join(clean_args))

class HapCutParams(click.ParamType):
    """A class for a click type that validates hapcut2 extra-params."""
    name = "hapcut2_params"
    def convert(self, value, param, ctx):
        harpy_options = "--fragments --vcf --out --nf --error_analysis_mode --call_homozygous --outvcf --threshold --no_prune".split() 
        valid_options = "--skip_prune --sp --discrete_pruning --dp --max_iter --mi --maxcut_iter --mc".split()
        opts = 0
        docs = "https://github.com/vibansal/HapCUT2"
        for i in shellsplit(value):
            if i.startswith("-"):
                opts += 1
                if i in harpy_options:
                    self.fail(f"{i} is already used by Harpy when calling hapcut2.", param, ctx)
                if i not in valid_options:
                    self.fail(f"{i} is not a valid hapcut2 option. See the hapcut2 documentation for a list of available options: {docs}.", param, ctx)
        if opts < 1:
            self.fail(f"No valid options recognized. Available hapcut2 options begin with two dashes (e.g. --dp). See the hapcut2 documentation for a list of available options: {docs}.", param, ctx)
        return sanitize_shell(value)

class LeviathanParams(click.ParamType):
    """A class for a click type that validates leviathan extra-params."""
    name = "leviathan_params"
    def convert(self, value, param, ctx):
        harpy_options = "-b -i -g -o -v --minVariantSize -c --minBarcodes -B --nbBins -t --threads -C --candidates -s --smallRate -m --mediumRate -l --largeRate -d --duplicates".split() 
        valid_options = "-r --regionSize -n --maxLinks -M --mediumSize -L --largeSize -s --skipTranslocations -p --poolSize".split()
        opts = 0
        docs = "https://github.com/morispi/LEVIATHAN?tab=readme-ov-file#options"
        for i in shellsplit(value):
            if i.startswith("-"):
                opts += 1
                if i in harpy_options:
                    self.fail(f"{i} is already used by Harpy when calling leviathan.", param, ctx)
                if i not in valid_options:
                    self.fail(f"{i} is not a valid leviathan option. See the leviathan documentation for a list of available options: {docs}.", param, ctx)
        if opts < 1:
            self.fail(f"No valid options recognized. Available leviathan options begin with one or two dashes (e.g. -m or -mediumRate). See the leviathan documentation for a list of available options: {docs}.", param, ctx)
        return sanitize_shell(value)

class NaibrParams(click.ParamType):
    """A class for a click type that validates naibr extra-params."""
    name = "naibr_params"
    def convert(self, value, param, ctx):
        harpy_options = "bam_file prefix outdir threads min_mapq d min_sv k".split() 
        valid_options = "blacklist candidates min_discs min_reads sd_mult".split()
        opts = 0
        docs = "https://github.com/pontushojer/NAIBR?tab=readme-ov-file#running-naibr"
        clean_args = []
        for idx,i in enumerate(shellsplit(value.replace("-", ""))):
            # if it's an even index, it's the argument name of an arg-val pair
            if idx % 2 == 0:
                opts += 1
                if i in harpy_options:
                    self.fail(f"{i} is already used by Harpy when calling naibr.", param, ctx)
                if i not in valid_options:
                    self.fail(f"{i} is not a valid naibr option. See the naibr documentation for a list of available options: {docs}.", param, ctx)
            clean_args.append(i.strip())
        if opts < 1:
            self.fail(f"No valid options recognized. Available naibr options begin without dashes in the form of ARG<space>VAL (e.g. blacklist inversions.txt). See the naibr documentation for a list of available options: {docs}.", param, ctx)
        return sanitize_shell(" ".join(clean_args))
    
class MpileupParams(click.ParamType):
    """A class for a click type that validates mpileup extra-params."""
    name = "mpileup_params"
    def convert(self, value, param, ctx):
        harpy_options = "--fasta-ref -f --bam-list -b --annotate -a --output-type -O -r --regions".split() 
        valid_options = "-6 --illumina1.3+ -A --count-orphans -B --no-BAQ -C --adjust-MQ -D --full-BAQ -d --max-depth -E --redo-BAQ -G --read-groups -q --min-MQ -Q --min-BQ --max-BQ INT --delta-BQ INT --ignore-RG --ls --skip-all-set --ns --skip-any-set --lu --skip-all-unset --nu --skip-any-unset -s --samples -S --samples-file -t --targets -T --targets-file -x --ignore-overlaps -g --gvcf -o -X --config -e --ext-prob -F --gap-frac -h --tandem-qual -I --skip-indels -L --max-idepth -m --min-ireads -M --max-read-len -o --open-prob -p --per-sample-mF -P --platforms --ar --ambig-reads --indel-bias --del-bias --score-vs-ref --indel-size --indels-2.0 --indels-cns --seqq-offset --no-indels-cns --poly-mqual".split()
        opts = 0
        docs = "https://samtools.github.io/bcftools/bcftools.html#mpileup"
        for i in shellsplit(value):
            if i.startswith("-"):
                opts += 1
                if i in harpy_options:
                    self.fail(f"{i} is already used by Harpy when calling mpileup.", param, ctx)
                if i not in valid_options:
                    self.fail(f"{i} is not a valid mpileup option. See the mpileup documentation for a list of available options: {docs}.", param, ctx)
        if opts < 1:
            self.fail(f"No valid options recognized. Available mpileup options begin one or with two dashes (e.g. -d or --max-depth). See the mpileup documentation for a list of available options: {docs}.", param, ctx)
        return sanitize_shell(value)
    
class FreebayesParams(click.ParamType):
    """A class for a click type that validates freebayes extra-params."""
    name = "freebayes_params"
    def convert(self, value, param, ctx):
        harpy_options = "-r --region -p --ploidy -f --fasta-reference -L --bam-list --populations".split() 
        valid_options = "-t --targets -s --samples -A --cnv-map -v --vcf --gvcf --gvcf-chunkUM -& --gvcf-dont-use-chunk -@ --variant-input -l --only-use-input-alleles --haplotype-basis-alleles --report-all-haplotype-alleles --report-monomorphic -P --pvar --strict-vcf -T --theta -J --pooled-discrete -K --pooled-continuous -Z --use-reference-allele --reference-quality -n --use-best-n-alleles -E --max-complex-gap --haplotype-length --min-repeat-size --min-repeat-entropy --no-partial-observations -O --dont-left-align-indels -4 --use-duplicate-reads -m --min-mapping-quality -q --min-base-quality -R --min-supporting-allele-qsum -Y --min-supporting-mapping-qsum -Q --mismatch-base-quality-threshold -U --read-mismatch-limit -z --read-max-mismatch-fraction -$ --read-snp-limit -e --read-indel-limit -0 --standard-filters -F --min-alternate-fraction -C --min-alternate-count -3 --min-alternate-qsum -G --min-alternate-total --min-coverage --limit-coverage -g --skip-coverage --trim-complex-tail -k --no-population-priors -w --hwe-priors-off -V --binomial-obs-priors-off -a --allele-balance-priors-off --observation-bias --base-quality-cap --prob-contamination --legacy-gls --contamination-estimates --report-genotype-likelihood-max -B --genotyping-max-iterations --genotyping-max-banddepth -W --posterior-integration-limits,M -N --exclude-unobserved-genotypes -S --genotype-variant-threshold -j --use-mapping-quality -H --harmonic-indel-quality -D --read-dependence-factor -= --genotype-qualities -d --debug".split()
        opts = 0
        docs = "https://github.com/freebayes/freebayes"
        for i in shellsplit(value):
            if i.startswith("-"):
                opts += 1
                if i in harpy_options:
                    self.fail(f"{i} is already used by Harpy when calling freebayes.", param, ctx)
                if i not in valid_options:
                    self.fail(f"{i} is not a valid freebayes option. See the freebayes documentation for a list of available options: {docs}.", param, ctx)
        if opts < 1:
            self.fail(f"No valid options recognized. Available freebayes options begin with one or two dashes (e.g. -t or --targets). See the freebayes documentation for a list of available options: {docs}.", param, ctx)
        return sanitize_shell(value)

class Barcodes(click.ParamType):
    """A class for a click type which accepts either a file or two integers, separated by a comma."""
    name = "barcodes"
    def convert(self, value, param, ctx):
        if os.path.isfile(value):
            return os.path.abspath(value)
        try:
            bp,count = value.split(",")
        except ValueError:
            self.fail(f"{value} is not a file, nor in int,int format", param, ctx)
        try:
            bp = int(bp)
        except ValueError:
            self.fail(f"{value} is not an integer.", param, ctx)
        try:
            count = int(count)
        except ValueError:
            self.fail(f"{value} is not an integer.", param, ctx)
        return f"{bp},{count}"

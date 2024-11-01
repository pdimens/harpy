"""Module with python-click types for command-line level validations of program-specific extra-params inputs"""

import click

class FastpParams(click.ParamType):
    """A class for a click type that validates fastp extra-params."""
    name = "fastp_params"
    def convert(self, value, param, ctx):
        harpy_options = "--length_required -l --max_len1 -b --detect_adapter_for_pe --disable_adapter_trimming -A --dedup -D".split() 
        valid_options = "--failed_out -6 --phred64 --reads_to_process --fix_mgi_id -a --adapter_sequence --adapter_sequence_r2 --adapter_fasta -f --trim_front1 -t --trim_tail1 -b --max_len1 -F --trim_front2 -T --trim_tail2 -B --max_len2 --dup_calc_accuracy --dont_eval_duplication --poly_g_min_len -x --trim_poly_x --poly_x_min_len -5 --cut_front -3 --cut_tail -W --cut_window_size -M --cut_mean_quality --cut_front_window_size --cut_front_mean_quality --cut_tail_window_size --cut_tail_mean_quality --cut_right_window_size --cut_right_mean_quality -Q --disable_quality_filtering -q --qualified_quality_phred -u --unqualified_percent_limit -n --n_base_limit -e --average_qual -L --disable_length_filtering -l --length_required --length_limit -y --low_complexity_filter -Y --complexity_threshold --filter_by_index1 --filter_by_index2 --filter_by_index_threshold -c --correction --overlap_len_require --overlap_diff_limit --overlap_diff_percent_limit -U --umi --umi_loc --umi_len --umi_prefix --umi_skip -p --overrepresentation_analysis -P --overrepresentation_sampling -s --split -S --split_by_lines -d --split_prefix_digits".split()
        opts = 0
        for i in value.split():
            if i.startswith("-"):
                opts += 1
                if i in harpy_options:
                    self.fail(f"{i} is already used by Harpy when calling fastp.", param, ctx)
                if i not in valid_options:
                    self.fail(f"{i} is not a valid fastp option. See the fastp documentation for a list of available options: https://github.com/OpenGene/fastp.", param, ctx)
        if opts < 1:
            self.fail("No valid options recognized. Available fastp options begin with one or two dashes (e.g. --phred64 or -a). See the fastp documentation for a list of available options: https://github.com/OpenGene/fastp.", param, ctx)
        return value

class BwaParams(click.ParamType):
    """A class for a click type that validates bwa extra-params."""
    name = "bwa_params"
    def convert(self, value, param, ctx):
        harpy_options = "-C -v -t -R".split() 
        valid_options = "-k -w -d -r -c -P -A -B -O -E -L -U -p -T -a -H -M ".split()
        opts = 0
        for i in value.split():
            if i.startswith("-"):
                opts += 1
                if i in harpy_options:
                    self.fail(f"{i} is already used by Harpy when calling bwa mem.", param, ctx)
                if i not in valid_options:
                    self.fail(f"{i} is not a valid bwa mem option. See the bwa documentation for a list of available options: https://bio-bwa.sourceforge.net/bwa.shtml.", param, ctx)
        if opts < 1:
            self.fail("No valid options recognized. Available bwa options begin with one dash (e.g. -M). See the bwa documentation for a list of available options: https://bio-bwa.sourceforge.net/bwa.shtml.", param, ctx)
        return value

class EmaParams(click.ParamType):
    """A class for a click type that validates ema extra-params."""
    name = "ema_params"
    def convert(self, value, param, ctx):
        harpy_options = "-t -p -d -r -R -x".split() 
        valid_options = "-i".split()
        opts = 0
        for i in value.split():
            if i.startswith("-"):
                opts += 1
                if i in harpy_options:
                    self.fail(f"{i} is already used by Harpy when calling ema.", param, ctx)
                if i not in valid_options:
                    self.fail(f"{i} is not a valid ema option. See the fastp documentation for a list of available options: https://github.com/arshajii/ema.", param, ctx)
        if opts < 1:
            self.fail("No valid options recognized. Available ema options begin with one dash (e.g. -i). See the ema documentation for a list of available options: https://github.com/arshajii/ema.", param, ctx)
        return value

class StrobeAlignParams(click.ParamType):
    """A class for a click type that validates strobealign extra-params."""
    name = "strobealign_params"
    def convert(self, value, param, ctx):
        harpy_options = "--use-index -i --create-index -r -N -t -U -C --rg-id --rg".split() 
        valid_options = "-x -A -B -O -E -L -f -S -M -R -m -k -l -u -c -s -b --aux-len --aemb --eqx --no-PG --details".split()
        opts = 0
        for i in value.split():
            if i.startswith("-"):
                opts += 1
                if i in harpy_options:
                    self.fail(f"{i} is already used by Harpy when calling strobealign.", param, ctx)
                if i not in valid_options:
                    self.fail(f"{i} is not a valid strobealign option. See the strobealign documentation for a list of available options: https://github.com/ksahlin/strobealign.", param, ctx)
        if opts < 1:
            self.fail("No valid options recognized. Available strobealign options begin with one or two dashes (e.g. --eqx or -L). See the strobealign documentation for a list of available options: https://github.com/ksahlin/strobealign.", param, ctx)
        return value

class SpadesParams(click.ParamType):
    """A class for a click type that validates spades extra-params."""
    name = "spades_params"
    def convert(self, value, param, ctx):
        harpy_options = "-t -m -k --gemcode1-1 --gemcode1-2 -o --isolate --pe1-1 --pe1-2".split() 
        valid_options = "--dataset --pacbio --nanopore --sanger --trusted-contigs --untrusted-contigs --assembly-graph --cov-cutoff --phred-offset --custom-hmms --gfa11".split()
        valid_options += [f"--mp{x}-{orient}" for x in range(1,10) for orient in ["1","2","12", "fr", "rf", "ff"]]
        valid_options += [f"--hqmp{x}-{orient}" for x in range(1,10) for orient in ["1","2","12", "s", "fr", "rf", "ff"]]
        opts = 0
        for i in value.split():
            if i.startswith("-"):
                opts += 1
                if i in harpy_options:
                    self.fail(f"{i} is already used by Harpy when calling spades.", param, ctx)
                if i not in valid_options:
                    self.fail(f"{i} is not a valid spades option. See the spades documentation for a list of available options: http://ablab.github.io/spades/running.html.", param, ctx)
        if opts < 1:
            self.fail("No valid options recognized. Available spades options begin with two dashes (e.g. --cov-cutoff). See the spades documentation for a list of available options: http://ablab.github.io/spades/running.html.", param, ctx)
        return value

class ArcsParams(click.ParamType):
    """A class for a click type that validates ARCS extra-params."""
    name = "arcs_params"
    def convert(self, value, param, ctx):
        harpy_options = "draft reads t mapq nm dist minsize span c z s l base_name".split() 
        valid_options = "G cut longmap window as trim ac u multfile g graph gap tsv barcodecounts m index_multiplicity d max_degree e end_length r error_percent k k_value j j_index B bin_sizeN D dist_est no_dist_est dist_median dist_upper dist_tsv samples_tsv P pair f d k o e a b r p x".split()
        opts = 0
        for i in value.split():
            if i.startswith("-"):
                self.fail(f"{i} begins with a dash, which would be interpreted as an argument to arcs-make rather than arcs. To avoid unexpected errors, arguments to arcs-make are disallowed. If this was inteded to be an argument to arcs, try using " + i.lstrip("-") + "=VAL instead", param, ctx)
            if "=" in i:
                opts += 1
                arg = i.split("=")[0].strip()
                if arg in harpy_options:
                    self.fail(f"{arg} is already used by Harpy when calling arcs.", param, ctx)
                if arg not in valid_options:
                    self.fail(f"{arg} is not a valid arcs option. See the documentation for a list of available options.\nTigmint: https://github.com/bcgsc/tigmint\nARCS: https://github.com/bcgsc/arcs\nLINKS: https://github.com/bcgsc/links", param, ctx)
        if opts < 1:
            self.fail("No valid options recognized. Available arcs options begin without dashes and must be in the form ARG=VAL, without spaces (e.g. k=15). See the documentation for a list of available options.\nTigmint: https://github.com/bcgsc/tigmint\nARCS: https://github.com/bcgsc/arcs\nLINKS: https://github.com/bcgsc/links", param, ctx)
        return value

class StitchParams(click.ParamType):
    """A class for a click type that validates stitch extra-params."""
    name = "stitch_params"
    def convert(self, value, param, ctx):
        harpy_options = "method posfile bamlist nCores nGen chr K S use_bx_tag bxTagUpperLimit outputdir output_filename tempdir".split() 
        valid_options = "nStarts sampleNames_file genfile B_bit_prob outputInputInVCFFormat downsampleToCov downsampleFraction readAware chrStart chrEnd regionStart regionEnd buffer maxDifferenceBetweenReads maxEmissionMatrixDifference alphaMatThreshold emissionThreshold iSizeUpperLimit bqFilter niterations shuffleHaplotypeIterations splitReadIterations expRate maxRate minRate Jmax regenerateInput originalRegionName keepInterimFiles keepTempDir outputHaplotypeProbabilities switchModelIteration generateInputOnly restartIterations refillIterations downsampleSamples downsampleSamplesKeepList subsetSNPsfile useSoftClippedBases outputBlockSize outputSNPBlockSize inputBundleBlockSize genetic_map_file reference_haplotype_file reference_legend_file reference_sample_file reference_populations reference_phred reference_iterations reference_shuffleHaplotypeIterations initial_min_hapProb initial_max_hapProb regenerateInputWithDefaultValues plotHapSumDuringIterations plot_shuffle_haplotype_attempts plotAfterImputation save_sampleReadsInfo gridWindowSize shuffle_bin_nSNPs shuffle_bin_radius keepSampleReadsInRAM useTempdirWhileWriting output_haplotype_dosages".split()
        opts = 0
        for i in value.split():
            if i.startswith("-"):
                self.fail(f"{i} begins with a dash, which is the wrong format. Try using " + i.lstrip("-") + "=VAL instead", param, ctx)
            if "=" in i:
                opts += 1
                arg = i.split("=")[0].strip()
                if arg in harpy_options:
                    self.fail(f"{arg} is already used by Harpy when calling stitch.", param, ctx)
                if arg not in valid_options:
                    self.fail(f"{arg} is not a valid stitch option. See the stitch documentation for a list of available options: https://github.com/rwdavies/STITCH/blob/master/Options.md", param, ctx)
        if opts < 1:
            self.fail("No valid options recognized. Available stitch options begin without dashes and must be in the form ARG=VAL (e.g. downsampleFraction=0.5). See the stitch documentation for a list of available options: https://github.com/rwdavies/STITCH/blob/master/Options.md.", param, ctx)
        return value

class HapCutParams(click.ParamType):
    """A class for a click type that validates hapcut2 extra-params."""
    name = "hapcut2_params"
    def convert(self, value, param, ctx):
        harpy_options = "--fragments --vcf --out --nf --error_analysis_mode --call_homozygous --outvcf --threshold --no_prune".split() 
        valid_options = "--skip_prune --sp --discrete_pruning --dp --max_iter --mi --maxcut_iter --mc".split()
        opts = 0
        for i in value.split():
            if i.startswith("-"):
                opts += 1
                if i in harpy_options:
                    self.fail(f"{i} is already used by Harpy when calling hapcut2.", param, ctx)
                if i not in valid_options:
                    self.fail(f"{i} is not a valid hapcut2 option. See the hapcut2 documentation for a list of available options: https://github.com/vibansal/HapCUT2.", param, ctx)
        if opts < 1:
            self.fail("No valid options recognized. Available hapcut2 options begin with two dashes (e.g. --dp). See the hapcut2 documentation for a list of available options: https://github.com/vibansal/HapCUT2.", param, ctx)
        return value

class LeviathanParams(click.ParamType):
    """A class for a click type that validates leviathan extra-params."""
    name = "leviathan_params"
    def convert(self, value, param, ctx):
        harpy_options = "-b -i -g -o -v --minVariantSize -c --minBarcodes -B --nbBins -t --threads -C --candidates".split() 
        valid_options = "-r --regionSize -n --maxLinks -M --mediumSize -L --largeSize -s --smallRate -m --mediumRate -l --largeRate -d --duplicates -s --skipTranslocations -p --poolSize".split()
        opts = 0
        for i in value.split():
            if i.startswith("-"):
                opts += 1
                if i in harpy_options:
                    self.fail(f"{i} is already used by Harpy when calling leviathan.", param, ctx)
                if i not in valid_options:
                    self.fail(f"{i} is not a valid leviathan option. See the leviathan documentation for a list of available options: https://github.com/morispi/LEVIATHAN.", param, ctx)
        if opts < 1:
            self.fail("No valid options recognized. Available leviathan options begin with one or two dashes (e.g. --mediumRate or -m). See the leviathan documentation for a list of available options: https://github.com/morispi/LEVIATHAN.", param, ctx)
        return value

class NaibrParams(click.ParamType):
    """A class for a click type that validates naibr extra-params."""
    name = "naibr_params"
    def convert(self, value, param, ctx):
        harpy_options = "bam_file prefix outdir threads min_mapq d min_sv k".split() 
        valid_options = "blacklist candidates".split()
        opts = 0
        for idx,i in enumerate(value.split()):
            if i.startswith("-"):
                self.fail(f"{i} begins with a dash, which is the wrong format. Try using " + i.lstrip("-") + " VAL instead", param, ctx)
            # if it's an odd index, it's the first (arg) of an arg val pair
            if idx % 2 != 0:
                opts += 1
                if i in harpy_options:
                    self.fail(f"{i} is already used by Harpy when calling naibr.", param, ctx)
                if i not in valid_options:
                    self.fail(f"{i} is not a valid naibr option. See the naibr documentation for a list of available options: https://github.com/pontushojer/NAIBR.", param, ctx)
        if opts < 1:
            self.fail("No valid options recognized. Available naibr options begin without dashes in the form of ARG VAL (e.g. blacklist inversions.ignore). See the naibr documentation for a list of available options: https://github.com/pontushojer/NAIBR.", param, ctx)
        return value
    
class MpileupParams(click.ParamType):
    """A class for a click type that validates mpileup extra-params."""
    name = "mpileup_params"
    def convert(self, value, param, ctx):
        harpy_options = "".split() 
        valid_options = "".split()
        opts = 0
        for i in value.split():
            if i.startswith("-"):
                opts += 1
                if i in harpy_options:
                    self.fail(f"{i} is already used by Harpy when calling mpileup.", param, ctx)
                if i not in valid_options:
                    self.fail(f"{i} is not a valid mpileup option. See the mpileup documentation for a list of available options: XXXX.", param, ctx)
        if opts < 1:
            self.fail("No valid options recognized. Available mpileup options begin with two dashes (e.g. --eqx or -L). See the mpileup documentation for a list of available options: XXXX.", param, ctx)
        return value
    
class FreebayesParams(click.ParamType):
    """A class for a click type that validates freebayes extra-params."""
    name = "freebayes_params"
    def convert(self, value, param, ctx):
        harpy_options = "".split() 
        valid_options = "".split()
        opts = 0
        for i in value.split():
            if i.startswith("-"):
                opts += 1
                if i in harpy_options:
                    self.fail(f"{i} is already used by Harpy when calling freebayes.", param, ctx)
                if i not in valid_options:
                    self.fail(f"{i} is not a valid freebayes option. See the freebayes documentation for a list of available options: XXXX.", param, ctx)
        if opts < 1:
            self.fail("No valid options recognized. Available freebayes options begin with two dashes (e.g. --eqx or -L). See the freebayes documentation for a list of available options: XXXX.", param, ctx)
        return value





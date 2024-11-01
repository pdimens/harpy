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
            self.fail("No valid options recognized. All arcs options begin with one or two dashes (e.g. --eqx or -L). See the documentation for a list of available options.\nTigmint: https://github.com/bcgsc/tigmint\nARCS: https://github.com/bcgsc/arcs\nLINKS: https://github.com/bcgsc/links", param, ctx)
        return value

class XXXParams(click.ParamType):
    """A class for a click type that validates XXX extra-params."""
    name = "XXX_params"
    def convert(self, value, param, ctx):
        harpy_options = "".split() 
        valid_options = "".split()
        opts = 0
        for i in value.split():
            if i.startswith("-"):
                opts += 1
                if i in harpy_options:
                    self.fail(f"{i} is already used by Harpy when calling XXX.", param, ctx)
                if i not in valid_options:
                    self.fail(f"{i} is not a valid XXX option. See the XXX documentation for a list of available options: XXXX.", param, ctx)
        if opts < 1:
            self.fail("No valid options recognized. All XXX options begin with two dashes (e.g. --eqx or -L). See the XXX documentation for a list of available options: XXXX.", param, ctx)
        return value



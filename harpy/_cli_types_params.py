"""Module with python-click types for command-line level validations of program-specific extra-params inputs"""

import click

class FastpParams(click.ParamType):
    """A class for a click type that validates fastp extra-params."""
    name = "fastp_params"
    def convert(self, value, param, ctx):
        harpy_fastp = "--length_required -l --max_len1 -b --detect_adapter_for_pe --disable_adapter_trimming -A --dedup -D".split() 
        faspt_options = " --failed_out -6 --phred64 --reads_to_process --fix_mgi_id -a --adapter_sequence --adapter_sequence_r2 --adapter_fasta -f --trim_front1 -t --trim_tail1 -b --max_len1 -F --trim_front2 -T --trim_tail2 -B --max_len2 --dup_calc_accuracy --dont_eval_duplication --poly_g_min_len -x --trim_poly_x --poly_x_min_len -5 --cut_front -3 --cut_tail -W --cut_window_size -M --cut_mean_quality --cut_front_window_size --cut_front_mean_quality --cut_tail_window_size --cut_tail_mean_quality --cut_right_window_size --cut_right_mean_quality -Q --disable_quality_filtering -q --qualified_quality_phred -u --unqualified_percent_limit -n --n_base_limit -e --average_qual -L --disable_length_filtering -l --length_required --length_limit -y --low_complexity_filter -Y --complexity_threshold --filter_by_index1 --filter_by_index2 --filter_by_index_threshold -c --correction --overlap_len_require --overlap_diff_limit --overlap_diff_percent_limit -U --umi --umi_loc --umi_len --umi_prefix --umi_skip -p --overrepresentation_analysis -P --overrepresentation_sampling -s --split -S --split_by_lines -d --split_prefix_digits".split()
        opts = 0
        for i in value.split():
            if i.startswith("-"):
                opts += 1
                if i in harpy_fastp:
                    self.fail(f"{i} is already used by Harpy when calling fastp.")
                if i not in faspt_options:
                    self.fail(f"{i} is not a valid fastp option. See the fastp documentation for a list of available fastp options: https://github.com/OpenGene/fastp.")
        if opts < 1:
            self.fail("No valid options recognized. All fastp options begin with one or two dashes (e.g. --phred64 or -a). See the fastp documentation for a list of available fastp options: https://github.com/OpenGene/fastp.")
        return value
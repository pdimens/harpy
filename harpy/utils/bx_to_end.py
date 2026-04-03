import re
import sys
import gzip
import click
import pysam

def format_bam(record):
    tags = record.get_tags()
    if "BX" not in tags:
        return record
    bx = record.get_tag("BX")
    # delete the BX tag
    tags = [i for i in tags if i[0] != "BX"]
    # add BX tag to end
    tags.append(("BX", bx))
    record.set_tags(tags)
    return record

def format_fastq(record):
    if "BX:Z" in record.comment:
        splitcomment = record.comment.split()
        # find the BX:Z tag, remove it, and add it to the end
        BX_idx = next((index for index, value in enumerate(splitcomment) if value.startswith("BX:Z")), None)
        bx_tag = splitcomment.pop(BX_idx)
        splitcomment += [bx_tag]
        comment = "\t".join(splitcomment)
    else:
        comment = record.comment
    fastq_req = [
        f"{record.name}\t" + comment,
        record.sequence,
        "+",
        record.quality
    ]
    return "\n".join(fastq_req) + "\n"

@click.command(no_args_is_help = True, epilog = "Documentation: https://pdimens.github.io/harpy/workflows/preprocess/")
@click.argument('input', required = True, type=click.Path(exists = True, dir_okay=False, resolve_path=True))
@click.help_option('--help', hidden = True)
def bx_to_end(input):
    """
    Move BX:Z tag to the end of records
    
    Input can be FASTQ or SAM/BAM. Writes to stdout
    """
    # VALIDATIONS
    fq_ext = re.compile(r"\.f(?:ast)?q(?:\.gz)?$", flags=re.IGNORECASE)
    if input.lower().endswith(".bam") or input.lower().endswith(".sam"):
        is_fastq = False
    elif fq_ext.search(input):
        is_fastq = True
    else:
        sys.stderr.write(f"Filetype not recognized as one of BAM or FASTQ for file {input}\n")
        sys.exit(1)

    if is_fastq:
        with (
            pysam.FastxFile(input, persist=False) as fq_in,
            gzip.GzipFile(fileobj= sys.stdout.buffer, mode= "wb", compresslevel=6) as fq_out
        ):
            for rec in fq_in:
                fq_out.write(format_fastq(rec).encode())
    else:
        try:
            bam_in = pysam.AlignmentFile(input, "rb", require_index=False)
        except (OSError, ValueError):
            try:
                bam_in = pysam.AlignmentFile(input, "r", require_index=False)
            except (OSError, ValueError) as e:
                print(f"Could not process {input} as a SAM/BAM file. See the error from pysam: {e}")
                sys.exit(1)
        with pysam.AlignmentFile(sys.stdout.buffer, "wb", template = bam_in) as bam_out:
            for aln_rec in bam_in.fetch(until_eof=True):
                bam_out.write(format_bam(aln_rec))
        bam_in.close()

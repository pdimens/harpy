#! /usr/bin/env python

import os
import re
import sys
import argparse
import subprocess
import pysam

def main():
    parser = argparse.ArgumentParser(
        prog='separate_singletons',
        description='Isolate singleton and non-singleton linked-read BAM records into separate files.',
        usage = "separate_singletons -t threads -b barcode_tag -s singletons.bam input.bam > output.bam",
        )
    parser.add_argument("-b", dest = "bx_tag", metavar = "barcode_tag", type=str, default = "BX", help="The header tag with the barcode (default: %(default)s)")
    parser.add_argument("-s", dest = "singletons", metavar = "singletons_file", type=str, default = "singletons.bam", help="Name of output singleton file (default: %(default)s)")
    parser.add_argument("-t", dest = "threads", metavar="threads", type=int, default = 4, help="Number of threads to use (default: %(default)s)")
    parser.add_argument('input', type = str, help = "Input bam file")
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()
    if args.threads <1:
        parser.error(f"Threads supplied to -t ({args.threads}) must be positive (e.g. >1)")
    if not os.path.exists(args.input):
        parser.error(f"{args.input} was not found")
    if len(args.bx_tag) != 2 or not args.bx_tag.isalnum():
        parser.error(f"The header tag supplied to -b ({args.bx_tag}) must be alphanumeric and exactly two characters long")

    invalid_pattern = re.compile(r'[AaBbCcDd]00')
    sorted_bam = f"{args.input[:-4]}.bxsort.bam"
    subprocess.run(f"samtools sort -@ {args.threads} -o {sorted_bam} -t {args.bx_tag} {args.input}".split(), stderr=sys.stderr)
    with (
        pysam.AlignmentFile(sorted_bam, "rb", check_sq=False) as infile,
        pysam.AlignmentFile(sys.stdout, "wb", template=infile) as nonsingleton,
        pysam.AlignmentFile(args.singletons, "wb", template=infile) as singleton,                      
    ):
        record_store = []
        read_count = 0
        last_barcode = None
        for record in infile:
            try:
                barcode = record.get_tag(args.bx_tag)
                if isinstance(barcode, int):
                    pass # an int from an MI-type tag
                elif invalid_pattern.search(barcode):
                    continue
            except KeyError:
                continue
            # write the stored records when the barcode changes
            if last_barcode and barcode != last_barcode:
                target_file = nonsingleton if read_count > 1 else singleton  
                for record in record_store:  
                    target_file.write(record)  

                # reset the record store and read count
                record_store = []
                read_count = 0

            record_store.append(record)
            if record.is_forward:
                # +1 for a forward read, whether it is paired or not
                read_count += 1
            elif record.is_reverse and not record.is_paired:
                # +1 for reverse only if it's unpaired, so the paired read doesn't count twice
                read_count += 1
            # update the last barcode with the current one
            last_barcode = barcode
        # After the for loop ends
        if record_store:
            target_file = nonsingleton if read_count > 1 else singleton
            for i in record_store:
                target_file.write(i)

    # final housekeeping to remove intermediate
    os.remove(sorted_bam)
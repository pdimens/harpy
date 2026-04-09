
import os
from pathlib import Path
import pysam
import pysam.bcftools
from shutil import which
import subprocess
from harpy.common.printing import HarpyPrint

class VCF():
    '''
    A class to contain and validate a VCF input file.
    '''
    def __init__(self, filename:str, workdir:str, quiet:int = 0):
        os.makedirs(workdir, exist_ok = True)

        self.file: str = filename
        self.workdir: str = workdir
        self.biallelic_file: str = ""
        self.contigs: dict[str,int] = {}
        self.print = HarpyPrint(quiet)

    def get_contigs(self):
        """reads the header of a vcf/bcf file and populate `self.contigs` with the contigs (keys) and their lengths (values)"""
        with pysam.VariantFile(self.file) as _vcf:
            for _contig,_info in _vcf.header.contigs.items():
                self.contigs[_contig] = _info.length

    def find_biallelic_contigs(self):
        """
        Identify which contigs have at least 5 biallelic SNPs and write them to `workdir/vcf.biallelic`
        Populates `self.biallelic` file and `self.biallelic_contigs`
        """
        self.print.log("Finding contigs with ≥ 5 biallelic SNPs", newline=False)

        self.biallelic_file = Path(os.path.join(self.workdir, os.path.basename(self.file) + ".biallelic")).resolve().as_posix()
        if not self.contigs:
            self.get_contigs()
        to_keep = []

        if self.file.lower().endswith("bcf") and not os.path.exists(f"{self.file}.csi"):
            pysam.bcftools.index(self.file)
        if self.file.lower().endswith("vcf.gz") and not os.path.exists(f"{self.file}.tbi"):
            pysam.bcftools.index("--tbi", self.file)

        bcftools = which("bcftools") or "bcftools"

        with open(self.biallelic_file, "w", encoding="utf-8") as f:
            for contig in list(self.contigs.keys()):
                # Use bcftools to count the number of biallelic SNPs in the contig
                viewcmd = subprocess.Popen(
                    [bcftools, 'view', '-H', '-r', str(contig), '-v', 'snps', '-m2', '-M2', '-c', '2', self.file],
                    stdout=subprocess.PIPE,
                    text = True
                )
                snpcount = 0
                keep = False
                while True:
                    # Read the next line of output
                    line = viewcmd.stdout.readline()
                    if not line:
                        break
                    snpcount += 1
                    # If there are at least 5 biallellic snps, terminate the process
                    if snpcount >= 5:
                        viewcmd.terminate()
                        try:
                            viewcmd.communicate(timeout=2)
                        except Exception:
                            viewcmd.kill()
                        f.write(f"{contig}\t{self.contigs[contig]}\n")
                        keep = True
                        break
                if not keep:
                    del self.contigs[contig]
 
        if not self.contigs:
            self.print.validation(False)
            self.print.error("insufficient data", "No contigs with at least 5 biallelic SNPs identified. Cannot continue with imputation.")

        self.print.validation(True)

    def check_phase(self):
        """Check to see if the input VCf file is phased or not, determined by the presence of ID=PS or ID=HP tags"""
        self.print.log("VCF file is phased ([green]PS[/] or [green]HP[/] tags)", newline=False)
        with pysam.VariantFile(self.file) as _vcf:
            formats = list(_vcf.header.formats)
        if 'PS' not in formats and 'HP' not in formats:
            bn = os.path.basename(self.file)
            self.print.validation(False)
            self.print.error(
                "vcf not phased",
                "The input variant file needs to be phased into haplotypes, but no [green]FORMAT/PS[/] or [green]FORMAT/HP[/] fields were found.",
                f"Phase [bold]{bn}[/] into haplotypes using [blue bold]harpy phase[/] or another manner of your choosing and use the phased vcf file as input. If you are confident this file is phased, then the phasing does not follow standard convention and you will need to make sure the phasing information appears as either [green]FORMAT/PS[/] or [green]FORMAT/HP[/] tags."
            )
        self.print.validation(True)

    def match_samples(self, bamlist: list[str], prioritize_vcf: bool) -> None:
        """Validate that the input VCF file and the samples in the list of BAM files. The directionality of this check is determined by 'prioritize_vcf', which prioritizes the sample list in the vcf file, rather than bamlist."""
        self.print.log("Alignment samples in VCF", newline=False)
        with pysam.VariantFile(self.file) as _vcf:
            vcfsamples = list(_vcf.header.samples)
        filesamples = [Path(i).stem for i in bamlist]
        if prioritize_vcf:
            fromthis, query = self.file, vcfsamples
            inthis, search = "the input files", filesamples
        else:
            fromthis, query = "the input files", filesamples
            inthis, search = self.file, vcfsamples
        missing_samples = [x for x in query if x not in search]
        # check that samples in VCF match input directory
        if len(missing_samples) > 0:
            self.print.validation(False)
            self.print.error(
                "mismatched inputs",
                f"There are [bold]{len(missing_samples)}[/] samples found in [blue]{fromthis}[/] that are not in [blue]{inthis}[/]. Terminating Harpy to avoid downstream errors.",
                f"[blue]{fromthis}[/] cannot contain samples that are absent in [blue]{inthis}[/]. Check the spelling or remove those samples from [blue]{fromthis}[/] or remake the vcf file to include/omit these samples. Alternatively, toggle [green]--vcf-samples[/] to aggregate the sample list from the input files or [blue]{self.file}[/].",
                "The samples causing this error are",
                ", ".join(sorted(missing_samples)) + "\n"
            )
        self.samples = query
        self.print.validation(True)

    def match_contigs(self, contigs: list[str]):
        """Check if the supplied contigs are present in the VCF"""
        self.print.log("Input contigs exist in VCF", newline=False)
        bad_names = []
        if not self.contigs:
            self.get_contigs()
        for i in contigs:
            if i not in self.contigs:
                bad_names.append(i)
        if bad_names:
            shortname = os.path.basename(self.file)
            self.print.validation(False)
            self.print.error(
                "contigs absent",
                f"Some of the provided contigs were not found in [blue]{shortname}[/]. This will definitely cause plotting errors in the workflow.",
                "Check that your contig names are correct, including uppercase and lowercase. It's possible that you listed a contig in the genome that isn't in the variant call file due to filtering.",
                f"Contigs absent in {shortname}",
                ",".join([i for i in bad_names])
            )
        self.print.validation(True)

    def validate_region(self, region: str) -> tuple[str,int,int]:
        """
        Validates the --region input of harpy impute to infer it's a properly formatted region
        Use the contigs and lengths of the vcf file to check that the region is valid. Returns
        a tuple of (contig, start, end).
        """
        self.print.log("Input regions exist in VCF", newline=False)
        if not self.contigs:
            self.get_contigs()
        contig, positions = region.split(":")
        parts = positions.split("-")
        startpos, endpos = (int(parts[0]), int(parts[1]))

        # check if the region is in the biallelic contig list
        if contig not in self.contigs:
            self.print.validation(False)
            self.print.error(
                "missing contig",
                f"The [bold yellow]{contig}[/] contig given in [blue]{region}[/] is not in the list of contigs identified to have at least 2 biallelic SNPs, therefore it cannot be processed.",
                f"Restrict the contig provided to [bold green]--region[/] to those with at least 2 biallelic SNPs, which can be reviewed in [blue]{self.biallelic_file}[/]."
            )
        if endpos > self.contigs[contig]:
            self.print.print("[yellow]𐄂[/]")
            self.print.notice(f"The region end position ({endpos}) is greater than the contig size ({self.contigs[contig]}), setting the end position to the end of the contig.")
            endpos = self.contigs[contig]
        else:
            self.print.validation(True)
        return contig, startpos, endpos

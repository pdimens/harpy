
import os
from pathlib import Path
import pysam
import pysam.bcftools
import subprocess
from harpy.common.printing import print_error

class VCF():
    '''
    A class to contain and validate a VCF input file.
    '''
    def __init__(self, filename:str, workdir:str):
        os.makedirs(workdir, exist_ok = True)

        self.file = filename
        self.workdir = workdir

    def find_biallelic_contigs(self):
        """
        Identify which contigs have at least 5 biallelic SNPs and write them to `workdir/vcf.biallelic`
        Populates `self.biallelic` file and `self.biallelic_contigs`
        """
        self.biallelic_file = Path(os.path.join(self.workdir, os.path.basename(self.file) + ".biallelic")).resolve().as_posix()
        self.biallelic_contigs = []
        with pysam.VariantFile(self.file) as _vcf:
            header_contigs = list(_vcf.header.contigs)
        if self.file.lower().endswith("bcf") and not os.path.exists(f"{self.file}.csi"):
            pysam.bcftools.index(self.file)
        if self.file.lower().endswith("vcf.gz") and not os.path.exists(f"{self.file}.tbi"):
            pysam.bcftools.index("--tbi", self.file)

        with open(self.biallelic_file, "w", encoding="utf-8") as f:
            for contig in header_contigs:
                # Use bcftools to count the number of biallelic SNPs in the contig
                viewcmd = subprocess.Popen(['bcftools', 'view', '-H', '-r', str(contig), '-v', 'snps', '-m2', '-M2', '-c', '2', self.file], stdout=subprocess.PIPE)
                snpcount = 0
                while True:
                    # Read the next line of output
                    line = viewcmd.stdout.readline().decode()
                    if not line:
                        break
                    snpcount += 1
                    # If there are at least 2 biallellic snps, terminate the process
                    if snpcount >= 5:
                        viewcmd.terminate()
                        f.write(f"{contig}\n")
                        self.biallelic_contigs.append(contig)
                        break
        if not self.biallelic_contigs:
            print_error("insufficient data", "No contigs with at least 2 biallelic SNPs identified. Cannot continue with imputation.")

    def check_phase(self):
        """Check to see if the input VCf file is phased or not, determined by the presence of ID=PS or ID=HP tags"""
        with pysam.VariantFile(self.file) as _vcf:
            formats = list(_vcf.header.formats)
        if 'PS' not in formats and 'HP' not in formats:
            bn = os.path.basename(self.file)
            print_error(
                "vcf not phased",
                "The input variant file needs to be phased into haplotypes, but no [green]FORMAT/PS[/] or [green]FORMAT/HP[/] fields were found.",
                f"Phase [bold]{bn}[/] into haplotypes using [blue bold]harpy phase[/] or another manner of your choosing and use the phased vcf file as input. If you are confident this file is phased, then the phasing does not follow standard convention and you will need to make sure the phasing information appears as either [green]FORMAT/PS[/] or [green]FORMAT/HP[/] tags."
            )

    def match_samples(self, bamlist: list[str], prioritize_vcf: bool) -> list[str]:
        """Validate that the input VCF file and the samples in the list of BAM files. The directionality of this check is determined by 'prioritize_vcf', which prioritizes the sample list in the vcf file, rather than bamlist."""
        vcfsamples = pysam.bcftools.head(self.file).split("\tINFO\tFORMAT\t")[-1].split()
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
            print_error(
                "mismatched inputs",
                f"There are [bold]{len(missing_samples)}[/] samples found in [blue]{fromthis}[/] that are not in [blue]{inthis}[/]. Terminating Harpy to avoid downstream errors.",
                f"[blue]{fromthis}[/] cannot contain samples that are absent in [blue]{inthis}[/]. Check the spelling or remove those samples from [blue]{fromthis}[/] or remake the vcf file to include/omit these samples. Alternatively, toggle [green]--vcf-samples[/] to aggregate the sample list from the input files or [blue]{self.file}[/].",
                "The samples causing this error are",
                ", ".join(sorted(missing_samples)) + "\n"
            )
        self.samples = query

    def match_contigs(self, contigs: list[str]):
        """Check if the supplied contigs are present in the VCF"""
        with pysam.VariantFile(self.file) as _vcf:
            vcf_contigs = list(_vcf.header.contigs)
        bad_names = []
        for i in contigs:
            if i not in vcf_contigs:
                bad_names.append(i)
        if bad_names:
            shortname = os.path.basename(self.file)
            print_error(
                "contigs absent",
                f"Some of the provided contigs were not found in [blue]{shortname}[/]. This will definitely cause plotting errors in the workflow.",
                "Check that your contig names are correct, including uppercase and lowercase. It's possible that you listed a contig in the genome that isn't in the variant call file due to filtering.",
                f"Contigs absent in {shortname}",
                ",".join([i for i in bad_names])
            )

    def get_contigs(self) -> dict:
        """reads the header of a vcf/bcf file and returns a dict of the contigs (keys) and their lengths (values)"""
        contigs = {}
        with pysam.VariantFile(self.file) as _vcf:
            for _contig,_info in _vcf.header.contigs.items():
                contigs[_contig] = _info.length
        return contigs

    def validate_region(self, region: str) -> list:
        """
        Validates the --region input of harpy impute to infer it's a properly formatted region
        Use the contigs and lengths of the vcf file to check that the region is valid. Returns
        a tuple of (contig, start, end).
        """
        contigs = self.get_contigs()
        contig, positions = region.split(":")
        startpos,endpos,buffer = [int(i) for i in positions.split("-")]
        # check if the region is in the genome
        if contig not in self.biallelic_contigs:
            print_error(
                "missing contig",
                f"The [bold yellow]{contig}[/] contig given in [blue]{region}[/] is not in the list of contigs identified to have at least 2 biallelic SNPs, therefore it cannot be processed.",
                f"Restrict the contig provided to [bold green]--region[/] to those with at least 2 biallelic SNPs. The contigs Harpy found with at least 2 biallelic can be reviewed in [blue]{self.biallelic_file}[/]."
            )
        if endpos > contigs[contig]:
            print_error(
                "invalid region",
                f"The region end position [yellow bold]({endpos})[/] is greater than the length of contig [yellow bold]{contig}[/] ({contigs[contig]})"
            )
        return contig, startpos, endpos, buffer

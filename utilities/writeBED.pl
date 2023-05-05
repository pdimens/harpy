#! /usr/bin/perl -W
use warnings;
use strict;

my $max_distance = 50000;
chomp($ARGV[0]);
my $outfile = $ARGV[1];
#$outfile=~s/\.bam/.all.bed/;
open(OUT, ">$outfile");
#my $header = `samtools view $ARGV[0] -H`;
my $last_bx = "";
my $last_read = "";
my $last_chr = "";
my $sam = "";

my @cache;
my %sort_order;
open (PIPE, "samtools view $ARGV[0] -F 1024|");
while (<PIPE>) {
    my @temp = split "\t";
    #$temp[0] = $1."|".$3."|".$2."//".$temp[0] if (/BX:Z:(\S+).+QX:Z:(\S+).+RX:Z:(\S+)/);
    my $bx = $1 if (/BX:Z:(\S+)/);
    next if ($temp[5] eq "*");
    next if ($bx=~/([ABCD]00){3,4}/);
    if (($bx ne $last_bx) || ($temp[2] ne $last_chr)) {
        my @mols;
        my $last_end = 0;
        if (@cache) {
            my $str = join("\t",@{$cache[0]});
            foreach my $ci (1..$#cache) {
                $str = join("\t",@{$cache[$ci]});
                #print "Comparing $cache[$ci][3] and ".($cache[$ci-1][3])." = ".($cache[$ci][3]-$cache[$ci-1][3])."\n";
                my @last_cigar = split(/(\d+[MIDNSHPX=])/, $cache[$ci-1][5]);
                my $last_a_len = 0;
                #my $length = length($cache[$order[$ix]][9]);
                foreach my $v (@last_cigar) {
                    $last_a_len += $1 if ($v=~/(\d+)[MDNX=]/);
                }
                if ($cache[$ci][3]-($cache[$ci-1][3]+$last_a_len) > $max_distance) {
                    #print "CACHE-CI-3 $cache[$ci][3] - CACHE-CI_LAST-3 $cache[$ci-1][3] + LAST_A_LEN $last_a_len => ".($cache[$ci][3]-($cache[$ci-1][3]+$last_a_len))." vs. MAX_DIST $max_distance\n";
                    push @mols, [ ($last_end, $ci-1) ];
                    $last_end = $ci;
                }
            }
            push @mols, [ ($last_end, $#cache) ];
            #print "NUM_MOL: ".(@mols)."\n";
            
            MOL: foreach my $m (0..$#mols) {
            my @starts;
            my @sizes;
            my @reads;
            my %bx;
            my %bx_name;
            my $mqual=0;
            my $mori=0;
            my %mchr;
            READ: foreach my $ix ($mols[$m][0]..$mols[$m][1]) {
            #print "INDEX $ix from : $mols[$m][0] to $mols[$m][1]\n";
            #Skip duplicates
            next READ if (substr(reverse(sprintf("%012b",$cache[$ix][1])),10,1));
            $mchr{$cache[$ix][2]}++;
            #Parse cigar string;
            my $line = join("\t", @{$cache[$ix]});
            #print $line."";
            my @cigar = split(/(\d+[MIDNSHPX=])/, $cache[$ix][5]);
            my $a_len = 0;
            #my $length = length($cache[$order[$ix]][9]);
            foreach my $v (@cigar) {
                $a_len += $1 if ($v=~/(\d+)[MDNX=]/);
            }
            #Get the starts into a sorted order as it is
            if (!(@starts)) {
                push @starts, $cache[$ix][3]-1;
                push @sizes, $a_len;
                push @reads, $cache[$ix][0];
                } elsif ($cache[$ix][3]-1 >= $starts[0]) {
                push @starts, $cache[$ix][3]-1;
                push @sizes, $a_len;
                push @reads, $cache[$ix][0];
                } else {
                unshift @starts, $cache[$ix][3]-1;
                unshift @sizes, $a_len;
                unshift @reads, $cache[$ix][0];
            }
            $mqual += $cache[$ix][4];
            #Count the orientation flags if it's the first read. -2*0.5*FLAG - revComp gives -1 for -ve and +1 for +ve
            $mori += 2*(0.5-(substr(reverse(sprintf("%012b",$cache[$ix][1])),4,1))) if (substr(reverse(sprintf("%012b",$cache[$ix][1])),6,1));
        }
        $mqual /= @starts if (@starts);
        my $mchr_out;
        CHROM: foreach my $c (sort {$mchr{$b} <=> $mchr{$a}} keys %mchr) {
        $mchr_out = $c;
        last CHROM;
    }
    my $mstart = $starts[0];
    my $mtstart = $starts[0];
    my $mend = $starts[$#starts] + $sizes[$#sizes];
    my $mtend = $starts[$#starts] + $sizes[$#sizes];
    my $mname=$last_bx;
    if ($mori >= 0) {$mori = "+"} else {$mori = "-"};
    my $mcol;
    my ($a, $b, $c, $d) = ($1, $2, $3, $4) if ($last_bx=~/^A(\d+)C(\d+)B(\d+)D(\d+)/);
    $mcol = 96**2 * $a + 96 * $b + $d;
    $mcol = "0,".int($mcol/(96**3-1)*253+0.5).",0";
    my $msizes_str = join(",", @sizes);
    map {$_-=$mstart} @starts;
    my $mstarts_str = join (",", @starts);
    my $mreads_str = join (",", @reads);
    #print "$mchr_out\t$mstart\t$mend\t$mname\t$mqual\t$mori\t$mtstart\t$mtend\t".(@starts)."\t$msizes_str\t$mstarts_str\t$mreads_str\n";
    print OUT "$mchr_out\t$mstart\t$mend\t$mname\t$mqual\t$mori\t$mtstart\t$mtend\t$mcol\t".(@starts)."\t$msizes_str\t$mstarts_str\t$mreads_str\n";
}
}
undef @cache;
};
if (!(@cache)) {
push @cache, [ @temp ];
} elsif (($last_read eq $temp[0])) {
#Checking to see if the read corresponds to first or second read;
if (substr(reverse(sprintf("%012b",$temp[1])),6,1)) {
    push @cache, [ @temp ];
    } else {
    push @cache, [ @temp ];
}
} elsif ($temp[3] >= $cache[0][3]) {
push @cache, [ @temp ];
} else {
unshift @cache, [ @temp ];
}
$last_bx = $bx;
$last_chr = $temp[2];
$last_read = $temp[0];
#print "\t\tPUSHED INTO CACHE - size: ".(@cache)."\n\n";
}
close (PIPE);
my @mols;
my $last_end = 0;
if (@cache) {
foreach my $ci (1..$#cache) {
#print "Comparing $cache[$ci][3] and ".($cache[$ci-1][3])." = ".($cache[$ci][3]-$cache[$ci-1][3])."\n";
my @last_cigar = split(/(\d+[MIDNSHPX=])/, $cache[$ci-1][5]);
my $last_a_len = 0;
#my $length = length($cache[$order[$ix]][9]);
foreach my $v (@last_cigar) {
    $last_a_len += $1 if ($v=~/(\d+)[MDNX=]/);
}
if ($cache[$ci][3]-($cache[$ci-1][3]+$last_a_len) > $max_distance) {
    push @mols, [ ($last_end, $ci-1) ];
    $last_end = $ci;
}
}
push @mols, [ ($last_end, $#cache) ];
#print "NUM_MOL: ".(@mols)."\n";

MOL: foreach my $m (0..$#mols) {
my @starts;
my @sizes;
my @reads;
my %bx;
my %bx_name;
my $mqual=0;
my $mori=0;
my %mchr;
READ: foreach my $ix ($mols[$m][0]..$mols[$m][1]) {
#print "INDEX $ix from : $mols[$m][0] to $mols[$m][1]\n";
#Skip duplicates
next READ if (substr(reverse(sprintf("%012b",$cache[$ix][1])),10,1));
$mchr{$cache[$ix][2]}++;
#Parse cigar string;
my $line = join("\t", @{$cache[$ix]});
#print $line."";
my @cigar = split(/(\d+[MIDNSHPX=])/, $cache[$ix][5]);
my $a_len = 0;
#my $length = length($cache[$order[$ix]][9]);
foreach my $v (@cigar) {
$a_len += $1 if ($v=~/(\d+)[MDNX=]/);
}
#Get the starts into a sorted order as it is
if (!(@starts)) {
push @starts, $cache[$ix][3]-1;
push @sizes, $a_len;
push @reads, $cache[$ix][0];
} elsif ($cache[$ix][3]-1 >= $starts[0]) {
push @starts, $cache[$ix][3]-1;
push @sizes, $a_len;
push @reads, $cache[$ix][0];
} else {
unshift @starts, $cache[$ix][3]-1;
unshift @sizes, $a_len;
unshift @reads, $cache[$ix][0];
}
$mqual += $cache[$ix][4];
#Count the orientation flags if it's the first read. -2*0.5*FLAG - revComp gives -1 for -ve and +1 for +ve
$mori += 2*(0.5-(substr(reverse(sprintf("%012b",$cache[$ix][1])),4,1))) if (substr(reverse(sprintf("%012b",$cache[$ix][1])),6,1));
}
$mqual /= @starts if (@starts);
my $mchr_out;
CHROM: foreach my $c (sort {$mchr{$b} <=> $mchr{$a}} keys %mchr) {
$mchr_out = $c;
last CHROM;
}
my $mstart = $starts[0];
my $mtstart = $starts[0];
my $mend = $starts[$#starts] + $sizes[$#sizes];
my $mtend = $starts[$#starts] + $sizes[$#sizes];
my $mname=$last_bx;
if ($mori >= 0) {$mori = "+"} else {$mori = "-"};
my $mcol;
my ($a, $b, $c, $d) = ($1, $2, $3, $4) if ($last_bx=~/^A(\d+)C(\d+)B(\d+)D(\d+)/);
$mcol = 96**2 * $a + 96 * $b + $d;
$mcol = "0,".int($mcol/(96**3-1)*253+0.5).",0";
my $msizes_str = join(",", @sizes);
map {$_-=$mstart} @starts;
my $mstarts_str = join (",", @starts);
my $mreads_str = join (",", @reads);
#print "$mchr_out\t$mstart\t$mend\t$mname\t$mqual\t$mori\t$mtstart\t$mtend\t".(@starts)."\t$msizes_str\t$mstarts_str\t$mreads_str\n";
print OUT "$mchr_out\t$mstart\t$mend\t$mname\t$mqual\t$mori\t$mtstart\t$mtend\t$mcol\t".(@starts)."\t$msizes_str\t$mstarts_str\t$mreads_str\n";
}
}
undef @cache;
close (OUT);

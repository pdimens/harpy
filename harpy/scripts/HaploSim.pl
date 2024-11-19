####
# LRSIM for Harpy
# Originally from https://github.com/aquaskyline/LRSIM
# This file was forked from commit 9ad7e54 for Harpy (https://github.com/pdimens/harpy)
# Summary:
# This fork removes everything from LRSIM except where it creates linked reads and modifies the inputs/outputs to be friendlier for Harpy
# Detailed explanation:
# LRSIM first created two haplotypes of a genome by introducing variants via SURVIVOR --REMOVED
# Sections that called external tools (samtools, dwgsim, SURVIVOR) --REMOVED
# The input now expects two haplotypes of a genome (-g) that were created separately, and HaploSim.pl will just
# create linked reads from them. No more, no less.
# NOTABLE CHANGES #
# The -g option now takes the outputs from DWGSIM
# the -a option is new and accepts the .fai files associated with the source genome haplotypes
# seq error insertion for barcodes has been skipped
# -r , -u removed
####

use File::Basename;
use strict;
use warnings;
use feature 'state';
use threads;
use threads::shared;
use IO::Handle;
use Getopt::Std;
use Data::Dumper;
use Cwd 'abs_path';
use Math::Random qw(random_poisson random_uniform_integer);
use Inline 'C';

my %fnToBeUnlinkAtExit = ();

&main;
0;

sub main {
    our %opts = (
        h => undef,
        c => undef,
        g => undef,
        a => undef,
        d => 2,
        l => 16,
        p => undef,
        b => undef,
        e => "0.0001,0.0016",
        E => "0.0001,0.0016",
        i => 350,
        s => 35,
        x => 600,
        f => 100,
        t => 1500,
        m => 10,
        z => 8,
        1 => 1000,
        2 => 1,
        3 => 50,
        4 => 1000,
        5 => 1000,
        6 => 10000,
        7 => 100,
        8 => 1000,
        9 => 10000,
        0 => 100
    );
    &usage( \%opts ) if ( @ARGV < 1 );
    getopts( 'hc:g:a:d:l:p:b:e:E:i:s:x:f:t:m:z:1:2:3:4:5:6:7:8:9:0:', \%opts );
    &usage( \%opts ) if ( defined $opts{h} );

    #Check options
    die "Please provide a output prefix with -p\n" if ( not defined $opts{p} );
    die "Output prefix (-p) cannot end with a /\n" if ( $opts{p} =~ /\/$/ );
    die "Please provide a barcodes file with -b\n" if ( not defined $opts{b} );
    die "Barcodes file $opts{b} not found\n"       if ( !-s "$opts{b}" );
    #Check options end

    #Global variables
    #Initialize Log routine
    our @haplotypes = split /,/, $opts{a};
    our @fastafai = split /,/, $opts{g};
    #Global variables end
    #Load barcodes
    our @barcodes                   = ();
    our $barcodesMutexLock : shared = 0;
    our $numBarcodes                = 0;
    &Log("Load barcodes: start");
    open my $fh, "$opts{b}"
        or &LogAndDie("Barcodes file $opts{b} not found");
    @barcodes = <$fh>;
    chomp(@barcodes);
    $numBarcodes = scalar(@barcodes);
    close $fh;
    &Log("Load barcodes: end");
    #Load barcodes end

    our @fragmentSizesList = ();
    our $sizesCount        = 0;
    our $readsPerMolecule  = int( 0.499 + ( $opts{x} * 1000000 ) / ( $opts{t} * 1000 / $opts{d} ) / $opts{m} / $opts{d} );
    &Log("readPairsPerMolecule: $readsPerMolecule");

    # For every Haplotype
    sub simReads {
        $SIG{'INT'} = $SIG{'TERM'} = $SIG{'KILL'} = sub { threads->exit(); };
        my $i = shift;
        &Log("Simulating: haplotype $i");
        if ( -e "$opts{p}.$i.manifest" ) {
            &Log("Haplotype $i simulation already completed");
            return;
        }

        &Log("Load read positions: haplotype $i");
        my @defaultBarcodeQualAry = split //, "AAAFFF" . "K" x ( $opts{l} - 6 );
        my $haplotype             = $haplotypes[$i];
        my %faidx                 = ();
        my @boundary              = ();
        my $genomeSize            = &LoadFaidx( \%faidx, \@boundary, $fastafai[$i] );
        &LogAndDie("Failed loading genome index " . $fastafai[$i])
            if ( $genomeSize == 0 );
        my $readPositionsInFile = mallocAry($genomeSize);
        initAryFF( $readPositionsInFile, $genomeSize );

        if ( -e "$opts{p}.$i.fp" ) {
            &Log("Importing $opts{p}.$i.fp");
            importAry( $readPositionsInFile, "$opts{p}.$i.fp", $genomeSize );
            &Log("Imported $opts{p}.$i.fp");
        }
        else {
            open my $fh, $haplotype
                or &LogAndDie("Error opening $haplotype");
            my $l1;
            my $l2;
            my $l3;
            my $l4;
            my $l5;
            my $l6;
            my $l7;
            my $l8;
            my $newFpos;
            my $fpos = tell($fh);
            &LogAndDie("Fail to tell file position") if $fpos == -1;
            my $failedRegistration = 0;
            my $rt;

            while ( $l1 = <$fh> ) {
                $l2      = <$fh>;
                $l3      = <$fh>;
                $l4      = <$fh>;
                $l5      = <$fh>;
                $l6      = <$fh>;
                $l7      = <$fh>;
                $l8      = <$fh>;
                $newFpos = tell($fh);
                unless ( $l1 =~ /@(\S+)_(\d+)_\d+_\d_\d_\d_\d_\d+:\d+:\d+_\d+:\d+:\d+_\S+\/1/)
                {
                    &LogAndDie("Cannot find correct chromosome and position in $l1.");
                }
                my $gCoord = &GenomeCoord2Idx( \%faidx, "$1", $2 );
                if ( $gCoord < 0 || $gCoord >= $genomeSize ) {
                    &LogAndDie("$1 $2 $gCoord $fpos");
                }
                $rt = writeToPos( $readPositionsInFile, $gCoord, $fpos );
                ++$failedRegistration if $rt == 0;
                $fpos = $newFpos;
            }
            close $fh;
            &Log("Reads failed to load: $failedRegistration ");
            &Log("Exporting: $opts{p}.$i.fp");
            ++$fnToBeUnlinkAtExit{"$opts{p}.$i.fp"};
            exportAry( $readPositionsInFile, "$opts{p}.$i.fp", $genomeSize );
            delete $fnToBeUnlinkAtExit{"$opts{p}.$i.fp"};
            &Log("Exported: $opts{p}.$i.fp");
        }

        open my $outputfh, "> $opts{p}.$i.manifest"
            or &LogAndDie("Error opening: $opts{p}.$i.manifest");
        ++$fnToBeUnlinkAtExit{"$opts{p}.$i.manifest"};

        my $readsCountDown = int( $opts{x} * 1000000 / $opts{d} );
        &Log("Reads remaining: $readsCountDown");

        while ( $readsCountDown > 0 ) {
            #Pick a barcode
            my $selectedBarcode;
            {
                my $idx = int( rand($numBarcodes) );
                lock($barcodesMutexLock);
                my $wentToZero = 0;
                while (1) {
                    if ( $barcodes[$idx] eq "" ) {
                        ++$idx;
                        if ( $idx == $numBarcodes && not($wentToZero) ) {
                            $idx        = 0;
                            $wentToZero = 1;
                        }
                        elsif ( $idx == $numBarcodes && $wentToZero ) {
                            &LogAndDie("Reached end of barcodes list. No more barcodes. Last read processed: $readsCountDown. Exiting.");
                        }
                        next;
                    }
                    $selectedBarcode = $barcodes[$idx];
                    $barcodes[$idx] = "";
                    last;
                }
            }
            my @precreatedSelectedBarcodeAry = split //, $selectedBarcode;
            my $numberOfMolecules = &PoissonMoleculePerPartition( $opts{m} );

            for ( my $j = 0 ; $j < $numberOfMolecules ; ++$j ) {
                #Pick a starting position
                my $startingPosition = int( rand($genomeSize) );
                #Pick a fragment size
                my $moleculeSize = ( $sizesCount == 0 ) ? ( &PoissonMoleculeSize( $opts{f} * 1000 ) ) : ( $fragmentSizesList[ rand($sizesCount) ] );
                my $readsToExtract = int($readsPerMolecule * $moleculeSize / ( $opts{f} * 1000 ) + 0.4999 );
                #Check and align to boundary
                my $lowerBoundary;
                my $upperBoundary;
                &bSearch($startingPosition, \@boundary, \$lowerBoundary, \$upperBoundary);
                if ( ( $startingPosition + $moleculeSize ) > $upperBoundary ){
                    my $newMoleculeSize = $upperBoundary - $startingPosition;
                    #skip molecule with length < 1000
                    if ( $newMoleculeSize < 1000 ){
                        --$j;
                        next;
                    }
                    $readsToExtract = int( $readsToExtract * $newMoleculeSize / $moleculeSize );
                    $moleculeSize = $newMoleculeSize;
                }

                #Get a list of read positions
                my @readPosToExtract = random_uniform_integer( $readsToExtract, $startingPosition, $startingPosition + $moleculeSize - 1 );
                foreach my $gCoord (@readPosToExtract) {
                    my $filePosToExtract = getFromPos( $readPositionsInFile, $gCoord, $genomeSize );
                    next if $filePosToExtract == -1;
                    #Introduce barcode mismatch
                    my @selectedBarcodeAry = @precreatedSelectedBarcodeAry;
                    my @barcodeQualAry     = @defaultBarcodeQualAry;
                    my $barcodeLength      = $opts{l};
                    #Output
                    print $outputfh "$filePosToExtract\t"
                        . ( join "", @selectedBarcodeAry ) . "\t"
                        . ( join "", @barcodeQualAry ) . "\n";

                    --$readsCountDown;
                    if ( $readsCountDown % 100000 == 0 ) {
                        &Log("Reads remaining: $readsCountDown");
                    }
                }
            }
        }
        close $outputfh;
        delete $fnToBeUnlinkAtExit{"$opts{p}.$i.manifest"};
        freeAry($readPositionsInFile);
        if ( !-s "$opts{p}.$i.manifest" ) {
            &LogAndDie("$opts{p}.$i.manifest empty");
        }
    }

    for ( my $i = 0 ; $i < $opts{d} ; ++$i ) {
        simReads($i);
        sleep( 2 + int( rand(3) ) );
    }
}

#Simulate reads end
0;

sub usage {
    my $opts = shift @_;
    die(
        qq/
    Usage:   $0 -r\/-g <reference\/haplotypes> -p <output prefix> [options]

    Reference genome and variants:
    -g STRING   Haploid FASTA .FAI files, separated by comma
    -a STRING   DWGSIM sequences, interleaved and separated by comma
    -p STRING   Output prefix

    Illumina reads characteristics:
    -e FLOAT    Per base error rate of the first read [$$opts{e}]
    -E FLOAT    Per base error rate of the second read [$$opts{E}]
    -i INT      Outer distance between the two ends for pairs [$$opts{i}]
    -s INT      Standard deviation of the distance for pairs [$$opts{s}]

    Linked reads parameters:
    -b STRING   Barcodes list
    -l INT      Barcode Length
    -x INT      # million reads pairs in total to simulated [$$opts{x}]
    -f INT      Mean molecule length in kbp [$$opts{f}]
    -c STRING   Input a list of fragment sizes
    -t INT      n*1000 partitions to generate [$$opts{t}]
    -m INT      Average # of molecules per partition [$$opts{m}]

    Miscellaneous:
    -h          Show this help

    /
    );
}

# Log routine
sub Log {
    my $time = localtime;
    print STDERR "$time: $_[0]\n";
}

sub LogAndDie {
    &Log(@_);
    die;
}

# Log routine end

sub LoadFaidx {
    my $faidx    = shift;
    my $boundary = shift;
    my $fai       = shift;
    open my $fh, "$fai" or &LogAndDie("Error opening faidx: $fai");
    my $accumulation = 0;
    while (<$fh>) {
        chomp;
        my @a = split;
        $$faidx{acc}{"$a[0]"}  = $accumulation;
        $$faidx{size}{"$a[0]"} = $a[1];
        push @$boundary, $accumulation;
        $accumulation += $a[1];
    }
    push @$boundary, $accumulation;
    close $fh;
    return $accumulation;
}

sub getChrSize  { return ${ $_[0] }{size}{ $_[1] }; }
sub getChrStart { return ${ $_[0] }{acc}{ $_[1] }; }

sub GenomeCoord2Idx {
    &LogAndDie("not defined $_[1]") unless defined ${ $_[0] }{acc}{ $_[1] };
    return ${ $_[0] }{acc}{ $_[1] } + $_[2];
}

sub bSearch {
    my ( $elem, $list, $lowerLimit, $upperLimit ) = @_;
    my $max = $#$list;
    my $min = 0;

    my $index;
    while ( $max >= $min ) {
        $index = int( ( $max + $min ) / 2 );
        if    ( $list->[$index] < $elem ) { $min = $index + 1; }
        elsif ( $list->[$index] > $elem ) { $max = $index - 1; }
        else                              { last; }
    }
    if ( $elem >= $list->[$index] ) {
        $$lowerLimit = $list->[$index];
        $$upperLimit = $list->[ $index + 1 ];
    }
    elsif ( $elem < $list->[$index] ) {
        $$lowerLimit = $list->[ $index - 1 ];
        $$upperLimit = $list->[$index];
    }
    else { die "bSearch: Should never reach here"; }
}

sub PoissonMoleculePerPartition {
    state $mu = $_[0];
    state $i  = 10000;
    state $pool;
    $i = 10000 if ( $mu != $_[0] );
    if ( $i == 10000 ) {
        @{$pool} = random_poisson( 10000, $_[0] );
        $i = 0;
    }
    return ${$pool}[ $i++ ];
}

sub PoissonMoleculeSize {
    state $mu = $_[0];
    state $i  = 10000;
    state $pool;
    $i = 10000 if ( $mu != $_[0] );
    if ( $i == 10000 ) {
        @{$pool} = random_poisson( 10000, $_[0] );
        $i = 0;
    }
    return ${$pool}[ $i++ ];
}

0;

__END__
__C__

#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>
#include<string.h>
#include<limits.h>

#define AMP_ON_SLOTS 1
long mallocAry(long size)
{
  long ptr = (long)malloc(size * sizeof(size_t) * AMP_ON_SLOTS);
  if(ptr == (long)NULL)
  {
    fprintf(stderr, "Error allocation, size %zu\n", (size_t)size);
    return 0;
  }
  return ptr;
}

void initAryFF(long pptr, long size)
{
  size_t* ptr = (size_t*)pptr;
  memset(ptr, 0xFF, size * sizeof(size_t) * AMP_ON_SLOTS);
}

void printAry(long pptr, long size)
{
  size_t* ptr = (size_t*)pptr;
  long i;
  for(i = 0; i < size*AMP_ON_SLOTS; ++i)
  {
    fprintf(stderr, "%l\t%lu\n", i, ptr[i]);
  }
}

void importAry(long pptr, char* fn, long size)
{
  void* ptr = (void*)pptr;
  FILE* fh = fopen(fn, "rb");
  fread(ptr, sizeof(size_t), size * AMP_ON_SLOTS, fh);
  fclose(fh);
}

void exportAry(long pptr, char* fn, long size)
{
  size_t* ptr = (size_t*)pptr;
  FILE* fh = fopen(fn, "wb");
  fwrite(ptr, sizeof(size_t), size * AMP_ON_SLOTS, fh);
  fclose(fh);
}

void freeAry(long pptr)
{
  size_t* ptr = (size_t*)pptr;
  free(ptr);
}

#define CHK_PREV_SLOT_LIMIT (3000*AMP_ON_SLOTS)
int writeToPos(long pptr, long pos, long toWrite)
{
  size_t* ptr = (size_t*)pptr;
  int limit = CHK_PREV_SLOT_LIMIT;
  pos = (pos + 1) * AMP_ON_SLOTS - 1;
  while(limit > 0)
  {
    if(ptr[pos] == ULLONG_MAX)
    {
      ptr[pos] = (size_t)toWrite;
      break;
    }
    --pos;
    if(pos < 0) { break; }
    --limit;
  }
  return limit;
}

long getFromPos(long pptr, long pos, long maxSize)
{
  size_t* ptr = (size_t*)pptr;
  int limit = 0;
  size_t result = ULLONG_MAX;
  if(pos >= maxSize) { pos = maxSize - 1; }
  if(pos < 0) { pos = 0; }
  pos = (pos + 1) * AMP_ON_SLOTS - 1;
  while(limit < CHK_PREV_SLOT_LIMIT)
  {
    if(ptr[pos] != ULLONG_MAX)
    {
      result = ptr[pos];
      ptr[pos] = ULLONG_MAX;
      break;
    }
    --pos;
    if(pos < 0) { break; }
    ++limit;
  }
  return (long)result;
}

#!/usr/bin/perl -w
use strict;
use Data::Dumper;
use Getopt::Long;

######################################################################################################################################################

# This script performs the analysis of a TiledC experiment. It depends on the combined bam file outputted by the CCseqBasic pipeline "COMBINED_reported_capture_reads_CS5.bam" (F6 folder), ran using an oligo file without proximity exclusion (see instructions here: http://userweb.molbiol.ox.ac.uk/public/telenius/captureManual/PipeSite.html), and converted to sam format. This file contains all mapped reads containing at least one interaction with a restriction fragment within the tile, post PCR duplicate and ploidy filtering. The scripts also requires a file with the coordinates of the restriction fragments (of correspondong enzyme) in the genome. 

# The script selects reads containing 2 or more interactions between restriction fragments within the tile and outputs a raw contact matrix. This matrix can be ICE normalised using the HiC-Pro package.The script also outputs a report file with some basic statistics and a bed file that specifies the coordinates of each bin number. 

# Example of run command:

# nohup perl /path/Tiled_sam2rawmatrix_MO.pl -sam /path/name_COMBINED_reported_capture_reads_CM5.sam -c 11 -s 29902951 -e 33226736 -r /path/mm9_dpnII_coordinates.txt -n name -b 2000 &

######################################################################################################################################################

&GetOptions
(
	"sam=s" => \my $sam_file, 	    # -sam          sam file
    "c=i" => \my $tile_chr,         # -c            chromosome of tiled region
    "s=i" => \my $tile_start,       # -s            start coordinate of tiled region
    "e=i" => \my $tile_stop,        # -e            end coordinate of tiled region
    "r=s" => \my $dig_genome,       # -r            file with restriction coordinates in genome 
    "n=s" => \my $name,             # -n            name of the experiment, can be whatever you like, will be used for names of output files
    "b=i" => \my $bin_size,         # -b            bin size
);

# Open filehandles and generate output directories and files

open (SAMFH, $sam_file) or die "can't open sam file $sam_file";
open (GENFH, $dig_genome) or die "can't open genome file with restriction coordinates $dig_genome";

my $path = "undefined"; 
if ($sam_file =~ /(.*\/)(\V++)/) {
    $path = $1;
    }
my $dir = "$path/matrix/$name/raw/$bin_size/";
unless (-d "$path/matrix/") {
    mkdir "$path/matrix/";
}
unless (-d "$path/matrix/$name") {
    mkdir "$path/matrix/$name";
}
unless (-d "$path/matrix/$name/raw") {
    mkdir "$path/matrix/$name/raw";
}
unless (-d "$dir") {
    mkdir "$dir";
}

my $report = "$path/Tiled_$name\_analysis_report.txt";
open (REP, ">$report") or die "can't open output report file $report";

# Store restriction fragment coordinates in RE_hash and sort in ascending order; to use for binary search

my %RE_hash;
while (my $line = <GENFH>) {
    chomp $line;
    my ($chr, $start, $stop) = split (/\W/, $line);
    push @{$RE_hash{$chr}}, $start;
}

my @chr_index = keys %RE_hash;
foreach my $chr (@chr_index) {
    @{$RE_hash{$chr}} = sort {$a <=> $b} @{$RE_hash{$chr}};
}

# Map reporter reads in tile to restriction fragments and store in data_hash; both as keys to ensure restriction fragments are only reported once per read, and in array to keep track of order and facilitate multi-way interaction analysis.

my %data_hash; my %counter; 
while (my $line = <SAMFH>) {
    chomp $line;
    unless ($line =~ /^@/) {
        $counter{"01 Number of data lines in sam file"}++;
        my ($name, $flag, $chr, $start, $map_qual, $cigar, $mate_ref, $mate_start, $mate_insert, $seq, $seq_qual, $options) = split (/\t/, $line, 12);
        my $read_name = "undefined";
        if ($name =~ /(.*):PE(.*)$/) {
            $read_name = $1;
        }
        $chr =~ s/chr//; 
        my $stop = 0;
        if ($cigar =~/(\d++)(.*)/) {
            $stop = $1 + $start;                # Use length of mapped read until first indel for mapping to restriction fragment
        }
        if ($options =~ /(.*)CO:Z:(.*)(_CISREP)/) {             # Select reporters in cis  
            if ($stop > $tile_start and $start < $tile_stop) {  # Select interactions in tiled region
                my ($start_frag, $end_frag) = binary_search(\@{$RE_hash{$chr}}, $start, $stop, \%counter);  # Map cis reporters to restriction fragment
                unless ($start_frag eq "error") {
                    unless (exists $data_hash{$read_name}{"reporters"}{"$chr:$start_frag-$end_frag"}) {     # Check mapped restriction fragments are only reported once per read
                        push (@{$data_hash{$read_name}{"reporter_array"}}, "$chr:$start_frag-$end_frag");
                        $data_hash{$read_name}{"reporters"}{"$chr:$start_frag-$end_frag"} = 1;
                    }
                }
            }
        }
    }
}

#print Dumper(\%data_hash);

# Count interactions between restriction fragments in all reads stored in data_hash and store these counts in frag_hash.

my %frag_hash;
foreach my $read_name (keys %data_hash) {
    if (exists $data_hash{$read_name}{"reporter_array"}) {
        my @frags_sorted = sort @{$data_hash{$read_name}{"reporter_array"}};
        if ($#frags_sorted == 0) {                                              # Singlets
            $counter{"02a Number of singlets in cis"}++;
        }
        if ($#frags_sorted == 1) {                                              # Doublets
            $counter{"02b Number of doublets in cis"}++;
            $counter{"02 Total number of interactions in cis"}++;
            unless (exists $frag_hash{$frags_sorted[0]}{$frags_sorted[1]}) {         
                $frag_hash{$frags_sorted[0]}{$frags_sorted[1]} = 1;
            }
            else {
                $frag_hash{$frags_sorted[0]}{$frags_sorted[1]} += 1;                                
            }
        }
        if ($#frags_sorted == 2) {                                              # Triplets  
            $counter{"02c Number of triplets in cis"}++;
            $counter{"02 Total number of interactions in cis"} += 3;
            unless (exists $frag_hash{$frags_sorted[0]}{$frags_sorted[1]}) {         
                $frag_hash{$frags_sorted[0]}{$frags_sorted[1]} = 1;
            }
            else {
                $frag_hash{$frags_sorted[0]}{$frags_sorted[1]} += 1;
            }
            unless (exists $frag_hash{$frags_sorted[0]}{$frags_sorted[2]}) {
                $frag_hash{$frags_sorted[0]}{$frags_sorted[2]} = 1;
            }
            else {
                $frag_hash{$frags_sorted[0]}{$frags_sorted[2]} += 1;
            }
            unless (exists $frag_hash{$frags_sorted[1]}{$frags_sorted[2]}) {
                $frag_hash{$frags_sorted[1]}{$frags_sorted[2]} = 1;
            }
            else {
                $frag_hash{$frags_sorted[1]}{$frags_sorted[2]} += 1;
            }
        }
        if ($#frags_sorted == 3) {                                              # Quadruplets 
            $counter{"02d Number of quadruplets in cis"}++;
            $counter{"02 Total number of interactions in cis"} += 6;
            unless (exists $frag_hash{$frags_sorted[0]}{$frags_sorted[1]}) {
                $frag_hash{$frags_sorted[0]}{$frags_sorted[1]} = 1;
            }
            else {
                $frag_hash{$frags_sorted[0]}{$frags_sorted[1]} += 1;
            }
            unless (exists $frag_hash{$frags_sorted[0]}{$frags_sorted[2]}) {
                $frag_hash{$frags_sorted[0]}{$frags_sorted[2]} = 1;
            }
            else {
                $frag_hash{$frags_sorted[0]}{$frags_sorted[2]} += 1;
            }
            unless (exists $frag_hash{$frags_sorted[0]}{$frags_sorted[3]}) {
                $frag_hash{$frags_sorted[0]}{$frags_sorted[3]} = 1;
            }
            else {
                $frag_hash{$frags_sorted[0]}{$frags_sorted[3]} += 1;
            }
            unless (exists $frag_hash{$frags_sorted[1]}{$frags_sorted[2]}) {
                $frag_hash{$frags_sorted[1]}{$frags_sorted[2]} = 1;
            }
            else {
                $frag_hash{$frags_sorted[1]}{$frags_sorted[2]} += 1;
            }
            unless (exists $frag_hash{$frags_sorted[1]}{$frags_sorted[3]}) {
                $frag_hash{$frags_sorted[1]}{$frags_sorted[3]} = 1;
            }
            else {
                $frag_hash{$frags_sorted[1]}{$frags_sorted[3]} += 1;
            }
            unless (exists $frag_hash{$frags_sorted[2]}{$frags_sorted[3]}) {
                $frag_hash{$frags_sorted[2]}{$frags_sorted[3]} = 1;
            }
            else {
                $frag_hash{$frags_sorted[2]}{$frags_sorted[3]} += 1;
            }
        }
        if ($#frags_sorted == 4) {                                              # Quintuplets
            $counter{"02e Number of quintuplets in cis"}++;
            $counter{"02 Total number of interactions in cis"} += 10;
            unless (exists $frag_hash{$frags_sorted[0]}{$frags_sorted[1]}) {
                $frag_hash{$frags_sorted[0]}{$frags_sorted[1]} = 1;
            }
            else {
                $frag_hash{$frags_sorted[0]}{$frags_sorted[1]} += 1;
            }
            unless (exists $frag_hash{$frags_sorted[0]}{$frags_sorted[2]}) {
                $frag_hash{$frags_sorted[0]}{$frags_sorted[2]} = 1;
            }
            else {
                $frag_hash{$frags_sorted[0]}{$frags_sorted[2]} += 1;
            }
            unless (exists $frag_hash{$frags_sorted[0]}{$frags_sorted[3]}) {
                $frag_hash{$frags_sorted[0]}{$frags_sorted[3]} = 1;
            }
            else {
                $frag_hash{$frags_sorted[0]}{$frags_sorted[3]} += 1;
            }
            unless (exists $frag_hash{$frags_sorted[0]}{$frags_sorted[4]}) {
                $frag_hash{$frags_sorted[0]}{$frags_sorted[4]} = 1;
            }
            else {
                $frag_hash{$frags_sorted[0]}{$frags_sorted[4]} += 1;
            }
            unless (exists $frag_hash{$frags_sorted[1]}{$frags_sorted[2]}) {
                $frag_hash{$frags_sorted[1]}{$frags_sorted[2]} = 1;
            }
            else {
                $frag_hash{$frags_sorted[1]}{$frags_sorted[2]} += 1;
            }
            unless (exists $frag_hash{$frags_sorted[1]}{$frags_sorted[3]}) {
                $frag_hash{$frags_sorted[1]}{$frags_sorted[3]} = 1;
            }
            else {
                $frag_hash{$frags_sorted[1]}{$frags_sorted[3]} += 1;
            }
            unless (exists $frag_hash{$frags_sorted[1]}{$frags_sorted[4]}) {
                $frag_hash{$frags_sorted[1]}{$frags_sorted[4]} = 1;
            }
            else {
                $frag_hash{$frags_sorted[1]}{$frags_sorted[4]} += 1;
            }
            unless (exists $frag_hash{$frags_sorted[2]}{$frags_sorted[3]}) {
                $frag_hash{$frags_sorted[2]}{$frags_sorted[3]} = 1;
            }
            else {
                $frag_hash{$frags_sorted[2]}{$frags_sorted[3]} += 1;
            }
            unless (exists $frag_hash{$frags_sorted[2]}{$frags_sorted[4]}) {
                $frag_hash{$frags_sorted[2]}{$frags_sorted[4]} = 1;
            }
            else {
                $frag_hash{$frags_sorted[2]}{$frags_sorted[4]} += 1;
            }
            unless (exists $frag_hash{$frags_sorted[3]}{$frags_sorted[4]}) {
                $frag_hash{$frags_sorted[3]}{$frags_sorted[4]} = 1;
            }
            else {
                $frag_hash{$frags_sorted[3]}{$frags_sorted[4]} += 1;
            }
        }
    }
}

# Bin intraction counts between individual restriction fragments.

my %bin_hash;
my %annot_hash;
foreach my $X (sort keys %frag_hash) {
    my ($chrX, $startX, $stopX) = split (/:|-/, $X);
    my $midX = ($startX + $stopX) / 2;
    my $rawX = int($midX / $bin_size);
    my $binX = $rawX - (int($tile_start / $bin_size));      # Calculating the bin in this way ensures that bins start at rounded counts, eg 29902000
    $annot_hash{$binX}{$midX} = 1;
    foreach my $Y (sort keys %{$frag_hash{$X}}) {
        my ($chrY, $startY, $stopY) = split (/:|-/, $Y);
        my $midY = ($startY + $stopY) / 2;
        my $rawY = int($midY / $bin_size);
        my $binY = $rawY - (int($tile_start / $bin_size));      
        $bin_hash{$binX}{$binY} += $frag_hash{$X}{$Y};
        $annot_hash{$binY}{$midY} = 1;
    }
}

#print Dumper (\%bin_hash);

# Print matrix

my $output_file = "$dir/tiled_$name\_$bin_size.matrix";
open (OUT, ">$output_file") or die "can't open output file $output_file";
foreach my $X (sort { $a <=> $b } keys %bin_hash ) {
    foreach my $Y (sort { $a <=> $b } keys %{$bin_hash{$X}} ) {
        print OUT "$X\t$Y\t$bin_hash{$X}{$Y}\n";
    }
}

# Print bed file with genomic coordinates of each bin number

my $bed_file = "$dir/tiled_$name\_$bin_size\_coordinates.bed";
open (BED, ">$bed_file") or die "can't open output report file $bed_file";

my $coord = $tile_start;
while ($coord < $tile_stop + $bin_size) {
    my $bin_start = int($coord/$bin_size) * $bin_size;
    my $bin_stop = $bin_start + $bin_size - 1;
    my $bin = int(($coord - $tile_start + 1) / $bin_size);
    print BED "chr$tile_chr\t$bin_start\t$bin_stop\t$bin\n";
    $coord += $bin_size; 
}
 
# Print counts to report file

foreach my $key (sort keys %counter) {
    printf REP "%-8s %s\n", $key, $counter{$key};
}

######################################################################################################################################################

sub binary_search {
    my ($chr_array, $start, $stop, $counter_hash) = @_;
    my $mid_value = ($start + $stop) / 2;
    my $first_pos = 0;
    my $last_pos = scalar @$chr_array - 1; 
    my $counter = 0;
    if (($mid_value < $$chr_array[$first_pos]) or ($mid_value > $$chr_array[$last_pos])) {
        #$$counter_hash{"Binary search error: search outside range of restriction enzyme coordinates"}++;
        return ("error", "error")
    }
    for (my $i = 0; $i < 99; $i++) {
        my $mid_search = int(($first_pos + $last_pos) / 2);
        if ($$chr_array[$mid_search] > $$chr_array[$mid_search + 1]) {
            #$$counter_hash{"Binary search error: restriction enzyme array coordinates not in ascending order"}++;
            return ("error", "error")
        }
        if (($$chr_array[$mid_search] <= $mid_value) and ($$chr_array[$mid_search + 1] > $mid_value)) {    
            if (($$chr_array[$mid_search] <= $start + 2) and ($$chr_array[$mid_search + 1] >= $stop - 2)) {    
                return ($$chr_array[$mid_search], $$chr_array[$mid_search + 1] - 1)
            }
            else {
                #$$counter_hash{"Binary search error: fragment overlaps multiple restriction sites"}++;
                return ("error", "error");
            }
        }       
        elsif ($$chr_array[$mid_search] > $mid_value) {
            $last_pos = $mid_search - 1;
        }    
        elsif ($$chr_array[$mid_search] < $mid_value) {
            $first_pos = $mid_search + 1;
        }
        else {
            #$$counter_hash{"Binary search error: end of loop reached"}++;
        }
    }
    #$$counter_hash{"Binary search error: couldn't map read to fragments"}++;
    return ("error", "error")
}


#!/usr/bin/env perl

use warnings;
use strict;

my $hashing_value = 10000;
my $idr_file = $ARGV[0];
my $file1 = $ARGV[1];
my $file2 = $ARGV[2];
my $file3 = $ARGV[3];
my $output_bed = $ARGV[4];
unless ($ARGV[0] && $ARGV[1] && $ARGV[2] && $ARGV[3]) {
    print STDERR "Usage: perl parse_idr_peaks.pl idr_peak_file peak_file1 peak_file2 [l2fc_cutoff] [l10p_cutoff]\n\n";
    exit;
}

my $l2fc_cutoff = 3;
if (exists $ARGV[4]) {
    $l2fc_cutoff = $ARGV[4];
}
my $l10p_cutoff = 3;
if (exists $ARGV[5]) {
    $l10p_cutoff = $ARGV[5];
}

my $idr_cutoff = 0.01;
if (exists $ARGV[6]) {
    $l10p_cutoff = $ARGV[6];
}

my %idr_cutoffs = ("0.001" => "1000", "0.005" => "955", "0.01" => "830", "0.02" => "705", "0.03" => "632",
    "0.04" => "580", "0.05" => "540", "0.06" => "507", "0.07" => "479", "0.08" => "455", "0.09" => "434",
    "0.1" => "415", "0.2" => "290", "0.3" => "217", "0.4" => "165", "0.5" => "125", "1" => "0");

my %idr_output;
&parse_idr_file($idr_file);

my %idr_region2peaks;
&parse_file($file1);
&parse_file($file2);
&parse_file($file3);

open(OUTBED, ">$output_bed") || die "Cannot open $output_bed for writing!";
for my $idr_region (keys %idr_region2peaks) {
    for my $peak (keys %{$idr_region2peaks{$idr_region}}) {
	my ($chrom, $position, $strand) = split(/\:/, $peak);
	my ($start, $stop) = split(/\-/, $position);
	print OUTBED "$chrom\t$start\t$stop\t.\t.\t$strand\n";
    }
}
close(OUTBED);


sub parse_file {
    my $file = shift;
    open(F,$file);
    for my $line (<F>) {
        chomp($line);

        my @tmp = split(/\t/,$line);
        my $chr = $tmp[0];
        my $start = $tmp[1];
        my $stop = $tmp[2];

        my ($chr2,$pos2,$str,$origpval) = split(/\:/,$tmp[3]);
        my $entropy = $tmp[15];
        my $l2fc = $tmp[11];
        my $l10p = $tmp[10];

        next unless ($l2fc >= 3 && $l10p >= 3);
#       print "chr $chr start $start stop $stop str $str l2 $l2fc l10p $l10p ent $entropy\n";


        my $x = int($start / $hashing_value);
        my $y = int($stop / $hashing_value);

        my %overlapping_idrs;


        for my $i ($x..$y) {
            for my $idr_peak (@{$idr_output{$chr}{$str}{$i}}) {
                my ($ichr,$ipos,$istr,$iidr) = split(/\:/,$idr_peak);
                my ($istart,$istop) = split(/\-/,$ipos);
                next if ($istart >= $stop);
                next if ($istop <= $start);

                $overlapping_idrs{$idr_peak} = $iidr;
            }
        }

	if (scalar(keys %overlapping_idrs) > 0) {

	    if (scalar(keys %overlapping_idrs) > 1) {
#		print STDERR "This should NEVER be hit - peak overlaps with more than one IDR region $line\n".join("\t",keys %overlapping_idrs)."\n---\n";
	    }

	    for my $overlapping_idrpeak (keys %overlapping_idrs) {
#	    my @sorted_idr = keys %overlapping_idrs;
#	    my $overlapping_idrpeak = $sorted_idr[0];
		my ($ichr,$ipos,$istr,$iidr) = split(/\:/,$overlapping_idrpeak);

		if ($iidr >= $idr_cutoffs{$idr_cutoff}) {
		    $idr_region2peaks{$overlapping_idrpeak}{$chr.":".$start."-".$stop.":".$str} = 1;
		}
	    }
	} else {
	    # peak not in IDR list

	}
    }
    close(F);
}

sub parse_idr_file {
    my $idr_file = shift;
    open(ID, $idr_file) || die "Cannot open $idr_file for reading!";
    for my $line (<ID>) {
        chomp($line);
        my @fields = split(/\t/,$line);
        my $chrom = $fields[0];
        my $start = $fields[1];
        my $stop = $fields[2];
        my $strand = $fields[5];
        # my $idr_score = $fields[4];
        my $x = int($start / $hashing_value);
        my $y = int($stop / $hashing_value);
        for my $i ($x..$y) {
            push @{$idr_output{$chrom}{$strand}{$i}},$chrom.":".$start."-".$stop.":".$strand.":"."1000";
        }
    }
    close(ID);
}


sub min {
    my $x = shift;
    my $y = shift;
    if ($x < $y) {
	    return($x);
    } else {
	    return($y);
    }
}
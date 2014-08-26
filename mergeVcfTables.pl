#!/usr/bin/env perl 

use warnings;
use strict;

# First GATK file then sift file

my $file = $ARGV[0];
my $file2 = $ARGV[1];
my @sift;
my @gatk;
my $i = 0;

open(file1, "<$file");
while (<file1>){
	if ($_ =~ m/^CHROM/){
		$_ =~ s/\n//;
		print $_,"\tGene\tEffect\tFunclass\tAA\n";
	}
       unless ($_ =~ m/^CHROM/){
		push (@gatk, $_);
	}
 }
open(file2, "<$file2");
while (<file2>){
        unless($_ =~ m/^#/){
                push (@sift, $_);
       }
}
close (file1);
close (file2);

foreach my $gatkline (@gatk){
	my @array = split ("\t",$sift[$i]);
	$gatkline =~ s/\n//;
	$gatkline = $gatkline . "\t" . $array[2] . "\t" . $array[3] . "\t" . $array[4] . "\t" . $array[5];
	print $gatkline;	
	++$i;
	}

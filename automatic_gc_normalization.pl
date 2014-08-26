#!/usr/bin/env perl 

use warnings;
use strict;

my $bam = shift;
my $fragmentlength = shift;
my $ref = shift;
my $line;
my $i = 0;
my $folder = $bam;
#$path =~ s/\.bam//; 
my @temp = split("/", $folder);
$folder = $temp[0];

my $output = $temp[2];
$output =~ s/.bam$/_gc-corrected\.bam/;

open(ST, "genomeCoverageBed -d -ibam $bam|");
        while ($line = <ST>){

         $line =~ m/(\w+)(\s+)(\d+)(\s+)(\d+)/;
               if ($5 != 0){
                       ++$i;
                        }
        }
system("computeGCBias --bamfile $bam --effectiveGenomeSize $i --genome $ref --fragmentLength $fragmentlength --GCbiasFrequenciesFile $folder/GC_bias.out --plotFileFormat pdf --biasPlot $folder/plots/computeGCBias_plot.pdf");

system("correctGCBias --bamfile $bam --effectiveGenomeSize $i --genome $ref --GCbiasFrequenciesFile $folder/GC_bias.out --correctedFile $folder/bam/$output");


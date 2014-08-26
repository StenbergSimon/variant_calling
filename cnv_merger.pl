#!/usr/bin/env perl


use warnings;
use strict;

open (FILE1, "<$ARGV[0]");
my @table = <FILE1>;
close (FILE1);
my @table_clean;
my @table_fixed;
my @matrix;
my $increment;
my $windowsize;
my @merger;

foreach (@table){
	if ($_ =~ m/(#Increment:)(\d+)/){
		$increment = $2;
		}
	elsif($_ =~ m/(#Windowsize:)(\d+)/){
		$windowsize = $2;
		}
	}

foreach (@table){
        unless ($_ =~ m/^#/){
                push (@table_clean, $_);
                }
	}
my $i = 0;


foreach my $line (@table_clean){
	my @array = split ("\t", $line);
	push(@matrix, \@array);
	}
my @lib; 
while (defined $matrix[$i]){
	push (@lib, $i);
	unless(!defined $matrix[($i-1)]){	
		if (($matrix[$i][0] eq $matrix[($i-1)][0]) && ((($matrix[$i][6])  == ($matrix[($i-1)][6] + $increment))) || (($matrix[$i][6])  == ($matrix[($i-1)][6] + $windowsize))){
			my $new = pop(@lib);
			my $last = pop(@lib);
			$new = $last . "," . $new;
			push(@lib, $new);
			}
		}
	++$i;
	}

$i = 0;
foreach my $line (@lib){
	my @log2;
	my @diff;
	my @start;
	my @stop;
	my @ref_avg;
	my @sam_avg;
	my $name;
		
	my @array = split (",", $line);
	foreach my $ref (@array){
		
		$name = $matrix[$ref][0];
		push (@log2, $matrix[$ref][1]);
		push (@ref_avg, $matrix[$ref][2]);
		push (@sam_avg, $matrix[$ref][3]);
		push (@diff, $matrix[$ref][4]);
		push (@start, $matrix[$ref][5]);
		push (@stop, $matrix[$ref][6]);
		}
	foreach (@stop){
		$_ =~ s/\n//g;
		}
		
	my @start_sort = sort {$a <=> $b} @start;
	my @stop_sort = sort{$a <=> $b} @stop;
	
	my $start = $start_sort[0];
	my $stop = pop(@stop_sort);
	my $length = $stop - $start;	

	my $new_entry = $name . "\t" . average(\@log2) . "\t" . average(\@ref_avg) . "\t" .  average(\@sam_avg) . "\t" . average(\@diff) . "\t" . $start . "\t" . $stop . "\t" . $length . "\n";
	push (@table_fixed, $new_entry);
	}

foreach my $line2 (@table){
        if ($line2 =~ m/^#Chrom\/Contig/){
                $line2 =~ s/\n//g;
                $line2 = $line2 . "\t" . "Region_length\n";
                print $line2;
                }
        elsif($line2 =~ m/^#/){
                print $line2;
                }
        }


foreach (@table_fixed){
	print $_;
	}	






sub average{

        my @array = @{$_[0]};
        my $sum = 0;
        my $mean = 0;

        foreach (@array){

                $sum = $sum + $_ ;

        }

        unless(scalar(@array) == 0){

                $mean = $sum / (scalar(@array));
                return $mean;
                }
        else{


        return 0;
        }
}
		

#!/usr/bin/env perl 

# Calculates how much of a gene that is covered by the CNV, if cnv contains 1 gene and its 100% then score = 1 if two then = 2, if two genes 50% covered score = 1

use warnings;
use strict;
 
open (FILE1, "<$ARGV[0]") || die ("COULD NOT OPEN FILE");
open (FILE2, "<$ARGV[1]") || die ("COULD NOT OPEN FILE");

my @table = <FILE1>;
my @gff = <FILE2>;
close(FILE1);
close(FILE2);
my @gff_clean;
my @table_clean;
my @scf;
my @pos1;
my @pos2;
my @ID;
my @means;
my @geness;

foreach (@gff){
        unless ($_ =~ m/^#/){
                push (@gff_clean, $_);
                }
}
foreach (@table){
        unless ($_ =~ m/^#/){
                push (@table_clean, $_);
                }
}

foreach (@gff_clean){

        chomp $_;

        my @temp = split('\t', $_);

        $temp[0] =~ s/\s//;
        push (@scf, $temp[0]);
        $temp[3] =~ s/\s//;
        push (@pos1, $temp[3]);
        $temp[4] =~ s/\s//;
        push (@pos2, $temp[4]);
        $temp[8] =~ s/ID=//;
        $temp[8] =~ s/\t+//g;
	$temp[8] =~ s/\s+//g;
	$temp[8] =~ s/\n+//g;
	chomp($temp[8]);
	push (@ID, $temp[8]);
}

my $i2 = 0;

foreach (@table_clean){
	my @genes = 0;
	$_ =~ s/\t/;/g;
	$_ =~ s/\s//;
	my @array = split(";", $_);
	my $numgenes = 0;
	my $sum = 0;
	my $cnv_start = $array[5];
	my $cnv_stop = $array[6];
	if(defined ($array[8])){
		@genes = split(',',$array[8]);
		foreach my $gene (@genes){
			$gene =~ s/\n+//g;
			$gene =~ s/\s+//g; 
			my $start;
                         my $stop;
			 my $i = 0;
			 my $amount = 0;
				until(($i == scalar(@ID)) ){
					if( $ID[$i] =~ m/^$gene$/){
						$start = $pos1[$i];
						$stop = $pos2[$i];
						}
					++$i;
					}
			if(defined($start) && defined($stop)){
			if (($start > $cnv_start) && ($stop < $cnv_stop)){
				$amount = 1;
				goto END;
				}
			
			elsif( ($start > $cnv_start) && ($stop > $cnv_stop)){
				my $diff = $stop - $cnv_stop;
				$amount = (($stop - $start) - $diff) / ($stop - $start);
				$amount = 1 - $amount;
				goto END;
				}
			elsif(( $start < $cnv_start) && ($stop < $cnv_stop)){
				my $diff = $cnv_start - $start;
				$amount = (($stop - $start) - $diff) / ($stop - $start);
				$amount = 1 - $amount;
				goto END;
				}
			elsif(( $start < $cnv_start) && ($stop > $cnv_stop)){
				my $diff = ($stop - $start) - ($cnv_stop - $cnv_start);
				$amount = $diff / ($stop - $start);
				$amount = 1 - $amount;
				goto END;
				}
				
			else{	
				die("\nSomething went wrong...\n");
				}
			END:
			$sum = $sum + $amount;
			}
		}
	
	}
	
	if ($sum > 0){
		 #my $mean = $sum/$numgenes;
			 $means[$i2] =  $sum;	
			
		}
	++$i2;

}


foreach my $line2 (@table){
        if ($line2 =~ m/^#Chrom\/Contig/){
                $line2 =~ s/\n//g;
                $line2 = $line2 . "\t" . "Score\n";
                print $line2;
                }
        elsif($line2 =~ m/^#/){
                print $line2;
                }
        }
my $i3 = 0;
foreach(@table_clean){
	$_ =~ s/;/\t/g;
	print $_;
	if (defined $means[$i3]){
		print "\t", $means[$i3];
		}
	print "\n";
	++$i3;
	}


#foreach (@means){
#	print $_;
#	}	

	


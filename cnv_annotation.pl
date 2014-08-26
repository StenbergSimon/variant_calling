#!/usr/bin/env perl 

use warnings;
use strict;

#use script.pl cnv.tsv gff.gff


open (FILE1, "<$ARGV[0]");
open (FILE2, "<$ARGV[1]");

my @table = <FILE1>;
my @gff = <FILE2>;
my @table_clean;
my @gff_clean;
my @scf;
my @pos1;
my @pos2;
my @ID;

close(FILE1);
close(FILE2);

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
        push (@ID, $temp[8]);

}


foreach my $line (@table_clean){
	my $i = 0;
	$line =~ s/\n//g;
	my @array = split('\t',$line);
	#print $array[0], "\n";
		foreach my $gffline (@scf){	
			if($gffline =~ m/$array[0]/){
				#print $array[0], "matched with", $gffline, "\n";
				if( ($array[5] > $pos1[$i] && $array[6] < $pos2[$i]) ||($array[1] < $pos1[$i] && $array[1] > $pos2[$i])){					
			#	print "scaffold match: $array[$i] \n\n";
				$line = $line . "," . $ID[$i];
				++$i;
				next;
				}
				++$i;
			}
			else { ++$i;}
	}
}



foreach my $line2 (@table){
	if ($line2 =~ m/^#Chrom\/Contig/){
		$line2 =~ s/\n//g;
		$line2 = $line2 . "\t" . "Gene_ID\n";
		print $line2;
		}
	elsif($line2 =~ m/^#/){
		print $line2;
		}
	}
foreach my $line1 (@table_clean){
	my $i2 =2;
	my @array = split(',',$line1);
	print $array[0], "\t";
	if (defined $array[1]){
		print $array[1];
		}
	while(defined $array[$i2]){
		print "," , $array[$i2];
		++$i2;
		}	
	print "\n";
	}	

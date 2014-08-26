#!/usr/bin/env perl 

use strict;
use warnings;

# This script calculates VARW variance of the width of the gap of the reads contributing to a variant in the vcf file.
# Use samtools view file.bam | VARW.pl file.vcf > output.vcf


#my $file = $ARGV[0];
my $file2 = $ARGV[0];
my @sam;
my @vcfclean;
my @vcfraw;


# Open and read into memory the important parts of the SAM file
#open(file1, "<$file");
while (<STDIN>){
	unless( $_ =~ m/^@/){
		my @array = split ("\t", $_);
		# Chro \t pos \t CIGAR \t length
		my $length = length ($array[9]);
		my $string = $array[2] . "\t" . $array[3] . "\t" .  $array[5] . "\t" . $length;
		push(@sam, $string);
	}
}

#Open and read VCF file containing all variants
open(file2, "<$file2"); 
while (<file2>){
	push (@vcfraw, $_);
	unless($_ =~ m/^#/){
		push (@vcfclean, $_);
		}
}

#close (file1);
close (file2);

my $addline = "##FORMAT=<ID=VARW,Number=1,Type=Integer,Description=\"Variance of the observed gap in all reads for given variant\">";

foreach (@vcfraw){
	if ($_ =~ m/^#/){
		if ($_ =~ m/^#CHROM/){
			print $addline;
			$_ =~ s/\s+/\t/;
			$_ = "\n" . $_;
		}
		print $_;
	}
}

foreach my $variantline (@vcfclean){
	my @gaps;
	# Foreach line in the vcf i e each variant, split up each line and store chr and pos
	my @array = split ("\t", $variantline);
	my $flag_field = $array[8];
	my $value_field = $array[9];
	my $chr_vcf = $array[0];
	my $pos_vcf = $array[1];
	$value_field =~ s/\n//g;
	$flag_field = $flag_field . ":VARW";
	foreach my $samline (@sam){
		# For each sam line do the same for this particular variant (in loop)
		my @array_sam = split ("\t", $samline);
		my $chr_sam = $array_sam[0];
		if ($chr_vcf eq $chr_sam){
			my $pos_sam = $array_sam[1];
			my $length_sam = $array_sam[3];
			# Check if variant overlaps with read
			if( ($pos_vcf >= $pos_sam ) && ($pos_vcf <= ($pos_sam + $length_sam))){
				my $cigar_sam = $array_sam[2];
				my @cigar_array;
				my $cigar_count = $pos_sam;
				my $gap;
				my $i = 0;
				# Store each cigar part into an array X[M/I/S...]
				until($cigar_sam eq ""){
					$cigar_sam =~ m/(^\d+[A-Z])/;
					my $n = length ($1);
					$1 =~ s/\s//g;
					push (@cigar_array, $1);
					substr($cigar_sam, 0, $n, "");
				}
			until($cigar_count >= $pos_vcf){
				unless(defined($cigar_array[$i])){
					goto SKIP;
					}
				$cigar_array[$i] =~ m/(^\d+)([MIDNSHP=X])/;
				unless($2 eq "I"){
					$cigar_count = $cigar_count + $1;
				}
				my $future_cigar_val = 0;
				my $future_cigar_type = 0;
				my $future_cigar_count = 0;
				if (defined $cigar_array[($i + 1)]){
					$cigar_array[($i + 1)] =~ m/(^\d+)([MIDNSHP=X])/;
					$future_cigar_val = $1;
					$future_cigar_type = $2;
					unless($future_cigar_type eq "I"){
						$future_cigar_count = $cigar_count + $future_cigar_val;
					}	
				}
				if (($cigar_count - 1) == $pos_vcf){
					if (defined($cigar_array[($i+1)])){
						$cigar_array[($i+1)] =~ m/(^\d+)([MIDNSHP=X])/;
						if(($2 eq "I") || ($2 eq "D")){
							$gap = $1;
							push (@gaps, $gap);
						}
						else{
							$gap = 0;
						}
					}
				}
				if (($future_cigar_type eq "D") || ($future_cigar_type eq "I")){
					if (($pos_vcf > ($cigar_count)) && ($pos_vcf < ($future_cigar_count))){
						$gap = $future_cigar_val;
						push(@gaps, $gap);
					}
				}
				++$i;
				
			}
			SKIP:	
		}
	}
}
	my $variance  = variance(\@gaps);
	$value_field = $value_field . ":$variance";
	$array[8] = $flag_field;
	$array[9] = $value_field;
	foreach (@array){
		print $_,"\t";
		}
	print "\n";
	
}


# Subroutines

sub variance{
	
	my @array = @{$_[0]};
	my $n = scalar(@array);
	my $sum1 = 0;
	my $sum2 = 0;
	my $var = "nan";
	my $mean = 0;
	foreach (@array){
		$sum1 = $sum1 + $_;
		}
	if (($sum1 > 0) && ($n > 1)){
		$mean = $sum1 / $n;
                foreach my $line (@array){
                        $sum2 = $sum2 + (($line - $mean) * ($line - $mean));
                        }
                $var = $sum2 / ($n - 1);
                return $var;	
	}
	else{
		return "nan";
	}
}

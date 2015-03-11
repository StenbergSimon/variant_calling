#!/usr/bin/env bash

# Usage:
#
#  ./freebays_pipeline.sh reads1.fastq reads2.fastq
#
#

#Variables:
filename=${1%.*}
pair1=$1
pair2=$2

########################################
##############  SETTINGS  ##############
########################################
#
#Reference fasta:
reference=/home/simons/iontorrent_vcalling/yps128_pacbio_assembly_final_only_ChrI_XVI.fasta
#Refernce bwa index:
bwaindex=/home/simons/iontorrent_vcalling/yps128_pacbio_assembly_final_only_ChrI_XVI.fasta
#Reference annotation:
gff=/home/simons/iontorrent_vcalling/YPS128_PacBio_S288C_LiftOver_BILS_Annotation.gff
#
#Fragment length:
fragmentlength=300
#
#2-bit format Reference genome:
regenome2bit=/home/simons/iontorrent_vcalling/yps128_pacbio_assembly_final_only_ChrI_XVI.2bit
#
#
########################################
########################################
########################################

# Create folders for each analysis

mkdir $filename
mkdir $filename/bam
mkdir $filename/logs
mkdir $filename/picard_temps
mkdir $filename/plots
touch $filename/logs/${filename}.log


########## Replace yps128_pacbio with reference genome name (indexed) (FOUR PLACES)


# Run aligner pipe to samtools for sam>bam that pipes to sort the bam file (1) (1b) (1c)
bwa mem -t 16 -M $bwaindex $pair1 $pair2 | samtools view -bSu - | samtools sort - ${filename}_sorted > $filename/logs/bwa.log 2>&1
mv ${filename}_sorted.bam $filename/bam


# Index .bam file

samtools index $filename/bam/${filename}_sorted.bam

# Output .bam stats
samtools flagstat $filename/bam/${filename}_sorted.bam >> $filename/logs/${filename}.log

# Remove PCR duplicates with Picard (2)

java -Xmx2g -jar ~/bin/MarkDuplicates.jar \
INPUT=$filename/bam/${filename}_sorted.bam \
OUTPUT=$filename/bam/${filename}_sorted_RMDUP.bam \
METRICS_FILE=$filename/logs/${filename}_picard.log \
AS=true \
TMP_DIR=$filename/picard_temps \
VALIDATION_STRINGENCY=LENIENT \
> $filename/logs/MarkDuplicates.log 2>&1

# Remove old file 
rm $filename/bam/${filename}_sorted.bam

# Index 
samtools index $filename/bam/${filename}_sorted_RMDUP.bam 

# Add read groups

java -Xmx2g -jar ~/bin/AddOrReplaceReadGroups.jar \
INPUT=$filename/bam/${filename}_sorted_RMDUP.bam \
OUTPUT=$filename/bam/${filename}_sorted_RMDUP_rg.bam \
TMP_DIR=$filename/picard_temps \
RGLB=$filename \
RGPL=illumina \
RGPU=$filename \
RGSM=$filename \
TMP_DIR=$filename/picard_temps \
> $filename/logs/AddOrReplaceReadGroups.log 2>&1

# Remove old file
rm $filename/bam/${filename}_sorted_RMDUP.bam 

# Index

samtools index $filename/bam/${filename}_sorted_RMDUP_rg.bam

# SNP Calling:

samtools calmd -Arb $filename/bam/${filename}_sorted_RMDUP_rg.bam $reference | \
freebayes \
--fasta-reference $reference \
--ploidy 1 \
-F 0.1 \
-no-complex \
--pooled-continuous \
--stdin \
1>$filename/${filename}.vcf 2>$filename/logs/freebays.error 


# Annotating SNPs
java -Xmx8g -jar ~/bin/snpEff.jar \
-v \
-noStats
yps128_bils \
-treatAllAsProteinCoding Auto \
-no-downstream \
-no-intergenic \
-no-intron \
-no-utr \
-o gatk \
$filename/${filename}.vcf \
> $filename/${filename}_snpeff.vcf

#Filter 

java -Xmx8g -jar ~/bin/SnpSift.jar \
filter \
'QUAL > 10' \
$filename/${filename}_snpeff.vcf \
> $filename/${filename}_filtered_snpeff.vcf

echo $filename/${filename}_filtered_snpeff.vcf >> vcf_list.txt
echo scp $filename/${filename}_filtered_snpeff.vcf simon@genmac33.gen.gu.se:~/ >> scp_vcf_list.txt

#Convert to table
java -Xmx8g -jar ~/bin/SnpSift.jar \
extractFields \
$filename/${filename}_filtered_snpeff.vcf \
CHROM \
POS \
"EFF[*].GENE" \
TYPE \
REF \
ALT \
QUAL \
FILTER \
DP \
AF \
AO \
RO \
QA \
QR \
DPRA \
SRP \
SAP \
AC \
AN \
"EFF[*].EFFECT" \
"EFF[*].FUNCLASS" \
"EFF[*].AA" \
> $filename/${filename}_filtered_variants.table 

echo $filename/${filename}_filtered_variants.table >> table_list.txt
echo scp $filename/${filename}_filtered_variants.table simon@genmac33.gen.gu.se:~/ >> scp_table_list.txt

# Collect GC bias plot

java -Xmx2g -jar ~/bin/CollectGcBiasMetrics.jar \
REFERENCE_SEQUENCE=$reference \
INPUT= $filename/bam/${filename}_sorted_RMDUP_rg.bam \
OUTPUT=$filename/${filename}_gcbias.out \
ASSUME_SORTED=true \
CHART_OUTPUT=$filename/plots/${filename}_gc_bias_plot.pdf \
TMP_DIR=$filename/picard_temps \
VALIDATION_STRINGENCY=LENIENT \

# Correct GC bias

automatic_gc_normalization.pl \
$filename/bam/${filename}_sorted_RMDUP_rg.bam \
$fragmentlength \
$regenome2bit \
1>$filename/logs/${filename}_gc_bias_correction.log \
2>>$filename/logs/${filename}_gc_bias_correction.log

# Collect corrected GC bias plot

java -Xmx2g -jar ~/bin/CollectGcBiasMetrics.jar \
REFERENCE_SEQUENCE=$reference \
INPUT=$filename/bam/${filename}_sorted_RMDUP_rg_gc-corrected.bam \
OUTPUT=$filename/${filename}_gcbias-corrected.out \
ASSUME_SORTED=true \
CHART_OUTPUT=$filename/plots/${filename}_gc_bias-corrected_plot.pdf \
TMP_DIR=$filename/picard_temps \
VALIDATION_STRINGENCY=LENIENT \

# CNV analysis
python ~/git/CNV_pipe/cnv_caller.py -f $filename/bam/${filename}_sorted_RMDUP_rg_gc-corrected.bam -m 1 -w 100 -o ${filename}/CNV_analysis -l ~/git/CNV_pipe/chr_only_I_XVI.list

#mkdir -p $filename/intermediate_files/
#mv $filename/${filename}.intervals $filename/intermediate_files/
#mv $filename/${filename}_filtered_snp.vcf $filename/intermediate_files/
#mv $filename/${filename}_gcbias-corrected.out $filename/intermediate_files/
#mv $filename/${filename}_gcbias.out $filename/intermediate_files/
#mv $filename/${filename}_indels.vcf.gz $filename/intermediate_files/
#mv $filename/${filename}_indels.vcf.gz.tbi $filename/intermediate_files/
#mv $filename/${filename}_indels_filtered.vcf $filename/intermediate_files/
#mv $filename/${filename}_indels_filtered_snpEff.vcf $filename/intermediate_files/
#mv $filename/${filename}_indels_rmf.vcf.gz $filename/intermediate_files/
#mv $filename/${filename}_indels_rmf.vcf.gz.tbi $filename/intermediate_files/
#mv $filename/${filename}_indels_snpEff_sift.table $filename/intermediate_files/
##mv $filename/${filename}_snp.vcf.gz $filename/intermediate_files/
##mv $filename/${filename}_snp.vcf.gz.tbi $filename/intermediate_files/
##mv $filename/${filename}_snp_rmf.vcf $filename/intermediate_files/
##mv $filename/${filename}_filtered_snp_snpEFF.vcf $filename/intermediate_files/
#mv $filename/${filename}_indels_snpEff_gatk.table $filename/intermediate_files/
#mv $filename/cnv_analysis_report.tsv $filename/intermediate_files/
#mv $filename/cnv_analysis_report_merged.tsv $filename/intermediate_files/
#mv $filename/cnv_analysis_report_merged_annotated.tsv $filename/intermediate_files/
#mv $filename/GC_bias.out $filename/intermediate_files/
#mv $filename/${filename}_indels_filtered_VARW.vcf.gz $filename/intermediate_files/
#mv $filename/${filename}_indels_filtered_VARW.vcf.gz.tbi $filename/intermediate_files/
#mv $filename/${filename}_indels_filtered_varw.vcf $filename/intermediate_files/
#mv $filename/${filename}__indels_filtered_snpEff.vcf.idx $filename/intermediate_files/
#mv $filename/${filename}_varscan_snp.vcf.gz $filename/intermediate_files/
#mv $filename/${filename}_varscan_snp_rmf.vcf $filename/intermediate_files/
#mv $filename/${filename}_varscan_snp_rmf_snpeff.vcf $filename/intermediate_files/
#mv $filename/${filename}_varscan_snp_rmf_snpeff_gatk.table $filename/intermediate_files/
#mv $filename/${filename}_varscan_snp_rmf_snpeff_sift.table $filename/intermediate_files/
#mv $filename/${filename}_varscan_snp.vcf.gz.tbi $filename/intermediate_files/
## Done :)
#

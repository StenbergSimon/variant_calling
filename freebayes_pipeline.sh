#!/usr/bin/env bash

# Usage:
#
#  ./freebays_pipeline.sh reads.fastq
#
#

#Variables:
filename=${1%.*}
filenamefull=$1

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
fragmentlength=200
#
#2-bit format Reference genome:
regenome2bit=/home/simons/iontorrent_vcalling/yps128_pacbio_assembly_final_only_ChrI_XVI.2bit
#
#Pacbio reference alignment:
pacbio_bam=/home/simons/iontorrent_vcalling/up_91_3.Founder_003_trimmed.sorted_rmdup_gc-corrected.bam
#
#Concordant SNP's for removal
conc=/home/simons/iontorrent_vcalling/complete_concordance.vcf.gz
#
conc_indel=/home/simons/iontorrent_vcalling/complete_concordance_indels.vcf.gz
#
compl=complete_concordance_varscan.vcf.gz
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
bwa mem -t 16 -M $bwaindex $filenamefull | samtools view -bSu - | samtools sort - ${filename}_sorted > $filename/logs/bwa.log 2>&1
mv ${filename}_sorted.bam $filename/bam


# Index .bam file

samtools index $filename/bam/${filename}_sorted.bam

# Output .bam stats
echo "Raw alignment metrics:" >> $filename/logs/${filename}.log
samtools flagstat $filename/bam/${filename}_sorted.bam >> $filename/logs/${filename}.log

# Remove PCR duplicates with Picard (2)

java -Xmx2g -jar ~/bin/MarkDuplicates.jar \
INPUT=$filename/bam/${filename}_sorted.bam \
OUTPUT=$filename/bam/${filename}_sorted_RMDUP.bam \
METRICS_FILE=$filename/logs/${filename}_picard.log \
AS=true \
REMOVE_DUPLICATES=true \
READ_NAME_REGEX='[a-zA-Z0-9]+:[0-9]:([0-9]+):([0-9]+):([0-9]+)' \
TMP_DIR=$filename/picard_temps \
VALIDATION_STRINGENCY=LENIENT \
> $filename/logs/MarkDuplicates.log 2>&1

samtools index $filename/bam/${filename}_sorted_RMDUP.bam 

# Add read groups
# Not really good for anything, GATK just wants it

java -Xmx2g -jar ~/bin/AddOrReplaceReadGroups.jar \
INPUT=$filename/bam/${filename}_sorted_RMDUP.bam \
OUTPUT=$filename/bam/${filename}_sorted_RMDUP_readgroup_realigned_BAQ.bam \
TMP_DIR=$filename/picard_temps \
RGLB=123 \
RGPL=ION_TORRENT \
RGPU=123 \
RGSM=123 \
TMP_DIR=$filename/picard_temps \
> $filename/logs/AddOrReplaceReadGroups.log 2>&1

# Index

samtools index $filename/bam/${filename}_sorted_RMDUP_readgroup_realigned_BAQ.bam

# Metrics

java -Xmx4g -jar ~/bin/BamIndexStats.jar \
I=$filename/bam/${filename}_sorted_RMDUP_readgroup_realigned_BAQ.bam \
#VALIDATION_STRINGENCY=LENIENT \
TMP_DIR=$filename/picard_temps \
> $filename/logs/rmdup_picard.log

# Index
samtools index $filename/bam/${filename}_sorted_RMDUP_realigned_BAQ.bam

# Indel Calling
samtools mpileup \
-f yps128_pacbio_assembly_final.fasta \
-h 100 \
$filename/bam/${filename}_sorted_RMDUP_readgroup_realigned_BAQ.bam \
| \
java -jar /home/simons/bin/VarScan.v2.3.6.jar mpileup2indel --min-reads2 10 --output-vcf 1> $filename/${filename}_indels.vcf


# SNP calling using VarScan:
samtools mpileup \
-f yps128_pacbio_assembly_final.fasta $filename/bam/${filename}_sorted_RMDUP_readgroup_realigned_BAQ.bam \
| java -jar ~/bin/VarScan.v2.3.6.jar mpileup2snp --min-coverage 5 --min-var-freq 0.01 --output-vcf > $filename/${filename}_varscan_snp.vcf
bgzip $filename/${filename}_varscan_snp.vcf
tabix -p vcf $filename/${filename}_varscan_snp.vcf.gz

vcf-isec -c -o $filename/${filename}_varscan_snp.vcf.gz $compl > $filename/${filename}_varscan_snp_rmf.vcf

java -Xmx8g -jar ~/bin/snpEff.jar \
-v \
yps128_bils \
-treatAllAsProteinCoding Auto \
-no-downstream \
-no-intergenic \
-no-intron \
-no-utr \
$filename/${filename}_varscan_snp_rmf.vcf \
> $filename/${filename}_varscan_snp_rmf_snpeff.vcf

java -jar ~/bin/GenomeAnalysisTK.jar \
-R $reference \
-T VariantsToTable \
-V $filename/${filename}_varscan_snp_rmf_snpeff.vcf \
-F CHROM \
-F POS \
-F ID \
-F REF \
-F ALT \
-GF DP \
-GF SDP \
-GF RD \
-GF AD \
-GF FREQ \
-GF GQ \
-GF PVAL \
-GF RBQ \
-GF ABQ \
-GF RDF \
-GF RDR \
-GF ADF \
-GF ADR \
-o $filename/${filename}_varscan_snp_rmf_snpeff_gatk.table

java -Xmx8g -jar ~/bin/SnpSift.jar \
extractFields \
-V $filename/${filename}_varscan_snp_rmf_snpeff.vcf \
CHROM \
POS \
"EFF[*].GENE" \
"EFF[*].EFFECT" \
"EFF[*].FUNCLASS" \
"EFF[*].AA" \
> $filename/${filename}_varscan_snp_rmf_snpeff_sift.table

mergeVcfTables.pl $filename/${filename}_varscan_snp_rmf_snpeff_gatk.table $filename/${filename}_varscan_snp_rmf_snpeff_sift.table > $filename/${filename}_varscan_snp_rmf_snpeff_final.table

# Preparing:
bgzip $filename/${filename}_indels.vcf
tabix -p vcf $filename/${filename}_indels.vcf.gz

# Remove founder variants from indels

vcf-isec -c -o $filename/${filename}_indels.vcf.gz $conc_indel | bgzip -c > $filename/${filename}_indels_rmf.vcf.gz

# Filtering indels

tabix -p vcf $filename/${filename}_indels_rmf.vcf.gz
~/bin/bcftools/bcftools filter -o $filename/${filename}_indels_filtered.vcf -i 'FMT/AD/(FMT/RD+FMT/AD)>0.7' $filename/${filename}_indels_rmf.vcf.gz

# Calculate and filter on VARW
samtools view $filename/bam/${filename}_sorted_RMDUP_readgroup_realigned_BAQ.bam | VARW.pl $filename/${filename}_indels_filtered.vcf > $filename/${filename}_indels_filtered_VARW.vcf
bgzip $filename/${filename}_indels_filtered_VARW.vcf
tabix -p vcf $filename/${filename}_indels_filtered_VARW.vcf.gz
~/bin/bcftools/bcftools filter -o $filename/${filename}_indels_filtered_varw.vcf -i 'FMT/VARW == 0' $filename/${filename}_indels_filtered_VARW.vcf.gz



mv snpEff_genes.txt $filename/logs
mv snpEff_summary.html $filename/logs

# Annotation of indel .vcf file

java -Xmx8g -jar ~/bin/snpEff.jar \
-v \
yps128_bils \
-treatAllAsProteinCoding Auto \
-no-downstream \
-no-intergenic \
-no-intron \
-no-upstream \
-no-utr \
$filename/${filename}_indels_filtered_varw.vcf \
> $filename/${filename}_indels_filtered_snpEff.vcf


# Convert indel vcf to table
java -jar ~/bin/GenomeAnalysisTK.jar \
-R $reference \
-T VariantsToTable \
-V $filename/${filename}_indels_filtered_snpEff.vcf \
-F CHROM \
-F POS \
-F ID \
-F REF \
-F ALT \
-GF DP \
-GF SDP \
-GF RD \
-GF AD \
-GF GQ \
-GF PVAL \
-GF RBQ \
-GF ABQ \
-GF RDF \
-GF RDR \
-GF ADF \
-GF ADR \
-GF VARW \
-o $filename/${filename}_indels_snpEff_gatk.table

java -Xmx8g -jar ~/bin/SnpSift.jar \
extractFields \
$filename/${filename}_indels_filtered_snpEff.vcf \
CHROM \
POS \
"EFF[*].GENE" \
"EFF[*].EFFECT" \
"EFF[*].FUNCLASS" \
"EFF[*].AA" \
> $filename/${filename}_indels_snpEff_sift.table

mergeVcfTables.pl $filename/${filename}_indels_snpEff_gatk.table $filename/${filename}_indels_snpEff_sift.table $filename/${filename}_indels_snpEff_gatk.table > $filename/${filename}_indels_filtered_snpEff_final.table


mv snpEff_genes.txt snpEff_genes_indel.txt
mv snpEff_summary.html snpEff_summary_indel.html

mv snpEff_genes_indel.txt $filename/logs
mv snpEff_summary_indel.html $filename/logs


# Collect GC bias plot

java -Xmx2g -jar ~/bin/CollectGcBiasMetrics.jar \
REFERENCE_SEQUENCE=$reference \
INPUT=$filename/bam/${filename}_sorted_RMDUP.bam \
OUTPUT=$filename/${filename}_gcbias.out \
ASSUME_SORTED=true \
CHART_OUTPUT=$filename/plots/${filename}_gc_bias_plot.pdf \
TMP_DIR=$filename/picard_temps \
VALIDATION_STRINGENCY=LENIENT \

# Correct GC bias

automatic_gc_normalization.pl \
$filename/bam/${filename}_sorted_RMDUP.bam \
$fragmentlength \
$regenome2bit \
1>$filename/logs/${filename}_gc_bias_correction.log \
2>>$filename/logs/${filename}_gc_bias_correction.log

# Collect corrected GC bias plot

java -Xmx2g -jar ~/bin/CollectGcBiasMetrics.jar \
REFERENCE_SEQUENCE=$reference \
INPUT=$filename/bam/${filename}_sorted_RMDUP_gc-corrected.bam \
OUTPUT=$filename/${filename}_gcbias-corrected.out \
ASSUME_SORTED=true \
CHART_OUTPUT=$filename/plots/${filename}_gc_bias-corrected_plot.pdf \
TMP_DIR=$filename/picard_temps \
VALIDATION_STRINGENCY=LENIENT \

# CNV analysis
python ~/git/CNV_pipe/cnv_caller.py -f ${filename}/bam/${filename}_sorted_RMDUP_gc-corrected.bam -m 1 -w 100 -o ${filename}/CNV_analysis -l ~/git/CNV_pipe/chr_only_I_XVI.list

mkdir -p $filename/intermediate_files/
mv $filename/${filename}.intervals $filename/intermediate_files/
mv $filename/${filename}_filtered_snp.vcf $filename/intermediate_files/
mv $filename/${filename}_gcbias-corrected.out $filename/intermediate_files/
mv $filename/${filename}_gcbias.out $filename/intermediate_files/
mv $filename/${filename}_indels.vcf.gz $filename/intermediate_files/
mv $filename/${filename}_indels.vcf.gz.tbi $filename/intermediate_files/
mv $filename/${filename}_indels_filtered.vcf $filename/intermediate_files/
mv $filename/${filename}_indels_filtered_snpEff.vcf $filename/intermediate_files/
mv $filename/${filename}_indels_rmf.vcf.gz $filename/intermediate_files/
mv $filename/${filename}_indels_rmf.vcf.gz.tbi $filename/intermediate_files/
mv $filename/${filename}_indels_snpEff_sift.table $filename/intermediate_files/
#mv $filename/${filename}_snp.vcf.gz $filename/intermediate_files/
#mv $filename/${filename}_snp.vcf.gz.tbi $filename/intermediate_files/
#mv $filename/${filename}_snp_rmf.vcf $filename/intermediate_files/
#mv $filename/${filename}_filtered_snp_snpEFF.vcf $filename/intermediate_files/
mv $filename/${filename}_indels_snpEff_gatk.table $filename/intermediate_files/
mv $filename/cnv_analysis_report.tsv $filename/intermediate_files/
mv $filename/cnv_analysis_report_merged.tsv $filename/intermediate_files/
mv $filename/cnv_analysis_report_merged_annotated.tsv $filename/intermediate_files/
mv $filename/GC_bias.out $filename/intermediate_files/
mv $filename/${filename}_indels_filtered_VARW.vcf.gz $filename/intermediate_files/
mv $filename/${filename}_indels_filtered_VARW.vcf.gz.tbi $filename/intermediate_files/
mv $filename/${filename}_indels_filtered_varw.vcf $filename/intermediate_files/
mv $filename/${filename}__indels_filtered_snpEff.vcf.idx $filename/intermediate_files/
mv $filename/${filename}_varscan_snp.vcf.gz $filename/intermediate_files/
mv $filename/${filename}_varscan_snp_rmf.vcf $filename/intermediate_files/
mv $filename/${filename}_varscan_snp_rmf_snpeff.vcf $filename/intermediate_files/
mv $filename/${filename}_varscan_snp_rmf_snpeff_gatk.table $filename/intermediate_files/
mv $filename/${filename}_varscan_snp_rmf_snpeff_sift.table $filename/intermediate_files/
mv $filename/${filename}_varscan_snp.vcf.gz.tbi $filename/intermediate_files/
# Done :)


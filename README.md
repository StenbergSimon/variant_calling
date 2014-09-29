variant_calling
===============

Scripts used to run variant calling in an automated way


Dependencies:
------------

* Samtools
* Freebayes
* Bamtools
* Deeptools
* Picard
* GATK
* snpEFF
* BWA
* vcftools
* bgzip
* tabix
* CNV_pipe

Known Issues:
------------

Conversion to table from indel .vcf does not work properly

Two realignments (Indel Realignment and BAQ) are currently in place, one will be removed after testing with reference dataset

#! /bin/bash

for i in vcflib samtools bcftools htslib bam2fastq freebayes; do
    cd $i
    make
    cd ..
done

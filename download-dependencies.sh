#! /bin/sh

git clone --recursive https://github.com/ekg/vcflib.git
git clone https://github.com/lh3/samtools.git
git clone https://github.com/samtools/bcftools.git
git clone https://github.com/samtools/htslib.git
git clone --recursive https://github.com/jts/bam2fastq.git
git clone --recursive https://github.com/ekg/freebayes.git

# Download dbsnp
mkdir dbsnp
cd dbsnp
wget ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b142_GRCh37p13/VCF/00-All.vcf.gz
wget ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b142_GRCh37p13/VCF/00-All.vcf.gz.tbi
cd ..

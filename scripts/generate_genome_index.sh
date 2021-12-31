#!/bin/bash

#STAR_Generate
#GRIP-seq pipeline
#Author: Hongijang Liu
#Email: hongjiang.liu@ucsf.edu

# Exit the script if an error happens
set -e

if [[ $# -eq 0 ]]; then
  echo "Error: GENOME_DIR must be provided as an input argument."
  exit 1
fi

GENOME_DIR="$1"
CPU_THREADS=${2:-"8"}  # Default using 8 threads
INDEX_LENGTH=${3:-"100"}  # Default using 100 as the index length

echo "Genome dir: ${GENOME_DIR} ..."
echo "--runThreadN: ${CPU_THREADS} ..."
echo "--sjdbOverhang: ${INDEX_LENGTH} ..."

# activate enviroment:
source ~/.bashrc

mkdir -p ${GENOME_DIR}/hg19

echo "Downloading Fasta File to ${GENOME_DIR} ..."
wget -P ${GENOME_DIR}/hg19 https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/GRCh37.p13.genome.fa.gz

echo "Unzipping Fasta File"
gzip -d ${GENOME_DIR}/hg19/GRCh37.p13.genome.fa.gz

echo "Downloading GTF File to ${GENOME_DIR} ..."
wget -P ${GENOME_DIR}/hg19 https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.chr_patch_hapl_scaff.annotation.gtf.gz

echo "Unzipping GTF File"
gzip -d ${GENOME_DIR}/hg19/gencode.v19.chr_patch_hapl_scaff.annotation.gtf.gz

echo "Running STAR to generate genome index ..."
mkdir -p ${GENOME_DIR}/STAR_index
STAR --runThreadN ${CPU_THREADS} \
--runMode genomeGenerate \
--genomeDir ${GENOME_DIR}/STAR_index \
--genomeFastaFiles ${GENOME_DIR}/hg19/GRCh37.p13.genome.fa \
--sjdbGTFfile ${GENOME_DIR}/hg19/gencode.v19.chr_patch_hapl_scaff.annotation.gtf \
--sjdbOverhang ${INDEX_LENGTH}

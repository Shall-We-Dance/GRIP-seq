#!/bin/bash

# Preprocess

# GRIP-seq pipeline
# Author: Hongijang Liu
# Email: hongjiang.liu@ucsf.edu

# Exit the script if an error happens
set -e

if [[ $# -eq 0 ]]; then
  echo "Error: GENOME_DIR must be provided as an input argument."
  exit 1
fi

REPEAT=$(cat ./repeat.txt)

ROOTDIR="$1"
GENOMEDIR_STAR="$2"
ID=${3:-"current"}  # Default ID is current
THREAD=${4:-"8"}  # Default using 8 threads


echo "ROOTDIR: ${ROOTDIR} ..."
echo "GENOMEDIR_STAR: ${GENOMEDIR_STAR} ..."
echo "REPEAT: ${REPEAT} ..."
echo "ID: ${ID} ..."
echo "Using CPU threads: ${THREAD} ..."

mkdir -p ${ROOTDIR}/fastp/${ID}
mkdir -p ${ROOTDIR}/STAR/${ID}

echo "Preprocessing the fastq files..."

for NAME in ${REPEAT};
do
# fastp
echo "Preprocessing the fastq files for ${NAME}..."
## Read length filter
fastp -i ${ROOTDIR}/raw_data/${NAME}_L003_R1_001.fastq.gz \
-I ${ROOTDIR}/raw_data/${NAME}_L003_R2_001.fastq.gz \
-o ${ROOTDIR}/fastp/${ID}/${NAME}_R1_length.fastq.gz \
-O ${ROOTDIR}/fastp/${ID}/${NAME}_R2_length.fastq.gz \
-h ${ROOTDIR}/fastp/${ID}/${NAME}_fastp_length.html \
-j ${ROOTDIR}/fastp/${ID}/${NAME}_fastp_length.json \
--thread ${THREAD} -Q -A \
--length_required 101 \
--length_limit 101
## quality control and 

illumina adapters
fastp -i ${ROOTDIR}/fastp/${ID}/${NAME}_R1_length.fastq.gz \
-I ${ROOTDIR}/fastp/${ID}/${NAME}_R2_length.fastq.gz \
-o ${ROOTDIR}/fastp/${ID}/${NAME}_R1_fastp_dedup_adapter.fastq.gz \
-O ${ROOTDIR}/fastp/${ID}/${NAME}_R2_fastp_dedup_adapter.fastq.gz \
-h ${ROOTDIR}/fastp/${ID}/${NAME}_fastp_dedup_adapter.html \
-j ${ROOTDIR}/fastp/${ID}/${NAME}_fastp_dedup_adapter.json \
--thread ${THREAD} \
--dedup

## cut GRIP-seq own adapters
fastp -i ${ROOTDIR}/fastp/${ID}/${NAME}_R1_fastp_dedup_adapter.fastq.gz \
-I ${ROOTDIR}/fastp/${ID}/${NAME}_R2_fastp_dedup_adapter.fastq.gz \
-o ${ROOTDIR}/fastp/${ID}/${NAME}_R1_fastp_final.fastq.gz \
-O ${ROOTDIR}/fastp/${ID}/${NAME}_R2_astp_final.fastq.gz \
-h ${ROOTDIR}/fastp/${ID}/${NAME}_fastp_final.html \
-j ${ROOTDIR}/fastp/${ID}/${NAME}_fastp_final.json \
--thread ${THREAD} \
-A \
-f 12 -t 10 -F 10 -T 12
done

#STAR
echo "Mapping reads to genome..."
for NAME in ${REPEAT};
do
echo "Mapping for ${NAME}..."
STAR  --runThreadN ${THREAD} \
--genomeDir  ${GENOMEDIR_STAR} \
--readFilesCommand zcat \
--readFilesIn  ${ROOTDIR}/fastp/${ID}/${NAME}_R2_astp_final.fastq.gz \
--outFileNamePrefix ${ROOTDIR}/STAR/${ID}/${NAME}/${NAME}.bam \
--outSAMtype BAM SortedByCoordinate ;
done

echo "Indexing and calculating depth..."
for NAME in ${REPEAT};
do
samtools index -b ${ROOTDIR}/STAR/${ID}/${NAME}/${NAME}.bamAligned.sortedByCoord.out.bam;
samtools depth ${ROOTDIR}/STAR/${ID}/${NAME}/${NAME}.bamAligned.sortedByCoord.out.bam > ${ROOTDIR}/STAR/${ID}/${NAME}/${NAME}.coverage
cat ${ROOTDIR}/STAR/${ID}/${NAME}/${NAME}.coverage | awk '$3 > 10 {print $0}' > ${ROOTDIR}/STAR/${ID}/${NAME}/${NAME}.positive.10.coverage
done
echo "Finished!"
exit 0

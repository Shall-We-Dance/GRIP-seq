#!/bin/bash

#Preprocess

#GRIP-seq pipeline
#Author: Hongijang Liu
#Email: hongjiang.liu@ucsf.edu

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
mkdir -p ${ROOTDIR}/cutadapt/${ID}
mkdir -p ${ROOTDIR}/STAR/${ID}

echo "Preprocessing the fastq files..."

for NAME in ${REPEAT};
do
#fastp
echo "Preprocessing the fastq files for ${NAME}..."
fastp -i ${ROOTDIR}/raw_data/${NAME}_L003_R1_001.fastq.gz \
-I ${ROOTDIR}/raw_data/${NAME}_L003_R2_001.fastq.gz \
-o ${ROOTDIR}/fastp/${ID}/${NAME}_R1_length.fastq.gz \
-O ${ROOTDIR}/fastp/${ID}/${NAME}_R2_length.fastq.gz \
-h ${ROOTDIR}/fastp/${ID}/${NAME}_fastp_length.html \
-j ${ROOTDIR}/fastp/${ID}/${NAME}_fastp_length.json \
--thread ${THREAD} -Q -A \
--length_required 101 \
--length_limit 101

fastp -i ${ROOTDIR}/fastp/${ID}/${NAME}_R1_length.fastq.gz \
-I ${ROOTDIR}/fastp/${ID}/${NAME}_R2_length.fastq.gz \
-o ${ROOTDIR}/fastp/${ID}/${NAME}_R1_fastp_dedup_adapter.fastq.gz \
-O ${ROOTDIR}/fastp/${ID}/${NAME}_R2_fastp_dedup_adapter.fastq.gz \
-h ${ROOTDIR}/fastp/${ID}/${NAME}_fastp_dedup_adapter.html \
-j ${ROOTDIR}/fastp/${ID}/${NAME}_fastp_dedup_adapter.json \
--thread ${THREAD} \
--dedup

#cutadapt
cutadapt -j ${THREAD} -q 10,10 -u 12 -u -10 -o ${ROOTDIR}/cutadapt/${ID}/${NAME}_R1_cutadapt_fastp_dedup_adapter.fastq.gz ${ROOTDIR}/fastp/${ID}/${NAME}_R1_fastp_dedup_adapter.fastq.gz;
cutadapt -j ${THREAD} -q 10,10 -u 10 -u -12 -o ${ROOTDIR}/cutadapt/${ID}/${NAME}_R2_cutadapt_fastp_dedup_adapter.fastq.gz ${ROOTDIR}/fastp/${ID}/${NAME}_R2_fastp_dedup_adapter.fastq.gz;
done

#STAR
echo "Mapping..."
ulimit -n 65535
for NAME in ${REPEAT};
do
echo "Mapping for ${NAME}..."
STAR  --runThreadN ${THREAD} \
--genomeDir  ${GENOMEDIR_STAR} \
--readFilesCommand zcat \
--readFilesIn  ${ROOTDIR}/cutadapt/${ID}/${NAME}_R2_cutadapt_fastp_dedup_adapter.fastq.gz \
--outFileNamePrefix ${ROOTDIR}/STAR/${ID}/${NAME}/${NAME}.bam \
--outSAMtype BAM SortedByCoordinate ;
done

echo "Indexing and calculating depth..."
for NAME in ${REPEAT};
do
samtools index -b ${ROOTDIR}/STAR/${ID}/${NAME}/${NAME}.bamAligned.sortedByCoord.out.bam;
samtools depth -l 10 ${ROOTDIR}/STAR/${ID}/${NAME}/${NAME}.bamAligned.sortedByCoord.out.bam > ${ROOTDIR}/STAR/${ID}/${NAME}/${NAME}.10.coverage
done
echo "Finished!"
exit 0

#!/bin/bash

#Preprocess

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

REPEAT = WS013 WS014 WS015 WS016

ROOTDIR="$1"
GENOMEDIR_STAR="$2"
ID=${3:-"current"}  # Default ID is current
THREAD=${4:-"8"}  # Default using 8 threads


echo "ROOTDIR: ${ROOTDIR} ..."
echo "GENOMEDIR_STAR: ${GENOMEDIR_STAR} ..."
echo "ID: ${ID} ..."
echo "Using CPU threads: ${THREAD} ..."

mkdir -p ${ROOTDIR}/fastp/${ID}
mkdir -p ${ROOTDIR}/cutadapt/${ID}
mkdir -p ${ROOTDIR}/STAR/${ID}
mkdir -p ${ROOTDIR}/bigWig/${ID}
mkdir -p ${ROOTDIR}/macs2/${ID}
mkdir -p ${ROOTDIR}/output/${ID}
mkdir -p ${ROOTDIR}/find_peak/${ID}
mkdir -p ${ROOTDIR}/upsetplot/${ID}
mkdir -p ${ROOTDIR}/metaPlotR/${ID}


print("Preprocessing the fastq files...")
for NAME in ${REPEAT};
do
#fastp
fastp -i ${ROOTDIR}/raw_data/${NAME}_L003_R1_001.fastq.gz \
-I ${ROOTDIR}/raw_data/${NAME}_L003_R2_001.fastq.gz \
-o ${ROOTDIR}/fastp/${ID}/${NAME}_R1_length.fastq.gz \
-O ${ROOTDIR}/fastp/${ID}/${NAME}_R2_length.fastq.gz \
-h ${ROOTDIR}/fastp/${ID}/${NAME}_fastp_length.html \
-j ${ROOTDIR}/fastp/${ID}/${NAME}_fastp_length.json \
--thread 16 -Q -A \
--length_required 101 \
--length_limit 101

fastp -i ${ROOTDIR}/fastp/${ID}/${NAME}_R1_length.fastq.gz \
-I ${ROOTDIR}/fastp/${ID}/${NAME}_R2_length.fastq.gz \
-o ${ROOTDIR}/fastp/${ID}/${NAME}_R1_fastp_dedup_adapter.fastq.gz \
-O ${ROOTDIR}/fastp/${ID}/${NAME}_R2_fastp_dedup_adapter.fastq.gz \
-h ${ROOTDIR}/fastp/${ID}/${NAME}_fastp_dedup_adapter.html \
-j ${ROOTDIR}/fastp/${ID}/${NAME}_fastp_dedup_adapter.json \
--thread 16 \
--length_required 40 \
--length_limit 40 \
--dedup

#cutadapt
cutadapt -j 40 -q 10,10 -u 12 -u -10 -o ${ROOTDIR}/cutadapt/${ID}/${NAME}_R1_cutadapt_fastp_dedup_adapter.fastq.gz ${ROOTDIR}/fastp/${ID}/${NAME}_R1_fastp_dedup_adapter.fastq.gz;
cutadapt -j 40 -q 10,10 -u 10 -u -12 -o ${ROOTDIR}/cutadapt/${ID}/${NAME}_R2_cutadapt_fastp_dedup_adapter.fastq.gz ${ROOTDIR}/fastp/${ID}/${NAME}_R2_fastp_dedup_adapter.fastq.gz;
done

#STAR
print("Mapping...")
ulimit -n 65535
for NAME in ${REPEAT};
do
STAR  --runThreadN 40 \
--genomeDir  ${GENOMEDIR_STAR} \
--readFilesCommand zcat \
--readFilesIn  ${ROOTDIR}/cutadapt/${ID}/${NAME}_R2_cutadapt_fastp_dedup_adapter.fastq.gz \
--outFileNamePrefix ${ROOTDIR}/STAR/${ID}/${NAME}/${NAME}.bam \
--outSAMtype BAM SortedByCoordinate ;
done

#samtools_index_bigwig

print("Indexing and calculating depth...")
for NAME in ${REPEAT};
do
samtools index -b ${ROOTDIR}/STAR/${ID}/${NAME}/${NAME}.bamAligned.sortedByCoord.out.bam;
samtools depth -l 10 ${ROOTDIR}/find_peak/${NAME}.bamAligned.sortedByCoord.out.bam > ${ROOTDIR}/find_peak/${NAME}.10.coverage
done

exit 0
